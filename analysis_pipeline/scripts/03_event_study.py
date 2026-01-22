"""
03_event_study.py - Event Study Regressions

This script runs weighted event study regressions replicating Danilo's approach.

Four specifications per outcome:
1. LEVELS: Y ~ treated×period FE + gene FE + semester FE
2. TREND: Y ~ treated×period FE + gene FE + treated×linear_trend
3. DY: ΔY ~ treated×period FE + gene FE + semester FE
4. DASINH: Δasinh(Y) ~ treated×period FE + gene FE + semester FE

The "good" specs (flat pre-trends) are typically DY and DASINH.

Output:
    - outputs/*.csv: Event study coefficients
    - figures/*.png: Event study plots
"""

import pandas as pd
import numpy as np
from pathlib import Path
import statsmodels.api as sm
from linearmodels.panel import PanelOLS
import warnings
warnings.filterwarnings('ignore')

# Paths
BASE_DIR = Path(__file__).parent.parent.parent
DATA_DIR = BASE_DIR / "analysis_pipeline" / "data"
OUTPUT_DIR = BASE_DIR / "analysis_pipeline" / "outputs"
FIGURE_DIR = BASE_DIR / "analysis_pipeline" / "figures"

# Configuration
OUTCOMES = ['n_papers', 'n_newcomer_papers', 'n_veteran_papers', 'n_top10_y']
BASE_PERIOD = -1  # Base period for event study (omitted category)


def load_semester_panel():
    """Load the semester-level matched panel."""
    path = DATA_DIR / "matched_panel_semester.parquet"
    print(f"Loading semester panel from {path}")
    df = pd.read_parquet(path)
    print(f"  Loaded {len(df):,} observations")
    return df


def run_event_study_ols(df, outcome, weight_col='cem_weight', base_period=BASE_PERIOD):
    """
    Run event study regression using OLS with clustered standard errors.

    Model: Y_it = Σ_k β_k (treated_i × 1[rel_sem=k]) + α_i + γ_t + ε_it

    Returns DataFrame with coefficients, SEs, and CIs for each period.
    """
    df = df.copy()

    # Drop missing outcomes
    df = df.dropna(subset=[outcome])

    if len(df) == 0:
        return None

    # Get all relative semester values
    periods = sorted(df['rel_semester'].unique())

    # Create interaction dummies (treated × period)
    for k in periods:
        df[f'treat_sem_{k}'] = ((df['treated'] == 1) & (df['rel_semester'] == k)).astype(int)

    # Remove base period (will be omitted)
    interaction_cols = [f'treat_sem_{k}' for k in periods if k != base_period]

    # Create gene and semester dummies
    df['gene_fe'] = pd.Categorical(df['gene_id'])
    df['sem_fe'] = pd.Categorical(df['semester'])

    # Prepare regression data
    # We'll use statsmodels with absorbed fixed effects approach
    # First demean by gene and time

    # Gene means
    gene_means = df.groupby('gene_id')[outcome].transform('mean')
    # Time means
    time_means = df.groupby('semester')[outcome].transform('mean')
    # Overall mean
    overall_mean = df[outcome].mean()

    # Demean outcome (within transformation for two-way FE)
    df['y_demeaned'] = df[outcome] - gene_means - time_means + overall_mean

    # Also demean the interaction terms
    for col in interaction_cols:
        gene_m = df.groupby('gene_id')[col].transform('mean')
        time_m = df.groupby('semester')[col].transform('mean')
        overall_m = df[col].mean()
        df[f'{col}_dm'] = df[col] - gene_m - time_m + overall_m

    interaction_cols_dm = [f'{col}_dm' for col in interaction_cols]

    # Weighted regression
    X = df[interaction_cols_dm].values
    y = df['y_demeaned'].values
    w = df[weight_col].values

    # Weighted least squares
    X_weighted = X * np.sqrt(w)[:, np.newaxis]
    y_weighted = y * np.sqrt(w)

    # Add constant (absorbed but needed for statsmodels)
    X_with_const = sm.add_constant(X_weighted)

    try:
        # Fit model
        model = sm.OLS(y_weighted, X_with_const)
        results = model.fit(cov_type='cluster', cov_kwds={'groups': df['gene_id']})

        # Extract coefficients (skip constant)
        coefs = results.params[1:]
        ses = results.bse[1:]

    except Exception as e:
        print(f"  Regression failed: {e}")
        return None

    # Build results DataFrame
    results_df = []
    for i, k in enumerate([p for p in periods if p != base_period]):
        results_df.append({
            'period': k,
            'coef': coefs[i],
            'se': ses[i],
            'ci_low': coefs[i] - 1.96 * ses[i],
            'ci_high': coefs[i] + 1.96 * ses[i]
        })

    # Add base period with zeros
    results_df.append({
        'period': base_period,
        'coef': 0.0,
        'se': 0.0,
        'ci_low': 0.0,
        'ci_high': 0.0
    })

    results_df = pd.DataFrame(results_df).sort_values('period')

    return results_df


def run_event_study_simple(df, outcome, weight_col='cem_weight', base_period=BASE_PERIOD):
    """
    Simpler event study using period-by-period DiD.

    For each period k, estimate:
        E[Y|treated, period=k] - E[Y|treated, period=base] -
        (E[Y|control, period=k] - E[Y|control, period=base])

    This is more robust and easier to interpret.
    """
    df = df.copy()
    df = df.dropna(subset=[outcome])

    if len(df) == 0:
        return None

    periods = sorted(df['rel_semester'].unique())

    # Get baseline means
    base_treated = df[(df['rel_semester'] == base_period) & (df['treated'] == 1)]
    base_control = df[(df['rel_semester'] == base_period) & (df['treated'] == 0)]

    if len(base_treated) == 0 or len(base_control) == 0:
        print(f"  No data in base period {base_period}")
        return None

    base_diff = (
        np.average(base_treated[outcome], weights=base_treated[weight_col]) -
        np.average(base_control[outcome], weights=base_control[weight_col])
    )

    results = []
    for k in periods:
        treated_k = df[(df['rel_semester'] == k) & (df['treated'] == 1)]
        control_k = df[(df['rel_semester'] == k) & (df['treated'] == 0)]

        if len(treated_k) == 0 or len(control_k) == 0:
            continue

        # Weighted means
        mean_treated = np.average(treated_k[outcome], weights=treated_k[weight_col])
        mean_control = np.average(control_k[outcome], weights=control_k[weight_col])

        # DiD coefficient
        diff_k = mean_treated - mean_control
        coef = diff_k - base_diff if k != base_period else 0.0

        # Bootstrap SE (simple)
        n_boot = 200
        boot_coefs = []
        for _ in range(n_boot):
            # Sample genes with replacement
            genes_t = treated_k['gene_id'].unique()
            genes_c = control_k['gene_id'].unique()

            boot_genes_t = np.random.choice(genes_t, size=len(genes_t), replace=True)
            boot_genes_c = np.random.choice(genes_c, size=len(genes_c), replace=True)

            boot_t = treated_k[treated_k['gene_id'].isin(boot_genes_t)]
            boot_c = control_k[control_k['gene_id'].isin(boot_genes_c)]

            if len(boot_t) > 0 and len(boot_c) > 0:
                m_t = np.average(boot_t[outcome], weights=boot_t[weight_col])
                m_c = np.average(boot_c[outcome], weights=boot_c[weight_col])
                boot_coefs.append((m_t - m_c) - base_diff)

        se = np.std(boot_coefs) if boot_coefs else 0.0

        results.append({
            'period': k,
            'coef': coef,
            'se': se,
            'ci_low': coef - 1.96 * se,
            'ci_high': coef + 1.96 * se,
            'n_treated': len(treated_k),
            'n_control': len(control_k)
        })

    return pd.DataFrame(results).sort_values('period')


def plot_event_study(results_df, title, output_path, base_period=BASE_PERIOD):
    """Create event study plot."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 6))

    periods = results_df['period'].values
    coefs = results_df['coef'].values
    ci_low = results_df['ci_low'].values
    ci_high = results_df['ci_high'].values

    # Plot confidence intervals
    ax.fill_between(periods, ci_low, ci_high, alpha=0.2, color='blue')

    # Plot coefficients
    ax.plot(periods, coefs, 'o-', color='blue', markersize=8)

    # Reference lines
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.axvline(x=base_period, color='red', linestyle='--', linewidth=1, alpha=0.7)

    # Labels
    ax.set_xlabel('Semesters Relative to AlphaFold2 (base = -1)', fontsize=12)
    ax.set_ylabel('Treatment Effect', fontsize=12)
    ax.set_title(title, fontsize=14)

    # Grid
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"  Saved figure: {output_path}")


def main():
    """Run event study analysis."""
    print("=" * 60)
    print("EVENT STUDY ANALYSIS")
    print("=" * 60)

    # Ensure output directories exist
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    # Load data
    df = load_semester_panel()

    # Run event studies for each outcome and specification
    specifications = [
        ('LEVELS', lambda y: y),
        ('DY', lambda y: f'D_{y}'),
        ('DASINH', lambda y: f'D_asinh_{y}'),
    ]

    all_results = {}

    for outcome in OUTCOMES:
        print(f"\n{'='*40}")
        print(f"Outcome: {outcome}")
        print('='*40)

        for spec_name, get_var in specifications:
            var = get_var(outcome)

            if var not in df.columns:
                print(f"  Skipping {spec_name}: {var} not found")
                continue

            print(f"\n  Running {spec_name} specification...")

            # Run event study
            results = run_event_study_simple(df, var)

            if results is None:
                print(f"  Failed to run {spec_name}")
                continue

            # Save coefficients
            csv_path = OUTPUT_DIR / f"eventstudy_{spec_name}_{outcome}.csv"
            results.to_csv(csv_path, index=False)
            print(f"  Saved: {csv_path}")

            # Create plot
            title = f"{outcome} - {spec_name} (Semester)"
            fig_path = FIGURE_DIR / f"eventstudy_{spec_name}_{outcome}.png"
            plot_event_study(results, title, fig_path)

            # Store results
            all_results[(outcome, spec_name)] = results

            # Print key coefficients
            pre_coefs = results[results['period'] < 0]['coef']
            post_coefs = results[results['period'] > 0]['coef']

            print(f"  Pre-period mean coef: {pre_coefs.mean():.4f}")
            print(f"  Post-period mean coef: {post_coefs.mean():.4f}")

    print("\n" + "=" * 60)
    print("EVENT STUDY ANALYSIS COMPLETE")
    print("=" * 60)

    return all_results


if __name__ == "__main__":
    main()
