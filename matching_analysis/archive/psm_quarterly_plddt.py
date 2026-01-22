"""
PSM Analysis: Quarterly Event Studies + PLDDT Heterogeneity
- Quarterly (3-month) aggregation for more granular time resolution
- Semester (6-month) for comparison
- PLDDT-based heterogeneity (does AlphaFold quality matter?)
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import cdist
import warnings
warnings.filterwarnings('ignore')

plt.style.use('seaborn-v0_8-whitegrid')

# =============================================================================
# 1. Load Data
# =============================================================================
print("Loading data...")
panel = pd.read_stata("../final/final_panel_CLEAN.dta")

TREATMENT_YM = 738
panel['ym_seq'] = (panel['ym'] - panel['ym'].min() + 1).astype(int)
TREATMENT_SEQ = int(TREATMENT_YM - panel['ym'].min() + 1)
panel['treated'] = (panel['num_deposits'] == 0).astype(int)

# Create quarter and semester
panel['quarter'] = ((panel['ym_seq'] - 1) // 3).astype(int) - 6  # 6 pre-treatment quarters
panel['semester'] = ((panel['ym_seq'] - 1) // 6).astype(int) - 3  # 3 pre-treatment semesters

print(f"Panel: {panel['gene_id'].nunique()} genes")
print(f"Quarters: {panel['quarter'].min()} to {panel['quarter'].max()}")
print(f"Semesters: {panel['semester'].min()} to {panel['semester'].max()}")

# Outcomes
OUTCOMES = ['n_papers', 'n_newcomer_papers', 'n_veteran_papers', 'n_top10_y']

# =============================================================================
# 2. Create Matching (Trajectory on n_papers)
# =============================================================================
print("\nCreating matches...")
pre = panel[panel['ym_seq'] < TREATMENT_SEQ].copy()

pre_agg = pre.groupby(['gene_id', 'ym_seq'])['n_papers'].sum().reset_index()
trajectory = pre_agg.pivot(index='gene_id', columns='ym_seq', values='n_papers').fillna(0)

genes = panel.groupby('gene_id').agg({
    'treated': 'first',
    'average_plddt': 'first'
}).reset_index()

treated = genes[genes['treated'] == 1].reset_index(drop=True)
control = genes[genes['treated'] == 0].reset_index(drop=True)

treated_traj = trajectory.loc[treated['gene_id']].values
control_traj = trajectory.loc[control['gene_id']].values

# Trajectory matching
print("Computing trajectory distances...")
dist_traj = cdist(treated_traj, control_traj, 'cityblock')
match_idx = np.argmin(dist_traj, axis=1)

psm_pairs = pd.DataFrame({
    'treated_id': treated['gene_id'].values,
    'control_id': control['gene_id'].values[match_idx],
    'treated_plddt': treated['average_plddt'].values,
})

# Add pre-mean for binning
pre_mean = pre.groupby('gene_id')['n_papers'].mean()
psm_pairs['treated_pre_mean'] = psm_pairs['treated_id'].map(pre_mean)

print(f"Pairs: {len(psm_pairs)}, Unique controls: {psm_pairs['control_id'].nunique()}")

# =============================================================================
# 3. Event Study Function (Generalized for any time aggregation)
# =============================================================================
def event_study(panel, pairs, outcome, time_var, ref_period):
    """
    Event study with explicit SE reporting.

    Parameters:
    - panel: full panel data
    - pairs: matched pairs DataFrame
    - outcome: outcome variable name
    - time_var: 'quarter' or 'semester'
    - ref_period: reference period (e.g., -1 for quarter before treatment)
    """
    # Aggregate to time level
    agg = panel.groupby(['gene_id', time_var])[outcome].sum().reset_index()

    # Merge treated
    t_out = agg.rename(columns={'gene_id': 'treated_id', outcome: 'y_t'})
    merged = pairs[['treated_id', 'control_id']].merge(t_out, on='treated_id')

    # Merge control
    c_out = agg.rename(columns={'gene_id': 'control_id', outcome: 'y_c'})
    merged = merged.merge(c_out, on=['control_id', time_var])

    # Pair-wise difference
    merged['diff'] = merged['y_t'] - merged['y_c']

    # Aggregate by time period
    results = merged.groupby(time_var)['diff'].agg(['mean', 'std', 'count']).reset_index()
    results.columns = [time_var, 'coef', 'std', 'n']
    results['se'] = results['std'] / np.sqrt(results['n'])

    # Normalize to reference
    ref_val = results[results[time_var] == ref_period]['coef'].values[0]
    results['coef_raw'] = results['coef']  # Keep raw for debugging
    results['coef'] = results['coef'] - ref_val
    results['ci_low'] = results['coef'] - 1.96 * results['se']
    results['ci_high'] = results['coef'] + 1.96 * results['se']

    # Pre/post stats
    pre_mask = results[time_var] < ref_period
    post_mask = results[time_var] > ref_period

    stats = {
        'pre_avg': results[pre_mask]['coef'].mean() if pre_mask.sum() > 0 else np.nan,
        'post_avg': results[post_mask]['coef'].mean() if post_mask.sum() > 0 else np.nan,
        'n_pairs': len(pairs),
        'n_controls': pairs['control_id'].nunique(),
        'ref_period': ref_period
    }

    return results, stats

# =============================================================================
# 4. Run Quarterly Event Studies
# =============================================================================
print("\n" + "="*80)
print("QUARTERLY EVENT STUDIES (3-month periods)")
print("="*80)

quarterly_results = {}
for outcome in OUTCOMES:
    coef, stats = event_study(panel, psm_pairs, outcome, 'quarter', ref_period=-1)
    quarterly_results[outcome] = {'coef': coef, 'stats': stats}
    print(f"{outcome:25s}: pre={stats['pre_avg']:+.3f}, post={stats['post_avg']:+.3f}")

# =============================================================================
# 5. Plot Quarterly Results with SE Table
# =============================================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()
colors = ['steelblue', 'orange', 'green', 'purple']
labels = ['Total Papers', 'Newcomer Papers', 'Veteran Papers', 'Top 10% Papers']

for i, outcome in enumerate(OUTCOMES):
    ax = axes[i]
    r = quarterly_results[outcome]
    coef = r['coef']
    stats = r['stats']

    # Plot with CI
    ax.fill_between(coef['quarter'], coef['ci_low'], coef['ci_high'],
                    alpha=0.3, color=colors[i])
    ax.plot(coef['quarter'], coef['coef'], 'o-', color=colors[i],
            markersize=5, linewidth=1.5)
    ax.axhline(0, color='black', linewidth=1)
    ax.axvline(-0.5, color='red', linestyle='--', linewidth=2, label='AlphaFold')

    ax.set_xlabel('Quarter (relative to Q3 2021)', fontsize=10)
    ax.set_ylabel('Avg(Treated - Control)', fontsize=10)
    ax.set_title(f'{labels[i]}\n(n={stats["n_pairs"]:,}, {stats["n_controls"]:,} controls)',
                 fontsize=11, fontweight='bold')

    # Stats box
    ax.text(0.02, 0.98, f"Pre: {stats['pre_avg']:+.2f}\nPost: {stats['post_avg']:+.2f}",
            transform=ax.transAxes, fontsize=9, va='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.suptitle('Quarterly Event Studies: Author Type Outcomes\n(PSM Trajectory Matching, normalized to Q2 2021)',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/psm_quarterly_author_outcomes.png', dpi=150, bbox_inches='tight')
print("\nSaved: figures/psm_quarterly_author_outcomes.png")

# =============================================================================
# 6. Print Detailed Coefficient Table with SE
# =============================================================================
print("\n" + "="*100)
print("QUARTERLY COEFFICIENTS WITH STANDARD ERRORS")
print("="*100)

# Get quarters
quarters = quarterly_results['n_papers']['coef']['quarter'].values

header = f"{'Quarter':>8}"
for label in labels:
    header += f" | {label[:12]:>20}"
print(header)
print("-"*100)

for q in quarters:
    row = f"{q:>8}"
    for outcome in OUTCOMES:
        coef_df = quarterly_results[outcome]['coef']
        c = coef_df[coef_df['quarter'] == q]
        if len(c) > 0:
            coef_val = c['coef'].values[0]
            se_val = c['se'].values[0]
            row += f" | {coef_val:>8.3f} ({se_val:.3f})"
        else:
            row += f" | {'N/A':>20}"
    print(row)

# =============================================================================
# 7. PLDDT Heterogeneity Analysis
# =============================================================================
print("\n" + "="*80)
print("PLDDT HETEROGENEITY ANALYSIS")
print("="*80)

# Create PLDDT bins
plddt_bins = [0, 70, 80, 90, 100]
plddt_labels = ['Low (<70)', 'Medium (70-80)', 'High (80-90)', 'Very High (90+)']

psm_pairs['plddt_bin'] = pd.cut(psm_pairs['treated_plddt'], bins=plddt_bins,
                                 labels=plddt_labels, include_lowest=True)

print("\nPairs by PLDDT bin:")
print(psm_pairs['plddt_bin'].value_counts().sort_index())

# Run event studies by PLDDT bin
plddt_results = {}
for plddt_bin in plddt_labels:
    pairs_bin = psm_pairs[psm_pairs['plddt_bin'] == plddt_bin]
    if len(pairs_bin) > 50:  # Only if enough pairs
        coef, stats = event_study(panel, pairs_bin, 'n_papers', 'semester', ref_period=-1)
        plddt_results[plddt_bin] = {'coef': coef, 'stats': stats}
        print(f"PLDDT {plddt_bin}: {stats['n_pairs']} pairs, pre={stats['pre_avg']:+.3f}, post={stats['post_avg']:+.3f}")

# =============================================================================
# 8. Plot PLDDT Heterogeneity
# =============================================================================
fig, axes = plt.subplots(1, 4, figsize=(16, 4))
colors_plddt = ['#d62728', '#ff7f0e', '#2ca02c', '#1f77b4']

for i, plddt_bin in enumerate(plddt_labels):
    ax = axes[i]
    if plddt_bin in plddt_results:
        r = plddt_results[plddt_bin]
        coef = r['coef']
        stats = r['stats']

        ax.fill_between(coef['semester'], coef['ci_low'], coef['ci_high'],
                        alpha=0.3, color=colors_plddt[i])
        ax.plot(coef['semester'], coef['coef'], 'o-', color=colors_plddt[i],
                markersize=6, linewidth=2)
        ax.axhline(0, color='black', linewidth=1)
        ax.axvline(-0.5, color='red', linestyle='--', linewidth=1.5)

        ax.set_xlabel('Semester', fontsize=10)
        ax.set_ylabel('Avg(T - C) Papers', fontsize=10)
        ax.set_title(f'pLDDT: {plddt_bin}\n(n={stats["n_pairs"]:,})', fontsize=11, fontweight='bold')

        ax.text(0.05, 0.95, f"Pre: {stats['pre_avg']:+.2f}\nPost: {stats['post_avg']:+.2f}",
                transform=ax.transAxes, fontsize=9, va='top',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

plt.suptitle('Treatment Effect by AlphaFold Prediction Quality (pLDDT)\n(Higher pLDDT = more confident structure prediction)',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/psm_plddt_heterogeneity.png', dpi=150, bbox_inches='tight')
print("\nSaved: figures/psm_plddt_heterogeneity.png")

# =============================================================================
# 9. PLDDT × Intensity Interaction
# =============================================================================
print("\n" + "="*80)
print("PLDDT × INTENSITY INTERACTION")
print("="*80)

# Create intensity bins
intensity_bins = [0, 3, 10, np.inf]
intensity_labels = ['Low (0-3)', 'Medium (3-10)', 'High (10+)']
psm_pairs['intensity_bin'] = pd.cut(psm_pairs['treated_pre_mean'], bins=intensity_bins,
                                     labels=intensity_labels, include_lowest=True)

# 2x2: High vs Low PLDDT × High vs Low Intensity
fig, axes = plt.subplots(2, 3, figsize=(15, 8))

for row, plddt_cat in enumerate(['Low (<70)', 'Very High (90+)']):
    for col, intensity_cat in enumerate(intensity_labels):
        ax = axes[row, col]

        pairs_subset = psm_pairs[(psm_pairs['plddt_bin'] == plddt_cat) &
                                  (psm_pairs['intensity_bin'] == intensity_cat)]

        if len(pairs_subset) > 30:
            coef, stats = event_study(panel, pairs_subset, 'n_papers', 'semester', ref_period=-1)

            color = colors_plddt[0] if row == 0 else colors_plddt[3]
            ax.fill_between(coef['semester'], coef['ci_low'], coef['ci_high'],
                           alpha=0.3, color=color)
            ax.plot(coef['semester'], coef['coef'], 'o-', color=color,
                   markersize=5, linewidth=1.5)
            ax.axhline(0, color='black', linewidth=0.8)
            ax.axvline(-0.5, color='red', linestyle='--', linewidth=1)

            ax.text(0.05, 0.95, f"n={stats['n_pairs']}\nPre:{stats['pre_avg']:+.2f}\nPost:{stats['post_avg']:+.2f}",
                   transform=ax.transAxes, fontsize=8, va='top',
                   bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
        else:
            ax.text(0.5, 0.5, f'n={len(pairs_subset)}\n(too few)', transform=ax.transAxes,
                   ha='center', va='center', fontsize=10)

        if row == 0:
            ax.set_title(f'{intensity_cat}', fontsize=11, fontweight='bold')
        if col == 0:
            ax.set_ylabel(f'pLDDT {plddt_cat}\nAvg(T-C)', fontweight='bold')
        if row == 1:
            ax.set_xlabel('Semester')

plt.suptitle('Treatment Effect: pLDDT × Pre-Treatment Intensity\n(Does AlphaFold quality matter more for certain genes?)',
             fontsize=13, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/psm_plddt_intensity_interaction.png', dpi=150, bbox_inches='tight')
print("Saved: figures/psm_plddt_intensity_interaction.png")

# =============================================================================
# 10. Summary Statistics
# =============================================================================
print("\n" + "="*80)
print("SUMMARY: PLDDT HETEROGENEITY")
print("="*80)
print(f"{'PLDDT Bin':<20} {'N Pairs':>10} {'Controls':>10} {'Pre-Avg':>10} {'Post-Avg':>10}")
print("-"*80)
for plddt_bin in plddt_labels:
    if plddt_bin in plddt_results:
        s = plddt_results[plddt_bin]['stats']
        print(f"{plddt_bin:<20} {s['n_pairs']:>10,} {s['n_controls']:>10,} {s['pre_avg']:>+10.3f} {s['post_avg']:>+10.3f}")

print("\n" + "="*80)
print("Done!")
