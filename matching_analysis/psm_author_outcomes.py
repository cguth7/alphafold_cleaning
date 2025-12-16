"""
PSM Event Studies for Author Type Outcomes
- n_newcomer_papers: Papers by authors new to this gene
- n_veteran_papers: Papers by authors who previously published on this gene
- n_top10_y: Papers by top 10% most productive authors (yearly measure)
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
panel['semester'] = ((panel['ym_seq'] - 1) // 6).astype(int) - 3

# Define outcomes of interest
OUTCOMES = {
    'n_papers': 'Total Papers',
    'n_newcomer_papers': 'Newcomer Papers',
    'n_veteran_papers': 'Veteran Papers',
    'n_top10_y': 'Top 10% Author Papers'
}

print(f"Panel: {panel['gene_id'].nunique()} genes")
print(f"Outcomes: {list(OUTCOMES.keys())}")

# Check data
for out in OUTCOMES:
    print(f"  {out}: mean={panel[out].mean():.2f}, max={panel[out].max()}")

# =============================================================================
# 2. Create Matching Variables (using n_papers trajectory)
# =============================================================================
print("\nCreating matching variables...")
pre = panel[panel['ym_seq'] < TREATMENT_SEQ].copy()

pre_agg = pre.groupby(['gene_id', 'ym_seq'])['n_papers'].sum().reset_index()
trajectory = pre_agg.pivot(index='gene_id', columns='ym_seq', values='n_papers').fillna(0)

genes = panel.groupby('gene_id')['treated'].first().reset_index()

treated = genes[genes['treated'] == 1].reset_index(drop=True)
control = genes[genes['treated'] == 0].reset_index(drop=True)

treated_traj = trajectory.loc[treated['gene_id']].values
control_traj = trajectory.loc[control['gene_id']].values

print(f"Treated: {len(treated)}, Control: {len(control)}")

# =============================================================================
# 3. PSM Trajectory Matching (using n_papers)
# =============================================================================
print("\nRunning PSM trajectory matching on n_papers...")
dist_traj = cdist(treated_traj, control_traj, 'cityblock')
match_idx_traj = np.argmin(dist_traj, axis=1)

psm_pairs = pd.DataFrame({
    'treated_id': treated['gene_id'].values,
    'control_id': control['gene_id'].values[match_idx_traj],
})

print(f"PSM Trajectory: {len(psm_pairs)} pairs, {psm_pairs['control_id'].nunique()} unique controls")

# =============================================================================
# 4. Event Study Function (generalized for any outcome)
# =============================================================================
def event_study_outcome(panel, pairs, outcome_var):
    """Event study for a specific outcome variable."""
    sem = panel.groupby(['gene_id', 'semester'])[outcome_var].sum().reset_index()

    treated_outcomes = sem.rename(columns={'gene_id': 'treated_id', outcome_var: 'y_treated'})
    merged = pairs[['treated_id', 'control_id']].merge(treated_outcomes, on='treated_id')

    control_outcomes = sem.rename(columns={'gene_id': 'control_id', outcome_var: 'y_control'})
    merged = merged.merge(control_outcomes, on=['control_id', 'semester'])

    merged['diff'] = merged['y_treated'] - merged['y_control']

    results = merged.groupby('semester')['diff'].agg(['mean', 'std', 'count']).reset_index()
    results.columns = ['semester', 'coef', 'std', 'n']
    results['se'] = results['std'] / np.sqrt(results['n'])

    ref = results[results['semester'] == -1]['coef'].values[0]
    results['coef'] = results['coef'] - ref
    results['ci_low'] = results['coef'] - 1.96 * results['se']
    results['ci_high'] = results['coef'] + 1.96 * results['se']

    pre_avg = results[results['semester'] < -1]['coef'].mean()
    post_avg = results[results['semester'] >= 0]['coef'].mean()

    stats = {
        'pre_avg': pre_avg,
        'post_avg': post_avg,
        'n_pairs': len(pairs),
        'n_control': pairs['control_id'].nunique()
    }

    return results, stats

# =============================================================================
# 5. Run Event Studies for All Outcomes
# =============================================================================
print("\nRunning event studies for all outcomes...")

results = {}
for outcome, label in OUTCOMES.items():
    coef, stats = event_study_outcome(panel, psm_pairs, outcome)
    results[outcome] = {'coef': coef, 'stats': stats, 'label': label}
    print(f"  {label}: pre={stats['pre_avg']:+.3f}, post={stats['post_avg']:+.3f}")

# =============================================================================
# 6. Plot: 2x2 Grid of Outcomes
# =============================================================================
fig, axes = plt.subplots(2, 2, figsize=(12, 10))
axes = axes.flatten()

colors = ['steelblue', 'orange', 'green', 'purple']

for i, (outcome, r) in enumerate(results.items()):
    ax = axes[i]
    coef = r['coef']
    stats = r['stats']
    label = r['label']

    ax.fill_between(coef['semester'], coef['ci_low'], coef['ci_high'],
                    alpha=0.3, color=colors[i])
    ax.plot(coef['semester'], coef['coef'], 'o-', color=colors[i],
            markersize=8, linewidth=2)
    ax.axhline(0, color='black', linewidth=1)
    ax.axvline(-0.5, color='red', linestyle='--', linewidth=2)

    ax.set_xlabel('Semester (relative to July 2021)', fontsize=11)
    ax.set_ylabel('Avg(Treated - Control)', fontsize=11)
    ax.set_title(f'{label}\n({stats["n_pairs"]:,} pairs)', fontsize=12, fontweight='bold')

    ax.text(0.05, 0.95, f"Pre: {stats['pre_avg']:+.2f}\nPost: {stats['post_avg']:+.2f}",
            transform=ax.transAxes, fontsize=10, va='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))

    ax.set_xticks(coef['semester'].values)

plt.suptitle('PSM Event Studies by Outcome Type\n(Trajectory Matching on n_papers)',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/psm_author_outcomes.png', dpi=150, bbox_inches='tight')
print("\nSaved: figures/psm_author_outcomes.png")

# =============================================================================
# 7. Plot: Overlay All Outcomes on Single Plot
# =============================================================================
fig, ax = plt.subplots(figsize=(10, 6))

for i, (outcome, r) in enumerate(results.items()):
    coef = r['coef']
    label = r['label']
    ax.plot(coef['semester'], coef['coef'], 'o-', color=colors[i],
            markersize=6, linewidth=2, label=label)

ax.axhline(0, color='black', linewidth=1)
ax.axvline(-0.5, color='red', linestyle='--', linewidth=2, label='AlphaFold Release')
ax.set_xlabel('Semester (relative to July 2021)', fontsize=12)
ax.set_ylabel('Avg(Treated - Control)', fontsize=12)
ax.set_title('PSM Event Studies: All Outcomes Compared', fontsize=14, fontweight='bold')
ax.legend(loc='lower left', fontsize=10)
ax.set_xticks(coef['semester'].values)

plt.tight_layout()
plt.savefig('figures/psm_outcomes_overlay.png', dpi=150, bbox_inches='tight')
print("Saved: figures/psm_outcomes_overlay.png")

# =============================================================================
# 8. Bucketed Analysis for Author Outcomes
# =============================================================================
print("\nRunning bucketed analysis for author outcomes...")

# Add pre-mean to pairs for bucketing
pre_mean = pre.groupby('gene_id')['n_papers'].mean()
psm_pairs['treated_pre_mean'] = psm_pairs['treated_id'].map(pre_mean)

BINS = [0, 1, 3, 5, 10, 20, np.inf]
BIN_LABELS = ['0-1', '1-3', '3-5', '5-10', '10-20', '20+']
psm_pairs['bin'] = pd.cut(psm_pairs['treated_pre_mean'], bins=BINS, labels=BIN_LABELS, include_lowest=True)

# Focus on key outcomes: newcomer vs veteran
KEY_OUTCOMES = ['n_newcomer_papers', 'n_veteran_papers']

fig, axes = plt.subplots(2, 6, figsize=(18, 8))

for row, outcome in enumerate(KEY_OUTCOMES):
    for col, bin_label in enumerate(BIN_LABELS):
        ax = axes[row, col]
        pairs_bin = psm_pairs[psm_pairs['bin'] == bin_label]

        if len(pairs_bin) > 0:
            coef, stats = event_study_outcome(panel, pairs_bin, outcome)

            color = 'orange' if outcome == 'n_newcomer_papers' else 'green'
            ax.fill_between(coef['semester'], coef['ci_low'], coef['ci_high'],
                           alpha=0.3, color=color)
            ax.plot(coef['semester'], coef['coef'], 'o-', color=color,
                   markersize=4, linewidth=1.5)
            ax.axhline(0, color='black', linewidth=0.8)
            ax.axvline(-0.5, color='red', linestyle='--', linewidth=1)

            ax.text(0.05, 0.95, f"Pre:{stats['pre_avg']:+.2f}\nPost:{stats['post_avg']:+.2f}",
                   transform=ax.transAxes, fontsize=7, va='top',
                   bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

            if row == 0:
                ax.set_title(f'{bin_label}\n({stats["n_pairs"]} pairs)', fontsize=10, fontweight='bold')
            if col == 0:
                label = 'Newcomer' if outcome == 'n_newcomer_papers' else 'Veteran'
                ax.set_ylabel(f'{label}\nAvg(T-C)', fontweight='bold')
            if row == 1:
                ax.set_xlabel('Semester', fontsize=9)

plt.suptitle('PSM Event Studies: Newcomer vs Veteran Papers by Intensity Bin',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/psm_newcomer_veteran_by_bin.png', dpi=150, bbox_inches='tight')
print("Saved: figures/psm_newcomer_veteran_by_bin.png")

# =============================================================================
# 9. Summary Table
# =============================================================================
print("\n" + "="*80)
print("SUMMARY: PSM AUTHOR OUTCOMES (AGGREGATE)")
print("="*80)
print(f"{'Outcome':<25} {'Pre-Avg':>10} {'Post-Avg':>10} {'Change':>10}")
print("-"*80)
for outcome, r in results.items():
    s = r['stats']
    change = s['post_avg'] - s['pre_avg']
    print(f"{r['label']:<25} {s['pre_avg']:>+10.3f} {s['post_avg']:>+10.3f} {change:>+10.3f}")
print("="*80)

print("\n" + "="*80)
print("INTERPRETATION")
print("="*80)
print("""
Key Questions:
1. Did AlphaFold bring NEW researchers to previously understudied genes? (newcomer effect)
2. Did existing researchers increase their output? (veteran effect)
3. Did top researchers shift their attention? (top 10% effect)

Compare newcomer vs veteran effects across bins to understand the mechanism.
""")
