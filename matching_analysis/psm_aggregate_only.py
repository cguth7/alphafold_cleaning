"""
PSM Aggregate Event Studies (No CEM)
Clean figure showing just PSM Pre-Mean and PSM Trajectory
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

print(f"Panel: {panel['gene_id'].nunique()} genes")

# =============================================================================
# 2. Create Matching Variables
# =============================================================================
print("Creating matching variables...")
pre = panel[panel['ym_seq'] < TREATMENT_SEQ].copy()
pre_mean = pre.groupby('gene_id')['n_papers'].mean()

pre_agg = pre.groupby(['gene_id', 'ym_seq'])['n_papers'].sum().reset_index()
trajectory = pre_agg.pivot(index='gene_id', columns='ym_seq', values='n_papers').fillna(0)

genes = panel.groupby('gene_id')['treated'].first().reset_index()
genes['pre_mean'] = genes['gene_id'].map(pre_mean)

treated = genes[genes['treated'] == 1].reset_index(drop=True)
control = genes[genes['treated'] == 0].reset_index(drop=True)

treated_traj = trajectory.loc[treated['gene_id']].values
control_traj = trajectory.loc[control['gene_id']].values

print(f"Treated: {len(treated)}, Control: {len(control)}")

# =============================================================================
# 3. PSM Matching
# =============================================================================
print("\nRunning PSM matching...")

# PSM Pre-Mean
t_vals = treated['pre_mean'].values.reshape(-1, 1)
c_vals = control['pre_mean'].values.reshape(-1, 1)
dist_mean = cdist(t_vals, c_vals, 'euclidean')
match_idx_mean = np.argmin(dist_mean, axis=1)

psm_mean = pd.DataFrame({
    'treated_id': treated['gene_id'].values,
    'control_id': control['gene_id'].values[match_idx_mean],
})

# PSM Trajectory
dist_traj = cdist(treated_traj, control_traj, 'cityblock')
match_idx_traj = np.argmin(dist_traj, axis=1)

psm_traj = pd.DataFrame({
    'treated_id': treated['gene_id'].values,
    'control_id': control['gene_id'].values[match_idx_traj],
})

print(f"PSM Pre-Mean: {len(psm_mean)} pairs, {psm_mean['control_id'].nunique()} unique controls")
print(f"PSM Trajectory: {len(psm_traj)} pairs, {psm_traj['control_id'].nunique()} unique controls")

# =============================================================================
# 4. Event Study Function
# =============================================================================
def event_study_fast(panel, pairs):
    """Fast event study using merges."""
    sem = panel.groupby(['gene_id', 'semester'])['n_papers'].sum().reset_index()

    treated_outcomes = sem.rename(columns={'gene_id': 'treated_id', 'n_papers': 'y_treated'})
    merged = pairs[['treated_id', 'control_id']].merge(treated_outcomes, on='treated_id')

    control_outcomes = sem.rename(columns={'gene_id': 'control_id', 'n_papers': 'y_control'})
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
# 5. Run Event Studies
# =============================================================================
print("\nRunning event studies...")
coef_mean, stats_mean = event_study_fast(panel, psm_mean)
coef_traj, stats_traj = event_study_fast(panel, psm_traj)

print(f"PSM Pre-Mean: pre={stats_mean['pre_avg']:+.3f}, post={stats_mean['post_avg']:+.3f}")
print(f"PSM Trajectory: pre={stats_traj['pre_avg']:+.3f}, post={stats_traj['post_avg']:+.3f}")

# =============================================================================
# 6. Plot: Side-by-side comparison
# =============================================================================
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# PSM Pre-Mean
ax = axes[0]
ax.fill_between(coef_mean['semester'], coef_mean['ci_low'], coef_mean['ci_high'],
                alpha=0.3, color='steelblue')
ax.plot(coef_mean['semester'], coef_mean['coef'], 'o-', color='steelblue',
        markersize=8, linewidth=2, label='Point estimate')
ax.axhline(0, color='black', linewidth=1)
ax.axvline(-0.5, color='red', linestyle='--', linewidth=2, label='AlphaFold Release')
ax.set_xlabel('Semester (relative to July 2021)', fontsize=11)
ax.set_ylabel('Avg(Treated - Control) Papers', fontsize=11)
ax.set_title(f'PSM: Pre-Mean Matching\n({stats_mean["n_pairs"]:,} pairs, {stats_mean["n_control"]} unique controls)',
             fontsize=12, fontweight='bold')
ax.text(0.05, 0.95, f"Pre-trend: {stats_mean['pre_avg']:+.2f}\nPost-effect: {stats_mean['post_avg']:+.2f}",
        transform=ax.transAxes, fontsize=10, va='top',
        bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.9))
ax.set_xticks(coef_mean['semester'].values)
ax.legend(loc='lower right', fontsize=9)

# PSM Trajectory
ax = axes[1]
ax.fill_between(coef_traj['semester'], coef_traj['ci_low'], coef_traj['ci_high'],
                alpha=0.3, color='green')
ax.plot(coef_traj['semester'], coef_traj['coef'], 'o-', color='green',
        markersize=8, linewidth=2, label='Point estimate')
ax.axhline(0, color='black', linewidth=1)
ax.axvline(-0.5, color='red', linestyle='--', linewidth=2, label='AlphaFold Release')
ax.set_xlabel('Semester (relative to July 2021)', fontsize=11)
ax.set_ylabel('Avg(Treated - Control) Papers', fontsize=11)
ax.set_title(f'PSM: Trajectory Matching\n({stats_traj["n_pairs"]:,} pairs, {stats_traj["n_control"]:,} unique controls)',
             fontsize=12, fontweight='bold')
ax.text(0.05, 0.95, f"Pre-trend: {stats_traj['pre_avg']:+.2f}\nPost-effect: {stats_traj['post_avg']:+.2f}",
        transform=ax.transAxes, fontsize=10, va='top',
        bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.9))
ax.set_xticks(coef_traj['semester'].values)
ax.legend(loc='lower right', fontsize=9)

plt.suptitle('Aggregate PSM Event Studies: Pair-wise Differences\n(Normalized to semester -1)',
             fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/psm_aggregate_comparison.png', dpi=150, bbox_inches='tight')
print("\nSaved: figures/psm_aggregate_comparison.png")

# =============================================================================
# 7. Summary stats
# =============================================================================
print("\n" + "="*70)
print("AGGREGATE PSM SUMMARY")
print("="*70)
print(f"{'Method':<20} {'Pairs':>8} {'Controls':>10} {'Pre-Avg':>10} {'Post-Avg':>10}")
print("-"*70)
print(f"{'PSM Pre-Mean':<20} {stats_mean['n_pairs']:>8,} {stats_mean['n_control']:>10} {stats_mean['pre_avg']:>+10.3f} {stats_mean['post_avg']:>+10.3f}")
print(f"{'PSM Trajectory':<20} {stats_traj['n_pairs']:>8,} {stats_traj['n_control']:>10,} {stats_traj['pre_avg']:>+10.3f} {stats_traj['post_avg']:>+10.3f}")
print("="*70)
