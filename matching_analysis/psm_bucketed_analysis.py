"""
PSM Bucketed Analysis: Event studies by pre-treatment intensity bins
FAST VERSION - vectorized operations
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

print(f"Panel: {panel['gene_id'].nunique()} genes, Treatment at ym_seq={TREATMENT_SEQ}")

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

# Store trajectories as arrays
treated_traj = trajectory.loc[treated['gene_id']].values
control_traj = trajectory.loc[control['gene_id']].values

print(f"Treated: {len(treated)}, Control: {len(control)}")

# =============================================================================
# 3. PSM Matching (Vectorized)
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
    'treated_pre_mean': treated['pre_mean'].values,
    'control_pre_mean': control['pre_mean'].values[match_idx_mean]
})

print(f"PSM Pre-Mean: {len(psm_mean)} pairs, {psm_mean['control_id'].nunique()} unique controls")

# PSM Trajectory
print("Computing trajectory distances (this takes a moment)...")
dist_traj = cdist(treated_traj, control_traj, 'cityblock')
match_idx_traj = np.argmin(dist_traj, axis=1)

psm_traj = pd.DataFrame({
    'treated_id': treated['gene_id'].values,
    'control_id': control['gene_id'].values[match_idx_traj],
    'treated_pre_mean': treated['pre_mean'].values,
    'control_pre_mean': control['pre_mean'].values[match_idx_traj]
})

print(f"PSM Trajectory: {len(psm_traj)} pairs, {psm_traj['control_id'].nunique()} unique controls")

# =============================================================================
# 4. Bucket the Pairs
# =============================================================================
BINS = [0, 1, 3, 5, 10, 20, np.inf]
BIN_LABELS = ['0-1', '1-3', '3-5', '5-10', '10-20', '20+']

psm_mean['bin'] = pd.cut(psm_mean['treated_pre_mean'], bins=BINS, labels=BIN_LABELS, include_lowest=True)
psm_traj['bin'] = pd.cut(psm_traj['treated_pre_mean'], bins=BINS, labels=BIN_LABELS, include_lowest=True)

print("\nPairs per bin:")
print(psm_mean['bin'].value_counts().sort_index())

# =============================================================================
# 5. Event Study Function (FAST - Vectorized)
# =============================================================================
def event_study_fast(panel, pairs):
    """Fast event study using merges instead of loops."""
    # Aggregate to semester level
    sem = panel.groupby(['gene_id', 'semester'])['n_papers'].sum().reset_index()

    # Merge treated outcomes
    treated_outcomes = sem.rename(columns={'gene_id': 'treated_id', 'n_papers': 'y_treated'})
    merged = pairs[['treated_id', 'control_id']].merge(treated_outcomes, on='treated_id')

    # Merge control outcomes
    control_outcomes = sem.rename(columns={'gene_id': 'control_id', 'n_papers': 'y_control'})
    merged = merged.merge(control_outcomes, on=['control_id', 'semester'])

    # Compute pair-wise difference
    merged['diff'] = merged['y_treated'] - merged['y_control']

    # Aggregate by semester
    results = merged.groupby('semester')['diff'].agg(['mean', 'std', 'count']).reset_index()
    results.columns = ['semester', 'coef', 'std', 'n']
    results['se'] = results['std'] / np.sqrt(results['n'])

    # Normalize to semester -1
    ref = results[results['semester'] == -1]['coef'].values[0]
    results['coef'] = results['coef'] - ref
    results['ci_low'] = results['coef'] - 1.96 * results['se']
    results['ci_high'] = results['coef'] + 1.96 * results['se']

    # Stats
    pre_avg = results[results['semester'] < -1]['coef'].mean()
    post_avg = results[results['semester'] >= 0]['coef'].mean()

    stats = {
        'pre_avg': pre_avg,
        'post_avg': post_avg,
        'n_pairs': len(pairs),
        'n_treated': pairs['treated_id'].nunique(),
        'n_control': pairs['control_id'].nunique()
    }

    return results, stats

# =============================================================================
# 6. Run Event Studies by Bin
# =============================================================================
print("\nRunning event studies by bin...")

results_psm_mean = {}
results_psm_traj = {}

for bin_label in BIN_LABELS:
    pairs_mean = psm_mean[psm_mean['bin'] == bin_label]
    pairs_traj = psm_traj[psm_traj['bin'] == bin_label]

    if len(pairs_mean) > 0:
        coef, stats = event_study_fast(panel, pairs_mean)
        results_psm_mean[bin_label] = {'coef': coef, 'stats': stats}
        print(f"  PSM Pre-Mean [{bin_label}]: {stats['n_pairs']} pairs, {stats['n_control']} controls, pre={stats['pre_avg']:+.3f}, post={stats['post_avg']:+.3f}")

    if len(pairs_traj) > 0:
        coef, stats = event_study_fast(panel, pairs_traj)
        results_psm_traj[bin_label] = {'coef': coef, 'stats': stats}
        print(f"  PSM Trajectory [{bin_label}]: {stats['n_pairs']} pairs, {stats['n_control']} controls, pre={stats['pre_avg']:+.3f}, post={stats['post_avg']:+.3f}")

# =============================================================================
# 7. Plot: PSM Pre-Mean by Bin (2x3 grid)
# =============================================================================
print("\nGenerating plots...")

fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()
colors = plt.cm.viridis(np.linspace(0.2, 0.8, len(BIN_LABELS)))

for i, bin_label in enumerate(BIN_LABELS):
    ax = axes[i]
    if bin_label in results_psm_mean:
        r = results_psm_mean[bin_label]
        coef, stats = r['coef'], r['stats']

        ax.fill_between(coef['semester'], coef['ci_low'], coef['ci_high'], alpha=0.3, color=colors[i])
        ax.plot(coef['semester'], coef['coef'], 'o-', color=colors[i], markersize=6, linewidth=1.5)
        ax.axhline(0, color='black', linewidth=0.8)
        ax.axvline(-0.5, color='red', linestyle='--', linewidth=1.5)

        ax.set_xlabel('Semester')
        ax.set_ylabel('Avg(T - C)')
        ax.set_title(f'Bin: {bin_label} papers/month\n(n={stats["n_pairs"]} pairs, {stats["n_control"]} controls)', fontweight='bold')

        ax.text(0.05, 0.95, f"Pre: {stats['pre_avg']:+.2f}\nPost: {stats['post_avg']:+.2f}",
                transform=ax.transAxes, fontsize=9, va='top',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

plt.suptitle('PSM (Pre-Mean Matching) Event Studies by Intensity Bin', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/psm_premean_by_bin.png', dpi=150, bbox_inches='tight')
print("Saved: figures/psm_premean_by_bin.png")

# =============================================================================
# 8. Plot: PSM Trajectory by Bin (2x3 grid)
# =============================================================================
fig, axes = plt.subplots(2, 3, figsize=(15, 10))
axes = axes.flatten()

for i, bin_label in enumerate(BIN_LABELS):
    ax = axes[i]
    if bin_label in results_psm_traj:
        r = results_psm_traj[bin_label]
        coef, stats = r['coef'], r['stats']

        ax.fill_between(coef['semester'], coef['ci_low'], coef['ci_high'], alpha=0.3, color=colors[i])
        ax.plot(coef['semester'], coef['coef'], 'o-', color=colors[i], markersize=6, linewidth=1.5)
        ax.axhline(0, color='black', linewidth=0.8)
        ax.axvline(-0.5, color='red', linestyle='--', linewidth=1.5)

        ax.set_xlabel('Semester')
        ax.set_ylabel('Avg(T - C)')
        ax.set_title(f'Bin: {bin_label} papers/month\n(n={stats["n_pairs"]} pairs, {stats["n_control"]} controls)', fontweight='bold')

        ax.text(0.05, 0.95, f"Pre: {stats['pre_avg']:+.2f}\nPost: {stats['post_avg']:+.2f}",
                transform=ax.transAxes, fontsize=9, va='top',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))

plt.suptitle('PSM (Trajectory Matching) Event Studies by Intensity Bin', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/psm_trajectory_by_bin.png', dpi=150, bbox_inches='tight')
print("Saved: figures/psm_trajectory_by_bin.png")

# =============================================================================
# 9. Combined Comparison Plot
# =============================================================================
fig, axes = plt.subplots(2, 6, figsize=(18, 8))

for i, bin_label in enumerate(BIN_LABELS):
    # Top row: PSM Pre-Mean
    ax = axes[0, i]
    if bin_label in results_psm_mean:
        r = results_psm_mean[bin_label]
        coef, stats = r['coef'], r['stats']
        ax.fill_between(coef['semester'], coef['ci_low'], coef['ci_high'], alpha=0.3, color='steelblue')
        ax.plot(coef['semester'], coef['coef'], 'o-', color='steelblue', markersize=4, linewidth=1)
        ax.axhline(0, color='black', linewidth=0.8)
        ax.axvline(-0.5, color='red', linestyle='--', linewidth=1)
        ax.set_title(f'{bin_label}\n({stats["n_pairs"]} pairs)', fontsize=10, fontweight='bold')
        ax.text(0.05, 0.95, f"Pre:{stats['pre_avg']:+.2f}\nPost:{stats['post_avg']:+.2f}",
                transform=ax.transAxes, fontsize=7, va='top',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))
        if i == 0:
            ax.set_ylabel('PSM Pre-Mean\nAvg(T-C)', fontweight='bold')

    # Bottom row: PSM Trajectory
    ax = axes[1, i]
    if bin_label in results_psm_traj:
        r = results_psm_traj[bin_label]
        coef, stats = r['coef'], r['stats']
        ax.fill_between(coef['semester'], coef['ci_low'], coef['ci_high'], alpha=0.3, color='green')
        ax.plot(coef['semester'], coef['coef'], 'o-', color='green', markersize=4, linewidth=1)
        ax.axhline(0, color='black', linewidth=0.8)
        ax.axvline(-0.5, color='red', linestyle='--', linewidth=1)
        ax.set_xlabel('Semester', fontsize=9)
        ax.text(0.05, 0.95, f"Pre:{stats['pre_avg']:+.2f}\nPost:{stats['post_avg']:+.2f}",
                transform=ax.transAxes, fontsize=7, va='top',
                bbox=dict(boxstyle='round', facecolor='lightgreen', alpha=0.8))
        if i == 0:
            ax.set_ylabel('PSM Trajectory\nAvg(T-C)', fontweight='bold')

plt.suptitle('PSM Event Studies: Pre-Mean vs Trajectory Matching by Bin', fontsize=14, fontweight='bold')
plt.tight_layout()
plt.savefig('figures/psm_comparison_by_bin.png', dpi=150, bbox_inches='tight')
print("Saved: figures/psm_comparison_by_bin.png")

# =============================================================================
# 10. Summary Table
# =============================================================================
print("\n" + "="*100)
print("SUMMARY: PSM BY INTENSITY BIN")
print("="*100)

print(f"\n{'Bin':<8} | {'PSM Pre-Mean':<45} | {'PSM Trajectory':<45}")
print(f"{'':8} | {'N':>6} {'Ctrls':>6} {'Pre':>8} {'Post':>8} {'Reuse':>8} | {'N':>6} {'Ctrls':>6} {'Pre':>8} {'Post':>8} {'Reuse':>8}")
print("-"*100)

for bin_label in BIN_LABELS:
    row = f"{bin_label:<8} |"

    if bin_label in results_psm_mean:
        s = results_psm_mean[bin_label]['stats']
        reuse = s['n_pairs'] / s['n_control'] if s['n_control'] > 0 else 0
        row += f" {s['n_pairs']:>6} {s['n_control']:>6} {s['pre_avg']:>+8.3f} {s['post_avg']:>+8.3f} {reuse:>8.1f}x |"
    else:
        row += " "*47 + "|"

    if bin_label in results_psm_traj:
        s = results_psm_traj[bin_label]['stats']
        reuse = s['n_pairs'] / s['n_control'] if s['n_control'] > 0 else 0
        row += f" {s['n_pairs']:>6} {s['n_control']:>6} {s['pre_avg']:>+8.3f} {s['post_avg']:>+8.3f} {reuse:>8.1f}x"

    print(row)

print("\n" + "="*100)
print("Done!")
