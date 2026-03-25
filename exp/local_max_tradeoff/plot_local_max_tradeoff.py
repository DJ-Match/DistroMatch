import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import FancyBboxPatch
from scipy.stats import gmean, gstd

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 30,
    "text.latex.preamble": r"""
        \usepackage{amsmath}
        \usepackage{bm}
        \renewcommand{\familydefault}{\sfdefault}  % optional: switch to sans-serif if you want
        \renewcommand{\seriesdefault}{\bfdefault}  % <-- this sets all text to bold by default
    """
})

# === Simulate / load data (same pre-processing as before) ===
file_path = 'results_run_local_max_tradeoff.txt'
data = pd.read_csv(file_path, sep=' ')

data = data[data['k'] != 2]

baseline_per_core = data[data['epsilon'] == 0].set_index(['graph', 'k', 'cores'])[['running_time', 'matching_weight']]
baseline_all_cores = data[(data['epsilon'] == 0) & (data['cores'] == 1)].set_index(['graph', 'k'])[['running_time', 'matching_weight']]

data = data.merge(baseline_per_core, on=['graph', 'k', 'cores'], suffixes=('', '_baseline_per_core'))
data = data.merge(baseline_all_cores, on=['graph', 'k'], suffixes=('', '_baseline_all_cores'))

data['relative_running_time'] = data['running_time'] / data['running_time_baseline_per_core']
data['relative_matching_weight'] = data['matching_weight'] / data['matching_weight_baseline_per_core']

# === Separate data for fb_ and rmat_ graphs ===
fb_data = data[data['graph'].str.startswith('fb_')]
rmat_data = data[data['graph'].str.startswith('rmat_')]

# === Group data for fb_ and rmat_ ===
fb_grouped = fb_data.groupby(['epsilon', 'cores']).agg(
    running_time_mean=('relative_running_time', gmean),
    running_time_std=('relative_running_time', gstd),
    matching_weight_mean=('relative_matching_weight', gmean),
    matching_weight_std=('relative_matching_weight', gstd)
).reset_index()

rmat_grouped = rmat_data.groupby(['epsilon', 'cores']).agg(
    running_time_mean=('relative_running_time', gmean),
    running_time_std=('relative_running_time', gstd),
    matching_weight_mean=('relative_matching_weight', gmean),
    matching_weight_std=('relative_matching_weight', gstd)
).reset_index()

# === Plotting function for subplots with FancyBboxPatch legends ===
def plot_grouped_data_subplot(ax, grouped_data, title, colors, base_markers):
    core_values = sorted(grouped_data['cores'].unique())

    # Assign markers
    marker_map = {
        cores: (base_markers[i * 2 % len(base_markers)], base_markers[(i * 2 + 1) % len(base_markers)])
        for i, cores in enumerate(core_values)
    }

    for i, cores in enumerate(core_values):
        subset = grouped_data[grouped_data['cores'] == cores]
        color = colors[i % len(colors)]
        rt_marker, mw_marker = marker_map[cores]

        # Running Time
        ax.plot(subset['epsilon'], subset['running_time_mean'],
                color=color, linestyle='-', marker=rt_marker,
                markersize=14, linewidth=4,
                label=f'RT p={cores}')
        
        # Matching Weight
        ax.plot(subset['epsilon'], subset['matching_weight_mean'],
                color=color, linestyle='--', marker=mw_marker,
                markersize=14, linewidth=4,
                label=f'MW p={cores}', markevery=0.2)

    # Axis setup
    if title == 'fb':
        ax.set_ylabel('Relative Value')

    ax.set_xlabel(r'$\epsilon$')
    ax.set_xscale('log')
    ax.grid(True, which='major', linestyle='--', alpha=0.5)
    ax.grid(True, which='minor', axis='y', linestyle=':', alpha=0.3)
    ax.set_title(title)

    if title == 'rmat':
        return
    
# === Create side-by-side subplots ===
fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)

# Colors and markers
colors = ['#004488', '#DDAA33', '#BB5566']  # Adjust colors if needed
base_markers = ['o', 's', 'D', '^', 'v', '*', 'X']

# Plot fb_ data on the first subplot
plot_grouped_data_subplot(axes[0], fb_grouped, "fb", colors, base_markers)

# Plot rmat_ data on the second subplot
plot_grouped_data_subplot(axes[1], rmat_grouped, "rmat", colors, base_markers)

# Add a shared legend
core_values = sorted(fb_grouped['cores'].unique())
marker_map = {
    cores: (base_markers[i * 2 % len(base_markers)], base_markers[(i * 2 + 1) % len(base_markers)])
    for i, cores in enumerate(core_values)
}

rt_handles = [
    Line2D([0], [0], color=colors[i % len(colors)], linestyle='-', marker=marker_map[cores][0],
            markersize=14, label=f'p={cores}')
    for i, cores in enumerate(core_values)
]
mw_handles = [
    Line2D([0], [0], color=colors[i % len(colors)], linestyle='--', marker=marker_map[cores][1],
            markersize=14, label=f'p={cores}')
    for i, cores in enumerate(core_values)
]

# Add FancyBboxPatch background
legend_x = 0.25
legend_y = -0.12
legend_width = 0.75
row_height = 0.09

box = FancyBboxPatch((legend_x -0.1, legend_y - row_height +0.1),
legend_width, 2 * row_height -0.02,
boxstyle="round,pad=0.02",
transform=fig.transFigure,
edgecolor='grey', facecolor='white', linewidth=1)
box.set_clip_on(False)
fig.patches.append(box)

# Add legends for RT and MW
rt_legend = fig.legend(rt_handles, [f'p={c}' for c in core_values],
                        loc='upper left', bbox_to_anchor=(legend_x + 0.1, legend_y + row_height +0.125),
                        ncol=len(core_values), frameon=False)
fig.add_artist(rt_legend)

mw_legend = fig.legend(mw_handles, [f'p={c}' for c in core_values],
                        loc='upper left', bbox_to_anchor=(legend_x + 0.1, legend_y +0.125),
                        ncol=len(core_values), frameon=False)
fig.add_artist(mw_legend)

# Add titles for the legend rows
fig.text(legend_x -0.1 , legend_y + 0.125, "Running Time:", transform=fig.transFigure,
ha='left', va='center')
fig.text(legend_x -0.1 , legend_y + 0.125 - row_height, "Solution Weight:", transform=fig.transFigure,
ha='left', va='center')


# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 0.95])  # Leave space for the legend at the top

# Save or Show
plt.savefig("epsilon_tradeoff.pdf", dpi=500, bbox_inches='tight')
# plt.show()

d_fb=data[data['graph'].str.startswith('fb')]
d_rmat=data[data['graph'].str.startswith('rmat')]
# data = data[data['k'] == 8]
# Print reductions
fb_rt_reduction = 1 - d_fb[d_fb['epsilon'] == 0.1]['relative_running_time'].mean()
fb_mw_reduction = 1 - d_fb[d_fb['epsilon'] == 0.1]['relative_matching_weight'].mean()
rmat_rt_reduction = 1 - d_rmat[d_rmat['epsilon'] == 0.1]['relative_running_time'].mean()
rmat_mw_reduction = 1 - d_rmat[d_rmat['epsilon'] == 0.1]['relative_matching_weight'].mean()

print(f"fb: epsilon=0.1: running_time_reduction: {fb_rt_reduction:.4f}, weight_reduction: {fb_mw_reduction:.4f}")
print(f"rmat: epsilon=0.1: running_time_reduction: {rmat_rt_reduction:.4f}, weight_reduction: {rmat_mw_reduction:.4f}")
