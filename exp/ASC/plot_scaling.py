import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gmean
from matplotlib.ticker import ScalarFormatter
import matplotlib.ticker as ticker

# matplotlib.use('tkagg')

# Use LaTeX for font rendering
# use latex for font rendering
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


# Load YAML config
with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)

graph_names = config.get("graph_names", {})
partition_config = config.get("partition_config", {})
algo_config = config.get("algo_config", {})

# Extract styling from config
palette = {v["display_name"]: v["color"] for v in graph_names.values()}
marker = {v["display_name"]: v["marker"] for v in graph_names.values()}





def plot_strong_scaling_runtime(k,zen_version =2, part='random', algo='CPal'):
    """
    Plot runtime vs number of cores for strong scaling (fixed problem size).
    """

    # Load the data
    data = pd.read_csv(f'results_zen{zen_version}.txt', sep=' ')
    df = pd.DataFrame(data)

    df['graph'] = df['graph'].map(lambda x: graph_names[x]['display_name'])
    df['algo'] = df['algo'].map(lambda x: algo_config[x]['display_name'])


    p = ['metis', 'none', 'random', 'umetis', 'uwmetis', 'ukahip', 'diststream']

    # Duplicate rows for one core
    for partitioning in p:
        if partitioning != 'random':
            copy = df[df['cores'] == 1].copy()
            copy['p'] = partitioning
            df = pd.concat([df, copy], ignore_index=True)

    # Set a consistent color palette for partitioning techniques
    # colors = ["#CC6677","#332288","#DDCC77","#117733","#88CCEE","#882255","#44AA99","#999933","#AA4499"]

    df_k = df[(df['k'] == k) & (df['p'] == part) & (df['algo'] == algo)]

    # Group by graph and cores
    grouped_df = df_k.groupby(['graph', 'cores']).agg(
        mean_runtime=('running_time', 'mean'),
        std_runtime=('running_time', 'std')
    ).reset_index()

    # Plotting
    plt.figure(figsize=(14,5))
    for graph, group in grouped_df.groupby('graph'):
        group = group.sort_values('cores')
        plt.plot(group['cores'], group['mean_runtime'], label=graph,
                 marker=marker.get(graph, 'o'), color=palette.get(graph, None),
                 linewidth=3, markersize=16)
        plt.fill_between(group['cores'],
                         group['mean_runtime'] - group['std_runtime'],
                         group['mean_runtime'] + group['std_runtime'],
                         alpha=0.2, color=palette.get(graph, None))
        


    plt.xlabel("\# Controllers")
    plt.ylabel("Runtime (s)")
    plt.legend(
        loc='center left',  # Place the legend to the left of the anchor point
        bbox_to_anchor=(1.05, 0.4),  # Anchor the legend outside the plot (right side, vertically centered)
        ncol=1,  # Single column for the legend
        fontsize=24  # Font size for the legend
    )
    # plt.xscale("log", base=2)  # optional: log-scale x-axis
    ax = plt.gca()
    ax.set_yscale("log")      # optional: log-scale y-axis

    desired_ticks = [0.25,0.5,1,2, 4, 8]
    ax.yaxis.set_major_locator(ticker.FixedLocator(desired_ticks))
    # ax.set_yticklabels([str(t) for t in desired_ticks], fontweight='normal')
    # ax.set_ylim(bottom=2)
    ax.set_yticklabels([r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_ticks], fontsize=24)  # Tick labels not bold


    # Disable minor ticks completely
    ax.yaxis.set_minor_locator(ticker.NullLocator())

    # If x-axis is log scale
    ax.set_xscale("log")

    # Fix x-axis ticks to specific values
    desired_xticks = [1, 2, 4, 8, 16, 32, 64,128]
    ax.xaxis.set_major_locator(ticker.FixedLocator(desired_xticks))
    # ax.set_xticklabels([str(t) for t in desired_xticks])
    ax.set_xticklabels([r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_xticks], fontsize=24)  # Tick labels not bold


    # Disable minor ticks on x-axis
    ax.xaxis.set_minor_locator(ticker.NullLocator())

    plt.grid()
    plt.tight_layout(rect=[0, 0, 0.85, 1])  # Leave space on the right for the legend

    plt.savefig(f"ASC_zen{zen_version}_k{k}_{algo}_{part}.pdf", dpi=500,bbox_inches='tight')
    plt.close()

def plot_strong_scaling_runtime_side_by_side(k, part='random', algo='CPal'):
    """
    Plot runtime vs number of cores for strong scaling (fixed problem size) for Zen2 and Zen3 side by side.
    """

    # Create subplots for Zen2 and Zen3
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for i, zen_version in enumerate([2, 3]):
        ax = axes[i]

        # Load the data for the current Zen version
        data = pd.read_csv(f'results_zen{zen_version}.txt', sep=' ')
        df = pd.DataFrame(data)

        df['graph'] = df['graph'].map(lambda x: graph_names[x]['display_name'])
        df['algo'] = df['algo'].map(lambda x: algo_config[x]['display_name'])

        p = ['metis', 'none', 'random', 'umetis', 'uwmetis', 'ukahip', 'diststream']

        # Duplicate rows for one core
        for partitioning in p:
            if partitioning != 'random':
                copy = df[df['cores'] == 1].copy()
                copy['p'] = partitioning
                df = pd.concat([df, copy], ignore_index=True)

        # Filter the data for the given parameters
        df_k = df[(df['k'] == k) & (df['p'] == part) & (df['algo'] == algo)]

        # Group by graph and cores
        grouped_df = df_k.groupby(['graph', 'cores']).agg(
            mean_runtime=('running_time', 'mean'),
            std_runtime=('running_time', 'std')
        ).reset_index()

        # Track which labels have already been added to the legend
        added_labels = set()

        # Plot runtime for each graph
        for graph, group in grouped_df.groupby('graph'):
            group = group.sort_values('cores')
            line, = ax.plot(group['cores'], group['mean_runtime'], label=graph,
                            marker=marker.get(graph, 'o'), color=palette.get(graph, None),
                            linewidth=3, markersize=16)
            ax.fill_between(group['cores'],
                            group['mean_runtime'] - group['std_runtime'],
                            group['mean_runtime'] + group['std_runtime'],
                            alpha=0.2, color=palette.get(graph, None))

            # Add to legend only if the label hasn't been added yet
            if graph not in added_labels:
                added_labels.add(graph)

        # Set axis labels and title
        ax.set_xlabel("\# Controllers")
        if i == 0:
            ax.set_ylabel("Runtime (s)")
        ax.set_title(f"Zen{zen_version}")

        # Set log scale for axes
        ax.set_yscale("log")
        ax.set_xscale("log")

        # Fix y-axis ticks
        desired_ticks = [0.25, 0.5, 1, 2, 4, 8]
        ax.yaxis.set_major_locator(ticker.FixedLocator(desired_ticks))
        ax.set_yticklabels([str(t) for t in desired_ticks], fontweight='normal')

        # Disable minor ticks on y-axis
        ax.yaxis.set_minor_locator(ticker.NullLocator())

        # Fix x-axis ticks
        desired_xticks = [1, 2, 4, 8, 16, 32, 64]
        ax.xaxis.set_major_locator(ticker.FixedLocator(desired_xticks))
        ax.set_xticklabels([str(t) for t in desired_xticks])

        # Disable minor ticks on x-axis
        ax.xaxis.set_minor_locator(ticker.NullLocator())

        ax.grid()

    added_labels = set()
    handles, labels = [], []
    for ax in axes:
        for handle, label in zip(*ax.get_legend_handles_labels()):
            if label not in added_labels:
                handles.append(handle)
                labels.append(label)
                added_labels.add(label)

    # Add a shared legend at the bottom
    fig.legend(
        handles, labels,
        loc='upper center',
        ncol=4,
        fontsize=20,
        bbox_to_anchor=(0.55, 0.05)
    )

    plt.tight_layout(rect=[0, 0, 1, 0.95])  # Adjust layout for legend
    plt.savefig(f"ASC_zen2_zen3_k{k}_{algo}_{part}.pdf", dpi=500, bbox_inches='tight')
    plt.close()

def print_runtime_ratios(part='random', zen_version=3, algo='kCS'):
    data = pd.read_csv(f'results_zen{zen_version}.txt', sep=' ')
    df = pd.DataFrame(data)
    df_filtered = df[(df['p'] == part) & (df['algo'] == algo)]
    print(f"Computing runtime ratios for partitioning='{part}', algo='{algo}'")

    # Compute runtime ratios for each pair of cores (e.g., 32/64, 16/32, ..., and 4/1)
    runtime_ratios = []
    unique_cores = sorted(df_filtered['cores'].unique())  # Sort cores in ascending order

    for i in range(1, len(unique_cores)):
        core_small = unique_cores[i - 1]
        core_large = unique_cores[i]

        # if core_large == 2 * core_small:  # Ensure the pair is in the form (x, 2x)
        cores_small = df_filtered[df_filtered['cores'] == core_small]
        cores_large = df_filtered[df_filtered['cores'] == core_large]

        if not cores_small.empty and not cores_large.empty:
            runtime_small = gmean(cores_small['running_time'])
            runtime_large = gmean(cores_large['running_time'])
            runtime_ratio = runtime_large/runtime_small
            runtime_ratios.append((core_small, core_large, runtime_ratio))
        else:
            runtime_ratios.append((core_small, core_large, None))  # Append None if data is insufficient

    # Print runtime ratios
    print("\nRuntime Ratios (for each pair of cores):")
    for core_small, core_large, ratio in runtime_ratios:
        if ratio is not None:
            print(f"Core {core_large}/{core_small}: Runtime ratio = {ratio:.2f}")
        else:
            print(f"Core {core_large}/{core_small}: Insufficient data to compute runtime ratio.")

    print("\n----------------------------------------------------------\n")

    return runtime_ratios


# plot_strong_scaling_runtime(8, 2, part='random',algo='kCS')
plot_strong_scaling_runtime(8, 3, part='random',algo='kCS')

# plot_strong_scaling_runtime_side_by_side(8, part='random', algo='kCS')

print_runtime_ratios()