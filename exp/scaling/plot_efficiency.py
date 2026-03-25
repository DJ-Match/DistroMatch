import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gmean
import matplotlib.ticker as ticker


# matplotlib.use('tkagg')

# Use LaTeX for font rendering
# use latex for font rendering
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 28,
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

palette_algo = {v["display_name"]: v["color"] for v in algo_config.values()}
markers_algo = {v["display_name"]: v["marker"] for v in algo_config.values()}
dashes_algo = {v["display_name"]: v.get("linestyle", "") for v in algo_config.values()}


# Load the data
data = pd.read_csv('results_run_tcp_scaling.txt', sep=' ')
df = pd.DataFrame(data)

df['graph'] = df['graph'].map(lambda x: graph_names[x]['display_name'])
df['algo'] = df['algo'].map(lambda x: algo_config[x]['display_name'])

p = ['metis', 'none', 'random', 'umetis', 'uwmetis', 'ukahip', 'diststream']

# Duplicate rows for one core
for partitioning in p:
    if partitioning != 'none':
        copy = df[df['cores'] == 1].copy()
        copy['p'] = partitioning
        df = pd.concat([df, copy], ignore_index=True)

# Set a consistent color palette for partitioning techniques
# colors = ["#CC6677","#332288","#DDCC77","#117733","#88CCEE","#882255","#44AA99","#999933","#AA4499"]


def print_geomean_speedup_per_core(part='random', algo='MRep'):
    df_filtered = df[(df['p'] == part) & (df['algo'] == algo)]
    print(f"Computing geometric mean speedup per core for partitioning='{part}', algo='{algo}'")

    # Compute T1 for each (graph, k)
    t1_df = df_filtered[df_filtered['cores'] == 1].groupby(['graph', 'k']).agg(
        T1=('running_time', gmean)
    ).reset_index()

    # Merge T1 back to full data
    df_with_t1 = df_filtered.merge(t1_df, on=['graph', 'k'], how='left')
    df_with_t1['speedup'] = df_with_t1['T1'] / df_with_t1['running_time']

    # Group by cores and compute geometric mean speedup over all k and graphs
    geomean_by_core = df_with_t1.groupby('cores').agg(
        geo_speedup=('speedup', gmean),
        std_log_speedup=('speedup', lambda x: np.std(np.log(x))),
        count=('speedup', 'count')
    ).reset_index()

    # Compute confidence intervals in log space
    geomean_by_core['speedup_lower'] = geomean_by_core['geo_speedup'] * np.exp(-geomean_by_core['std_log_speedup'])
    geomean_by_core['speedup_upper'] = geomean_by_core['geo_speedup'] * np.exp(geomean_by_core['std_log_speedup'])

    # Print it
    print("\nGeometric Mean Speedup per Core (aggregated over all k):")
    print(geomean_by_core.to_string(index=False, float_format="{:.2f}".format))
    print("\n----------------------------------------------------------\n")

    return geomean_by_core

def print_runtime_ratios(part='random', algo='MRep'):
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


def plot_strong_scaling_runtime_geo_mean_algorithms(k_values, algos, partition, plot_legend=True):
    """
    Plot runtime vs number of cores for strong scaling (geometric mean over multiple k values)
    for four algorithms and a single partitioning method, in a 2x2 subplot layout.
    """
    # Create 2x2 subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharey=True, sharex=True)
    axes = axes.flatten()  # Flatten the 2D array of axes for easier iteration

    # To collect handles and labels for a shared legend
    handles = []
    labels = []

    def gmean(x):
        """Compute geometric mean."""
        return np.exp(np.mean(np.log(x)))

    for i, algo in enumerate(algos):
        ax = axes[i]
        ax.grid()
        # Filter the DataFrame for the current algorithm, partition, and the specified k values
        df_k = df[(df['algo'] == algo) & (df['p'] == partition) & (df['k'].isin(k_values))]

        # Group by graph and cores, and compute geometric mean and standard deviation of runtime
        grouped_df = df_k.groupby(['graph', 'cores']).agg(
            mean_runtime=('running_time', gmean),
            std_runtime=('running_time', lambda x: np.std(np.log(x)))
        ).reset_index()

        # Plot runtime for each graph
        for graph, group in grouped_df.groupby('graph'):
            group = group.sort_values('cores')  # Ensure cores are sorted for proper plotting
            line, = ax.plot(
                group['cores'], group['mean_runtime'], label=graph,
                marker=marker.get(graph, 'o'), color=palette.get(graph, None),
                linewidth=3, markersize=16
            )
            # ax.fill_between(
            #     group['cores'],
            #     group['mean_runtime'] * np.exp(-group['std_runtime']),
            #     group['mean_runtime'] * np.exp(group['std_runtime']),
            #     alpha=0.2, color=palette.get(graph, None)
            # )
            handles.append(line)  # Collect handles for the legend
            labels.append(graph)  # Collect labels for the legend

        # Set title and labels
        ax.set_title(f"{algo}")
        ax.set_xlabel("\# Controllers")
        if i % 2 == 0:  # Set y-label only for the leftmost plots
            ax.set_ylabel("Runtime (s)")

        # Optional: log-scale axes
        ax.set_yscale("log")  # Log scale for runtime
        ax.set_xscale("log")  # Log scale for number of cores

        # Fix y-axis ticks
        desired_ticks = [0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128]
        ax.yaxis.set_major_locator(ticker.FixedLocator(desired_ticks))
        # ax.set_yticklabels([str(t) for t in desired_ticks], fontweight='normal', fontsize=20)
        ax.set_yticklabels([r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_ticks], fontsize=24)  # Tick labels not bold

        # Fix x-axis ticks
        desired_xticks = [1, 2, 4, 8, 16, 32, 64, 128]
        ax.xaxis.set_major_locator(ticker.FixedLocator(desired_xticks))
        # ax.set_xticklabels([str(t) for t in desired_xticks], fontsize=20)
        ax.set_xticklabels([r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_xticks], fontsize=24)  # Tick labels not bold


        # Disable minor ticks
        ax.yaxis.set_minor_locator(ticker.NullLocator())
        ax.xaxis.set_minor_locator(ticker.NullLocator())

    # Hide any unused subplots (if algos < 4)
    for j in range(len(algos), len(axes)):
        fig.delaxes(axes[j])

    # Add a shared legend at the bottom
    if plot_legend:
        fig.legend(
            handles[:len(set(labels))],  # Unique handles
            labels[:len(set(labels))],  # Unique labels
            loc='upper center',
            ncol=4,
            fontsize=19,
            bbox_to_anchor=(0.5, 0.1)
        )

    plt.tight_layout(rect=[0, 0.1, 1, 1])  # Adjust layout for legend
    plt.savefig(f"strong_scaling.pdf", bbox_inches='tight', dpi=500)
    plt.close()


def plot_strong_scaling_runtime_geo_mean_algorithms_per_graph_type(k_values, algos, partition, plot_legend=True):
    """
    Plot runtime vs number of cores for strong scaling (geometric mean over multiple k values)
    in two subplots:
      - Left: graphs starting with 'fb'
      - Right: graphs starting with 'rmat'
    Lines correspond to algorithms.
    """

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True, sharex=True)

    def gmean(x):
        return np.exp(np.mean(np.log(x)))

    graph_groups = {
        "fb": lambda g: g.startswith("fb"),
        "rmat": lambda g: g.startswith("rmat"),
    }

    handles = []
    labels = []

    for ax, (group_name, graph_filter) in zip(axes, graph_groups.items()):
        ax.grid()

        # Filter relevant graphs
        relevant_graphs = [g for g in df['graph'].unique() if graph_filter(g)]

        for algo in algos:
            df_k = df[
                (df['algo'] == algo) &
                (df['p'] == partition) &
                (df['k'].isin(k_values)) &
                (df['graph'].isin(relevant_graphs))
            ]

            if df_k.empty:
                continue

            # Geometric mean over graphs and k values per core count
            grouped_df = df_k.groupby(['cores']).agg(
                mean_runtime=('running_time', gmean),
                std_runtime=('running_time', lambda x: np.std(np.log(x)))
            ).reset_index()

            grouped_df = grouped_df.sort_values('cores')

            line, = ax.plot(
                grouped_df['cores'],
                grouped_df['mean_runtime'],
                marker=markers_algo[algo],
                color=palette_algo[algo],
                linewidth=3,
                markersize=10,
                label=algo
            )

            if ax == axes[0]:
                handles.append(line)
                labels.append(algo)

        ax.set_title(f"{group_name.upper()} Graphs")
        ax.set_xlabel("\# Controllers")
        if ax == axes[0]:
            ax.set_ylabel("Runtime (s)")

        ax.set_xscale("log")
        ax.set_yscale("log")

        desired_ticks = [0.125, 0.25, 0.5, 1, 2, 4, 8, 16, 32, 64, 128]
        ax.yaxis.set_major_locator(ticker.FixedLocator(desired_ticks))
        ax.set_yticklabels(
            [r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_ticks],
            fontsize=20
        )

        desired_xticks = [1, 2, 4, 8, 16, 32, 64, 128]
        ax.xaxis.set_major_locator(ticker.FixedLocator(desired_xticks))
        ax.set_xticklabels(
            [r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_xticks],
            fontsize=20
        )

        ax.yaxis.set_minor_locator(ticker.NullLocator())
        ax.xaxis.set_minor_locator(ticker.NullLocator())

    if plot_legend:
        fig.legend(
            handles,
            labels,
            loc='upper center',
            ncol=len(algos),
            fontsize=16,
            bbox_to_anchor=(0.5, 0.15)
        )

    plt.tight_layout(rect=[0, 0.1, 1, 1])
    plt.savefig("strong_scaling.pdf", bbox_inches='tight', dpi=500)
    plt.close()

def plot_normalized_runtime_with_std(k_values, algo, partitions=['random', 'diststream']):
    """
    Plot runtime normalized by the 'random' partition runtime, with standard deviation,
    for a single algorithm and multiple k values.
    """
    # Ensure 'random' is in the partitions list
    if 'random' not in partitions:
        raise ValueError("'random' partition must be included in the partitions list.")

    # Create subplots
    fig, ax = plt.subplots(figsize=(8, 6))

    def gmean(x):
        """Compute geometric mean."""
        return np.exp(np.mean(np.log(x)))

    # Filter the DataFrame for the current algorithm and the specified k values
    df_k = df[(df['algo'] == algo) & (df['k'].isin(k_values))]

    # Group by graph, cores, and partition, and compute geometric mean and standard deviation of runtime
    grouped_df = df_k.groupby(['graph', 'cores', 'p']).agg(
        mean_runtime=('running_time', gmean),
        std_runtime=('running_time', lambda x: np.std(np.log(x)))
    ).reset_index()

    # Pivot the data to get 'random' runtime for normalization
    pivoted_df = grouped_df.pivot(index=['graph', 'cores'], columns='p', values=['mean_runtime', 'std_runtime'])
    pivoted_df.columns = ['_'.join(col).strip() for col in pivoted_df.columns.values]  # Flatten multi-index columns
    pivoted_df = pivoted_df.reset_index()

    # Plot normalized runtime for each graph
    for graph in pivoted_df['graph'].unique():
        graph_data = pivoted_df[pivoted_df['graph'] == graph]
        graph_data = graph_data.sort_values('cores')  # Ensure cores are sorted for proper plotting

        # Compute normalized runtime and standard deviation
        graph_data['normalized_diststream'] = graph_data['mean_runtime_diststream'] / graph_data['mean_runtime_random']
        graph_data['normalized_std'] = np.sqrt(
            (graph_data['std_runtime_diststream'] / graph_data['mean_runtime_diststream'])**2 +
            (graph_data['std_runtime_random'] / graph_data['mean_runtime_random'])**2
        ) * graph_data['normalized_diststream']

        # Plot normalized runtime
        line, = ax.plot(
            graph_data['cores'], graph_data['normalized_diststream'], label=graph,
            marker=marker.get(graph, 'o'), color=palette.get(graph, None),
            linewidth=3, markersize=16
        )

        # Add error bars for standard deviation
        ax.fill_between(
            graph_data['cores'],
            graph_data['normalized_diststream'] - graph_data['normalized_std'],
            graph_data['normalized_diststream'] + graph_data['normalized_std'],
            alpha=0.2, color=palette.get(graph, None)
        )

    # Set title and labels
    ax.set_title(f"Normalized Runtime (diststream / random) for {algo}")
    ax.set_xlabel("\# Controllers")
    ax.set_ylabel("Normalized Runtime")

    # Optional: log-scale axes
    ax.set_xscale("log")  # Log scale for number of cores

    # Fix x-axis ticks
    desired_xticks = [1, 2, 4, 8, 16, 32, 64, 128]
    ax.xaxis.set_major_locator(ticker.FixedLocator(desired_xticks))
    ax.set_xticklabels([str(t) for t in desired_xticks])

    # Add a legend
    ax.legend(loc='upper left', fontsize=12)

    plt.tight_layout()
    plt.savefig(f"normalized_runtime_{algo}.pdf", bbox_inches='tight', dpi=500)
    plt.close()


# plot_strong_scaling_runtime_geo_mean_algorithms(k_values=[4, 8, 16, 32], algos=['kMM','kCS', 'MRepLM','MRepS'], partition='random')
plot_strong_scaling_runtime_geo_mean_algorithms_per_graph_type(k_values=[4, 8, 16, 32], algos=['kMM','kCS', 'MRepLM','MRepS'], partition='random')

# plot_strong_scaling_runtime_geo_mean_algorithms(k_values=[4, 8, 16, 32], algos=['MRepLM','MRepS'], partition='random')

# Example usage
# plot_normalized_runtime_with_std(k_values=[32], algo='kMM')
# plot_normalized_runtime_with_std(k_values=[32], algo='MRep')
# plot_normalized_runtime_with_std(k_values=[32], algo='kCS')


print_runtime_ratios(part='random', algo='kMM')
print_runtime_ratios(part='random', algo='kCS')
print_runtime_ratios(part='random', algo='MRepLM')
print_runtime_ratios(part='random', algo='MRepS')
