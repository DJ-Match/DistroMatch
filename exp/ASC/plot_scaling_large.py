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





def plot_strong_scaling_runtime(k, part='random', algo='CPal'):
    """
    Plot runtime vs number of cores for strong scaling (fixed problem size).
    """

    # Load the data
    data = pd.read_csv(f'results_strong_scaling_large.txt', sep=' ')
    df = pd.DataFrame(data)

    df['graph'] = df['graph'].str.replace(r'_(exp|uni)$', '', regex=True)
    df['algo'] = df['algo'].map(lambda x: algo_config[x]['display_name'])
    df['graph'] = df['graph'].map(lambda x: graph_names[x]['display_name'])

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
    grouped_df = (
    df_k.groupby(['graph', 'cores'])
        .agg(
            geomean_runtime=('running_time', lambda x: gmean(x)),
            std_runtime=('running_time', lambda x: np.std(np.log(x)))

        )
        .reset_index()
    )
    # Plotting
    plt.figure(figsize=(13,12))
    for graph, group in grouped_df.groupby('graph'):
        group = group.sort_values('cores')
        plt.plot(group['cores'], group['geomean_runtime'], label=graph,
                 marker=marker.get(graph, 'o'), color=palette.get(graph, None),
                 linewidth=3, markersize=16)
        plt.fill_between(group['cores'],
                         group['geomean_runtime'] - group['std_runtime'],
                         group['geomean_runtime'] + group['std_runtime'],
                         alpha=0.2, color=palette.get(graph, None))
        


    plt.xlabel("\# Controllers")
    plt.ylabel("Runtime (s)")
    plt.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, 1.45),
        ncol=3
        # fontsize=18
    )
    # plt.xscale("log", base=2)  # optional: log-scale x-axis
    ax = plt.gca()
    ax.set_yscale("log")      # optional: log-scale y-axis

    desired_ticks = [0.25,0.5,1,2,4,8,16,32,64,128]
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
    plt.tight_layout(rect=[0, 0, 1, 0.9])

    plt.savefig(f"ASC_BIG_strong.pdf", dpi=500,bbox_inches='tight')
    plt.close()

def print_runtime_ratios(part='random', zen_version=3, algo='kCS'):
    data = pd.read_csv(f'results_strong_scaling_large.txt', sep=' ')
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
plot_strong_scaling_runtime(4, part='random',algo='kCS')

# plot_strong_scaling_runtime_side_by_side(8, part='random', algo='kCS')

print_runtime_ratios()