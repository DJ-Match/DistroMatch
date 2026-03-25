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





def plot_weak_scaling_runtime(algo='GKM'):
    """
    Plot runtime vs number of cores for weak scaling.
    """

    # Load the data
    data = pd.read_csv(f'results_weak_scaling_large.txt', sep=' ')
    df = pd.DataFrame(data)

    df['algo'] = df['algo'].map(lambda x: algo_config[x]['display_name'])


    # Set a consistent color palette for partitioning techniques
    # colors = ["#CC6677","#332288","#DDCC77","#117733","#88CCEE","#882255","#44AA99","#999933","#AA4499"]

    df_k = df[(df['algo'] == algo)]

    # Extract "type" from graph name
    # e.g. rmat_s21_ef24_0.25_0.25_0.25_0.25_exp -> 0.25_0.25_0.25_0.25
    df_k['graph'] = df_k['graph'].str.replace(r'_(exp|uni)$', '', regex=True)
    print(df_k['graph'].unique())
    df_k['type'] = df_k['graph'].str.extract(r'(\d+\.\d+_\d+\.\d+_\d+\.\d+_\d+\.\d+)')[0]

    df_k['type'] = df_k['type'].map(lambda x: graph_names[x]['display_name'])


    #  Group by type and cores
    grouped_df = (
    df_k.groupby(['type', 'cores'])
        .agg(
            mean_runtime=('running_time', lambda x: gmean(x)),
            std_runtime=('running_time', lambda x: np.std(np.log(x)))

        )
        .reset_index()
    )

    # Plotting
    plt.figure(figsize=(11, 8))
    # plt.figure(figsize=(13,12))
    for t, group in grouped_df.groupby('type'):
        group = group.sort_values('cores')
        plt.plot(group['cores'], group['mean_runtime'], label=t,
                 marker=marker.get(t, 'o'),
                 color=palette.get(t, None),
                 linewidth=3, markersize=12)
        plt.fill_between(group['cores'],
                         group['mean_runtime'] - group['std_runtime'],
                         group['mean_runtime'] + group['std_runtime'],
                         alpha=0.2, color=palette.get(t, None))

        


    plt.xlabel("\# Controllers")
    plt.ylabel("Runtime (s)")
    plt.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, 1.3),  # centered above the axes
        ncol=3
        # fontsize=18
    )
    # plt.xscale("log", base=2)  # optional: log-scale x-axis
    ax = plt.gca()
    ax.set_yscale("log")      # optional: log-scale y-axis

    desired_ticks = [0.25,0.5,1,2,4,8,16,32]
    ax.yaxis.set_major_locator(ticker.FixedLocator(desired_ticks))
    # ax.set_yticklabels([str(t) for t in desired_ticks], fontweight='normal')
    ax.set_ylim(bottom=2, top=32)
    ax.set_yticklabels([r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_ticks], fontsize=24)  # Tick labels not bold


    # Disable minor ticks completely
    ax.yaxis.set_minor_locator(ticker.NullLocator())

    # If x-axis is log scale
    ax.set_xscale("log")

    # Fix x-axis ticks to specific values
    desired_xticks = [1, 2, 4, 8, 16, 32, 64,128]
    ax.xaxis.set_major_locator(ticker.FixedLocator(desired_xticks))
    ax.set_xlim(left=30, right=135)
    # ax.set_xticklabels([str(t) for t in desired_xticks])
    ax.set_xticklabels([r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_xticks], fontsize=24)  # Tick labels not bold


    # Disable minor ticks on x-axis
    ax.xaxis.set_minor_locator(ticker.NullLocator())

    plt.grid()
    plt.tight_layout(rect=[0, 0, 1, 0.9])  # Leave space on the right for the legend

    plt.savefig(f"ASC_BIG_weak.pdf", dpi=500,bbox_inches='tight')
    plt.close()


plot_weak_scaling_runtime(algo='kCS')
