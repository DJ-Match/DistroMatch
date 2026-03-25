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


def plot_scaling_combined(k=4, part='random', algo='CPal'):

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(22,10))

    # =====================================================
    # Strong Scaling
    # =====================================================

    data = pd.read_csv('results_strong_scaling_large.txt', sep=' ')
    df = pd.DataFrame(data)

    df['graph'] = df['graph'].str.replace(r'_(exp|uni)$', '', regex=True)
    df['algo'] = df['algo'].map(lambda x: algo_config[x]['display_name'])
    df['graph'] = df['graph'].map(lambda x: graph_names[x]['display_name'])

    p = ['metis','none','random','umetis','uwmetis','ukahip','diststream']

    for partitioning in p:
        if partitioning != 'random':
            copy = df[df['cores'] == 1].copy()
            copy['p'] = partitioning
            df = pd.concat([df, copy], ignore_index=True)

    df_k = df[(df['k'] == k) & (df['p'] == part) & (df['algo'] == algo)]

    grouped_df = (
        df_k.groupby(['graph', 'cores'])
        .agg(
            geomean_runtime=('running_time', lambda x: gmean(x)),
            std_runtime=('running_time', lambda x: np.std(np.log(x)))
        )
        .reset_index()
    )



    for graph, group in grouped_df.groupby('graph'):
        group = group.sort_values('cores')

        ax1.plot(group['cores'], group['geomean_runtime'],
                 label=graph,
                 marker=marker.get(graph, 'o'),
                 color=palette.get(graph, None),
                 linewidth=3, markersize=12)

        ax1.fill_between(group['cores'],
                         group['geomean_runtime'] - group['std_runtime'],
                         group['geomean_runtime'] + group['std_runtime'],
                         alpha=0.2,
                         color=palette.get(graph, None))

    ax1.set_xlabel("\# Controllers")
    ax1.set_ylabel("Runtime (s)")

    ax1.set_xscale("log")
    ax1.set_yscale("log")

    desired_ticks = [0.25,0.5,1,2,4,8,16,32,64,128]
    ax1.yaxis.set_major_locator(ticker.FixedLocator(desired_ticks))
    ax1.set_yticklabels([r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_ticks], fontsize=24)
    ax1.yaxis.set_minor_locator(ticker.NullLocator())

    desired_xticks = [1, 2, 4, 8, 16, 32, 64,128]
    ax1.xaxis.set_major_locator(ticker.FixedLocator(desired_xticks))
    ax1.set_xticklabels([r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_xticks], fontsize=24)
    ax1.xaxis.set_minor_locator(ticker.NullLocator())

    ax1.grid()

    ax1.set_title("(a) Strong scaling")

    ax1.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.25),
        ncol=3,
        fontsize=24
    )

    # =====================================================
    # Weak Scaling
    # =====================================================

    data = pd.read_csv('results_weak_scaling_large.txt', sep=' ')
    df = pd.DataFrame(data)

    df['algo'] = df['algo'].map(lambda x: algo_config[x]['display_name'])

    df_k = df[(df['algo'] == algo)]

    df_k['graph'] = df_k['graph'].str.replace(r'_(exp|uni)$', '', regex=True)
    df_k['type'] = df_k['graph'].str.extract(r'(\d+\.\d+_\d+\.\d+_\d+\.\d+_\d+\.\d+)')[0]
    df_k['type'] = df_k['type'].map(lambda x: graph_names[x]['display_name'])

    grouped_df = (
        df_k.groupby(['type', 'cores'])
        .agg(
            mean_runtime=('running_time', lambda x: gmean(x)),
            std_runtime=('running_time', lambda x: np.std(np.log(x)))
        )
        .reset_index()
    )

    for t, group in grouped_df.groupby('type'):
        group = group.sort_values('cores')

        ax2.plot(group['cores'], group['mean_runtime'],
                 label=t,
                 marker=marker.get(t, 'o'),
                 color=palette.get(t, None),
                 linewidth=3, markersize=12)

        ax2.fill_between(group['cores'],
                         group['mean_runtime'] - group['std_runtime'],
                         group['mean_runtime'] + group['std_runtime'],
                         alpha=0.2,
                         color=palette.get(t, None))

    ax2.set_xlabel("\# Controllers")
    ax2.set_ylabel("Runtime (s)")

    ax2.set_xscale("log")
    ax2.set_yscale("log")

    ax2.set_ylim(bottom=2, top=32)

    desired_ticks = [0.25,0.5,1,2,4,8,16,32]
    ax2.yaxis.set_major_locator(ticker.FixedLocator(desired_ticks))
    ax2.set_yticklabels([r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_ticks], fontsize=24)  # Tick labels not bold

    desired_xticks = [1, 2, 4, 8, 16, 32, 64,128]
    ax2.xaxis.set_major_locator(ticker.FixedLocator(desired_xticks))
    ax2.set_xlim(left=30, right=135)
    ax2.set_xticklabels([r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_xticks], fontsize=24)  # Tick labels not bold


    # Disable minor ticks on x-axis
    ax2.xaxis.set_minor_locator(ticker.NullLocator())

    # Disable minor ticks completely
    ax2.yaxis.set_minor_locator(ticker.NullLocator())


    ax2.grid()

    ax2.set_title("(b) Weak scaling")

    ax2.legend(
        loc='upper center',
        bbox_to_anchor=(0.5, -0.25),
        ncol=3,
        fontsize=24
    )


    plt.tight_layout(rect=[0,0.08,1,1])

    # plt.tight_layout()

    plt.savefig("ASC_BIG_scaling_combined.pdf",
                dpi=500,
                bbox_inches='tight')

    plt.close()

plot_scaling_combined(algo='kCS')