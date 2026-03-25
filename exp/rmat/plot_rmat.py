import yaml
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import gmean
import numpy as np
import itertools


# Plot settings
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.size": 24,
    "text.latex.preamble": r"""
        \usepackage{amsmath}
        \usepackage{bm}
        \renewcommand{\familydefault}{\sfdefault}  % optional: switch to sans-serif if you want
        \renewcommand{\seriesdefault}{\bfdefault}  % <-- this sets all text to bold by default
    """
})
# === Load Configuration ===
with open('../config.yaml', 'r') as f:
    config = yaml.safe_load(f)

algo_config = config.get("algo_config", {})
param_config = config.get("short_param_config", {})

# === Load and Preprocess Data ===
data = pd.read_csv("results_run_rmat.csv")

# Filter out specific algorithms
data = data[~data['algo'].isin(['CCM_1', 'CCM_32', 'NCb'])]

# Map display names from config
data['algo'] = data['algo'].map(lambda x: algo_config[x]['display_name'])

# Extract plotting styles from config
palette = {v["display_name"]: v["color"] for _, v in algo_config.items()}
markers = {v["display_name"]: v["marker"] for _, v in algo_config.items()}
dashes = {v["display_name"]: v.get("linestyle", "") for _, v in algo_config.items()}

# Helper to compute geometric mean safely
def safe_gmean(x):
    x = np.array(x)
    x = x[x > 0]
    return gmean(x) if len(x) > 0 else np.nan

# === Function 1: Plot Parameter over k ===
def plot_per_algo(data, param, baseline_algo="", out_file="rmat"):
    fig, ax1 = plt.subplots(figsize=(10, 6))

    plot_data = data.copy()

    # Normalize by baseline if specified
    if baseline_algo:
        plot_data = plot_data.pivot(index=['k', 'graph', 'seed'], columns='algo', values=param)
        plot_data = plot_data.apply(lambda row: row / row[algo_config[baseline_algo]["display_name"]], axis=1)
        plot_data = plot_data.reset_index().melt(id_vars=['k', 'graph', 'seed'], var_name='algo', value_name=param)

    sns.boxplot(data=plot_data, x='algo', hue='algo', legend=False, y=param, palette=palette, ax=ax1)

    ax1.set_xlabel("")
    ax1.set_ylabel(("Relative " if baseline_algo else "") + param_config[param])
    if param == "running_time":
        ax1.set_yscale("log")
    ax1.grid(True)

    # === Overlay average matching weight as stars (secondary axis) ===
    if param == "running_time":
        ax2 = ax1.twinx()
        avg_weights = data.groupby('algo')['matching_weight'].agg(safe_gmean)

        algo_positions = {label.get_text(): i for i, label in enumerate(ax1.get_xticklabels())}

        for algo_name, weight in avg_weights.items():
            if algo_name in algo_positions:
                xpos = algo_positions[algo_name]
                ax2.plot(
                    xpos, weight,
                    marker='*',
                    markersize=15,
                    color=palette[algo_name],
                    markeredgecolor='black',
                    markeredgewidth=1.2,
                    linestyle='None'
                )

        ax2.set_ylabel("Average Matching Weight")
        ax2.tick_params(axis='y', labelsize=12)

        ax1.plot(
            [0.005], [-0.18],
            transform=ax1.transAxes,
            marker='*',
            markersize=12,
            color='black',
            linestyle='None'
        )

    plt.tight_layout()
    plt.savefig(f"{out_file}.pdf")
    plt.close()

def plot_per_algo_weight_baseline(data, param, baseline_algo="", out_file="rmat"):
    fig, ax1 = plt.subplots(figsize=(10, 7))

    display_baseline = algo_config[baseline_algo]["display_name"] if baseline_algo else ""

    plot_data = data.copy()

    # Normalize by baseline if specified
    if baseline_algo:
        plot_data = plot_data.pivot(index=['k', 'graph', 'seed'], columns='algo', values=param)
        plot_data = plot_data.apply(lambda row: row / row[display_baseline], axis=1)
        plot_data = plot_data.reset_index().melt(id_vars=['k', 'graph', 'seed'], var_name='algo', value_name=param)

    # Boxplot
    sns.boxplot(data=plot_data, x='algo', hue='algo', legend=False, y=param, palette=palette, ax=ax1)

    ax1.set_xlabel("")
    ax1.set_ylabel(("Rel. " if baseline_algo else "") + param_config[param])
    plt.setp(ax1.get_xticklabels(), rotation=45, ha='right')

    if param == "running_time":
        ax1.set_yscale("log")
    ax1.grid(True)

    # === Overlay average relative matching weight as stars ===
    if param == "running_time" and baseline_algo:
        ax2 = ax1.twinx()

        # Step 1: Filter valid matching weights
        incomplete_algos = data[data['matching_weight'].isna()]['algo'].unique()


        matching_data = data[['graph', 'algo', 'matching_weight']].dropna()
        matching_data = matching_data[matching_data['matching_weight'] > 0]

        # Step 2: Compute geo mean of raw values per algo
        geo_means = matching_data.groupby('algo')['matching_weight'].agg(gmean)

        avg_relative_weights = {}
        baseline_val = geo_means.get(display_baseline, np.nan)

        for algo_name, algo_val in geo_means.items():
                avg_relative_weights[algo_name] = algo_val / baseline_val

        print("Incomplete results:",incomplete_algos)


        # Step 3: For incomplete algos, fallback to row-wise relative values
        all_algos = data['algo'].unique()
        for algo in all_algos:
            if algo in incomplete_algos:
                # Find graphs where both algo and baseline exist
                subset = matching_data[matching_data['algo'].isin([algo, display_baseline])]
                valid_graphs = subset.groupby('graph')['algo'].nunique() == 2
                valid_graphs = valid_graphs[valid_graphs].index

                # print(valid_graphs)

                rels = []
                for g in valid_graphs:
                    vals = subset[subset['graph'] == g].set_index('algo')['matching_weight']

                    # print(vals)
                    if display_baseline in vals and algo in vals:
                        # print(vals[algo]/ vals[display_baseline])
                        rels.append(safe_gmean(vals[algo]) / safe_gmean(vals[display_baseline]))
                # print(rels)
                if rels:
                    avg_relative_weights[algo] = safe_gmean(rels)


        # Plot stars
        algo_positions = {label.get_text(): i for i, label in enumerate(ax1.get_xticklabels())}
        for algo_name, rel_weight in avg_relative_weights.items():
            if algo_name in algo_positions and not pd.isna(rel_weight):
                xpos = algo_positions[algo_name]
                ax2.plot(
                    xpos, rel_weight,
                    marker='*',
                    markersize=15,
                    color=palette.get(algo_name, 'gray'),
                    markeredgecolor='black',
                    markeredgewidth=1.2,
                    linestyle='None'
                )

    ax2.set_ylabel("Rel. Weight")
    ax2.tick_params(axis='y', labelsize=20)

    # Legend marker
    ax1.plot(
        [0.005], [-0.18],
        transform=ax1.transAxes,
        marker='*',
        markersize=12,
        color='black',
        linestyle='None'
    )

    plt.tight_layout()
    plt.savefig(f"{out_file}.pdf",bbox_inches='tight')
    plt.close()


# === Function 3: Plot per Graph for Fixed k ===
def plot_param_per_graph(data, param, k=32, baseline_algo=""):
    fig, ax = plt.subplots(figsize=(12, 8))

    data_k = data[data['k'] == k]

    # Aggregate using geometric mean
    avg_data = data_k.groupby(['graph', 'algo'])[param].agg(safe_gmean).reset_index()

    if baseline_algo:
        avg_data = avg_data.pivot(index='graph', columns='algo', values=param)
        avg_data = avg_data.apply(lambda row: row / row[algo_config[baseline_algo]["display_name"]], axis=1)
        avg_data = avg_data.reset_index().melt(id_vars='graph', var_name='algo', value_name=param)

    sns.barplot(data=avg_data, x='graph', y=param, hue='algo', palette=palette, ax=ax)
    ax.set_xlabel("Graph")
    ax.set_ylabel(("Relative " if baseline_algo else "") + param)
    ax.set_yscale("log")
    ax.set_title(f"{'Relative ' if baseline_algo else ''}{param} for each graph (k={k})")
    plt.xticks(rotation=45, ha='right')
    ax.legend(title='Algorithm')
    plt.tight_layout()

    plt.savefig(f"rmat_{param}_for_k{k}_per_graph.pdf")
    plt.close()

def create_comparison_table(data, param, k_val):
    # Step 1: Filter data to the specified k
    df_k = data[data['k'] == k_val]

    # Step 2: Compute geo mean for each algorithm
    geo_means = (
        df_k[df_k[param] > 0]
        .groupby('algo')[param]
        .agg(gmean)
    )

    # Step 3: Create ratio table
    algos = geo_means.index.tolist()

    ratio_table = pd.DataFrame(index=algos, columns=algos)
    print(k_val)
    print(geo_means)

    for a1, a2 in itertools.product(algos, repeat=2):
        # print(a1, "/", a2, geo_means[a1],'/', geo_means[a2], geo_means[a1] / geo_means[a2])
        ratio_table.loc[a1, a2] = geo_means[a1] / geo_means[a2]

    ratio_table = ratio_table.astype(float)
    return ratio_table

def print_incomplete_graphs(data, param_col='matching_weight'):
    # Count number of algos per graph with valid entries
    non_null_counts = (
        data.dropna(subset=[param_col])
            .groupby('graph')['algo']
            .nunique()
    )

    total_algos = data['algo'].nunique()
    incomplete_graphs = non_null_counts[non_null_counts < total_algos].index.tolist()

    print(incomplete_graphs)


# === Usage ===
param_list = ['running_time', 'matching_weight']


plot_per_algo_weight_baseline(data, param="running_time", baseline_algo="kEC")


# table = create_comparison_table(data, param="running_time", k_val=16)
# print(table.round(3))

# table_w = create_comparison_table(data, param="matching_weight", k_val=16)
# print(table_w.round(3))

# print_incomplete_graphs(data)