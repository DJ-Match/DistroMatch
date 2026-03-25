import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gmean
import yaml
from matplotlib.ticker import FuncFormatter

# matplotlib.use('tkagg')

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
palette = {k: v["color"] for k, v in partition_config.items()}
markers = {k: v["marker"] for k, v in partition_config.items()}
dashes = {k: v.get("linestyle", "") for k, v in partition_config.items()}


# Load the data
data = pd.read_csv('results_run_partitioning.csv', sep=' ')
df = pd.DataFrame(data)

df = df[df['p'] != 'metis']

# Data partitioning times
data_pt = pd.read_csv('graph_partition_times.csv')
df_pt = pd.DataFrame(data_pt)
df_pt['none'] =0

df['algo'] = df['algo'].map(lambda x: algo_config[x]['display_name'])
# df['p'] = df['p'].map(lambda x: partition_config[x]['display_name'])


# p = 

parts = ['random', 'diststream', 'umetis', 'uwmetis', 'ukahip',]
# df['p'].unique()

# Duplicate rows for one core
for partitioning in parts:
    if partitioning != 'random':
        copy = df[df['cores'] == 1].copy()
        copy['p'] = partitioning
        df = pd.concat([df, copy], ignore_index=True)

# Set a consistent color palette for partitioning techniques
# colors = sns.color_palette("tab10", len(p))


grouped = df[(df['p'] == 'random') &(df['cores'] == 16) & (df['k'] == 4)].groupby(['algo', 'k', 'cores'])
result = grouped['rounds'].agg(gmean).reset_index()
result2 = data_pt[data_pt['cores']==16][["metis","umetis","uwmetis","ukahip","random","diststream", "cores"]].agg(gmean).reset_index()
print(result.to_string())
print(result2.to_string())

def plot_parameter_mean_std(k, param='running_time', use_baseline=False, add_part_time=False):
    # df_k = df[df['k'] == k]



    if add_part_time:
        # Merge partitioning data with the main dataframe based on graph, cores, and seed
        df_k = df_k.merge(df_pt, on=['graph', 'cores', 'seed'], how='left')

        df_k.fillna(0,inplace=True)

        # Add the partitioning time to the running time based on the corresponding partitioning algorithm
        df_k['running_time'] += df_k.apply(lambda row: row[row['p']], axis=1)

    # Group by 'graph', 'algo', 'p', 'cores' and calculate mean and std for the selected parameter
    grouped_df = df_k.groupby(['graph', 'algo', 'p', 'cores']).agg(
        mean=(param, 'mean'),
        std=(param, 'std')
    ).reset_index()

    print(grouped_df)
    
    # Apply baseline normalization if use_baseline is True
    if use_baseline:
        baseline = grouped_df[grouped_df['p'] == 'random'].set_index(['graph', 'algo', 'cores'])['mean']
        grouped_df['mean'] = grouped_df.apply(
            lambda row: row['mean'] / baseline.loc[(row['graph'], row['algo'], row['cores'])], axis=1
        )
        grouped_df['std'] = grouped_df.apply(
            lambda row: row['std'] / baseline.loc[(row['graph'], row['algo'], row['cores'])], axis=1
        )

    # Get unique graphs and algorithms
    graphs = grouped_df['graph'].unique()
    algos = grouped_df['algo'].unique()

    # Set up the subplots grid |algo| x |graph| + additional partitioning subplot
    fig, axes = plt.subplots(len(algos) + 1, len(graphs), figsize=(15, 12), sharex=False, sharey=False)
    axes = axes.flatten()  # Flatten axes array for easier indexing

    # Iterate over the combinations of graph and algo to create subplots
    for i, (algo, graph) in enumerate([(a, g) for a in algos for g in graphs]):
        ax = axes[i]
        data = grouped_df[(grouped_df['graph'] == graph) & (grouped_df['algo'] == algo)]

        # if(algo == 'CCM'):
            # print(data)

        for p_value in data['p'].unique():
            subset = data[data['p'] == p_value]

            # Plot the parameter mean and std deviation
            ax.plot(subset['cores'], subset['mean'], label=f'{p_value}', color=palette.get(p_value))
            ax.fill_between(subset['cores'], subset['mean'] - subset['std'], subset['mean'] + subset['std'],
                            color=palette.get(p_value), alpha=0.3)

        ax.set_xlabel('\# Controllers')
        # if add_part_time:
        #     ax.set_yscale('log')
            
    handles, labels = ax.get_legend_handles_labels()

    # Plot partitioning times in the last row
    for j, graph in enumerate(graphs):
        ax = axes[len(algos) * len(graphs) + j]  # Last row for partitioning times
        partition_data = df_pt[df_pt['graph'] == graph]


        for partitioning in ['metis', 'umetis', 'uwmetis', 'ukahip', 'random']:
            partition_group = partition_data.groupby('cores')[partitioning].agg(['mean', 'std']).reset_index()

            ax.plot(partition_group['cores'], partition_group['mean'], label=f'{partitioning}', color=   palette.get(partitioning))
            ax.fill_between(partition_group['cores'],
                            partition_group['mean'] - partition_group['std'],
                            partition_group['mean'] + partition_group['std'],
                            color=palette.get(partitioning), alpha=0.3)

        ax.set_xlabel('\# Controllers')
        ax.set_ylabel('Partitioning Time')
        # ax.set_yscale('log')

    # Ensure axes are 2D for consistent plotting
    axes = np.array(axes).reshape(len(algos) + 1, len(graphs))

    # Add titles for each column (graph) and row (algo)
    for ax, graph in zip(axes[0], graphs):
        ax.set_title(f"{graph}", fontsize=16)
    for ax, algo in zip(axes[:, 0], list(algos) + ["Partitioning"]):
        ax.set_ylabel(f"Algorithm: {"\nPartitioning + " if add_part_time and algo!= 'Partitioning'  else ""} {algo}\n\n {"Rel " if add_part_time and algo!= 'Partitioning'  else ""}{param.capitalize()}", rotation=90, size='large', labelpad=15)

    # Add legend and title
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.1), ncol=3)
    fig.suptitle(f"Partitioning Comparison (k={k}, {param})", fontsize=30)

    # Adjust layout and save plot
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig(f"pc_{k}_{param}" + ("_diff" if use_baseline else "") + ".pdf", bbox_inches='tight', dpi=500)

def plot_overall_geomean_comparison_by_algo(param='running_time', use_baseline=False, plot_std=False):
    df_all = df[df[param] > 0].copy()
    df_all['add_part_time'] = False
    
    df_with = df[df[param] > 0].copy()
    df_with = df_with.merge(df_pt, on=['graph', 'cores', 'seed'], how='left')
    df_with.fillna(0, inplace=True)
    df_with['running_time'] += df_with.apply(lambda row: row[row['p']], axis=1)
    df_with['add_part_time'] = True

    combined_df = pd.concat([df_all, df_with], ignore_index=True)
    combined_df = combined_df[combined_df[param] > 0]

    grouped = combined_df.groupby(['add_part_time', 'algo', 'p', 'cores'])

    def geo_stats(x):
        logs = np.log(x)
        return pd.Series({
            'geo_mean': np.exp(logs.mean()),
            'geo_std': np.exp(logs.std())
        })

    stats_df = grouped[param].apply(geo_stats).unstack().reset_index()

    if use_baseline:
        baseline = stats_df[stats_df['p'] == 'random'].set_index(['add_part_time', 'algo', 'cores'])['geo_mean']
        stats_df['geo_mean'] = stats_df.apply(
            lambda row: row['geo_mean'] / baseline.loc[(row['add_part_time'], row['algo'], row['cores'])], axis=1
        )

    # print(stats_df.to_string())

    algos = stats_df['algo'].unique()
    nrows = len(algos)
    fig, axes = plt.subplots(nrows, 2, figsize=(14, 4 * nrows), sharex=True)
    axes = np.array(axes)

    # print(stats_df['algo'] .unique())
    for i, algo in [(0,'MRepLM'), (1, 'MRepS'), (2, 'kMM'), (3, 'kCS')]:
        for j, add_part_time in enumerate([False, True]):
            ax = axes[i, j]
            for method in parts:
                subset = stats_df[(stats_df['algo'] == algo) & 
                                  (stats_df['p'] == method) & 
                                  (stats_df['add_part_time'] == add_part_time)]
                # print(method, subset)
                if not subset.empty:
                    style = (0, dashes.get(method, '-'))

                    label = method if i == 0 and j == 0 else None  # only first plot
                    ax.plot(subset['cores'], subset['geo_mean'], label=partition_config[method]['display_name'], color=palette.get(method), linestyle=style, linewidth=4)
                    if plot_std:
                        lower = subset['geo_mean'] / subset['geo_std']
                        upper = subset['geo_mean'] * subset['geo_std']
                        ax.fill_between(subset['cores'], lower, upper, color=palette.get(method), alpha=0.2)
            if add_part_time:
                ax.set_yscale('log')

            if j == 0:
                ax.set_ylabel(f"{'Rel. ' if use_baseline else ''} Time {'' if use_baseline else '(s)'}")
                ax.set_title(f"{algo}", loc='left')
            if i == nrows - 1:
                ax.set_xlabel("\# Controllers")

            ax.grid(True)

    # Add a single column-level title above the columns
    fig.text(0.26, 0.95, "w/o Partition Time", ha='center')
    fig.text(0.75, 0.95, "w/ Partition Time", ha='center')

    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 0.14), ncol=3)
    # fig.suptitle(f"GeoMean Comparison with and without Partition Time ({param})", fontsize=20)

    plt.tight_layout(rect=[0, 0.125, 1, 0.95])
    plt.savefig(f"geomean_comparison_{param}" + ("_rel" if use_baseline else "") + ".pdf", dpi=500)



# plot_parameter_mean_std(4, 'running_time', False)
# plot_parameter_mean_std(4, 'running_time', True)
# plot_parameter_mean_std(32, 'running_time', False)
# plot_parameter_mean_std(16, 'running_time', False)
# plot_parameter_mean_std(8, 'running_time', False)
# plot_parameter_mean_std(4, 'running_time', False)
# plot_parameter_mean_std(16, 'rounds', False)
# plot_parameter_mean_std(16, 'running_time', True)

# plot_parameter_mean_std(4, 'matching_weight', False)
# plot_parameter_mean_std(4, 'matching_weight', True)
# plot_parameter_mean_std(16, 'matching_weight', False)
# plot_parameter_mean_std(16, 'matching_weight', True)


# plot_overall_geomean_comparison_by_algo('running_time', True)
plot_overall_geomean_comparison_by_algo('running_time', False)

