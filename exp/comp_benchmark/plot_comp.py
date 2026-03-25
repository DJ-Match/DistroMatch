import argparse
import yaml
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import gmean
from itertools import product
import itertools
from matplotlib.ticker import LogLocator, NullFormatter
import numpy as np

# Global matplotlib settings for better visuals
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

algo_config = config.get("algo_config", {})
param_config = config.get("param_config", {})

# Toggle flags
use_large = False
use_only_rmat = False
use_only_fb = False

# Extract styling from config
palette = {v["display_name"]: v["color"] for v in algo_config.values()}
markers = {v["display_name"]: v["marker"] for v in algo_config.values()}
dashes = {v["display_name"]: v.get("linestyle", "") for v in algo_config.values()}


# --- Plot 1: Parameter over k ---
def plot_param_over_k(data, param, use_baseline='', plot_std=False):
    plt.figure(figsize=(9.6, 6.9))

    grouped = data.groupby(['algo', 'k'])[param].agg(lambda x: gmean(x[x > 0])).reset_index()

    if use_baseline != '':
        pivot = grouped.pivot(index='k', columns='algo', values=param)
        relative = pivot.div(pivot[use_baseline], axis=0)
        grouped = relative.reset_index().melt(id_vars='k', var_name='algo', value_name=param)

    sns.lineplot(
        data=grouped,
        x='k', y=param,
        hue='algo', style='algo',
        palette=palette, markers=markers,
        dashes=dashes, markersize=16
    )

    plt.xlabel("k")
    plt.ylabel(("Relative " if use_baseline else "") + param_config[param])
    if param == "running_time" and not use_large:
        plt.ylim(bottom=0, top=3)
    plt.grid(True, linestyle='--', alpha=0.5)

    leg_pos = 'lower right' if param == 'matching_weight' else 'upper ' + ('right' if use_only_rmat else 'center')
    if use_large and param == 'running_time':
        leg_pos = 'center right'

    plt.legend(title='Algorithm', ncol=2, loc=leg_pos)
    plt.tight_layout()

    filename = f"{param}_per_k{'_rmat' if use_only_rmat else ('_fb' if use_only_fb else ('_large' if use_large else ''))}.pdf"
    plt.savefig(filename)
    plt.close()


    data_used = ('large ' if use_large else '') + ('rmat' if use_only_rmat else ('fb' if use_only_fb else 'all'))

    table_k2 = create_comparison_table(data, param=param, k_val=2)
    print(f"Comparison Table {param} {data_used} (row/col) for k=2:")
    print(table_k2.round(3))

    table_k8 = create_comparison_table(data, param=param, k_val=8)
    print(f"Comparison Table {param} {data_used} (row/col) for k=8:")
    print(table_k8.round(3))

    table_k32 = create_comparison_table(data, param=param, k_val=32)
    print(f"Comparison Table {param} {data_used} (row/col) for k=32:")
    print(table_k32.round(3))

    table_k64 = create_comparison_table(data, param=param, k_val=64)
    print(f"\nComparison Table {param} {data_used} (row/col) for k=64:")
    print(table_k64.round(3))



def plot_params_side_by_side(data, use_baseline='', plot_std=False):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharex=True)

    params = ['running_time', 'matching_weight']
    params_short = ['Time (s)', 'Weight']
    legend_handles, legend_labels = None, None

    for i,(ax, param) in enumerate(zip(axes, params)):
        grouped = data.groupby(['algo', 'k'])[param].agg(lambda x: gmean(x[x > 0])).reset_index()

        if use_baseline != '' and param == 'matching_weight':
            pivot = grouped.pivot(index='k', columns='algo', values=param)
            # print(pivot[use_baseline])
            relative = pivot.div(pivot[use_baseline], axis=0)
            grouped = relative.reset_index().melt(id_vars='k', var_name='algo', value_name=param)
        print(grouped)

        # Reorder the DataFrame: move rows with 'algo_to_front' to the end
        grouped = grouped[grouped['algo'] != 'kCS_1']._append(
            grouped[grouped['algo'] == 'kCS_1']
        )
        grouped = grouped[grouped['algo'] != 'kCS_32']._append(
            grouped[grouped['algo'] == 'kCS_32']
        )

        show_legend = (i == 0)
        sns.lineplot(
            data=grouped,
            x='k', y=param,
            hue='algo', style='algo',
            palette=palette, markers=markers,
            dashes=dashes, markersize=16,
            linewidth=3,
            ax=ax,
            legend=show_legend  # no legend in individual plots
        )

        if show_legend:
            legend_handles, legend_labels = ax.get_legend_handles_labels()
            ax.get_legend().remove()  # remove local legend after capturing


        ax.set_xlabel("k")
        desired_ticks = [2,8,32,64]
        ax.set_xticks(desired_ticks)
        ax.set_xticklabels([r"{\fontseries{m}\selectfont " + str(t) + r"}" for t in desired_ticks])  # Tick labels not bold


    
        if param != 'matching_weight':
            ax.set_yscale('log')
        # if param == 'matching_weight':
        #     ax.grid(True, which="both")
        #     ax.invert_yaxis()
        #     # Add minor ticks to Y-axis
        #     ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10)*0.1, numticks=100))
        #     ax.yaxis.set_minor_formatter(NullFormatter())
        # else:
        
        ax.grid(True, linestyle='--', alpha=0.5)

    

        ax.set_ylabel(("Rel. " if use_baseline and param == 'matching_weight' else "") + params_short[i])

        # if param == "running_time" and not use_large:
        #     ax.set_ylim(bottom=0, top=3)
        #     if use_only_rmat:
        #         ax.set_ylim(bottom=0, top=1.5)


    # Create one shared legend just below the plots
    if use_only_rmat:
        handles, labels = plt.gca().get_legend_handles_labels()
        order = [0,1,4,2,3,7,5,6,8,9,10]
        plt.legend()

        fig.legend([legend_handles[idx] for idx in order],[legend_labels[idx] for idx in order],handlelength=3.5,
                loc='upper center', ncol=4, fontsize=24, bbox_to_anchor=(0.53, 0.08))

    plt.tight_layout(rect=[0, 0, 1, 1])  # leave space at bottom for legend


    filename = f"combined_params_per_k{'_rmat' if use_only_rmat else ('_fb' if use_only_fb else ('_large' if use_large else ''))}.pdf"
    plt.savefig(filename, bbox_inches='tight')
    plt.close()

# --- Plot 2: Parameter per graph at fixed k ---
def plot_param_per_graph(data, param, k=32, use_baseline=''):
    plt.figure(figsize=(10, 8))

    data_k = data[data['k'] == k]
    agg_df = data_k.groupby(['graph', 'algo'])[param].agg(lambda x: gmean(x[x > 0])).reset_index()

    if use_baseline != '':
        pivot = agg_df.pivot(index='graph', columns='algo', values=param)
        rel = pivot.apply(lambda row: row / row[use_baseline], axis=1)
        agg_df = rel.reset_index().melt(id_vars='graph', var_name='algo', value_name=param)

    sns.barplot(
        data=agg_df,
        x='graph', y=param,
        hue='algo',
        palette=palette
    )

    plt.xlabel("Graph")
    plt.ylabel(("Relative " if use_baseline else "") + param_config[param])
    plt.yscale('log')
    plt.xticks(rotation=45, ha='right')
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.legend(title='Algorithm', ncol=2, loc='best')
    plt.tight_layout()

    filename = f"{param}_for_k{k}_per_graph{'_rmat' if use_only_rmat else ('_fb' if use_only_fb else ('_large' if use_large else ''))}.pdf"
    plt.savefig(filename)
    plt.close()


from matplotlib.lines import Line2D
from matplotlib.patches import Patch, FancyBboxPatch

def plot_runtime_with_weight_stars(data, use_baseline='stk', real_time=True):
    plt.rcParams.update({
        "font.size": 24,
    })

    fig, ax = plt.subplots(figsize=(10, 5))

    # --- Aggregate ---
    grouped_time = data.groupby(['algo', 'k'])['running_time'].agg(
        lambda x: gmean(x[x > 0])
    ).reset_index()

    grouped_weight = data.groupby(['algo', 'k'])['matching_weight'].agg(
        lambda x: gmean(x[x > 0])
    ).reset_index(name='matching_weight')

    # --- Relative computations ---
    if use_baseline:
        pivot_time = grouped_time.pivot(index='k', columns='algo', values='running_time')
        pivot_weight = grouped_weight.pivot(index='k', columns='algo', values='matching_weight')
        print(pivot_time)

        rel_time = pivot_time.div(pivot_time[use_baseline], axis=0) \
                             .reset_index().melt(id_vars='k', var_name='algo', value_name='running_time')

        rel_weight = pivot_weight.div(pivot_weight[use_baseline], axis=0) \
                                 .reset_index().melt(id_vars='k', var_name='algo', value_name='matching_weight')
    else:
        rel_time = grouped_time
        rel_weight = grouped_weight

    plot_df = pd.merge(rel_time, rel_weight, on=['algo', 'k'])

    # --- Second axis if real time ---
    if real_time:
        ax2 = ax.twinx()
    else:
        ax2 = ax

    rt_handles = []
    mw_handles = []

    for algo in plot_df['algo'].unique():
        algo_data = plot_df[plot_df['algo'] == algo]
        style = (0, dashes.get(algo, '-'))

        # --- Running Time ---
        time_values = grouped_time[grouped_time['algo'] == algo]['running_time'] if real_time else algo_data['running_time']

        line, = ax.plot(
            algo_data['k'],
            time_values,
            color=palette.get(algo),
            marker=markers.get(algo, 'o'),
            linestyle=style,
            linewidth=3,
            markersize=16,
            zorder=2
        )

        rt_handles.append(Line2D(
            [0], [0],
            color=palette[algo],
            marker=markers.get(algo, 'o'),
            label=algo,
            linestyle=style
        ))

        # --- Matching Weight (stars) ---
        ax2.scatter(
            algo_data['k'],
            algo_data['matching_weight'],
            color=palette.get(algo),
            marker='*',
            linewidths=2,
            s=300,
            edgecolors='black',
            zorder=3
        )

        mw_handles.append(Line2D(
            [0], [0],
            color=palette[algo],
            marker='*',
            linestyle='None',
            markeredgecolor='black',
            markersize=16,
            label=algo
        ))

    # --- Labels ---
    ax.set_xlabel("k")

    if real_time:
        ax.set_ylabel("Time (s)")
        ax.set_yscale('log')
        ax2.set_ylabel("Relative Weight")
    else:
        ax.set_ylabel("Relative Value")

    ax.grid(True, linestyle='--', alpha=0.5)

    # --- Custom Legend Box ---
    legend_x = 0.5
    legend_y = 0.6
    legend_width = 0.49
    row_height = 0.06

    # Background box
    box = FancyBboxPatch((legend_x - 0.03, legend_y - 0.01 - 3*row_height),
                         legend_width, 5 * row_height,
                         boxstyle="round,pad=0.02",
                         transform=ax.transAxes,
                         facecolor='white', 
                         edgecolor='gray',
                         alpha=0.95,
                         zorder=5,
                         linewidth=1)
    # box.set_clip_on(False)
    ax.add_patch(box)

    # Running time legend
    rt_legend = ax.legend(rt_handles,  [''] * len(rt_handles),
                          loc='upper left',
                          bbox_to_anchor=(legend_x +0.12 , legend_y + 1.5*row_height),
                          ncol=1, frameon=False, handletextpad=0.0, columnspacing=0.0, markerscale=2)
    ax.add_artist(rt_legend)

    # Matching weight legend
    mw_legend = ax.legend(mw_handles,  [''] * len(mw_handles),
                          loc='upper left',
                          bbox_to_anchor=(legend_x + 0.32, legend_y + 1.5*row_height),
                          ncol=1, frameon=False, handletextpad=0.0, columnspacing=0.0,markerscale=1)
    ax.add_artist(mw_legend)

    for i,h in enumerate(mw_handles):
        label = h.get_label()
        ax.text(legend_x-0.03 , legend_y -0.03- 2*row_height+ (row_height+0.05 if i==0 else 0), label,fontsize=18, transform=ax.transAxes,zorder=6,
        ha='left', va='center')

    # Titles as legend headers
    ax.text(legend_x + 0.15  , legend_y + row_height, "Time", fontsize=18, transform=ax.transAxes,zorder=6,
            ha='left', va='center')
    ax.text(legend_x + 0.34, legend_y + row_height , "Weight", fontsize=18, transform=ax.transAxes,zorder=6,
            ha='left', va='center')

    plt.tight_layout()

    # Save
    filename = f"runtime_weight_line_star_over_k{'_large' if use_large else ''}{'_rmat' if use_only_rmat else '_fb' if use_only_fb else ''}.pdf"
    plt.savefig(filename, bbox_inches='tight')
    plt.close()


    # for k in algo_data['k']:
    #     for param in ['running_time', 'matching_weight']:

    #         table = create_comparison_table(data, param=param, k_val=k)
    #         print(f"Comparison Table {param} large (row/col) for k={k}:")
    #         print(table.round(3))

    # generate_combined_latex_table(data)


def generate_combined_latex_table(data, target_algo='stk', baseline_algo='DCCM_32'):
    output_lines = []
    print(data['algo'].unique())

    for graph in sorted(data['graph'].unique()):
        graph_df = data[data['graph'] == graph]
        k_vals = sorted(graph_df['k'].unique())

        # Filter out k-values where either target or baseline is missing
        valid_k_data = []
        for k_val in k_vals:
            k_df = graph_df[graph_df['k'] == k_val]
            target_df = k_df[k_df['algo'] == target_algo]
            baseline_df = k_df[k_df['algo'] == baseline_algo]
            if target_df.empty or baseline_df.empty:
                continue

            t_runtime = gmean(target_df['running_time'][target_df['running_time'] > 0])
            t_weight = gmean(target_df['matching_weight'][target_df['matching_weight'] > 0])
            b_runtime = gmean(baseline_df['running_time'][baseline_df['running_time'] > 0])
            b_weight = gmean(baseline_df['matching_weight'][baseline_df['matching_weight'] > 0])

            rel_runtime = t_runtime / b_runtime if b_runtime > 0 else float('inf')
            rel_weight = t_weight / b_weight if b_weight > 0 else float('inf')

            valid_k_data.append({
                'k': k_val,
                'Runtime': f"{t_runtime:.2f}",
                'RelRuntime': f"{rel_runtime:.2f}",
                'Weight': f"{t_weight:.2e}",
                'RelWeight': f"{rel_weight:.2f}"
            })

        if not valid_k_data:
            continue

        # Start multirow
        first_row = True
        for row in valid_k_data:
            if first_row:
                output_lines.append(
                    f"\\multirow{{{len(valid_k_data)}}}{{2.5cm}}{{{graph}}} & {row['k']} & {row['Runtime']} & {row['RelRuntime']} & {row['Weight']} & {row['RelWeight']} \\\\"
                )
                first_row = False
            else:
                output_lines.append(
                    f"& {row['k']} & {row['Runtime']} & {row['RelRuntime']} & {row['Weight']} & {row['RelWeight']} \\\\"
                )
        output_lines.append("\\hline")

    if output_lines:
        latex_table = "\n".join(output_lines)
        print("\nFormatted LaTeX Table:\n")
        print(latex_table)
        return latex_table
    else:
        print("No data to generate table.")
        return ""

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

def create_algo_comparison_table(data, algo1, algo2):
    """
    Create a comparison table for two algorithms across all k values for running_time and matching_weight.
    
    Parameters:
        data (pd.DataFrame): The input dataset with columns ['k', 'algo', 'running_time', 'matching_weight'].
        algo1 (str): The name of the first algorithm.
        algo2 (str): The name of the second algorithm.
    
    Returns:
        pd.DataFrame: A table where columns are k values and rows are ratios for the parameters.
    """
    # Get all unique k values
    k_values = sorted(data['k'].unique())
    
    # Initialize a dictionary to store the results
    results = {
        f'running_time: {algo1} / {algo2}': [],
        f'running_time: {algo2} / {algo1}': [],
        f'matching_weight: {algo1} / {algo2}': [],
        f'matching_weight: {algo2} / {algo1}': []
    }
    
    # Iterate over all k values
    for k_val in k_values:
        # Filter data for the current k value
        df_k = data[data['k'] == k_val]
        
        # Compute geometric means for algo1 and algo2 for both parameters
        geo_means = (
            df_k[df_k['algo'].isin([algo1, algo2])]
            .groupby('algo')[['running_time', 'matching_weight']]
            .agg(gmean)
        )
        
        # Ensure both algorithms are present for this k value
        if algo1 in geo_means.index and algo2 in geo_means.index:
            # Compute ratios for running_time
            rt_algo1 = geo_means.loc[algo1, 'running_time']
            rt_algo2 = geo_means.loc[algo2, 'running_time']
            results[f'running_time: {algo1} / {algo2}'].append(rt_algo1 / rt_algo2)
            results[f'running_time: {algo2} / {algo1}'].append(rt_algo2 / rt_algo1)
            
            # Compute ratios for matching_weight
            mw_algo1 = geo_means.loc[algo1, 'matching_weight']
            mw_algo2 = geo_means.loc[algo2, 'matching_weight']
            results[f'matching_weight: {algo1} / {algo2}'].append(mw_algo1 / mw_algo2)
            results[f'matching_weight: {algo2} / {algo1}'].append(mw_algo2 / mw_algo1)
        else:
            # If one of the algorithms is missing, append NaN
            results[f'running_time: {algo1} / {algo2}'].append(float('nan'))
            results[f'running_time: {algo2} / {algo1}'].append(float('nan'))
            results[f'matching_weight: {algo1} / {algo2}'].append(float('nan'))
            results[f'matching_weight: {algo2} / {algo1}'].append(float('nan'))
    
    # Create a DataFrame with k values as columns and the results as rows
    comparison_table = pd.DataFrame(results, columns=results.keys(), index=k_values).T
    comparison_table.columns = [f'k={k}' for k in k_values]
    
    print(comparison_table.round(3).applymap(lambda x: f"${x} \\times$,"))
    return comparison_table

def main():
    global use_large, use_only_rmat, use_only_fb, data

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Run experiments with different configurations.")
    parser.add_argument("--size", choices=["small", "large"], required=True, help="Specify the graph size: small or large.")
    args = parser.parse_args()

    use_large = args.size == "large"

    flag_options = {
        "use_only_rmat": [False, True],
        "use_only_fb": [False, True],
    }

    combinations = list(product(*flag_options.values()))
    if use_large:
        data = pd.read_csv(f"results_run_seq_comp_large.csv")
        data['algo'] = data['algo'].map(lambda x: algo_config[x]['display_name'])

        plot_runtime_with_weight_stars(data, use_baseline='stk')

        # generate_combined_latex_table(data, target_algo='stk', baseline_algo='kCS_128')
        
        # plot_param_over_k(data, param="matching_weight", use_baseline='stk')


    else:
        for (only_rmat, only_fb) in combinations:
            # Skip invalid combinations
            if (use_large and only_fb)  or (use_large and only_rmat or (only_fb and only_rmat)):
                continue

            # Set flags
            use_only_rmat = only_rmat
            use_only_fb = only_fb

            # Reload data
            data = pd.read_csv(f"results_run_seq_comp.csv")
            data = data[data['algo'] != 'NCb']

            graph_prefix = '^rmat' if use_only_rmat else ('^fb' if use_only_fb else '')
            data = data[data["graph"].str.match(graph_prefix)]

            data['algo'] = data['algo'].map(lambda x: algo_config[x]['display_name'])

            # Plot
            plot_params_side_by_side(data, use_baseline='kEC')


            # Print comparisons
            # print("\n-------------------------------------------------\n")
            # print(data["graph"].unique())
            # print(f"Processing config: rmat={only_rmat}, fb={only_fb}")
            # create_algo_comparison_table(data, 'MRepLM_1', 'MRepS_1')
            # create_algo_comparison_table(data, 'MRepLM_32', 'MRepS_32')
            # create_algo_comparison_table(data, 'kMM_1', 'kEC')
            # create_algo_comparison_table(data, 'kMM_32', 'stk')
            # create_algo_comparison_table(data, 'kMM_32', 'kEC')
            # create_algo_comparison_table(data, 'kCS_1', 'kEC')
            # create_algo_comparison_table(data, 'kCS_1', 'stk')
            # create_algo_comparison_table(data, 'kCS_32', 'kEC')
            # create_algo_comparison_table(data, 'kCS_32', 'stk')
            # create_algo_comparison_table(data, 'kCS_32', 'kMM_32')
            # print("\n-------------------------------------------------\n")

if __name__ == "__main__":
    main()

