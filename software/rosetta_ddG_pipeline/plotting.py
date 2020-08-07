"""plotting.py plots heatmaps and other things.

Author: Johanna K.S. Tiemann
Contribution: Magnus Haraldson Høie

Date of last major changes: 2020-05-04

"""
# Standard library imports
import logging as logger
import os


# Third party imports
import matplotlib
matplotlib.use('Agg')
from matplotlib.collections import PathCollection
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns


# Local application imports
from helper import AttrDict, drop_numerical_outliers
from prism_rosetta_parser import read_from_prism, write_prism
import rosetta_paths


# Definitions
aa_order_alphabetical = pd.Series(["A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M",
                                   "N", "P", "Q", "R", "S", "T", "V", "W", "Y"])
aa_order_alphabetical_index_series = pd.Series(
    range(len(aa_order_alphabetical)), index=aa_order_alphabetical)
aa_order_group = pd.Series(["R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "G", "P",
                            "A", "V", "I", "L", "M", "F", "Y", "W"])
convert_aa_order_alph_to_group = {}
for residue, pos_alph in zip(aa_order_alphabetical, range(len(aa_order_alphabetical))):
    pos_group = np.where(aa_order_group == residue)[0]
    convert_aa_order_alph_to_group[pos_alph] = pos_group

# Generate sequence map function


def generate_seqmap(df_variant, df_score, verbose=0):
    # Plot sequence map of scores across a protein sequence
    # Usage: plot_sequence_heatmap(prism_dataframe, variant_col = "variant", score_col = "score")
    # By Magnus Haraldson Høie
    # maghoi@bio.ku.dk
    # May 2020
    arr = np.array(df_variant).reshape(-1, 1)
    index_positions = np.apply_along_axis(
        lambda x: int(str(x[0])[1:-1]), 1, arr)
    mut_positions = np.apply_along_axis(lambda x: str(x[0][-1:]), 1, arr)
    score_positions = df_score

# Extract start end
    start = int(np.sort(index_positions)[0])
    end = int(np.sort(index_positions)[-1])
    protein_sequence_positions = pd.Series(list(range(start, end + 1)))

# Plot
    seq_map = np.zeros(shape=(protein_sequence_positions.iloc[-1], 20))
    seq_map[:] = np.nan

    for i, pos, mut, score in zip(range(len(arr)), index_positions, mut_positions, score_positions):
        if mut in '*=~X':
            continue
        seq_map[pos - 1, aa_order_alphabetical_index_series[mut]] = score

# Remove all empty positions before the first MAVE position
    seq_map = seq_map[start - 1: end + 1]

# Transpose
    seq_map = seq_map.transpose()

    position_means = np.nanmedian(seq_map, axis=0).reshape(1, -1)
    residue_means = np.nanmedian(seq_map, axis=1).reshape(-1, 1)

    if verbose >= 1:
        logger.info('Returning seq_map, position_means, residue_means')
        logger.info(seq_map.shape, protein_sequence_positions.shape,
                    residue_means.shape)
    return(seq_map, protein_sequence_positions, position_means, residue_means)


# Plot sequence map function


def plot_sequence_heatmap(df, out_folder, variant_col='variant', score_col='predicted_ddG', name='protein', title='Title',
                          figsize=(22, 6), font_scale=1.7, reverse_colors=False, return_values=False,
                          verbose=0):
    # Prepare data
    seq_map, positions, position_means, residue_means = generate_seqmap(
        df[variant_col], df[score_col], verbose=verbose)
    position_means = pd.DataFrame(
        position_means.flatten(), index=positions).transpose()

# Stats string
    values = seq_map.flatten()
    values = values[~np.isnan(values)]
    stats_df = np.round(pd.DataFrame(values, columns=['Stats']).describe(), 3)
    stats_string = str(stats_df)

# Re-order alphabetical to chemical groups
    seq_map_aa_group = seq_map.copy()
    for key in convert_aa_order_alph_to_group.keys():
        value = convert_aa_order_alph_to_group[key]
        seq_map_aa_group[value] = seq_map[key]

    residue_means_aa_group_order = residue_means.copy()
    for key in convert_aa_order_alph_to_group.keys():
        value = convert_aa_order_alph_to_group[key]
        residue_means_aa_group_order[value] = residue_means[key]

# Set colors
    color_map = sns.color_palette('coolwarm', 400)
    if reverse_colors:
        color_map = sns.color_palette('coolwarm_r', 400)

# Plot heatmap and position, variant means on same figure
    sns.set_context(context='paper', font_scale=font_scale)
    fig = plt.figure(constrained_layout=True, figsize=figsize)
    widths, heights = [48, 2], [6, 1]
    spec = fig.add_gridspec(ncols=2, nrows=2, width_ratios=widths,
                            height_ratios=heights)

    for i_row, row in enumerate(range(2)):
        for i_col, col in enumerate(range(2)):
            ax = fig.add_subplot(spec[row, col])

    if i_row == 0 and i_col == 0:
        sns.heatmap(seq_map_aa_group, xticklabels=False,
                    yticklabels=aa_order_group, cmap=color_map)
        ax.set_title(f'Sequence map for {name} \n' +
                     f'(Mean value {str(stats_df.loc["mean"][0])} )')
        ax.set_ylabel('Mutation')

    if i_row == 1 and i_col == 0:
        sns.heatmap(position_means, yticklabels=False,
                    cmap=color_map, cbar=False)
        ax.set_xlabel('Sequence position')
        ax.set_ylabel('Median')

    if i_row == 0 and i_col == 1:
        sns.heatmap(residue_means_aa_group_order, yticklabels=aa_order_group, cmap=color_map, cbar=False,
                    annot=True, annot_kws={"size": 10, "weight": "bold"}, fmt='.2g')
        ax.set_ylabel("Median")

# Save
    plt.savefig(os.path.join(out_folder, 'sequence_heatmap.png'),
                facecolor='w', dpi=300, bbox_inches='tight')

    if return_values == True:
        logger.info(
            "Returning sequence map, position means and residue means ordered in aa alphabetical order")
        return(seq_map, position_means, residue_means)


def simple_plot_heatmap(pd_matrix, out_folder, sys_name='', ranges=[], title=""):
    sns.set()
    sns.set_context(context="paper", font_scale=1.5)
    a4_dims = (20, 15)
    fig = plt.figure(constrained_layout=True, figsize=a4_dims)
    if ranges != []:
        ax = sns.heatmap(pd_matrix, mask=pd_matrix.isnull(),
                         cmap="afmhot", vmin=ranges[0], vmax=ranges[1])
    else:
        ax = sns.heatmap(pd_matrix, mask=pd_matrix.isnull(), cmap="afmhot")
    ax.set_title(title, fontsize=60)
    plt.savefig(os.path.join(out_folder, f'{sys_name}_simple_heatmap.png'),
                facecolor='w', dpi=300, bbox_inches='tight')


def joint_plot(pd_dataframe, out_folder, variant_col='variant', exp_score_col='experimental_ddG', ros_score_col='predicted_ddG', sys_name=''):
    g = sns.jointplot(x=exp_score_col, y=ros_score_col,
                      data=pd_dataframe,
                      hue=variant_col,  # "REF",
                      col=variant_col,  # "REF",
                      col_wrap=3,
                      kind="kde", height=7, space=0)

    x0, x1 = g.ax_joint.get_xlim()
    y0, y1 = g.ax_joint.get_ylim()
    lims = [max(x0, y0), min(x1, y1)]
    g.ax_joint.plot(lims, lims, ':k')
    plt.savefig(os.path.join(out_folder, f'{sys_name}_joint.png'),
                facecolor='w', dpi=300, bbox_inches='tight')


def violin_plot(pd_dataframe, out_folder, variant_col='variant', diff_score_col='diff', sys_name=''):
    plt.figure(figsize=(15, 4))
    ax = sns.violinplot(x=variant_col, y=diff_score_col,
                        data=pd_dataframe,
                        color=".8"
                        )
    for artist in ax.lines:
        artist.set_zorder(10)
    for artist in ax.findobj(PathCollection):
        artist.set_zorder(11)

    sns.stripplot(x=variant_col, y=diff_score_col, data=pd_dataframe,
                  jitter=True, ax=ax)

    #ax = (ax.set(ylim=(-5,5)))
    plt.title("Deviation of ddG prediction from experiment")
    plt.savefig(os.path.join(out_folder, f'{sys_name}_violin.png'),
                facecolor='w', dpi=300, bbox_inches='tight')


def volcano_plot(pd_dataframe, out_folder, variant_col='variant', exp_score_col='experimental_ddG', ros_score_col='predicted_ddG', sys_name=''):
    def const_line(*args, **kwargs):
        mini = min(pd_dataframe.min(axis=0)[
                   exp_score_col], pd_dataframe.min(axis=0)[ros_score_col])
        maxi = max(pd_dataframe.max(axis=0)[
                   ros_score_col], pd_dataframe.max(axis=0)[exp_score_col])
        plt.plot([mini, maxi], [mini, maxi], color='grey', linestyle='--')

    plt.figure(figsize=(10, 10))
    markers = ['o', 'v', '^', '<', '>',
               '8', 's', 'p', '*', 'h',
               'H', 'D', 'd', 'P', 'X']
    g = sns.scatterplot(
        x=ros_score_col,  # "predicted_ddG",
        y=exp_score_col,  # "experimental_ddG",
        hue=variant_col,
        data=pd_dataframe,
        markers=markers,
        # palette=['green','orange','brown','dodgerblue','red'],
        legend='full')

    # Plot each individual point separately
    # for i,row in enumerate(pandas_data_no_outlier.predicted_ddG):
    #g.set_prop_cycle( marker=markers)

    g.legend(loc='center right', bbox_to_anchor=(1.16, 0.5), ncol=1)
    const_line()
    plt.savefig(os.path.join(out_folder, f'{sys_name}_volcano.png'),
                facecolor='w', dpi=300, bbox_inches='tight')


def plot_all(folder, sys_name='', drop_outliers=True):

    rosetta_prims_file = os.path.join(folder.output, f'prims_rosetta_XXX_{sys_name}.txt')
    rosetta_metadata, rosetta_dataframe = read_from_prism(rosetta_prims_file)

    experimental_prims_file = os.path.join(
        folder.input, 'prism_mave_input.txt')
    experimental_metadata, experimental_dataframe = read_from_prism(
        experimental_prims_file)

    merged_dataframe = pd.concat([rosetta_dataframe[['variant', 'norm_ddG']],
                                  experimental_dataframe[['score']]], axis=1)
    merged_dataframe = merged_dataframe.rename(
        columns={'norm_ddG': 'predicted_ddG', 'score': 'experimental_ddG'})

    merged_dataframe['diff'] = merged_dataframe['experimental_ddG'].sub(
        merged_dataframe['predicted_ddG'], axis=0)

    merged_metadata = rosetta_metadata.copy()
    merged_metadata['mave'] = experimental_metadata['mave']
    merged_metadata['columns'] = {
        'experimental_ddG': 'Experimental ddG',
        'diff': 'Difference of experimental vs predicted ddG',
        'predicted_ddG': merged_metadata[
            'columns'].pop('norm_ddG')}

    merged_prism_file = os.path.join(folder.analysis, f'prims_rosetta_mave_XXX_{sys_name}_merged.txt')
    comment = ['Automatically generated file from stability_pipeline']
    write_prism(merged_metadata, merged_dataframe,
                merged_prism_file, comment='')

    if drop_outliers:
        drop_numerical_outliers(
            merged_dataframe, variant_col="variant", score_col="predicted_ddG", z_thresh=3)

    plot_sequence_heatmap(merged_dataframe, folder.analysis,
                          variant_col="variant", score_col="predicted_ddG")

    ranges = [min(merged_dataframe[['predicted_ddG']].min(axis=1)), max(
        merged_dataframe[['predicted_ddG']].mean(axis=1)) * 1.75]
    simple_plot_heatmap(merged_dataframe[['predicted_ddG']], folder.analysis,
                        sys_name=sys_name, ranges=ranges, title='Normalized Rosetta run')

    joint_plot(merged_dataframe, folder.analysis, variant_col='variant',
               exp_score_col='experimental_ddG', ros_score_col='predicted_ddG', sys_name=sys_name)

    violin_plot(merged_dataframe, folder.analysis,
                variant_col='variant', diff_score_col='diff', sys_name=sys_name)

    volcano_plot(merged_dataframe, folder.analysis, variant_col='variant',
                 exp_score_col='experimental_ddG', ros_score_col='predicted_ddG', sys_name=sys_name)


if __name__ == '__main__':
    folder = AttrDict()
    folder.update({'input': sys.argv[1], 'output': sys.argv[
                  2], 'analysis': sys.argv[3]})

    if len(sys.argv) > 4:
        sys_name = sys.argv[4]
    else:
        sys_name = ''

    plot_all(folder, sys_name=sys_name)
