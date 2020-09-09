"""analysis and statistics of results

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-06-15

"""
# Standard library imports
import logging as logger
import os


# Third party imports
import numpy as np
import pandas as pd
import json
from scipy import stats
import subprocess


# Local application imports
from helper import AttrDict, drop_numerical_outliers
from prism_rosetta_parser import read_from_prism, write_prism
from PrismData import PrismParser
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


def calc_cc(prismfile):
    parser = PrismParser()
    data = parser.read(prismfile)
    data_frame_raw = data.dataframe
    data_frame1 = data_frame_raw.dropna()

    spearman = stats.spearmanr(np.array(data_frame1)[
                               :, 1], np.array(data_frame1)[:, 2])
    pearson = stats.pearsonr(np.array(data_frame1)[
                             :, 1], np.array(data_frame1)[:, 2])
    mannwhitney = stats.mannwhitneyu(
        np.array(data_frame1)[:, 1], np.array(data_frame1)[:, 2])
    linregress = stats.linregress(np.array(data_frame1)[:, 1].astype(
        float), np.array(data_frame1)[:, 2].astype(float))
    best_fit = np.polyfit(np.array(data_frame1)[:, 1].astype(
        float), np.array(data_frame1)[:, 2].astype(float), 1)
    statistics = {
        'spearman': {
            'correlation': spearman[0],
            'p-value': spearman[1]
        },
        'pearson': {
            'correlation': pearson[0],
            'p-value': pearson[1]
        },
        'mannwhitney': {
            'correlation': mannwhitney[0],
            'p-value': mannwhitney[1]
        },
        'linregress': {
            'r-square': linregress[2]**2,
            'r-value': linregress[2],
            'p-value': linregress[3]
        },
        'polyfit': {
            'slope': best_fit[0],
            'intercept': best_fit[1]
        },
    }
    return statistics


def calc_all(folder, sys_name='', drop_outliers=True, drop_pro=True):

    rosetta_prims_file = os.path.join(folder.output, f'prims_rosetta_XXX_{sys_name}.txt')
    rosetta_metadata, rosetta_dataframe = read_from_prism(rosetta_prims_file)

    experimental_prims_file = os.path.join(
        folder.input, 'prism_mave_input.txt')
    experimental_metadata, experimental_dataframe = read_from_prism(
        experimental_prims_file)
    experimental_dataframe = experimental_dataframe[~experimental_dataframe.variant.str.contains(r".\*", case=False, na=False)]
    experimental_dataframe = experimental_dataframe[~experimental_dataframe.variant.str.contains(r".\=", case=False, na=False)]

    merged_prism_file_raw = os.path.join(folder.analysis, f'prims_merged_XXX_{sys_name}_merged_raw.txt')

    shell_command = (f'python {os.path.join(rosetta_paths.prims_parser, "PrismData.py")} '
                     f'{rosetta_prims_file} {experimental_prims_file} '
                     f'--merge {merged_prism_file_raw}')
    subprocess.call(shell_command, shell=True)
    merged_metadata, merged_dataframe = read_from_prism(merged_prism_file_raw)

    merged_dataframe = merged_dataframe.rename(
        columns={merged_dataframe.keys()[1]: 'predicted_ddG', merged_dataframe.keys()[2]: 'std_ddG', merged_dataframe.keys()[3]: 'experimental_ddG'})

    merged_dataframe['diff'] = merged_dataframe['predicted_ddG'].sub(
        merged_dataframe['experimental_ddG'], axis=0)

    merged_dataframe.sort_values('variant', inplace=True)
    merged_dataframe.reset_index(drop=True, inplace=True)
    merged_dataframe = merged_dataframe[['variant', 'predicted_ddG', 'experimental_ddG', 'diff']]
    merged_metadata = merged_metadata.copy()
   # merged_metadata['mave'] = experimental_metadata['mave']
    merged_metadata['columns'] = {
        'predicted_ddG': rosetta_metadata[
            'columns'].pop('norm_ddG'),
        'experimental_ddG': 'Experimental ddG',
        'diff': 'Difference of experimental vs predicted ddG',
    }

    comment = ['Automatically generated file from stability_pipeline']
    merged_prism_file = os.path.join(folder.analysis, f'prims_merged_XXX_{sys_name}_merged.txt')

    write_prism(merged_metadata, merged_dataframe,
                merged_prism_file, comment=comment)
    statistics = {}
    stats = calc_cc(merged_prism_file)
    statistics['correlation'] = stats

    if drop_pro:
        merged_dataframe_pro = merged_dataframe.copy()
        merged_dataframe_pro = merged_dataframe_pro[~merged_dataframe_pro.variant.str.contains(r".P", case=False, na=False, regex=True)]
        merged_dataframe_pro.reset_index(drop=True, inplace=True)
        merged_prism_file3 = os.path.join(folder.analysis, f'prims_merged_XXX_{sys_name}_merged_no-pro.txt')
        write_prism(merged_metadata, merged_dataframe_pro,
                    merged_prism_file3, comment=comment)

        stats_no_pro = calc_cc(
            merged_prism_file3)
        statistics['correlation_no_pro'] = stats_no_pro

    if drop_outliers:
        merged_dataframe_no_outlier = drop_numerical_outliers(
            merged_dataframe, variant_col="variant", score_col="predicted_ddG", z_thresh=3)
        merged_prism_file2 = os.path.join(folder.analysis, f'prims_merged_XXX_{sys_name}_merged_no-outlier.txt')
        write_prism(merged_metadata, merged_dataframe_no_outlier,
                    merged_prism_file2, comment=comment)

        stats_no_outlier = calc_cc(
            merged_prism_file2)
        statistics['correlation_no_outlier'] = stats_no_outlier

    statistic_outfile = os.path.join(folder.analysis, 'statistic.json')
    with open(statistic_outfile, 'w') as fp:
        json.dump(statistics, fp, indent=1, sort_keys=True)

if __name__ == '__main__':
    folder = AttrDict()
    folder.update({'input': sys.argv[1], 'output': sys.argv[
                  2], 'analysis': sys.argv[3]})

    if len(sys.argv) > 4:
        sys_name = sys.argv[4]
    else:
        sys_name = ''

    calc_all(folder, sys_name=sys_name)
