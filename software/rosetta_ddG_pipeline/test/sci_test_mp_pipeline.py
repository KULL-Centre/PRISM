"""sci_test_mp_pipeline.py tests the complete integrity of the mp pipeline.

Author: Johanna K.S. Tiemann

Date of last major changes: 2021-05-20

How to run all tests:
=======
>>> python -m unittest sci_test_mp_pipeline


consider changing to component test or similar instead of unittest
"""

# Standard library imports
from argparse import Namespace
from datetime import timedelta
import os
import shutil
import sys
import time
import unittest

# 3rd party library imports
import numpy as np
from scipy.stats import linregress
import pandas as pd

DIR = os.path.split(os.path.abspath(__file__))[0]
PARENT_DIR = os.path.split(DIR)[0]
sys.path.insert(1, sys.path.insert(0, PARENT_DIR))
os.environ['ddG_pipeline'] = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

# Local application imports
import run_pipeline
import rosetta_paths
from prism_rosetta_parser import read_prism


DATA_DIR = os.path.join(PARENT_DIR, 'data', 'test')
TMP_DIR = os.path.join(DIR, 'tmp')


def data(file_name):
    return os.path.join(DATA_DIR, file_name)


def tmp(*dir_name):
    return os.path.join(TMP_DIR, 'sci-mp-pipeline', *dir_name)


class MPpipelineFullrunGlpGTestCase(unittest.TestCase):
    """
    Description:
    =======
    Unittest for membrane protein pipeline

    Commands:
    =======
    - test_fullrun_rosettamutfile_prov_flag

    Example:
    =======
    Running individual classes or methods
    >>> python -m unittest sci_test_mp_pipeline.MPpipelineFullrunGlpGTestCase
    >>> python -m unittest sci_test_mp_pipeline.MPpipelineFullrunGlpGTestCase.test_fullrun_rosettamutfile_prov_flag
    """

    def setUp(self):
        # Creates the testrun directory
        self.global_test_dir = tmp()

    def tearDown(self):
        # Remove the directory after the test
        #shutil.rmtree(self.global_test_dir, ignore_errors=True)
        print('not automatic deletion - please perform manually')

    def test_a_fullrun_rosettamutfile_prov_flag(self):
        # initiates output and reference directory
        self.output_dir = tmp('fullrun_rosettamutfile_GlpG_prov_flag')
        self.reference_dir = data('mp-pipeline/output/MPpipelineFullrunGlpGTestCase_rosettamutfile')

        # initiates the test
        inputdict = {
            'STRUC_FILE': data('mp-pipeline/input/6xro_final_renum.pdb'),
            'OUTPUT_FILE': self.output_dir,
            'MODE': 'fullrun',
            'MUT_MODE': 'mut_file',
            'MUTATION_INPUT': data('mp-pipeline/input/mutfile_all'),
            'CHAIN': 'A',
            'RUN_STRUC': None,
            'LIGAND': None,
            'OVERWRITE_PATH': True,
            'SLURM_PARTITION': 'sbinlab',
            'GAPS_OUTPUT': False,
            'DUMP_PDB': 0,
            'DO_CHECKING': True,
            'NO_ZIP': False, 
            'ZIP_FILES': True,
            'VERBOSE': False,
            'IS_MP': True,
            'MP_SPAN_INPUT': None,
            'MP_CALC_SPAN_MODE': 'deepTMHMM', #'DSSP',
            'MP_ALIGN_REF': '-', #'6xro_A',
            'MP_ALIGN_MODE': 'span', #'OPM',
            'SUPERPOSE_ONTM': True,
            'DDG_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'cartesian_ddg_flagfile'),
            'RELAX_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'relax_flagfile'),
            'RELAX_XML_INPUT': os.path.join(rosetta_paths.path_to_data, 'mp', 'mp_cart_relax.xml'),
            'UNIPROT_ID': '',
            'SCALE_FACTOR': 1,
            'MP_THICKNESS': 15,
            'MP_LIPIDS': 'DLPC',
            'MP_TEMPERATURE': 20.0,
            'MP_PH': -1.0,
            'BENCH_MP_REPACK': 8.0,
            'BENCH_MP_REPEAT': 5,
            'BENCH_MP_RELAX_REPEAT': 5,
            'BENCH_MP_RELAX_STRUCS': 20,
            'MP_IGNORE_RELAX_MP_FLAGS': False,
            'MP_ENERGY_FUNC': 'franklin2019',
            'MP_ENERGY_FUNC_WEIGHTS': os.path.join(rosetta_paths.path_to_data, 'mp', 'f19_cart_1.5.wts'),
            'MP_REPACK_PROTOCOL': 'MP_flex_relax_ddG',
            'MP_MULTISTRUC_PROTOCOL': 0 ,
            'MP_CART_DDG': 1,
        }
        input_args = Namespace(**inputdict)
        self.create = run_pipeline.predict_stability(input_args)

        # runs the test
        self.create

    def test_b_analyze_exp_rosettamutfile_prov_flag(self):
        self.output_dir = tmp('fullrun_rosettamutfile_GlpG_prov_flag')
        self.reference_dir_in = data('mp-pipeline/input')
        self.reference_dir_out = data('mp-pipeline/output/MPpipelineFullrunGlpGTestCase_rosettamutfile')


        prism_file_test = os.path.join(self.output_dir, 'output', 'prism_rosetta_XXX_6xro_final_renum.txt')
        date_time = pd.Timestamp.today()
        date_time += timedelta(hours = 2)
        if not os.path.isfile(prism_file_test):
            print(f"Please wait for ~ 2h ({date_time.strftime('%Y-%m-%d %H:%M')}) for the calculation to finish.")
        while not os.path.isfile(prism_file_test):
            time.sleep(20)


        correlation_list_file = os.path.join(self.reference_dir_out, 'GlpG_cutoffs.txt')
        if os.path.isfile(correlation_list_file):
            correlation_all_list = pd.read_csv(correlation_list_file, header=[0])
        else:
            correlation_all_list = pd.DataFrame(data=[], columns=['date', 'exp', 'corr', 'slope', 'intercept', 'count'])


        cutoffs_ref = dict()
        for exp in correlation_all_list['exp'].unique():
            tmp_df = correlation_all_list.copy()
            tmp_df = tmp_df.loc[(tmp_df['exp']==exp)].reset_index(drop=True)
            cutoffs_ref[exp] = tmp_df['corr'].min()-(tmp_df['corr'].min()*0.1)
            #cutoffs_ref[exp] = tmp_df['corr'].median()-(tmp_df['corr'].median()*0.25)
            #cutoffs_ref[exp] = tmp_df['corr'].median()-0.2


        dataframe_test = read_prism(prism_file_test).dataframe[['variant', 'mean_ddG']]
        prism_file_exp = os.path.join(self.reference_dir_in, 'prism_merged_files.txt')
        dataframe_exp = read_prism(prism_file_exp).dataframe

        dataframe_test = dataframe_test.rename(columns={'mean_ddG':'new'})
        ddG_column = [column for column in dataframe_exp.columns if column.startswith('ddG')]
        dataframe_exp = pd.merge(dataframe_test[['variant', 'new']], dataframe_exp[['variant']+ddG_column], on='variant', suffixes=('',''))


        correlations_list = [['date', 'exp', 'corr', 'slope', 'intercept', 'count']]
        date_time = pd.Timestamp.today().strftime('%Y-%m-%d-%H-%M')
        exp_list = []
        ref_list = []
        cutoffs_test = dict()
        for exp in ddG_column:
            tmp_df = dataframe_exp.copy()
            tmp_df = tmp_df[['new', exp]].dropna(how='any').reset_index(drop=True)
            exp_list += list(tmp_df[exp])
            ref_list += list(tmp_df['new'])

            linreg = linregress(tmp_df[exp], tmp_df['new'])
            corr = round( linreg.rvalue, 3 )
            slope = round( linreg.slope, 3 )
            intercept = round( linreg.intercept, 3 )
            correlations_list.append([date_time, exp, corr, slope, intercept, len(tmp_df[exp])])
            cutoffs_test[exp] = corr
        linreg = linregress(exp_list, ref_list)
        corr = round( linreg.rvalue, 3 )
        slope = round( linreg.slope, 3 )
        intercept = round( linreg.intercept, 3 )
        correlations_list.append([date_time, 'all', corr, slope, intercept, len(exp_list)])
        cutoffs_test['all'] = corr
        correlations_list = pd.DataFrame(data=correlations_list[1:], columns=correlations_list[0])
        correlation_merged = pd.concat([correlation_all_list, correlations_list], ignore_index=True)
        correlation_merged.to_csv(correlation_list_file, index=False)
        print(correlation_merged)


        result_dic = dict()
        for key in cutoffs_ref.keys():
            if abs(cutoffs_test[key]) > cutoffs_ref[key]:
                result_dic[key] = True
            else:
                result_dic[key] = False
        print(f"{result_dic}, {cutoffs_test}, {cutoffs_ref}")
        
        for key in cutoffs_ref.keys():
            self.assertTrue(abs(cutoffs_test[key]) > cutoffs_ref[key])

    def test_c_analyze_comp_rosettamutfile_prov_flag(self):
        self.output_dir = tmp('fullrun_rosettamutfile_GlpG_prov_flag')
        self.reference_dir_out = data('mp-pipeline/output/MPpipelineFullrunGlpGTestCase_rosettamutfile')

        ref_ddG_file = os.path.join(self.reference_dir_out , 'ref_ddG.txt')
        ref_ddG_df = pd.read_csv(ref_ddG_file, header=[0])
        ref_ddG_df['median_refs'] = ref_ddG_df.median(axis = 1)

        prism_file1 = os.path.join(self.output_dir, 'output', 'prism_rosetta_XXX_6xro_final_renum.txt')
        dataframe1 = read_prism(prism_file1).dataframe[['variant', 'mean_ddG']]
        dataframe1['mean_ddG'] = dataframe1['mean_ddG'].round(3)
        date_time = pd.Timestamp.today().strftime('%Y-%m-%d-%H-%M')
        df = pd.merge(ref_ddG_df, dataframe1[['variant', 'mean_ddG']], on='variant')
        df = df.rename(columns={'mean_ddG':date_time})

        linreg = linregress(df['median_refs'], df[date_time])
        corr = round( linreg.rvalue, 3 )
        slope = round( linreg.slope, 3 )
        intercept = round( linreg.intercept, 3 )
        
        df = df.drop(columns=['median_refs'])
        df.to_csv(ref_ddG_file, index=False)
        print(f"corr: {corr}, slope: {slope}, intercept: {intercept}")
        
        self.assertTrue(abs(corr) > 0.85)

    def test_d_deleteFolder(self):
        # Remove the directory after the test
        shutil.rmtree(self.global_test_dir, ignore_errors=True)
        print('Folder deleted')

