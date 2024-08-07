"""test_sp_pipeline.py tests the integrity of the sp pipeline.

Author: Johanna K.S. Tiemann

Date of last major changes: 2022-07-14

How to run all tests:
=======
>>> python -m unittest test_sp_pipeline
"""

# Standard library imports
from argparse import Namespace
from datetime import date
import os
import re
import shutil
import sys
import unittest

DIR = os.path.split(os.path.abspath(__file__))[0]
PARENT_DIR = os.path.split(DIR)[0]
sys.path.insert(1, sys.path.insert(0, PARENT_DIR))
os.environ['ddG_pipeline'] = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

# Local application imports
import run_pipeline
import rosetta_paths


DATA_DIR = os.path.join(PARENT_DIR, 'data', 'test')
TMP_DIR = os.path.join(DIR, 'tmp')


def data(file_name):
    return os.path.join(DATA_DIR, file_name)


def tmp(*dir_name):
    return os.path.join(TMP_DIR, 'sp-pipeline', *dir_name)


def clean_reference_from_local_path(dir_name, local_path):
    for dname, dirs, files in os.walk(dir_name):
        for fname in files:
            if not fname.startswith('.'):
                fpath = os.path.join(dname, fname)
                s = ''
                with open(fpath, 'r') as fp:
                    s = fp.read()
                if type(local_path)==list:
                    for elem in local_path:
                        s = s.replace(elem, '')
                else:
                    s = s.replace(local_path, '')
                with open(fpath, 'w') as fp:
                    fp.write(s)


def clean_version(dir_name, new_val=' XXXtagvXXX ', old_val=r'\s?.\S+tagv.*?\s'):
    for dname, dirs, files in os.walk(dir_name):
        for fname in files:
            if not fname.startswith('.'):
                fpath = os.path.join(dname, fname)
                s = ''
                with open(fpath, 'r') as fp:
                    s = fp.read()
                s = re.sub(old_val, new_val, s)
                with open(fpath, 'w') as fp:
                    fp.write(s)


class SPpipelineCreateDHFRTestCase(unittest.TestCase):
    """
    Description:
    =======
    Unittest for soluble protein pipeline

    Commands:
    =======
    - test_create_prov_flag

    Example:
    =======
    Running individual classes or methods
    >>> python -m unittest test_mp_pipeline.SPpipelineCreateDHFRTestCase
    >>> python -m unittest test_mp_pipeline.SPpipelineCreateDHFRTestCase.test_create_prov_flag
    """

    def setUp(self):
        # Creates the testrun directory
        self.global_test_dir = tmp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.global_test_dir, ignore_errors=True)
        # pass

    def test_create_rosettamutfile_prov_flag(self):
        # initiates output and reference directory
        self.output_dir = tmp('create_rosettamutfile_DHFR_prov_flag')
        self.reference_dir = data('sp-pipeline/output/SPpipelineCreateDHFRTestCase_rosettamutfile')

        # initiates the test
        inputdict = {
            'STRUC_FILE': data('sp-pipeline/input/4M6J.pdb'),
            'OUTPUT_FILE': self.output_dir,
            'MODE': 'create',
            'MUT_MODE': 'mut_file',
            'MUTATION_INPUT': data('sp-pipeline/input/mutfile_short'),
            'CHAIN': 'A',
            'RUN_STRUC': 'A',
            'LIGAND': None,
            'OVERWRITE_PATH': True,
            'SLURM_PARTITION': 'sbinlab',
            'GAPS_OUTPUT': False,
            'DUMP_PDB': 0,
            'DO_CHECKING': True,
            'ZIP_FILES': True,
            'NO_ZIP': False, 
            'VERBOSE': False,
            'IS_MP': False,
            'MP_SPAN_INPUT': None,
            'MP_CALC_SPAN_MODE': 'False',
            'MP_ALIGN_REF': '',
            'MP_ALIGN_MODE': 'False',
            'SUPERPOSE_ONTM': False,
            'DDG_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'cartesian_ddg_flagfile'),
            'RELAX_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'relax_flagfile'),
            'RELAX_XML_INPUT': os.path.join(rosetta_paths.path_to_data, 'mp', 'mp_relax.xml'),
            'UNIPROT_ID': '',
            'SCALE_FACTOR': 2.9,
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

        # evaluates the test
        skip_files = [
            'logs/info.log', 
            'span.log', 
        ]

        today = date.today()
        clean_version(self.output_dir)
        clean_version(self.output_dir, new_val='', old_val=today.strftime("%d-%b-%y").upper())
        output_dic = {}
        for r, d, f in os.walk(self.output_dir):
            for file in f:
                file_join = os.path.join(r, file)
                # skip specific files
                skip=False
                for skip_file in skip_files:
                    if file_join.endswith(skip_file):
                        skip=True
                if skip==False:
                    parent_directory = r.split(self.output_dir)[1]
                    # remove local paths within files
                    if not file.startswith('.'):
                        with open(file_join, 'r') as fp:
                            file_read = fp.read()
                            for elem in [self.output_dir, PARENT_DIR, rosetta_paths.ddG_pipeline,
                                         rosetta_paths.Rosetta_main_path, rosetta_paths.Rosetta_tools_path,
                                         rosetta_paths.Rosetta_database_path, rosetta_paths.Rosetta_extension]:
                                file_read = file_read.replace(elem, '')
                            output_dic[os.path.join(
                                parent_directory, file)] = file_read

        reference_dic = {}
        elem  = [self.output_dir, PARENT_DIR, rosetta_paths.ddG_pipeline,
                     rosetta_paths.Rosetta_main_path, rosetta_paths.Rosetta_tools_path,
                     rosetta_paths.Rosetta_database_path, rosetta_paths.Rosetta_extension]
        clean_reference_from_local_path(
            self.reference_dir, elem)
        clean_version(self.reference_dir)
        clean_version(self.reference_dir, new_val='', old_val=today.strftime("%d-%b-%y").upper())

        for r, d, f in os.walk(self.reference_dir):
            for file in f:
                file_join = os.path.join(r, file)
                # skip specific files
                skip=False
                for skip_file in skip_files:
                    if file_join.endswith(skip_file):
                        skip=True
                if skip==False:
                    parent_directory = r.split(self.reference_dir)[1]
                    # remove local paths within files
                    if not file.startswith('.'):
                        with open(file_join, 'r') as fp:
                            reference_dic[os.path.join(
                                parent_directory, file)] = fp.read()
        self.maxDiff = None
        # print('test_create_rosettamutfile_prov_flag')
        # for key in output_dic.keys():
        #     if output_dic[key]!=reference_dic[key]:
        #         print(key)
        self.assertDictEqual(output_dic, reference_dic)


    def test_create_ligand_prov_flag(self):
        # initiates output and reference directory
        self.output_dir = tmp('create_ligand_DHFR_prov_flag')
        self.reference_dir = data('sp-pipeline/output/SPpipelineCreateDHFRTestCase_ligand')

        # initiates the test
        inputdict = {
            'STRUC_FILE': data('sp-pipeline/input/4M6J.pdb'),
            'OUTPUT_FILE': self.output_dir,
            'MODE': 'create',
            'MUT_MODE': 'mut_file',
            'MUTATION_INPUT': data('sp-pipeline/input/mutfile_short'),
            'CHAIN': 'A',
            'RUN_STRUC': 'A',
            'LIGAND': True,
            'OVERWRITE_PATH': True,
            'SLURM_PARTITION': 'sbinlab',
            'GAPS_OUTPUT': False,
            'DUMP_PDB': 0,
            'DO_CHECKING': True,
            'ZIP_FILES': True,
            'NO_ZIP': False, 
            'VERBOSE': False,
            'IS_MP': False,
            'MP_SPAN_INPUT': None,
            'MP_CALC_SPAN_MODE': 'False',
            'MP_ALIGN_REF': '',
            'MP_ALIGN_MODE': 'False',
            'SUPERPOSE_ONTM': False,
            'DDG_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'cartesian_ddg_flagfile'),
            'RELAX_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'relax_flagfile'),
            'RELAX_XML_INPUT': os.path.join(rosetta_paths.path_to_data, 'mp', 'mp_relax.xml'),
            'UNIPROT_ID': '',
            'SCALE_FACTOR': 2.9,
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

        # evaluates the test
        skip_files = [
            'logs/info.log', 
            'span.log', 
        ]

        today = date.today()
        clean_version(self.output_dir)
        clean_version(self.output_dir, new_val='', old_val=today.strftime("%d-%b-%y").upper())
        output_dic = {}
        for r, d, f in os.walk(self.output_dir):
            for file in f:
                file_join = os.path.join(r, file)
                # skip specific files
                skip=False
                for skip_file in skip_files:
                    if file_join.endswith(skip_file):
                        skip=True
                if skip==False:
                    parent_directory = r.split(self.output_dir)[1]
                    # remove local paths within files
                    if not file.startswith('.'):
                        with open(file_join, 'r') as fp:
                            file_read = fp.read()
                            for elem in [self.output_dir, PARENT_DIR, rosetta_paths.ddG_pipeline,
                                         rosetta_paths.Rosetta_main_path, rosetta_paths.Rosetta_tools_path,
                                         rosetta_paths.Rosetta_database_path, rosetta_paths.Rosetta_extension]:
                                file_read = file_read.replace(elem, '')
                            output_dic[os.path.join(
                                parent_directory, file)] = file_read

        reference_dic = {}
        elem  = [self.output_dir, PARENT_DIR, rosetta_paths.ddG_pipeline,
                     rosetta_paths.Rosetta_main_path, rosetta_paths.Rosetta_tools_path,
                     rosetta_paths.Rosetta_database_path, rosetta_paths.Rosetta_extension]
        clean_reference_from_local_path(
            self.reference_dir, elem)
        clean_version(self.reference_dir)
        clean_version(self.reference_dir, new_val='\n', old_val=today.strftime("%d-%b-%y").upper())

        for r, d, f in os.walk(self.reference_dir):
            for file in f:
                file_join = os.path.join(r, file)
                # skip specific files
                skip=False
                for skip_file in skip_files:
                    if file_join.endswith(skip_file):
                        skip=True
                if skip==False:
                    parent_directory = r.split(self.reference_dir)[1]
                    # remove local paths within files
                    if not file.startswith('.'):
                        with open(file_join, 'r') as fp:
                            reference_dic[os.path.join(
                                parent_directory, file)] = fp.read()
        self.maxDiff = None
        # print('test_create_ligand_prov_flag')
        # for key in output_dic.keys():
        #     if output_dic[key]!=reference_dic[key]:
        #         print(key)
        self.assertDictEqual(output_dic, reference_dic)

        
class SPpipelineCreateMTSP1TestCase(unittest.TestCase):
    """
    Description:
    =======
    Unittest for soluble protein pipeline

    Commands:
    =======
    - test_create_prov_flag

    Example:
    =======
    Running individual classes or methods
    >>> python -m unittest test_mp_pipeline.SPpipelineCreateMTSP1TestCase
    >>> python -m unittest test_mp_pipeline.SPpipelineCreateMTSP1TestCase.test_create_insertioncode_prov_flag
    """

    def setUp(self):
        # Creates the testrun directory
        self.global_test_dir = tmp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.global_test_dir, ignore_errors=True)
        # pass

    def test_create_insertioncode_prov_flag(self):
        # initiates output and reference directory
        self.output_dir = tmp('create_rosettamutfile_MTSP1_prov_flag')
        self.reference_dir = data('sp-pipeline/output/SPpipelineCreateMTSP1TestCase_insertioncode')

        # initiates the test
        inputdict = {
            'STRUC_FILE': data('sp-pipeline/input/1EAW.pdb'),
            'OUTPUT_FILE': self.output_dir,
            'MODE': 'create',
            'MUT_MODE': 'mut_file',
            'MUTATION_INPUT': data('sp-pipeline/input/mutfile_short_1EAW_AB_AB'),
            'CHAIN': 'AB',
            'RUN_STRUC': 'AB',
            'LIGAND': None,
            'OVERWRITE_PATH': True,
            'SLURM_PARTITION': 'sbinlab',
            'GAPS_OUTPUT': False,
            'DUMP_PDB': 0,
            'DO_CHECKING': True,
            'ZIP_FILES': True,
            'NO_ZIP': False, 
            'VERBOSE': False,
            'IS_MP': False,
            'MP_SPAN_INPUT': None,
            'MP_CALC_SPAN_MODE': 'False',
            'MP_ALIGN_REF': '',
            'MP_ALIGN_MODE': 'False',
            'SUPERPOSE_ONTM': False,
            'DDG_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'cartesian_ddg_flagfile'),
            'RELAX_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'relax_flagfile'),
            'RELAX_XML_INPUT': os.path.join(rosetta_paths.path_to_data, 'mp', 'mp_relax.xml'),
            'UNIPROT_ID': '',
            'SCALE_FACTOR': 2.9,
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

        # evaluates the test
        skip_files = [
            'logs/info.log', 
            'span.log', 
        ]

        today = date.today()
        clean_version(self.output_dir)
        clean_version(self.output_dir, new_val='', old_val=today.strftime("%d-%b-%y").upper())
        output_dic = {}
        for r, d, f in os.walk(self.output_dir):
            for file in f:
                file_join = os.path.join(r, file)
                # skip specific files
                skip=False
                for skip_file in skip_files:
                    if file_join.endswith(skip_file):
                        skip=True
                if skip==False:
                    parent_directory = r.split(self.output_dir)[1]
                    # remove local paths within files
                    if not file.startswith('.'):
                        with open(file_join, 'r') as fp:
                            file_read = fp.read()
                            for elem in [self.output_dir, PARENT_DIR, rosetta_paths.ddG_pipeline,
                                         rosetta_paths.Rosetta_main_path, rosetta_paths.Rosetta_tools_path,
                                         rosetta_paths.Rosetta_database_path, rosetta_paths.Rosetta_extension]:
                                file_read = file_read.replace(elem, '')
                            output_dic[os.path.join(
                                parent_directory, file)] = file_read

        reference_dic = {}
        elem  = [self.output_dir, PARENT_DIR, rosetta_paths.ddG_pipeline,
                     rosetta_paths.Rosetta_main_path, rosetta_paths.Rosetta_tools_path,
                     rosetta_paths.Rosetta_database_path, rosetta_paths.Rosetta_extension]
        clean_reference_from_local_path(
            self.reference_dir, elem)
        clean_version(self.reference_dir)
        clean_version(self.reference_dir, new_val='', old_val=today.strftime("%d-%b-%y").upper())

        for r, d, f in os.walk(self.reference_dir):
            for file in f:
                file_join = os.path.join(r, file)
                # skip specific files
                skip=False
                for skip_file in skip_files:
                    if file_join.endswith(skip_file):
                        skip=True
                if skip==False:
                    parent_directory = r.split(self.reference_dir)[1]
                    # remove local paths within files
                    if not file.startswith('.'):
                        with open(file_join, 'r') as fp:
                            reference_dic[os.path.join(
                                parent_directory, file)] = fp.read()
        self.maxDiff = None
        # print('test_create_insertioncode_prov_flag')
        # for key in output_dic.keys():
        #     if (key in output_dic.keys()) and (key in reference_dic.keys()):
        #         if output_dic[key]!=reference_dic[key]:
        #             print(key)
        #     else:
        #         print('Missing', key)
        self.assertDictEqual(output_dic, reference_dic)
