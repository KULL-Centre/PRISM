"""test_mp_pipeline.py tests the integrity of the mp pipeline.

Author: Johanna K.S. Tiemann

Date of last major changes: 2020-04-29

How to run all tests:
=======
>>> python -m unittest test_mp_pipeline
"""

# Standard library imports
from argparse import Namespace
import os
import shutil
import sys
import unittest

DIR = os.path.split(os.path.abspath(__file__))[0]
PARENT_DIR = os.path.split(DIR)[0]
sys.path.insert(1, sys.path.insert(0, PARENT_DIR))

# Local application imports
import run_pipeline
import rosetta_paths


DATA_DIR = os.path.join(PARENT_DIR, 'data', 'test')
TMP_DIR = os.path.join(DIR, 'tmp')


def data(file_name):
    return os.path.join(DATA_DIR, file_name)


def tmp(*dir_name):
    return os.path.join(TMP_DIR, 'mp-pipeline', *dir_name)


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


class MPpipelineCreateGlpGTestCase(unittest.TestCase):
    """
    Description:
    =======
    Unittest for membrane protein pipeline

    Commands:
    =======
    - test_create_prov_flag

    Example:
    =======
    Running individual classes or methods
    >>> python -m unittest test_mp_pipeline.MPpipelineCreateGlpGTestCase
    >>> python -m unittest test_mp_pipeline.MPpipelineCreateGlpGTestCase.test_create_prov_flag
    """

    def setUp(self):
        # Creates the testrun directory
        self.global_test_dir = tmp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.global_test_dir, ignore_errors=True)

    def test_create_rosettamutfile_prov_flag(self):
        # initiates output and reference directory
        self.output_dir = tmp('create_rosettamutfile_GlpG_prov_flag')
        self.reference_dir = data('mp-pipeline/output/MPpipelineCreateGlpGTestCase_rosettamutfile')

        # initiates the test
        inputdict = {
            'STRUC_FILE': data('mp-pipeline/input/6xro_final_renum.pdb'),
            'OUTPUT_FILE': self.output_dir,
            'MODE': 'create',
            'MUT_MODE': 'mut_file',
            'MUTATION_INPUT': data('mp-pipeline/input/mutfile_short'),
            'CHAIN': 'A',
            'RUN_STRUC': None,
            'LIGAND': None,
            'OVERWRITE_PATH': True,
            'SLURM_PARTITION': 'sbinlab',
            'GAPS_OUTPUT': False,
            'DUMP_PDB': 0,
            'VERBOSE': False,
            'IS_MP': True,
            'MP_SPAN_INPUT': None,
            'MP_CALC_SPAN_MODE': 'DSSP',
            'MP_ALIGN_REF': '6xro_A',
            'MP_ALIGN_MODE': 'OPM',
            'DDG_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'cartesian_ddg_flagfile'),
            'RELAX_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'relax_flagfile'),
            'RELAX_XML_INPUT': os.path.join(rosetta_paths.path_to_data, 'mp', 'mp_relax.xml'),
            'UNIPROT_ID': '',
            'MP_THICKNESS': 15,
            'MP_LIPIDS': 'DLPC',
            'MP_TEMPERATURE': 37.0,
            'MP_PH': -1.0,
            'BENCH_MP_REPACK': 8.0,
            'BENCH_MP_REPEAT': 5,
            'BENCH_MP_RELAX_REPEAT': 5,
            'BENCH_MP_RELAX_STRUCS': 20,
            'MP_IGNORE_RELAX_MP_FLAGS': False,
            'MP_ENERGY_FUNC': 'franklin2019',
            'MP_REPACK_PROTOCOL': 'MP_flex_relax_ddG',
            'MP_MULTISTRUC_PROTOCOL': 0 ,
        }
        input_args = Namespace(**inputdict)
        self.create = run_pipeline.predict_stability(input_args)

        # runs the test
        self.create

        # evaluates the test
        skip_files = ['logs/info.log', 'span.log', 
            'mp_lipid_acc/input_A_0001.pdb',
        ]

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
        
        self.assertDictEqual(output_dic, reference_dic)



    def test_create_pipemut_prov_flag(self):
        # initiates output and reference directory
        self.output_dir = tmp('create_pipemutfile_GlpG_prov_flag')
        self.reference_dir = data('mp-pipeline/output/MPpipelineCreateGlpGTestCase_pipemutfile')

        # initiates the test
        inputdict = {
            'STRUC_FILE': data('mp-pipeline/input/6xro_final_renum.pdb'),
            'OUTPUT_FILE': self.output_dir,
            'MODE': 'create',
            'MUT_MODE': 'mut_file',
            'MUTATION_INPUT': data('mp-pipeline/input/pipeline_mutfile_short.txt'),
            'CHAIN': 'A',
            'RUN_STRUC': None,
            'LIGAND': None,
            'OVERWRITE_PATH': True,
            'SLURM_PARTITION': 'sbinlab',
            'GAPS_OUTPUT': False,
            'DUMP_PDB': 0,
            'VERBOSE': False,
            'IS_MP': True,
            'MP_SPAN_INPUT': None,
            'MP_CALC_SPAN_MODE': 'DSSP',
            'MP_ALIGN_REF': '6xro_A',
            'MP_ALIGN_MODE': 'OPM',
            'DDG_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'cartesian_ddg_flagfile'),
            'RELAX_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'relax_flagfile'),
            'RELAX_XML_INPUT': os.path.join(rosetta_paths.path_to_data, 'mp', 'mp_relax.xml'),
            'UNIPROT_ID': '',
            'MP_THICKNESS': 15,
            'MP_LIPIDS': 'DLPC',
            'MP_TEMPERATURE': 37.0,
            'MP_PH': -1.0,
            'BENCH_MP_REPACK': 8.0,
            'BENCH_MP_REPEAT': 5,
            'BENCH_MP_RELAX_REPEAT': 5,
            'BENCH_MP_RELAX_STRUCS': 20,
            'MP_IGNORE_RELAX_MP_FLAGS': False,
            'MP_ENERGY_FUNC': 'franklin2019',
            'MP_REPACK_PROTOCOL': 'MP_flex_relax_ddG',
            'MP_MULTISTRUC_PROTOCOL': 0 ,
        }
        input_args = Namespace(**inputdict)
        self.create = run_pipeline.predict_stability(input_args)

        # runs the test
        self.create

        # evaluates the test
        skip_files = ['logs/info.log', 'span.log', 
            'mp_lipid_acc/input_A_0001.pdb',
        ]

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
        
        self.assertDictEqual(output_dic, reference_dic)


    def test_create_mutdir_prov_flag(self):
        # initiates output and reference directory
        self.output_dir = tmp('create_mutdir_GlpG_prov_flag')
        self.reference_dir = data('mp-pipeline/output/MPpipelineCreateGlpGTestCase_mutdir')

        # initiates the test
        inputdict = {
            'STRUC_FILE': data('mp-pipeline/input/6xro_final_renum.pdb'),
            'OUTPUT_FILE': self.output_dir,
            'MODE': 'create',
            'MUT_MODE': 'mut_file',
            'MUTATION_INPUT': data('mp-pipeline/input/mut_files'),
            'CHAIN': 'A',
            'RUN_STRUC': None,
            'LIGAND': None,
            'OVERWRITE_PATH': True,
            'SLURM_PARTITION': 'sbinlab',
            'GAPS_OUTPUT': False,
            'DUMP_PDB': 0,
            'VERBOSE': False,
            'IS_MP': True,
            'MP_SPAN_INPUT': None,
            'MP_CALC_SPAN_MODE': 'DSSP',
            'MP_ALIGN_REF': '6xro_A',
            'MP_ALIGN_MODE': 'OPM',
            'DDG_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'cartesian_ddg_flagfile'),
            'RELAX_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'relax_flagfile'),
            'RELAX_XML_INPUT': os.path.join(rosetta_paths.path_to_data, 'mp', 'mp_relax.xml'),
            'UNIPROT_ID': '',
            'MP_THICKNESS': 15,
            'MP_LIPIDS': 'DLPC',
            'MP_TEMPERATURE': 37.0,
            'MP_PH': -1.0,
            'BENCH_MP_REPACK': 8.0,
            'BENCH_MP_REPEAT': 5,
            'BENCH_MP_RELAX_REPEAT': 5,
            'BENCH_MP_RELAX_STRUCS': 20,
            'MP_IGNORE_RELAX_MP_FLAGS': False,
            'MP_ENERGY_FUNC': 'franklin2019',
            'MP_REPACK_PROTOCOL': 'MP_flex_relax_ddG',
            'MP_MULTISTRUC_PROTOCOL': 0 ,
        }
        input_args = Namespace(**inputdict)
        self.create = run_pipeline.predict_stability(input_args)

        # runs the test
        self.create

        # evaluates the test
        skip_files = ['logs/info.log', 'span.log', 
            'mp_lipid_acc/input_A_0001.pdb', 'mutation_clean.txt',
        ]

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
        
        self.assertDictEqual(output_dic, reference_dic)



    def test_create_homodimer_prov_flag(self):
        # initiates output and reference directory
        self.output_dir = tmp('create_homodimer_GlpG_prov_flag')
        self.reference_dir = data('mp-pipeline/output/MPpipelineCreateGlpGTestCase_homodimer')

        # initiates the test
        inputdict = {
            'STRUC_FILE': data('mp-pipeline/input/homodimer-1afo_renum_uniquechain.pdb'),
            'OUTPUT_FILE': self.output_dir,
            'MODE': 'create',
            'MUT_MODE': 'mut_file',
            'MUTATION_INPUT': data('mp-pipeline/input/homodimer-mutations'),
            'CHAIN': 'A',
            'RUN_STRUC': None,
            'LIGAND': None,
            'OVERWRITE_PATH': True,
            'SLURM_PARTITION': 'sbinlab',
            'GAPS_OUTPUT': False,
            'DUMP_PDB': 0,
            'VERBOSE': False,
            'IS_MP': True,
            'MP_SPAN_INPUT': None,
            'MP_CALC_SPAN_MODE': 'DSSP',
            'MP_ALIGN_REF': '1AFO_A',
            'MP_ALIGN_MODE': 'OPM',
            'DDG_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'cartesian_ddg_flagfile'),
            'RELAX_FLAG_FILE': os.path.join(rosetta_paths.path_to_data, 'sp', 'relax_flagfile'),
            'RELAX_XML_INPUT': os.path.join(rosetta_paths.path_to_data, 'mp', 'mp_relax.xml'),
            'UNIPROT_ID': '',
            'MP_THICKNESS': 15,
            'MP_LIPIDS': 'DLPC',
            'MP_TEMPERATURE': 37.0,
            'MP_PH': -1.0,
            'BENCH_MP_REPACK': 8.0,
            'BENCH_MP_REPEAT': 5,
            'BENCH_MP_RELAX_REPEAT': 5,
            'BENCH_MP_RELAX_STRUCS': 20,
            'MP_IGNORE_RELAX_MP_FLAGS': False,
            'MP_ENERGY_FUNC': 'franklin2019',
            'MP_REPACK_PROTOCOL': 'MP_flex_relax_ddG',
            'MP_MULTISTRUC_PROTOCOL': 0 ,
        }
        input_args = Namespace(**inputdict)
        self.create = run_pipeline.predict_stability(input_args)

        # runs the test
        self.create

        # evaluates the test
        skip_files = ['logs/info.log', 'span.log', 
            'mp_lipid_acc/input_A_0001.pdb', 'mutation_clean.txt',
        ]

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
        
        self.assertDictEqual(output_dic, reference_dic)

