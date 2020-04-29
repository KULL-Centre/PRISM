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


DATA_DIR = os.path.join(PARENT_DIR, "data", "test")
TMP_DIR = os.path.join(DIR, "tmp")


def data(file_name):
    return os.path.join(DATA_DIR, file_name)


def tmp(*dir_name):
    return os.path.join(TMP_DIR, "mp-pipeline", *dir_name)


def clean_reference_from_local_path(dir_name, local_path):
    for dname, dirs, files in os.walk(dir_name):
        for fname in files:
            fpath = os.path.join(dname, fname)
            s = open(fpath, 'r').read()
            s = s.replace(local_path, '')
            open(fpath, "w").write(s)


class MPpipelineCreateTestCase(unittest.TestCase):
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
    >>> python -m unittest test_mp_pipeline.MPpipelineCreateTestCase
    >>> python -m unittest test_mp_pipeline.MPpipelineCreateTestCase.test_create_prov_flag
    """

    def setUp(self):
        # Creates the testrun directory
        self.global_test_dir = tmp()

    def tearDown(self):
        # Remove the directory after the test
        shutil.rmtree(self.global_test_dir, ignore_errors=True)

    def test_create_prov_flag(self):
        # initiates output and reference directory
        self.output_dir = tmp('create_prov_flag')
        self.reference_dir = data('mp-pipeline/output')

        # initiates the test
        inputdict = {
            'STRUC_FILE': data('mp-pipeline/input/PagP.pdb'),
            'UNIPROT_ID': '',
            'MUTATION_INPUT': data('mp-pipeline/input/PagP_short.mut'),
            'OUTPUT_FILE': self.output_dir,
            'DDG_FLAG_FILE': data('mp-pipeline/input/ddg_mp_flagfile'),
            'RELAX_FLAG_FILE': data('mp-pipeline/input/relax_mp_flagfile'),
            'MODE': 'create',
            'CHAIN': 'A',
            'RUN_STRUC': None,
            'LIGAND': None,
            'IS_MP': True,
            'MP_SPAN_INPUT': data('mp-pipeline/input/PagP.span'),
            'MP_CALC_SPAN_MODE': 'DSSP',
            'MP_THICKNESS': 15,
            'MP_ALIGN_REF': '3gp6',
            'MP_ALIGN_MODE': 'OPM',
            'RELAX_XML_INPUT': os.path.join(PARENT_DIR, 'data', 'mp', 'membrane_relax.xml'),
            'OVERWRITE_PATH': True,
            'SLURM_PARTITION': 'sbinlab',
        }
        input_args = Namespace(**inputdict)
        self.create = run_pipeline.predict_stability(input_args)

        # runs the test
        self.create

        # evaluates the test
        output_dic = {}
        for r, d, f in os.walk(self.output_dir):
            for file in f:
                parent_directory = r.split(self.output_dir)[1]
                # remove local paths within files
                file_read = open(os.path.join(r, file), 'r').read()
                file_read = file_read.replace(self.output_dir, '')
                file_read = file_read.replace(PARENT_DIR, '')

                output_dic[os.path.join(parent_directory, file)] = file_read

        reference_dic = {}
        for r, d, f in os.walk(self.reference_dir):
            for file in f:
                parent_directory = r.split(self.reference_dir)[1]
                # remove local paths within files
                clean_reference_from_local_path(
                    self.reference_dir, self.output_dir)
                clean_reference_from_local_path(self.reference_dir, PARENT_DIR)

                reference_dic[os.path.join(parent_directory, file)] = open(
                    os.path.join(r, file), 'r').read()
        self.maxDiff = None
        self.assertDictEqual(output_dic, reference_dic)
