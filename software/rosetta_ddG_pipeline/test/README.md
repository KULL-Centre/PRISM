# Description:
=======
Unittest for stability protein pipeline

# Consistency test for soluble pipeline (only prepare is tested):
=======
Running all tests 
> cd ./test

> python -m unittest

> python -m unittest discover

Running individual files (all tests of this script)
> python -m unittest test_sp_pipeline

Running individual classes or methods
> python -m unittest test_sp_pipeline.SPpipelineCreateDHFRTestCase

> python -m unittest test_sp_pipeline.SPpipelineCreateDHFRTestCase.test_create_rosettamutfile_prov_flag

> python -m unittest test_sp_pipeline.SPpipelineCreateDHFRTestCase.test_create_ligand_prov_flag

Please expect some waiting time for the last test.


# Consistency test for membrane pipeline:
=======
Note: make sure to activate the python environment before running
> conda activate /groups/sbinlab/software/PRISM_tools/py3_ros_ddG_env

Running all tests 
> cd ./test

> python -m unittest

> python -m unittest discover

Running individual files (all tests of this script)
> python -m unittest test_mp_pipeline

Running individual classes or methods
> python -m unittest test_mp_pipeline.MPpipelineCreateGlpGTestCase

> python -m unittest test_mp_pipeline.MPpipelineCreateGlpGTestCase.test_create_rosettamutfile_prov_flag

> python -m unittest test_mp_pipeline.MPpipelineCreateGlpGTestCase.test_create_pipemut_prov_flag

> python -m unittest test_mp_pipeline.MPpipelineCreateGlpGTestCase.test_create_mutfileall_prov_flag

> python -m unittest test_mp_pipeline.MPpipelineCreateGlpGTestCase.test_create_mutdir_prov_flag

> python -m unittest test_mp_pipeline.MPpipelineCreateGlpGTestCase.test_create_homodimer_prov_flag

> python -m unittest test_mp_pipeline.MPpipelineCreateGlpGTestCase.test_create_deepTMHMM_prov_flag

Please expect some waiting time for the last test.


# Scientific tests mp (full runs - cost time and resources):
=======
Run all and wait
> python -m unittest sci_test_mp_pipeline.MPpipelineFullrunGlpGTestCase

Or run in this order
> python -m unittest sci_test_mp_pipeline.MPpipelineFullrunGlpGTestCase.test_a_fullrun_rosettamutfile_prov_flag

wait until calculations have finished (ca. 2h)


This algorithm runs on 154 CPU (any) for ~ 0 h and 37 min each, for a total runtime of 89 h and 10 min, with a memory usage of 78.13 GB, drawing 1.82 kWh. Based in Denmark, this program produces 280.56 g of CO2e, which is equivalent to 0.31 tree-months carbon sequestration, a train ride of 6.84 km, 1.60 km of driving a passenger car in Europe or 3.26 h of netflix streaming.

> python -m unittest sci_test_mp_pipeline.MPpipelineFullrunGlpGTestCase.test_b_analyze_exp_rosettamutfile_prov_flag

> python -m unittest sci_test_mp_pipeline.MPpipelineFullrunGlpGTestCase.test_c_analyze_comp_rosettamutfile_prov_flag

> python -m unittest sci_test_mp_pipeline.MPpipelineFullrunGlpGTestCase.test_d_deleteFolder


# Scientific tests sp (full runs - cost time and resources):
=======
Run all and wait
> python -m unittest sci_test_sp_pipeline.SPpipelineFullrunDHFRTestCase

Or run in this order

> python -m unittest sci_test_sp_pipeline.SPpipelineFullrunDHFRTestCase.test_a_fullrun_rosettamutfile_prov_flag

wait until calculations have finished (ca. 2h)


This algorithm runs on 60 CPU (any) for ~ 0 h and 27 min each, for a total runtime of 20 h and 58 min, with a memory usage of 10.05 GB, drawing 0.42 kWh. Based in Denmark, this program produces 65.30 g of CO2e, which is equivalent to 0.07 tree-months carbon sequestration, a train ride of 1.59 km, 0.37 km of driving a passenger car in Europe or 0.76 h of netflix streaming.

> python -m unittest sci_test_sp_pipeline.SPpipelineFullrunDHFRTestCase.test_c_analyze_comp_rosettamutfile_prov_flag

> python -m unittest sci_test_sp_pipeline.SPpipelineFullrunDHFRTestCase.test_d_deleteFolder


# Trouble shooting:
=======
Sometimes the pybiolib library needs to be updated. Execute one (or multiple) of the following commands:

> pip3 install -U pybiolib

> conda install -y -c conda-forge -c anaconda -c defaults -c bioconda biolib