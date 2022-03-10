# Description:
=======
Unittest for stability protein pipeline

# Example:
=======
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

> python -m unittest test_mp_pipeline.MPpipelineCreateGlpGTestCase.test_create_mutdir_prov_flag

> python -m unittest test_mp_pipeline.MPpipelineCreateGlpGTestCase.test_create_homodimer_prov_flag

Please expect some waiting time for the last test.


# Scientific tests mp (full runs - cost time and resources):
=======
Run all and wait
> python -m unittest sci_test_mp_pipeline.MPpipelineFullrunGlpGTestCase

Or run in this order
> python -m unittest sci_test_mp_pipeline.MPpipelineFullrunGlpGTestCase.test_a_fullrun_rosettamutfile_prov_flag

wait until calculations have finished (~2h)

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

wait until calculations have finished (~2h)

> python -m unittest sci_test_sp_pipeline.SPpipelineFullrunDHFRTestCase.test_c_analyze_comp_rosettamutfile_prov_flag

> python -m unittest sci_test_sp_pipeline.SPpipelineFullrunDHFRTestCase.test_d_deleteFolder
