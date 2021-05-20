#!/bin/sh 
#SBATCH --job-name=collect_rosetta_ddgs_6xro_final_renum 
#SBATCH --array=1 
#SBATCH --nodes=1 
#SBATCH --time=0:10:00 
#SBATCH --partition=sbinlab 

#This sbatch script launches the parse postprocess_rosetta_ddg_mp_pyrosetta function, from the mp_ddG 
python3 /mp_ddG.py /prepare/checking /ddG/run /ddG/output /ddG/input /output A 6xro_final_renum ddG.out False