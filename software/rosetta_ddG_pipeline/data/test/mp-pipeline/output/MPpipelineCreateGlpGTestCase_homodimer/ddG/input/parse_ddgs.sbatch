#!/bin/sh 
#SBATCH --job-name=homodimer-1afo_renum_uniquechain_p-ddg
#SBATCH --array=1 
#SBATCH --nodes=1 
#SBATCH --time=0:10:00 
#SBATCH --partition=sbinlab 

#This sbatch script launches the parse parse_rosetta_ddgs function, from the parse_cartesian_ddgs 
python3XXXtagvXXX1