#!/bin/sh 
#SBATCH --job-name=cartesian_PagP
#SBATCH --array=0-1
#SBATCH --time=32:00:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab
#SBATCH --nice 
LST=(`ls /ddG/input/mutfiles/mutfile*`)
OFFSET=0 
INDEX=$((OFFSET+SLURM_ARRAY_TASK_ID))
echo $INDEX

# launching rosetta 
/home/j0kaso/repos/Rosetta/main/source/bin/cartesian_ddg.default.linuxgccdebug -s /ddG/input/input.pdb -ddg:mut_file ${LST[$INDEX]} -out:prefix ddg-$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID @/ddG/input/ddg_flagfile