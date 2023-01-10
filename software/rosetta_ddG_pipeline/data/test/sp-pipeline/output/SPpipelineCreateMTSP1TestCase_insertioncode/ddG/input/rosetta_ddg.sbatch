#!/bin/sh 
#SBATCH --job-name=1EAW_ddg
#SBATCH --array=0-27
#SBATCH --time=48:00:00
#SBATCH --mem 2000
#SBATCH --partition=sbinlab
#SBATCH --nice 
LST=(`ls /ddG/input/mutfiles/mutfile*`)
OFFSET=0 
INDEX=$((OFFSET+SLURM_ARRAY_TASK_ID))
echo $INDEX

# launching rosetta 
if test -f "/ddG/input/input.pdb"; then
	bin/cartesian_ddg. -s /ddG/input/input.pdb -ddg:mut_file ${LST[$INDEX]}  -out:prefix ddg-$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID @/ddG/input/ddg_flagfile -ddg::dump_pdbs 0 
else
	echo "/ddG/input/input.pdb does not exist - exiting the call"
fi
