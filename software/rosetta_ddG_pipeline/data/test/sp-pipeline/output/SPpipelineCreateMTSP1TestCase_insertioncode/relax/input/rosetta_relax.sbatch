#!/bin/sh
#SBATCH --job-name=1EAW_relax
#SBATCH --time=24:00:00
#SBATCH --mem 5000
#SBATCH --array=0-19
#SBATCH --partition=sbinlab

# launching rosetta relax 
bin/relax. -s /relax/input/input.pdb -relax:script /data/sp/cart2.script @/relax/input/relax_flagfile -out:prefix $SLURM_ARRAY_TASK_ID-