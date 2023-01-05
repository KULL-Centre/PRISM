#!/bin/sh 
#SBATCH --job-name=6xro_final_renum_MPcartddg
#SBATCH --array=0-130
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
	bin/cartesian_ddg. -database  -s /ddG/input/input.pdb -ddg:mut_file ${LST[$INDEX]} -score:weights /data/mp/f19_cart_1.5.wts -mp:lipids:composition DLPC -ddg:iterations 5 -ddg::dump_pdbs 0 -mp:setup:spanfiles /ddG/input/spanfiles/input_A.span -in:membrane -fa_max_dis 9 -missing_density_to_jump -has_pore 0 -ddg:legacy true -ddg:optimize_proline 1 -ddg:frag_nbrs 4 -ddg:bbnbrs 1 -ddg::cartesian -out:prefix ddg-$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID 
else
	echo "/ddG/input/input.pdb does not exist - exiting the call"
fi
