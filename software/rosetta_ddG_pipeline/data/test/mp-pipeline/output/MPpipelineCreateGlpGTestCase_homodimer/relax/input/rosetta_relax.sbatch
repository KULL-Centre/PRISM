#!/bin/sh
#SBATCH --job-name=relax_homodimer-1afo_renum_uniquechain
#SBATCH --time=24:00:00
#SBATCH --mem 5000
#SBATCH --array=0-19
#SBATCH --partition=sbinlab

# launching rosetta relax 
bin/rosetta_scripts. -parser:protocol /relax/input/relax.xml -parser:script_vars repeats=5 energy_func=/data/mp/f19_cart_1.5.wts energy_fawtb=0 -in:file:s /relax/input/input.pdb -mp:setup:spanfiles /relax/input/spanfiles/input_A.span -mp:scoring:hbond -mp:lipids:composition DLPC -mp::thickness 15 -relax:jump_move true -nstruct 1 -packing:pack_missing_sidechains 0 -out:pdb -out:file:scorefile /relax/run/relax_scores.sc -out:prefix $SLURM_ARRAY_TASK_ID- -database  -ignore_unrecognized_res true -score:weights /data/mp/f19_cart_1.5.wts -fa_max_dis 9 -ex1 -ex2 -flip_HNQ -missing_density_to_jump -relax:coord_constrain_sidechains -relax:constrain_relax_to_start_coords  