#!/bin/sh 
#SBATCH --job-name=6xro_final_renum_MPddg
#SBATCH --array=0-2
#SBATCH --time=48:00:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab
#SBATCH --nice 
LST=(`ls /ddG/input/mutfiles/mutfile*`)
OFFSET=0 
INDEX=$((OFFSET+SLURM_ARRAY_TASK_ID))
echo $INDEX

# launching rosetta 
python3 /pyrosetta_ddG.py --in_pdb /ddG/input/input.pdb --in_span /ddG/input/spanfiles/input_A.span --outdir /ddG/run --outfile ddG.out --out_add ddG_additional.out --repack_radius 8.0 --include_pH 0 --pH_value 7 --repeats 5 --lowest 1 --lipids DLPC --temperature 20.0 --score_function franklin2019 --repack_protocol MP_flex_relax_ddG --lip_has_pore /prepare/mp_files/mp_lipid_acc/mp_lipid_acc_dic.json --dump_pdb 0 --mutfile ${LST[$INDEX]} 