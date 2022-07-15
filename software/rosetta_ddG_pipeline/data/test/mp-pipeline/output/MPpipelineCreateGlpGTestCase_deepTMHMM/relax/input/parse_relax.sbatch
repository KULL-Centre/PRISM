#!/bin/sh
#SBATCH --job-name=6xro_final_renum_p-relax
#SBATCH --time=00:20:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab

# launching parsing script 
python /relax_parse_results.py True True - /relax/input /relax/run /relax/output /ddG/input relax_scores