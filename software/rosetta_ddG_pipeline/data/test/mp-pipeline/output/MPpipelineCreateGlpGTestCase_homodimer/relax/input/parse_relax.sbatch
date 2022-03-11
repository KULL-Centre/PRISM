#!/bin/sh
#SBATCH --job-name=homodimer-1afo_renum_uniquechain_p-relax
#SBATCH --time=00:20:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab

# launching parsing script 
python /relax_parse_results.py True 1AFO_A /input /relax/run /relax/output /ddG/input relax_scores