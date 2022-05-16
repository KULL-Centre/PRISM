#!/bin/sh
#SBATCH --job-name=4M6J_p-relax
#SBATCH --time=00:20:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab

# launching parsing script 
python /relax_parse_results.py True False - /relax/input /relax/run /relax/output /ddG/input score_bn15_calibrated