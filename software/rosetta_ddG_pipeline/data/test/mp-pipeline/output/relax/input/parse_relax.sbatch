#!/bin/sh
#SBATCH --job-name=parse_relax_rosetta_relax_scores
#SBATCH --time=00:20:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab

# launching parsing script 
python /relax_parse_results.py /relax/run /relax/output /ddG/input