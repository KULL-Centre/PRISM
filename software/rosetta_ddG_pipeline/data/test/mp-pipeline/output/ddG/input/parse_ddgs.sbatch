#!/bin/sh 
#SBATCH --job-name=collect_rosetta_ddgs_PagP 
#SBATCH --array=1 
#SBATCH --nodes=1 
#SBATCH --time=0:10:00 
#SBATCH --partition=sbinlab 

#This sbatch script launches the parse parse_rosetta_ddgs function, from the parse_cartesian_ddgs 
python3 /parse_rosetta_ddgs.py PagP A MNADEWMTTFRENIAQTWQQPEHYDLYIPAITWHARFAYNERPWGGGFGLSRWDEKGNWHGLYAMAFKDSWNKWEPIAGYGWESTWRPLADENFHLGLGFTAGVTARDNWNYIPLPVLLPLASVGYGPVTFQMTYIPGTYNNGNVYFAWMRFQFL /ddG/run /ddG/output