#!/bin/sh 
#SBATCH --job-name=homodimer-1afo_renum_uniquechain_p-ddg
#SBATCH --array=1 
#SBATCH --nodes=1 
#SBATCH --time=0:10:00 
#SBATCH --partition=sbinlab 

#This sbatch script launches the parse parse_rosetta_ddgs function, from the parse_cartesian_ddgs 
python3 /parser_ddg_v2.py homodimer-1afo_renum_uniquechain A VQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKKVQLAHHFSEPEITLIIFGVMAGVIGTILLISYGIRRLIKK /ddG/run /ddG/output /prepare/checking/structure_input.json /ddG/input /output /prepare/checking False True 94714b9891abbb08e5b64c5f67e5d1afad815edbtagv0.1.1-19-g94714b9 /prepare/cleaning True 1