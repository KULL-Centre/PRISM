#!/bin/sh 
#SBATCH --job-name=4M6J_p-ddg
#SBATCH --array=1 
#SBATCH --nodes=1 
#SBATCH --time=10:00:00 
#SBATCH --partition=sbinlab 

#This sbatch script launches the parse parse_rosetta_ddgs function, from the parse_cartesian_ddgs 
python3 /parser_ddg_v2.py 4M6J A A SLNCIVAVSQNMGIGKNGDLPWPPLRNEFRYFQRMTTTSSVEGKQNLVIMGKKTWFSIPEKNRPLKGRINLVLSRELKEPPQGAHFLSRSLDDALKLTEQPELANKVDMVWIVGGSSVYKEAMNHPGHLKLFVTRIMQDFESDTFFPEIDLEKYKLLPEYPGVLSDVQEEKGIKYKFEVYEKND /ddG/run /ddG/output /prepare/checking/structure_input.json /ddG/input /output /prepare/checking False True XXXtagvXXX /prepare/cleaning False 2.9