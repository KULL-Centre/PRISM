#!/bin/sh 
#SBATCH --job-name=4M6J_p-ddg
#SBATCH --array=1 
#SBATCH --nodes=1 
#SBATCH --time=0:10:00 
#SBATCH --partition=sbinlab 

#This sbatch script launches the parse parse_rosetta_ddgs function, from the parse_cartesian_ddgs 
python3 /parser_ddg_v2.py 4M6J A SLNCIVAVSQNMGIGKNGDLPWPPLRNEFRYFQRMTTTSSVEGKQNLVIMGKKTWFSIPEKNRPLKGRINLVLSRELKEPPQGAHFLSRSLDDALKLTEQPELANKVDMVWIVGGSSVYKEAMNHPGHLKLFVTRIMQDFESDTFFPEIDLEKYKLLPEYPGVLSDVQEEKGIKYKFEVYEKND /ddG/run /ddG/output /prepare/checking/structure_input.json /ddG/input /output /prepare/checking False True 351877c0dacd5079cbde563e4c1e1c93f621a253tagv0.2.2-9-g351877c /prepare/cleaning False 2.9