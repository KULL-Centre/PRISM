import os
import subprocess
import time

# yes, I know global vars are bad...
default_path = {
    'muscle_exec': '/groups/sbinlab/software/muscle/muscle3.8.31_i86linux64',
    'TMalign_exec': '/groups/sbinlab/software/TMalign/TMalign',
    'MMLigner_exec': '/groups/sbinlab/software/MMLigner/mmligner',
    # 'ddG_pipeline': '/groups/sbinlab/software/PRISM_tools/rosetta_stability-dev/software/rosetta_ddG_pipeline',
    'ddG_pipeline': '/lustre/hpc/sbinlab/panf/PRISM/software/rosetta_ddG_pipeline/',
    # 'Rosetta_main_path': '/sbinlab/software/Rosetta_2021_Aug_c7009b3/source/',
    'Rosetta_main_path': '/sbinlab/software/Rosetta_2025_Jan_e5e4b27/source/',
    # 'Rosetta_tools_path': '/sbinlab/software/Rosetta_tools/tools/',
    'Rosetta_tools_path': '/sbinlab/software/Rosetta_2025_Jan_e5e4b27/tools/',
    # 'Rosetta_database_path': '/sbinlab/software/Rosetta_2021_Aug_c7009b3/database/',
    'Rosetta_database_path': '/sbinlab/software/Rosetta_2025_Jan_e5e4b27/database/',
    'Rosetta_extension': 'linuxgccrelease',
#    'prism_parser': '/groups/sbinlab/software/PRISM_tools/prism_parser/scripts',
    'prism_parser': '/groups/sbinlab/software/PRISM_tools/prism_parser/scripts',
    'BIOLIB_TOKEN': 'wEogNuisIGYhUmKe8T3KBaPnWYb6vtTRElSSZRFmRkE'
}


def load_env(env):
    # load exec if within environment
    if os.getenv(env) != None:
        return os.getenv(env)
    else:
        os.environ[env] = default_path[env]
    return default_path[env]


muscle_exec = load_env('muscle_exec')
TMalign_exec = load_env('TMalign_exec')
MMLigner_exec = load_env('MMLigner_exec')
ddG_pipeline = load_env('ddG_pipeline')
Rosetta_main_path = load_env('Rosetta_main_path')
Rosetta_tools_path = load_env('Rosetta_tools_path')
Rosetta_database_path = load_env('Rosetta_database_path')
Rosetta_extension = load_env('Rosetta_extension')
prism_parser = load_env('prism_parser')
BIOLIB_TOKEN = load_env('BIOLIB_TOKEN')

try:
    pipes2 = subprocess.Popen("git describe --tags", shell=True, cwd=ddG_pipeline, stdout=subprocess.PIPE,stderr=subprocess.PIPE,)
    std_out2, std_err = pipes2.communicate()
    pipes = subprocess.Popen("git log -n 1 | grep -i commit", shell=True, cwd=ddG_pipeline, stdout=subprocess.PIPE,stderr=subprocess.PIPE,)
    std_out, std_err = pipes.communicate()
    time.sleep(3)
    sha = std_out.strip().decode('UTF-8').split()[1]
    tag = std_out2.strip().decode('UTF-8').split()[0]
except:
    sha = 'dev'
    tag = 'XXXXX'


print('current env paths & exec:', TMalign_exec, muscle_exec, ddG_pipeline, Rosetta_main_path,
      Rosetta_tools_path, Rosetta_database_path, Rosetta_extension, prism_parser, sha, tag, BIOLIB_TOKEN)

# Rosetta paths
path_to_rosetta = Rosetta_main_path
path_to_clean_pdb = os.path.join(
    Rosetta_tools_path, 'protein_tools', 'scripts', 'clean_pdb.py')
#path_to_clean_pdb = 'Cartesian_pipeline_runs/448_1EAW_new/clean_pdb.py'
path_to_clean_keep_ligand = os.path.join(
    path_to_rosetta, 'src', 'apps', 'public', 'relax_w_allatom_cst', 'clean_pdb_keep_ligand.py')
# Personal paths
path_to_stability_pipeline = ddG_pipeline
default_output_path = os.path.join(path_to_stability_pipeline, 'output')
path_to_data = os.path.join(
    path_to_stability_pipeline, 'data')

# Muscle exec
path_to_muscle = muscle_exec

# TMalign exec
path_to_TMalign = TMalign_exec
