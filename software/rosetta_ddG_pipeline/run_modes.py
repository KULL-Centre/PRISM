import subprocess
import os
def relaxation(relax_input,relax_output):

    relax_call = subprocess.Popen(f'sbatch {os.path.join(relax_input,"rosetta_relax.sbatch")}', stdout=subprocess.PIPE, shell=True, cwd=relax_output)
    
    relax_process_id_info = relax_call.communicate()
    relax_process_id = str(relax_process_id_info[0]).split()[3][0:-3]
    parse_relaxation_call = subprocess.Popen(f'sbatch --dependency=afterany:{relax_process_id} {os.path.join(relax_input,"parse_relax.sbatch")}', stdout=subprocess.PIPE, shell=True, cwd=relax_output)
    
    parse_relax_process_id_info = parse_relaxation_call.communicate()
    parse_relax_process_id = str(parse_relax_process_id_info[0]).split()[3][0:-3]

    return parse_relax_process_id



def ddg_calculation(ddG_input,ddG_output,parse_relax_process_id=None):
    dependency = '--dependency=afterany:'
    if parse_relax_process_id == None:
        dependency = ''
    
    cart_ddg_call = subprocess.Popen(f'sbatch {dependency}{parse_relax_process_id} {os.path.join(ddG_input,"rosetta_cartesian_saturation_mutagenesis.sbatch")}', stdout=subprocess.PIPE, shell=True, cwd=ddG_output)

    cart_ddg_process_id_info = cart_ddg_call.communicate()
    print(cart_ddg_process_id_info)      
    cart_ddg_process_id = str(cart_ddg_process_id_info[0]).split()[3][0:-3]

    parse_results_call = subprocess.Popen(f'sbatch --dependency=afterany:{cart_ddg_process_id} {os.path.join(ddG_input,"parse_ddgs.sbatch")}', stdout=subprocess.PIPE, shell=True, cwd=ddG_output) 

    return


