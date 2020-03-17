import subprocess

def relaxation(runfolder):
    
    relax_call = subprocess.Popen('sbatch rosetta_relax.sbatch', stdout=subprocess.PIPE, shell=True, cwd=runfolder)
 
    relax_process_id_info = relax_call.communicate()
    relax_process_id = str(relax_process_id_info[0]).split()[3][0:-3]
    
    parse_relaxation_call = subprocess.Popen('sbatch --dependency=afterany:{} parse_relax.sbatch'.format(relax_process_id), stdout=subprocess.PIPE, shell=True, cwd=runfolder)
 

    parse_relax_process_id_info = parse_relaxation_call.communicate()
    parse_relax_process_id = str(parse_relax_process_id_info[0]).split()[3][0:-3]

    return parse_relax_process_id



def ddg_calculation(runfolder,parse_relax_process_id=None):
    dependency = '--dependency=afterany:'
    if parse_relax_process_id == None:
        dependency = ''
    
    cart_ddg_call = subprocess.Popen('sbatch {}{} rosetta_cartesian_saturation_mutagenesis.sbatch'.format(dependency,parse_relax_process_id), stdout=subprocess.PIPE, shell=True, cwd=runfolder)

    cart_ddg_process_id_info = cart_ddg_call.communicate()
    print(cart_ddg_process_id_info)      
    cart_ddg_process_id = str(cart_ddg_process_id_info[0]).split()[3][0:-3]

    parse_results_call = subprocess.Popen('sbatch --dependency=afterany:{} parse_ddgs.sbatch'.format(cart_ddg_process_id), stdout=subprocess.PIPE, shell=True, cwd=runfolder) 

    return


