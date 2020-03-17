import os
import re
import sys, getopt
import shutil

from structure_input import structure
from pdb_to_fasta_seq import pdb_to_fasta_seq
import rosetta_paths
from checks import compare_mutfile
from Args_pipeline import parse_args2
import run_modes
import storeinputs 


def predict_stability(argv):

    os.chdir(os.getcwd())
    print(os.getcwd())
    
    args = parse_args2()
    
    structure_list = args.STRUC_FILE; uniprot_accesion = args.UNIPROT_FILE;outpath=args.OUTPUT_FILE; mutation_input= args.MUTATION_INPUT; relaxfile= args.RELAX_FLAG_FILE; cartesianfile= args.CART_FLAG_FILE; mode= args.MODE; chain_id=args.CHAIN
    

    ## System name
    name = structure_list.split('/')[-1].split('.')[0] 
    path_to_input = '{}/{}/{}/'.format(outpath,name,"inputs")
        
########################################################################################################        
#                                        PRINT AND PROCEED MODE        
########################################################################################################
    
    if mode == "print":
        open("relax_flag_file_copy", "w").writelines(open(relax_flag_file).readlines())     
        open("ddg_flag_file_copy", "w").writelines(open(ddg_flag_file).readlines())        
    if mode=="proceed" or mode =="relax" or mode == "ddg_calculation":
        mutation_input == "proceed"
    
########################################################################################################        
#                                         HANDLING THE STRUCTURES       
########################################################################################################
        
    ## Defining structure parameters
    structure_instance = structure(name, output_path=outpath)    
    structure_instance.sys_name = name; structure_instance.chain_id = chain_id; structure_instance.path = structure_list;
    
    ## Cleaning pdb and making fasta based on pdb
    structure_instance.clean_up_and_isolate(structure_instance.path, structure_instance.chain_id)
    structure_instance.fasta_seq = pdb_to_fasta_seq(structure_instance.path_to_cleaned_pdb) 
    

    
    
    
    
    if mode == "create" or mode == "fullrun":
        ## Making mutfiles and checks
        check2, resids = structure_instance.make_mutfiles(mutation_input) 
         
        
        path_to_relax_sbatch = structure_instance.rosetta_sbatch_relax(relaxfile=relaxfile)
        path_to_parse_relax_results_sbatch = structure_instance.parse_relax_sbatch(structure_instance.path_to_run_folder + '/score_bn15_calibrated.sc', structure_instance.path_to_run_folder)
        path_to_ddg_calc_sbatch = structure_instance.write_rosetta_cartesian_sbatch(cartesianfile=cartesianfile,resids=resids)
        path_to_parse_ddgs_sbatch = structure_instance.write_parse_ddg_sbatch()
        if mode == "create":    
            print("Creating all required files \n Inspect before continuing the calculation \n continue calculation by using same command and -i proceed")
     
    
    storeinputs.storeinputfuc(name,chain_id,structure_list,outpath,uniprot_accesion,relaxfile,cartesianfile,path_to_input)   
    check1 = compare_mutfile(structure_instance.fasta_seq,structure_instance.path_to_run_folder,mutation_input)  

    if check1 == True or check2 == True:
        print("ERROR: STOPPING SCRIPT DUE TO RESIDUE MISMATCH BETWEEN MUTFILE AND PDB SEQUENCE")
        sys.exit()
    
    if mode == "relax":
        parse_relax_process_id = run_modes.relaxation(structure_instance.path_to_run_folder)
        
    
    if mode == "fullrun":
        parse_relax_process_id = run_modes.relaxation(structure_instance.path_to_run_folder)
        run_modes.ddg_calculation(structure_instance.path_to_run_folder,parse_relax_process_id)
    
    if mode == "ddg_calculation":
        run_modes.ddg_calculation(structure_instance.path_to_run_folder)
    
    
    if mode == "proceed":
        parse_relax_process_id = run_modes.relaxation(structure_instance.path_to_run_folder)
        run_modes.ddg_calculation(structure_instance.path_to_run_folder,parse_relax_process_id) 
    

########################################################################################################        
#                                          Initiate scripts    
########################################################################################################
if __name__ == '__main__':
    predict_stability(sys.argv[1:])    