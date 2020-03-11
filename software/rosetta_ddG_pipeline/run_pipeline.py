from structure_input import structure
from pdb_to_fasta_seq import pdb_to_fasta_seq
import subprocess
import sys, getopt
import os
import rosetta_paths
from checks import compare_mutfile


#Hello
def predict_stability(argv):

	os.chdir(os.getcwd())
	print(os.getcwd())
	relax_flag_file=rosetta_paths.path_to_stability_pipeline+'/rosetta_parameters/relax_flagfile'
	ddg_flag_file=rosetta_paths.path_to_stability_pipeline+'/rosetta_parameters/cartesian_ddg_flagfile'  
    
	inputfile=''; indexfile=''; uniprot_accesion=''; mutation_input=None; chain_id='nochain'; relaxfile=relax_flag_file; cartesianfile=ddg_flag_file; mode="create"

	try:
		opts, args = getopt.getopt(argv,"hr:o:n:u:m:s:c:i:p:",["rfile=","ofile=","nfile=","ufile=","mfile=","sfile=","cfile=","mode="])
	except getopt.GetoptError:
		print (' -o output_path \n -s Structure.pdb \n -m Mutation_Input (optional) \n -r Relax_flags (optional) \n -c cartesian_ddg_flags (optional) \n -i modes [print,create,proceed,fullrun] \n        print: prints default flag files \n        create: Creates all run files \n        proceed: Starts calculations with created run files \n        fullrun: runs full pipeline')
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print (' -o output_path \n -s Structure.pdb \n -m Mutation_Input (optional) \n -r Relax_flags (optional) \n -c cartesian_ddg_flags (optional) \n -i modes [print,create,proceed,fullrun] \n        print: prints default flag files \n        create: Creates all run files \n        proceed: Starts calculations with created run files \n        fullrun: runs full pipeline')
			sys.exit()
		if opt in ("-r", "--rfile"):
			relaxfile = arg
			print ('Relax flags are ', relaxfile)
		if opt in ("-o", "--ofile"):
			out_path = arg
			print ('Output folder is ', out_path)
		if opt in ("-n", "--nfile"):
			indexfile = arg
			print ('Index file is ', indexfile)
		if opt in ("-u", "--ufile"):
			uniprot_accesion = arg
			print ('Uniprot file is ', uniprot_accesion)
		if opt in ("-m", "--mfile"):
			mutation_input = arg
			print ('Mutational file is ', mutation_input)  
		if opt in ("-s", "--sfile"):
			structure_list = arg
			print ('structure file is ', structure_list)
		if opt in ("-c", "--cfile"):
			cartesianfile = arg
			print ('Cartesian flags are ', cartesianfile)
		if opt in ("-i", "--mode"):
			mode = arg
			print('The selected mode is ', mode)            

            


########################################################################################################        
#                                        PRINT AND PROCEED MODE        
########################################################################################################
    
	if mode == "print":
		open("relax_flag_file_copy", "w").writelines(open(relaxfile).readlines())     
		open("ddg_flag_file_copy", "w").writelines(open(cartesianfile).readlines())        
	if mode=="proceed":
		mutation_input == "proceed"
        
########################################################################################################        
#                                         HANDLING THE STRUCTURES       
########################################################################################################
        
	name = structure_list.split('/')[-1].split('.')[0] 
	structure_instance = structure(name, output_path=out_path)    
	structure_instance.sys_name = name
	structure_instance.chain_id = chain_id
	structure_instance.path = structure_list
	structure_instance.clean_up_and_isolate(structure_instance.path, structure_instance.chain_id)
	structure_instance.fasta_seq = pdb_to_fasta_seq(structure_instance.path_to_cleaned_pdb)        
	check2, resids = structure_instance.make_mutfiles(mutation_input) 
	print(resids)
    
########################################################################################################        
#                                 CREATING SBATCH AND INPUT FILES        
########################################################################################################    
   

        
	if mode == "create" or mode == "fullrun":
		path_to_relax_sbatch = structure_instance.rosetta_sbatch_relax(relaxfile=relaxfile)
		path_to_parse_relax_results_sbatch = structure_instance.parse_relax_sbatch(structure_instance.path_to_run_folder + '/score_bn15_calibrated.sc', structure_instance.path_to_run_folder)
		path_to_ddg_calc_sbatch = structure_instance.write_rosetta_cartesian_sbatch(cartesianfile=cartesianfile,resids=resids)
		path_to_parse_ddgs_sbatch = structure_instance.write_parse_ddg_sbatch()
        
		print("Creating all required files \n Inspect before continuing the calculation \n continue calculation by using same command and -i proceed")
 

   
	path_to_input = '{}/{}/{}/'.format(out_path,name,"inputs")
	with open(path_to_input+'inputs','w') as file:
		file.write('NAME ='+name +'\n' + 'CHAIN =' +chain_id +'\n' +'STRUCTURE ='+structure_list +'\n'+inputfile+'\n'+'OUTPATH ='+out_path+'\n'+indexfile+'\n'+uniprot_accesion+'\n')
	file.close()
	open(path_to_input+"relax_flag_file", "w").writelines(open(relaxfile).readlines())     
	open(path_to_input+"ddg_flag_file", "w").writelines(open(cartesianfile).readlines()) 
	print("Input file can be found at",path_to_input) 
    
########################################################################################################        
#                                           Check phase 1        
########################################################################################################    
    
	check1 = compare_mutfile(structure_instance.fasta_seq,structure_instance.path_to_run_folder,mutation_input)        

    
	if check1 == True or check2 == True:
		print("ERROR: STOPPING SCRIPT DUE TO RESIDUE MISMATCH BETWEEN MUTFILE AND PDB SEQUENCE")
		sys.exit()

########################################################################################################        
#                                 THIS PART STARTS THE CALCULATIONS        
########################################################################################################
        
	if mode == "fullrun" or mode == "proceed":
        
		print('submitting sbatch jobs')
		relax_call = subprocess.Popen('sbatch rosetta_relax.sbatch', stdout=subprocess.PIPE, shell=True, cwd=structure_instance.path_to_run_folder)
 
		relax_process_id_info = relax_call.communicate()
		relax_process_id = str(relax_process_id_info[0]).split()[3][0:-3]
		parse_relaxation_call = subprocess.Popen('sbatch --dependency=afterany:{} parse_relax.sbatch'.format(relax_process_id), stdout=subprocess.PIPE, shell=True, cwd=structure_instance.path_to_run_folder)
 

		parse_relax_process_id_info = parse_relaxation_call.communicate()
		parse_relax_process_id = str(parse_relax_process_id_info[0]).split()[3][0:-3]
       

		cart_ddg_call = subprocess.Popen('sbatch --dependency=afterany:{} rosetta_cartesian_saturation_mutagenesis.sbatch'.format(parse_relax_process_id), stdout=subprocess.PIPE, shell=True, cwd=structure_instance.path_to_run_folder)
#		cart_ddg_call = subprocess.Popen('sbatch rosetta_cartesian_saturation_mutagenesis.sbatch'.format(parse_relax_process_id), stdout=subprocess.PIPE, shell=True, cwd=structure_instance.path_to_run_folder)

		cart_ddg_process_id_info = cart_ddg_call.communicate()
		print(cart_ddg_process_id_info)      
		cart_ddg_process_id = str(cart_ddg_process_id_info[0]).split()[3][0:-3]

		parse_results_call = subprocess.Popen('sbatch --dependency=afterany:{} parse_ddgs.sbatch'.format(cart_ddg_process_id), stdout=subprocess.PIPE, shell=True, cwd=structure_instance.path_to_run_folder) 

        
########################################################################################################        
#                                          Initiate scripts    
########################################################################################################
if __name__ == '__main__':
    predict_stability(sys.argv[1:])

