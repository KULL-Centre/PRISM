#! python3
import sys
import subprocess
import os

import rosetta_paths
import pdb_to_fasta_seq



class structure:

    def __init__(self, name, output_path=rosetta_paths.default_output_path):

        self.name = name

        if output_path[-1] == '/':
            self.out_path = output_path[:-1]
        else:
            self.out_path = output_path
            
        if not os.path.isdir('{}/{}/'.format(self.out_path, self.name)):
            os.mkdir('{}/{}/'.format(self.out_path, self.name))
     
        if not os.path.isdir('{}/{}/cleaned_structures'.format(self.out_path, self.name)):
            os.mkdir('{}/{}/cleaned_structures'.format(self.out_path, self.name))

        if not os.path.isdir('{}/{}/rosetta_runs'.format(self.out_path, self.name)):
            os.mkdir('{}/{}/rosetta_runs'.format(self.out_path, self.name))

        if not os.path.isdir('{}/{}/predictions'.format(self.out_path, self.name)):
            os.mkdir('{}/{}/predictions'.format(self.out_path, self.name))
            
        if not os.path.isdir('{}/{}/inputs'.format(self.out_path, self.name)):
            os.mkdir('{}/{}/inputs'.format(self.out_path, self.name))

        print(rosetta_paths.path_to_rosetta)

################################################################################################
#                                   Cleaning and isolating pdb
################################################################################################
        
    def clean_up_and_isolate(self, path_to_pdb, chains):

        path_to_clean_pdb = rosetta_paths.path_to_clean_pdb

        shell_command = 'python2 {} {} {}'.format(path_to_clean_pdb, path_to_pdb, chains)
        print('here is some output from the clean_pdb.py script')
        subprocess.call(shell_command, cwd='{}/{}/cleaned_structures/'.format(self.out_path, self.name), shell=True)
        print('end of output from clean_pdb.py')
        
        self.path_to_cleaned_pdb = '{}/{}/cleaned_structures/{}_{}.pdb'.format(self.out_path, self.name, self.sys_name, chains)
        path_to_cleaned_fasta = '{}/{}/cleaned_structures/{}_{}.fasta'.format(self.out_path, self.name, self.sys_name, chains)
        fasta_file = open(path_to_cleaned_fasta, 'r')
        fasta_lines = fasta_file.readlines()
        fasta_file.close()
        self.fasta_seq = ''

        for line in fasta_lines[1:]:
            self.fasta_seq = self.fasta_seq + line.strip()

        return(self.path_to_cleaned_pdb)

################################################################################################
#                                     Making mutfiles
################################################################################################    

    def make_mutfiles(self,mutation_input):
        check2 = False
        resids = []
        self.path_to_run_folder = '{}/{}/rosetta_runs/{}_{}'.format(self.out_path, self.name, self.sys_name, self.chain_id)

        if not os.path.isdir(self.path_to_run_folder):
            os.mkdir(self.path_to_run_folder)
        path_to_mutfiles = '{}/mutfiles/'.format(self.path_to_run_folder)
        if not os.path.isdir(path_to_mutfiles):
            os.mkdir(path_to_mutfiles)
       
        if mutation_input != None and mutation_input != "proceed":
            mutate = {}
            print("Printing point mutations [WT, ResID, Mutations]")
            with open(mutation_input) as f:
                for line in f:
                    (wt,key, val) = line.split()
                    mutate[int(key)] = wt, val
                    c= mutate[int(key)][0] in list(mutate[int(key)][1])
                    if c == True:
                        mutate[int(key)] = wt, val      
                    else: 
                        val = wt + val
                        mutate[int(key)] = wt, val       
            
            for residue_number in mutate:
                check = self.fasta_seq[residue_number-1] in list(mutate[residue_number][0])
                resids += [residue_number]
                if check == False:
                    check2=True
                final_list = [] 
                for num in mutate[residue_number][1]: 
                    if num not in final_list: 
                        final_list.append(num) 
                print(mutate[residue_number][0],residue_number,''.join(final_list))
                
                mutfile = open(path_to_mutfiles+'mutfile{:0>5}'.format(str(residue_number)), 'w')
                mutfile.write('total '+ str(len(final_list)))
                for AAtype in final_list:
                    mutfile.write('\n1\n')
                    mutfile.write(self.fasta_seq[residue_number-1] + ' ' + str(residue_number) + ' ' + AAtype )
                mutfile.close()
                
        if mutation_input == None:            
            for residue_number in range(1, len(self.fasta_seq)+1):
                mutfile = open(path_to_mutfiles+'mutfile{:0>5}'.format(str(residue_number)), 'w')
                mutfile.write('total 20')
                resids += [residue_number]
                # and then a line for each type of AA
                for AAtype in 'ACDEFGHIKLMNPQRSTVWY':
                    mutfile.write('\n1\n')
                    mutfile.write(self.fasta_seq[residue_number-1] + ' ' + str(residue_number) + ' ' + AAtype )
                mutfile.close()
        return(check2,resids)            

################################################################################################
#                                     Creating sbatch relax
################################################################################################

    def rosetta_sbatch_relax(self, structure_path='defaults to self.path_to_cleaned',relaxfile=None):
        structure_path = self.path_to_cleaned_pdb

        self.path_to_run_folder = '{}/{}/rosetta_runs/{}_{}'.format(self.out_path, self.name, self.sys_name, self.chain_id)
        if not os.path.isdir(self.path_to_run_folder):
            os.makedirs(self.path_to_run_folder)

        path_to_sbatch = '{}/rosetta_relax.sbatch'.format(self.path_to_run_folder)
        if relaxfile==None:
            path_to_relaxflags = rosetta_paths.path_to_parameters + '/relax_flagfile'
        else: 
            path_to_relaxflags = relaxfile
            
        sbatch = open(path_to_sbatch, 'w')
        sbatch.write('''#!/bin/sh
#SBATCH --job-name=pre_relax_rosetta
#SBATCH --time=10:00:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab

# launching rosetta relax
{}bin/relax.linuxgccrelease -s {} -relax:script {}/cart2.script @{}'''.format(rosetta_paths.path_to_rosetta, structure_path, rosetta_paths.path_to_parameters, path_to_relaxflags))
        sbatch.close()
        print(path_to_sbatch)
        return(path_to_sbatch)

################################################################################################
#                                     Creating sbatch parse relax
################################################################################################    
    
    def parse_relax_sbatch(self, path_to_scorefile, path_to_run_folder):
        path_to_parse_relax_script = rosetta_paths.path_to_stability_pipeline + '/relax_parse_results.py'

        path_to_sbatch = '{}/parse_relax.sbatch'.format(self.path_to_run_folder)
        sbatch = open(path_to_sbatch, 'w')
        sbatch.write('''#!/bin/sh
#SBATCH --job-name=parse_relax_rosetta
#SBATCH --time=00:10:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab

# launching parsing script
python {} {} {}'''.format(path_to_parse_relax_script, path_to_scorefile, path_to_run_folder))
        sbatch.close()
        print(path_to_sbatch)
        return(path_to_sbatch)

################################################################################################
#                                     Creating sbatch Rosetta cartesian
################################################################################################    
    
    def write_rosetta_cartesian_sbatch(self,cartesianfile=None,resids=None):
        path_to_sbatch = '{}/rosetta_cartesian_saturation_mutagenesis.sbatch'.format(self.path_to_run_folder)
        resids=resids
        if cartesianfile==None:
            path_to_cartesianflags = rosetta_paths.path_to_parameters + '/cartesian_ddg_flagfile'
        else: 
            path_to_cartesianflags = cartesianfile
        muts=os.listdir(self.path_to_run_folder + '/mutfiles/')
        
        sbatch = open(path_to_sbatch, 'w')
        sbatch.write('''#!/bin/sh
#SBATCH --job-name=Rosetta_cartesian_ddg
#SBATCH --array=0-{}
#SBATCH --time=24:00:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab
#SBATCH --nice
LST=(`ls mutfiles/mutfile*`)
OFFSET=0
INDEX=$((OFFSET+SLURM_ARRAY_TASK_ID))
echo $INDEX

# launching rosetta
{}/bin/cartesian_ddg.linuxgccrelease -s {} -ddg:mut_file ${{LST[$INDEX]}} -out:prefix ddg-$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID @{}'''.format(len(muts), rosetta_paths.path_to_rosetta, '*_bn15_calibrated*.pdb', path_to_cartesianflags))
        sbatch.close()
        print(path_to_sbatch)
        return path_to_sbatch

################################################################################################
#                                     Creating sbatch parse ddg
################################################################################################    
    
    
    def write_parse_ddg_sbatch(self):
        score_sbatch_path = '{}/parse_ddgs.sbatch'.format(self.path_to_run_folder)
        score_sbatch = open(score_sbatch_path, 'w')
        score_sbatch.write('''#!/bin/sh
#SBATCH --job-name=collect_rosetta_ddgs
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --time=0:20:00
#SBATCH --partition=sbinlab

# This sbatch script launches the parse parse_rosetta_ddgs function, from the parse_cartesian_ddgs
# it will output a file in the prediction_files/ folder.
python3 {}/parse_rosetta_ddgs.py {} {} {} {}'''.format(rosetta_paths.path_to_stability_pipeline, self.sys_name, self.chain_id, self.fasta_seq, self.out_path))
        score_sbatch.close()
        return(score_sbatch_path)


