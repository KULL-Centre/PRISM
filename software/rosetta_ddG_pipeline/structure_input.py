#! python3
import sys
import subprocess
import os
import numpy as np
import rosetta_paths
import pdb_to_fasta_seq



class structure:



################################################################################################
#                                   Cleaning and isolating pdb
################################################################################################
        
    def clean_up_and_isolate(self,input_cleaning, path_to_pdb, chains,name):

        path_to_clean_pdb = rosetta_paths.path_to_clean_pdb

        shell_command = 'python2 {} {} {}'.format(path_to_clean_pdb, path_to_pdb, chains)
        print('here is some output from the clean_pdb.py script')
        subprocess.call(shell_command, cwd='{}'.format(input_cleaning), shell=True)
        print('end of output from clean_pdb.py')
        
        self.path_to_cleaned_pdb = '{}{}_{}.pdb'.format(input_cleaning,self.sys_name, chains)
        self.path_to_cleaned_fasta = '{}{}_{}.fasta'.format(input_cleaning,self.sys_name, chains)
        fasta_file = open(self.path_to_cleaned_fasta, 'r')
        fasta_lines = fasta_file.readlines()
        fasta_file.close()
        self.fasta_seq = ''

        for line in fasta_lines[1:]:
            self.fasta_seq = self.fasta_seq + line.strip()

        return(self.path_to_cleaned_pdb)
    
    
    def read_fasta(self, uniprot_accesion):


        fastau_file = open(uniprot_accesion, 'r')
        fastau_lines = fastau_file.readlines()
        fastau_file.close()
        self.uniprot_seq = ''

        for line in fastau_lines[1:]:
            self.uniprot_seq = self.uniprot_seq + line.strip()

        return(self.uniprot_seq)

################################################################################################
#                                   Alignment to uniprot
################################################################################################    
    
    
    def muscle_align_to_uniprot(self,input_checking,uniprot_sequence):

        path_to_muscle= rosetta_paths.path_to_muscle
        self.path_to_fasta = input_checking+'fasta_file.fasta'
        self.path_to_alignment = input_checking+'alignment.txt'
        with open(self.path_to_fasta, 'w') as fasta_file:
            fasta_file.write('>{}_structure_sequence\n'.format(self.sys_name))
            fasta_file.write('{}\n'.format(self.fasta_seq))
            fasta_file.write('>{}_uniprot_sequence\n'.format(self.sys_name))
            fasta_file.write('{}\n'.format(uniprot_sequence))
   
        shell_call = '{} -in {} -out {}'.format(path_to_muscle, self.path_to_fasta, self.path_to_alignment)
        subprocess.call(shell_call, shell=True)

        alignment_sequences = {}
        current_seq = ''    
        current_name = 'first'
        with open(self.path_to_alignment, 'r') as alignment_file:
            for line in alignment_file.readlines():
                if line[0] == '>':
                    if current_name != 'first':
                        alignment_sequences[current_name] = current_seq
    
                    current_name = line[1:].strip()
                    current_seq = ''
                else:
                    current_seq = current_seq + line.strip()
            alignment_sequences[current_name] = current_seq
    
        alignment_index_numbers = []
        current_index = 0
        for key in alignment_sequences:
            if 'uniprot_sequence' in key:
                for letter in alignment_sequences[key]:
                    current_index = current_index + 1
                    if letter != '-':
                        alignment_index_numbers.append(current_index)
    

        current_index = 0
    
        structure_index_numbers = []
        for key in alignment_sequences:
            if 'uniprot_sequence' not in key:
                for letter in alignment_sequences[key]:
                    if letter != '-':
                        structure_index_numbers.append(str(alignment_index_numbers[current_index]))
                    current_index = current_index + 1
    
        self.structure_index_numbers = structure_index_numbers

        self.path_to_index_string = input_checking + 'uniprot_index_list.txt'
        with open(self.path_to_index_string, 'w') as index_file:
            index_list_as_string = '\n'.join(structure_index_numbers)
            index_file.write(str(index_list_as_string))
    
        path_to_index_string = self.path_to_index_string
        return(path_to_index_string)    
    
################################################################################################
#                                     Making mutfiles
################################################################################################    

    def make_mutfiles(self,mutation_input,input_mutfiles):
        check2 = False
        resids = []
        
        path_to_alignment = self.path_to_index_string
        alignment = np.loadtxt(path_to_alignment)

        path_to_mutfiles = input_mutfiles

        alignment_dic={}
        for n in enumerate(alignment):
            alignment_dic[int(n[1])] = int(n[0]+1)

        
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
                residue_number_ros = alignment_dic[residue_number]
                print(residue_number_ros,residue_number)

                check = self.fasta_seq[residue_number_ros] in list(mutate[residue_number][0])
                resids += [residue_number]
                
                
                if check == False:
                    check2=True
                final_list = [] 
                for num in mutate[residue_number][1]: 
                    if num not in final_list: 
                        final_list.append(num) 
                print(mutate[residue_number][0],residue_number,''.join(final_list))
                
                mutfile = open(path_to_mutfiles+'mutfile{:0>5}'.format(str(residue_number_ros)), 'w')
                mutfile.write('total '+ str(len(final_list)))
                for AAtype in final_list:
                    mutfile.write('\n1\n')
                    mutfile.write(self.fasta_seq[residue_number_ros-1] + ' ' + str(residue_number_ros) + ' ' + AAtype )
                mutfile.close()
                
        if mutation_input == None:            
            for residue_number in alignment:
                residue_number_ros = alignment_dic[residue_number]
                mutfile = open(path_to_mutfiles+'mutfile{:0>5}'.format(str(residue_number_ros)), 'w')
                mutfile.write('total 20')
                resids += [residue_number]
                # and then a line for each type of AA
                for AAtype in 'ACDEFGHIKLMNPQRSTVWY':
                    mutfile.write('\n1\n')
                    mutfile.write(self.fasta_seq[residue_number_ros-1] + ' ' + str(residue_number_ros) + ' ' + AAtype )
                mutfile.close()
        return(check2,resids)            

################################################################################################
#                                     Creating sbatch relax
################################################################################################

    def rosetta_sbatch_relax(self,relax_input, structure_path='defaults to self.path_to_cleaned',relaxfile=None):
        structure_path = self.path_to_cleaned_pdb


        path_to_sbatch = relax_input+'rosetta_relax.sbatch'
        if relaxfile==None:
            path_to_relaxflags = rosetta_paths.path_to_parameters + '/relax_flagfile'
        else: 
            path_to_relaxflags = relaxfile
            
        sbatch = open(path_to_sbatch, 'w')
        sbatch.write('''#!/bin/sh
#SBATCH --job-name=relax_{}
#SBATCH --time=10:00:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab

# launching rosetta relax
{}bin/relax.linuxgccrelease -s {} -relax:script {}/cart2.script @{}'''.format(self.sys_name,rosetta_paths.path_to_rosetta, structure_path, rosetta_paths.path_to_parameters, path_to_relaxflags))
        sbatch.close()
        print(path_to_sbatch)
        return(path_to_sbatch)

################################################################################################
#                                     Creating sbatch parse relax
################################################################################################    
    
    def parse_relax_sbatch(self, relax_input,relax_output):
        path_to_parse_relax_script = rosetta_paths.path_to_stability_pipeline + '/relax_parse_results.py'
        path_to_scorefile= relax_output + 'score_bn15_calibrated.sc'
        path_to_sbatch = '{}/parse_relax.sbatch'.format(relax_input)
        sbatch = open(path_to_sbatch, 'w')
        sbatch.write('''#!/bin/sh
#SBATCH --job-name=parse_relax_rosetta
#SBATCH --time=00:20:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab

# launching parsing script
python {} {} {}'''.format(path_to_parse_relax_script, path_to_scorefile, relax_output))
        sbatch.close()
        print(path_to_sbatch)
        return(path_to_sbatch)

################################################################################################
#                                     Creating sbatch Rosetta cartesian
################################################################################################    
    
    def write_rosetta_cartesian_sbatch(self,relax_input,ddG_input,input_mutfiles,cartesianfile=None,resids=None):
        path_to_sbatch = ddG_input + 'rosetta_cartesian_saturation_mutagenesis.sbatch'
        resids=resids
        if cartesianfile==None:
            path_to_cartesianflags = rosetta_paths.path_to_parameters + '/cartesian_ddg_flagfile'
        else: 
            path_to_cartesianflags = cartesianfile
        muts=os.listdir(input_mutfiles)
        
        sbatch = open(path_to_sbatch, 'w')
        sbatch.write('''#!/bin/sh
#SBATCH --job-name=cartesian_{}
#SBATCH --array=0-{}
#SBATCH --time=32:00:00
#SBATCH --mem 5000
#SBATCH --partition=sbinlab
#SBATCH --nice
LST=(`ls {}mutfile*`)
OFFSET=0
INDEX=$((OFFSET+SLURM_ARRAY_TASK_ID))
echo $INDEX

# launching rosetta
{}/bin/cartesian_ddg.linuxgccrelease -s {} -ddg:mut_file ${{LST[$INDEX]}} -out:prefix ddg-$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID @{}'''.format(self.sys_name,len(muts),input_mutfiles, rosetta_paths.path_to_rosetta, relax_input+'*_bn15_calibrated*.pdb', path_to_cartesianflags))
        sbatch.close()
        print(path_to_sbatch)
        return path_to_sbatch

################################################################################################
#                                     Creating sbatch parse ddg
################################################################################################    
    
    
    def write_parse_ddg_sbatch(self,ddG_input,ddG_output):
        score_sbatch_path = ddG_input+'parse_ddgs.sbatch'
        score_sbatch = open(score_sbatch_path, 'w')
        score_sbatch.write('''#!/bin/sh
#SBATCH --job-name=collect_rosetta_ddgs
#SBATCH --array=1
#SBATCH --nodes=1
#SBATCH --time=0:10:00
#SBATCH --partition=sbinlab

# This sbatch script launches the parse parse_rosetta_ddgs function, from the parse_cartesian_ddgs
python3 {}/parse_rosetta_ddgs.py {} {} {} {}'''.format(rosetta_paths.path_to_stability_pipeline, self.sys_name, self.chain_id, self.fasta_seq, ddG_input, ddG_output))
        score_sbatch.close()
        return(score_sbatch_path)


