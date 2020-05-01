"""structure_inputs.py includes several preprocessing and script generating functions.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-04-15

"""

# Standard library imports
import logging as logger
import os
import subprocess
import sys

# Third party imports
import numpy as np

# Local application imports
import pdb_to_fasta_seq
import rosetta_paths




class structure:

    ##########################################################################
    #                                   Cleaning and isolating pdb
    ##########################################################################

    def clean_up_and_isolate(self, input_cleaning, path_to_pdb, chains, run_struc, name='input',ligand=None):
        if  ligand == None:
            
            path_to_clean_pdb = rosetta_paths.path_to_clean_pdb
    
            shell_command = f'python2 {path_to_clean_pdb} {path_to_pdb} {run_struc}'
            print('here is some output from the clean_pdb.py script')
            subprocess.call(shell_command, cwd=input_cleaning, shell=True)
            print('end of output from clean_pdb.py')
    
            self.path_to_cleaned_pdb = os.path.join(input_cleaning, f'{name}_{run_struc}.pdb')
            print(str(run_struc))
            if type(run_struc) != list:
                run_struc = [run_struc]
            for chain in run_struc:
                self.path_to_cleaned_fasta = os.path.join(input_cleaning, f'{name}_{chain}.fasta')
                fasta_lines = open(self.path_to_cleaned_fasta, 'r').readlines()
                self.fasta_seq = ''
        
                for line in fasta_lines[1:]:
                    self.fasta_seq = self.fasta_seq + line.strip()
                    
        if ligand == True:
            path_to_clean_pdb = rosetta_paths.path_to_clean_keep_ligand
            
            shell_command = f'python2 {path_to_clean_pdb} {path_to_pdb} {chains}'
            print('here is some output from the clean_pdb_keep_ligand.py script')
            subprocess.call(shell_command, cwd=input_cleaning, shell=True)
            print('end of output from clean_pdb_keep_ligand.py')
            self.path_to_cleaned_pdb = os.path.join(input_cleaning, f'{name}_{run_struc}.pdb{chains}.pdb')
        return(self.path_to_cleaned_pdb)

    def read_fasta(self, uniprot_accesion):

        fastau_file = open(uniprot_accesion, 'r')
        fastau_lines = fastau_file.readlines()
        fastau_file.close()
        self.uniprot_seq = ''

        for line in fastau_lines[1:]:
            self.uniprot_seq = self.uniprot_seq + line.strip()

        return(self.uniprot_seq)

##########################################################################
#                                   Alignment to uniprot
##########################################################################

    def muscle_align_to_uniprot(self, input_checking, uniprot_sequence):

        path_to_muscle = rosetta_paths.path_to_muscle
        self.path_to_fasta = os.path.join(input_checking, 'fasta_file.fasta')
        self.path_to_alignment = os.path.join(input_checking, 'alignment.txt')
        with open(self.path_to_fasta, 'w') as fasta_file:
            fasta_file.write('>{}_structure_sequence\n'.format(self.sys_name))
            fasta_file.write('{}\n'.format(self.fasta_seq))
            fasta_file.write('>{}_uniprot_sequence\n'.format(self.sys_name))
            fasta_file.write('{}\n'.format(uniprot_sequence))

        shell_call = '{} -in {} -out {}'.format(
            path_to_muscle, self.path_to_fasta, self.path_to_alignment)
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
                        structure_index_numbers.append(
                            str(alignment_index_numbers[current_index]))
                    current_index = current_index + 1

        self.structure_index_numbers = structure_index_numbers

        self.path_to_index_string = os.path.join(input_checking, 'uniprot_index_list.txt')
        with open(self.path_to_index_string, 'w') as index_file:
            index_list_as_string = '\n'.join(structure_index_numbers)
            index_file.write(str(index_list_as_string))

        path_to_index_string = self.path_to_index_string
        return(path_to_index_string)

##########################################################################
#                                     Making mutfiles
##########################################################################

    def make_mutfiles(self, mutation_input, path_to_mutfiles,structure_dic,chain_id):
        check2 = False
        resdata = structure_dic["resdata"]
        strucdata = structure_dic["strucdata"]
        
        path_to_alignment = self.path_to_index_string
        alignment = np.loadtxt(path_to_alignment)
        
        alignment_dic = {}
        for n in enumerate(alignment):
            alignment_dic[int(n[1])] = int(n[0] + 1)

        if mutation_input != None and mutation_input != "proceed":
            mutate = {}
            print("Printing point mutations [WT, ResID, Mutations]")
            with open(mutation_input) as f:
                for line in f:
                    (wt, key, val) = line.split()
                    mutate[int(key)] = wt, val
                    c = mutate[int(key)][0] in list(mutate[int(key)][1])
                    if c == True:
                        mutate[int(key)] = wt, val
                    else:
                        val = wt + val
                        mutate[int(key)] = wt, val

            for residue_number in mutate:
                residue_number_ros = alignment_dic[residue_number]
                print(residue_number_ros, residue_number)

                check = self.fasta_seq[residue_number_ros] in list(
                    mutate[residue_number][0])
                

                if check == False:
                    check2 = True
                final_list = []
                for num in mutate[residue_number][1]:
                    if num not in final_list:
                        final_list.append(num)
                print(mutate[residue_number][0],
                      residue_number, ''.join(final_list))

                mutfile = open(os.path.join(path_to_mutfiles, f'mutfile{str(residue_number_ros):0>5}'), 'w')
                
                mutfile.write('total ' + str(len(final_list)))
                for AAtype in final_list:
                    mutfile.write('\n1\n')
                    mutfile.write(self.fasta_seq[
                                  residue_number_ros - 1] + ' ' + str(residue_number_ros) + ' ' + AAtype)
                mutfile.close()

        if mutation_input == None:
            print("Printing mutfiles!")
            for residue_number_ros in resdata:
                if resdata[residue_number_ros][2] == chain_id:
                    mutfile = open(os.path.join(path_to_mutfiles, f'mutfile{str(residue_number_ros):0>5}'), 'w')
                    mutfile.write('total 20')
    
                    # and then a line for each type of AA
                    for AAtype in 'ACDEFGHIKLMNPQRSTVWY':
                        mutfile.write('\n1\n')
                        mutfile.write(strucdata[chain_id][residue_number_ros-1] + ' ' + str(residue_number_ros) + ' ' + AAtype )
                    mutfile.close()
        return(check2)

##########################################################################
#                                     Creating sbatch relax
##########################################################################

    def rosetta_sbatch_relax(self, folder, relaxfile='', sys_name='', partition='sbinlab'):
        structure_path = os.path.join(folder.relax_input, 'input.pdb')

        path_to_sbatch = os.path.join(
            folder.relax_input, 'rosetta_relax.sbatch')
        if relaxfile == '':
            path_to_relaxflags = os.path.join(
                folder.relax_input, 'relax_flagfile')
        else:
            path_to_relaxflags = relaxfile

        with open(path_to_sbatch, 'w') as fp:
            fp.write(f'''#!/bin/sh
#SBATCH --job-name=relax_{sys_name}
#SBATCH --time=10:00:00
#SBATCH --mem 5000
#SBATCH --partition={partition}

# launching rosetta relax 
''')
            fp.write((f'{os.path.join(rosetta_paths.path_to_rosetta, f"bin/relax.{rosetta_paths.Rosetta_extension}")} '
                      f'-s {structure_path} '
                      f'-relax:script {rosetta_paths.path_to_parameters}/cart2.script @{path_to_relaxflags}'))
        logger.info(path_to_sbatch)
        return(path_to_sbatch)

##########################################################################
#                                     Creating sbatch parse relax
##########################################################################

    def parse_relax_sbatch(self, folder, sys_name='', partition='sbinlab'):
        path_to_parse_relax_script = os.path.join(
            rosetta_paths.path_to_stability_pipeline, 'relax_parse_results.py')

        path_to_sbatch = os.path.join(folder.relax_input, 'parse_relax.sbatch')
        with open(path_to_sbatch, 'w') as fp:
            fp.write(f'''#!/bin/sh
#SBATCH --job-name=parse_relax_rosetta_{sys_name}
#SBATCH --time=00:20:00
#SBATCH --mem 5000
#SBATCH --partition={partition}

# launching parsing script 
''')
            fp.write(f'python {path_to_parse_relax_script} {folder.relax_run} {folder.relax_output} {folder.ddG_input}')
        logger.info(path_to_sbatch)
        return path_to_sbatch

##########################################################################
#                                     Creating sbatch Rosetta cartesian
##########################################################################

    def write_rosetta_cartesian_ddg_sbatch(self, folder, input_mutfiles='', ddgfile='', sys_name='', partition='sbinlab'):
        path_to_sbatch = os.path.join(folder.ddG_input, 'rosetta_ddg.sbatch')
        structure_path = os.path.join(folder.ddG_input, 'input.pdb')
        relax_input = os.path.join(folder.ddG_input, 'input.pdb')
        if input_mutfiles == '':
            input_mutfiles = os.path.join(folder.ddG_input, 'mutfiles')
        if ddgfile == '':
            # path_to_ddgflags = os.path.join(
            # rosetta_paths.path_to_parameters, 'cartesian_ddg_flagfile')
            path_to_ddgflags = os.path.join(folder.ddG_input, 'ddg_flagfile')
        else:
            path_to_ddgflags = ddgfile

        muts = os.listdir(input_mutfiles)

        with open(path_to_sbatch, 'w') as fp:
            fp.write(f'''#!/bin/sh 
#SBATCH --job-name=cartesian_{sys_name}
#SBATCH --array=0-{len(muts)}
#SBATCH --time=32:00:00
#SBATCH --mem 5000
#SBATCH --partition={partition}
#SBATCH --nice 
LST=(`ls {input_mutfiles}/mutfile*`)
OFFSET=0 
INDEX=$((OFFSET+SLURM_ARRAY_TASK_ID))
echo $INDEX

# launching rosetta 
''')
            fp.write((f'{os.path.join(rosetta_paths.path_to_rosetta, f"bin/cartesian_ddg.{rosetta_paths.Rosetta_extension}")} '
                      f'-s {structure_path} -ddg:mut_file ${{LST[$INDEX]}} '
                      f'-out:prefix ddg-$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID @{path_to_ddgflags}'))
        logger.info(path_to_sbatch)
        return path_to_sbatch

##########################################################################
#                                     Creating sbatch parse ddg
##########################################################################

    def write_parse_cartesian_ddg_sbatch(self, folder, partition='sbinlab'):
        score_sbatch_path = os.path.join(folder.ddG_input, 'parse_ddgs.sbatch')
        with open(score_sbatch_path, 'w') as fp:
            fp.write(f'''#!/bin/sh 
#SBATCH --job-name=collect_rosetta_ddgs_{self.sys_name} 
#SBATCH --array=1 
#SBATCH --nodes=1 
#SBATCH --time=0:10:00 
#SBATCH --partition={partition} 

#This sbatch script launches the parse parse_rosetta_ddgs function, from the parse_cartesian_ddgs 
''')
            fp.write((f'python3 {rosetta_paths.path_to_stability_pipeline}/parse_rosetta_ddgs.py '
                      f'{self.sys_name} {self.chain_id} {self.fasta_seq} {folder.ddG_run} {folder.ddG_output}'))
        return score_sbatch_path