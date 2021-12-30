"""structure_inputs.py includes several preprocessing and script generating functions.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-04-15

"""

# Standard library imports
import glob
import logging 
import os
import subprocess
import sys

# Third party imports
import numpy as np

# Local application imports
import pdb_to_fasta_seq
import rosetta_paths
from AnalyseStruc import get_structure_parameters
from helper import read_fasta

class structure:
    
    def __init__(self,chain_id,name,folder,prep_struc,run_struc,logger,input_dict,uniprot_accesion=''):
        #Initiating parameters for internal use in the script
        try:
            self.chain_id=chain_id
            self.sys_name=name
            self.prep_struc=prep_struc
            self.struc_dic= get_structure_parameters(
                folder.prepare_checking, self.prep_struc)
            self.run_struc=run_struc
            self.folder=folder
            self.logger=logger
            self.input_dict=input_dict
            if uniprot_accesion != '':
                if os.path.isfile(uniprot_accesion):
                    self.uniprot_seq = read_fasta(uniprot_accesion)
                else:
                    self.uniprot_seq = extract_by_uniprot_fasta(uniprot_accesion)[1][1]
            self.name='input'
        except:
            print('Failed initiation. Please check input parameters')
            sys.exit()
            
    def clean_up_and_isolate(self, name='input',ligand=None):
        """This script is for cleaning the pdb file from unwanted entities"""
        
        #This cleans the protein and removes ligands
        if  ligand == None:           
            self.path_to_clean_pdb = rosetta_paths.path_to_clean_pdb
            #Runs shell script
            shell_command = f'python2 {self.path_to_clean_pdb} {self.prep_struc} {self.run_struc}'
            self.logger.info('Running clean_pdb.py script')
            subprocess.call(shell_command, cwd=self.folder.prepare_cleaning, shell=True)
            self.logger.info('end of output from clean_pdb.py')
    
            self.path_to_cleaned_pdb = os.path.join(self.folder.prepare_cleaning, f'{name}_{self.run_struc}.pdb')
            path_to_cleaned_pdb=self.path_to_cleaned_pdb
            #Gets the fasta_seq
            for chain in list(str(self.run_struc)):
                self.path_to_cleaned_fasta = os.path.join(self.folder.prepare_cleaning, f'{name}_{chain}.fasta')
                with open(self.path_to_cleaned_fasta, 'r') as fp:
                    fasta_lines = fp.readlines()
                    self.fasta_seq = ''
            
                    for line in fasta_lines[1:]:
                        self.fasta_seq = self.fasta_seq + line.strip()
                    
        #This cleans the protein but keeps the ligands          
        if ligand == True:
            self.path_to_clean_pdb = rosetta_paths.path_to_clean_keep_ligand
            #Runs shell script
            shell_command = f'python2 {self.path_to_clean_pdb} {self.prep_struc} {self.chain_id}'
            self.logger.info('Running clean_pdb_keep_ligand.py script')
            subprocess.call(shell_command, cwd=self.folder.prepare_cleaning, shell=True)
            self.logger.info('end of output from clean_pdb_keep_ligand.py')
            path_to_cleaned_pdb = os.path.join(self.folder.prepare_cleaning, f'{name}.pdb{self.chain_id}.pdb')
        #Creates struc.json from cleaned pdb file           
        struc_dic_cleaned= get_structure_parameters(
            self.folder.prepare_cleaning, path_to_cleaned_pdb,printing=False)
        self.struc_dic_cleaned= struc_dic_cleaned
        return(path_to_cleaned_pdb,struc_dic_cleaned)



    def muscle_align_to_uniprot(self, uniprot_sequence,name='input'):
        """This script aligns the uniprot sequence to the structure sequence using muscle. This should be changed from muscle to biopairwise2 eventially  """
        
        #Making fasta file with both structure sequence and uniprot sequence
        path_to_muscle = rosetta_paths.path_to_muscle
        self.path_to_fasta = os.path.join(self.folder.prepare_checking, 'fasta_file.fasta')
        self.path_to_alignment = os.path.join(self.folder.prepare_checking, 'alignment.txt')
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

        self.path_to_index_string = os.path.join(self.folder.prepare_checking, 'uniprot_index_list.txt')
        with open(self.path_to_index_string, 'w') as index_file:
            index_list_as_string = '\n'.join(structure_index_numbers)
            index_file.write(str(index_list_as_string))

        path_to_index_string = self.path_to_index_string
        return(path_to_index_string)

    def make_mut_input(self, mutation_input=None, struc_dic=None, alignment=None, fasta=None, chain_id=None, AA='ACDEFGHIKLMNPQRSTVWY'):
        """This script creates a mutation input file, which contains all variants for all position. If a mutation input file is provided, then no file will be made"""
        
        #Creating mutation_input file for desired chain
        chain_id = self.chain_id
        if mutation_input == None:
            path_to_mutation_input = self.folder.prepare_input    
            if struc_dic == None:
                resdata = self.struc_dic_cleaned["resdata"]
                strucdata = self.struc_dic_cleaned["strucdata"]
            else:
                resdata = struc_dic["resdata"]
                strucdata = struc_dic["strucdata"]
            with open(path_to_mutation_input + "/mutation_input.txt", 'w') as mutation_input:
                for n in range(1, len(resdata)+1):                   
                    if resdata[n][2] == str(chain_id):
                        wt = resdata[n][0]
                        pos = resdata[n][1]
                        mutation_input.write(f'{wt} {pos}  {AA} \n')
    
    def make_mutfiles(self, mutation_input, mutation_mode):
        """This scripts generation all the mutfiles for the Rosetta run. It uses a mutational input file to dictate which files to create"""
  

        check2 = False
        mut_dic = {}
        resdata = self.struc_dic_cleaned["resdata"]
        strucdata = self.struc_dic_cleaned["strucdata"]        
        path_to_alignment = self.path_to_index_string
        alignment = np.loadtxt(path_to_alignment)
        
        alignment_dic = {}
        for n in enumerate(alignment):
            alignment_dic[int(n[1])] = int(n[0] + 1)

        if mutation_input != "proceed":

            if mutation_mode == 'all':
                self.make_mut_input()
                path_to_mutation_input = self.folder.prepare_input    
                mutation_input=os.path.join(path_to_mutation_input,'mutation_input.txt') 

            mutate = {}
            self.logger.debug("Printing point mutations [WT, ResID, Mutations]")
            if os.path.isdir(mutation_input):
                for file in glob.glob(os.path.join(mutation_input, '*')):
                    with open(file, 'r') as f:
                        total_muts = int(next(f).split()[1].split('\\n')[0])
                        lines = f.readlines()
                        key = file.split('/')[-1].split("mutfile")[1].split('_')
                        res_nums = len(key)
                        mutsNum = int(total_muts/res_nums)
                        wt=[]
                        val=[]            
                        for resii in range(res_nums):
                            single_var = []
                            for singlemuts in range(mutsNum):
                                target = lines[(1+1*singlemuts)+(resii+(singlemuts*res_nums))].split()
                                single_var.append(target[2])
                                single_wt = target[0]
                            val.append("".join(single_var))
                            wt.append(single_wt)
                        mutate["_".join(key)] = "_".join(wt), "_".join(val)
            else:
                with open(mutation_input) as f:
                    first_line = next(f).split()[0]
                with open(mutation_input, 'r') as f:
                    if first_line == 'total':
                        save = False
                        next(f)
                        for lines in f:
                            lines = lines.split()
                            if len(lines) == 1:
                                if save:
                                    if "_".join(key) in mutate.keys():
                                        old_wt, old_val = mutate["_".join(key)]
                                        new_val = []
                                        for indl, valu in enumerate(old_val.split("_")):
                                            new_val.append(f"{valu}{val[indl]}")
                                        mutate["_".join(key)] = "_".join(wt), "_".join(new_val)
                                    else:
                                        mutate["_".join(key)] = "_".join(wt), "_".join(val)
                                else:
                                    save = True
                                key=[]
                                wt=[]
                                val=[]
                            else:
                                key.append(lines[1])
                                wt.append(lines[0])
                                val.append(lines[2])
                        if save:
                            if "_".join(key) in mutate.keys():
                                old_wt, old_val = mutate["_".join(key)]
                                new_val = []
                                for indl, valu in enumerate(old_val.split("_")):
                                    new_val.append(f"{valu}{val[indl]}")
                                mutate["_".join(key)] = "_".join(wt), "_".join(new_val)
                            else:
                                mutate["_".join(key)] = "_".join(wt), "_".join(val)
                        else:
                            save = True
                    else:
                        for line in f:
                            lines = line.split()
                            if len(lines) == 3:
                                (wt, key, val) = line.split()
                                mutate[int(key)] = wt, val
                                c = mutate[int(key)][0] in list(mutate[int(key)][1])
                                if c == True:
                                    mutate[int(key)] = wt, val
                                else:
                                    val = wt + val
                                    mutate[int(key)] = wt, val
                            else:
                                key=[]
                                wt=[]
                                val=[]
                                for multi_muts in range(int(len(lines)/3)):
                                    wt_single = lines[0+(3*multi_muts)]
                                    var_single = lines[2+(3*multi_muts)]
                                    if not wt_single in var_single:
                                        var_single += wt_single
                                    #var_single = ''.join(set(var_single))
                                    wt.append(wt_single)
                                    key.append(lines[1+(3*multi_muts)])
                                    val.append(var_single)
                                mutate["_".join(key)] = "_".join(wt), "_".join(val)
                    # remove duplicate mutants and add WT
                    for key in mutate.keys():
                        wt_list = mutate[key][0].split('_')
                        mut_list = mutate[key][1].split('_')
                        if type(key) == str:
                            resi_list = key.split('_')
                        else:
                            resi_list = [key]
                        joined_muts = []
                        for mut_ind in range(0,len(mut_list[0])):
                            tmp_joined = []
                            for resi_indi, resi in enumerate(resi_list):
                                tmp_joined.append(mut_list[resi_indi][mut_ind])
                            joined_muts.append(":".join(tmp_joined))
                        tmp_joined = []
                        for resi_indi, resi in enumerate(resi_list):
                            tmp_joined.append(wt_list[resi_indi])
                        joined_muts.append(":".join(tmp_joined))
                        joined_muts = sorted(list(set(joined_muts)))
                        #split back to mutate dic format
                        mut_dic = {}
                        for muti_ind, muti in enumerate(joined_muts):
                            for mut_ind, mut in enumerate(joined_muts[muti_ind].split(':')):
                                if mut_ind in mut_dic.keys():
                                    mut_dic[mut_ind].append(mut)
                                else:
                                    mut_dic[mut_ind] = [mut]
                        mut_list = []
                        for sub_key in mut_dic.keys():
                            mut_list.append("".join(mut_dic[sub_key]))
                        mutate[key] = (mutate[key][0], "_".join(mut_list))

            with open(os.path.join(self.folder.prepare_cleaning, 'mutation_clean.txt'), 'w') as fp:
                i = 10000
                for residue_number in mutate:
                    if mutation_mode == 'all':
                        residue_number_ros = alignment_dic[residue_number]
                    else:
                        try:
                            residue_number_ros = []
                            for res_num in str(residue_number).split('_'):
                                residue_number_ros.append(self.struc_dic['resdata_reverse'][int(res_num)])
                        except Exception as e:
                            self.logger.warning(f"Residue {residue_number} present in the mutation file is not resolved in the structure. No output generated for it.")
                            continue
                    self.logger.debug(f"{residue_number_ros}  {residue_number}")
                    #This checks that the position in the mutfile is the correct one in the fasta sequence
                    if type(residue_number_ros)!=list:
                        residue_number_ros = [residue_number_ros]
                    if len(residue_number_ros)==1:
                        residue_number_ros = int(residue_number_ros[0])
                        check = self.fasta_seq[residue_number_ros-1] in list(
                            mutate[residue_number][0])
                        if check == False:
                            check2 = True
                            self.logger.warning(f'Missmatch{self.fasta_seq[residue_number_ros-1]}, {residue_number},{mutate[residue_number][0]}')
                        final_list = []
                        for num in mutate[residue_number][1]:                    
                            if num not in final_list:
                                final_list.append(num)                
                        self.logger.debug(f"{mutate[residue_number][0]} {residue_number}  {''.join(final_list)}")

                        with open(os.path.join(self.folder.prepare_mutfiles, f'mutfile{str(residue_number_ros):0>5}'), 'w') as mutfile:
                            mutfile.write('total ' + str(len(final_list)))
                            mut_dic[str(residue_number_ros)] = "".join(final_list)
                            fp.write(f'{self.fasta_seq[residue_number_ros - 1]} {residue_number_ros} {"".join(final_list)}\n')
                            for AAtype in final_list:
                                mutfile.write('\n1\n')
                                mutfile.write(self.fasta_seq[
                                              residue_number_ros - 1] + ' ' + str(residue_number_ros) + ' ' + AAtype)
                    else:
                        for res_index, residue in enumerate(residue_number_ros):
                            check = self.fasta_seq[int(residue)-1] in list(
                                mutate[residue_number][0].split("_")[res_index])
                            if check == False:
                                check2 = True
                                self.logger.warning(f'Missmatch{self.fasta_seq[int(residue)-1]}, {residue},{mutate[int(residue)][0]}')
                            fp.write(f'{self.fasta_seq[int(residue) - 1]} {int(residue)} {mutate[residue_number][1].split("_")[res_index]} ')
                        fp.write('\n')

                        with open(os.path.join(self.folder.prepare_mutfiles, f'mutfile{str(residue_number)}'), 'w') as mutfile:#f'mutfile{str(i)}'), 'w') as mutfile:
                            mutfile.write('total ' + str(len(residue_number_ros)*len(mutate[residue_number][1].split('_')[0]))+'\n')
                            mut_dic[residue_number] = "".join(mutate[residue_number][1])
                            for res_rounds in range(len(mutate[residue_number][1].split('_')[0])):
                                mutfile.write(f'{len(residue_number_ros)}\n')
                                
                                for indi, res in enumerate(residue_number_ros):
                                    towrite = self.fasta_seq[int(res) - 1] + ' ' + str(res) + ' ' + mutate[residue_number][1].split("_")[indi].split()[0][res_rounds]+ '\n'
                                    mutfile.write(towrite)
                        i += 1

        return(check2, mut_dic)


    def rosetta_sbatch_relax(self, folder, relaxfile='', sys_name='', partition='sbinlab'):
        """This script creates the rosetta_relax.sbatch script"""
                                  
        structure_path = os.path.abspath(os.path.join(self.folder.relax_input, 'input.pdb'))
        path_to_sbatch = os.path.join(
            self.folder.relax_input, 'rosetta_relax.sbatch')
        if relaxfile == '':
            path_to_relaxflags = os.path.join(
                folder.relax_input, 'relax_flagfile')
        else:
            path_to_relaxflags = relaxfile
                                  
        with open(path_to_sbatch, 'w') as fp:
            fp.write(f'''#!/bin/sh
#SBATCH --job-name={self.sys_name}_relax
#SBATCH --time=24:00:00
#SBATCH --mem 5000
#SBATCH --array=0-19
#SBATCH --partition={partition}

# launching rosetta relax 
''')
            fp.write((f'{os.path.join(rosetta_paths.path_to_rosetta, f"bin/relax.{rosetta_paths.Rosetta_extension}")} '
                      f'-s {structure_path} '
                      f'-relax:script {rosetta_paths.path_to_data}/sp/cart2.script @{path_to_relaxflags}'
                     ' -out:prefix $SLURM_ARRAY_TASK_ID-'))
        self.logger.info(path_to_sbatch)
        return(path_to_sbatch)


    def parse_relax_sbatch(self, folder, sys_name='', partition='sbinlab', sc_name='score_bn15_calibrated', mp_multistruc=0):
        """This script creates the parse_relax.sbatch script"""
                                  
        path_to_parse_relax_script = os.path.join(
        rosetta_paths.path_to_stability_pipeline, 'relax_parse_results.py')
        path_to_sbatch = os.path.join(self.folder.relax_input, 'parse_relax.sbatch')
        with open(path_to_sbatch, 'w') as fp:
            fp.write(f'''#!/bin/sh
#SBATCH --job-name={self.sys_name}_p-relax
#SBATCH --time=00:20:00
#SBATCH --mem 5000
#SBATCH --partition={partition}

# launching parsing script 
''')
            if mp_multistruc == 0:
                fp.write(f'python {path_to_parse_relax_script} {folder.relax_run} {folder.relax_output} {folder.ddG_input} {sc_name}')
            else:
                fp.write(f'python {path_to_parse_relax_script} {folder.relax_run} {folder.relax_output}' )
                for ddg_subfolder in folder.ddG_input:
                    fp.write(f' {ddg_subfolder}')
                fp.write(f' {sc_name} {mp_multistruc}')
        self.logger.info(path_to_sbatch)
        return path_to_sbatch


    def write_rosetta_cartesian_ddg_sbatch(self, folder, input_mutfiles='', ddgfile='', sys_name='', partition='sbinlab'):
        """This script creates the rosetta_dddg.sbatch script"""
                                  
        path_to_sbatch = os.path.join(self.folder.ddG_input, 'rosetta_ddg.sbatch')
        structure_path = os.path.abspath(os.path.join(self.folder.ddG_input, 'input.pdb'))
        relax_input = os.path.join(self.folder.ddG_input, 'input.pdb')
        if input_mutfiles == '':
            input_mutfiles = os.path.join(self.folder.ddG_input, 'mutfiles')
        if ddgfile == '':
            # path_to_ddgflags = os.path.join(
            # rosetta_paths.path_to_data, 'sp', 'cartesian_ddg_flagfile')
            path_to_ddgflags = os.path.join(self.folder.ddG_input, 'ddg_flagfile')
        else:
            path_to_ddgflags = ddgfile

        muts = os.listdir(input_mutfiles)

        with open(path_to_sbatch, 'w') as fp:
            fp.write(f'''#!/bin/sh 
#SBATCH --job-name={sys_name}_ddg
#SBATCH --array=0-{len(muts)-1}
#SBATCH --time=48:00:00
#SBATCH --mem 2000
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
                      f' -out:prefix ddg-$SLURM_ARRAY_JOB_ID-$SLURM_ARRAY_TASK_ID @{path_to_ddgflags}'))
        self.logger.info(path_to_sbatch)
        return path_to_sbatch


    def write_parse_cartesian_ddg_sbatch(self, folder, partition='sbinlab', output_gaps=False, zip_files=True, sha_tag='', is_MP=False, scale_factor=2.9):
        """This script creates the parse_ddgs.sbatch script"""
                                  
        score_sbatch_path = os.path.join(self.folder.ddG_input, 'parse_ddgs.sbatch')
        structure_input = os.path.join(self.folder.prepare_checking,'structure_input.json')
        with open(score_sbatch_path, 'w') as fp:
            fp.write(f'''#!/bin/sh 
#SBATCH --job-name={self.sys_name}_p-ddg
#SBATCH --array=1 
#SBATCH --nodes=1 
#SBATCH --time=0:10:00 
#SBATCH --partition={partition} 

#This sbatch script launches the parse parse_rosetta_ddgs function, from the parse_cartesian_ddgs 
''')
            fp.write((f'python3 {rosetta_paths.path_to_stability_pipeline}/parser_ddg_v2.py '
                      f'{self.sys_name} {self.chain_id} {self.fasta_seq} {folder.ddG_run} {folder.ddG_output} {structure_input}'
                      f' {folder.ddG_input} {folder.output} {folder.prepare_checking} {output_gaps} {zip_files} {sha_tag} {folder.prepare_cleaning} {is_MP} {scale_factor}'))
        return score_sbatch_path