"""helper.py contains helpful functions.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-05-04

"""

# Standard library imports
import logging as logger
import os
import shutil
import subprocess
import urllib.request
import sys
from os import listdir
from os.path import isfile, join
import json

# Third party imports
import numpy as np
from scipy import stats

# Local application imports
from get_memory_stats import check_memory
from prism_rosetta_parser import rosetta_to_prism


def create_symlinks(source_file, dist_folder, name=''):
    # WIP
    # import stat
    if name == '':
        basename = os.path.basename(source_file)
    else:
        basename = name
    dist_file = os.path.join(dist_folder, basename)
    # st = os.stat(source_file)
    # os.chown(source_file, st[stat.ST_UID], st[stat.ST_GID])
    os.symlink(source_file, dist_file, True)
    logger.info(f'Symbolic link created: {dist_file} --> {source_file}')
    return dist_file


def create_copy(source_file, dist_folder, name='', directory=False):
    if name == '':
        basename = os.path.basename(source_file)
    else:
        basename = name
    dist_file = os.path.join(dist_folder, basename)
    if directory:
        try:
            shutil.copytree(source_file, dist_file)
        except FileExistsError:
            logger.warn(f'Directory {dist_file} already exists. No files will be copied.')
    else:
        shutil.copy(source_file, dist_file)
    logger.info(f'Hardcopy created: {dist_file} --> {source_file}')

    return dist_file


def find_copy(searchfolder, searchword, resultfolder, resultname):
    for file in os.listdir(searchfolder):
        if file.endswith(searchword):
            return create_copy(os.path.join(searchfolder, file), resultfolder, name=resultname)


class AttrDict(dict):

    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self


def get_mut_dict(mutfile):
    # extracts mutations & generates the mutation dictionary
    logger.info(
        'Extract information from mutfile and generate dictionary')
    mut_dic = {}
    with open(mutfile, 'r') as fp:
        for line in fp:
            muts = line.split()
            mut_dic[muts[1]] = muts[0] + muts[2]
    return mut_dic


def longest_common_subsequence(str1, str2):
    # adapted from
    # https://www.codespeedy.com/find-longest-common-subsequence-in-python/
    a = len(str1)
    b = len(str2)
    string_matrix = [[0 for i in range(b + 1)] for i in range(a + 1)]
    for i in range(1, a + 1):
        for j in range(1, b + 1):
            if i == 0 or j == 0:
                string_matrix[i][j] = 0
            elif str1[i - 1] == str2[j - 1]:
                string_matrix[i][j] = 1 + string_matrix[i - 1][j - 1]
            else:
                string_matrix[i][j] = max(
                    string_matrix[i - 1][j], string_matrix[i][j - 1])
    index = string_matrix[a][b]
    res = [''] * index
    i = a
    j = b
    while i > 0 and j > 0:
        if str1[i - 1] == str2[j - 1]:
            res[index - 1] = str1[i - 1]
            i -= 1
            j -= 1
            index -= 1
        elif string_matrix[i - 1][j] > string_matrix[i][j - 1]:
            i -= 1
        else:
            j -= 1
    return ''.join(res)


def extract_by_uniprot_fasta(keyword):
    # extract information from uniprot
    url_base = "https://www.uniprot.org/uniprot/"
    search_params = "?query=reviewed:yes" +\
        "+AND+accession:" + keyword
    return_params = "+&format=tab&columns=id,sequence,organism,entry%20name"
    url = url_base + search_params + return_params

    data_array = []
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8')
        unprocessed = line[:-1].split('\t')
        data_array.append(unprocessed)

    return data_array


def extract_uniprot_fasta(keyword):
    # extract information from uniprot
    url_base = "https://www.uniprot.org/uniprot/"
    search_params = "?query=reviewed:yes" +\
        "+AND+gene:" + keyword +\
        "+AND+organism:human"
    return_params = "+&format=tab&columns=id,sequence,organism,entry%20name"
    url = url_base + search_params + return_params

    data_array = []
    for line in urllib.request.urlopen(url):
        line = line.decode('utf-8')
        unprocessed = line[:-1].split('\t')
        data_array.append(unprocessed)

    return data_array


def drop_numerical_outliers(df, variant_col='variant', score_col='score', z_thresh=3):
    # Constrains will contain `True` or `False` depending on if it is a value
    # below the threshold.
    constrains = df[[score_col]].select_dtypes(include=[np.number]) \
        .apply(lambda x: np.abs(stats.zscore(x)) < z_thresh) \
        .all(axis=1)
    # Drop (inplace) values set to be rejected
    logger.info(df[variant_col][~constrains].tolist())
    logger.info(df[score_col][~constrains].tolist())
    df.drop(df.index[~constrains], inplace=True)

    df.sort_values('variant', inplace=True)
    df.reset_index(drop=True, inplace=True)
    return df


def read_fasta(fasta_file):
    fastau_file = open(fasta_file, 'r')
    fastau_lines = fastau_file.readlines()
    fastau_file.close()
    uniprot_seq = ''
    for line in fastau_lines[1:]:
        uniprot_seq = uniprot_seq + line.strip()
    return(uniprot_seq)

def check_path(path):
    if path and (os.path.isfile(path) or os.path.isdir(path)):
        path=os.path.abspath(path)
    return(path)

def read_slurms(path, printing=False):
    files = [f for f in listdir(path) if isfile(join(path, f))]
    mypath=path
    warn = []
    warnings = 0
    canc = []
    cancels = 0
    slurm_file_canc = []
    slurm_file_warn = []
    nums_warn = []
    nums_canc = []
    for file in files:
        if str(file[0:5]) == "slurm":
            path = join(mypath, file)
            with open(path) as f:
                for line in f:
                    if "WARNING" in line:
                        num_warn = file.split('_')
                        num_warn = num_warn[1].split('.')
                        num_warn = num_warn[0]

                        nums_warn.append(num_warn)
                        warn.append(line)
                        warnings += 1
                        slurm_file_warn.append(file)
                    if "CANCELLED" in line:

                        num_canc = file.split('_')
                        num_canc = num_canc[1].split('.')
                        num_canc = num_canc[0]

                        nums_canc.append(num_canc)
                        slurm_file_canc.append(file)
                        canc.append(line)
                        cancels += 1

    nums_canc = list(dict.fromkeys(nums_canc))
    nums_warn = list(dict.fromkeys(nums_warn))
    print("Number of cancels = ", cancels)
    print(nums_canc)
    #print("Number of warnings = ",warnings)
    # print(nums_warn)
    if printing == True:
        if warnings != 0:
            with open(join(mypath, "warningfile"), 'w') as errorfile:
                errorfile.write("Number of warnings = "+str(warnings)+'\n')
                errorfile.write("Files= "+str(nums_warn)+'\n')
                for n, m in zip(warn, slurm_file_warn):
                    errorfile.write(m+'\n')
                    errorfile.write(n+'\n')
        if cancels != 0:
            with open(join(mypath, "cancelfile"), 'w') as cancelfile:
                cancelfile.write("Number of cancels = "+str(cancels)+'\n')
                cancelfile.write("Files = "+str(nums_canc)+'\n')
                for n, m in zip(canc, slurm_file_canc):
                    cancelfile.write(m+'\n')
                    cancelfile.write(n+'\n')


def ddG_postprocessing(in_ddg, out_ddg, sec_all=None, startnr=1):
    """does sort the variants and optionally shifts the numbering, depending on the sequence"""
    result_list = []
    start_resi = 0
    if sec_all:
        seqdic = sec_all['resdata']
        if startnr ==1:
            minkey = min(sec_all['resdata_reverse'], key=sec_all['resdata_reverse'].get)
            start_resi = int(minkey)-1
    with open(in_ddg, 'r') as fp2, open(out_ddg, 'w') as fp3:
        for line in fp2:
            line_str = line.split(',')
            if sec_all:
                new_number = str(seqdic[line_str[0][1:-1]][1]-start_resi)
                new_variant = f'{line_str[0][0]}{new_number}{line_str[0][-1]},{line_str[1]},{line_str[2]}'
                result_list.append([new_number, line_str[0][-1], new_variant])
            else:
                result_list.append([line_str[0][1:-1], line_str[0][-1], line])
                
        sorted_list = sorted(result_list, key=lambda x: (int(x[0]), x[1]))
        for elem in sorted_list:
            fp3.write(f'{elem[-1]}')


def shift_pdb_numbering(in_pdb, out_pdb, sec_all, startnr=1):
    result_list = []
    start_resi = 0
    if sec_all:
        seqdic = sec_all['resdata']
        if startnr ==1:
            minkey = min(sec_all['resdata_reverse'], key=sec_all['resdata_reverse'].get)
            start_resi = int(minkey)-1
    with open(in_pdb, 'r') as fp2, open(out_pdb, 'w') as fp3:
        for line in fp2:
            if line.startswith('ATOM'):
                new_number = str(seqdic[line[22:26].strip()][1]-start_resi)
                line_str = ' '*(4-len(new_number))
                fp3.write(f'{line[:22]}{line_str}{new_number}{line[26:]}')
            else:
                fp3.write(line)


def add_ddG_to_bfactor(in_pdb, in_ddg, out_pdb, threshold=0.0):
    with open(in_pdb, 'r') as fp, open(in_ddg, 'r') as fp2, open(out_pdb, 'w') as fp3:
        ddG_resi_dic = {}
        for line in fp2:
            line_str = line.split(',')
            resid = line_str[0][1:-1]
            if resid in ddG_resi_dic.keys():
                ddG_resi_arr = ddG_resi_dic[resid]
                ddG_resi_arr.append(float(line_str[1]))
                ddG_resi_dic[resid] = ddG_resi_arr
            else:
                ddG_resi_dic[resid] = [float(line_str[1])]
        for line in fp:
            if line.startswith('ATOM'):
                resid = line[22:26].strip()
                if resid in ddG_resi_dic.keys():
                    bfactor = round(len(np.where( np.array(ddG_resi_dic[resid]) > threshold ))/len(ddG_resi_dic[resid]),2)
                else:
                    bfactor = 0.00
                bfac_str = '%.2f'%(bfactor)
                fp3.write(f'{line[:62]}{bfac_str}{line[66:]}')
            else:
                fp3.write(line)


def generate_output(folder, output_name='ddG.out', sys_name='', version=1, prims_nr='XXX', chain_id='A', output_gaps=False, bfac=True):
    ddg_file = os.path.join(folder.ddG_run, output_name)
    pdb_file_raw = os.path.join(folder.ddG_input, 'input.pdb')
    if bfac:
        pdb_file = os.path.join(folder.ddG_run, 'input+bfac.pdb')
        add_ddG_to_bfactor(pdb_file_raw, ddg_file, pdb_file, threshold=3.5)
    else:
        pdb_file = pdb_file_raw
    seqdicfile = os.path.join(folder.prepare_checking, 'structure_input.json')
    with open(seqdicfile, 'r') as fp:
        sec_all = json.load(fp)
        rosetta_seq = sec_all['strucdata'][chain_id][0]
        sequence_pdbnbr = sec_all['strucdata'][chain_id][2]
        seqdic = sec_all['resdata']
        minkey = min(sec_all['resdata_reverse'], key=sec_all['resdata_reverse'].get)
        first_residue_number = int(minkey)

    ddg_sorted_file = os.path.join(folder.ddG_run, f'{output_name[:-4]}_sorted_continuous{output_name[-4:]}')
    ddG_postprocessing(ddg_file, ddg_sorted_file, sec_all=None, startnr=1)
    prims_file = os.path.join(folder.ddG_output, f'prims_rosetta_{prims_nr}_{sys_name}.txt')
    rosetta_to_prism(ddg_sorted_file, prims_file, rosetta_seq, rosetta_info=None,
                     version=version, sys_name=sys_name, first_residue_number=1)
    create_copy(prims_file, folder.output)
    create_copy(pdb_file, folder.output, name=f'{sys_name}_final.pdb')

    if output_gaps:
        ini_d = True
        sequence = ''
        for elem in sequence_pdbnbr:
          if elem == '-' and ini_d:
            pass
          else:
            ini_d = False
            sequence += elem
        if first_residue_number >1:
            ddg_shifted_gap_file = os.path.join(folder.ddG_run, f'{output_name[:-4]}_gap-shifted{output_name[-4:]}')
            prims_gap_shifted_file = os.path.join(folder.ddG_output, f'prims_rosetta_{prims_nr}_{sys_name}_gap-shifted.txt')
            ddG_postprocessing(ddg_file, ddg_shifted_gap_file, sec_all=sec_all, startnr=first_residue_number)
            rosetta_to_prism(ddg_shifted_gap_file, prims_gap_shifted_file, sequence, rosetta_info=None,
                             version=version, sys_name=sys_name, first_residue_number=first_residue_number)
            create_copy(prims_gap_shifted_file, folder.output)

            pdb_gap_shifted_file = os.path.join(folder.ddG_output, 'relaxed_gap_shifted.pdb')
            shift_pdb_numbering(pdb_file, pdb_gap_shifted_file, sec_all, startnr=first_residue_number)
            create_copy(pdb_gap_shifted_file, folder.output, name=f'{sys_name}_final_gap_shifted.pdb')
        else:
            logger.warn(f'original pdb-numbering starts with residue-nr <=0 ({first_residue_number}) - no fitting file generated')
        ddg_gap_file = os.path.join(folder.ddG_run, f'{output_name[:-4]}_gap{output_name[-4:]}')
        print(sequence)
        ddG_postprocessing(ddg_file, ddg_gap_file, sec_all=sec_all, startnr=1)
        prims_gap_file = os.path.join(folder.ddG_output, f'prims_rosetta_{prims_nr}_{sys_name}-gap.txt')
        rosetta_to_prism(ddg_gap_file, prims_gap_file, sequence, rosetta_info=None,
                         version=version, sys_name=sys_name, first_residue_number=1)
        create_copy(prims_gap_file, folder.output)

        pdb_gap_file = os.path.join(folder.ddG_output, 'relaxed_gap.pdb')
        shift_pdb_numbering(pdb_file, pdb_gap_file, sec_all, startnr=1)
        create_copy(pdb_gap_file, folder.output, name=f'{sys_name}_final_gap.pdb')


def runtime_memory_stats(ddG_run_folder):
    #This part collects informations about runtime and memory use
    try:
        job_id_pos= join(ddG_run_folder, 'job_id_ddg.txt')
        with open(job_id_pos, 'r') as job_id_file:
            ddg_process_id=str(job_id_file.readlines()[-1])
            print(ddg_process_id)
         
        #Get stats 
        shell_command = f'sacct --format="JobID,Start,End,CPUTime,ReqMem,MaxRSS,MaxVMSize,AveVMSize,JobName" > memory_usage_{ddg_process_id}.log'
        subprocess.call(shell_command, cwd=ddG_run_folder, shell=True)
                                                
        check_memory(ddg_process_id,ddG_run_folder) 
        
    except:     
        print('no memory file')
    
