"""helper.py contains helpful functions.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-05-04

"""

# Standard library imports
import io
from datetime import timedelta
import logging as logger
import glob
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
import pandas as pd
from scipy import stats

# Local application imports
from get_memory_stats import check_memory
from prism_rosetta_parser import rosetta_to_prism
sys.path.insert(1, os.path.join(os.path.dirname(__file__), '../scripts/'))
from pdb_to_prism import rosetta_energy_to_prism


def remove_hetatms(input, output):
    with open(input) as fp, open(output, 'w') as fp2:
        for line in fp:
            if not line.startswith('HETATM'):
                if line.startswith('ATOM'):
                    if line[16] in [' ', 'A']:
                        fp2.write(line)
                else:
                    fp2.write(line)
    return output

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
    logger.info(f'Hardcopy created: {source_file} --> {dist_file}')

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
        path_new=os.path.abspath(path)
        return(path_new)
    else:
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


def ddG_postprocessing(in_ddg, out_ddg, sec_all=None, startnr=1, chain_id='A'):
    """does sort the variants and optionally shifts the numbering, depending on the sequence"""
    result_list = []
    start_resi = 0
    if sec_all:
        seqdic = sec_all['resdata']
        if startnr ==1:
            resis = [ int(key.split(';')[0]) for key in sec_all['resdata_reverse2'].keys() if key.split(';')[1] == chain_id]
            start_resi = min(resis)-1
    with open(in_ddg, 'r') as fp2, open(out_ddg, 'w') as fp3:
        for line in fp2:
            line_str = line.split(',')
            if sec_all:
                new_number = []
                new_variant = []
                for resi in line_str[0].split(':'):
                    new_numb = str(seqdic[resi[1:-1]][1]-start_resi)
                    new_number.append(new_numb)
                    new_variant.append(f'{resi[0]}{new_numb}{resi[-1]}')
                new_variants = f'{":".join(new_variant)},{line_str[1]},{line_str[2]}'
                result_list.append([new_number[-1], line_str[0][-1], new_variants])
            else:
                resinumber = line_str[0].split(':')[0][1:-1]
                resivari = line_str[0].split(':')[0][-1]
                result_list.append([resinumber, resivari, line])
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
    rr=1
    with open(in_pdb, 'r') as fp, open(in_ddg, 'r') as fp2, open(out_pdb, 'w') as fp3:
        ddG_resi_dic = {}
        ddG_resi_array = []
        for line in fp2:
            line_str = line.split(',')
            resid = line_str[0][1:-1]
            if resid in ddG_resi_dic.keys():
                ddG_resi_arr = ddG_resi_dic[resid]
                ddG_resi_arr.append(float(line_str[1]))
                ddG_resi_dic[resid] = ddG_resi_arr
            else:
                ddG_resi_dic[resid] = [float(line_str[1])]
            ddG_resi_array.append(float(line_str[1]))

        for line in fp:
            if line.startswith('ATOM'):
                resid = line[22:26].strip()
                if resid in ddG_resi_dic.keys():
                    #bfactor = round(len(np.where( np.array(ddG_resi_dic[resid]) > threshold ))/len(ddG_resi_dic[resid]),2)
                    #bfactor = round(np.average(ddG_resi_dic[resid]),2)
                    df = ddG_resi_dic[resid]
                    df2 = np.average(df)
                    normalized_df=(df2-np.min(df))/(np.max(df)-np.min(df))
                    bfactor = round(normalized_df,2)
                else:
                    bfactor = -1.00
                bfac_str = '%.2f'%(bfactor)
                fp3.write(f'{line[:61]}{" "*(5-len(bfac_str))}{bfac_str}{line[66:]}')
            elif line.startswith('HETATM'):
                if line.split()[-1]=='X':
                    if rr==1:
                        rr=0
                        bfac_str = '%.2f'%(-1.00)
                    else:
                        bfac_str = '%.2f'%(1.00)
                    fp3.write(f'{line[:61]}{" "*(5-len(bfac_str))}{bfac_str}{line[66:]}')
                else:
                    fp3.write(line)
            else:
                fp3.write(line)


def get_slurm_ids(path):
    slurm_ids = []
    slurms = glob.glob(os.path.join(path, 'relax', 'run', 'slurm-*')) + glob.glob(os.path.join(path, 'ddG', 'run', 'slurm-*'))
    for slurm in slurms:
        slurm_ids.append(slurm.split('-')[-1].split('_')[0].split('.')[0])
    return list(set(slurm_ids))


def convert_mem(x):
    if str(x)[-1]=='K':
        return float(x[:-1])*1024
    elif str(x)[-1]=='M':
        return float(x[:-1])*(1024*1024)
    elif str(x)[-1]=='G':
        return float(x[:-1])*(1024*1024*1024)
    else:
        logger.info(x)
        return float(x)


def convert_mem_to_GB(x):
    return float(x)/(1024*1024*1024)


def calculate_emissions(memory, runTime_hours, runTime_min, runTime_sec, n_cores=1, coreType='CPU', coreModel='any', log=False):
    runTime = runTime_hours + runTime_min/60 + runTime_sec/(60*60)
    location = 'Denmark'
    usage = 1
    PUE = 1.67
    PSF = 1
    selected_platform = 'localServer'
    selected_provider = ''
    existing_state = ''

    ###
    carbonIntensity = 154.44 #gCO2e/kWh - Denmark
    PUE_used = PUE
    corePower = 12.0# TDP_per_core, any=coreModel 
    #
    ref_dic = {
        'memoryPower': 0.3725,
        'passengerCar_EU_perkm': 175,
        'passengerCar_US_perkm': 251,
        'train_perkm': 41,
        'flight_economy_perkm': 171,
        'treeYear': 11000,
        'flight_NY-SF': 570000,
        'flight_PAR-LON': 50000,
        'flight_NYC-MEL': 2310000,
        'streaming_netflix_perhour': 86,
        'google_search': 10,
        'tree_month': 917,
    }


    # Power needed, in Watt
    powerNeeded_core = PUE_used * (n_cores * corePower) * usage
    refValues_dict_memoryPower = ref_dic['memoryPower']
    powerNeeded_memory = PUE_used * (memory * refValues_dict_memoryPower)
    powerNeeded = powerNeeded_core + powerNeeded_memory

    # Energy needed, in kWh (so dividing by 1000 to convert to kW)
    energyNeeded_core = runTime * powerNeeded_core * PSF / 1000
    eneregyNeeded_memory = runTime * powerNeeded_memory * PSF / 1000
    energyNeeded = runTime * powerNeeded * PSF / 1000

    # Carbon emissions: carbonIntensity is in g per kWh, so results in gCO2
    CE_core = energyNeeded_core * carbonIntensity
    CE_memory  = eneregyNeeded_memory * carbonIntensity
    carbonEmissions = energyNeeded * carbonIntensity

    output = dict()
    
    output['coreType'] = coreType
    output['coreModel'] = coreModel
    output['n_cores'] = n_cores
    output['corePower'] = corePower
    output['memory'] = memory
    output['runTime_hours'] = runTime_hours
    output['runTime_min'] = runTime_min
    output['runTime'] = runTime
    output['location'] = location
    output['carbonIntensity'] = carbonIntensity
    output['PUE'] = PUE_used
    output['PSF'] = PSF
    output['selected_platform'] = selected_platform
    output['carbonEmissions'] = carbonEmissions
    output['CE_core'] = CE_core
    output['CE_memory'] = CE_memory
    output['energy_needed'] = energyNeeded
    output['power_needed'] = powerNeeded

    output['n_treeMonths'] = carbonEmissions / ref_dic['treeYear'] * 12

    output['nkm_drivingUS'] = carbonEmissions / ref_dic['passengerCar_US_perkm']
    output['nkm_drivingEU'] = carbonEmissions / ref_dic['passengerCar_EU_perkm']
    output['nkm_train'] = carbonEmissions / ref_dic['train_perkm']

    if carbonEmissions < 0.5 * ref_dic['flight_NY-SF']:
        output['flying_context'] = carbonEmissions / ref_dic['flight_PAR-LON']
        output['flying_text'] = "Paris-London"
    elif carbonEmissions < 0.5 * ref_dic['flight_NYC-MEL']:
        output['flying_context'] = carbonEmissions / ref_dic['flight_NY-SF']
        output['flying_text'] = "NYC-San Francisco"
    else:
        output['flying_context'] = carbonEmissions / ref_dic['flight_NYC-MEL']
        output['flying_text'] = "NYC-Melbourne"
    output['streaming_netflix_perhour'] = carbonEmissions / ref_dic['streaming_netflix_perhour']
    output['google_search'] = carbonEmissions / ref_dic['google_search']

    if log:
        logger.info((f"This algorithm runs in {output['runTime_hours']} h and {output['runTime_min']} min "
                 f"on {output['n_cores']} {output['coreType']} ({output['coreModel']}), "
                 f"which draws {output['energy_needed']:.2f} kWh. Based in {location}, "
                 f"this program produces {output['carbonEmissions']:.2f} g of CO2e, "
                 f"which is equivalent to {output['n_treeMonths']:.2f} tree-months carbon sequestration, "
                 f"{output['nkm_train']:.2f} km in a train, "
                 f"{output['nkm_drivingEU']:.2f} km of driving a passenger car in Europe or "
                 f"{output['streaming_netflix_perhour']:.2f} h of netflix streaming "
                 f"(calculated using green-algorithms.org v1.1)."))

    return output

def generate_emission_stats(test_dir):
    slurm_ids = ",".join(get_slurm_ids(test_dir))
    logger.info(slurm_ids)
    shell_command = f'sacct -j {slurm_ids} --format="JobID,Start,End,CPUTime,NCPUS,Elapsed,ReqMem,MaxRSS,MaxVMSize,AveVMSize,JobName"'
    pipes = subprocess.Popen(shell_command, shell=True, cwd=test_dir,stdout=subprocess.PIPE,stderr=subprocess.PIPE,)
    std_out, std_err = pipes.communicate()
    iostringout = io.StringIO(std_out.decode('utf-8'))
    iostringout.seek(0)
    merged_memory = pd.read_fwf(iostringout, colspecs=[(0,12), (12,32), (32,52), (52,65), (65,75), (75,87), (87,98), (98,107), (107,118), (118,130), (130,150)])

    merged_memory = merged_memory.iloc[1: , :].reset_index(drop=True)
    merged_memory['JobID_parent'] = merged_memory['JobID'].apply(lambda x: x.split('_')[0])
    merged_memory = merged_memory.loc[(merged_memory['JobID'].str.endswith('+'))].reset_index(drop=True)

    def make_datetimes(x):
        s = x.split('-')
        if len(s)>1:
            return timedelta(days=int(x.split('-')[0]), hours=int(x.split('-')[1].split(':')[0]), minutes=int(x.split('-')[1].split(':')[1]), seconds=int(x.split('-')[1].split(':')[2]))
        else:
            return timedelta(hours=int(x.split('-')[0].split(':')[0]), minutes=int(x.split('-')[0].split(':')[1]), seconds=int(x.split('-')[0].split(':')[2]))

    merged_memory['Elapsed'] = merged_memory['Elapsed'].apply(lambda x: make_datetimes(x))
    merged_memory['MaxRSS_byte'] = merged_memory['MaxRSS'].apply(lambda x: convert_mem(x))

    emission_df = []
    sum_mem =[]
    sum_time = []
    sum_cores = []
    for index, row in merged_memory.iterrows():
        memory = convert_mem_to_GB(row['MaxRSS_byte'])
        df = pd.DataFrame(data=[row['Elapsed']], columns=['Elapsed'])
        d = df['Elapsed'].dt.components[['days', 'hours', 'minutes', 'seconds']]
        d['hours'] = d['hours'].add(d.pop('days') * 24)
        df['Elapsed'] = d.astype(str).agg(lambda s: ':'.join(s.str.zfill(2)), axis=1)
        runTime_hours = int(df['Elapsed'][0].split(':')[0])
        runTime_min = int(df['Elapsed'][0].split(':')[1])
        runTime_sec = int(df['Elapsed'][0].split(':')[2])
        n_cores = int(row['NCPUS'])
        if (runTime_hours <0) or (runTime_min <0) or(runTime_sec <0):
            print(runTime_hours, runTime_min, runTime_sec)
        emission_dic = calculate_emissions(memory, runTime_hours, runTime_min, runTime_sec, n_cores=n_cores, coreType='CPU', coreModel='any')
        tmp_df = pd.DataFrame([emission_dic])#data=list(emission_dic.items()), columns = emission_dic.keys())
        emission_df.append(tmp_df)

    emission_df = pd.concat(emission_df)
    output_text = ((f"This algorithm runs on {emission_df['n_cores'].sum()} {emission_df['coreType'].iloc[0]} ({emission_df['coreModel'].iloc[0]}) "
                 f"for ~ {int(emission_df['runTime'].median())} h and {(emission_df['runTime'].median() % 1 ) * 60:.0f} min each, "
                 f"for a total runtime of {int(emission_df['runTime'].sum())} h and {(emission_df['runTime'].sum() % 1 ) * 60:.0f} min, "
                 f"with a memory usage of {emission_df['memory'].sum():.2f} GB, "
                 f"drawing {emission_df['energy_needed'].sum():.2f} kWh. Based in {emission_df['location'].iloc[0]}, "
                 f"this program produces {emission_df['carbonEmissions'].sum():.2f} g of CO2e, "
                 f"which is equivalent to {emission_df['n_treeMonths'].sum():.2f} tree-months carbon sequestration, "
                 f"a train ride of {emission_df['nkm_train'].sum():.2f} km, "
                 f"{emission_df['nkm_drivingEU'].sum():.2f} km of driving a passenger car in Europe or "
                 f"{emission_df['streaming_netflix_perhour'].sum():.2f} h of netflix streaming."))

    memory = convert_mem_to_GB(merged_memory['MaxRSS_byte'].sum())
    df = pd.DataFrame(data=[merged_memory['Elapsed'].sum()], columns=['Elapsed'])
    d = df['Elapsed'].dt.components[['days', 'hours', 'minutes', 'seconds']]
    d['hours'] = d['hours'].add(d.pop('days') * 24)
    df['Elapsed'] = d.astype(str).agg(lambda s: ':'.join(s.str.zfill(2)), axis=1)
    runTime_hours = int(df['Elapsed'][0].split(':')[0])
    runTime_min = int(df['Elapsed'][0].split(':')[1])
    runTime_sec = int(df['Elapsed'][0].split(':')[2])
    output = calculate_emissions(memory, runTime_hours, runTime_min, runTime_sec, n_cores=1, coreType='CPU', coreModel='any')
    output_text2 = ((f"If this would be run on 1 CPU instead of {emission_df['n_cores'].sum()} CPUs, it would "
                     f"draw {output['energy_needed']:.2f} kWh, producing {output['carbonEmissions']:.2f} g of CO2e, "
                     f"which is equivalent to {output['n_treeMonths']:.2f} tree-months carbon sequestration, "
                     f"a train ride of {output['nkm_train']:.2f} km, "
                     f"{output['nkm_drivingEU']:.2f} km of driving a passenger car in Europe or "
                     f"{output['streaming_netflix_perhour']:.2f} h of netflix streaming. (calculated using https://doi.org/10.1002/advs.202100707 v1.1)"))
    emission_out_file = os.path.join(test_dir, 'output', 'emission_stats.txt')
    with open(emission_out_file, 'w') as fp:
        fp.write(output_text)
        fp.write('\n')
        fp.write(output_text2)
        fp.write('\n')
        fp.write(emission_df.to_string())
    logger.info(output_text)
    logger.info(output_text2)

def parse_rosetta_Eres_to_mut(folder, sys_name, prism_nr='XXX', chain='A', uniprot_id='', organism='', version=1, count=True):
    # generate Eres file
    prism_file = os.path.join(folder.ddG_output, f'prism_rosettapdb_{prism_nr}_{sys_name}.txt')
    infile = os.path.join(folder.ddG_input, 'input.pdb')
    rosetta_energy_to_prism(infile, prism_file, sys_name, chain, folder.ddG_run, 
        uniprot_id=uniprot_id, organism=organism, version=version, count=count)
    create_copy(prism_file, folder.output)

def generate_output(folder, output_name='ddG.out', sys_name='', version=1, prism_nr='XXX', chain_id='A', output_gaps=False, bfac=True, zip_files=True, sha_tag='', MP=False, scale=2.9):
    # generate emission stats
    generate_emission_stats(folder.output[:-7])
    if MP:
        span_file = glob.glob(os.path.join(folder.prepare_checking[:-9], 'mp_files', 'membrane_span', '*.span'))[0]
        lipid_file = glob.glob(os.path.join(folder.prepare_checking[:-9], 'mp_files', 'mp_lipid_acc', '*.json'))[0]
    else:
        span_file = ''
        lipid_file = ''

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
        first_residue_number = 1
        for elem in sec_all['strucdata'][chain_id][2]:
            if elem in ['-', 'x', 'X']:
                first_residue_number = first_residue_number+1
            else:
                break

    ddg_sorted_file = os.path.join(folder.ddG_run, f'{output_name[:-4]}_sorted_continuous{output_name[-4:]}')
    ddG_postprocessing(ddg_file, ddg_sorted_file, sec_all=None, startnr=1, chain_id=chain_id)
    prism_file = os.path.join(folder.ddG_output, f'prism_rosetta_{prism_nr}_{sys_name}.txt')
    rosetta_to_prism(ddg_sorted_file, prism_file, rosetta_seq, rosetta_info=None,
                     version=version, sys_name=sys_name, first_residue_number=1, sha_tag=sha_tag, MP=MP, 
                     span_file=span_file, lipid_file=lipid_file, scale=scale)
    create_copy(prism_file, folder.output)
    create_copy(pdb_file, folder.output, name=f'{sys_name}_final.pdb')

    parse_rosetta_Eres_to_mut(folder, sys_name, prism_nr=prism_nr, chain=chain_id, version=version)

    if output_gaps:
        ini_d = True
        sequence = ''
        for elem in sequence_pdbnbr:
          if elem in ['-', 'X', 'x'] and ini_d:
            pass
          else:
            ini_d = False
            sequence += elem
        if first_residue_number >1:
            ddg_shifted_gap_file = os.path.join(folder.ddG_run, f'{output_name[:-4]}_gap-shifted{output_name[-4:]}')
            prism_gap_shifted_file = os.path.join(folder.ddG_output, f'prism_rosetta_{prism_nr}_{sys_name}_gap-shifted.txt')
            ddG_postprocessing(ddg_file, ddg_shifted_gap_file, sec_all=sec_all, startnr=first_residue_number, chain_id=chain_id)
            rosetta_to_prism(ddg_shifted_gap_file, prism_gap_shifted_file, sequence, rosetta_info=None,
                             version=version, sys_name=sys_name, first_residue_number=first_residue_number, sha_tag=sha_tag, MP=MP, 
                             span_file=span_file, lipid_file=lipid_file, scale=scale)
            create_copy(prism_gap_shifted_file, folder.output)

            pdb_gap_shifted_file = os.path.join(folder.ddG_output, 'relaxed_gap_shifted.pdb')
            shift_pdb_numbering(pdb_file, pdb_gap_shifted_file, sec_all, startnr=first_residue_number)
            create_copy(pdb_gap_shifted_file, folder.output, name=f'{sys_name}_final_gap_shifted.pdb')
        else:
            logger.warn(f'original pdb-numbering starts with residue-nr <=1 ({first_residue_number}) - no fitting file generated')
        ddg_gap_file = os.path.join(folder.ddG_run, f'{output_name[:-4]}_gap{output_name[-4:]}')
        ddG_postprocessing(ddg_file, ddg_gap_file, sec_all=sec_all, startnr=1, chain_id=chain_id)
        prism_gap_file = os.path.join(folder.ddG_output, f'prism_rosetta_{prism_nr}_{sys_name}-gap.txt')
        rosetta_to_prism(ddg_gap_file, prism_gap_file, sequence, rosetta_info=None,
                         version=version, sys_name=sys_name, first_residue_number=1, sha_tag=sha_tag, MP=MP, 
                         span_file=span_file, lipid_file=lipid_file, scale=scale)
        create_copy(prism_gap_file, folder.output)

        pdb_gap_file = os.path.join(folder.ddG_output, 'relaxed_gap.pdb')
        shift_pdb_numbering(pdb_file, pdb_gap_file, sec_all, startnr=1)
        create_copy(pdb_gap_file, folder.output, name=f'{sys_name}_final_gap.pdb')

    if zip_files:
        logger.info('Zip files...')
        output_filename = os.path.join(folder.output[:-7], 'calculation')
        shutil.make_archive(output_filename, 'tar', folder.output[:-7])
        for r, d, f in os.walk(folder.output[:-7]):
            if not r in [folder.output, folder.output[:-7]]:
                shutil.rmtree(r, ignore_errors=True)


def runtime_memory_stats(ddG_run_folder):
    #This part collects informations about runtime and memory use
    try:
        job_id_pos= join(ddG_run_folder, 'job_id_ddg.txt')
        with open(job_id_pos, 'r') as job_id_file:
            ddg_process_id=str(job_id_file.readlines()[-1])
            print(ddg_process_id)
         
        #Get stats 
        shell_command = f'sacct -j {relax_process_id} --format="JobID,Start,End,CPUTime,NCPUS,Elapsed,ReqMem,MaxRSS,MaxVMSize,AveVMSize,JobName" > memory_usage_{ddg_process_id}.log'
        subprocess.call(shell_command, cwd=ddG_run_folder, shell=True)
                                                
        check_memory(ddg_process_id,ddG_run_folder) 
        
    except:     
        print('no memory file')
    
