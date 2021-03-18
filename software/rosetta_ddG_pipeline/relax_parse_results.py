"""relax_parse_results.py 
This function will parse the results of a rosetta pre_relax run. it will be callable from an 
sbatch script, that can wait for the relaxation to finish, and then select the best one. The 
point is to implement the 20x pre relaxation for the stability pipeline, that Amelie requested.

Author: Anders Frederiksen

Date of last major changes: 2020-04

"""

# Standard library imports
import logging as logger
import os
import sys
from os.path import join
import subprocess

# Local application imports
from helper import AttrDict, create_symlinks, create_copy, find_copy
from get_memory_stats import check_memory

def parse_relax_results(folder, sc_name='score_bn15_calibrated', logger_mode='info', keep=1):
    '''This function parses the scorefile from a rosetta
    pre-relaxation, and selects the lowest scoring one'''
    print(sc_name, keep, logger_mode, folder)
    try:
        job_id_pos= join(folder.relax_run, 'job_id_relax.txt')
        with open(job_id_pos, 'r') as job_id_file:
            relax_process_id=str(job_id_file.readlines()[-1])
            print(relax_process_id)
         
        #Get stats 
        shell_command = f'sacct --format="JobID,Start,End,CPUTime,ReqMem,MaxRSS,MaxVMSize,AveVMSize,JobName" > memory_usage_{relax_process_id}.log'
        subprocess.call(shell_command, cwd=folder.relax_run, shell=True)
                                                
        check_memory(relax_process_id,folder.relax_run) 
        
    except: 
        print('no memory file')
    
    
    relax_scores = {}
    relax_models = {}
    for n in range(20):
        if sc_name=='score_bn15_calibrated':
            path_to_scorefile = os.path.join(
                folder.relax_run, f'{n}-{sc_name}.sc')
        else:
            path_to_scorefile = os.path.join(
                folder.relax_run, f'{sc_name}.sc')
        with open(path_to_scorefile) as scorefile:
            scorelines = scorefile.readlines()

        # Put the relaxation scores in this dict
        for line in scorelines:
            score_fields = line.split()
            name = score_fields[-1].strip()
            if score_fields[0].strip() == 'SCORE:' and name != 'description':
                score = float(score_fields[1])
                relax_scores[score] = name
                relax_models[name] = score
    print(relax_scores, keep)
    most_relaxed_array = sorted(list(relax_scores.keys()), reverse=False)[:int(keep)]
    most_relaxed = relax_scores[most_relaxed_array[0]]
    print(most_relaxed_array)
    most_relaxed_models = [relax_scores[score] for score in most_relaxed_array]
    logger.info(f'most relaxed structure is {most_relaxed}.')
    for key in relax_models:
        if not key in most_relaxed_models:
            path_to_tense = os.path.join(folder.relax_run, f'{key}.pdb')
            print('deleting', path_to_tense)
            os.remove(path_to_tense)
    if int(keep) > 1:
        for ind, ddg_subfolder in enumerate(folder.ddG_input):
            relax_output_strucfile = find_copy(
                folder.relax_run, f'{relax_scores[most_relaxed_array[ind]]}.pdb', folder.relax_output, f'output_{ind}.pdb')
            create_copy(
                os.path.join(folder.relax_output, f'output_{ind}.pdb'), ddg_subfolder, name=f'input.pdb')
    else:
        relax_output_strucfile = find_copy(
            folder.relax_run, f'{most_relaxed}.pdb', folder.relax_output, 'output.pdb')
        create_copy(
            os.path.join(folder.relax_output, 'output.pdb'), folder.ddG_input, name='input.pdb')

    return os.path.join(folder.relax_output, f'{most_relaxed}.pdb')


if __name__ == '__main__':
    folder = AttrDict()
    print(sys.argv)
    if len(sys.argv) <= 4:
        folder.update({'relax_run': sys.argv[1], 'relax_output': sys.argv[2], 
            'ddG_input': sys.argv[3]})
    else:
        folder.update({'relax_run': sys.argv[1], 'relax_output': sys.argv[2], 
            'ddG_input': [sys.argv[x] for x in range(3, len(sys.argv)-2)]})
    print(sys.argv)
    print(folder)
    if len(sys.argv) == 5:
        parse_relax_results(folder, sc_name=sys.argv[4])
    elif len(sys.argv) > 5:
        parse_relax_results(folder, sc_name=sys.argv[-2], keep=int(sys.argv[-1]))
    else:
        parse_relax_results(folder)
