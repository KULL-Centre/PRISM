import sys
import os
from helper import create_symlinks, create_copy, find_copy
# this function will parse the results of a rosetta pre_relax run.
# it will be callable from an sbatch script, that can wait for the relaxation to finish,
# and then select the best one.
# the point is to implement the 20x pre relaxation for the stability pipeline, that Amelie requested


def parse_relax_results(relax_run, path_to_run_folder, ddG_input):
    '''This function parses the scorefile from a rosetta
    pre-relaxation, and selects the lowest scoring one'''
    path_to_scorefile = os.path.join(
            relax_run, 'score_bn15_calibrated.sc')
    with open(path_to_scorefile) as scorefile:
        scorelines = scorefile.readlines()

    # Put the relaxation scores in this dict
    relax_scores = {}

    for line in scorelines:
        score_fields = line.split()
        name = score_fields[-1].strip()
        # check if this is an actual scoring line
        if score_fields[0].strip() == 'SCORE:' and name != 'description':
            score = float(score_fields[1])
            relax_scores[name] = score

    # now find the lowest scoring one.
    # set it to the last one, just to start with something
    most_relaxed = name
    for key in relax_scores:
        # it the key is more relaxed than the previous best, update the best.
        if relax_scores[key] < relax_scores[most_relaxed]:
            most_relaxed = key

    print('most relaxed structure is ', most_relaxed)
    print('deleting the rest')
    # and now delete all the structures that are not the best scoring one.
    # it seems a little crude, but whatever.
    for key in relax_scores:
        if key != most_relaxed:
            path_to_tense = os.path.join(relax_run, f'{key}.pdb')
            print('deleting', path_to_tense)
            os.remove(path_to_tense)
    relax_output_strucfile = find_copy(
            relax_run, '.pdb', path_to_run_folder, 'output.pdb')
    
    create_copy(
            os.path.join(path_to_run_folder, 'output.pdb'), ddG_input, name='input.pdb')
    
    return os.path.join(path_to_run_folder, f'{most_relaxed}.pdb')


if __name__ == '__main__':
    parse_relax_results(relax_run=sys.argv[1], path_to_run_folder=sys.argv[2] ,ddG_input=sys.argv[3])
