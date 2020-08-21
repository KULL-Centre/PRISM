import json
import sys
from os.path import join
import subprocess
from parse_cartesian_functions import rosetta_cartesian_read, ddgs_from_dg


def parse_rosetta_ddgs(sys_name, chain_id, fasta_seq, ddG_input, ddG_output):

    path_to_run_folder = ddG_input
    print('the path to run folder is')
    print(path_to_run_folder)

    rosetta_summary_file = f'{sys_name}_{chain_id}.rosetta_cartesian.dgs'
    shell_command = f'cat *.ddg | grep -v WT > {rosetta_summary_file}'
    print('calling to the shell:')
    print(shell_command)
    subprocess.call(shell_command, cwd=path_to_run_folder, shell=True)

    rosetta_cartesian_ddgs_dict = ddgs_from_dg(rosetta_cartesian_read(
        join(path_to_run_folder, rosetta_summary_file), fasta_seq))
    line = []
    list_keys = list(rosetta_cartesian_ddgs_dict.keys())
    uniprot_numbering_ddgs_dict = {}
    scorefile = open(join(ddG_output, f'{sys_name}_ddg.txt'), 'w')
    scorefile.close
    scorefile = open(join(ddG_output, f'{sys_name}_ddg.txt'), 'w')

    for key in rosetta_cartesian_ddgs_dict:

        position = int(key[1:-1])
        uniprot_position = position
        uniprot_numbering_ddgs_dict[
            key[0] + str(uniprot_position) + key[-1]] = rosetta_cartesian_ddgs_dict[key]

    scorefile.write(f'#Rosetta cartesian_ddg stability predictions for {sys_name}\n')
    scorefile.write(f'#sequence is {fasta_seq}\n')
    scorefile.write(
        'UAC_pos\t A \t C \t D \t E \t F \t G \t H \t I \t K \t L \t M \t N \t P \t Q \t R \t S \t T \t V \t W \t Y \n')
    scorefile_line = '{}' + '\t {:.3}' * 20 + '\n'

    for i in range(1, len(fasta_seq.strip()) + 1):
        refs = 'ACDEFGHIKLMNPQRSTVWY'

        for AA in refs:
            a = str(fasta_seq[i - 1]) + str(i) + str(AA)
            if a in str(list_keys):

                line.append(rosetta_cartesian_ddgs_dict[str(a)])
                
            else:
                line.append('-')
        scorefile.write(scorefile_line.format(i, line[0 + (i - 1) * 20], line[1 + (i - 1) * 20], line[2 + (i - 1) * 20], line[3 + (i - 1) * 20], line[4 + (i - 1) * 20], line[5 + (i - 1) * 20], line[6 + (i - 1) * 20], line[7 + (i - 1) * 20], line[8 + (i - 1) * 20], line[9 + (
            i - 1) * 20], line[10 + (i - 1) * 20], line[11 + (i - 1) * 20], line[12 + (i - 1) * 20], line[13 + (i - 1) * 20], line[14 + (i - 1) * 20], line[15 + (i - 1) * 20], line[16 + (i - 1) * 20], line[17 + (i - 1) * 20], line[18 + (i - 1) * 20], line[19 + (i - 1) * 20]))

    scorefile.close
    line = []


if __name__ == '__main__':
    parse_rosetta_ddgs(sys_name=sys.argv[1], chain_id=sys.argv[
                       2], fasta_seq=sys.argv[3], ddG_input=sys.argv[4], ddG_output=sys.argv[5])
