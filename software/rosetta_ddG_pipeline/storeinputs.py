from os.path import join, basename
from shutil import copyfile


def storeinputfuc(name, chain_id, structure, out_path, uniprot_accesion, relax_flag_file, ddg_flag_file, path_to_input, mutations=None):

    with open(join(path_to_input, 'inputs'), 'w') as file:
        file.write(
            (f'NAME = {name}\n'
             f'CHAIN = {chain_id}\n'
             f'STRUCTURE = {structure}\n'
             f'OUTPATH = {out_path}\n'
             f'{uniprot_accesion}\n')
        )
    copyfile(relax_flag_file, join(path_to_input, basename(relax_flag_file)))
    copyfile(ddg_flag_file, join(path_to_input, basename(ddg_flag_file)))
    print('Input file can be found at', path_to_input)
    # Copy structure file
    copyfile(structure, join(path_to_input, basename(structure)))
    # Copy mutation file if present
    if mutations:
        copyfile(mutations, join(path_to_input, 'mutations'))
    return
