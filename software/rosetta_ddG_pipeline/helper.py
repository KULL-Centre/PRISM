"""helper.py contains helpful functions.

Author: Anders Frederiksen
Contribution: Johanna K.S. Tiemann

Date of last major changes: 2020-04-15

"""

# Standard library imports
import logging as logger
import os
import shutil


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
