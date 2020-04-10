import sys
import os
import logging as logger
from args_pipeline import parse_args2


def check_paths(path):
    args = parse_args2()
    overwrite_path = args.OVERWRITE_PATH
    try:
        os.makedirs(path)
        logger.info(f"Global working directory ({path}) created.")
    except FileExistsError:
        if overwrite_path:
            logger.warn(f"Directory {path} already exists. overwrite_path set to {overwrite_path}, so we will use this path.")
        else:
            logger.error(f"Directory {path} already exists. overwrite_path set to {overwrite_path}, so we stop the execution. Please provide a different output directory.")
            sys.exit() #will terminate the complete script
    return(path)        

class folder2:
    
    def __init__(self,output_path):
        ##Create
        #Main folders
        if output_path[-1] == '/':
            output_path = output_path[:-1]
        else:
            output_path = output_path
            
        self.output_path = output_path+'/'
        self.input = check_paths(self.output_path + 'input/') 
        self.prepare = check_paths(self.output_path + 'prepare/') 
        self.relax = check_paths(self.output_path + 'relax/') 
        self.ddG = check_paths(self.output_path + 'ddG/') 
        self.output = check_paths(self.output_path + 'output/') 
        self.analysis = check_paths(self.output_path + 'analysis/') 
        
        #Subfolders
        self.input_mutfiles = check_paths(self.input +'mutfiles/')
        self.input_cleaning = check_paths(self.input +'cleaning/')
        self.input_checking = check_paths(self.input +'checking/')
        self.input_mp_files = check_paths(self.input +'mp_files/')
        self.input_output = check_paths(self.input +'output/')
        
        self.relax_input = check_paths(self.relax +'input/')
        self.relax_run = check_paths(self.relax +'run/')
        self.relax_output = check_paths(self.relax +'output/')        
    
        self.ddG_input = check_paths(self.ddG +'input/')
        self.ddG_run = check_paths(self.ddG +'run/')
        self.ddG_output = check_paths(self.ddG +'output/') 
        
        
        return;