"""make_logs.py controls all the logging in the pipeline.

Author: Anders Frederiksen


Date of last major changes: 2020-05-05

"""
from os.path import join
import logging


class make_log:
    def __init__(self, folder, verbose=False):
        #First we define the different handlers 
        self.path_to_logs = folder.logs
        self.verbose = verbose
        self.logger = logging.getLogger('Pipeline_logger')
        self.logger.setLevel(logging.INFO)

        # This is the formats for the logger
        formatter = logging.Formatter(
            '%(asctime)s - %(levelname)s - %(message)s', datefmt='%Y-%m-%d %H:%M:%S')
        
        #This is the filehandler that handles the info.log
        fh = logging.FileHandler(join(self.path_to_logs, 'info.log'))
        fh.setLevel(logging.INFO)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        
        #This is the stream handler that handles what is displayed in the console
        ch = logging.StreamHandler()
        ch.setLevel(logging.INFO)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
        
        if self.verbose == True:
            #If verbose = True then we make a handler for the debug file
            vh = logging.FileHandler(join(self.path_to_logs, 'debug.log'))
            vh.setLevel(logging.DEBUG)
            vh.setFormatter(formatter)
            self.logger.addHandler(vh)

    def debug(self, string):
        self.logger.debug(string)

    def info(self, string):
        self.logger.info(string)

    def error(self, string):
        self.logger.error(string)

    def warning(self, string):
        self.logger.warning(string)

    def critical(self, string):
        self.logger.critical(string)