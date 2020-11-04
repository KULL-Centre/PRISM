
import sys
import os
from os import listdir
from os.path import isfile, join
import pandas as pd
import subprocess
import numpy as np

def check_memory(job_id,path):

    memory_file=f'memory_usage_{job_id}.log'
    memory=join(path,memory_file)   
    memory = pd.read_fwf(memory)   
    hours=0;minutes=0;seconds=0;jobs=0;max_memo=[];max_time=0
    
    for n in range(len(memory)):
        jobname=memory['JobID'][n]
        job=jobname.split('_')
        if str(job[0]) == str(job_id)  and memory['JobName'][n] != 'batch':
            cpu=str(memory['CPUTime'][n])            
            time=cpu.split(':')
            time_split=time[0].split('-')

            if len(time) > 0 and len(time_split) < 2:
                hour=int(time[0]);minute=int(time[1]);second=int(time[2])
            if len(time_split) > 1:
                
                hour=int(time_split[1])+int(time_split[0])*int(24);minute=int(time[1]);second=int(time[2]) 
            time=hour*60+minute
            if time > max_time:
                max_time=time
            if minute > 1:
                jobs+=1;hours+=hour;minutes+=minute;seconds+=second
        if str(job[0]) == str(job_id) and memory['JobName'][n] == 'batch':       
            maxmem=str(memory['MaxVMSize'][n])
            
            maxmems=float(maxmem[:-1])
            
            #if hour > 1:
            max_memo=np.hstack((max_memo,maxmems/1000000))
          
    total_time=hours + minutes/60 + seconds/3600
    total_time_per_residue=total_time/jobs
    
    print('Total_time=',total_time,'   Time_per_residue=',total_time_per_residue,'   Time_per_residue_min=',total_time_per_residue*60,'   Max_time=',max_time/60,'   #jobs=',jobs,'   Max_memory=',np.max(max_memo))
        
    
    with open(join(path,f'memory_stats_{job_id}.txt'), 'w') as memory_data:       
  
        memory_data.write('# ------------------------------------\n')
        memory_data.write(f'# Job_ID = {job_id}\n')
        memory_data.write(f'#\t Total_time = {total_time}\n')
        memory_data.write(f'#\t Time_per_job = {total_time_per_residue}\n')
        memory_data.write(f'#\t Time_per_job_minutes = {total_time_per_residue*60}\n')
        memory_data.write(f'#\t Max_time = {max_time/60}\n')
        memory_data.write(f'#\t Max_memory = {np.max(max_memo)}\n')
        memory_data.write(f'#\t #Jobs = {jobs}\n')
        
    
    
if __name__ == '__main__':
    job_id = sys.argv[1]
    path = sys.argv[2]
    check_memory(job_id,path)