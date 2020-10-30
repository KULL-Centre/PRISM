"""plotting includes functions for plotting the output of UniprotToPDB.

Author: Anders Frederiksen

Date of last major changes: 2020-10-09

"""

# Standard library imports
import sys
import os
from os import listdir
from os.path import isfile,join
import warnings
warnings.filterwarnings('ignore')

# Third party imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure


def plot_all_seqs(df, uniprot_accession='Protein', save_fig=False, location=False):
    f = plt.figure(figsize=(20, 6))
    sum_list = np.zeros(len(df['binary_sequence']))
    num = 1
    ax = plt.subplot(131)
    temp = []

    for p in range(len(df['binary_sequence'])):
        x = df['binary_sequence'][p]
        if x == 'NaN':
            x = 0
        if num == 1:
            temp = np.zeros(len(df['binary_sequence'][p]))
            num += 1
        if num != 1:
            if x == 'NaN':
                x = 0

            temp = x+temp
            plt.plot(x*num, '|')
            num += 1

    ax.tick_params(labelleft=True)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.suptitle(uniprot_accession, fontsize=20)
    plt.xlabel('Uniprot positions', fontsize=20)
    plt.ylabel('Sequences', fontsize=20)

    ax = plt.subplot(132)
    ax.tick_params(labelleft=True)
    ax.tick_params(axis='both', which='major', labelsize=20)
    plt.suptitle(uniprot_accession, fontsize=20)
    plt.xlabel('Uniprot positions', fontsize=20)
    plt.ylabel('#PDBs', fontsize=20)
    plt.plot(temp)

    #####

    data_temp = df

    temp_seq = np.zeros(len(data_temp['binary_sequence'][0]))
    for n in range(len(data_temp)):
        dat = data_temp['binary_sequence'][n]
        try:
            for m in range(len(dat)):
                if int(dat[m]) == 1:
                    temp_seq[m] = 1
        except:
            continue
    temp_seq
    ax = plt.subplot(133)
    for m in range(len(data_temp['binary_sequence'])):
        try:
            if data_temp['Method'][m] == 'X-ray diffraction':
                seq1 = data_temp['binary_sequence'][m]
                seq2 = temp_seq
                plt.plot(sum(seq1)/sum(seq2),
                         data_temp['resolution'][m], '.', MarkerSize=20)
            elif data_temp['Method'][m] == 'Solution NMR':
                seq1 = data_temp['binary_sequence'][m]
                seq2 = temp_seq
                plt.plot(sum(seq1)/sum(seq2), -1, '.', MarkerSize=20)
        except:
            continue
    plt.suptitle(uniprot_accession, fontsize=20)
    plt.xlabel('Coverage', fontsize=20)
    plt.ylabel('Resolution', fontsize=20)
    ax.axhspan(-1.3, -0.7, facecolor='grey', alpha=0.5)
    locs, labels = plt.yticks()  # Get the current locations and labels.
    # Set text labels and properties.
    plt.yticks([-1, 0, 1, 2, 3, 4], ['NMR', 0, 1, 2, 3, 4], rotation=0)
    if save_fig == True and location != False:
        location = location
        figure_name = f'{uniprot_accession}_overview.png'
        output = join(location, figure_name)
        plt.savefig(output)
        
def res_x_cover(df,uniprot_accession='Protein'):
    data_temp = df

    temp_seq = np.zeros(len(data_temp['binary_sequence'][0]))
    for n in range(len(data_temp)):
        dat = data_temp['binary_sequence'][n]
        try:
            for m in range(len(dat)):
                if int(dat[m]) == 1:
                    temp_seq[m] = 1
        except:
            continue
    temp_seq
    fig, ax = plt.subplots()
    for m in range(len(data_temp['binary_sequence'])):
        try:
            if data_temp['Method'][m] == 'X-ray diffraction':
                seq1 = data_temp['binary_sequence'][m]
                seq2 = temp_seq
                plt.plot(sum(seq1)/sum(seq2),
                         data_temp['resolution'][m], '.', MarkerSize=20)
            elif data_temp['Method'][m] == 'Solution NMR':
                seq1 = data_temp['binary_sequence'][m]
                seq2 = temp_seq
                plt.plot(sum(seq1)/sum(seq2), -1, '.', MarkerSize=20)
        except:
            continue
    plt.suptitle(uniprot_accession, fontsize=20)
    plt.xlabel('Coverage', fontsize=20)
    plt.ylabel('Resolution', fontsize=20)
    ax.axhspan(-1.3, -0.7, facecolor='grey', alpha=0.5)
    locs, labels = plt.yticks()  # Get the current locations and labels.
    # Set text labels and properties.
    plt.yticks([-1, 0, 1, 2, 3], ['NMR', 0, 1, 2, 3], rotation=0)
    
    
    
    
def plot_pdb_seq(df,uni,pdbid_list):
    f, ax = plt.subplots(figsize=(10, 2.5))  
    colors=['limegreen','goldenrod','purple']
    numX=0
    for pdbid in pdbid_list:
        numX+=1
        for n in range(len(df['pdb'])):
            num=0
            num2=0
    
            for m in df['binary_sequence'][n]:
                if m == 1 and num2 ==0:
                    start=num
                    num2+=1                
                if m == 1:
                    end=num                
                num+=1
            pdb=df['pdb'][n]
    
        
            if pdb==pdbid:
                        
                if len(pdbid_list) > 3:          
                    f=plt.plot(df['binary_sequence'][n]*numX,'|',MarkerSize=14, solid_capstyle='butt')            
                else:
                    f=plt.plot(df['binary_sequence'][n]*numX,'|',MarkerSize=14, solid_capstyle='butt', color=colors[numX-1]) 
                ax.tick_params(labelleft=True)
                ax.tick_params(axis='both', which='major', labelsize=20)
                ax.tick_params(axis='y', which='major', labelsize=20,grid_alpha=0.5)
                ax.spines['left'].set_visible(False)
                ax.spines['top'].set_visible(False)
                ax.spines['bottom'].set_visible(False)
                ax.spines['right'].set_visible(False)
                #plt.tick_params(axis='y', which='major', labelsize=__)
                plt.suptitle(uni, fontsize=20)
                plt.xlabel('Uniprot positions', fontsize=20)
                #plt.ylabel('#PDBs', fontsize=20)
    num_var=0
    arr1=[0]
    arr2=[ 'Not covered']
    for var in pdbid_list: 
        num_var+=1
        arr1.append(num_var)      
        arr2.append(var)
    plt.yticks(arr1[0:(len(pdbid_list)+1)], arr2[0:(len(pdbid_list)+1)], rotation=0)
    plt.ylim((0.5, numX+0.5))
    
    
    
def calc_overlaps(df,uniprot_accession):
    df_coords=df
    add=np.zeros(len(df_coords['binary_sequence'][0]))
    ns=[]
    empty_arr=np.zeros(len(df_coords['binary_sequence'][0]))
    arr_dic={}
    
    for n in range(len(df_coords)):   
        add=add+df_coords['binary_sequence'][n]
        empty_arr=np.zeros(len(df_coords['binary_sequence'][0]))
        for m in range(len(df_coords)):
            namen=df_coords['pdb'][n]
            namem=df_coords['pdb'][m]
            if n != m:
                empty_arr=df_coords['binary_sequence'][n]+df_coords['binary_sequence'][m]
                arr_dic[f'{namen}_{namem}']=[empty_arr,(df_coords['resolution'][n]+df_coords['resolution'][m])/2]          
            if n == m:
                empty_arr=df_coords['binary_sequence'][n]
                arr_dic[f'{namen}_{namem}']=[empty_arr,(df_coords['resolution'][n]+df_coords['resolution'][m])/2]
    
    #########

    plottings_single_x=[]
    plottings_single_y=[]
    plottings_resgood_x=[]
    plottings_resgood_y=[]
    plottings_resbad_x=[]
    plottings_resbad_y=[]
    
    for p in arr_dic:
        q=p.split('_')
        
        if q[0]!=q[1]:
            zeros=list(arr_dic[p][0]).count(0)
            ones=list(arr_dic[p][0]).count(1)
            twos=list(arr_dic[p][0]).count(2)
            cover=len(arr_dic[p][0])-zeros
            
            if arr_dic[p][1] < 2.5:
                plottings_resgood_x=np.hstack((plottings_resgood_x,twos))
                plottings_resgood_y=np.hstack((plottings_resgood_y,cover))
            else:
                plottings_resbad_x=np.hstack((plottings_resbad_x,twos))
                plottings_resbad_y=np.hstack((plottings_resbad_y,cover))
        if q[0]==q[1]:
            zeros=list(arr_dic[p][0]).count(0)
            ones=list(arr_dic[p][0]).count(1)
            twos=list(arr_dic[p][0]).count(2)
            cover=len(arr_dic[p][0])-zeros
            plottings_single_x=np.hstack((plottings_single_x,-30))
            plottings_single_y=np.hstack((plottings_single_y,cover))

    plottings_single=[plottings_single_x,plottings_single_y]
    plottings_resbad=[plottings_resbad_x,plottings_resbad_y]        
    plottings_resgood=[plottings_resgood_x,plottings_resgood_y] 
    plot_dic={'add':add,'plottings_resgood':plottings_resgood,'plottings_resbad':plottings_resbad,'plottings_single':plottings_single}

    return(arr_dic,plot_dic)
 
    
def plot_overlaps(arr_dic,plot_dic,uniprot_accession):
    plottings_single=plot_dic['plottings_single']
    plottings_resbad=plot_dic['plottings_resbad']        
    plottings_resgood=plot_dic['plottings_resgood']  
    add=plot_dic['add']
    
    zeros_cover=list(add).count(0)
    coverage=len(add)-zeros_cover    
    
    f = plt.figure(figsize=(6, 6))  
    f = plt.plot(coverage,coverage,'*',Markersize=15)       
    f = plt.plot(plottings_resgood[0],plottings_resgood[1],'v',Markersize=12) 
    f = plt.plot(plottings_single[0],plottings_single[1],'X',Markersize=10)
    f = plt.plot(plottings_resbad[0],plottings_resbad[1],'o',Markersize=10,alpha=0.5)
                
    plt.ylim(-40,len(add)+10)
    plt.xlim(-40,len(add)+10) 
    plt.suptitle(f'Overlaps x Coverage \n {uniprot_accession}', fontsize=20)
    plt.xlabel('Multi-PDBs', fontsize=20)
    plt.ylabel('Coverage', fontsize=20)

    
    plt.axhspan((coverage/100) * 95, (coverage/100) * 100, facecolor='darkgreen', alpha=0.5)
    plt.axhspan((coverage/100) * 90, (coverage/100) * 95, facecolor='darkgreen', alpha=0.4)
    plt.axhspan((coverage/100) * 85, (coverage/100) * 90, facecolor='darkgreen', alpha=0.3)
    plt.axhspan((coverage/100) * 80, (coverage/100) * 85, facecolor='darkgreen', alpha=0.2)
    plt.axhspan((coverage/100) * 75, (coverage/100) * 80, facecolor='darkgreen', alpha=0.1)
    plt.axhspan((coverage/100) * 70, (coverage/100) * 75, facecolor='turquoise', alpha=0.1)
    plt.axhspan((coverage/100) * 65, (coverage/100) * 70, facecolor='turquoise', alpha=0.2)
    plt.axhspan((coverage/100) * 60, (coverage/100) * 65, facecolor='turquoise', alpha=0.3)
    plt.axhspan((coverage/100) * 55, (coverage/100) * 60, facecolor='turquoise', alpha=0.4)
    plt.axhspan((coverage/100) * 50, (coverage/100) * 55, facecolor='turquoise', alpha=0.5)
    plt.axhspan((coverage/100) * 45, (coverage/100) * 50, facecolor='orange', alpha=0.5)
    plt.axhspan((coverage/100) * 40, (coverage/100) * 45, facecolor='orange', alpha=0.4)
    plt.axhspan((coverage/100) * 35, (coverage/100) * 40, facecolor='orange', alpha=0.3)
    plt.axhspan((coverage/100) * 30, (coverage/100) * 35, facecolor='orange', alpha=0.2)
    plt.axhspan((coverage/100) * 25, (coverage/100) * 30, facecolor='orange', alpha=0.1)
    plt.axhspan((coverage/100) * 20, (coverage/100) * 25, facecolor='red', alpha=0.1)
    plt.axhspan((coverage/100) * 15, (coverage/100) * 20, facecolor='red', alpha=0.2)
    plt.axhspan((coverage/100) * 10, (coverage/100) * 15, facecolor='red', alpha=0.3)
    plt.axhspan((coverage/100) * 5,  (coverage/100) * 10, facecolor='red', alpha=0.4)
    plt.axhspan((coverage/100) * 0,  (coverage/100) * 5, facecolor='red', alpha=0.5)
    plt.axhspan((coverage/100) * 0,  (coverage/100) * -15, facecolor='grey', alpha=0.5)
    plt.axhspan((coverage/100) * 100, (len(add)+10), facecolor='grey', alpha=0.5)
    plt.axhspan((coverage/100) * 100, len(add), facecolor='purple', alpha=0.2)  
    
    f = plt.plot([], [], color='black', marker='o',markersize=15, label='Overlaps')
    f = plt.plot([], [], color='black', marker='X',markersize=15, label='Single-PDB')
    f = plt.plot([], [], color='black', marker='*',markersize=15, label='Max Cover')
    f = plt.legend( loc='lower right')

    
def findbest(df,uniprot_accession):
    arr_dic,plot_dic=calc_overlaps(df,uniprot_accession)
    add=plot_dic['add']
    zeros_add=list(add).count(0)
    coverage=len(add)-zeros_add
    good_cover=[]
    for n in arr_dic:       
        zeros=list(arr_dic[n][0]).count(0)
        ones=list(arr_dic[n][0]).count(1)
        twos=list(arr_dic[n][0]).count(2)
        cover=len(add)-zeros    
        
        if cover >((coverage/100) * 95):
            good_cover.append([n,twos,cover,arr_dic[n][1]])
        else:
            continue  

    if good_cover==[]:
        for n in arr_dic:       
            zeros=list(arr_dic[n][0]).count(0)
            ones=list(arr_dic[n][0]).count(1)
            twos=list(arr_dic[n][0]).count(2)
            cover=len(add)-zeros    
            
            if cover >((coverage/100) * 85):
                good_cover.append([n,twos,cover,arr_dic[n][1]])
            else:
                continue 
    res=10
    final=[]
    for n in good_cover:   
        try:
            if float(n[3]) < res:
                res=n[3]
                final=n
        except:
            continue 
    if final==[]:
        for n in good_cover:   
            try:
                res=n[3]
                final=n
            except:
                continue  
            
    #print('Choose_pair:  ',final[0],'\nCoverage of possible: ',f'{round(final[2]/coverage*100,1)}%','\nCoverage of uniprot: ',f'{round(final[2]/len(add)*100,1)}%','\nOverlaps: ',final[1],'\nResolution: ',final[3])
    return(final,coverage,add)

    
    