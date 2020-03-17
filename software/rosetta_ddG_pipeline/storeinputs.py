

def storeinputfuc(name,chain_id,structure,out_path,uniprot_accesion,relax_flag_file,ddg_flag_file,path_to_input):
    
    with open(path_to_input+'inputs','w') as file:
        file.write('NAME ='+name +'\n' + 'CHAIN =' +chain_id +'\n' +'STRUCTURE ='+structure +'\n'+'OUTPATH ='+out_path+'\n'+uniprot_accesion+'\n')
    file.close()
    open(path_to_input+"relax_flag_file", "w").writelines(open(relax_flag_file).readlines())     
    open(path_to_input+"ddg_flag_file", "w").writelines(open(ddg_flag_file).readlines()) 
    print("Input file can be found at",path_to_input) 
    return