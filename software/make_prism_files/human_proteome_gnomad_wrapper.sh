#!/bin/bash

while read name ID seq
do
	#echo $name
	echo $ID
	#echo $seq
	 
	outfile=$(python3 /storage1/hezscha/src/gnomad_to_prism/scripts/make_prism_gnomad_file_seq.py -e mis -uniprot $ID -seq $seq)
	#outfile=$(python3 /storage1/hezscha/genome_proteome_map/scripts/parse_human_proteome_step2_v01.py -e mis -uniprot $ID -seq $seq -out_folder /storage1/hezscha/genome_proteome_map/results/test)
	
	echo "$name $ID $outfile" >> /storage1/hezscha/genome_proteome_map/data/UP000005640_9606_gnomad_files_Mar2021.list
	
done < /storage1/hezscha/genome_proteome_map/data/UP000005640_9606_ID_seq.list


