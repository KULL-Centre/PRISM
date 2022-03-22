
#v04: rolled back disqualify, kept nhomalt
#v03: Extract also the number of homozygots: nhomalt=
#rolled back: v03: disqualify variants called in less than half the samples, i.e. where AN is less than the number of samples (since humans are diploid)
#v02: change the way I write the files. Header is overwrite mode, var is append mode
#v02: only write an index if extraction mode is 'all' or it was specifically requested
#parse the vcf file for the requested chromosome. Put exome and genome data in separate files and also separate by var type in case we change our mind about which complex variants we want to look at. So in total produce 6 tmp files per transcript

import argparse
import sys
import cyvcf2
import re
import os

from vcf_parser_func_v07_revep import parse_HGVSp

#parse arguments
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-chrom', dest="chromosome", type = str, required=True, help="The chromosome for which to parse the gnomad_v2 exomes vcf and gnomad_v3 genomes vcf")
parser.add_argument('-filter', dest="FILTER", default='sg,sl,idel,iins', help='The complex variant types to extract, given as a comma separeted list (no whitespace!). Missense and synonymous variants are always extracted (unless omitted in -e option). The currently available types are: stop gained : sg, start lost : sl, inframe deletion : idel, inframe insertion : iins. Default is "sg,sl,idel,iins"')
parser.add_argument('-e', dest="extract", choices=['mis', 'syn', 'complex', 'all'], default = 'mis', help="Type of variants to write files for. Choose one of the list. The current default is 'mis'. You probably only want missense ('mis') and synonymous ('syn') since the notation for complex variants is still experimental.")
parser.add_argument('-tmp_folder', dest='tmp_folder', help="Where to store the output files. Default location is /storage1/hezscha/gnomad_to_prism_parser/step1_files/+args.chromosome/")
parser.add_argument('-v', dest="v", action='store_true', help="Verbose")
parser.add_argument('-index', dest="index", action='store_true', help="Write an index of prefixes into the tmp dir. Automatically on if -e is 'all'. Default is off (as for all arguments with the action store_true)")
parser.add_argument('-d', dest="data_source", choices=['exomes', 'genomes', 'both'], default = 'both', help="Choose whether to extract data from 'exomes', 'genomes' or 'both'(default).")
args = parser.parse_args()
################################################################################

#step 1: setup

#always write an index if all variant types are extracted. Otherwise, only if requested. The index is used by step2 to know which files it should attempt to open because I want to open corresponding exome and genome files directly after each other and if I just go through all file in dir I'm not sure it's guaranteed they are after each other. I also don't want to go through either .exomes or .genomes and open the other one because either of these sets might have variants for a transcript the other one doesn't have variants in (more individuals in exomes, but deep seq in genomes). And I don't want to write double prism files. So an index of the existing transcripts seems a good way. 
if args.extract == 'all':
	args.index = True

filter_d = {'sg':'stop_gained','sl':'start_lost','idel':'inframe_deletion', 'iins':'inframe_insertion'}

#which var types are we looking to extract?
if args.extract == 'all':
	filter_type=set()
	for item in args.FILTER.split(','):
		if item in filter_d:
			filter_type.add(filter_d[item])
		else:
			print('Variant type "', item, 'not recognized!', file = sys.stderr)
			sys.exit()
	#add missense and synonymous to filter_type
	filter_type.add('missense_variant')
	filter_type.add('synonymous_variant')
elif args.extract == 'mis':
	filter_type=set()
	filter_type.add('missense_variant')
elif args.extract == 'syn':
	filter_type=set()	
	filter_type.add('synonymous_variant')	
else:
	filter_type=set()
	for item in args.FILTER.split(','):
		if item in filter_d:
			filter_type.add(filter_d[item])
		else:
			print('Variant type "', item, 'not recognized!', file = sys.stderr)
			sys.exit()

#where to write the tmp files
tmp_folder = args.tmp_folder if args.tmp_folder else '/storage1/hezscha/gnomad_to_prism_parser/step1_files/'+args.chromosome+'/'

#if we are using the default tmp folder make sure it exists
if not args.tmp_folder:
	if not os.path.exists('/storage1/hezscha/gnomad_to_prism_parser/step1_files/'+args.chromosome+'/'):
		os.makedirs('/storage1/hezscha/gnomad_to_prism_parser/step1_files/'+args.chromosome+'/')
#if we are using the custom tmp folder that doesn't exist complain and exit		
else:
	if not os.path.exists(args.tmp_folder):
		print('Specified folder for temporary files', args.tmp_folder, 'does not exist. Exiting now.')
		sys.exit()

#print some general infos (can be saved into a log file f.x.)
print('Chromosome is', args.chromosome)
print('Parsed tmp files go to', tmp_folder)

if args.data_source == 'both':
	print('Parsing exome and genome vcfs')
elif args.data_source == 'exomes':
	print('Parsing exome vcf')
else:
	print('Parsing genome vcf')	

#v04: a file that collects instances of high allel frequencies so I can inspect them more easily
alertOUT = open(tmp_folder+'/high_AF.txt' ,'w')

#step 2: open exomes vcf file
#this is a list of files we have opened before. Used to know if I should write a header to the tmp file. Since the index writer attempts to merge this with the seen genome vars both of them need to exist regardless of whether we have chosen to parse exomes, genomes or both (but can be empty)
already_seen = set()
if args.data_source == 'both' or args.data_source == 'exomes':
	gnomad_version = '2.1.1'

	#2.1. open vcf file
	vcffile = '/storage1/shared/data/gnomAD/gnomad.exomes.r'+gnomad_version+'.sites.'+args.chromosome+'.liftover_grch38.vcf.VEP.cache100.bgz'
	print('Parsing exome vcf', vcffile, '...')
	try:
		vcf_reader = cyvcf2.VCF(vcffile)
	except:
		raise
		sys.exit()

	#2.2. find vep annotation: Which fields exist and what are their positions?
	for line in vcf_reader.raw_header.split('\n'):
		if line.startswith('##INFO=<ID=CSQ'):
			#get this part:"... Format: Allele|..."
			fields = line.split('Format: ')[1].split('|')
			pos = {} #a dict storing the position of each field we search for, or False if the field is not available
			#find out where certain fields are or if they are not available
			#the fields I'm looking for:
			target_fields = ['Gene', 'SYMBOL', 'BIOTYPE', 'CANONICAL', 'Consequence', 'HGVSp', 'Codons', 'SWISSPROT', 'DOMAINS', 'Feature', 'Protein_position', 'Amino_acids']
			#the information I want printed per variant. F.x. if a transcript is canon does not change and neither does the uniprot ID so no reason to print them for every var 
			print_fields = ['Consequence', 'HGVSp', 'Codons', 'SWISSPROT', 'DOMAINS', 'Protein_position', 'Amino_acids']
			for target_field in target_fields:
				pos[target_field] = False
				for (position, field) in enumerate(fields):
					if field == target_field: 
						pos[target_field] = position
						break
			print('Found vep field information')
			#debug
			if args.v:
				print('Field positions:')
				for target_field in pos:
					print(target_field, pos[target_field]) 
			
			#debug
			break

		else:
			continue
	else:
		#if we didn't break we didn't find the vep line
		print('Did not find annotation metadata line ##INFO=<ID=CSQ, please check vcf file')
		sys.exit()

	#2.3. parse var lines
	for record in vcf_reader:

		#first of, skip SNPs that did not PASS quality control. In vcf reader, PASS evaluates to nothing, so if there is something in the filter you want to skip the entry
		if record.FILTER:
			continue

		# ~ #v03: disqualify variants at sites with coverage < Nr of samples which roughly corresponds to sites called in less than half the samples since there are two allels per site.
		# ~ if record.INFO.get('AN') < nr_samples:
			# ~ continue

		#in this vcf readers, the entire vep entry is one string so split it on ','
		#for re-annotation w online vep, change the tag to 'CSQ'. For orig vep annotation in gnomAD files use 'vep' 
		for transcript in record.INFO.get('CSQ').split(','): 
			fields = transcript.split('|')
			
			#we are only interested in coding transcripts and the variant types requested by argv (filter set is adjusted accordingly in step 1)
			if fields[pos['BIOTYPE']] == 'protein_coding' and set(fields[pos['Consequence']].split('&')).intersection(filter_type):

				tID = fields[pos['Feature']]
				gID = fields[pos['Gene']]
				uniprot = fields[pos['SWISSPROT']] if (pos['SWISSPROT'] and fields[pos['SWISSPROT']] ) else 'NA'
				#define the type of variant as missense, syn or complex
				if 'missense_variant' in fields[pos['Consequence']].split('&'):
					vartype = 'mis'
				elif 'synonymous_variant' in fields[pos['Consequence']].split('&'):
					vartype = 'syn'
				else:
					vartype = 'complex'
					
				pref = tID+'_'+gID+'_'+uniprot+'.'+vartype
				#the first time we see this combination of transcript and vartype, write a header. Open file in 'w' mode instead of 'a', this is will also overwrite any existing files.
				if not pref in already_seen:
					with open(tmp_folder+'/'+pref +'.exomes.tmp', 'w') as OUT:
						#meta data for this transcript:
						print('#gnomad_version:', gnomad_version, file = OUT)
						print('#transcript_ID:', tID, file = OUT)
						print('#gene_ID:', gID, file = OUT)
						print('#symbol:', fields[pos['SYMBOL']], file = OUT)
						print('#uniprot_ID:', uniprot, file = OUT)
						if fields[pos['CANONICAL']] == '1':
							print('#canonical: 1', file = OUT)
						else:
							print('#canonical: 0', file = OUT)
						print('#extracted_variants:', ','.join(filter_type), file = OUT)	
						
						print('#prism_var AF AC AN Homozygotes basepos', ' '.join(print_fields), file = OUT)
						already_seen.add(pref)
					
				#write a line for the variant, no matter if we have written a header or not and this is always in append mode (if we didn't write a header it simply means it exists already)
				with open(tmp_folder+'/'+pref +'.exomes.tmp', 'a') as OUT:
					var_name = parse_HGVSp(transcript, HGSV_field_pos=pos['HGVSp'], prot_pos_field_pos=pos['Protein_position'], aa_field_pos=pos['Amino_acids'], v=args.v)
					print(var_name, record.INFO.get('AF'), record.INFO.get('AC'), record.INFO.get('AN'), record.INFO.get('nhomalt'), record.POS, sep = ' ', end = '', file = OUT)
					#v04: print an alert if the allel freq is higher than 0.5
					if float(record.INFO.get('AF')) > 0.5:
						print(pref, record.POS, var_name ,record.INFO.get('AF'), record.INFO.get('AN'), file = alertOUT)
					
					for field in print_fields:
						if (pos[field] and fields[pos[field]]):
							print(' ', fields[pos[field]], sep = '', end = '', file = OUT)
						else:
							print(' NA', end = '', file = OUT)
					#lastly print a newline
					print('\n', end = '', file = OUT)

############################################################################
#step 3: open genomes vcf file
#this is a list of files we have opened before. Used to know if I should write a header to the tmp file.
already_seen_g = set()
if args.data_source == 'both' or args.data_source == 'genomes':
	gnomad_version = '3.0'

	#do what we did above
	#3.1 open vcf
	vcffile = '/storage1/shared/data/gnomAD/gnomad.genomes.r'+gnomad_version+'.sites.chr'+args.chromosome+'.vcf.VEP.cache100.bgz'
	print('Parsing genome vcf', vcffile, '...')
	try:
		vcf_reader = cyvcf2.VCF(vcffile)
	except:
		raise
		sys.exit()

	#3.2. find vep annotation: Which fields exist and what are their positions?
	for line in vcf_reader.raw_header.split('\n'):
		if line.startswith('##INFO=<ID=CSQ'):
			#get this part:"... Format: Allele|..."
			fields = line.split('Format: ')[1].split('|')
			pos = {} #a dict storing the position of each field we search for, or False if the field is not available
			#find out where certain fields are or if they are not available
			#the fields I'm looking for:
			target_fields = ['Gene', 'SYMBOL', 'BIOTYPE', 'CANONICAL', 'Consequence', 'HGVSp', 'Codons', 'SWISSPROT', 'DOMAINS', 'Feature', 'Protein_position', 'Amino_acids']
			#the information I want printed per variant. F.x. if a transcript is canon does not change and neither does the uniprot ID so no reason to print them for every var 
			print_fields = ['Consequence', 'HGVSp', 'Codons', 'SWISSPROT', 'DOMAINS', 'Protein_position', 'Amino_acids']
			for target_field in target_fields:
				pos[target_field] = False
				for (position, field) in enumerate(fields):
					if field == target_field: 
						pos[target_field] = position
						break
			print('Found vep field information')
			#debug
			if args.v:
				print('Field positions:')
				for target_field in pos:
					print(target_field, pos[target_field]) 
			
			#debug
			break

		else:
			continue
	else:
		#if we didn't break we didn't find the vep line
		print('Did not find annotation metadata line ##INFO=<ID=CSQ, please check vcf file')
		sys.exit()

	#3.3. parse var lines
	for record in vcf_reader:

		#first of, skip SNPs that did not PASS quality control. In vcf reader, PASS evaluates to nothing, so if there is something in the filter you want to skip the entry
		if record.FILTER:
			continue

		#in this vcf readers, the entire vep entry is one string so split it on ','
		#for re-annotation w online vep, change the tag to 'CSQ'. For orig vep annotation in gnomAD files use 'vep' 
		for transcript in record.INFO.get('CSQ').split(','): 
			fields = transcript.split('|')
			
			#we are only interested in coding transcripts and the variant types requested by argv (filter set is adjusted accordingly in step 1)
			if fields[pos['BIOTYPE']] == 'protein_coding' and set(fields[pos['Consequence']].split('&')).intersection(filter_type):

				tID = fields[pos['Feature']]
				gID = fields[pos['Gene']]
				uniprot = fields[pos['SWISSPROT']] if (pos['SWISSPROT'] and fields[pos['SWISSPROT']] ) else 'NA'
				#define the type of variant as missense, syn or complex
				if 'missense_variant' in fields[pos['Consequence']].split('&'):
					vartype = 'mis'
				elif 'synonymous_variant' in fields[pos['Consequence']].split('&'):
					vartype = 'syn'
				else:
					vartype = 'complex'
					
				pref = tID+'_'+gID+'_'+uniprot+'.'+vartype
				#the first time we see this combination of transcript and vartype, write a header. Open file in 'w' mode instead of 'a', this is will also overwrite any existing files.
				if not pref in already_seen_g:
					with open(tmp_folder+'/'+pref +'.genomes.tmp', 'w') as OUT:
						#meta data for this transcript:
						print('#gnomad_version:', gnomad_version, file = OUT)
						print('#transcript_ID:', tID, file = OUT)
						print('#gene_ID:', gID, file = OUT)
						print('#symbol:', fields[pos['SYMBOL']], file = OUT)
						print('#uniprot_ID:', uniprot, file = OUT)
						if fields[pos['CANONICAL']] == '1':
							print('#canonical: 1', file = OUT)
						else:
							print('#canonical: 0', file = OUT)
						print('#extracted_variants:', ','.join(filter_type), file = OUT)	
						
						print('#prism_var AF AC AN Homozygotes basepos', ' '.join(print_fields), file = OUT)
						already_seen_g.add(pref)
					
				#write a line for the variant, no matter if we have written a header or not and this is always in append mod
				with open(tmp_folder+'/'+pref +'.genomes.tmp', 'a') as OUT:
					var_name = parse_HGVSp(transcript, HGSV_field_pos=pos['HGVSp'], prot_pos_field_pos=pos['Protein_position'], aa_field_pos=pos['Amino_acids'], v=args.v)
					print(var_name, record.INFO.get('AF'), record.INFO.get('AC'), record.INFO.get('AN'), record.INFO.get('nhomalt'), record.POS, sep = ' ', end = '', file = OUT)				
					#v04: print an alert if the allel freq is higher than 0.5
					if float(record.INFO.get('AF')) > 0.5:
						print(pref, record.POS, var_name ,record.INFO.get('AF'), record.INFO.get('AN'), file = alertOUT)
					
					for field in print_fields:
						if (pos[field] and fields[pos[field]]):
							print(' ', fields[pos[field]], sep = '', end = '', file = OUT)
						else:
							print(' NA', end = '', file = OUT)
					#lastly print a newline
					print('\n', end = '', file = OUT)



#2.4. print an index using the already_seen list
if args.index:
	print('Writing prefix index to ', tmp_folder, 'idx.txt', sep = '')
	pref_set = already_seen.union(already_seen_g)
	with open(tmp_folder + 'idx.txt','w') as idx_out:
		for item in pref_set:
			print(item, file = idx_out)

#close remaining file handles
alertOUT.close()
