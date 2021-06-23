
#v09: fix again clinvar lookup to repeat if there is a empty return dict
#v08: return not only the significance but also clinvar IDs 
#v07: added splitting the protein change during clinvar lookup because some entries have several consequences (I think it names the consequences for all transcripts).
#v07: try to deal with server errors by returning an error code and writing a log instead of exiting
#v07: changed prism var notation to use '.' instead of '~' consistently
#v07: enable get_uniprot to get all ids if there are several matching to a transcript
#revep: looking for my vep re-annotation in the CSQ field instead of the original gnomAD one in the vep field
#v06: updating the clinvar lookup to catch HTTPError: 429 Client Error: Too Many Requests for url . This happens when you exceed the the Transactions per second (TPS) limit and should result in a sleep + retry later. Added my NCBI API key so I can have more requests per second


#v04 to v05 switch to extracting aa changes to the transcript requested instead of gene ID+canonical field
#v03 to v04 switch to using the ensembl canonical transcript ID to find out if a transcript annotated in vep is the canonical, since I found a case where they disagree and I'm more inclined to believe ensembl 

#import cyvcf2
import re
import sys
import requests
from time import sleep
from requests import HTTPError

def get_uniprot(t_ID, no_isoform = False):
	with open('/storage1/hezscha/data/uniprot_idmapping/HUMAN_9606_idmapping.dat', 'r') as ID_IN:
		ids = []
		for line in ID_IN:
			fields = line.rstrip().split('\t')
			if fields[1] == 'Ensembl_TRS' and fields[2] == t_ID:
				if no_isoform: #split of the isoform number and dash: Q07157-1 to Q07157
					ids.append(fields[0].split('-')[0])					
				else:
					ids.append(fields[0])
		
		if ids:
			uniprot_ID = '|'.join(ids)			
			return(uniprot_ID)
		else:	
			print('Did not find uniprot ID to ', t_ID)
			return("NA") #return the string 'None' instead of actual None since we can still proceed with this

def get_uniprot_seq(uID):
	r = requests.get('https://www.uniprot.org/uniprot/'+uID+'.fasta', headers={ "Content-Type" : "text/plain"})
	if not r.ok:
		r.raise_for_status()
		sys.exit()
	else:
		seq = ''.join(r.text.split('\n')[1:])
		return seq
		#head = r.text.split('\n')[0]
		#return(head,seq)

#v07 deal with server errors (added filehandle and file prefix to complain about)
# ~ def get_prot_seq(transcript, server, fh, prefix):
	# ~ ext = "/sequence/id/"+transcript+"?type=protein"
	# ~ r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})   

	# ~ #v07 deal w server errors
	# ~ if not r.ok:
	  # ~ print(prefix, ':', transcript, r.raise_for_status(), file = fh)
	  # ~ return('NA due to server error')
	  # ~ #sys.exit()
	# ~ else:
		# ~ prot_seq = r.text
		# ~ sleep(0.5)
		# ~ return(prot_seq)

def get_prot_seq(transcript, server, prefix):
	ext = "/sequence/id/"+transcript+"?type=protein"
	
	wait = 10 #the initial wait is 10 seconds
	n_iter = 1
	while True:
		#give up after 3 tries otherwise we'll wait forever
		if n_iter > 3:
			print(prefix, ':', var_name, 'lookup failed at esearch', file = fh)
			return('lookup failed at esearch')
			break
		try:
			r = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
			#r.raise_for_status()
			#v07 deal w server errors
			if not r.ok:
			  print(prefix, ':', transcript, r.raise_for_status(), file = sys.stderr)
			  return('NA due to server error')
			  #sys.exit()
			else:
				prot_seq = r.text
				sleep(0.5)
				return(prot_seq)
		#if something went wrong, try again up to 3 times
		except:
			sleep(wait)
			n_iter += 1
			continue #go back to the start of the while loop and try again		

#this function is no longer used or maintained since we now make prism_clinvar files from a downloaded flat file.
#v08 change the return to a dict?
#v07 deal with server errors (added filehandle and file prefix to complain about)
def clinvar_lookup_dict(assembly, chromosome, basepos, var_name, fh, prefix, server = "https://eutils.ncbi.nlm.nih.gov", give_prot_changes=False):
	
	#step 1: esearch. get the relevant ids
	ext = "/entrez/eutils/esearch.fcgi?db=clinvar&term="+chromosome+"[Chromosome]+"+str(basepos)+"[Base Position for Assembly "+assembly+"]+AND+rev_at_least_one_star[Filter]&retmode=json&retmax=100&api_key=fe42ab6b5e8d76a7de4d795a849a30c54f08"
		
	#v07 deal with server errors
	#request server with exponential backoff if rate is exceeded (code 429)
	#actually I'll just request every 10 seconds for 3 times
	wait = 10 #the initial wait is 10 seconds
	n_iter = 1
	while True:
		#give up after 3 tries otherwise we'll wait forever
		if n_iter > 3:
			print(prefix, ':', var_name, 'lookup failed at esearch', file = fh)
			return('lookup failed at esearch')
			break
		try:
			r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
			r.raise_for_status()
			break
		except HTTPError as e:
			#if it's 429 sleep and try again
			status_code = e.response.status_code
			if status_code == '429':
				sleep(wait)
				n_iter += 1
				continue #go back to the start of the while loop and try again

	decoded = r.json() 
	sleep(1)
		
	#check how many results there were
	#if there were 0 results there is no clinvar entry with at least 1 star at this base position
	if decoded['esearchresult']['count'] == '0':
		#debug
		#print('No clinvar entries for', var_name, 'at', basepos, 'in', chromosome, file = fh)
		#debug
		#changed: return None
		return(None)
		# ~ if give_prot_changes:
			# ~ return('NA','NA','NA','NA')
		# ~ else:	
			# ~ return('NA','NA','NA')
	
	#step 2: get esummaries
	#clinsig = 'NA'
	
	#there can be several entries, i.e. at pos 43060472 in chrom 21, there are G>A (ID:643284) and G>T (ID:212864). So go through them and find the correct one
	for ID in decoded['esearchresult']['idlist']:
		ext = "/entrez/eutils/esummary.fcgi?db=clinvar&id="+ID+"&retmode=json&api_key=fe42ab6b5e8d76a7de4d795a849a30c54f08"
		
		#debug
		print('For', var_name, 'at', basepos, 'looking up', ID, file = fh)
		#debug
		
		#try 3 times
		wait = 10 #the initial wait is 10 seconds
		n_iter = 1
		while True:
			#give up after 3 tries otherwise we'll wait forever
			if n_iter > 3:
				return('lookup failed at esummary to '+ID)
				break
			
			else:
				try:
					r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
					report_error = r.raise_for_status()
				
					if r.ok:
						decoded_esum = r.json()
						if len(decoded_esum['result']['uids']) > 0:
							prot_change = decoded_esum['result'][ID]['protein_change'].split(', ')
							break
					else:
						print('Server error at esummary lookup:', report_error, file = fh)		
				
				except HTTPError as e:
					#if it's 429 sleep and try again
					status_code = e.response.status_code
					if status_code == '429':
						sleep(wait)
						n_iter += 1
						continue #go back to the start of the while loop and try again
					#different HTTP errors:
					else:
						print('Server error at esummary lookup:', e.response, file = fh)
						n_iter += 1
						pass
				
				#other errors, tell me what they are
				except:
					print('Server error at esummary lookup:', e.response, file = fh)
					n_iter += 1
					pass

		
		
		#debug: I got a strange error when looking up entry 421617 where it complained that no entry with key 421617 exists, but when I re-did it manually it was fine. The one thing I can think of are clinvar hickups where the server sometimes doesn't respond. I thought I catch these but seems not. 
		#update: see doc/clinvar 'Clinvar and esearch/esummay/entrez' . There is a kind of random fail where a dict is returned but the key to the accession/ID doesn't exist and decoded_esum['result']['uids'] is an empty list
		#if len(decoded_esum['result']['uids']) == 0:
			
		
		#is this the protein change mentioned in gnomad? Otherwise, look the the next entry
		#we have discovered that 'protein_change' can have several entries like this for K367N in chromosome='21', basepos='45182607': 
		#print('protein change: *', decoded_esum['result'][ID]['protein_change'], '*', sep = '')
		#protein change: *K367N, K416N*
		#therefore, we need to split on comma and see if any of the elements is the change we're looking at
		
		# ~ if ID in decoded_esum['result']:
			# ~ prot_change = decoded_esum['result'][ID]['protein_change'].split(', ')
		# ~ else:
			# ~ print('ID', ID, 'failed', file = fh)
			# ~ print('clinvar ID list', decoded['esearchresult']['idlist'], file = fh)
			# ~ print('decode esummary entry', decoded_esum, file = fh)
			# ~ print(decoded_esum['result'], file = fh)
			# ~ print(decoded_esum['result'].keys(), file = fh)
			# ~ for key in decoded_esum['result']:
				# ~ print(key, decoded_esum['result'][key], file = fh)
			# ~ print(decoded_esum['result'][ID].keys(), file = fh)
			# ~ sys.exit()
		#debug	
		
		if var_name in prot_change:
			rev_stat = re.sub(' ', '_', decoded_esum['result'][ID]['clinical_significance']['review_status'])
			clinsig = re.sub(' ', '_', decoded_esum['result'][ID]['clinical_significance']['description'])
			#break and return this info (and the current ID)
			break
		#bascially redundant since we are at the end of the loop anyway
		else:
			continue	
	#if we didn't break we didn't find an entry that is about the gnomad protein change
	else:
		# ~ ID = 'NA'
		# ~ clinsig = 'NA'
		# ~ rev_stat = 'NA'
		#changed to return 'Entry for different var'. The reason is that if the CV entry is about a diff protein change it is not the right entry for the currently queried var and nothing should be stored so that another var at same pos but w different protein change can be queried again and not wrongly get the info that there is no entry or that it is NA.
		#this is different from above where we know that there are no entries for pos X. Here we know there is an entry for pos X, but it does not concern this variant but another one
		return('Entry for different var') 
	
	if give_prot_changes:
		return(ID,clinsig,rev_stat,prot_change)
	else:	
		return(ID,clinsig,rev_stat)

#the simple version for when I'm not using a dict to store the lookups
def clinvar_lookup(assembly, chromosome, basepos, var_name, fh, prefix, server = "https://eutils.ncbi.nlm.nih.gov", give_prot_changes=False):
	
	#step 1: esearch. get the relevant ids
	ext = "/entrez/eutils/esearch.fcgi?db=clinvar&term="+chromosome+"[Chromosome]+"+str(basepos)+"[Base Position for Assembly "+assembly+"]+AND+rev_at_least_one_star[Filter]&retmode=json&retmax=100&api_key=fe42ab6b5e8d76a7de4d795a849a30c54f08"
		
	#v07 deal with server errors
	#request server with exponential backoff if rate is exceeded (code 429)
	#actually I'll just request every 10 seconds for 3 times
	wait = 10 #the initial wait is 10 seconds
	n_iter = 1
	while True:
		#give up after 3 tries otherwise we'll wait forever
		if n_iter > 3:
			print(prefix, ':', var_name, 'lookup failed at esearch', file = fh)
			return('lookup failed at esearch')
			break
		try:
			r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
			r.raise_for_status()
			break
		except HTTPError as e:
			#if it's 429 sleep and try again
			status_code = e.response.status_code
			if status_code == '429':
				sleep(wait)
				n_iter += 1
				continue #go back to the start of the while loop and try again

	decoded = r.json() 
	sleep(1)
		
	#check how many results there were
	#if there were 0 results there is no clinvar entry with at least 1 star at this base position
	if decoded['esearchresult']['count'] == '0':
		if give_prot_changes:
			return('NA','NA','NA','NA')
		else:	
			return('NA','NA','NA')
	
	#step 2: get esummaries
	#clinsig = 'NA'
	
	#there can be several entries, i.e. at pos 43060472 in chrom 21, there are G>A (ID:643284) and G>T (ID:212864). So go through them and find the correct one
	for ID in decoded['esearchresult']['idlist']:
		ext = "/entrez/eutils/esummary.fcgi?db=clinvar&id="+ID+"&retmode=json&api_key=fe42ab6b5e8d76a7de4d795a849a30c54f08"
		
		#debug
		print('For', var_name, 'at', basepos, 'looking up', ID, file = fh)
		#debug
		
		#try 3 times
		wait = 10 #the initial wait is 10 seconds
		n_iter = 1
		while True:
			#give up after 3 tries otherwise we'll wait forever
			if n_iter > 3:
				return('lookup failed at esummary to '+ID)
				break
			
			else:
				try:
					r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
					report_error = r.raise_for_status()
				
					if r.ok:
						decoded_esum = r.json()
						if len(decoded_esum['result']['uids']) > 0:
							prot_change = decoded_esum['result'][ID]['protein_change'].split(', ')
							break
					else:
						print('Server error at esummary lookup:', report_error, file = fh)		
				
				except HTTPError as e:
					#if it's 429 sleep and try again
					status_code = e.r.status_code
					if status_code == '429':
						sleep(wait)
						n_iter += 1
						continue #go back to the start of the while loop and try again
					#different HTTP errors:
					else:
						print('Server error at esummary lookup:', e.r.status_code, file = fh)
						n_iter += 1
						pass
				
				#other errors, tell me what they are
				except:
					print('Server error at esummary lookup:', e.r.status_code, file = fh)
					n_iter += 1
					pass

		
		
		#debug: I got a strange error when looking up entry 421617 where it complained that no entry with key 421617 exists, but when I re-did it manually it was fine. The one thing I can think of are clinvar hickups where the server sometimes doesn't respond. I thought I catch these but seems not. 
		#update: see doc/clinvar 'Clinvar and esearch/esummay/entrez' . There is a kind of random fail where a dict is returned but the key to the accession/ID doesn't exist and decoded_esum['result']['uids'] is an empty list
		#if len(decoded_esum['result']['uids']) == 0:
			
		
		#is this the protein change mentioned in gnomad? Otherwise, look the the next entry
		#we have discovered that 'protein_change' can have several entries like this for K367N in chromosome='21', basepos='45182607': 
		#print('protein change: *', decoded_esum['result'][ID]['protein_change'], '*', sep = '')
		#protein change: *K367N, K416N*
		#therefore, we need to split on comma and see if any of the elements is the change we're looking at
		
		# ~ if ID in decoded_esum['result']:
			# ~ prot_change = decoded_esum['result'][ID]['protein_change'].split(', ')
		# ~ else:
			# ~ print('ID', ID, 'failed', file = fh)
			# ~ print('clinvar ID list', decoded['esearchresult']['idlist'], file = fh)
			# ~ print('decode esummary entry', decoded_esum, file = fh)
			# ~ print(decoded_esum['result'], file = fh)
			# ~ print(decoded_esum['result'].keys(), file = fh)
			# ~ for key in decoded_esum['result']:
				# ~ print(key, decoded_esum['result'][key], file = fh)
			# ~ print(decoded_esum['result'][ID].keys(), file = fh)
			# ~ sys.exit()
		#debug	
		
		if var_name in prot_change:
			rev_stat = re.sub(' ', '_', decoded_esum['result'][ID]['clinical_significance']['review_status'])
			clinsig = re.sub(' ', '_', decoded_esum['result'][ID]['clinical_significance']['description'])
			#break and return this info (and the current ID)
			break
		#bascially redundant since we are at the end of the loop anyway
		else:
			continue	
	#if we didn't break we didn't find an entry that is about the gnomad protein change
	else:
		ID = 'NA'
		clinsig = 'NA'
		rev_stat = 'NA'
		prot_change = 'NA'
		
	if give_prot_changes:
		return(ID,clinsig,rev_stat,prot_change)
	else:	
		return(ID,clinsig,rev_stat)
	

#clean parse function
#changed to looking for a set of transcript IDs instead of a certain gene ID
def process_exome_vcf(transcript_IDs, chromosome, assembly, vcf_file_ls = [], 
				filter_type=set(['synonymous_variant', 'missense_variant', 'start_lost', 'stop_gained', 'inframe_deletion']), v = False):
	#needs imports of cyvcf2, re
	#attempts to open the vcf file at path 'vcf_file' and extract 
	
	#trying a dict instead of the list to quickly access variants by their variant name f.x. F801E
	list_of_records = {}
	#to get the number of mutated positions. != the number of variants since several variants can be in the same position. We also want to know how many times a mutated position was hit for the depth calc
	mut_pos_d = {}
	
	for tID in transcript_IDs:
		list_of_records[tID] = {}
		mut_pos_d[tID] = {}
	
	#since Kristoffer told me that the ordering doesn't matter there is no need to impose one or remember the order of variants, which will also become messy since we potentially have vars for more than one transcript
	#record_index = [] #order of the variants for printing later
	#find out where certain fields are or if they are not available
	#the fields I'm looking for: 
	#The 'Feature' field has the transcript ID. I need that for the protein sequence and Uniprot ID.
	#v04 EDIT: removed canonical field
	#v05 EDIT: removed Gene field, removed CLIN_SIG field
	target_fields = ['BIOTYPE', 'Consequence', 'HGVSp', 'Codons', 'SWISSPROT', 'DOMAINS', 'Condel', 'CADD_PHRED', 'Feature', 'Protein_position', 'Amino_acids']
	
	#go through the list of vcf files passed and try to open and parse each of them
	for vcf_file in vcf_file_ls:
		try:
			vcf_reader = cyvcf2.VCF(vcf_file)
		except:
			raise	
			#print("Can't open file: ", vcf_file)
			sys.exit()
		
		#try to get info out of the header: I want to order and position of fields in the vep annotation
		#https://github.com/brentp/cyvcf2/issues/94
		#certain header types exist like fileformat, but for example FILTER and INFO don't so it might be necessary to just parse the whole header string
		for line in vcf_reader.raw_header.split('\n'):
			if line.startswith('##INFO=<ID=CSQ'):
			#if line.startswith('##INFO=<ID=vep'):
				#get this part:"... Format: Allele|..."
				fields = line.split('Format: ')[1].split('|')
				pos = {} #a dict storing the position of each field we search for, or False if the field is not available
				for target_field in target_fields:
					pos[target_field] = False
					for (position, field) in enumerate(fields):
						if field == target_field: pos[target_field] = position
				
				print('Found vep field information')
				#debug
				if v:
					print('Field positions:')
					for target_field in pos:
						print(target_field, pos[target_field]) 
				
				#debug
				break

			else:
				continue
		else:
			#if we didn't break we didn't find the vep line
			print('Did not find annotation metadata line ##INFO=<ID=CSQ, please check vcf file', file = sys.stderr)
			sys.exit()
		
		for record in vcf_reader:

			#first of, skip SNPs that did not PASS quality control. In vcf reader, PASS evaluates to nothing, so if there is something in the filter you want to skip the entry
			if record.FILTER:
				#print(record.FILTER)
				continue
			
			#in this vcf readers, the entire vep entry is one string so split it on ','
			#for re-annotation w online vep, change the tag to 'CSQ'. For orig vep annotation in gnomAD files use 'vep' 
			for transcript in record.INFO.get('CSQ').split(','): 
			#for transcript in record.INFO.get('vep').split(','): 
				
				fields = transcript.split('|')
				#v05 get only changes to the transcripts requested (don't need to check gene ID
				if fields[pos['BIOTYPE']] == 'protein_coding' and fields[pos['Feature']] in transcript_IDs and set(fields[pos['Consequence']].split('&')).intersection(filter_type):
					
					#Update: gnomAD stores only one ALT allel per line, probably to be clear about the allel counts and frequencies.
					var_name = parse_HGVSp(transcript, HGSV_field_pos=pos['HGVSp'], prot_pos_field_pos=pos['Protein_position'], aa_field_pos=pos['Amino_acids'], v=v)
					#I have seen cases of empty HGVSp fields (thought the nucl level field is not empty, see example_VEP_annotation). In those case skip to the next transcript. I won't try to guess the protein consequence if VEP didn't find one
					if not var_name:
						continue
					tID = fields[pos['Feature']] 

					#if this var is already in the list, we need to aggregate freqs
					if var_name in list_of_records[tID]:
						#we're not really super using the entire transcript, so just leave that alone for now
						#add the allel count
						list_of_records[tID][var_name]['AC'] += record.INFO.get('AC')
						#find the max of the allel numbers (either one already saved or the one in the current line
						list_of_records[tID][var_name]['AN'] = max(int(list_of_records[tID][var_name]['AN']), int(record.INFO.get('AN')))
						#calc allel freq
						list_of_records[tID][var_name]['AF'] = list_of_records[tID][var_name]['AC'] / list_of_records[tID][var_name]['AN']
						#add the codon
						if pos['Codons']:
							#add the codon if not the same. Since I load both exome and genome files a variant can be mentioned in both with the same codon and then I only want that codon once
							known_codons = list_of_records[tID][var_name]['Codons'].split(',')
							#if current codon does not match the already saved codon, add current
							if not fields[pos['Codons']] in known_codons:
								list_of_records[tID][var_name]['Codons'] += ',' + fields[pos['Codons']]		
						else:
							list_of_records[tID][var_name]['Codons'] = 'NA'
						
					else:
						list_of_records[tID][var_name] = {}
						
						list_of_records[tID][var_name]['vep'] = transcript
						list_of_records[tID][var_name]['AC'] = record.INFO.get('AC')
						list_of_records[tID][var_name]['AN'] = record.INFO.get('AN')
						list_of_records[tID][var_name]['AF'] = record.INFO.get('AF')
						list_of_records[tID][var_name]['ID'] = record.ID
						
						#the same variant, i.e. protein change should have the same clinvar entry since I don't think it's based on the nucl change? I'm trying to say we only to look this up once, even if subsequent nucl changes lead to the same variant
						list_of_records[tID][var_name]['clinsig'] = clinvar_lookup(assembly=assembly, chromosome=chromosome, basepos=record.POS, var_name=var_name)

						for target_field in target_fields:
							list_of_records[tID][var_name][target_field] = fields[pos[target_field]] if (pos[target_field] and fields[pos[target_field]]) else 'NA'

						prot_pos = fields[pos['Protein_position']]
						if not prot_pos in mut_pos_d[tID]:
							mut_pos_d[tID][prot_pos] = 1
						else:
							mut_pos_d[tID][prot_pos] += 1  

	return(list_of_records, filter_type, mut_pos_d)

def parse_HGVSp(transcript, HGSV_field_pos, prot_pos_field_pos, aa_field_pos, v = False):	
	text = transcript.split('|')[HGSV_field_pos]
	
	#include end of string '$' in the match to distinguish single muts from longer insertions
	#the ? is to match M1?, start lost. The = was for silent mutations but I think those are coded %3D in gnomad VEP annotation
	#change the pattern to only match Upper_case{1}lower_case{2} so things like 'dup' don't get matched and only amino acid 3-letter names are caught
	#single_sub_pattern = r'([A-Za-z]{3})(\d{1,4})([A-Za-z]{3}|=|\?)$'
	single_sub_pattern = r'([A-Z]{1}[a-z]{2})(\d+)([A-Z]{1}[a-z]{2}|=|\?)$'

	#new case: when parsing spliceAI data sometimes the HGVSp field is just empty
	if not text:
		if v:
			print('No protein change is recorded for this DNA level variant.')
		return(None)

	#substitution case:
	elif re.search(single_sub_pattern, text):
		single_mutation_re = re.search(single_sub_pattern, text)
		wt_re, pos_re, mut_re = single_mutation_re.groups()
		if v:
			print(single_mutation_re)
			print(wt_re)
			print(pos_re)
			print(mut_re)

			#debug
			try:
				 wt = aa_threeletter_to_letter[wt_re]
			except:
				print(wt_re)
				
			try:
				mut = aa_threeletter_to_letter[mut_re]
			except:
				print(mut_re)
				
		else:
			wt = aa_threeletter_to_letter[wt_re]
			mut = aa_threeletter_to_letter[mut_re]			 

		return(wt+pos_re+mut)

	#if the expression did not match, try to match a silent mutation
	#changed regex from '\(p.%3D\)' to just matching %3D since the notation changed
	elif re.search('%3D', text):
		#now we need to get the position wrt the amino acid and the (unchanged) amino acid. For all other consequence this info is encoded in the HGVSp but not here
		aa = transcript.split('|')[aa_field_pos]
		pos = transcript.split('|')[prot_pos_field_pos]
		return(aa+pos+'=')
	
	
	#delins case first so I can differentiate it from the only del and only ins cases
	#examples:
	#Cys28delinsTrpVal			- C28.:.28T:29V
	#Asp224_Gln225delinsGlu		- D224.:Q225.:.224E
	#Pro578_Lys579delinsLeuTer	- P578.:K579.:.578L:.579*
	elif re.search('(\w{3})(\d+)(?:_(\w{3})(\d+))*delins(.*)', text):
		bla = re.search('(\w{3})(\d+)(?:_(\w{3})(\d+))*delins(.*)', text)
		#everything before the delins (and everything between it) is deleted, everything after it is added
		#deletion
		#if there is a second aa before the delins we're deleting a range
		if bla.groups()[2]:
			start_pos = int(bla.groups()[1])
			end_pos = int(bla.groups()[3])
			prism_expr = aa_threeletter_to_letter[bla.groups()[0]] + str(start_pos) + '.'
			for i in range(1,end_pos - start_pos):
				#I guess I should try to get the actual aa's from the protein sequence of the transcript (in the info_d) but let's see if we actually need that
				prism_expr += ':X' + str(start_pos+i) + '.'
			prism_expr += ':' + aa_threeletter_to_letter[bla.groups()[2]] + str(end_pos) + '.'
			
		
		else:
			prism_expr = aa_threeletter_to_letter[bla.groups()[0]] + bla.groups()[1] + '.'
	
		#now insertion. Just create one position for every three letters after the delins, i.e. group 4
		#there is always something being deleted before since it is a delins, so start with the delimiter ':'
		prism_expr += ':'
		pos = int(bla.groups()[1])
		insert_aas = []
		for i in range(0,len((bla.groups()[4])),3):
			aa=bla.groups()[4][i:i+3]
			#insert_aas.append('.'+str(pos)+aa_threeletter_to_letter[aa])
			insert_aas.append('.'+str(pos)+aa_threeletter_to_letter[aa])
			#the first inserted aa is in the same place as the first deleted one, so we count up after having added to the prism expression
			pos += 1
		
		prism_expr += ':'.join(insert_aas)	
		return(prism_expr)
	
	#simple deletion case: I noticed I didn't catch this one: Glu28del which should be E28.
	elif re.search('(\w{3})(\d+)del', text):
		bla = re.search('(\w{3})(\d+)del', text)
		prism_expr = aa_threeletter_to_letter[bla.groups()[0]] + bla.groups()[1] + '.'
		return(prism_expr)
	
	#make a match case for long insertions and deletions.
	#ins: [proteinID]:p.Gln597_Gly598insPheCysTerSerTyrSerLys
	#ins: ENSP00000264554.4:p.Gly129_Ser130insIleArgAsp to .130I:.131R:.132D
	#del: NP_003997.1:p.Lys23_Val25del deletes all positions from Lys23 to Val25, so in total 3 amino acid, translate to prism expression: L23.:X24.:V25.
	#we're dealing with delins further up so no case here
	elif re.search('(\w{3})(\d+)_(\w{3})(\d+)(ins|del)(.*)', text):
		bla = re.search('(\w{3})(\d+)_(\w{3})(\d+)(ins|del)(.*)', text)
		prism_expr = ''
		
		if bla.groups()[4] == 'ins':
			#first, split up the last group by 3 to get all the inserted aa
			pos = int(bla.groups()[1])
			insert_aas = []
			for i in range(0,len((bla.groups()[5])),3):
				pos += 1 #we start one after the position of the first named aa (see example). Then each subsequently inserted aa needs to be one position later
				aa=bla.groups()[5][i:i+3]
				#insert_aas.append('.'+str(pos)+aa_threeletter_to_letter[aa])
				insert_aas.append('.'+str(pos)+aa_threeletter_to_letter[aa])
			
			prism_expr = ':'.join(insert_aas)

		#else it's a long deletion. The aa's between the two flanking ones (that are also deleted) are not named in the HGVSp, but never fear, I can get them from the amino acids field. F.x. for E26.:X27.:E28. which is ENSP00000482968.1:p.Glu26_Glu28del, the amino acids field tells us that we lost 3 E's: LEEE/L	
		else:
			#upon reading the reference for protein deletion again I come to the conclusion that the deleted AA inbetween the flanking are never named, so we don't need to read groups()[5] 
			start_pos = int(bla.groups()[1])
			del_aas = []
			#add the amino acids named as deleted in the amino acids field. skip the first one, that one stays.
			for ppos,aa in enumerate(transcript.split('|')[aa_field_pos].split('/')[0][1:]):
				del_aas.append(aa + str(start_pos + ppos) + '.')
			
			prism_expr = ":".join(del_aas)
			
		return(prism_expr)		
		
	#and a dup case
	#ENSP00000269812.1:p.Arg161_Leu162dup
	#ENSP00000482968.1:p.Glu28dup		1 extra E	L/LE	.29E
	#ENSP00000482968.1:p.Glu27_Glu28dup	2 extra Es	L/LEE	.29E:.30E
	#ENSP00000482968.1:p.Glu26_Glu28dup	3 extra Es	-/-EEE	.29E:.30E:.31E
	elif re.search('(?:(\w{3})(\d+)_)*(\w{3})(\d+)dup$', text):		
		bla = re.search('(?:(\w{3})(\d+)_)*(\w{3})(\d+)dup$', text)
		#if the first two capturing groups are empty it is a simple duplication
		if bla.groups()[0] is None and bla.groups()[1] is None:
				prism_expr = '.' + aa_threeletter_to_letter[bla.groups()[2]] + bla.groups()[3]
		
		#else the range between the two flanking aa including them has been duplicated.
		#use the same strategy as for long deletions
		else:
			start_pos = int(bla.groups()[1])
			insert_aas = []
			#add the amino acids named as inserted in the amino acids field. skip the first one, that one was there already and is sort of the anchor, BUT doesn't need to relate to where the insertion is named to be by the HGVSp (read debug/ALM1_ENTS... file to understand)
			for ppos,aa in enumerate(transcript.split('|')[aa_field_pos].split('/')[1][1:]):
				insert_aas.append('.' + str(start_pos + ppos) + aa)

			prism_expr = ":".join(insert_aas)

		return(prism_expr)

############
#???pending stop loss case
#	elif 'extTer' in text:
		
##############		

	#TODO
	#in the final version the else case should be to complain. Until I find out how I want to treat long insertions and deletions, the else case will be to just return the HGVSp field
	#~ else:
		#~ print('Could not parse variant', text)
		#~ return(None)

	#tmp else case!
	else:
		print('HGVSp field:', text, 'could not be parsed and will be skipped', file = sys.stderr)
		return(None)
		#return(text)
		
#~ def parse_HGVSp(transcript, v = False):
	#~ #edit: to take the entire vep field for a transcript, cut out the HGVSp field then
	#~ #parse it into single letter pos single letter variant naming
	#~ text = transcript.split('|')[11]

	#~ #new 
	#~ #include end of string '$' in the match to distinguish single muts from longer insertions
	#~ single_sub_pattern = r'([A-Za-z]{3})(\d{1,4})([A-Za-z]{3}|=|\?)$'

	#~ #me from here:
	#~ if re.search(single_sub_pattern, text):
		#~ single_mutation_re = re.search(single_sub_pattern, text)
		#~ wt_re, pos_re, mut_re = single_mutation_re.groups()
		#~ if v:
			#~ print(single_mutation_re)
			#~ print(wt_re)
			#~ print(pos_re)
			#~ print(mut_re)

		#~ #debug
		#~ try:
			 #~ wt = aa_threeletter_to_letter[wt_re]
		#~ except:
			#~ print(wt_re)
			
		#~ try:
			#~ mut = aa_threeletter_to_letter[mut_re]
		#~ except:
			#~ print(mut_re)		 

		#~ return(wt+pos_re+mut)

	#~ #if the expression did not match, try to match a silent mutation
	#~ #changed regex from '\(p.%3D\)' to just matching %3D since the notation changed
	#~ elif re.search('%3D', text):
		#~ #now we need to get the position wrt the amino acid and the (unchanged) amino acid. For all other consequence this info is encoded in the HGVSp but not here
		#~ aa = transcript.split('|')[15]
		#~ pos = transcript.split('|')[14]
		#~ return(aa+pos+aa)
		
	#~ #else this is not a simple substitution, try to match a long insertion/deletion pattern	
	#~ #elif ????:
				

	#~ #TODO
	#~ #in the final version the else case should be to complain. Until I find out how I want to treat long insertions and deletions, the else case will be to just return the HGVSp field
	#else:
	#	print('Could not parse variant', text)
	#	return(None)

	#~ #tmp else case!
	#~ else:
		#~ return(text)


#from Magnus' mavedb_to_prism.py
# HVGS encoding
# http://www.hgvs.org/mutnomen/references.html
aa_letter_to_threeletter = {
	"A": "Ala",
	"V": "Val",
	"I": "Ile",
	"L": "Leu",
	"M": "Met",
	"F": "Phe",
	"Y": "Tyr",
	"W": "Trp",
	"S": "Ser",
	"T": "Thr",
	"N": "Asn",
	"Q": "Gln",
	"C": "Cys",
	"G": "Gly",
	"P": "Pro",
	"R": "Arg",
	"H": "His",
	"K": "Lys",
	"D": "Asp",
	"E": "Glu",
	"B": "Asx",
	"U": "Sec",
	"X": "Xaa",
	"Z": "Glx",
	"*": "*", #Termination
	"*": "Ter", # Termination
	#".": "del", # Deletion
	"~": "?" #start lost
}

aa_threeletter_to_letter = {v: k for k, v in aa_letter_to_threeletter.items()}


#the func that def works but has a lot of shit in it
def process_exome_vcf_legacy(vcf_file, gene_ID, 
				#changed from list to set so I can intersect with the split list of all consequences 
				filter_type=set(['synonymous_variant', 'missense_variant', 'start_lost', 'stop_gained', 'inframe_deletion']), 
				canonical = True):
	#needs imports of cyvcf2, re
	#attempts to open the vcf file at path 'vcf_file' and extract 
	
	#trying a dict instead of the list to quickly access variants by their variant name f.x. F801E
	list_of_records = {}
	record_index = [] #order of the variants for printing later
	
	#debug
	count = 0
	c_trans = 0
	tmp_OUT = open('tmp.txt', 'w')
	#first = True
	#debug
	
	try:
		vcf_reader = VCF(vcf_file)
	except:	
		print("Can't open file: ", vcf_file)
		#sys.exit()
	
	#try to get info out of the header: I want to order and position of fields in the vep annotation
	#https://github.com/brentp/cyvcf2/issues/94
	#certain header types exist like fileformat, but for example FILTER and INFO don't so it might be necessary to just parse the whole header string
	for line in vcf_reader.raw_header.split('\n'):
		#debug
		print(line, file = tmp_OUT)
		#debug
		if line.startswith('##INFO=<ID=CSQ'):
			#get this part:"... Format: Allele|..."
			fields = line.split('Format: ')[1].split('|')
			pos = {} #a dict storing the position of each field we search for, or False if the field is not available
			#find out where certain fields are or if they are not available
			#the fields I'm looking for: 
			for target_field in ['Gene', 'BIOTYPE', 'CANONICAL', 'Consequence', 'HGVSp', 'Codons', 'SWISSPROT', 'DOMAINS', 'Condel', 'CADD_PHRED', 'CLIN_SIG']:
				pos[target_field] = False
				for (position, field) in enumerate(fields):
					if field == target_field: pos[target_field] = position
			
			print('Found vep field information')
			break

		else:
			continue
	else:
		#if we didn't break we didn't find the vep line
		print('Did not find ##INFO=<ID=CSQ, please check vcf file')
		
	#~ print(vcf_reader.get_header_type('fileformat'))
	#~ print(vcf_reader.get_header_type('fileformat')['fileformat'])
	#~ print(vcf_reader.get_header_type('FILTER'))
	#~ #print(vcf_reader.raw_header, file = tmp_OUT)
	#~ print(vcf_reader.get_header_type('INFO'))
	
	for record in vcf_reader:
		#debug
		#print(record)
		#print('#######################')
		#print(record.INFO.get('vep'))
		#debug
		
		
		#first of, skip SNPs that did not PASS quality control. In vcf reader, PASS evaluates to nothing, so if there is something in the filter you want to skip the entry
		if record.FILTER:
			#print(record.FILTER)
			continue
		
		#in this vcf readers, the entire vep entry is one string so split it on ','
		#for transcript in record.INFO.get('vep').split(','): 
		#for re-annotation w online vep, change the tag to CSQ
		for transcript in record.INFO.get('CSQ').split(','): 
			#debug
			#print(transcript, file = OUT)
			#if first:
			#	print(transcript, file = OUT)
			#	first = False
			#debug
			
			fields = transcript.split('|')
			#the first thing we check: is this line related to the protein/gene that was queried? And is the 
			#currently read transcript the canonical and coding and the protein consequence among those filtered 
			#for in 'filter_type'? canonical is field nr 26 and Ensemble gene ID should be in field 4 (0-indexed). 
			#field 7 Biotype should be 'protein_coding'
			#I don't think it's necessary right now to make a special case for wanting only canonical transcripts and therefore there is no reason to check it every record
			#if canonical:
			
			#the consequence field can have several entries separated by & and we want to know if any of them are in the list of consequence we're looking for (filter_type)
			#adjusted field numbers for online vep annotation
			if fields[4] == gene_ID and fields[7] == 'protein_coding' and fields[23] == 'YES' and set(fields[1].split('&')).intersection(filter_type): 
			#if fields[4] == gene_ID and fields[7] == 'protein_coding' and fields[26] == 'YES' and set(fields[1].split('&')).intersection(filter_type): 
				
				#Update: gnomAD stores only one ALT allel per line, probably to be clear about the allel counts and frequencies.
				c_trans += 1
				#I want to try storing the info in dict that uses the variant notation as keys, so get that first
				var_name = parse_HGVSp(transcript)

				#if this var is already in the list, we need to aggregate freqs
				if var_name in list_of_records:
					#we're not really super using the entire transcript, so just leave that alone for now
					#add the allel count
					list_of_records[var_name]['AC'] += record.INFO.get('AC')
					#find the max of the allel numbers (either one already saved or the one in the current line
					list_of_records[var_name]['AN'] = max(int(list_of_records[var_name]['AN']), int(record.INFO.get('AN')))
					#calc allel freq
					list_of_records[var_name]['AF'] = list_of_records[var_name]['AC'] / list_of_records[var_name]['AN']
					#add the codon
					if pos['Codons']:
						list_of_records[var_name]['Codons'] += ',' + fields[pos['Codons']]		
					else:
						list_of_records[var_name]['Codons'] = 'NA'
					
					#debug
					print(var_name)
					#debug
					
				else:
					list_of_records[var_name] = {}
					
					list_of_records[var_name]['vep'] = transcript
					list_of_records[var_name]['AC'] = record.INFO.get('AC')
					list_of_records[var_name]['AN'] = record.INFO.get('AN')
					list_of_records[var_name]['AF'] = record.INFO.get('AF')
					list_of_records[var_name]['ID'] = record.ID

					for target_field in  ['HGVSp', 'Codons', 'SWISSPROT', 'DOMAINS', 'Condel', 'CADD_PHRED', 'CLIN_SIG']:
						list_of_records[var_name][target_field] = fields[pos[target_field]] if (pos[target_field] and fields[pos[target_field]]) else 'NA'

					record_index.append(var_name)
					count += 1				  

	#debug
	print('Eligble transcripts', c_trans)
	print('protein products', count)
	tmp_OUT.close()
	#debug
	
	return(list_of_records, record_index)

