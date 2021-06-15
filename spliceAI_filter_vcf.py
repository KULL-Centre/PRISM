
#script to go through the spliceAI vcf files and pick out the lines where any delta is >= 0.5 in accordance with their recommended cutoff. 
#can make a second version later where we pick out anything >= 0.2 for recall
#the thing is that most lines in the vcf have no splice altering and we're not interested in these or in annotation them with vep. So we should select the snvs we're interested in first and then run annotation

import argparse
import sys
import cyvcf2
import re
import os

#parse arguments
################################################################################
parser = argparse.ArgumentParser()
parser.add_argument('-chrom', dest="chromosome", type = str, required=True, help="The chromosome for which to parse the gnomad_v2 exomes vcf and gnomad_v3 genomes vcf")
parser.add_argument('-filter', dest="filter", type = float, required=True, help="SpliceAI delta value to filter on. If one of the four deltas is greater equal the filter, the line passed to the filtered file, otherwise not.")
parser.add_argument('-v', dest="v", action='store_true', help="Verbose")
args = parser.parse_args()
################################################################################

data_dir = '/storage0/shared/data/spliceAI/'
# ~ out_dir = '/storage1/hezscha/data/spliceAI_step11_files/'
out_dir = '/storage0/shared/data/spliceAI/'

#vcffile = data_dir+'spliceai_scores.MLH1.vcf'
vcffile = data_dir+'spliceai_scores.masked.snv.hg38.chromosome'+args.chromosome+'.vcf.gz'

#outfile = open(out_dir+'spliceai_scores.MLH1.filter.spliceAI'+str(args.filter)+'.vcf', 'w')
outfile = open(out_dir+'spliceai_scores.masked.snv.hg38.chromosome'+args.chromosome+'.filter.spliceAI'+str(args.filter)+'.vcf', 'w')

#open input vcf
print('Parsing vcf', vcffile, '...')
try:
	vcf_reader = cyvcf2.VCF(vcffile)
except:
	raise
	sys.exit()

#just write the header to the new file:
print(vcf_reader.raw_header, end = '', file=outfile)

#go through the body and pick out the positions that fullfill the filter criteria
for record in vcf_reader:
	#first of, skip SNPs that did not PASS quality control. In vcf reader, PASS evaluates to nothing, so if there is something in the filter you want to skip the entry
	if record.FILTER:
		continue
		
	#we are only interested in DNA level snvs that are splice altering, so only those that have one of the deltas >= 0.5:
	#https://github.com/Illumina/SpliceAI
	sfields = record.INFO.get('SpliceAI').split('|')
	if any(float(x) >= args.filter for x in sfields[2:6]):
		print(record, end = '', file = outfile)







outfile.close()
