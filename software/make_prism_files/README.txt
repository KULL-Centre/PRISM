1. Creating prism_uniprot files
------

Files are created by calling make_prism_uniprot_file.py on either a list of uniprot IDs or a file containing that list.

./make_prism_uniprot_file.py -uniprot P04637 
./make_prism_uniprot_file.py -uniprot 'P04637,P35520' 
./make_prism_uniprot_file.py -fromfile my_uniprot_list.txt

arguments:
  -uniprot UNIPROT      Comma separated list of uniprot IDs (no white space)
  -fromfile FROMFILE    A file from which to read uniprot IDs (one per line)
  -outdir OUTDIR        Output directory for prism_uniprot files. Will be
                        current working directory if not given.
  -version VERSION      Version of the files we're creating. Right now we are
                        on version 002 since the addition of the mobidb
                        disorder priority column (named mobi_db_consensus).                        
  -m {overwrite,leave,check}  
                        What do when the output file already exists. Leave
                        (default), overwrite the file or check the version and 
                        overwrite if the version of the new file would be 
                        higher.
  -v, --verbose         Level of output.
  -fail FAIL            A file listing uniprot IDs for which extraction failed
                        before (i.e. because they have no features of interest) 
                        so we don't waste time querying for them.

At least -uniprot or -fromfile are required.
The default destination folder is /storage1/shared/data/prism_uniprot/[uniprotID[0:2]]/[uniprotID[2:4]]/[uniprotID[4:6]]/. 

For the purpose of merging with other prism files you need to first expand it to all variants with:

./FillVariants.py prism_uniprot_XXX_P00374.txt

The merge wrapper detailed below does that automatically.

File Format
------

Header
------

The header is in machine readable YAML format and contains information on the protein and the available info. It is divided into fields and subfields, noticable by indent.

The fields are
**version** [number] : version of the data set used to check if the data has changed

**Protein** sub fields give information on the target protein

- name [text] : Typically the gene name or a common name of the protein, possible followed by the domain considered. Need not match the file names and is mostly meant for human readability.
- organism [text] : The source of the target protein. Mostly ment for human readability.
- uniprot [text] : The uniprot ID. This will be the uniprot ID requested at file creation.
- sequence [word] : The sequence the information in the data section applies to. All amino acid given in the data section need to match this.

**Column** fields describe the recorded features and measured quantities given in the file. The first column (variant, residue or similar) is not described in the header. Other columns are named and described here and must appear in order of appearence in the data section. A column name must be a single word without white space since the columns are white-spece separated.

**uniprot** sub fields give information on the uniprot entry the data derives from
- reviewed [text] : review status
- uncertain_placement [text] : indicates features where one or both of the positions are uncertain, i.e. DISULFID_265..?274 in A4D0S4. These features are recorded in the data section.
- missing_placement [text] : indicates features who's start or end position is unknown. These features are not recorded in the data section.

Mandatory header fields:
- version (always)
- protein (always) -> sequence
- columns (always) -> all named columns
- uniprot (since this file type is prism_uniprot)


Data Section
------

The first line of the Data Section (no preceeding #) are the column headers. All other lines contain info on the position named in the beginning of the line. 
Note that the recorded features correspond to the wild type as given in the uniprot entry (the sequence in the header). For the purpose of merging with other, variant based prism files we assume them to always pertain to that position, irrespective of whether it is the WT or a variant which is of course an approximiation since some substitutions can destroy certain features like secondary structure elements or active sites. 

Columns are separated by a single space.
There are currently a maximum of 29 features available, however for most uniprot entries not all features apply and prism_uniprot files therefore have varying numbers of columns. All present columns are always named in the header section.

White space in fields has been replaced with '_' since that is the column separator, i.e. Transcription_activation_(acidic). Each field also contains a number separated by '|' to count up the feature instance. For example, the first disulfide bond in a protein is DISULFID|1, the second one is DISULFID|2 and so on.

------
2. Creating prism_gnomad files
------
A python parser to take a protein idenitfier and a gnomAD vcf file and produce a prism file for the protein

2.0. Setup step1 files
This was done already on the binf system, so just info unless things need redoing. Skip to 2.1. Prerequists or 2.2. Prism file generation

2.0.1. downloaded vcf files from gnomad:
https://gnomad.broadinstitute.org/downloads
https://gnomad.broadinstitute.org/faq

The FAQ details that: 
"The gnomAD v2.1 data set contains data from 125,748 exomes and 15,708 whole genomes, all mapped to the GRCh37/hg19 reference sequence. The gnomAD v3.1 data set contains 76,156 whole genomes (and no exomes), all mapped to the GRCh38 reference sequence. Most of the genomes from v2 are included in v3.1." 

I therefore used exome data from gnomad v2 and genome data from gnomad v3 to get the best available coverage in terms of both variants seen and their allel frequencies in the population.

2.0.2. annotation with Ensembl's Variant Effect Predictor

However, gnomad v2 vcfs are only annotated for protein variants for the transcript set of GRCh37. Since we aim to combine the v2 and v3 data, we used the GRCh38 liftover of the v2 exome vcf (also downloaded from gnomad) and then annotated with Variant Effect Predictor (VEP) for GRCh38. At the time the project was started there were also vital fields missing from the VEP annotation of the gnomad v3 genome data, so I re-annotated those vcfs as well for GRCh38, same VEP version and same transcript set as used for the v2 exome liftover.

I managed to install vep on deic using the following procedure:

Prerequisite: conda
#add common conda install to your path. Perhaps check with deic/ the sbinlab wiki if this is still the recommended way to do it
$ echo ". /lustre/hpc/sbinlab/software/anaconda2/etc/profile.d/conda.sh" >> ~/.bashrc
#or install your own mini conda

This instruction installs a conda package which is version 95 of vep. No other method worked when I tried. The current version is 100 and some transcripts have changed between releases 95 and 100, see https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html

#make env:
conda create -n conda_vep
environment location: /groups/sbinlab/hezscha/.conda/envs/conda_vep

#activate env
#source .bashrc #if it can't find the conda command, source the .bashrc
conda activate conda_vep

#install the vep package. Might want to search for a newer one than version 95
conda install -c bioconda/label/cf201901 ensembl-vep

#find where this was actually installed
which vep_install
~/.conda/envs/conda_vep/bin/vep_install
#the executable for running vep is in the same dir: /groups/sbinlab/hezscha/.conda/envs/conda_vep/bin/vep

#test
cd /groups/sbinlab/hezscha/.conda/envs/conda_vep/bin
./vep

#get the latest cache you can find, even if it is newer than the vep version. The later scripts are querying Ensembl database for info with gene and transcript IDs from the annotation and those occasionally change between cache versions. The Ensembl database always be on the latest cache version which is why you should also use the latest cache for the anno, otherwise there will be mismatches. The below example is for cache version 100.

#get the cache manually:
#check documentation here: https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache

mkdir /groups/sbinlab/hezscha/software/vep/cache_GRCh38_100/
cd /groups/sbinlab/hezscha/software/vep/cache_GRCh38_100/
curl -O ftp://ftp.ensembl.org/pub/release-100/variation/indexed_vep_cache/homo_sapiens_vep_100_GRCh38.tar.gz
#start interactive node because the untar may be computionally heavy
srun -p sbinlab_ib --pty bash
#load conda env
conda activate conda_vep
cd /groups/sbinlab/hezscha/software/vep/cache_GRCh38_100
tar xzf homo_sapiens_vep_100_GRCh38.tar.gz

#I checked and there is no fasta for the HGVS notation, so try to get that
vep_install -a f -s homo_sapiens -y GRCh38 --CACHE_VERSION [put version here] -c /groups/sbinlab/hezscha/software/vep/cache_GRCh38_100

#now test the manual cache
/groups/sbinlab/hezscha/.conda/envs/conda_vep/bin/vep -i /groups/sbinlab/hezscha/test/ensembl-vep/examples/homo_sapiens_GRCh38.vcf --cache --dir_cache /groups/sbinlab/hezscha/software/vep/cache_GRCh38_100 --fasta /groups/sbinlab/hezscha/software/vep/cache_GRCh38_100/homo_sapiens/100_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz -o /groups/sbinlab/hezscha/test/vep_test_deic/homo_sapiens_GRCh38.100.vcf --vcf --offline --variant_class --sift b --polyphen b --distance 0 --hgvs --uniprot --domains --canonical --gene_phenotype --force_overwrite --CACHE_VERSION [put version here]

2.0.3. temporary files 

To avoid reading in the huge vcf files every time, I made a directory with sub directories per chromosome that contains temporary files for each transcript. 

Specifically per transcript ID there are:
1 file with missense mutations in exome data: .mis.exomes.tmp
1 file with missense mutations in genome data: .mis.genomes.tmp
1 file with synonymous mutations in exome data: .syn.exomes.tmp
1 file with synonymous mutations in genome data: .syn.genomes.tmp
1 file with complex mutations in exome data: .complex.exomes.tmp
1 file with complex mutations in genome data: .complex.genomes.tmp

Files may not exist if there was no such data. I.e. if missense mutations were found for a certain transcript in the genome data but not in the exome data (i.e. due to the genomic region not being in the exome data), only the .mis.genomes.tmp will exist, not the .mis.exomes.tmp . The complex files are experimental and only exist for chrom Y.

These tmp files are made with: 
python3 parse_whole_chromosome_split_vartypes_v04.py

required arguments:
  -chrom CHROMOSOME     The chromosome for which to parse the gnomad_v2 exomes
                        vcf and gnomad_v3 genomes vcf

optional arguments:
  -filter FILTER        The complex variant types to extract, given as a comma
                        separeted list (no whitespace!). Missense and
                        synonymous variants are always extracted (unless
                        omitted in -e option). The currently available types
                        are: stop gained : sg, start lost : sl, inframe
                        deletion : idel, inframe insertion : iins. Default is
                        "sg,sl,idel,iins"
  -e {mis,syn,complex,all}
                        Type of variants to write files for. Choose one of the
                        list. The current default is 'mis'. You probably only 
                        want missense ('mis') and synonymous ('syn') since the 
                        notation for complex variants is still experimental.
                    
  -tmp_folder TMP_FOLDER
                        Where to store the output files. Default location
                        is /storage1/hezscha/gnomad_to_prism_parser/step1_file
                        s/+args.chromosome/
  -v                    Verbose
  -index                Write an index of prefixes into the tmp dir.
                        Automatically on if -e is 'all'. Default is off (as
                        for all arguments with the action store_true)
  -d {exomes,genomes,both}
                        Choose whether to extract data from 'exomes',
                        'genomes' or 'both'(default).

The files are in:
/storage1/hezscha/gnomad_to_prism_parser/step1_files/
and in a tar archive:
/storage1/hezscha/gnomad_to_prism_parser/gnomad_step1_files.tgz

The tmp_folder option of the gnomad file generating scripts below points there automatically, change it if you untar the tar file somewhere else.

2.1. Prerequisits 

python libraries: 
cyvcf2
argparse
sys
re
os
requests
datetime
numpy
time

2.2. Prism file generation 

To generate prism compatible files with gnomad data, run either make_prism_gnomad_file_seq.py or make_prism_gnomad_file_all.py. 

usage: make_prism_gnomad_file_all.py [-h] [-e {mis,syn,complex,all}]
                                     [-m {overwrite,leave}]
                                     [-tmp_folder TMP_FOLDER]
                                     [-out_folder OUT_FOLDER] [-d] -uniprot
                                     UNIPROT [-v]

Minimally requires a uniprot ID. Will make all possible prism gnomad files to that ID, one per unique protein product.

required arguments:
  -uniprot UNIPROT      Uniprot ID to get transcripts and make prism_gnomad
                        files for

optional arguments:
  -h, --help            show this help message and exit
  -e {mis,syn,complex,all}
                        Type of variants to write files for.
  -m {overwrite,leave}  What do when the output file already exists. Leave
                        (default) or overwrite
  -tmp_folder TMP_FOLDER
                        Where the temporary step 1 files are stored. Default
                        location is /storage1/hezscha/gnomad_to_prism_parser/s
                        tep1_files/+args.chromosome/
  -out_folder OUT_FOLDER
                        Where the output files should be written. Default
                        location is /storage1/shared/data/prism_gnomad/uniprot
                        [0:2]/unipro[2:4]/uniprot[4:6]/
  -d                    Debug mode.
  -v, --verbose         Level of output, default zero is no output

usage: make_prism_gnomad_file_seq.py [-h] [-e {mis,syn,complex,all}]
                                     [-m {overwrite,leave}]
                                     [-tmp_folder TMP_FOLDER]
                                     [-out_folder OUT_FOLDER] [-d] -uniprot
                                     UNIPROT -seq SEQ [-v]

Minimally requires a uniprot ID and a seq. Will make one prism file for the Ensembl transcript that matches this uniprot and sequence, if it exists. This is done to specifically make a file to a certain isoform/transcript.

required arguments:
  -uniprot UNIPROT      Uniprot ID to get transcripts and make prism_gnomad
                        files for
  -seq SEQ              Protein sequence of interest. Extract variants for the
                        first transcript that matches this.

optional arugments:
same as make_prism_gnomad_file_all

The reason it has to be this sucky and you can't just specific isoform 1 is that while you can get all Ensembl transcripts associated to a uniprot ID there is no guaranteed way to match a specific isoform to its Ensembl transcripts. So we compare the sequences and the transcript with the same seq as has been requested must be that isoform.

2.3. proteome wide data

now:
I later discovered that the approach below caused issues during merge because sometimes we are more interested in gnomad data on isoforms other than isoform 1, i.e. when there is clinvar data on another isoform. Therefore, I decided to instead make all possible files per uniprot ID in the human proteome. Re-ran make_prism_gnomad_file_all.py with only the uniprot ID as a bash oneliner:

while read line; do echo $line; python3 make_prism_gnomad_file_all.py -uniprot $line -e mis -v; done </storage1/shared/data/uniprot_datasets/human_proteome_UP000005640_9606.list

make_prism_gnomad_file_all.py includes a check if the output file already exists and skips if so (change this behavior by changing the mode to overwrite with -m overwrite). Therefore, you can just run this and it will only make the files that are missing.

--
first approach (has been superseeded!):
You can use this to make gnomad file only for isoform 1 but it turned out that was not so useful.

For generating this data I have downloaded the human proteome with one sequence per protein from uniprot from here: https://www.uniprot.org/proteomes/UP000005640 (click on "Download one protein sequence per gene")
resulting in: UP000005640_9606.fasta.gz

The fasta file was then converted into a white space separated file listing the full identifier, uniprotID and sequence like so:
>tr|A0A024R1R8|A0A024R1R8_HUMAN A0A024R1R8 MSSHEGGKKKALKQPKKQAKEMDEEEKAFKQKQKEEQKKLEVLKAKVVGKGPLATGGIKKSGKK

The following wrapper script reads in the space separated list and attempts to produce one prism_gnomad file per entry:
human_proteome_gnomad_wrapper.sh

Some proteins have no gnomad variants or are not on the main assembly. This means their genomic location is on what is called a patch. We do not have data for those locations since gnomad maps only to the main assembly. Therefore these protein do not have corresponding prism_gnomad files.
--


------
3. Creating prism clinvar files
------

3.1. Generating prism_clinvar files

You should not have to unless there are large scale changes.

3.1.0.
download a flat file with clinvar summary information to parse:
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz

unzip and sort the file by te symbol column (important!):
gunzip variant_summary.txt.gz
sort variant_summary.txt -k5 -t$'\t' > variant_summary.sort.symbol

other pages or interest
website:
https://www.ncbi.nlm.nih.gov/clinvar/
Significance levels:
https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/
Review status:
https://www.ncbi.nlm.nih.gov/clinvar/docs/review_status/
FTP directory:
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/

3.1.1.
run: 
make_prism_clinvar.py
, making sure line 526 is pointed to the flat file sorted by symbol! This is very important, otherwise it will not work properly. The file I've used is /storage1/hezscha/data/clinvar/variant_summary.sort.symbol . The script will go through the file and process it symbol by symbol, producing one prism_clinvar file per named refseq transcript (within the same symbol) and per uniprot ID associated to the refseq transcript, or to the HGNC if there are no uniprot IDs directly cross-referenced to the transcript.

------
4. Creating prism spliceAI files
------

4.1. Generating prism spliceAI files

You should not have to unless there are large scale changes. The files are created by parsing an entire chromosome vcf.

4.1.0.
go to https://github.com/Illumina/SpliceAI and obtain the vcf file spliceai_scores.masked.snv.hg38.vcf or similar

create separate file per one chromosome:
cd /storage0/shared/data/spliceAI
#if unzipped:
grep -P '^21\t' spliceai_scores.masked.snv.hg38.vcf > spliceai_scores.masked.snv.hg38.chromosome21.vcf
#otherwise:
gunzip -c spliceai_scores.masked.snv.hg38.vcf | grep -P '^Y\t' > spliceai_scores.masked.snv.hg38.chromosomeY.vcf
#add header
cat spliceai_scores.masked.snv.hg38.vcf.head spliceai_scores.masked.snv.hg38.chromosome21.vcf > tmp 
mv tmp spliceai_scores.masked.snv.hg38.chromosome21.vcf

2. filter for only positions predicted to be splice altering
https://github.com/Illumina/SpliceAI
according to their documentation the cutoffs are: 0.2 (high recall), 0.5 (recommended), and 0.8 (high precision). We used 0.2 to get the most info out.

I made this script to do the filtering:
python3 spliceAI_filter_vcf.py -chrom 7 -filter 0.2
#bash loop
for i in {1..5..1}; do; python3 spliceAI_filter_vcf.py -filter 0.2 -chrom $i; done

gzip the resulting file again

4.1.1. annotate the filtered vcf files with vep
if vep is not installed take a look at the gnomad section, it explains how to do it
annotated each filtered chromosome vcf with vep

#the command I used. I used sbatch on deic, an example job script can be found under vep_spliceai_chr21.sh
srun /groups/sbinlab/hezscha/.conda/envs/conda_vep/bin/vep -i /groups/sbinlab/hezscha/data/spliceAI/spliceai_scores.masked.snv.hg38.chromosome2.filter.spliceAI0.2.vcf.gz --cache --dir_cache /groups/sbinlab/hezscha/software/vep/cache_GRCh38_100 --fasta /groups/sbinlab/hezscha/software/vep/cache_GRCh38_100/homo_sapiens/100_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz -o /groups/sbinlab/hezscha/results/vep/spliceai_scores.masked.snv.hg38.chromosome2.filter.spliceAI0.2.vcf.VEP.cache100.bgz --vcf --offline --variant_class --sift b --polyphen b --distance 0 --hgvs --uniprot --domains --canonical --gene_phenotype --compress_output bgzip --force_overwrite --CACHE_VERSION 100 --gencode_basic

4.1.2. make prism files

python3 make_prism_spliceAI_files_v04.py -chrom Y -cutoff 0.2

------
5. Creating prism netsurfp files
------

Netsurfp predicts secondary structure in q3 and q8 categories as well as relative surface accessible area and disorder from sequence.

5.0.1 installing netsurfp
downloaded from:
https://services.healthtech.dtu.dk/cgi-bin/sw_request
tar -xvzf netsurfp-2.0.Any.tar.gz
prerequisits:
https://github.com/soedinglab/hh-suite
https://github.com/soedinglab/MMseqs2

#more infos, manuals, dbs
https://github.com/soedinglab/hh-suite/wiki#hh-suite-databases
https://wwwuser.gwdg.de/~compbiol/data/hhsuite/databases/hhsuite_dbs/
https://github.com/soedinglab/mmseqs2/wiki#downloading-databases

#I kept having dependency issues. What worked in the end was this order of installing things:
conda create -n netsurfp_tensor 
conda activate netsurfp_tensor
conda install -c conda-forge tensorflow=1.14

conda install -c conda-forge -c bioconda hhsuite 
#install mmseqs
conda install -c conda-forge -c bioconda mmseqs2
#for using the prism parser inside the env we need biopython and pandas
conda install -c anaconda biopython
conda install -c anaconda pandas
pip install pyyaml

cd netsurfp
mkdir build
#I removed the install_requires from setup.py (since both tensorflow and numpy(prereq of tensorflow) have been installed into the conda env manually) and re-ran installation
python setup.py develop --install-dir /storage1/hezscha/src/netsurfp/build

5.1. Running netsurfp

calling netsurfp only works after activating the conda env and standing in the directory it was installed into (full path to executable doesn't help):
conda activate netsurfp_tensor
cd /storage1/hezscha/src/netsurfp
netsurfp2 --help

#example run:
netsurfp2 --csv example/dhfr_human.csv --hhdb /storage1/shared/data/HH_databases/UniRef30_2020/UniRef30_2020_02 hhblits /storage1/hezscha/src/netsurfp/models/hhsuite.pb ../dhfr_human.fasta example/

#run over a list of uniprot IDs:
#given a a list of uniprot IDs and their isoform1 seq like so:
>tr|A0A024R1R8|A0A024R1R8_HUMAN A0A024R1R8 MSSHEGGKKKALKQPKKQAKEMDEEEKAFKQKQKEEQKKLEVLKAKVVGKGPLATGGIKKSGKK
>...
#this script generates temporary fasta files and submits them to netsurfp
python3 netsurfp_preds.py -in data/UP000005640_chromosome21_ID_seq.list -v -spliceAI -proteome
#the -proteome flag is important because it limits the prediction to IDs on the proteome list, /storage1/shared/data/uniprot_datasets/human_proteome_UP000005640_9606.list. The single chromosome lists contain more proteins than are in the 1 seq per gene list. They contain in total about 74K as oposed to 20K on the 1 seq per gene list!

#when done, run this to create prism files from the csv output files of netsurfp
python3 make_prism_netsurfp_files.py

------
6. Merging prism files
------

Merging is done with a wrapper around the prism parser. The procedure is roughly:
1. start with uniprot ID
2. choose target seq
3. pick files
4. attempt to resolve SNPs

Choosing the target seq happens automatically.Since each prism file can only have one sequence WT and the sequences in clinvar, uniprot, gnomad and spliceAI file can differ due to coming from different data sources as well as potentially being different isoforms, we need to choose one target sequence to merge all other sequences to. Variants that do not conform to the WT stated in the target sequence (the seq that will be in the metadata of the output prism file) are dropped during merge.

The target sequence is selected with the following hierarchy:
1. the sequence with the most non-VUS clinvar variants
2. the sequence with the most gnomad variants
3. the sequence with the most spliceAI variants
4. uniprot isoform 1

usage: merge_wrapper_v005.py [-h] [--typeslog] [-uniprot UNIPROT] [-ff FF]
                             [-out_folder OUT_FOLDER] [-swiss]
                             [-m {overwrite,leave}] [-v] [-allow]
                             [-check_uniprot] [-cleanup] [-fail]
                             [-min_id MIN_ID] [-min_cov MIN_COV]
                             [-snp_af SNP_AF] [--fill_invert]

optional arguments:
  -h, --help            show this help message and exit
  --typeslog            Write a log file listing the components of every merge
                        file for stats later.
  -uniprot UNIPROT      The uniprot ID for which to merge files. The
                        authorative sequence for merging is selected with
                        following hierarchy: 1. transcript with most clinvar
                        vars 2. transcript with most gnomad vars 3. transcript
                        with most spliceAI vars 4. uniprot isoform 1, unless
                        otherwise specified with -seq
  -ff FF, --fromfile FF
                        A file from which to read uniprot IDs (one per line).
  -out_folder OUT_FOLDER
                        Where the output files should be written. Default
                        location is /storage1/shared/data/prism_merge/uniprot[
                        0:2]/uniprot[2:4]/uniprot[4:6]/
  -swiss                Put a folder substructure into the designated output
                        folder (only for -out_folder, it's done automatically
                        for the default output folder).
  -m {overwrite,leave}  What do when the output file already exists. Leave
                        (default) or overwrite
  -v, --verbose         Level of output, default zero is no output
  -allow                Allow merging with files that do not have the exact
                        same sequence (but still the same uniprot ID).
  -check_uniprot        If there is no prism_uniprot file, check if this ID is
                        on the fail list (if it's not we can try to make a
                        prism_uniprot file to this ID).
  -cleanup              Reduce the resulting prism file to only rows that have
                        values in the gnomad/spliceAI/clinvar columns since
                        those are actual variants we have data for. Remove
                        rows, i.e. variants for which we only have uniprot
                        feature or netsurfp data.
  -fail                 Consult a file listing uniprot IDs for which we tried
                        to make this type of merge file before and failed.
                        Used to limit the amount of seq lookups we do to
                        uniprot since this is failing a lot now
  -min_id MIN_ID        Minimum identity (1.00 = 100 procent) to authoritative
                        sequence needed to qualify file for the merge. Default
                        0.8.
  -min_cov MIN_COV      Minimum coverage (1.00 = 100 procent) of authoritative
                        sequence needed to qualify file for the merge. Default
                        0.1.
  -snp_af SNP_AF        Variants with allel frequencies greater than or equal
                        to this will be regarded as SNPs, i.e. common
                        variants. Default 0.01.
  --fill_invert         If there are known SNP positions in the gnomad file
                        this will attempt to resolve mismatches with the
                        target seq by inverting those SNPs, also for uniprot
                        and netsurfp files.

additional information: The options -uniprot and -ff are mutually exlusive.
Use only one of them.


You should always use --fill_invert unless you have a specific reason not to. It will not affect the merging of files that do not have SNPs. It will attempt to resolved mismatches in WT that originate from common variants between the target seq and other sequences by matching the target seq WT. In the case of gnomad data the variant line is inverted (i.e. from T546A to A546T) to conform to the target seq WT and allel frequency estimates are calc as 1 - orig AF. In the case of uniprot nad netsurfp files which carry information per position, not per variant, the WT is simply switchen to conform to the target seq WT (i.e. from T546= to A546=). This only occurs if the variant has an allel freq of at least 1% in the gnomad data.

some typical ways of using the merge wrapper are: 

#whole proteome: 
output goes into /storage1/shared/data/prism_merge
!remember to cd into the directory where you want the log files before running since they are created in current working dir!

cd /storage1/hezscha/genome_proteome_map/results/proteome_merge_log_files
python3 /storage1/hezscha/genome_proteome_map/scripts/merge_wrapper_v005.py -ff /storage1/shared/data/uniprot_datasets/human_proteome_UP000005640_9606.list -allow -cleanup --fill_invert --typeslog

#marks genes: specific output folder w swiss structure

cd /storage1/hezscha/marks_disease_genes/data
(summary file and typeslog will be saved here)
python3 /storage1/hezscha/genome_proteome_map/scripts/merge_wrapper_v005.py -ff /storage1/shared/data/uniprot_datasets/marks_disease_genes.list -out_folder /storage1/hezscha/marks_disease_genes/data/clinvar_priority_merge_mismatch_match -allow -cleanup --fill_invert --typeslog -swiss -v

#single protein run
python3 /storage1/hezscha/genome_proteome_map/scripts/merge_wrapper_v005.py -uniprot Q9Y6B7 -out_folder /storage1/hezscha/marks_disease_genes/data/clinvar_priority_merge/ -allow -cleanup -vv -m overwrite --fill_invert -swiss


