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
  -m {overwrite,leave}  What do when the output file already exists. Leave
                        (default) or overwrite
  -v, --verbose         Level of output.
  -fail FAIL            A file listing uniprot IDs for which extraction failed
                        before (i.e. because they have no features of interest) 
                        so we don't waste time querying for them.

At least -uniprot or -fromfile are required.
The default destination folder is /storage1/shared/data/prism_uniprot/[uniprotID[0:2]]/[uniprotID[2:4]]/[uniprotID[4:6]]/. 

For the purpose of merging with other prism files you need to first expand it to all variants with:

./FillVariants.py prism_uniprot_XXX_P00374.txt

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

2. Creating prism_gnomad files
------
A python parser to take a protein idenitfier and a gnomAD vcf file and produce a prism file for the protein

0. Setup step1 files
This was done already on the binf system, so just info unless things need redoing.

0.1. downloaded vcf files from gnomad:
https://gnomad.broadinstitute.org/downloads
https://gnomad.broadinstitute.org/faq

The FAQ details that: 
"The gnomAD v2.1 data set contains data from 125,748 exomes and 15,708 whole genomes, all mapped to the GRCh37/hg19 reference sequence. The gnomAD v3.1 data set contains 76,156 whole genomes (and no exomes), all mapped to the GRCh38 reference sequence. Most of the genomes from v2 are included in v3.1." 

I therefore used exome data from gnomad v2 and genome data from gnomad v3 to get the best available coverage in terms of both variants seen and their allel frequencies in the population.

0.2. annotation with Ensembl's Variant Effect Predictor

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
mkdir /groups/sbinlab/hezscha/software/vep/cache_GRCh38_100/
https://www.ensembl.org/info/docs/tools/vep/script/vep_cache.html#cache
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

0.3. temporary files 

To avoid reading in the huge vcf files every time, I made a directory with sub directories per chromosome that contains temporary files for each transcript. 

Specifically per transcript ID there are:
1 file with missense mutations in exome data: .mis.exomes.tmp
1 file with missense mutations in genome data: .mis.genomes.tmp
1 file with synonymous mutations in exome data: .syn.exomes.tmp
1 file with synonymous mutations in genome data: .syn.genomes.tmp
1 file with complex mutations in exome data: .complex.exomes.tmp
1 file with complex mutations in genome data: .complex.genomes.tmp

Files may not exist if there was no such data. I.e. if missense mutations were found for a certain transcript in the genome data but not in the exome data (i.e. due to the genomic region not being in the exome data), only the .mis.genomes.tmp will exist, not the .mis.exomes.tmp . The complex files are experimental and only exist for chrom Y.

These files are made with: 
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

1. Prerequisits 

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

2. File generation 

To generate prism compatible files with gnomad data, run make_prism_gnomad_file_seq.py

usage: make_prism_gnomad_file_seq.py [-h] [-e {mis,syn,complex,all}]
                                     [-m {overwrite,leave}]
                                     [-tmp_folder TMP_FOLDER]
                                     [-out_folder OUT_FOLDER] [-d] -uniprot
                                     UNIPROT -seq SEQ [-v]

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
  -uniprot UNIPROT      Uniprot ID to get transcripts and make prism_gnomad
                        files for
  -seq SEQ              Protein sequence of interest. Extract variants for the
                        first transcript that matches this.
  -v, --verbose         Level of output, default zero is no output


The reason it has to be this sucky and you can't just specific isoform 1 is that while you can get all Ensembl transcripts associated to a uniprot ID there is no guaranteed way to match a specific isoform to its Ensembl transcripts. So we compare the sequences and the transcript with the same seq as has been requested must be that isoform.

3. proteome wide data

For generating this data I have downloaded the human proteome with one sequence per protein from uniprot from here: https://www.uniprot.org/proteomes/UP000005640 (click on "Download one protein sequence per gene")
resulting in: UP000005640_9606.fasta.gz

The fasta file was then converted into a white space separated file listing the full identifier, uniprotID and sequence like so:
>tr|A0A024R1R8|A0A024R1R8_HUMAN A0A024R1R8 MSSHEGGKKKALKQPKKQAKEMDEEEKAFKQKQKEEQKKLEVLKAKVVGKGPLATGGIKKSGKK

The following wrapper script reads in the space separated list and attempts to produce one prism_gnomad file per entry:
human_proteome_gnomad_wrapper.sh

Some proteins have no gnomad variants or are not on the main assembly. This means their genomic location is on what is called a patch. We do not have data for those locations since gnomad maps only to the main assembly. Therefore these protein do not have corresponding prism_gnomad files.
