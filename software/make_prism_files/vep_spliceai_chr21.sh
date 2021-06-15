#!/bin/bash
#SBATCH --job-name=21
#SBATCH --nodes=1
#SBATCH --time=2:00:00
#SBATCH --partition=sbinlab_ib
#SBATCH --mem-per-cpu=5G
#SBATCH --nodelist=node158

#load the conda env vep is installed into
#I'm trying to get it to be able to execute the conda command. Maybe I can try to source. Otherwise, put srun in front?
cd 
source .bashrc
conda activate conda_vep

#I'm setting distance to 0 since I don't unse the upstream or downstream modifications anyway
#omitted --species arg because Default = "homo_sapiens"
#I didn't specify the assembly or transcript set to use but since I put the path to the cache and used --offline it should use that one  
srun /groups/sbinlab/hezscha/.conda/envs/conda_vep/bin/vep -i /groups/sbinlab/hezscha/data/spliceAI/spliceai_scores.masked.snv.hg38.chromosome21.filter.spliceAI0.2.vcf.gz --cache --dir_cache /groups/sbinlab/hezscha/software/vep/cache_GRCh38_100 --fasta /groups/sbinlab/hezscha/software/vep/cache_GRCh38_100/homo_sapiens/100_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz -o /groups/sbinlab/hezscha/results/vep/spliceai_scores.masked.snv.hg38.chromosome21.filter.spliceAI0.2.vcf.VEP.cache100.bgz --vcf --offline --variant_class --sift b --polyphen b --distance 0 --hgvs --uniprot --domains --canonical --gene_phenotype --compress_output bgzip --force_overwrite --CACHE_VERSION 100 --gencode_basic

