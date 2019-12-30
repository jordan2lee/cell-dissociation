#!/bin/bash

# Man https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
echo '####################'
echo 'Align - STAR version'
../../apps/STAR-2.7.3a/bin/MacOSX_x86_64/STAR --version
echo '####################'

ref_fa=${1} #ensembl reference fasta file, Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gtf=${2} #ensembl reference gtf file, Homo_sapiens.GRCh38.98.gtf.gz
outdir=${3}

###############
# Create genome indexces
###############
# Create a temp unzipped file of the fasta files (STAR requires unzipped), create where running script
gzcat ${ref_fa} > ./TEMP-ref.fa
# create indexes
echo '1. Generate genome indexes files'
../../apps/STAR-2.7.3a/bin/MacOSX_x86_64/STAR \
    --runMode genomeGenerate \
    --runThreadN 4 \
    --genomeDir ${outdir} \
    --genomeFastaFiles TEMP-ref.fa \
    --sjdbGTFfile ${gtf} \
    --sjdbOverhang 100
# Clean up workspace
rm TEMP-ref.fa

###############
# Map reads to reference
###############
#tbc
