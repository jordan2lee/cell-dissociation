#!/bin/bash

# Man https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

echo '####################'
echo 'Align - STAR version'
../../apps/STAR-2.7.3a/bin/MacOSX_x86_64/STAR --version
echo '####################'

echo '1. Generate genome indexes files'
../../apps/STAR-2.7.3a/bin/MacOSX_x86_64/STAR \
    --runMode genomeGenerate \
    --runThreadN 4 \
    --genomeDir ../../data/01_process-bulkrna/04_align/genome_indexes \
    --genomeFastaFiles ../../data/01_process-bulkrna/03_readtrim_and_filter/example-trim_39_39_1P.fq.gz ../../data/01_process-bulkrna/03_readtrim_and_filter/example-trim_39_39_2P.fq.gz \
    --sjdbGTPfile <GTP file> \
    --sjdbOverhang 100
