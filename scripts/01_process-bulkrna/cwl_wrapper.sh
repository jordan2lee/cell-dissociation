#!/bin/bash

###########
# Readqc - run fastqc, raw reads
##########
# Read 1
cwltool --outdir ../../data/01_process-bulkrna/readqc/ \
    fastqc.cwl \
    --f1 ../../data/01_process-bulkrna/test-data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz
# Read 2
cwltool --outdir ../../data/01_process-bulkrna/readqc/ \
    fastqc.cwl \
    --f1 ../../data/01_process-bulkrna/test-data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
