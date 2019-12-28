#!/bin/bash

###############
# Quality filtering and trimming
###############
# Create raw read quality report
    #input1, input2, outpath
bash fastqc.sh \
    fastq \
    ../../raw-data/test-data/example_r1.fq.gz \
    ../../raw-data/test-data/example_r2.fq.gz \
    ../../data/01_process-bulkrna/02_readqc/

# Filter and/or trim as per above quality report
    #input1, input2, min_qual_leading, min_qual_trailing,out_basename
bash trimmomatic.sh \
    ../../data/01_process-bulkrna/01_demultiplex_and_adpttrim/example_r1-da.fq.gz \
    ../../data/01_process-bulkrna/01_demultiplex_and_adpttrim/example_r2-da.fq.gz \
    39 39 \
    ../../data/01_process-bulkrna/03_readtrim_and_filter/example-trim_39_39.fq.gz

# See changes to fastqc
    #input1, input2, outpath
bash fastqc.sh \
    fastq \
    ../../data/01_process-bulkrna/03_readtrim_and_filter/example-trim_39_39_1P.fq.gz \
    ../../data/01_process-bulkrna/03_readtrim_and_filter/example-trim_39_39_2P.fq.gz \
    ../../data/01_process-bulkrna/03_readtrim_and_filter
