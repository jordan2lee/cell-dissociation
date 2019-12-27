#!/bin/bash

##########
# FastQC
#   App, input, input-format, output-dir
##########
# Read 1
../../apps/FastQC/fastqc \
    ../../raw-data/test-data/example_r1.fq -f fastq \
    -o ../../data/01_process-bulkrna/01_readqc/
# Read 2
../../apps/FastQC/fastqc \
    ../../raw-data/test-data/example_r2.fq -f fastq \
    -o ../../data/01_process-bulkrna/01_readqc/
