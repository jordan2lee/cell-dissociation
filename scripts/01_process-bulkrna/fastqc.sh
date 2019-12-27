#!/bin/bash

##########
# Raw read QC - FastQC
#   App, input, input-format, output-dir
##########

raw_seqfile1=${1}
raw_seqfile2=${2}

echo '####################'
echo 'Raw read QC - FastQC version'
../../apps/FastQC/fastqc -version
echo '####################'

# Read 1
../../apps/FastQC/fastqc \
    ../../raw-data/test-data/${raw_seqfile1} -f fastq \
    -o ../../data/01_process-bulkrna/02_readqc/
# Read 2
../../apps/FastQC/fastqc \
    ../../raw-data/test-data/${raw_seqfile2} -f fastq \
    -o ../../data/01_process-bulkrna/02_readqc/
