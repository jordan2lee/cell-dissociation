#!/bin/bash

##########
# Raw read QC - FastQC
#   App, input, input-format, output-dir
##########
input_type=${1}
raw_seqfile1=${2}
raw_seqfile2=${3}
outpath=${4}

echo '####################'
echo 'Raw read QC - FastQC version'
../../apps/FastQC/fastqc -version
echo '####################'

# Read 1
../../apps/FastQC/fastqc \
    ${raw_seqfile1} -f ${input_type} \
    -o ${outpath}
# Read 2
../../apps/FastQC/fastqc \
    ${raw_seqfile2} -f ${input_type} \
    -o ${outpath}
