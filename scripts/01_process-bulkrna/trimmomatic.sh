#!/bin/bash

############
# Filter/trim raw reads
#   Man: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
############
input_f1=${1}
input_f2=${2}
qual_min_leading=${3}
qual_min_trailing=${4}
out_basename=${5}

echo '####################'
echo 'Filter/trim raw reads - Trimmomatic version'
# java -jar ../../apps/Trimmomatic-0.39/trimmomatic-0.39.jar -version
echo '####################'

# Trim for paired end reads
java -jar ../../apps/Trimmomatic-0.39/trimmomatic-0.39.jar \
    PE ${input_f1} ${input_f2} \
    -baseout ${out_basename} \
    -phred33 -threads 4 \
    LEADING:${qual_min_leading} TRAILING:${qual_min_trailing}
