#!/bin/bash

############
# Filter/trim raw reads
#   Man: http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/TrimmomaticManual_V0.32.pdf
############
inpath=${1}
outpath=${2}

echo '####################'
echo 'Filter/trim raw reads - Trimmomatic version'
# java -jar ../../apps/Trimmomatic-0.39/trimmomatic-0.39.jar -version
echo '####################'

java -jar ../../apps/Trimmomatic-0.39/trimmomatic-0.39.jar \
    PE ${inpath}/example_r1-da.fq.gz ${inpath}/example_r2-da.fq.gz \
    -baseout ${outpath}/example-trim.fq.gz \
    LEADING:3 TRAILING:3 MINLEN:36
