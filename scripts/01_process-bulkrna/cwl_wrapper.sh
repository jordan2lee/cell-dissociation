#!/bin/bash

###########
# Readqc - run fastqc, raw reads
##########
# # Read 1
# cwltool --outdir ../../data/01_process-bulkrna/readqc/ \
#     fastqc.cwl \
#     --input_file ../../data/01_process-bulkrna/test-data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz \
#     --type 'fastq'
# # Read 2
# cwltool --outdir ../../data/01_process-bulkrna/readqc/ \
#     fastqc.cwl \
#     --input_file ../../data/01_process-bulkrna/test-data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz \
#     --type 'fastq'
###########
# Trim reads
###########
# Filter and/or trim as per above quality report
    #input1, input2, min_qual_leading, min_qual_trailing,out_basename
cwltool --outdir ../../data/01_process-bulkrna/readtrim_and_filter/ \
    trimmomatic.cwl \
    --fq1 ../../data/01_process-bulkrna/test-data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz \
    --fq2 ../../data/01_process-bulkrna/test-data/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz 
