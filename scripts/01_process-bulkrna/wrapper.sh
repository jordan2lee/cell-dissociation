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
    ../../data/01_process-bulkrna/readqc/

# # Filter and/or trim as per above quality report
#     #input1, input2, min_qual_leading, min_qual_trailing,out_basename
# bash trimmomatic.sh \
#     ../../data/01_process-bulkrna/demultiplex_and_adpttrim/example_r1-da.fq.gz \
#     ../../data/01_process-bulkrna/demultiplex_and_adpttrim/example_r2-da.fq.gz \
#     39 39 \
#     ../../data/01_process-bulkrna/readtrim_and_filter/example-trim_39_39.fq.gz
#
# # See changes to fastqc
#     #input1, input2, outpath
# bash fastqc.sh \
#     fastq \
#     ../../data/01_process-bulkrna/readtrim_and_filter/example-trim_39_39_1P.fq.gz \
#     ../../data/01_process-bulkrna/readtrim_and_filter/example-trim_39_39_2P.fq.gz \
#     ../../data/01_process-bulkrna/readtrim_and_filter
#
# # Create genome indexes
#     # reference fa, reference gtf, existing outdir
# bash star.sh \
#     ../../data/01_process-bulkrna/align/metadata/ref_fa/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz \
#     ../../data/01_process-bulkrna/align/metadata/Homo_sapiens.GRCh38.98.gtf.gz \
#     ../../data/01_process-bulkrna/align/genome_indexes
