#!/bin/bash
read1=UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz
read2=UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
trimlog=UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22

# 1. Demultiplex reads
    #

# # # 2a. QC raw reads - fastqc
# cwl-runner --outdir ../../data/01_process-bulkrna/data_dump/02a_fastqc \
#     workflows/fastqc_raw-workflow.cwl \
#     tools/fastqc_raw-inputs.yml


# Trim low qual reads - trimmomatic
docker pull quay.io/biocontainers/trimmomatic:0.39--1
docker run --rm -v /home/ubuntu/cell-dissociation:/tmp \
    quay.io/biocontainers/trimmomatic:0.39--1 \
    trimmomatic PE \
    -threads 11 -phred33 \
    -trimlog /tmp/data/01_process-bulkrna/data_dump/02b_trimmomatic/trimlog_${trimlog}.txt\
    /tmp/data/01_process-bulkrna/test-data/${read1} \
    /tmp/data/01_process-bulkrna/test-data/${read2} \
    /tmp/data/01_process-bulkrna/data_dump/02b_trimmomatic/${read1}_trimmed.fastq.gz \
    /tmp/data/01_process-bulkrna/data_dump/02b_trimmomatic/${read2}_trimmed.fastq.gz \
    LEADING:28 TRAILING:28\
    SLIDINGWINDOW:4:28 \
    MINLEN:90
