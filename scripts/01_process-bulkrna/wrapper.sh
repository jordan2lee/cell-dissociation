#!/bin/bash

###############
# Quality filtering and trimming
###############
# Create raw read quality report
bash fastqc.sh example_r1.fq.gz example_r2.fq.gz

# # Filter and/or trim as per above quality report
# bash trimmomatic.sh ../../data/01_process-bulkrna/01_demultiplex_and_adpttrim ../../data/01_process-bulkrna/03_readtrim_and_filter
