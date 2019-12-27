#!/bin/bash

###############
# Quality filtering and trimming
###############
# Create raw read quality report
bash fastqc.sh example_r1.fq example_r2.fq
# Filter and/or trim as per above quality report
bash trimmomatic.sh
