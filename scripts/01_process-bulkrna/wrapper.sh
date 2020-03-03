#!/bin/bash

###########
# PURPOSE: run preprocessing of bulk RNA
###########

# Patient 7319 - bulk rna
sbatch RUN.sh /home/groups/EllrottLab/cell-dissociation patient_7319 P2000996_02292020/FASTQ/SCC_7319_B2_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7319_B2_R2.fastq.gz SCC_7319_B2 5

# Patient 7320 - bulk rna
sbatch RUN.sh /home/groups/EllrottLab/cell-dissociation patient_7320 P2000996_02292020/FASTQ/SCC_7320_B2_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7320_B2_R2.fastq.gz SCC_7320_B2 5
