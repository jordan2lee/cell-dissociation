#!/bin/bash

############
# Purpose: run pre-processing, kallisto quant
# 	   for bulk-RNA data
############

# 1. Build genome index for downstream alignment
sbatch RUN_build-genom_kallisto.sh

# 2A. Run analytical pipeline (uses Kallisto) - patient 7320
sbatch RUN_pipeline.sh /home/groups/EllrottLab/cell-dissociation patient_7320_B2 P2000996_02292020/FASTQ/SCC_7320_B2_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7320_B2_R2.fastq.gz SCC_7320_B2 5

# 2B. Run analytical pipeline (uses Kallisto) - patient 7319
sbatch RUN_pipeline.sh /home/groups/EllrottLab/cell-dissociation patient_7319_B2 P2000996_02292020/FASTQ/SCC_7319_B2_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7319_B2_R2.fastq.gz SCC_7319_B2 5

# 3A. Run analytical pipeline (uses Kallisto) - GT patient 7320 - jobID 13019693
sbatch RUN_pipeline.sh /home/groups/EllrottLab/cell-dissociation patient_7320_GT P2000996_02292020/FASTQ/SCC_7320_GT_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7320_GT_R2.fastq.gz SCC_7320_GT 5

# 3B. Run analytical pipeline (uses Kallisto) - GT patient 7319 - jobID 13019694
sbatch RUN_pipeline.sh /home/groups/EllrottLab/cell-dissociation patient_7319_GT P2000996_02292020/FASTQ/SCC_7319_GT_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7319_GT_R2.fastq.gz SCC_7319_GT 5

# 4A. Run analytical pipeline (uses Kallisto) - normal adj patient 7320 - jobID 13019695
sbatch RUN_pipeline.sh /home/groups/EllrottLab/cell-dissociation patient_7320_N P2000996_02292020/FASTQ/SCC_7320_N_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7320_N_R2.fastq.gz SCC_7320_N 5
