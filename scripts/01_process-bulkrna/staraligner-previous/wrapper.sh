###########
# PURPOSE: run preprocessing of bulk RNA
###########

# # Patient 7319 - bulk rna
# sbatch RUN.sh /home/groups/EllrottLab/cell-dissociation patient_7319_B2 P2000996_02292020/FASTQ/SCC_7319_B2_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7319_B2_R2.fastq.gz SCC_7319_B2 5
# # Patient 7320 - bulk rna
# sbatch RUN.sh /home/groups/EllrottLab/cell-dissociation patient_7320_B2 P2000996_02292020/FASTQ/SCC_7320_B2_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7320_B2_R2.fastq.gz SCC_7320_B2 5
#
#
# # Patient 7319 - ground truth
# sbatch RUN.sh /home/groups/EllrottLab/cell-dissociation patient_7319_GT P2000996_02292020/FASTQ/SCC_7319_GT_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7319_GT_R2.fastq.gz SCC_7319_GT 5
# # Patient  7320 - ground truth
# sbatch RUN.sh /home/groups/EllrottLab/cell-dissociation patient_7320_GT P2000996_02292020/FASTQ/SCC_7320_GT_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7320_GT_R2.fastq.gz SCC_7320_GT 5
#

# Patient normal
sbatch RUN.sh /home/groups/EllrottLab/cell-dissociation patient_7319_N P2000996_02292020/FASTQ/SCC_7319_N_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7319_N_R2.fastq.gz SCC_7319_N 5
# Patient normal
sbatch RUN.sh /home/groups/EllrottLab/cell-dissociation patient_7320_N P2000996_02292020/FASTQ/SCC_7320_N_R1.fastq.gz P2000996_02292020/FASTQ/SCC_7320_N_R2.fastq.gz SCC_7320_N 5
