#!/bin/bash
dir=/home/groups/EllrottLab/cell-dissociation
ses=round1
read1=TESTDATA--UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz
read2=TESTDATA--UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
basename=TESTDATA--UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22

################
# Set up workspace
################
# Move to main dir
cd ${dir}
# move files to scratch, computing cluster specific
mkdir /mnt/scratch/5420
mkdir /mnt/scratch/5420/${ses}
cp -r raw-data /mnt/scratch/5420/${ses} #move raw seq reads dir
# make output structure
mkdir /mnt/scratch/5420/${ses}/output
mkdir /mnt/scratch/5420/${ses}/output/02a_fastqc
mkdir /mnt/scratch/5420/${ses}/output/03_trimmomatic

#########################
# QC raw reads - fastqc
#########################
# Read 1
sudo /opt/acc/sbin/exadocker pull biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1 \
    fastqc \
    /tmp/raw-data/${read1} \
    -f fastq \
    -o /tmp/output/02a_fastqc

# Read 2
sudo /opt/acc/sbin/exadocker pull biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1 \
    fastqc \
    /tmp/raw-data/${read2} \
    -f fastq \
    -o /tmp/output/02a_fastqc


# # ####################################
# # # Trim low qual reads - trimmomatic
# # ####################################
# sudo /opt/acc/sbin/exadocker pull quay.io/biocontainers/trimmomatic:0.39--1
# sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
#     quay.io/biocontainers/trimmomatic:0.39--1 \
#     trimmomatic PE \
#     -threads 11 -phred33 \
#     -trimlog /tmp/output/LOG_${basename}.txt\
#     /tmp/raw-data/${read1} \
#     /tmp/raw-data/${read2} \
#     /tmp/output/03_trimmomatic/${read1}_trimmed.fastq.gz \
#     /tmp/output/03_trimmomatic/${read1}_trimmedUNPAIRED.fastq.gz \
#     /tmp/output/03_trimmomatic/${read2}_trimmed.fastq.gz \
#     /tmp/output/03_trimmomatic/${read2}_trimmedUNPAIRED.fastq.gz \
#     LEADING:30 TRAILING:30\
#     SLIDINGWINDOW:4:30 \
#     MINLEN:90
#
# # want to automate this step
# # QC raw reads - fastqc
# cwl-runner --outdir ../../data/01_process-bulkrna/data_dump/04_fastqc-trimmed \
#     workflows/fastqc_raw-workflow.cwl \
#     tools/fastqc_trimmedinputs.yml
#
# # Build genome indexes for alignment
# docker pull alexdobin/star:2.6.1d
# zcat /home/groups/EllrottLab/cell-dissociation/data/01_process-bulkrna/data_dump/06_star/metadata/ref_fa/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > /home/groups/EllrottLab/cell-dissociation/data/01_process-bulkrna/data_dump/06_star/metadata/ref_fa/TEMP-ref.fa
# zcat /home/groups/EllrottLab/cell-dissociation/data/01_process-bulkrna/data_dump/06_star/metadata/gtf/Homo_sapiens.GRCh38.98.gtf.gz > /home/groups/EllrottLab/cell-dissociation/data/01_process-bulkrna/data_dump/06_star/metadata/gtf/TEMP-ref.gtf
# docker run --rm -v /home/groups/EllrottLab/cell-dissociation:/tmp \
#     alexdobin/star:2.6.1d \
#     STAR \
#     --runThreadN 11 \
#     --runMode genomeGenerate \
#     --genomeDir /tmp/data/01_process-bulkrna/data_dump/06_star/genome_indexes \
#     --genomeFastaFiles /tmp/data/01_process-bulkrna/data_dump/06_star/metadata/ref_fa/TEMP-ref.fa \
#     --sjdbGTFfile /tmp/data/01_process-bulkrna/data_dump/06_star/metadata/gtf/TEMP-ref.gtf \
#     --sjdbOverhang 100 \
#     --genomeSAsparseD 10
# rm ../../data/01_process-bulkrna/data_dump/06_star/metadata/ref_fa/TEMP-ref.fa
#
# # Map Reads - output: unsorted BAM
# docker pull alexdobin/star:2.6.1d
# docker run --rm -v /home/groups/EllrottLab/cell-dissociation:/tmp \
#     alexdobin/star:2.6.1d \
#     STAR \
#     --runThreadN 11 \
#     --readFilesCommand zcat \
#     --genomeDir /tmp/data/01_process-bulkrna/data_dump/06_star/genome_indexes \
#     --readFilesIn /tmp/data/01_process-bulkrna/data_dump/03_trimmomatic/${read1}_trimmed.fastq.gz /tmp/data/01_process-bulkrna/data_dump/03_trimmomatic/${read2}_trimmed.fastq.gz \
#     --outFileNamePrefix /tmp/data/01_process-bulkrna/data_dump/06_star/mapped/${basename}_mapped \
#     --outSAMtype BAM Unsorted
#
# # BAM QC
#
#     # TBA
#
#
# # Quantify reads: featureCounts
# docker pull alexgilgal/featurecount:latest
# docker run --rm -v /home/groups/EllrottLab/cell-dissociation:/tmp \
#     alexgilgal/featurecount:latest \
#     featureCounts \
#     -T 10 -F GTF \
#     -g gene_id -t exon \
#     -R \
#     -s 1 \
#     -a /tmp/data/01_process-bulkrna/data_dump/06_star/metadata/gtf/TEMP-ref.gtf \
#     -o /tmp/data/01_process-bulkrna/data_dump/07_ct_matrices/${basename}_CTmatrix.txt  \
#     /tmp/data/01_process-bulkrna/data_dump/06_star/mapped/${basename}_mappedAligned.out.bam
# rm ../../data/01_process-bulkrna/data_dump/06_star/metadata/gtf/TEMP-ref.gtf


# copy files from scratch and clean up scratch
cp -r /mnt/scratch/5420/${ses}/output/* /home/groups/EllrottLab/cell-dissociation/data/01_process-bulkrna/
rm -rf /mnt/scratch/5420/${ses}
