#!/bin/bash
read1=UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz
read2=UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2.fastq.gz
trimlog=UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22
# trim_read1=UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1_trimmed.fastq.gz
# trim_read2=UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2_trimmed.fastq.gz


# # QC raw reads - fastqc
# cwl-runner --outdir ../../data/01_process-bulkrna/data_dump/02a_fastqc \
#     workflows/fastqc_raw-workflow.cwl \
#     tools/fastqc_raw-inputs.yml
#
#
# # Trim low qual reads - trimmomatic
# docker pull quay.io/biocontainers/trimmomatic:0.39--1
# docker run --rm -v /home/ubuntu/cell-dissociation:/tmp \
#     quay.io/biocontainers/trimmomatic:0.39--1 \
#     trimmomatic PE \
#     -threads 11 -phred33 \
#     -trimlog /tmp/data/01_process-bulkrna/data_dump/03_trimmomatic/LOG_${trimlog}.txt\
#     /tmp/data/01_process-bulkrna/test-data/${read1} \
#     /tmp/data/01_process-bulkrna/test-data/${read2} \
#     /tmp/data/01_process-bulkrna/data_dump/03_trimmomatic/${read1}_trimmed.fastq.gz \
#     /tmp/data/01_process-bulkrna/data_dump/03_trimmomatic/${read1}_trimmedUNPAIRED.fastq.gz \
#     /tmp/data/01_process-bulkrna/data_dump/03_trimmomatic/${read2}_trimmed.fastq.gz \
#     /tmp/data/01_process-bulkrna/data_dump/03_trimmomatic/${read2}_trimmedUNPAIRED.fastq.gz \
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
# zcat /home/ubuntu/cell-dissociation/data/01_process-bulkrna/data_dump/06_star/metadata/ref_fa/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > /home/ubuntu/cell-dissociation/data/01_process-bulkrna/data_dump/06_star/metadata/ref_fa/TEMP-ref.fa
# zcat /home/ubuntu/cell-dissociation/data/01_process-bulkrna/data_dump/06_star/metadata/gtf/Homo_sapiens.GRCh38.98.gtf.gz > /home/ubuntu/cell-dissociation/data/01_process-bulkrna/data_dump/06_star/metadata/gtf/TEMP-ref.gtf
# docker run --rm -v /home/ubuntu/cell-dissociation:/tmp \
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
# rm ../../data/01_process-bulkrna/data_dump/06_star/metadata/gtf/TEMP-ref.gtf
#
# Map Reads - output: unsorted BAM
docker pull alexdobin/star:2.6.1d
docker run --rm -v /home/ubuntu/cell-dissociation:/tmp \
    alexdobin/star:2.6.1d \
    STAR \
    --runThreadN 11 \
    --readFilesCommand zcat \
    --genomeDir /tmp/data/01_process-bulkrna/data_dump/06_star/genome_indexes \
    --readFilesIn /tmp/data/01_process-bulkrna/data_dump/03_trimmomatic/${read1}_trimmed.fastq.gz /tmp/data/01_process-bulkrna/data_dump/03_trimmomatic/${read2}_trimmed.fastq.gz \
    --outFileNamePrefix /tmp/data/01_process-bulkrna/data_dump/06_star/mapped \
    --outSAMtype BAM Unsorted
