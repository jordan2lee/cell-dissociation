#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=docker_pipeline
#SBATCH --time=0-12:00:00
#SBATCH --partition=exacloud
#SBATCH --ntasks=12
#SBATCH --cpus-per-task=1
#SBATCH --output=../../logs/slurm.%N.%j.out
#SBATCH --error=../../logs/slurm.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=leejor@ohsu.edu


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
# clean up any remaining dirs
rm -rf /mnt/scratch/5420
# move files to scratch, computing cluster specific
mkdir /mnt/scratch/5420
mkdir /mnt/scratch/5420/${ses}
cp -r raw-data /mnt/scratch/5420/${ses} #move raw seq reads dir
# make output structure
mkdir /mnt/scratch/5420/${ses}/output
mkdir /mnt/scratch/5420/${ses}/output/02a_fastqc
mkdir /mnt/scratch/5420/${ses}/output/03_trimmomatic
mkdir /mnt/scratch/5420/${ses}/output/04_fastqc-trimmed
mkdir /mnt/scratch/5420/${ses}/output/06_star/
mkdir /mnt/scratch/5420/${ses}/output/06_star/genome_ref
mkdir /mnt/scratch/5420/${ses}/output/06_star/mapped
mkdir /mnt/scratch/5420/${ses}/output/07_ct_matrices

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
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1 \
    fastqc \
    /tmp/raw-data/${read2} \
    -f fastq \
    -o /tmp/output/02a_fastqc

####################################
# Trim low qual reads - trimmomatic
####################################
sudo /opt/acc/sbin/exadocker pull quay.io/biocontainers/trimmomatic:0.39--1
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    quay.io/biocontainers/trimmomatic:0.39--1 \
    trimmomatic PE \
    -threads 11 -phred33 \
    -trimlog /tmp/output/LOG_${basename}.txt\
    /tmp/raw-data/${read1} \
    /tmp/raw-data/${read2} \
    /tmp/output/03_trimmomatic/${read1}_trimmed.fastq.gz \
    /tmp/output/03_trimmomatic/${read1}_trimmedUNPAIRED.fastq.gz \
    /tmp/output/03_trimmomatic/${read2}_trimmed.fastq.gz \
    /tmp/output/03_trimmomatic/${read2}_trimmedUNPAIRED.fastq.gz \
    LEADING:30 TRAILING:30\
    SLIDINGWINDOW:4:30 \
    MINLEN:90

#########################
# QC raw reads - fastqc
#########################
# Read 1
sudo /opt/acc/sbin/exadocker pull biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1 \
    fastqc \
    /tmp/output/03_trimmomatic/${read1}_trimmed.fastq.gz \
    -f fastq \
    -o /tmp/output/04_fastqc-trimmed
# Read 2
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1 \
    fastqc \
    /tmp/output/03_trimmomatic/${read2}_trimmed.fastq.gz \
    -f fastq \
    -o /tmp/output/04_fastqc-trimmed

##########################
# Build genome indexes for alignment
#       genomeSAsparseD=1 (default)
##########################
sudo /opt/acc/sbin/exadocker pull alexdobin/star:2.6.1d
zcat /mnt/scratch/5420/${ses}/raw-data/metadata-STAR/ref_fa/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > /mnt/scratch/5420/${ses}/raw-data/metadata-STAR/ref_fa/TEMP-ref.fa
zcat /mnt/scratch/5420/${ses}/raw-data/metadata-STAR/gtf/Homo_sapiens.GRCh38.98.gtf.gz > /mnt/scratch/5420/${ses}/raw-data/metadata-STAR/gtf/TEMP-ref.gtf
echo '### finished unzipping files ###'
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    alexdobin/star:2.6.1d \
    STAR \
    --runThreadN 11 \
    --runMode genomeGenerate \
    --genomeDir /tmp/output/06_star/genome_ref/ \
    --genomeFastaFiles /tmp/raw-data/metadata-STAR/ref_fa/TEMP-ref.fa \
    --sjdbGTFfile /tmp/raw-data/metadata-STAR/gtf/TEMP-ref.gtf \
    --sjdbOverhang 100 \
    --outFileNamePrefix /tmp/output/06_star/genome_ref/Log \
    --genomeSAsparseD 1

####################################
# Map Reads - output: unsorted BAM
####################################
##########
# testing
# cp -r /home/groups/EllrottLab/cell-dissociation/data/01_process-bulkrna/03_trimmomatic/* /mnt/scratch/5420/${ses}/output/03_trimmomatic/
# cp -r /home/groups/EllrottLab/cell-dissociation/data/01_process-bulkrna/06_star/* /mnt/scratch/5420/${ses}/output/06_star/
##########
sudo /opt/acc/sbin/exadocker pull alexdobin/star:2.6.1d
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    alexdobin/star:2.6.1d \
    STAR \
    --runThreadN 11 \
    --readFilesCommand zcat \
    --genomeDir /tmp/output/06_star/genome_ref/ \
    --readFilesIn /tmp/output/03_trimmomatic/${read1}_trimmed.fastq.gz /tmp/output/03_trimmomatic/${read2}_trimmed.fastq.gz \
    --outFileNamePrefix /tmp/output/06_star/mapped/${basename}_mapped-unsorted.bam \
    --outSAMtype BAM Unsorted

###############
# BAM QC
###############
# TBA

##############################
# Quantify reads: featureCounts
##############################
sudo /opt/acc/sbin/exadocker pull alexgilgal/featurecount:latest
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    alexgilgal/featurecount:latest \
    featureCounts \
    -T 10 -F GTF \
    -g gene_id -t exon \
    -R \
    -s 1 \
    -a /tmp/raw-data/metadata-STAR/gtf/TEMP-ref.gtf \
    -o /tmp/output/07_ct_matrices/${basename}_CTmatrix.txt  \
    /tmp/output/06_star/mapped/${basename}_mapped-unsorted.bam

###############
# Clean up workspace
###############
# copy files from scratch and clean up scratch
cp -r /mnt/scratch/5420/${ses}/output/* /home/groups/EllrottLab/cell-dissociation/data/01_process-bulkrna/
rm -rf /mnt/scratch/5420/${ses}
