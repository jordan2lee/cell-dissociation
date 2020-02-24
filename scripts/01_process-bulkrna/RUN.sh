#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=docker_bulkpipeline
#SBATCH --time=0-00:10:00
#SBATCH --partition=exacloud
#SBATCH --ntasks=6
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
threads=1

################
# Set up workspace
################
# make output structure
mkdir /mnt/scratch/5420
mkdir /mnt/scratch/5420/${ses}
mkdir /mnt/scratch/5420/${ses}/output
mkdir /mnt/scratch/5420/${ses}/output/02a_fastqc
mkdir /mnt/scratch/5420/${ses}/output/03_trimmomatic
mkdir /mnt/scratch/5420/${ses}/output/04_fastqc-trimmed
mkdir /mnt/scratch/5420/${ses}/output/06_star
mkdir /mnt/scratch/5420/${ses}/output/07_ct_matrices
# move files to scratch, computing cluster specific
cp -r ${dir}/raw-data /mnt/scratch/5420/${ses} #move raw seq reads dir

#########################
# QC raw reads - fastqc
#########################
echo '##### QC raw reads #####'
# Read 1
echo 'starting read 1'
sudo /opt/acc/sbin/exadocker pull biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1 \
    fastqc \
    /tmp/raw-data/${read1} \
    -f fastq \
    -o /tmp/output/02a_fastqc
# Read 2
echo 'starting read 2'
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1 \
    fastqc \
    /tmp/raw-data/${read2} \
    -f fastq \
    -o /tmp/output/02a_fastqc

####################################
# Trim low qual reads - trimmomatic
####################################
echo '##### Trim low qual raw reads #####'
sudo /opt/acc/sbin/exadocker pull quay.io/biocontainers/trimmomatic:0.39--1
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    quay.io/biocontainers/trimmomatic:0.39--1 \
    trimmomatic PE \
    -threads ${threads} -phred33 \
    -trimlog /tmp/output/03_trimmomatic/LOG_${basename}.txt\
    /tmp/raw-data/${read1} \
    /tmp/raw-data/${read2} \
    /tmp/output/03_trimmomatic/${basename}_R1_trimmed.fastq.gz \
    /tmp/output/03_trimmomatic/${basename}_R1_trimmedUNPAIRED.fastq.gz \
    /tmp/output/03_trimmomatic/${basename}_R2_trimmed.fastq.gz \
    /tmp/output/03_trimmomatic/${basename}_R2_trimmedUNPAIRED.fastq.gz \
    LEADING:30 TRAILING:30\
    SLIDINGWINDOW:4:30 \
    MINLEN:90

#########################
# QC trimmed reads - fastqc
#########################
echo '##### QC trimmed reads #####'
# Read 1
echo 'starting read 1'
sudo /opt/acc/sbin/exadocker pull biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1 \
    fastqc \
    /tmp/output/03_trimmomatic/${basename}_R1_trimmed.fastq.gz \
    -f fastq \
    -o /tmp/output/04_fastqc-trimmed
# Read 2
echo 'starting read 2'
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1 \
    fastqc \
    /tmp/output/03_trimmomatic/${basename}_R2_trimmed.fastq.gz \
    -f fastq \
    -o /tmp/output/04_fastqc-trimmed

####################################
# Map Reads - output: unsorted BAM
####################################
echo '##### Mapping reads #####'
sudo /opt/acc/sbin/exadocker pull alexdobin/star:2.6.1d
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    alexdobin/star:2.6.1d \
    STAR \
    --runThreadN ${threads} \
    --readFilesCommand zcat \
    --genomeDir /tmp/raw-data/genome_ref \
    --readFilesIn /tmp/output/03_trimmomatic/${basename}_R1_trimmed.fastq.gz /tmp/output/03_trimmomatic/${basename}_R2_trimmed.fastq.gz \
    --outFileNamePrefix /tmp/output/06_star/${basename}_mapped-unsorted \
    --outSAMtype BAM Unsorted

###############
# BAM QC
###############
echo '##### QC BAM file #####'
# TBA
echo 'tba'

##############################
# Quantify reads: featureCounts
##############################
# featureCounts requires decompressed input files
zcat /mnt/scratch/5420/${ses}/raw-data/reference_files/Homo_sapiens.GRCh38.98.gtf.gz > /mnt/scratch/5420/${ses}/raw-data/reference_files/TEMP-ref.gtf
echo '##### Create count matrix #####'
sudo /opt/acc/sbin/exadocker pull alexgilgal/featurecount:latest
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    alexgilgal/featurecount:latest \
    featureCounts \
    -T ${threads} -F GTF \
    -g gene_id -t exon \
    -R \
    -s 1 \
    -a /tmp/raw-data/reference_files/TEMP-ref.gtf \
    -o /tmp/output/07_ct_matrices/${basename}_prelim_COUNTmatrix.txt /tmp/output/06_star/${basename}_mapped-unsortedAligned.out.bam

##############################
# Extract gene counts from count matrix
###############################
cut -f 1,7,8,9,10,11,12 /mnt/scratch/5420/${ses}/output/07_ct_matrices/${basename}_prelim_COUNTmatrix.txt > /mnt/scratch/5420/${ses}/output/07_ct_matrices/${basename}_COUNTmatrix.txt

###############
# Clean up workspace
###############
# copy files from scratch and clean up scratch
echo '##### Cleaning up workspace #####'
cp -r /mnt/scratch/5420/${ses}/output/* /home/groups/EllrottLab/cell-dissociation/data/01_process-bulkrna/
rm -rf /mnt/scratch/5420/${ses}
rm -rf /mnt/scratch/5420 #run only if not currently running other jobs
