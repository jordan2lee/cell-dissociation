#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=docker-quant-kallisto
#SBATCH --time=0-05:30:00
#SBATCH --partition=exacloud
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=1
#SBATCH --output=../../logs/slurm.%N.%j.out
#SBATCH --error=../../logs/slurm.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=leejor@ohsu.edu


dir=${1}
ses=${2}
read1=${3}
read2=${4}
basename=${5}
threads=${6}

# dir=/home/groups/EllrottLab/cell-dissociation
# ses=patient_7320_B2
# read1=P2000996_02292020/FASTQ/SCC_7320_B2_R1.fastq.gz
# read2=P2000996_02292020/FASTQ/SCC_7320_B2_R2.fastq.gz
# basename=SCC_7320_B2
# threads=5

# dir=/home/groups/EllrottLab/cell-dissociation
# ses=patient_7319_B2
# read1=P2000996_02292020/FASTQ/SCC_7319_B2_R1.fastq.gz
# read2=P2000996_02292020/FASTQ/SCC_7319_B2_R2.fastq.gz
# basename=SCC_7319_B2
# threads=5


################
# Set up workspace
################
# make output structure
echo '##### Creating output structure #####'
mkdir /mnt/scratch/5420
mkdir /mnt/scratch/5420/${ses}
mkdir /mnt/scratch/5420/${ses}/output
mkdir /mnt/scratch/5420/${ses}/output/02a_fastqc
mkdir /mnt/scratch/5420/${ses}/output/03_trimmomatic
mkdir /mnt/scratch/5420/${ses}/output/04_fastqc-trimmed
mkdir /mnt/scratch/5420/${ses}/output/05_kallisto
mkdir /mnt/scratch/5420/${ses}/output/05_kallisto/${basename}

# move files to scratch, computing cluster specific
echo '##### Transfering input files to scratch #####'
cp -r ${dir}/raw-data /mnt/scratch/5420/${ses} #move raw seq reads dir

#############
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


#########################
# Quant
##########################
echo '##### Quantify abundances of the transcripts - Kallisto #####'
# see https://hub.docker.com/r/zlskidmore/kallisto, version==0.46.0
sudo /opt/acc/sbin/exadocker pull zlskidmore/kallisto
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    zlskidmore/kallisto \
    kallisto quant \
    -i /tmp/raw-data/kallisto_index/transcripts_kallisto.idx \
	--rf-stranded \
	-t ${threads} \
	-o /tmp/output/05_kallisto/${basename} \
	/tmp/output/03_trimmomatic/${basename}_R1_trimmed.fastq.gz /tmp/output/03_trimmomatic/${basename}_R2_trimmed.fastq.gz


###############
# BAM QC
###############
echo '##### QC BAM file #####'
# TBA
echo 'tba'


###############
# Clean up workspace
###############
# copy files from scratch and clean up scratch
# echo '##### Cleaning up workspace #####'
cp -r /mnt/scratch/5420/${ses}/output/* /home/groups/EllrottLab/cell-dissociation/data/01_process-bulkrna/
# rm -rf /mnt/scratch/5420/${ses}
# rm -rf /mnt/scratch/5420 #run only if not currently running other jobs
