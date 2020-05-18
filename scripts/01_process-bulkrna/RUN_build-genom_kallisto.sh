#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name=docker-index-kallisto
#SBATCH --time=0-04:00:00
#SBATCH --partition=exacloud
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --output=../../logs/slurm.%N.%j.out
#SBATCH --error=../../logs/slurm.%N.%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=leejor@ohsu.edu


dir=/home/groups/EllrottLab/cell-dissociation
ses=build-ref-genom

################
# Set up workspace
################
# Move to main dir
cd ${dir}
# make output structure (for this script and RUN.sh)
mkdir /mnt/scratch/5420
mkdir /mnt/scratch/5420/${ses}
mkdir /mnt/scratch/5420/${ses}/input
mkdir /mnt/scratch/5420/${ses}/input/reference_files
mkdir /mnt/scratch/5420/${ses}/output
mkdir /mnt/scratch/5420/${ses}/output/kallisto_index


# if not already there, get metadata for kallisto
cd /mnt/scratch/5420/${ses}/input/reference_files
# get reference file
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.transcripts.fa.gz

# Move to main dir
cd ${dir}


#########################
# RUN ONLY ONCE
# Build genome index
##########################
echo '##### Building kallisto genome index #####'
# kallisto index -i transcripts.idx /mnt/scratch/5420/${ses}/input/reference_files/gencode.v34.transcripts.fa.gz

# see https://hub.docker.com/r/zlskidmore/kallisto, version==0.46.0
sudo /opt/acc/sbin/exadocker pull zlskidmore/kallisto
sudo /opt/acc/sbin/exadocker run --rm -v /mnt/scratch/5420/${ses}:/tmp \
    zlskidmore/kallisto \
    kallisto index \
    -i /tmp/output/kallisto_index/transcripts_kallisto.idx /tmp/input/reference_files/gencode.v34.transcripts.fa.gz


###############
# Clean up workspace
###############
# copy files from scratch and clean up scratch
echo '##### Cleaning up workspace #####'
cp -r /mnt/scratch/5420/${ses}/output/* /home/groups/EllrottLab/cell-dissociation/raw-data
cp -r /mnt/scratch/5420/${ses}/input/reference_files /home/groups/EllrottLab/cell-dissociation/raw-data
rm -rf /mnt/scratch/5420/${ses}
