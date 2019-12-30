# General

This dir contains the scripts for processing raw sequence reads of bulk RNA-seq samples through generating count matrices and inspecting count matrices for odd occurances

Primary output of pipeline is in `data/01_process-bulkrna/`

# Prerequisites

Pipeline assumes that user has properly installed python3 and Java. And additional installations as indicated in `requirements.txt`

Development environment will occur in a virtual environment and finalized scripts will be deployed using Docker

Reference files for alignment (fasta and GTF) need to be manually downloaded prior to running pipeline


# Demultiplex reads & Adaptor trimming

# Quality filtering and trimming

### Raw read QC Report

Generate QC summary report of reads `fastqc.sh`

Output: 2 files in `data/readqc/` for each sequence file

### Filter and/or trim poor reads

Remove bp as per above QC summary report `trimmomatic.sh`

Output: `data/readtrim_and_filter/`

### Filtered/trimmed read QC report

Run again FastQC on the filtered/trimmed read files to see the new quality of the reads

Output: 2 files in `data/readtrim_and_filter/` for each sequence file

# Potential intermediate steps

K-mer filtering or K-mer normalization ?

# Alignment

STAR aligner

### Building genome indexes

Reference files downloaded: `Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz` and `Homo_sapiens.GRCh38.98.gtf.gz`

Ensembl homo sapien reference fasta and GTF file stored in `data/align/metadata/gtf/` and `data/align/metadata/ref_fa`

Reminder: Reference files need to be manually downloaded prior to running pipeline

Output: `data/align/genome_indexes/`

### Read mapping

Output: `data/align/mapped`

# QC BAM file

# Approximate strandedness

# Generate count matrices

# Extended

The analysis of this section will not be used for input of downstream cross project analysis. Instead this section allows us to examine our count matrices - and consider if hyperparameters in the above steps should be altered to improve quality of count matrices

### Normalize

### Dimensionality reduction

### QC for odd data
