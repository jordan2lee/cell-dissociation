# General

This dir contains the scripts for processing raw sequence reads of bulk RNA-seq samples through generating count matrices and inspecting count matrices for odd occurances

Primary output of pipeline is in `data/01_process-bulkrna/`

# Prerequisites

Pipeline assumes that user has properly installed python3 and Java. And additional installations as indicated in `requirements.txt`

Development environment will occur in a virtual environment and finalized scripts will be deployed using Docker

# Demultiplex reads & Adaptor trimming

# Quality filtering and trimming

### Raw read QC

Generate QC summary report of reads `fastqc.sh`

Output: 2 files in `data/02_readqc/` for each sequence file

### Filter and/or trim poor reads

Remove bp as per above QC summary report

Output: `03_readtrim_and_filter`

# K-mer filtering

# K-mer normalization

# Alignment

# QC BAM file

# Approximate strandedness

# Generate count matrices

# Extended

The analysis of this section will not be used for input of downstream cross project analysis. Instead this section allows us to examine our count matrices - and consider if hyperparameters in the above steps should be altered to improve quality of count matrices

### Normalize

### Dimensionality reduction

### QC for odd data
