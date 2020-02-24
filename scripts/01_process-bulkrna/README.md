# PREPROCESSING BULK RNA-SEQ WORKFLOW

This dir contains the scripts for processing raw sequence reads of bulk RNA-seq samples through generating count matrices and inspecting count matrices for odd occurances

### Overview

1. Demultiplex reads
2. Read quality control
3. Adapter and poor quality read trimming
4. Check read quality
5. Genome indexing (create alignment indices)
6. Alignment
7. Alignment File QC
8. Generate count matrices - from aligned reads (featureCounts) -> called `${basename}_COUNTmatrix.txt`


This workflow contains the following steps:

# Demultiplex Reads

*Will be added once have seq results from the Core*

# Read Quality Control (Docker)

The CWL tool `fastqc.cwl` runs **FastQC** in a docker container for read 1 and read 2 independently. For each input file, there are two output files (`.html` `.zip`)

Requires manual inspection of two output summary files to determine input parameters for read trimming.

# Adapter and Poor Quality Read Trimming (Docker)

**Adapter sequence trimming will be added in future iterations**

Trimmomatic 0.39. Paired trimming. Although virtually all adapter sequences should already have been trimmed, we will conduct a secondary pass to remove any remaining adapter sequences. Then low quality reads will be trimmed.

# Check Read Quality after trimming (Docker)

Run same CWL tool on output file from previous section. Currently input file names are hardcoded.

# Build Genome Indexes (Docker)

STAR aligner

# Align (Docker)

STAR aligner

Each FASTQ file has 5 output files, including unsorted by name BAM file

# Alignment File QC

**TBA**

Samtools

# Generate Count Matrices (Docker)

featureCounts (gene-level counting) and produces final count matrix `${basename}_COUNTmatrix.txt`

+ Parameters set for stranded, ignore multi-mapping reads, not current for paired end data

*Note* Transcript-level quantification is less accurate than gene-level quantification (e.g. salmon, RSEM, kallisto). Also transcript-level quantification has less clear biological interpretability. Less statistical power if split counts between isoforms.

# Additional Notes

Example CWL included in this scripts dir and will be implemented for publication. For now it is simply placed there as an example.
