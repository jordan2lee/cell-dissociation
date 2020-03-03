# Purpose
Investigate artifacts produced by cell dissociation methods that confound scRNA-seq analysis

# Repository Structure
Set up workspace by running `./setup-dir.sh` which will create the following nested structure:

`scripts/`
+ scripts/01_process-bulkrna
+ scripts/02_process-scrna
+ scripts/03_compare-bulk_scrna

`data/`
+ data/01_process-bulkrna
+ data/02_process-scrna
+ data/03_compare-bulk_scrna

`raw-data/`

`sandbox/`

`notes/`


# Getting Started

### Cwltool

Pipeline largely relies on docker and cwltools to be installed. Make sure to have these installed prior to running this pipeline. Some examples include `pip install cwltool` and rather than `sudo apt install cwltool` in case you have multiple CWL implementations

### Reference meta files

Building a reference genome requires reference Ensembl files and this pipeline assumes you have already retrieved these. The existing pipeline uses the following files (but if you would like to use different files then be sure to update the scripts directly - as these files have been hardcoded)

+ Fasta for STAR `Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`

+ GTF for STAR and featureCounts `wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz`

# 1. Preprocessing Bulk RNA-seq Samples

**PURPOSE: to take raw sequencing reads and convert to a gene count matrix**


# PREPROCESSING BULK RNA-SEQ WORKFLOW

<details><summary>Click for more info</summary><p>

The dir `scripts/01_process-bulkrna/` contains the scripts for processing raw sequence reads of bulk RNA-seq samples through generating count matrices and inspecting count matrices for odd occurances

**The wrapper script allows for batch scripting tracking** `wrapper.sh`

*Step 1:*

1. Genome indexing (create alignment indices)

*Step 2:*

1. Demultiplex reads
2. Read quality control
3. Adapter and poor quality read trimming
4. Check read quality
5. Alignment
6. Alignment File QC
7. Generate count matrices - from aligned reads (featureCounts) -> called `${basename}_COUNTmatrix.txt`

## Step 1: `RUN_build-genom.sh`

### Build Genome Indexes (Docker)

STAR aligner. Job ID 12460222

## Step 2: `RUN.sh` via `wrapper.sh`

Job ID 12467117

### Demultiplex Reads

**[TODO] Will be added once have seq results from the Core**

### Read Quality Control (Docker)

Produce read quality reports of raw seq files (that have been demultiplexed).

Requires manual inspection of two output summary files to determine input parameters for read trimming.

### Adapter and Poor Quality Read Trimming (Docker)

**[TODO] Adapter sequence trimming will be added once have seq results from the Core**

Trimmomatic 0.39. Paired trimming. Although virtually all adapter sequences should already have been trimmed, we will conduct a secondary pass to remove any remaining adapter sequences. Then low quality reads will be trimmed.

### Check Read Quality after trimming (Docker)

Produce read quality reports of trimmed seq files.

### Align (Docker)

STAR aligner

Each FASTQ file has 5 output files, including unsorted by name BAM file

### Alignment File QC

**[TODO] Adapter sequence trimming will be added once have seq results from the Core**

Samtools

### Generate Count Matrices (Docker)

featureCounts (gene-level counting) and produces final count matrix `${basename}_COUNTmatrix.txt`

+ Parameters set for stranded, ignore multi-mapping reads, not current for paired end data

*Note* Transcript-level quantification is less accurate than gene-level quantification (e.g. salmon, RSEM, kallisto). Also transcript-level quantification has less clear biological interpretability. Less statistical power if split counts between isoforms.

### Additional Notes

Example CWL included in this scripts dir and will be implemented for publication. For now it is simply placed there as an example.

</details>
