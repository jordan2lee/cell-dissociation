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

# 1. DOWNLOAD RAW READ DATA

**PURPOSE: explanation of raw sequencing files location**

<details><summary>Click for details</summary><p>

# Where did the data come from?

Medgenome ftp. Files were moved into raw-data dir using `wget -r --no-parent --password <pwd>  https://portal.us.medgenome.com/<dir>/ --no-check-certificate --user <user>`

# Directory Structure - Patient Samples

### 1A. P2000997_02242020 - scRNA patient samples

*Last updated 2/26/20*

Samples paired with P2000996_02292020

2 patient samples. Recieved seq files on Feb. 25, 2020. This does NOT include the bulk RNA-seq data.
That data will be in a separate dir.

Fq files:

```
7319-Enz_S5_L002_I1.fastq.gz
7319-Enz_S5_L002_R1.fastq.gz
7319-Enz_S5_L002_R2.fastq.gz
7319-Mech_S4_L002_I1.fastq.gz
7319-Mech_S4_L002_R1.fastq.gz
7319-Mech_S4_L002_R2.fastq.gz
7320E_S7_L002_I1.fastq.gz
7320E_S7_L002_R1.fastq.gz
7320E_S7_L002_R2.fastq.gz
7320M_S6_L002_I1.fastq.gz
7320M_S6_L002_R1.fastq.gz
7320M_S6_L002_R2.fastq.gz
```

### 1B. P2000996_02292020 - bulk RNA patient samples

*Last updated 3/2/20*

Samples that are paired with P2000997_02242020. Received on March 2, 2020.

Fq files:

```
# Samples bulk rna
SCC_7319_B2_R1.fastq.gz     #patient 7319
SCC_7319_B2_R2.fastq.gz     #patient 7319
SCC_7320_B2_R1.fastq.gz     #patient 7320
SCC_7320_B2_R2.fastq.gz     #patient 7320

# Samples ground truth
SCC_7319_GT_R1.fastq.gz     #patient 7319
SCC_7319_GT_R2.fastq.gz     #patient 7319
SCC_7320_GT_R1.fastq.gz     #patient 7319
SCC_7320_GT_R2.fastq.gz     #patient 7319

# N fq
SCC_7320_N_R1.fastq.gz
SCC_7320_N_R2.fastq.gz
```

### external-u54

*Last updated 2/26/20*

All additional data outside of our project. Will be used for deconvolution models

</details>


# 2. PREPROCESSING BULK RNA-SEQ ANALYSIS

**PURPOSE: to take raw sequencing reads and convert to a gene count matrix**

Base script dir `scripts/01_process-bulkrna/`

<details><summary>Click for details</summary><p>

The dir `scripts/01_process-bulkrna/` contains the scripts for processing raw sequence reads of bulk RNA-seq samples through generating count matrices and inspecting count matrices for odd occurances

**The wrapper script allows for batch scripting tracking** `wrapper.sh`

*Step 1:*

1. Genome indexing (create alignment indices)

*Step 2:*

1. Demultiplex reads - *skip because seq core did this for us*
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

Job IDs:

+ SCC_7319_B2_R1.fastq.gz – job id 12515579
+ SCC_7319_B2_R2.fastq.gz – job id 12515580
+ SCC_7319_GT_R1.fastq.gz – job id 12515709
+ SCC_7320_GT_R1.fastq.gz – job id 12515710


### Demultiplex Reads

Sequencing core completed this for us. Therefore not included in our built pipeline.

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

# 3. PREPROCESSING SINGLE CELL RNA-SEQ ANALYSIS

**PURPOSE: to take raw sequencing reads and convert to a gene-cell count matrix**

Base script dir `scripts/02_process-scrna/`

<details><summary>Click for details</summary><p>

Scripts here

</details>


# 4. IDENTIFY DISSOCIATION SIGNATURES IN SC-RNA DATA ALONE (WIP)

**PURPOSE: to identify cell populations through gene marker analysis in single cell data. In otherwords, comparing profile of mechanical dissociation scRNA-seq to enzymatic scRNA-seq**

Base script dir `scripts/03_compare-bulk_scrna/signatures/`

<details><summary>Click for details</summary><p>

Scripts here

</details>

# 5. COMPARING PSEUDO-BULK TO BULK SAMPLE ANALYSIS (WIP)

**PURPOSE: to convert scRNA gene-cell count matrices into a pseudo-bulk count matrix. Then compare this pseduo-bulk data to the bulk RNA count matrix**

Base script dir `scripts/03_compare-bulk_scrna/scrna2bulk/`

<details><summary>Click for details</summary><p>

Analysis hardcoded in `dge.Rmd` and the rendered version `dge.html`

</details>

# 6. COMPARING PREDICTED CELL POPULATIONS TO SINGLE CELL SAMPLE ANALYSIS (WIP)

**PURPOSE: to convert bulk RNA count matrices into a pseudo-scRNA count matrix - this will be achieved by the cell-wise predictions of several deconvolution algorithms. Then compare this pseduo-scRNA data to the scRNA gene-cell count matrix**

Base script dir `scripts/03_compare-bulk_scrna/deconv/`

<details><summary>Click for details</summary><p>

Scripts here

</details>
