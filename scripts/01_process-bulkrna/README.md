# PREPROCESSING BULK RNA-SEQ WORKFLOW

This dir contains the scripts for processing raw sequence reads of bulk RNA-seq samples through generating count matrices and inspecting count matrices for odd occurances

### Overview

1. Demultiplex reads
2. Read quality control
3. Adapter and poor quality read trimming
4. Check read quality
5. Genome indexing (create alignment indices)
6. Alignment
7. Generate count matrices - from aligned reads (samtools + HTSeqcount)


Primary output of pipeline is in `data/01_process-bulkrna/`

Run workflow `RUN.sh`
```
cwl-runner --outdir ../../data/01_process-bulkrna/data_dump workflows/bulk_prepro-workflow.cwl tools/bulk_prepro-workflow-job.yml
```

This workflow contains the following steps:

# Read Quality Control (CWL, Docker)

The CWL tool `fastqc.cwl` runs **FastQC** in a docker container for read 1 and read 2 independently. For each input file, there are two output files (`.html` `.zip`)

Requires manual inspection of two output summary files to determine input parameters for read trimming.

# Adapter and Poor Quality Read Trimming (Docker)

**Adapter sequence trimming will be added in future iterations**

Trimmomatic 0.39. Paired trimming. Although virtually all adapter sequences should already have been trimmed, we will conduct a secondary pass to remove any remaining adapter sequences. Then low quality reads will be trimmed.

# Check Read Quality after trimming

**Rewrite to automate yaml input file**

Run same CWL tool on output file from previous section. Currently input file names are hardcoded.

# Build Genome Indexes

STAR
