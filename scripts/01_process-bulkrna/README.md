# PREPROCESSING BULK RNA-SEQ WORKFLOW

This dir contains the scripts for processing raw sequence reads of bulk RNA-seq samples through generating count matrices and inspecting count matrices for odd occurances

Primary output of pipeline is in `data/01_process-bulkrna/`

Run workflow `RUN.sh`
```
cwl-runner --outdir ../../data/01_process-bulkrna/data_dump workflows/bulk_prepro-workflow.cwl tools/bulk_prepro-workflow-job.yml
```

This workflow contains the following steps:

# Read Quality Control

The CWL tool `fastqc.cwl` runs **FastQC** in a docker container for read 1 and read 2 independently. For each input file, there are two output files (`.html` `.zip`)

Requires manual inspection of two output summary files to determine input parameters for read trimming.

# Trim Reads

`trimmomatic.cwl` pulled from https://github.com/Duke-GCB/GGR-cwl/tree/master/v1.0/trimmomatic
