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

### Virtual Environment
```
. venv3/bin/activate
```

### Cwltool Install
```
pip install cwltool
```

Rather than `sudo apt install cwltool` in case you have multiple CWL implementations

### Additional Installs
Install all other requirements in `requirements.txt`





# PREPROCESSING BULK RNA-SEQ WORKFLOW

Run workflow `cwl-runner --outdir ../../data/01_process-bulkrna/readqc/ bulk_prepro-workflow.cwl bulk_prepro-workflow-job.yml`

This workflow contains the following steps:

### Read Quality Control

The CWL tool `fastqc.cwl` runs **FastQC** in a docker container for read 1 and read 2 independently. For each input file, there are two output files (`.html` `.zip`)

Requires manual inspection of two output summary files to determine input parameters for read trimming.
