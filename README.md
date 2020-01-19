# Purpose
Investigate artifacts produced by cell dissociation methods that confound scRNA-seq analysis

# Repository structure
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

# Getting started

### Virtual Environment
```
. venv3/bin/activate
```

### Cwltool Install
`pip install cwltool`

### Additional Installs
Install all other requirements in `requirements.txt`

### Get Bulk RNA-seq Example Data
```
curl -o data/test/HBR_UHR_ERCC_ds_5pc.tar http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
tar -C data/test -xvf data/test/HBR_UHR_ERCC_ds_5pc.tar
```

# Preprocessing: Bulk RNA-seq

`scripts/01_process-bulkrna/cwl_wrapper.sh`

This analysis includes:

1. Read Quality Control

The CWL tool `scripts/01_process-bulkrna/fastqc.cwl` runs **FastQC** in a docker container for read 1 and read 2 independently.

Requires manual inspection of two output summary files (`.html` `.zip`) to determine input parameters for read trimming.

2. Read Quality Trimming

 [WIP] The CWL tool `scripts/01_process-bulkrna/trimmomatic.cwl` runs **Trimomatic** in a docker container for read 1 and read 2.
