# cell-dissociation
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



# Testing

## Virtual environment
```
virtual venv
. venv/bin/activate
```

## Cwltool Install
`pip install cwltool`

## Get test data
```
curl -o data/test/HBR_UHR_ERCC_ds_5pc.tar http://genomedata.org/rnaseq-tutorial/HBR_UHR_ERCC_ds_5pc.tar
tar -C data/test -xvf data/test/HBR_UHR_ERCC_ds_5pc.tar
```

## Run Fastqc Tool
```
cwltool --outdir data/01_process-bulkrna/ workflows/fastqc.cwl.yaml --fastq_file data/test/UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1.fastq.gz
```
