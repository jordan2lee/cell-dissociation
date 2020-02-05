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

### Meta data

Metadata for STAR
```
Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
```

Metadata for STAR and featureCounts
```
wget ftp://ftp.ensembl.org/pub/release-99/gtf/homo_sapiens/Homo_sapiens.GRCh38.98.gtf.gz
```
