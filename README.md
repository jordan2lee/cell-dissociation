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
