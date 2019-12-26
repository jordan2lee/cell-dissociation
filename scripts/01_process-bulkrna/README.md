# General
---

This dir contains the scripts for processing raw sequence reads of bulk RNA-seq samples through generating count matrices and inspecting count matrices for odd occurances

Primary output of pipeline is in `data/01_process-bulkrna/`

# Prerequisites
---

Pipeline assumes that user has properly installed python3 and Java. And additional installations as indicated in `requirements.txt`

Development environment will occur in a virtual environment and finalized scripts will be deployed using Docker


# Raw read QC
---

Generate summary of reads using FastQC GUI

Output: 2 files in `01_readqc` for each sequence file

# Filter and/or trim poor reads


