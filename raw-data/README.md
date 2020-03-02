Original unchanged data stored here (example sequencing data)

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

Samples that are paired with P2000997_02242020. Recieved on March 2, 2020.

Fq files:

```
SCC_7319_B2_R1.fastq.gz
SCC_7319_B2_R2.fastq.gz
SCC_7319_GT_R1.fastq.gz
SCC_7319_GT_R2.fastq.gz
SCC_7320_B2_R1.fastq.gz
SCC_7320_B2_R2.fastq.gz
SCC_7320_GT_R1.fastq.gz
SCC_7320_GT_R2.fastq.gz
SCC_7320_N_R1.fastq.gz
SCC_7320_N_R2.fastq.gz
```

### external-u54

*Last updated 2/26/20*

All additional data outside of our project. Will be used for deconvolution models












# ftp files from medgenome

wget -r --no-parent --password <pwd>  https://portal.us.medgenome.com/<dir>/ --no-check-certificate --user <user>
