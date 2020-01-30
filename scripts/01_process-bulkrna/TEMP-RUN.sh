trim_read1=UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read1_trimmed.fastq.gz
trim_read2=UHR_Rep1_ERCC-Mix1_Build37-ErccTranscripts-chr22.read2_trimmed.fastq.gz

# Build genome indexes for alignment - output: unsorted BAM
docker pull alexdobin/star:2.6.1d
docker run --rm -v /home/ubuntu/cell-dissociation:/tmp \
    alexdobin/star:2.6.1d \
    STAR \
    --runThreadN 11 \
    --readFilesCommand zcat \
    --genomeDir /tmp/data/01_process-bulkrna/data_dump/06_star/genome_indexes \
    --readFilesIn /tmp/data/01_process-bulkrna/data_dump/04_fastqc-trimmed/${trim_read1} /tmp/data/01_process-bulkrna/data_dump/04_fastqc-trimmed/${trim_read2} \
    --outFileNamePrefix /tmp/data/01_process-bulkrna/data_dump/06_star/mapped \
    --outSAMtype BAM Unsorted

