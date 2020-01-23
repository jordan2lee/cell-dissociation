#!/bin/bash

cwl-runner --outdir ../../data/01_process-bulkrna/data_dump \
    workflows/bulk_prepro-workflow.cwl \
    tools/bulk_prepro-workflow-job.yml
