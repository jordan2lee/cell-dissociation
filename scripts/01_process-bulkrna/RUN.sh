#!/bin/bash

cwl-runner --outdir ../../data/01_process-bulkrna/data_dump \
    bulk_prepro-workflow.cwl \
    bulk_prepro-workflow-job.yml
