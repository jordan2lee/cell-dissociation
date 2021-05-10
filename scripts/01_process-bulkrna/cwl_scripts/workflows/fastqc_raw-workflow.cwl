#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
requirements:
  - class: DockerRequirement
    dockerPull: "biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1"





inputs:
  input_file_r1:
    doc: Input read1 fastq file of raw sequencing
    type: File
  input_file_r2:
    doc: Input read2 fastq file of raw sequencing
    type: File





outputs:
  # RAW SEQ - QC
  rawreadqc_r1_REPORT:
    doc: HTML file of FastQC report of raw sequencing read1
    type: File
    outputBinding:
      glob: "*.html"
    outputSource: rawreadqc_r1/report
  rawreadqc_r1_SUMMARY:
    doc: ZIP of FastQC report of raw sequencing read1
    type: File
    outputBinding:
      glob: "*.zip"
    outputSource: rawreadqc_r1/summary

  rawreadqc_r2_REPORT:
    doc: HTML file of FastQC report of raw sequencing read2
    type: File
    outputBinding:
      glob: "*.html"
    outputSource: rawreadqc_r2/report
  rawreadqc_r2_SUMMARY:
    doc: ZIP of FastQC report of raw sequencing read2
    type: File
    outputBinding:
      glob: "*.zip"
    outputSource: rawreadqc_r2/summary




steps:
  # RAW SEQ - QC
  rawreadqc_r1:
    run: ../tools/fastqc.cwl
    in:
      input_file: input_file_r1
    out: [report, summary]

  rawreadqc_r2:
    run: ../tools/fastqc.cwl
    in:
      input_file: input_file_r2
    out: [report, summary]
