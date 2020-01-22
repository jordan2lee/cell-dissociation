#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: Workflow
requirements:
  - class: DockerRequirement
    dockerPull: "biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1"




inputs:
  input_file_r1:
    type: File
  input_file_r2:
    type: File
outputs:
  rawreadqc_r1_REPORT:
    type: File
    outputBinding:
      glob: "*.html"
    outputSource: rawreadqc_r1/report
  rawreadqc_r1_SUMMARY:
    type: File
    outputBinding:
      glob: "*.zip"
    outputSource: rawreadqc_r1/summary
    
  rawreadqc_r2_REPORT:
    type: File
    outputBinding:
      glob: "*.html"
    outputSource: rawreadqc_r2/report
  rawreadqc_r2_SUMMARY:
    type: File
    outputBinding:
      glob: "*.zip"
    outputSource: rawreadqc_r2/summary



steps:
  rawreadqc_r1:
    run: fastqc.cwl
    in:
      input_file: input_file_r1
    out: [report, summary]

  rawreadqc_r2:
    run: fastqc.cwl
    in:
      input_file: input_file_r2
    out: [report, summary]
