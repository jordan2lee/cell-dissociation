#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool
baseCommand: [fastqc, -o, ./] # name of command to run
requirements:
  InlineJavascriptRequirement: {}
hints:
  DockerRequirement:
    dockerPull: biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1


inputs:
  # Read 1 file
  input_file:
    type: File
    inputBinding:
      position: 1

  # format type, default passes 'fastq' to baseCommand
  type:
    type: string
    default: fastq
    inputBinding:
      position: 2
      prefix: "-f"


outputs:
  report:
    type: File
    outputBinding:
      glob: "*.html"
  summary:
    type: File
    outputBinding:
      glob: "*.zip"
