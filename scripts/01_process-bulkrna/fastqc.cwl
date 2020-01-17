cwlVersion: v1.0
class: CommandLineTool

requirements:
  - class: DockerRequirement
    dockerPull: "biocontainers/fastqc:v0.11.8dfsg-2-deb_cv1"

baseCommand: [fastqc, -o, ./] # name of command to run

inputs:
  f1:
    type: File
    inputBinding:
      position: 1

  input_type:
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
