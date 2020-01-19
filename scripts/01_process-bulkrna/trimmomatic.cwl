cwlVersion: v1.0
id: trimmomatic-0.38
label: trimmomatic-0.38
class: CommandLineTool

requirements:
  - class: DockerRequirement
    dockerPull: "quay.io/biocontainers/trimmomatic:0.38--1"

baseCommand: [java, -jar]

inputs:
  nthreads:
    type: int?
    default: 2
    inputBinding:
      prefix: -threads
      position: 1
    doc: number of cpu cores to be used
  fq1:
    type: File
    inputBinding:
      position: 2
    doc: FastQ file from next-generation sequencers
  fq2:
    type: File
    inputBinding:
      position: 3
    doc: FastQ file from next-generation sequencers

outputs: []
