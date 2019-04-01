#!/usr/bin/env cwl-runner
cwlVersion: v1.0
class: CommandLineTool

label: BAMscale-scale
doc: Scale one or multiple BAM files

requirements:
  InlineJavascriptRequirement: {}

hints:
  - $import: bamscale.yml

inputs:
  l:
    type: string?
    inputBinding:
      position: 1
      prefix: -l
    doc: |
      Sequencing type to be used. Can be: single, paired, and auto (default: autodetect)
  f:
    type: string?
    inputBinding:
      position: 1
      prefix: -f
    doc: |
      Compute coverage using fragments instead of reads (default: no)
  a:
    type: int?
    inputBinding:
      position: 1
      prefix: -a
    doc: |
      Fragment size to be used to extend single-end library reads
  y:
    type: string?
    inputBinding:
      position: 2
      prefix: -y
    doc: |
      Type of normalization. (default: base)
      If no normalization is needed, set '--scale no' argument, the program will disregard this option.
      Options:
        1) reads: No. of mapped reads/fragments
        2) base: Sum of per-base coverage of reads/fragments
  k:
    type: string?
    inputBinding:
      position: 2
      prefix: -k
    doc: |
      Method to scale samples together. (default: genome)
      Options are:
        1) no: no scaling, just calculate coverage
        2) smallest: scale reads to smallest library (multiple-samples only)
        3) genome: scale samples to 1x genome coverage (only possible with 'base' normalization type)
  r:
    type: string?
    inputBinding:
      position: 1
      prefix: -r
    doc: |
      Operation to perform when scaling samples. Default: scaled
      Options are:
        1) scaled: output scaled tracks.
        2) unscaled: do not scale files in any way.
        2) log2: log2 transform against first BAM file.
        3) ratio: coverage ratio against first BAM file.
        4) subtract: subtract coverage against first BAM file.
        5) rfd: OK-seq RFD calculation
  z:
    type: int?
    inputBinding:
      position: 2
      prefix: -z
    doc: |
      Size of bins for output bigWig/bedgraph generation (default: 5)
  e:
    type: int?
    inputBinding:
      position: 3
      prefix: -e
    doc: |
      Compute sequencing coverage from BAM file. (default: '1', count reads while parsing BAM)
      Options are:
        1) 0: use reads in index (only if normalization is set to 'reads')
        2) 1: count reads while parsing BAM(s)
      WARNING: this option is only useful when 'reads' are used for normalization
  c:
    type: File?
    inputBinding:
      position: 3
      prefix: -c
    doc: |
      Input file with list of chromosomes to blacklist when computing coverage for normalization
  u:
    type: int?
    inputBinding:
      position: 3
      prefix: -u
    doc: |
      BED file with regions to subtract when computing coverage for normalization
      These coordinates should not overlap so reads are not counted multiple times
  j:
    type: int?
    inputBinding:
      position: 3
      prefix: -j
    doc: |
      Smoothen signal by calculating mean of N bins flanking both sides of each bin (default: 0)
      If set to '0', the signal is not smoothened. To turn on specify a value greater than '0'.
      For replication timing, a good value is to smoothen to 100k bases. If binSize is 100bp, this would be '1000'
  b:
    type: int?
    inputBinding:
      position: 3
      prefix: -b
    doc: |
      Which tracks should be smoothened when performing smoothening (default: '1' meaning only binned track).
      Options are:
        1) 0: Smoothen scaled and transformed tracks (log2, ratio or subtracted)
        2) 1: Smoothen only the scaled sequencing track
        3) 2: Smoothen only the transformed (log2, ratio or subtract) track
  q:
    type: int?
    inputBinding:
      position: 4
      prefix: -q
    doc: |
      Minimum (at least) mapping quality (default: 0)
  d:
    type: string?
    inputBinding:
      position: 4
      prefix: -d
    doc: |
      Keep duplicated reads (default: no)
  p:
    type: string?
    inputBinding:
      position: 4
      prefix: -p
    doc: |
      Do not filter un-proper alignments (default: filter)
  m:
    type: string?
    inputBinding:
      position: 4
      prefix: -m
    doc: |
      Do not remove reads with unmapped pairs
  g:
    type: int?
    inputBinding:
      position: 4
      prefix: -g
    doc: |
      Minimum fragment size for read pairs (default: 0)
  x:
    type: int?
    inputBinding:
      position: 4
      prefix: -x
    doc: |
      Maximum fragment size for read pairs (default: 2000)
  w:
    type: int?
    inputBinding:
      position: 4
      prefix: -w
    doc: |
      Filter reads based on fragment size (default: no)
  t:
    type: int?
    inputBinding:
      position: 5
      prefix: -t
    doc: |
      No. of threads to use (default: 1)
  bam:
    type:
      type: array
      items: File
      inputBinding:
        prefix: --bam
        separate: true
    secondaryFiles: .bai
    inputBinding:
      position: 6
    doc: |
      Input BAM file. This can be specified multiple times in case of multiple BAM files

outputs:
  output:
    type: File[]
    outputBinding:
       glob: "*.bw"

baseCommand: ["BAMscale", "scale"]

s:author:
  - class: s:Person
    s:identifier: https://orcid.org/0000-0002-4108-5982
    s:email: mailto:r78v10a07@gmail.com
    s:name: Roberto Vera Alvarez

s:codeRepository: https://github.com/ncbi/BAMscale
s:license: https://spdx.org/licenses/OPL-1.0

$namespaces:
  s: http://schema.org/

$schemas:
  - http://schema.org/docs/schema_org_rdfa.html

