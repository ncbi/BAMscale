BAMscale
===

<img src="https://github.com/pongorlorinc/BAMscale/blob/master/doc/images/MAIN.png"  width="600" height="200" />

BAMscale is a one-step tool for either 1) quantifying and normalizing the coverage of peaks or 2) generated scaled BigWig files for easy visualization of commonly used DNA-seq capture based methods.

## Reference
\<NA\>

## Requirements

### samtools
http://www.htslib.org/

### libBigWig
Clone the libBigWig repository from GitHub: https://github.com/dpryan79/libBigWig

## Usage

### Peak quantification

Usage: ./bin/BAMscale cov --bed \<BED_FILE\> --bam \<BAM1\> --bam \<BAM2\> --bam \<BAM3\> ... --bam \<BAMn\>

### Generating scaled coverage tracks

Usage ./bin/BAMscale scale --bam \<BAM_FILE\> [--bam \<BAM2\> .. --bam \<BAMn\>]
