BAMscale
===

<img src="https://github.com/pongorlorinc/BAMscale/blob/master/doc/images/MAIN.png"  width="600" height="200" />

BAMscale is a one-step tool for either 1) quantifying and normalizing the coverage of peaks or 2) generated scaled BigWig files for easy visualization of commonly used DNA-seq capture based methods.

For detailed information, visit the [wiki](https://github.com/pongorlorinc/BAMscale/wiki) page

## Reference



## Requirements

### samtools
http://www.htslib.org/

### libBigWig
Clone the libBigWig repository from GitHub: https://github.com/dpryan79/libBigWig

## Usage

### Peak quantification

    BAMscale cov --bed <BED_FILE> --bam <BAM1> --bam <BAM2> --bam <BAM3> ... --bam <BAMn>

### Generating scaled coverage tracks

    BAMscale scale --bam <BAM_FILE> [--bam <BAM2> .. --bam <BAMn>]


## Docker

### Build docker image

    docker build -t bamscale https://raw.githubusercontent.com/pongorlorinc/BAMscale/master/Dockerfile

### Peak quantification with Docker

    docker run -v `pwd`:/data bamscale BAMscale cov --bed <BED_FILE> --bam <BAM1> --bam <BAM2> --bam <BAM3> ... --bam <BAMn>

### Generating scaled coverage tracks with Docker

    docker run -v `pwd`:/data bamscale BAMscale scale --bam <BAM_FILE> [--bam <BAM2> .. --bam <BAMn>]

# Public Domain notice

National Center for Biotechnology Information.

This software is a "United States Government Work" under the terms of the United States
Copyright Act. It was written as part of the authors' official duties as United States
Government employees and thus cannot be copyrighted. This software is freely available
to the public for use. The National Library of Medicine and the U.S. Government have not
 placed any restriction on its use or reproduction.

Although all reasonable efforts have been taken to ensure the accuracy and reliability
of the software and data, the NLM and the U.S. Government do not and cannot warrant the
performance or results that may be obtained by using this software or data. The NLM and
the U.S. Government disclaim all warranties, express or implied, including warranties
of performance, merchantability or fitness for any particular purpose.

Please cite NCBI in any work or product based on this material.
