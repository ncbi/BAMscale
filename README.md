BAMscale
===

<p align="center">
<img src="https://github.com/ncbi/BAMscale/blob/master/doc/images/MAIN.png"  width="800" height="250" />
</p>

BAMscale is a one-step tool for either 1) quantifying and normalizing the coverage of peaks or 2) generated scaled BigWig files for easy visualization of commonly used DNA-seq capture based methods.


In the [wiki](https://github.com/ncbi/BAMscale/wiki) page we have more detailed tutorials for:

1. [OK-seq and RFD Track Generation](https://github.com/ncbi/BAMscale/wiki/Detailed-Use:-OKseq-RFD-(Replication-Fork-Directionality)-Track-Generation)
2. [Quantifying Peaks](https://github.com/ncbi/BAMscale/wiki/Detailed-Use:-Quantifying-Peak-Coverages-from-Multiple-BAM-Files#comparing-atac-seq-changes-induced-from-treatment)
3. [Generating Scaled Coverage Tracks](https://github.com/ncbi/BAMscale/wiki/Detailed-Use:-Generating-Scaled-Coverage-Tracks#preparing-input-data-for-bamscale)
4. [END-seq data](https://github.com/ncbi/BAMscale/wiki/Detailed-Use:-Processing-END-seq-Data)
5. [Log2 Coverage Tracks for Replication Timing Data](https://github.com/ncbi/BAMscale/wiki/Detailed-Use:-Replication-Timing-log2-Coverage-Ratio-from-Two-BAM-Files)
6. [Smoothening Function for Coverage Tracks](https://github.com/ncbi/BAMscale/wiki/Detailed-Use:-Smooth-Coverage-Tracks)


For additional information, visit the [wiki](https://github.com/ncbi/BAMscale/wiki) page.

## Reference

BAMscale can be found at **bioR&chi;iv** ([https://doi.org/10.1101/669275](https://www.biorxiv.org/content/10.1101/669275v1))


## Requirements

We have a detailed installation for [Linux](https://github.com/ncbi/BAMscale/wiki/Installation#detailed-installation-for-linux-based-os) and [MAC](https://github.com/ncbi/BAMscale/wiki/Installation#detailed-installation-for-mac-os-with-homebrew) (with homebrew) based systems or through [conda](https://github.com/ncbi/BAMscale/wiki/Installation#detailed-installation-for-mac-os-with-conda). There is also a precompiled version for linux ready for usage available at the [releases](https://github.com/ncbi/BAMscale/releases).

### samtools
http://www.htslib.org/

### libBigWig
Clone the libBigWig repository from GitHub: https://github.com/dpryan79/libBigWig

    git clone https://github.com/dpryan79/libBigWig.git

Compile it and set the environment variables for BAMscale

    cd libBigWig/
    make
    export LIBBIGWIG_DIR=`pwd`
    export CPPFLAGS="-I $LIBBIGWIG_DIR"
    export LDFLAGS="-L $LIBBIGWIG_DIR -Wl,-rpath,$LIBBIGWIG_DIR"
    
Optionally (and if you have permission), the libbigwig can also be installed

    make install
    
In this case, the flags don't have to be set in the terminal.

## Installation

After compiling the libBigWig library and samtools (if not already installed) clone the BAMscale from GitHub

    git clone https://github.com/ncbi/BAMscale.git
    
and go to the BAMscale folder to compile the program:

    cd BAMscale/
    make
    
A bin folder will be created with the BAMscale executable.

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
