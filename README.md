BAMscale
===

**Overview of BAMscale applications**


<p align="center">
<img src="https://github.com/ncbi/BAMscale/blob/master/doc/images/MAIN_figure.png"  width="800" height="220" />
</p>


BAMscale is a one-step tool to 

    1) quantify/normalize peak coverages from multiple BAM files 
    2) Create scaled BigWig files for easy visualization

In the [wiki](https://github.com/ncbi/BAMscale/wiki) pages we have more detailed tutorials for creating bigWig files and quantifying peaks

## Update

20190821: We recently added support for [RNA-seq](https://github.com/ncbi/BAMscale/wiki/Detailed-usage:-RNA-seq-coverage-tracks) data as well to create coverage tracks. The new method enables accurate representations of exon-intron boundaries (splicing). 

## Manuals

In the [wiki](https://github.com/ncbi/BAMscale/wiki) page we have more detailed tutorials for creating bigWig files and quantifying peaks:

1. [OK-seq and RFD Track Generation](https://github.com/ncbi/BAMscale/wiki/Detailed-Use:-OKseq-RFD-(Replication-Fork-Directionality)-Track-Generation)
2. [Quantifying Peaks](https://github.com/ncbi/BAMscale/wiki/Detailed-Use:-Quantifying-Peak-Coverages-from-Multiple-BAM-Files#comparing-atac-seq-changes-induced-from-treatment)
3. [Generating Scaled Coverage Tracks](https://github.com/ncbi/BAMscale/wiki/Detailed-Use:-Generating-Scaled-Coverage-Tracks#preparing-input-data-for-bamscale)
4. [END-seq data](https://github.com/ncbi/BAMscale/wiki/Detailed-Use:-Processing-END-seq-Data)
5. [Log2 Coverage Tracks for Replication Timing Data](https://github.com/ncbi/BAMscale/wiki/Detailed-Use:-Replication-Timing-log2-Coverage-Ratio-from-Two-BAM-Files)
6. [Smoothening Function for Coverage Tracks](https://github.com/ncbi/BAMscale/wiki/Detailed-Use:-Smooth-Coverage-Tracks)


We also added a few R scripts that might be helpful for basic visualizations:

1. [Segmenting replication timing bigwigs](https://github.com/ncbi/BAMscale/wiki/Replication-timing-BED-segments-from-bigwig)
2. [Identifying OK-seq strand switches](https://github.com/ncbi/BAMscale/wiki/Finding-OK-seq-strand-switched-from-the-RFD-track)

For additional information, visit the [wiki](https://github.com/ncbi/BAMscale/wiki) page.

For any other requests, or if you need help either open an issue, or feel free to email me: *pongorlorinc@gmail.com*


## Usage for the impatient

These examples assume you have 4 processing threads, so we set '-t 4' for multithreading.

#### Peak quantification

    BAMscale cov -t 4 --bed <BED_FILE> --bam <BAM1> --bam <BAM2> --bam <BAM3> ... --bam <BAMn>

#### Generating scaled coverage tracks

***Creating scaled coverage tracks***

    BAMscale scale -t 4 --bam <BAM_FILE> [--bam <BAM2> .. --bam <BAMn>]

***Creating stranded RNA-seq coverage tracks***

    BAMscale scale --operation strandrna --bam <RNAseq.bam>
    
***Creating unstranded coverage from RNA-seq***

    BAMscale scale --operation rna --bam <RNAseq.bam>

***Getting RFD score from OKseq data***

    BAMscale scale -t 4 --operation rfd --binsize 1000 --bam <BAM_FILE>
    
***Processing replication timing and Repli-seq data***

    BAMscale scale -t 4 --operation reptime --bam <G1_phase.bam> --bam <S_phase.bam>
    
***Creating stranded END-seq coverages***

    BAMscale scale -t 4 --operation endseq --bam <ENDseq.bam>


## Reference

BAMscale can be found at **bioR&chi;iv** ([https://doi.org/10.1101/669275](https://www.biorxiv.org/content/10.1101/669275v1))

## Bioconda instalation

[BAMscale](https://bioconda.github.io/recipes/bamscale/README.html) is available through [Bioconda](https://bioconda.github.io/). Read the Bioconda [Getting Started](https://bioconda.github.io/user/install.html#install-conda) page for a detailed description on how to get Bioconda installed.

Once Bioconda is available you can install BAMscale using this command.

    conda install bamscale

## Docker

BAMscale docker image is available in [quay.io/biocontainers/bamscale](https://quay.io/repository/biocontainers/bamscale).

### Pulling the image

    docker pull quay.io/biocontainers/bamscale:0.0.5--ha85820d_0
    
### Using the Docker image

#### Peak quantification with Docker

    docker run -v `pwd`:/data bamscale BAMscale cov --bed <BED_FILE> --bam <BAM1> --bam <BAM2> --bam <BAM3> ... --bam <BAMn>

#### Generating scaled coverage tracks with Docker

    docker run -v `pwd`:/data bamscale BAMscale scale --bam <BAM_FILE> [--bam <BAM2> .. --bam <BAMn>]

### Creating a custom docker image

    docker build -t bamscale https://raw.githubusercontent.com/pongorlorinc/BAMscale/master/Dockerfile

## Local compilation

### Requirements

We have a detailed installation for [Linux](https://github.com/ncbi/BAMscale/wiki/Installation#detailed-installation-for-linux-based-os) and [MAC](https://github.com/ncbi/BAMscale/wiki/Installation#detailed-installation-for-mac-os-with-homebrew) (with homebrew) based systems or through [conda](https://github.com/ncbi/BAMscale/wiki/Installation#detailed-installation-for-mac-os-with-conda). There is also a precompiled version for linux ready for usage available at the [releases](https://github.com/ncbi/BAMscale/releases).

#### samtools
http://www.htslib.org/

#### libBigWig
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

### Installation

After compiling the libBigWig library and samtools (if not already installed) clone the BAMscale from GitHub

    git clone https://github.com/ncbi/BAMscale.git
    
and go to the BAMscale folder to compile the program:

    cd BAMscale/
    make
    
A bin folder will be created with the BAMscale executable.

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
