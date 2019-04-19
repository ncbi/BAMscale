# Base Image
FROM continuumio/anaconda3

# Metadata
LABEL base.image="continuumio/anaconda3"
LABEL version="1"
LABEL software="BAMscale"
LABEL software.version="0.0.1"
LABEL description="BAMscale is a one-step tool for either 1) quantifying and normalizing the coverage of peaks or 2) generated scaled BigWig files for easy visualization of commonly used DNA-seq capture based methods."
LABEL tags="BAM"
LABEL website="https://github.com/ncbi/BAMscale"

# Maintainer
MAINTAINER Roberto Vera Alvarez <r78v10a07@gmail.com>

USER root

RUN apt-get update && \
    apt-get install -y apt-utils wget bzip2 sudo gcc make && \
    apt-get clean && \
    apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# Updating Anaconda packages
RUN /opt/conda/bin/conda update conda
RUN /opt/conda/bin/conda update anaconda
RUN /opt/conda/bin/conda update --all
RUN /opt/conda/bin/conda install -c bioconda htslib libbigwig

# Add user ubuntu with no password, add to sudo group
RUN adduser --disabled-password --gecos '' ubuntu
RUN chmod a+rwx /home/ubuntu/
RUN mkdir /home/ubuntu/bin
RUN chown -R ubuntu /home/ubuntu
USER ubuntu

ENV URL=https://github.com/ncbi/BAMscale
ENV FOLDER=BAMscale
ENV PATH="/home/ubuntu/bin:${PATH}"
ENV CONDA_DIR="/opt/conda/"
ENV CPPFLAGS="-I $CONDA_DIR/include"
ENV LDFLAGS="-L $CONDA_DIR/lib -Wl,-rpath,$CONDA_DIR/lib"

RUN cd /home/ubuntu/ && \
        git clone $URL && \
        cd $FOLDER && \
	make && \
	mv bin/BAMscale /home/ubuntu/bin/ && \
        cd .. && \
        rm -rf $FOLDER

WORKDIR /data

CMD ["/home/ubuntu/bin/BAMscale"]
