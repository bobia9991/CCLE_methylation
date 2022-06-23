FROM rocker/r-ubuntu:18.04

ENV DEBIAN_FRONTEND="noninteractive" TZ="America/New_York"

RUN apt-get update && apt-get -y upgrade && \
    apt-get install -y build-essential curl git language-pack-en libbz2-dev libcurl3-dev liblzma-dev libncurses5-dev libtbb-dev libxml2-dev software-properties-common wget zip zlib1g-dev && \
    apt-get clean && apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* 

RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('methylKit')"

#Set language
ENV LANGUAGE en_US.UTF-8
ENV LC_ALL en_US.UTF-8

WORKDIR /usr/src

# Bowtie 2
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.4.2/bowtie2-2.4.2-linux-x86_64.zip && \
    unzip bowtie2-2.4.2-linux-x86_64.zip && \
    rm bowtie2-2.4.2-linux-x86_64.zip

# Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2 && \
    tar jxf samtools-1.12.tar.bz2 && \
    rm samtools-1.12.tar.bz2 && \
    cd samtools-1.12 && \
    ./configure --prefix $(pwd) && \
    make

# Bismark
RUN wget https://github.com/FelixKrueger/Bismark/archive/refs/tags/0.23.0.tar.gz && \
    tar xvfz 0.23.0.tar.gz && \
    rm 0.23.0.tar.gz

COPY run_methylkit.R /usr/src/

ENV PATH="${PATH}:/usr/src/bowtie2-2.4.2-linux-x86_64:/usr/src/Bismark-0.23.0:/usr/src/samtools-1.12"
