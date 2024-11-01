FROM broadinstitute/gatk:4.6.1.0

ENV DEBIAN_FRONTEND=noninteractive

USER root

## Install essential
RUN apt-get update && apt-get install -y --no-install-recommends build-essential python3 python3-pip python3-setuptools python3-dev wget bedtools bcftools imagemagick less ghostscript gawk nano openjdk-8-jdk vcftools liblzma-dev libbz2-dev zlib1g-dev libncurses5-dev libncursesw5-dev libcurl4-openssl-dev perl

RUN apt-get install -y libz-dev gcc zlib-devel

RUN git clone https://github.com/lh3/seqtk.git && \
    cd seqtk && \
    make

WORKDIR /workspace