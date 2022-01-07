FROM docker pull quay.io/hdc-workflows/ubuntu

LABEL maintainer "Jared Galloway <jgallowa@fredhutch.rg>" \
      version "1.0" \
      description "Common PhIP-Seq Workflows"

# install needed tools
RUN apt-get update --fix-missing -qq && apt-get install -y -q \
    git \
    curl \
    locales \
    libncurses5-dev  \
    libncursesw5-dev \
    build-essential \
    pkg-config \
    zlib1g-dev \
    python3 \
    python3-pip \ 
    python3-venv \
    bowtie

ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# install SAMtools
RUN curl -fksSL https://github.com/samtools/samtools/releases/download/1.3.1/samtools-1.3.1.tar.bz2 | tar xj && \
    cd samtools-1.3.1 && \
    make all all-htslib && make install install-htslib

# install phippery
RUN git clone https://github.com/matsengrp/phippery.git \
    && (cd phippery && pip install .)
#COPY ../../phippery/ .
#RUN (cd phippery && pip install .)
