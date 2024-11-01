FROM broadinstitute/gatk:4.6.1.0

RUN git clone https://github.com/lh3/seqtk.git && \
    cd seqtk && \
    make
