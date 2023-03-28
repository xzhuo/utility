FROM conda/miniconda3

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -y htslib samtools bedtools pysam last

RUN mkdir /utility/
ADD ./ /utility