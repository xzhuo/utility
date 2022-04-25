FROM conda/miniconda3

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -y htslib samtools bedtools pysam

RUN mkdir /utility/
ADD ./ /utility