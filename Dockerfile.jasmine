FROM condaforge/miniforge3

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && mamba install -y wget vim snakemake pbjasmine pbmm2 pbtk samtools bedtools hiphase

RUN wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.3.2/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz
RUN tar -xzf pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu.tar.gz
ENV PATH="$PATH:/pb-CpG-tools-v2.3.2-x86_64-unknown-linux-gnu/bin"