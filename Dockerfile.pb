FROM snakemake/snakemake

RUN conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install -y wget htslib samtools bedtools minimap2 pbccs pbmm2 pbbam pbtk primrose

RUN wget https://github.com/PacificBiosciences/pb-CpG-tools/releases/download/v2.2.0/pb-CpG-tools-v2.2.0-x86_64-unknown-linux-gnu.tar.gz
RUN tar -xzf pb-CpG-tools-v2.2.0-x86_64-unknown-linux-gnu.tar.gz
ENV PATH="$PATH:/tmp/repo/pb-CpG-tools-v2.2.0-x86_64-unknown-linux-gnu/bin"