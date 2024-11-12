FROM nfcore/base
LABEL authors="Urmo Võsa" \
      description="Docker image containing tools for GWAS-eQTL LDSC analyses"

COPY environment.yml /
RUN apt-get update && apt install -y libgmp-dev && apt install -y build-essential
RUN conda install -n base -c conda-forge mamba && \
    mamba env create -f environment.yml && \
    conda clean -a
ENV PATH /opt/conda/envs/eqtlgenldsc/bin:$PATH
ENV TAR="/bin/tar"
RUN ln -s /bin/tar /bin/gtar
