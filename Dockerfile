# Use a lightweight miniconda base
FROM mambaorg/micromamba:1.5.3

# Install system libraries required by R spatial packages and hdf5
USER root
RUN apt-get update && apt-get install -y \
    libhdf5-dev libcurl4-openssl-dev libssl-dev libxml2-dev \
    fonts-liberation build-essential && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

USER $MAMBA_USER

# Create a single unified environment with both the R and Python spatial stacks
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    "r-base>=4.3.0" \
    "r-seurat>=5.0.0" \
    "r-hdf5r" \
    "r-remotes" \
    "r-tidyverse" \
    "r-optparse" \
    "r-patchwork" \
    "r-fs" \
    "r-glue" \
    "python>=3.10" \
    "scanpy>=1.9.0" \
    "squidpy>=1.3.0" \
    "matplotlib" \
    && micromamba clean --all --yes

# Install SeuratDisk directly from GitHub (it is not on CRAN or Conda)
RUN micromamba run -n base Rscript -e "remotes::install_github('mojaveazure/seurat-disk', upgrade = 'never')"

ENV PATH="/opt/conda/bin:$PATH"
CMD ["/bin/bash"]