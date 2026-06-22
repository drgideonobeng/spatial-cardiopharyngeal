# Spatial Transcriptomics: Pharyngeal Pipeline (E8-E11)

A modular, reproducible Nextflow (DSL2) pipeline designed to process spatial transcriptomics datasets with a specific focus on the developing mouse pharyngeal and cardiopharyngeal regions (E8-E11).

## Architecture
This pipeline combines two distinct analytical environments into a unified workflow using Docker:
1. **R (Seurat v5)**: Converts `.h5ad` files, computes Spatially Variable Features (SVFs), and maps target cardiopharyngeal markers (e.g., Ttn, Foxf1).
2. **Python (Squidpy/Scanpy)**: Constructs spatial neighborhood graphs and calculates spatial autocorrelation (Moran's I metrics) to detect morphological patterning.

## Dataset Target
* **Source:** Slide-seq Mouse Organogenesis Atlas (Chen et al., 2023) via CZ CELLxGENE.
* **Stages:** E8.5, E9.0, E9.5.

## Requirements
* Nextflow (>=24.04.0)
* Docker
