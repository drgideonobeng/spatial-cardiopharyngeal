#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * =============================================================================
 * Spatial Cardiopharyngeal Pipeline
 * =============================================================================
 */

// ============================ PROCESSES ====================================

process SEURAT_SPATIAL {
  tag "${meta.id}"
  publishDir "${params.results_dir}/${meta.id}/seurat", mode: 'copy'
  
  input:
    tuple val(meta), path(h5ad_file)
    
  output:
    tuple val(meta), path("${meta.id}_spatial.rds"), emit: rds
    path "*.pdf"                                   , emit: plots
    
  script:
  """
  01_spatial_seurat.R \\
    --in_h5ad        ${h5ad_file} \\
    --sample_id      ${meta.id} \\
    --timepoint      ${meta.timepoint} \\
    --replicate      ${meta.replicate} \\
    --puck_id        ${meta.puck_id} \\
    --genotype       ${meta.genotype} \\
    --markers        '${params.spatial_genes}' \\
    --out_rds        ${meta.id}_spatial.rds \\
    --out_pdf        ${meta.id}_featureplots.pdf
  """
}

process SQUIDPY_TOPOLOGY {
  tag "${meta.id}"
  publishDir "${params.results_dir}/${meta.id}/squidpy", mode: 'copy'
  
  input:
    tuple val(meta), path(h5ad_file)
    
  output:
    path "*.csv" , emit: metrics
    path "*.png" , emit: plots
    
  script:
  """
  02_topology_squidpy.py \\
    --in_h5ad   ${h5ad_file} \\
    --sample_id ${meta.id} \\
    --timepoint ${meta.timepoint} \\
    --replicate ${meta.replicate} \\
    --puck_id   ${meta.puck_id} \\
    --genotype  ${meta.genotype}
  """
}

// ============================ WORKFLOW =====================================
workflow {
  def cfg = new org.yaml.snakeyaml.Yaml().load(file(params.input_samples).text)

  samples_ch = Channel
    .fromList(cfg.samples)
    .map { s ->
        def meta = [
          id          : s.id,
          timepoint   : s.timepoint,
          replicate   : s.replicate,
          puck_id     : s.puck_id,
          genotype    : s.genotype,
          tissue      : s.tissue,
          assay       : s.assay,
          organism    : s.organism
        ]
        def h5ad_file = file("${params.raw_data_dir}/${s.id}.h5ad")
        tuple(meta, h5ad_file)
    }

  SEURAT_SPATIAL(samples_ch)
  SQUIDPY_TOPOLOGY(samples_ch)
}