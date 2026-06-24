#!/usr/bin/env python3
# =============================================================================
# Step 02: SQUIDPY_TOPOLOGY
# =============================================================================

import argparse
import scanpy as sc
import squidpy as sq
import os

def main():
    parser = argparse.ArgumentParser(description="Squidpy Spatial Topology Module")
    parser.add_argument("--in_h5ad", required=True, help="Input .h5ad file")
    parser.add_argument("--sample_id", required=True, help="Sample ID")
    parser.add_argument("--timepoint", required=True, help="Timepoint")
    parser.add_argument("--replicate", required=True, help="Replicate run")
    parser.add_argument("--puck_id", required=True, help="Puck ID")
    parser.add_argument("--genotype", required=True, help="Genotype")
    args = parser.parse_args()

    print(f"[02] Loading AnnData object for {args.sample_id}...")
    adata = sc.read_h5ad(args.in_h5ad)
    
    # Inject Metadata into observation matrix (adata.obs)
    adata.obs["sample_id"] = args.sample_id
    adata.obs["timepoint"] = args.timepoint
    adata.obs["replicate"] = args.replicate
    adata.obs["puck_id"] = args.puck_id
    adata.obs["genotype"] = args.genotype
    adata.obs["pipeline_version"] = "v1.0-spatial-cardiopharyngeal"

    # Fix spatial coordinate mappings if necessary
    if 'spatial' not in adata.obsm:
        if 'X_spatial' in adata.obsm:
            adata.obsm['spatial'] = adata.obsm['X_spatial']
        else:
            raise ValueError(f"Spatial coordinates not found in {args.in_h5ad}")

    # Calculate Neighborhood Topology Graphs
    print(f"[02] Constructing spatial neighborhood graph...")
    sq.gr.spatial_neighbors(adata)

    if 'highly_variable' not in adata.var:
        sc.pp.highly_variable_genes(adata, n_top_genes=1000, flavor='seurat_v3')
    
    hvg_list = adata.var_names[adata.var['highly_variable']].tolist()[:1000]

    # Compute spatial structure parameters
    print(f"[02] Computing Moran's I for top spatially variable genes...")
    sq.gr.spatial_autocorr(adata, mode="moran", genes=hvg_list)
    
    csv_out = f"{args.sample_id}_morans_I.csv"
    adata.uns['moranI'].to_csv(csv_out)
    
    # Plot spatial maps
    top_genes = adata.uns['moranI'].head(3).index.tolist()
    fig_out = f"{args.sample_id}_top_moranI_genes.png"
    
    sq.pl.spatial_scatter(
        adata,
        color=top_genes,
        cmap="Reds",
        size=1.5,
        save=os.path.basename(fig_out)
    )
    
    default_path = f"figures/spatial_scatter{os.path.basename(fig_out)}"
    if os.path.exists(default_path):
        os.rename(default_path, fig_out)

    print(f"[02] Process completed. Wrote {csv_out} and {fig_out}")

if __name__ == "__main__":
    main()