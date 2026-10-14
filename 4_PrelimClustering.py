import os
import pandas as pd
import numpy as np
import scanpy as sc
from sklearn.metrics import adjusted_rand_score
# ############################################################
# 1. INPUT AND OUTPUT PATHS
# Set the Harmony and scVI input objects and the clustering output directory.
# ############################################################
BASE_DIR = "/home/cv1u24/GIST/4_Preprocessing_and_Integration"
HARMONY_H5AD = os.path.join(
    BASE_DIR,
    "Outputs",
    "BatchToolComparisons",
    "Harmony",
    "Harmony_30PC_n15_md0p5",
    "h5ad",
    "merged_singlets_harmony_30PC_n15_md0p5.h5ad"
)
SCVI_H5AD = os.path.join(
    BASE_DIR,
    "Outputs",
    "BatchToolComparisons",
    "scVI",
    "scVI_30latent",
    "h5ad",
    "merged_singlets_scvi_integrated_30latent.h5ad"
)
OUT_BASE = os.path.join(
    BASE_DIR,
    "Outputs",
    "Clustering_Comparison"
)
os.makedirs(OUT_BASE, exist_ok=True)

# ############################################################
# 2. CLUSTERING PARAMETERS
# Define the neighbour counts, Leiden resolutions, and metadata keys to test.
# ############################################################
# These are the parameters we agreed to test
N_NEIGHBORS_OPTIONS = [15, 25, 50]
LEIDEN_RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0]
N_PCS = 30
BATCH_KEY = "sample"
GENOTYPE_KEY = "genotype_group"
sc.settings.verbosity = 1

# ############################################################
# 3. DEFINE CLUSTERING WORKFLOW
# Run clustering grids for a supplied representation and save figures and summaries.
# ############################################################
def run_clustering(
    adata_base,
    method,
    representation_key,
    n_neighbors_options,
    resolutions,
    output_base
):
    print("\n")
    print("============================================================")
    print(f"RUNNING LEIDEN CLUSTERING: {method}")
    print("============================================================")
    if representation_key not in adata_base.obsm:
        raise ValueError(
            f"Representation '{representation_key}' not found in adata.obsm.\n"
            f"Available representations:\n{list(adata_base.obsm.keys())}"
        )
    print(f"Using representation: {representation_key}")
    print(f"Cells: {adata_base.n_obs}")
    print(f"Genes: {adata_base.n_vars}")
    method_out = os.path.join(output_base, method)
    h5ad_out = os.path.join(method_out, "h5ad")
    fig_out = os.path.join(method_out, "Figures")
    csv_out = os.path.join(method_out, "CSVs")
    os.makedirs(h5ad_out, exist_ok=True)
    os.makedirs(fig_out, exist_ok=True)
    os.makedirs(csv_out, exist_ok=True)
    summary_rows = []
    composition_rows = []
    percentage_rows = []
    # --------------------------------------------------------
    # Loop through neighbour numbers
    # --------------------------------------------------------
    for n_neighbors in n_neighbors_options:
        print("\n------------------------------------------------------------")
        print(f"{method}: n_neighbors = {n_neighbors}")
        print("------------------------------------------------------------")
        # Fresh copy for each neighbour setting
        adata = adata_base.copy()
        neighbors_key = f"neighbors_{method}_n{n_neighbors}"
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            use_rep=representation_key,
            key_added=neighbors_key
        )
        # ----------------------------------------------------
        # Generate UMAP
        # ----------------------------------------------------
        umap_key = f"X_umap_{method}_n{n_neighbors}"
        sc.tl.umap(
            adata,
            neighbors_key=neighbors_key,
            min_dist=0.3
        )
        adata.obsm[umap_key] = adata.obsm["X_umap"].copy()
        # ----------------------------------------------------
        # Test Leiden resolutions
        # ----------------------------------------------------
        for resolution in resolutions:
            resolution_label = str(resolution).replace(".", "p")
            cluster_key = (
                f"leiden_{method}_n{n_neighbors}_r{resolution_label}"
            )
            print(
                f"Running Leiden: "
                f"method={method}, "
                f"n_neighbors={n_neighbors}, "
                f"resolution={resolution}"
            )
            sc.tl.leiden(
                adata,
                resolution=resolution,
                neighbors_key=neighbors_key,
                key_added=cluster_key,
                directed=False,
                flavor="igraph",
                n_iterations=2
            )
            # ------------------------------------------------
            # Cluster counts
            # ------------------------------------------------
            n_clusters = adata.obs[cluster_key].nunique()
            cluster_counts = (
                adata.obs[cluster_key]
                .value_counts()
                .sort_index()
            )
            for cluster, count in cluster_counts.items():
                summary_rows.append({
                    "method": method,
                    "n_neighbors": n_neighbors,
                    "resolution": resolution,
                    "cluster": cluster,
                    "cluster_size": count,
                    "total_cells": adata.n_obs,
                    "n_clusters": n_clusters
                })
            # ------------------------------------------------
            # Sample composition
            # ------------------------------------------------
            if BATCH_KEY in adata.obs.columns:
                sample_table = pd.crosstab(
                    adata.obs[cluster_key],
                    adata.obs[BATCH_KEY]
                )
                for cluster in sample_table.index:
                    row = {
                        "method": method,
                        "n_neighbors": n_neighbors,
                        "resolution": resolution,
                        "cluster": cluster
                    }
                    for sample in sample_table.columns:
                        row[sample] = sample_table.loc[
                           cluster, sample
                        ]
                    composition_rows.append(row)
                # --------------------------------------------
                # Sample percentages
                # --------------------------------------------
                sample_percentage = sample_table.div(
                    sample_table.sum(axis=1),
                    axis=0
                ) * 100
                for cluster in sample_percentage.index:
                    row = {
                        "method": method,
                        "n_neighbors": n_neighbors,
                        "resolution": resolution,
                        "cluster": cluster
                    }
                    for sample in sample_percentage.columns:
                        row[sample] = sample_percentage.loc[
                           cluster, sample
                        ]
                    percentage_rows.append(row)
            # ------------------------------------------------
            # Save cluster UMAP by cluster
            # ------------------------------------------------
            sc.pl.embedding(
                adata,
                basis=umap_key.replace("X_", ""),
                color=cluster_key,
                show=False,
                title=(
                    f"{method} Leiden | "
                    f"30 dimensions | "
                    f"n={n_neighbors} | "
                    f"resolution={resolution}"
                ),
                legend_loc="right margin"
            )
            fig_file = os.path.join(
                fig_out,
                f"{method}_leiden_n{n_neighbors}_"
                f"r{resolution_label}.png"
            )
            import matplotlib.pyplot as plt
            plt.savefig(
                fig_file,
                dpi=300,
                bbox_inches="tight"
            )
            plt.close()
            # ------------------------------------------------
            # Save UMAP coloured by sample
            # ------------------------------------------------
            if BATCH_KEY in adata.obs.columns:
                sc.pl.embedding(
                    adata,
                    basis=umap_key.replace("X_", ""),
                    color=BATCH_KEY,
                    show=False,
                    title=(
                        f"{method} sample | "
                        f"30 dimensions | "
                        f"n={n_neighbors}"
                    )
                )
                fig_file = os.path.join(
                    fig_out,
                    f"{method}_sample_n{n_neighbors}.png"
                )
                plt.savefig(
                    fig_file,
                    dpi=300,
                    bbox_inches="tight"
                )
                plt.close()
            # ------------------------------------------------
            # Save UMAP coloured by genotype
            # ------------------------------------------------
            if GENOTYPE_KEY in adata.obs.columns:
                sc.pl.embedding(
                    adata,
                    basis=umap_key.replace("X_", ""),
                    color=GENOTYPE_KEY,
                    show=False,
                    title=(
                        f"{method} genotype | "
                        f"30 dimensions | "
                        f"n={n_neighbors}"
                    )
                )
                fig_file = os.path.join(
                    fig_out,
                    f"{method}_genotype_n{n_neighbors}.png"
                )
                plt.savefig(
                    fig_file,
                    dpi=300,
                    bbox_inches="tight"
                )
                plt.close()
        # ----------------------------------------------------
        # Save object for this neighbour number
        # ----------------------------------------------------
        output_h5ad = os.path.join(
            h5ad_out,
            f"{method}_30dim_n{n_neighbors}_leiden.h5ad"
        )
        adata.write_h5ad(output_h5ad)
        print(f"Saved: {output_h5ad}")
    # ############################################################
    # 3A. SAVE CLUSTERING SUMMARIES
    # Write cluster counts and sample composition tables for the current representation.
    # ############################################################
    summary_df = pd.DataFrame(summary_rows)
    summary_file = os.path.join(
        csv_out,
        f"{method}_cluster_counts.csv"
    )
    summary_df.to_csv(
        summary_file,
        index=False
    )
    composition_df = pd.DataFrame(composition_rows)
    composition_file = os.path.join(
        csv_out,
        f"{method}_cluster_sample_composition.csv"
    )
    composition_df.to_csv(
        composition_file,
        index=False
    )
    percentage_df = pd.DataFrame(percentage_rows)
    percentage_file = os.path.join(
        csv_out,
        f"{method}_cluster_sample_percentage.csv"
    )
    percentage_df.to_csv(
        percentage_file,
        index=False
    )
    print("\nSaved CSVs:")
    print(summary_file)
    print(composition_file)
    print(percentage_file)
    return summary_df

# ############################################################
# 4. LOAD AND CLUSTER HARMONY
# Apply the clustering workflow to the selected Harmony representation.
# ############################################################
print("\n============================================================")
print("LOADING HARMONY OBJECT")
print("============================================================")
if not os.path.exists(HARMONY_H5AD):
    raise FileNotFoundError(
        f"Harmony file not found:\n{HARMONY_H5AD}"
    )
adata_harmony = sc.read_h5ad(HARMONY_H5AD)
print(adata_harmony)
harmony_summary = run_clustering(
    adata_base=adata_harmony,
    method="Harmony",
    representation_key="X_pca_harmony_30PC",
    n_neighbors_options=N_NEIGHBORS_OPTIONS,
    resolutions=LEIDEN_RESOLUTIONS,
    output_base=OUT_BASE
)

# ############################################################
# 5. LOAD AND CLUSTER scVI
# Apply the clustering workflow to the selected scVI latent representation.
# ############################################################
print("\n============================================================")
print("LOADING scVI OBJECT")
print("============================================================")
if not os.path.exists(SCVI_H5AD):
    raise FileNotFoundError(
        f"scVI file not found:\n{SCVI_H5AD}"
    )
adata_scvi = sc.read_h5ad(SCVI_H5AD)
print(adata_scvi)
scvi_summary = run_clustering(
    adata_base=adata_scvi,
    method="scVI",
    representation_key="X_scVI_30latent",
    n_neighbors_options=N_NEIGHBORS_OPTIONS,
    resolutions=LEIDEN_RESOLUTIONS,
    output_base=OUT_BASE
)

# ############################################################
# 6. COMBINE SUMMARY RESULTS
# Combine the Harmony and scVI cluster-count summaries into one CSV.
# ############################################################
combined_summary = pd.concat(
    [harmony_summary, scvi_summary],
    ignore_index=True
)
combined_file = os.path.join(
    OUT_BASE,
    "combined_cluster_counts.csv"
)
combined_summary.to_csv(
    combined_file,
    index=False
)

# ############################################################
# 7. PRINT COMPLETION SUMMARY
# Report the settings tested and output location.
# ############################################################
print("\n")
print("============================================================")
print("LEIDEN CLUSTERING COMPLETE")
print("============================================================")
print(f"\nOutput directory:")
print(OUT_BASE)
print("\nMethods tested:")
print("  Harmony - 30 PCs")
print("  scVI    - 30 latent dimensions")
print("\nNeighbour numbers:")
print(N_NEIGHBORS_OPTIONS)
print("\nLeiden resolutions:")
print(LEIDEN_RESOLUTIONS)
print("\nMain combined output:")
print(combined_file)
print("\nDone.")
