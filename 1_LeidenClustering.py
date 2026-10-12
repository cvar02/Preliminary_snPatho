#!/usr/bin/env python
# ############################################################
# 1. SCRIPT SETUP
# Load the packages used by this analysis.
# ############################################################
import os
import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# ############################################################
# 2. PATHS
# Define the project root and the two integrated scVI input files.
# ############################################################
BASE_DIR = "/home/cv1u24/GIST/5_Clustering"
INPUTS = {
   25: (
       "/iridisfs/home/cv1u24/GIST/4_Preprocessing_and_Integration/"
       "Outputs/BatchToolComparisons/scVI/h5ad/"
       "merged_scVI_25PCs.h5ad"
   ),
   30: (
       "/iridisfs/home/cv1u24/GIST/4_Preprocessing_and_Integration/"
       "Outputs/BatchToolComparisons/scVI/h5ad/"
       "merged_scVI_30PCs.h5ad"
   )
}

# ############################################################
# 3. OUTPUT DIRECTORIES
# Create folders for AnnData objects, figures and tables.
# ############################################################
OUT_DIR = os.path.join(
   BASE_DIR,
   "Outputs",
   "Leiden"
)
H5AD_OUT = os.path.join(
   OUT_DIR,
   "h5ad"
)
FIG_OUT = os.path.join(
   OUT_DIR,
   "Figures"
)
CSV_OUT = os.path.join(
   OUT_DIR,
   "CSVs"
)
os.makedirs(H5AD_OUT, exist_ok=True)
os.makedirs(FIG_OUT, exist_ok=True)
os.makedirs(CSV_OUT, exist_ok=True)

# ############################################################
# 4. PARAMETERS
# Set the neighbour counts, Leiden resolutions and metadata keys to screen.
# ############################################################
N_NEIGHBORS_LIST = [
   10,
   15,
   20,
   30
]
RESOLUTIONS = [
   0.2,
   0.4,
   0.6,
   0.8,
   1.0
]
RANDOM_STATE = 0
LATENT_KEY = "X_scVI"
SAMPLE_KEY = "sample"
SUBTYPE_KEY = "gist_subtype"

# ############################################################
# 5. STORAGE FOR MASTER RESULTS
# Initialise the table that will collect one row per parameter combination.
# ############################################################
parameter_results = []

# ============================================================
# PROCESS scVI 25 AND scVI 30
# ============================================================
for N_LATENT, INPUT_H5AD in INPUTS.items():
   print("\n============================================================")
   print(f"PROCESSING scVI {N_LATENT}")
   print("============================================================")
   # --------------------------------------------------------
   # Check input exists
   # --------------------------------------------------------
   if not os.path.exists(INPUT_H5AD):
       raise FileNotFoundError(
           f"Input file not found:\n{INPUT_H5AD}"
       )
   # --------------------------------------------------------
   # Load AnnData
   # --------------------------------------------------------
   print("\nLoading:")
   print(INPUT_H5AD)
   adata = sc.read_h5ad(INPUT_H5AD)
   print("\nLoaded:")
   print(adata)

   # ========================================================
   # CHECK REQUIRED DATA
   # ========================================================
   if LATENT_KEY not in adata.obsm:
       raise ValueError(
           f"{LATENT_KEY} not found in adata.obsm."
       )
   actual_latent_dimensions = (
       adata.obsm[LATENT_KEY].shape[1]
   )
   if actual_latent_dimensions != N_LATENT:
       raise ValueError(
           f"Expected {N_LATENT} scVI dimensions, "
           f"but {INPUT_H5AD} contains "
           f"{actual_latent_dimensions} dimensions."
       )
   if SAMPLE_KEY not in adata.obs.columns:
       raise ValueError(
           f"{SAMPLE_KEY} not found in adata.obs."
       )
   if SUBTYPE_KEY not in adata.obs.columns:
       print(
           f"WARNING: {SUBTYPE_KEY} not found. "
           "Subtype composition will be skipped."
       )

   print(
       f"\nConfirmed scVI latent space: "
       f"{LATENT_KEY} = "
       f"{adata.obsm[LATENT_KEY].shape}"
   )
   print(
       f"Samples: "
       f"{adata.obs[SAMPLE_KEY].unique().tolist()}"
   )
   if SUBTYPE_KEY in adata.obs.columns:
       print(
           f"GIST subtypes: "
           f"{adata.obs[SUBTYPE_KEY].unique().tolist()}"
       )

   # ========================================================
   # LOOP THROUGH NEIGHBOUR SETTINGS
   # ========================================================
   for N_NEIGHBORS in N_NEIGHBORS_LIST:
       print("\n------------------------------------------------------------")
       print(
           f"scVI {N_LATENT} | "
           f"n_neighbors = {N_NEIGHBORS}"
       )
       print("------------------------------------------------------------")

       # ----------------------------------------------------
       # Compute neighbour graph
       # ----------------------------------------------------
       neighbors_key = (
           f"neighbors_scVI{N_LATENT}"
           f"_n{N_NEIGHBORS}"
       )
       sc.pp.neighbors(
           adata,
           use_rep=LATENT_KEY,
           n_neighbors=N_NEIGHBORS,
           key_added=neighbors_key,
           random_state=RANDOM_STATE
       )

       # ----------------------------------------------------
       # Compute UMAP
       # ----------------------------------------------------
       print("Computing UMAP...")
       sc.tl.umap(
           adata,
           neighbors_key=neighbors_key,
           random_state=RANDOM_STATE
       )
       umap_key = (
           f"X_umap_scVI{N_LATENT}"
           f"_n{N_NEIGHBORS}"
       )
       adata.obsm[umap_key] = (
           adata.obsm["X_umap"].copy()
       )

       # ====================================================
       # LOOP THROUGH LEIDEN RESOLUTIONS
       # ====================================================
       for resolution in RESOLUTIONS:
           res_label = str(
               resolution
           ).replace(".", "_")
           cluster_key = (
               f"leiden_scVI{N_LATENT}"
               f"_n{N_NEIGHBORS}"
               f"_r{res_label}"
           )
           print("\n--------------------------------------------")
           print(
               f"Leiden | "
               f"scVI {N_LATENT} | "
               f"n={N_NEIGHBORS} | "
               f"resolution={resolution}"
           )
           print("--------------------------------------------")

           # ------------------------------------------------
           # Leiden clustering
           # ------------------------------------------------
           sc.tl.leiden(
               adata,
               resolution=resolution,
               key_added=cluster_key,
               neighbors_key=neighbors_key,
               random_state=RANDOM_STATE
           )
           n_clusters = (
               adata.obs[cluster_key]
               .nunique()
           )
           print(
               f"Clusters identified: "
               f"{n_clusters}"
           )

           # =================================================
           # CLUSTER SIZE INFORMATION
           # =================================================
           cluster_counts = (
               adata.obs[cluster_key]
               .value_counts()
               .sort_index()
           )
           cluster_counts_df = (
               cluster_counts
               .reset_index()
           )
           cluster_counts_df.columns = [
               "cluster",
               "cell_count"
           ]
           min_cluster_size = (
               int(cluster_counts.min())
           )
           max_cluster_size = (
               int(cluster_counts.max())
           )
           median_cluster_size = (
               float(cluster_counts.median())
           )
           clusters_under_50 = (
               int((cluster_counts < 50).sum())
           )
           clusters_under_100 = (
               int((cluster_counts < 100).sum())
           )
           clusters_over_500 = (
               int((cluster_counts > 500).sum())
           )

           # ------------------------------------------------
           # Save cluster counts
           # ------------------------------------------------
           counts_csv = os.path.join(
               CSV_OUT,
               f"cluster_counts_scVI{N_LATENT}"
               f"_n{N_NEIGHBORS}"
               f"_resolution_{res_label}.csv"
           )
           cluster_counts_df.to_csv(
               counts_csv,
               index=False
           )

           # =================================================
           # SAMPLE COMPOSITION
           # =================================================
           sample_comp = pd.crosstab(
               adata.obs[cluster_key],
               adata.obs[SAMPLE_KEY]
           )
           sample_csv = os.path.join(
               CSV_OUT,
               f"sample_composition_scVI{N_LATENT}"
               f"_n{N_NEIGHBORS}"
               f"_resolution_{res_label}.csv"
           )
           sample_comp.to_csv(
               sample_csv
           )

           # ------------------------------------------------
           # Sample dominance
           # ------------------------------------------------
           sample_proportions = (
               sample_comp.div(
                   sample_comp.sum(axis=1),
                   axis=0
               )
           )
           dominant_sample_fraction = (
               sample_proportions.max(axis=1)
           )
           mean_dominant_sample_fraction = (
               float(
                   dominant_sample_fraction.mean()
               )
           )
           clusters_over_75pct_single_sample = (
               int(
                   (
                       dominant_sample_fraction >= 0.75
                   ).sum()
               )
           )
           clusters_over_90pct_single_sample = (
               int(
                   (
                       dominant_sample_fraction >= 0.90
                   ).sum()
               )
           )

           # =================================================
           # SUBTYPE COMPOSITION
           # =================================================
           subtype_csv = None
           subtype_dominance = np.nan
           clusters_over_75pct_single_subtype = np.nan
           clusters_over_90pct_single_subtype = np.nan
           if SUBTYPE_KEY in adata.obs.columns:
               subtype_comp = pd.crosstab(
                   adata.obs[cluster_key],
                   adata.obs[SUBTYPE_KEY]
               )
               subtype_csv = os.path.join(
                   CSV_OUT,
                   f"subtype_composition_scVI{N_LATENT}"
                   f"_n{N_NEIGHBORS}"
                   f"_resolution_{res_label}.csv"
               )
               subtype_comp.to_csv(
                   subtype_csv
               )

               subtype_proportions = (
                   subtype_comp.div(
                       subtype_comp.sum(axis=1),
                       axis=0
                   )
               )
               dominant_subtype_fraction = (
                   subtype_proportions.max(axis=1)
               )
               subtype_dominance = (
                   float(
                       dominant_subtype_fraction.mean()
                   )
               )
               clusters_over_75pct_single_subtype = (
                   int(
                       (
                           dominant_subtype_fraction >= 0.75
                       ).sum()
                   )
               )
               clusters_over_90pct_single_subtype = (
                   int(
                       (
                           dominant_subtype_fraction >= 0.90
                       ).sum()
                   )
               )

           # =================================================
           # UMAP: LEIDEN
           # =================================================
           print("Creating Leiden UMAP...")
           sc.pl.embedding(
               adata,
               basis=umap_key.replace("X_", ""),
               color=cluster_key,
               legend_loc="on data",
               title=(
                   f"scVI {N_LATENT} | "
                   f"n={N_NEIGHBORS} | "
                   f"Leiden {resolution}"
               ),
               show=False
           )
           plt.savefig(
               os.path.join(
                   FIG_OUT,
                   f"umap_leiden_scVI{N_LATENT}"
                   f"_n{N_NEIGHBORS}"
                   f"_resolution_{res_label}.png"
               ),
               dpi=300,
               bbox_inches="tight"
           )
           plt.close()

           # =================================================
           # UMAP: SAMPLE + LEIDEN
           # =================================================
           sc.pl.embedding(
               adata,
               basis=umap_key.replace("X_", ""),
               color=[
                   SAMPLE_KEY,
                   cluster_key
               ],
               wspace=0.4,
               show=False
           )
           plt.savefig(
               os.path.join(
                   FIG_OUT,
                   f"umap_sample_vs_leiden_scVI{N_LATENT}"
                   f"_n{N_NEIGHBORS}"
                   f"_resolution_{res_label}.png"
               ),
               dpi=300,
               bbox_inches="tight"
           )
           plt.close()

           # =================================================
           # UMAP: SUBTYPE + LEIDEN
           # =================================================
           if SUBTYPE_KEY in adata.obs.columns:
               sc.pl.embedding(
                   adata,
                   basis=umap_key.replace("X_", ""),
                   color=[
                       SUBTYPE_KEY,
                       cluster_key
                   ],
                   wspace=0.4,
                   show=False
               )
               plt.savefig(
                   os.path.join(
                       FIG_OUT,
                       f"umap_subtype_vs_leiden_scVI{N_LATENT}"
                       f"_n{N_NEIGHBORS}"
                       f"_resolution_{res_label}.png"
                   ),
                   dpi=300,
                   bbox_inches="tight"
               )
               plt.close()

           # =================================================
           # STORE MASTER SUMMARY
           # =================================================
           parameter_results.append({
               "scVI_latent": N_LATENT,
               "n_neighbors": N_NEIGHBORS,
               "resolution": resolution,
               "n_clusters": n_clusters,
               "min_cluster_size":
                   min_cluster_size,
               "median_cluster_size":
                   median_cluster_size,
               "max_cluster_size":
                   max_cluster_size,
               "clusters_under_50":
                   clusters_under_50,
               "clusters_under_100":
                   clusters_under_100,
               "clusters_over_500":
                   clusters_over_500,
               "mean_dominant_sample_fraction":
                   mean_dominant_sample_fraction,
               "clusters_over_75pct_single_sample":
                   clusters_over_75pct_single_sample,
               "clusters_over_90pct_single_sample":
                   clusters_over_90pct_single_sample,
               "mean_dominant_subtype_fraction":
                   subtype_dominance,
               "clusters_over_75pct_single_subtype":
                   clusters_over_75pct_single_subtype,
               "clusters_over_90pct_single_subtype":
                   clusters_over_90pct_single_subtype,
               "cluster_key":
                   cluster_key,
               "cluster_counts_csv":
                   counts_csv,
               "sample_composition_csv":
                   sample_csv,
               "subtype_composition_csv":
                   subtype_csv
           })

   # ========================================================
   # SAVE H5AD
   # ========================================================
   output_h5ad = os.path.join(
       H5AD_OUT,
       f"scVI{N_LATENT}"
       f"_leiden_parameter_screen.h5ad"
   )
   print("\nSaving AnnData:")
   adata.write_h5ad(
       output_h5ad
   )
   print(output_h5ad)

# ############################################################
# 6. MASTER PARAMETER TABLE
# Write the combined metrics for all tested configurations.
# ############################################################
parameter_df = pd.DataFrame(
   parameter_results
)
parameter_output = os.path.join(
   CSV_OUT,
   "parameter_comparison_scVI25_vs_scVI30.csv"
)
parameter_df.to_csv(
   parameter_output,
   index=False
)

# ############################################################
# 7. DETAILED TEXT SUMMARY
# Write a short text record of the screening design.
# ############################################################
summary_output = os.path.join(
   CSV_OUT,
   "parameter_comparison_scVI25_vs_scVI30_summary.txt"
)
with open(summary_output, "w") as f:
   f.write(
       "scVI 25 vs scVI 30 Leiden Parameter Screen\n"
   )
   f.write(
       "============================================\n\n"
   )
   f.write(
       f"Total configurations tested: "
       f"{len(parameter_df)}\n\n"
   )
   f.write(
       "scVI latent dimensions: 25, 30\n"
   )
   f.write(
       "n_neighbors: 10, 15, 20, 30\n"
   )
   f.write(
       "Leiden resolutions: "
       f"{RESOLUTIONS}\n\n"
   )
   f.write(
       "The main comparison CSV contains one row "
       "per parameter combination.\n"
   )
   f.write(
       "Use this CSV to compare cluster numbers, "
       "cluster sizes, sample dominance and "
       "GIST subtype composition.\n"
   )

# ############################################################
# 8. FINISH
# Print the locations and scope of the completed screening outputs.
# ############################################################
print("\n============================================================")
print("scVI PARAMETER SCREENING COMPLETE")
print("============================================================")
print("\nConfigurations tested:")
print(
   f"  scVI latent dimensions: "
   f"{list(INPUTS.keys())}"
)
print(
   f"  n_neighbors: "
   f"{N_NEIGHBORS_LIST}"
)
print(
   f"  Leiden resolutions: "
   f"{RESOLUTIONS}"
)
print(
   f"\nTotal configurations: "
   f"{len(parameter_df)}"
)
print("\nMain comparison CSV:")
print(parameter_output)
print("\nSummary:")
print(summary_output)
print("\nFigures:")
print(FIG_OUT)
print("\nH5AD:")
print(H5AD_OUT)
print("\nCSVs:")
print(CSV_OUT)
print("\nDone.")