#!/usr/bin/env python
# ############################################################
# 1. SCRIPT SETUP
# Load the packages used by this analysis.
# ############################################################
import os
import scanpy as sc
import pandas as pd
import numpy as np
from sklearn.metrics import (
   adjusted_rand_score,
   normalized_mutual_info_score
)

# ############################################################
# 2. PATHS
# Set the two screened AnnData inputs and the validation output folder.
# ############################################################
BASE_DIR = "/home/cv1u24/GIST/5_Clustering"
SCVI25_H5AD = os.path.join(
   BASE_DIR,
   "Outputs",
   "Leiden",
   "h5ad",
   "scVI25_leiden_parameter_screen.h5ad"
)
SCVI30_H5AD = os.path.join(
   BASE_DIR,
   "Outputs",
   "Leiden",
   "h5ad",
   "scVI30_leiden_parameter_screen.h5ad"
)
OUT_DIR = os.path.join(
   BASE_DIR,
   "Outputs",
   "Leiden",
   "CSVs",
   "Marker_Analysis"
)
os.makedirs(
   OUT_DIR,
   exist_ok=True
)

# ############################################################
# 3. PARAMETERS
# Select the cluster labels and metadata used in comparisons.
# ############################################################
SCVI25_CLUSTER = "leiden_scVI25_n30_r0_6"
SCVI30_N20_CLUSTER = "leiden_scVI30_n20_r0_6"
SCVI30_CLUSTER = "leiden_scVI30_n30_r0_6"
SAMPLE_KEY = "sample"
SUBTYPE_KEY = "gist_subtype"
N_TOP_GENES = 50

# ############################################################
# 4. CHECK INPUT FILES
# Confirm both expected screened objects are available.
# ############################################################
print("\n============================================================")
print("CHECKING INPUT FILES")
print("============================================================")
for path in [
   SCVI25_H5AD,
   SCVI30_H5AD
]:
   if not os.path.exists(path):
       raise FileNotFoundError(
           f"\nFile not found:\n{path}"
       )
   print(
       f"Found: {path}"
   )

# ============================================================
# LOAD scVI25
# ============================================================
print("\n============================================================")
print("LOADING scVI25")
print("============================================================")
adata25 = sc.read_h5ad(
   SCVI25_H5AD
)
print(adata25)
print(
   f"Cells: {adata25.n_obs}"
)

# ============================================================
# LOAD scVI30
# ============================================================
print("\n============================================================")
print("LOADING scVI30")
print("============================================================")
adata30 = sc.read_h5ad(
   SCVI30_H5AD
)
print(adata30)
print(
   f"Cells: {adata30.n_obs}"
)

# ############################################################
# 5. CHECK CELL BARCODE OVERLAP
# Identify barcodes shared by both objects before comparing labels.
# ############################################################
print("\n============================================================")
print("CHECKING CELL OVERLAP")
print("============================================================")
cells25 = set(
   adata25.obs_names
)
cells30 = set(
   adata30.obs_names
)
common_cells = sorted(
   cells25.intersection(
       cells30
   )
)
print(
   f"scVI25 cells: "
   f"{len(cells25)}"
)
print(
   f"scVI30 cells: "
   f"{len(cells30)}"
)
print(
   f"Common cells: "
   f"{len(common_cells)}"
)
if len(common_cells) == 0:
   raise ValueError(
       "No overlapping cell barcodes were found."
   )

if len(common_cells) != len(cells25) or \
  len(common_cells) != len(cells30):
   print(
       "\nWARNING: The two datasets do not contain "
       "exactly the same cells."
   )
   print(
       "ARI/NMI will be calculated using only "
       "the common cells."
   )

# ############################################################
# 6. SUBSET AND ALIGN
# Subset both objects to the same cells and verify their order.
# ############################################################
adata25_common = adata25[
   common_cells
].copy()
adata30_common = adata30[
   common_cells
].copy()
print(
   f"\nAligned cells: "
   f"{adata25_common.n_obs}"
)
if not all(
   adata25_common.obs_names
   ==
   adata30_common.obs_names
):
   raise ValueError(
       "Cell order could not be aligned."
   )

# ============================================================
# PART 1
# ARI / NMI
# ============================================================
print("\n============================================================")
print("PART 1: ARI / NMI STABILITY")
print("============================================================")

clusterings = {
   "scVI25_n30_r0.6":
       adata25_common.obs[
           SCVI25_CLUSTER
       ].astype(str),
   "scVI30_n20_r0.6":
       adata30_common.obs[
           SCVI30_N20_CLUSTER
       ].astype(str),
   "scVI30_n30_r0.6":
       adata30_common.obs[
           SCVI30_CLUSTER
       ].astype(str)
}

stability_results = []

clustering_names = list(
   clusterings.keys()
)

for i in range(
   len(clustering_names)
):
   for j in range(
       i + 1,
       len(clustering_names)
   ):
       name1 = clustering_names[i]
       name2 = clustering_names[j]
       labels1 = clusterings[name1]
       labels2 = clusterings[name2]

       ari = adjusted_rand_score(
           labels1,
           labels2
       )
       nmi = normalized_mutual_info_score(
           labels1,
           labels2
       )

       print("\n--------------------------------------------")
       print(
           f"{name1}"
       )
       print(
           "vs"
       )
       print(
           f"{name2}"
       )
       print(
           f"ARI: {ari:.4f}"
       )
       print(
           f"NMI: {nmi:.4f}"
       )

       stability_results.append({
           "clustering_1":
               name1,
           "clustering_2":
               name2,
           "ARI":
               ari,
           "NMI":
               nmi,
           "clusters_1":
               len(
                   labels1.unique()
               ),
           "clusters_2":
               len(
                   labels2.unique()
               )
       })

# ############################################################
# 7. SAVE ARI/NMI
# Export pairwise adjusted Rand and normalized mutual information scores.
# ############################################################
stability_df = pd.DataFrame(
   stability_results
)
stability_csv = os.path.join(
   OUT_DIR,
   "scVI_ARI_NMI_stability_comparison.csv"
)
stability_df.to_csv(
   stability_csv,
   index=False
)

# ============================================================
# PART 2
# CLUSTER OVERLAP TABLES
# ============================================================
print("\n============================================================")
print("PART 2: CLUSTER OVERLAP")
print("============================================================")

for i in range(
   len(clustering_names)
):
   for j in range(
       i + 1,
       len(clustering_names)
   ):
       name1 = clustering_names[i]
       name2 = clustering_names[j]
       labels1 = clusterings[name1]
       labels2 = clusterings[name2]

       overlap = pd.crosstab(
           labels1,
           labels2
       )

       filename = (
           f"cluster_overlap_"
           f"{name1}_vs_{name2}.csv"
       )
       overlap.to_csv(
           os.path.join(
               OUT_DIR,
               filename
           )
       )

# ============================================================
# PART 3
# SAMPLE COMPOSITION OF LEADING CLUSTERING
# ============================================================
print("\n============================================================")
print("PART 3: LEADING CLUSTER SAMPLE COMPOSITION")
print("============================================================")

leading_sample = pd.crosstab(
   adata30_common.obs[
       SCVI30_CLUSTER
   ],
   adata30_common.obs[
       SAMPLE_KEY
   ]
)
leading_sample_csv = os.path.join(
   OUT_DIR,
   "scVI30_n30_r0_6_sample_composition.csv"
)
leading_sample.to_csv(
   leading_sample_csv
)

# ############################################################
# 8. SUBTYPE COMPOSITION
# Summarise subtype composition when the metadata column is present.
# ############################################################
if SUBTYPE_KEY in adata30_common.obs.columns:
   leading_subtype = pd.crosstab(
       adata30_common.obs[
           SCVI30_CLUSTER
       ],
       adata30_common.obs[
           SUBTYPE_KEY
       ]
   )
   leading_subtype_csv = os.path.join(
       OUT_DIR,
       "scVI30_n30_r0_6_subtype_composition.csv"
   )
   leading_subtype.to_csv(
       leading_subtype_csv
   )
else:
   leading_subtype_csv = None
   print(
       "\nWARNING: gist_subtype not found."
   )

# ============================================================
# PART 4
# MARKER GENE ANALYSIS
# ============================================================
print("\n============================================================")
print("PART 4: MARKER GENE ANALYSIS")
print("============================================================")
print(
   "\nLeading clustering:"
)
print(
   SCVI30_CLUSTER
)

# ============================================================
# COPY scVI30 OBJECT
# ============================================================
marker_adata = adata30_common.copy()

# ############################################################
# 9. CHECK EXPRESSION LAYER
# Use the log-normalised layer when available, otherwise retain the current matrix.
# ############################################################
print(
   "\nAvailable layers:"
)
print(
   list(
       marker_adata.layers.keys()
   )
)

if "lognorm" in marker_adata.layers:
   print(
       "\nUsing lognorm layer "
       "for marker analysis."
   )
   marker_adata.X = (
       marker_adata.layers[
           "lognorm"
       ].copy()
   )
else:
   print(
       "\nWARNING: lognorm layer not found."
   )
   print(
       "Using current adata.X."
   )

# ############################################################
# 10. RUN WILCOXON
# Rank genes by cluster using Scanpy’s Wilcoxon method.
# ############################################################
print(
   "\nRunning Wilcoxon "
   "rank-gene analysis..."
)
sc.tl.rank_genes_groups(
   marker_adata,
   groupby=SCVI30_CLUSTER,
   method="wilcoxon",
   key_added="markers_wilcoxon",
   pts=True
)

# ############################################################
# 11. EXTRACT ALL MARKERS
# Convert the ranked marker results into a table.
# ############################################################
marker_df = sc.get.rank_genes_groups_df(
   marker_adata,
   group=None,
   key="markers_wilcoxon"
)

# ############################################################
# 12. SAVE ALL MARKERS
# Export the complete marker table.
# ############################################################
all_markers_csv = os.path.join(
   OUT_DIR,
   "scVI30_n30_r0_6_all_markers.csv"
)
marker_df.to_csv(
   all_markers_csv,
   index=False
)

# ############################################################
# 13. TOP 50 MARKERS
# Select the 50 highest-scoring genes per cluster.
# ############################################################
top50 = (
   marker_df
   .sort_values(
       [
           "group",
           "scores"
       ],
       ascending=[
           True,
           False
       ]
   )
   .groupby(
       "group",
       observed=True
   )
   .head(
       N_TOP_GENES
   )
)
top50_csv = os.path.join(
   OUT_DIR,
   "scVI30_n30_r0_6_top50_markers.csv"
)
top50.to_csv(
   top50_csv,
   index=False
)

# ############################################################
# 14. SIGNIFICANT MARKERS
# Filter marker results using adjusted p-value below 0.05.
# ############################################################
significant = marker_df[
   marker_df[
       "pvals_adj"
   ] < 0.05
].copy()
significant_csv = os.path.join(
   OUT_DIR,
   "scVI30_n30_r0_6_significant_markers.csv"
)
significant.to_csv(
   significant_csv,
   index=False
)

# ############################################################
# 15. TOP 10 MARKERS
# Select the 10 highest-scoring genes per cluster.
# ############################################################
top10 = (
   marker_df
   .sort_values(
       [
           "group",
           "scores"
       ],
       ascending=[
           True,
           False
       ]
   )
   .groupby(
       "group",
       observed=True
   )
   .head(10)
)
top10_csv = os.path.join(
   OUT_DIR,
   "scVI30_n30_r0_6_top10_markers.csv"
)
top10.to_csv(
   top10_csv,
   index=False
)

# ############################################################
# 16. MARKER SUMMARY
# Summarise significant marker counts and top genes by cluster.
# ############################################################
summary_rows = []

for cluster in sorted(
   marker_df[
       "group"
   ].astype(str).unique(),
   key=lambda x: int(x)
):
   cluster_markers = marker_df[
       marker_df[
           "group"
       ].astype(str)
       == cluster
   ].copy()

   significant_cluster = (
       cluster_markers[
           cluster_markers[
               "pvals_adj"
           ] < 0.05
       ]
   )

   top_genes = (
       cluster_markers
       .sort_values(
           "scores",
           ascending=False
       )
       .head(10)[
           "names"
       ]
       .astype(str)
       .tolist()
   )

   summary_rows.append({
       "cluster":
           cluster,
       "n_significant_markers":
           len(
               significant_cluster
           ),
       "top_10_markers":
           "; ".join(
               top_genes
           )
   })

marker_summary = pd.DataFrame(
   summary_rows
)

marker_summary_csv = os.path.join(
   OUT_DIR,
   "scVI30_n30_r0_6_marker_summary.csv"
)
marker_summary.to_csv(
   marker_summary_csv,
   index=False
)

# ############################################################
# 17. FINAL COMBINED SUMMARY
# Write a text summary of cell overlap and marker-analysis settings.
# ############################################################
final_summary = os.path.join(
   OUT_DIR,
   "ANALYSIS_README.txt"
)
with open(
   final_summary,
   "w"
) as f:
   f.write(
       "scVI clustering stability and marker analysis\n"
   )
   f.write(
       "=============================================\n\n"
   )
   f.write(
       f"Common cells analysed: "
       f"{len(common_cells)}\n\n"
   )
   f.write(
       "ARI/NMI comparisons:\n"
   )
   f.write(
       "  scVI25 n30 r0.6\n"
   )
   f.write(
       "  scVI30 n20 r0.6\n"
   )
   f.write(
       "  scVI30 n30 r0.6\n\n"
   )
   f.write(
       "Leading clustering for marker analysis:\n"
   )
   f.write(
       "  scVI30 n30 r0.6\n\n"
   )
   f.write(
       "Marker method:\n"
   )
   f.write(
       "  Wilcoxon rank-sum test\n"
   )
   f.write(
       "  Expression layer: lognorm\n"
   )

# ############################################################
# 18. FINISH
# Print the main output paths.
# ############################################################
print("\n============================================================")
print("ANALYSIS COMPLETE")
print("============================================================")
print(
   "\nCommon cells analysed:",
   len(common_cells)
)
print("\nARI/NMI results:")
print(stability_csv)
print("\nCluster overlap tables:")
print(OUT_DIR)
print("\nLeading sample composition:")
print(leading_sample_csv)
if leading_subtype_csv:
   print("\nLeading subtype composition:")
   print(leading_subtype_csv)
print("\nAll marker genes:")
print(all_markers_csv)
print("\nTop 50 markers:")
print(top50_csv)
print("\nSignificant markers:")
print(significant_csv)
print("\nTop 10 markers:")
print(top10_csv)
print("\nMarker summary:")
print(marker_summary_csv)
print("\nDone.")
