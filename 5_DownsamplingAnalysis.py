#!/usr/bin/env python3
# ############################################################
# 1. PURPOSE AND PACKAGES
# Balance subtype cell counts, recluster the selected cells, and summarise every cluster.
# ############################################################
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
# ############################################################
# 2. INPUT AND OUTPUT LOCATIONS
# Read the integrated scVI30 object and create folders for tables, plots, and the saved object.
# ############################################################
INPUT_H5AD = (
   "/iridisfs/home/cv1u24/GIST/5_Clustering/"
   "Outputs/Leiden/h5ad/"
   "scVI30_leiden_parameter_screen.h5ad"
)
BASE_OUT = (
   "/iridisfs/home/cv1u24/GIST/5_Clustering/"
   "Outputs/Downsampling_Sensitivity"
)
FIG_OUT = os.path.join(BASE_OUT, "Figures")
CSV_OUT = os.path.join(BASE_OUT, "CSVs")
H5AD_OUT = os.path.join(BASE_OUT, "h5ad")
os.makedirs(FIG_OUT, exist_ok=True)
os.makedirs(CSV_OUT, exist_ok=True)
os.makedirs(H5AD_OUT, exist_ok=True)
# ############################################################
# 3. ANALYSIS SETTINGS
# Define one reproducible random draw, equal cells per subtype, and the graph clustering settings.
# ############################################################
SUBTYPE_KEY = "gist_subtype"
# Fixed target used for this one random draw; each subtype contributes this many nuclei.
TARGET_N = 839
# Use a 30-neighbour graph and Leiden resolution 0.6 on scVI30 coordinates.
N_NEIGHBORS = 30
LEIDEN_RESOLUTION = 0.6
# Use the existing scVI representation.
SCVI_KEY = "X_scVI"
# Dimension count corresponding to the existing 30-dimensional scVI representation.
N_DIMS = 30
RANDOM_SEED = 1
# ############################################################
# 4. LOAD AND CHECK THE INPUT
# Confirm subtype, sample, and scVI fields are present before analysis.
# ############################################################
print("")
print("Loading AnnData...")
adata = sc.read_h5ad(INPUT_H5AD)
print(adata)
if SUBTYPE_KEY not in adata.obs:
   raise ValueError(
       "gist_subtype not found in adata.obs"
   )
if SCVI_KEY not in adata.obsm:
   raise ValueError(
       "X_scVI not found in adata.obsm"
   )
if "sample" not in adata.obs.columns:
   raise ValueError(
       "sample not found in adata.obs"
   )
# ############################################################
# 5. COUNT NUCLEI BEFORE BALANCING
# Show how many nuclei each subtype contributes to the input object.
# ############################################################
print("")
print("Original subtype counts:")
original_counts = (
   adata.obs[SUBTYPE_KEY]
   .value_counts()
   .sort_index()
)
print(original_counts)
# ############################################################
# 6. DRAW AN EQUAL SUBTYPE SAMPLE
# Randomly select TARGET_N nuclei from each subtype; samples are pooled, not balanced separately.
# ############################################################
print("")
print("Performing balanced downsampling...")
rng = np.random.default_rng(RANDOM_SEED)
selected_cells = []
for subtype in original_counts.index:
   subtype_cells = adata.obs_names[
       adata.obs[SUBTYPE_KEY] == subtype
   ].to_numpy()
   if len(subtype_cells) < TARGET_N:
       raise ValueError(
           f"{subtype} has only "
           f"{len(subtype_cells)} nuclei. "
           f"Cannot sample {TARGET_N}."
       )
   selected = rng.choice(
       subtype_cells,
       size=TARGET_N,
       replace=False
   )
   selected_cells.extend(selected)
selected_cells = np.array(selected_cells)
adata_balanced = adata[selected_cells].copy()
print("")
print("Balanced subtype counts:")
balanced_counts = (
   adata_balanced.obs[SUBTYPE_KEY]
   .value_counts()
   .sort_index()
)
print(balanced_counts)
# ############################################################
# 7. RECORD THE SELECTED NUCLEI
# Save barcodes and their subtype/sample labels so this single seeded draw can be reproduced.
# ############################################################
selection_df = pd.DataFrame({
   "barcode": adata_balanced.obs_names,
   "gist_subtype": adata_balanced.obs[SUBTYPE_KEY].values,
   "sample": adata_balanced.obs["sample"].values
})
selection_csv = os.path.join(
   CSV_OUT,
   "balanced_downsampling_cell_selection.csv"
)
selection_df.to_csv(
   selection_csv,
   index=False
)
# ############################################################
# 8. REBUILD THE NEIGHBOUR GRAPH
# Recalculate cell-to-cell neighbours in the existing scVI latent representation.
# ############################################################
print("")
print("Building balanced neighbour graph...")
sc.pp.neighbors(
   adata_balanced,
   use_rep=SCVI_KEY,
   n_neighbors=N_NEIGHBORS,
   n_pcs=N_DIMS,
   random_state=RANDOM_SEED
)
# ############################################################
# 9. CLUSTER ALL SELECTED NUCLEI
# Run Leiden on the complete balanced subset, not on one sample or one cluster.
# ############################################################
print("")
print("Running balanced Leiden clustering...")
balanced_key = "leiden_balanced_n30_r0_6"
sc.tl.leiden(
   adata_balanced,
   resolution=LEIDEN_RESOLUTION,
   key_added=balanced_key,
   random_state=RANDOM_SEED
)
# ############################################################
# 10. CALCULATE A UMAP FOR THE FULL BALANCED SET
# Create one embedding that will be coloured by cluster, subtype, sample, and candidate genes.
# ############################################################
print("")
print("Calculating UMAP...")
sc.tl.umap(
   adata_balanced,
   random_state=RANDOM_SEED
)
# ############################################################
# 11. SUMMARISE EVERY CLUSTER
# Report the nucleus count and percentage for each Leiden cluster.
# ############################################################
cluster_counts = (
   adata_balanced.obs[balanced_key]
   .value_counts()
   .sort_index()
)
cluster_summary = pd.DataFrame({
   "cluster": cluster_counts.index,
   "n_nuclei": cluster_counts.values,
   "percent": (
       cluster_counts.values /
       len(adata_balanced) *
       100
   )
})
cluster_csv = os.path.join(
   CSV_OUT,
   "balanced_cluster_sizes.csv"
)
cluster_summary.to_csv(
   cluster_csv,
   index=False
)
print("")
print("Balanced cluster sizes:")
print(cluster_summary)
# ############################################################
# 12. CHECK SUBTYPE MIXTURE IN EACH CLUSTER
# Calculate the percentage of each subtype represented within every cluster.
# ############################################################
composition = pd.crosstab(
   adata_balanced.obs[balanced_key],
   adata_balanced.obs[SUBTYPE_KEY],
   normalize="index"
) * 100
composition_csv = os.path.join(
   CSV_OUT,
   "balanced_cluster_subtype_composition.csv"
)
composition.to_csv(composition_csv)
print("")
print("Cluster subtype composition:")
print(composition)
# ############################################################
# 13. CHECK SAMPLE MIXTURE IN EACH CLUSTER
# Show which samples contribute nuclei to every cluster; this is not sample-level balancing.
# ############################################################
sample_composition = pd.crosstab(
   adata_balanced.obs[balanced_key],
   adata_balanced.obs["sample"],
   normalize="index"
) * 100
sample_csv = os.path.join(
   CSV_OUT,
   "balanced_cluster_sample_composition.csv"
)
sample_composition.to_csv(sample_csv)
# ############################################################
# 14. RANK MARKERS FOR ALL CLUSTERS
# Use Wilcoxon ranking to identify genes enriched in each cluster of the balanced subset.
# ############################################################
print("")
print("Running differential expression analysis...")
marker_key = "balanced_cluster_markers"
sc.tl.rank_genes_groups(
   adata_balanced,
   groupby=balanced_key,
   method="wilcoxon",
   key_added=marker_key
)
markers_df = sc.get.rank_genes_groups_df(
   adata_balanced,
   group=None,
   key=marker_key
).rename(
   columns={"group": "cluster"}
)
# ############################################################
# 15. SAVE FULL MARKER RANKINGS FOR ALL CLUSTERS
# Export every ranked gene for every cluster in the balanced dataset.
# ############################################################
all_markers_csv = os.path.join(
   CSV_OUT,
   "balanced_all_cluster_ranked_markers.csv"
)
markers_df.to_csv(
   all_markers_csv,
   index=False
)
# ############################################################
# 16. SELECT THE TOP 10 MARKERS FOR EACH CLUSTER
# Keep the ten highest-ranked genes for every cluster, with statistics.
# ############################################################
print("")
print("Extracting top 10 DEGs for each cluster...")
marker_rows = []
groups = (
   adata_balanced.obs[balanced_key]
   .cat.categories
)
for cluster in groups:
   genes = adata_balanced.uns[
       marker_key
   ]["names"][cluster]
   scores = adata_balanced.uns[
       marker_key
   ]["scores"][cluster]
   pvals = adata_balanced.uns[
       marker_key
   ]["pvals_adj"][cluster]
   # logfoldchanges may not always be available
   if "logfoldchanges" in adata_balanced.uns[marker_key]:
       logfc = adata_balanced.uns[
           marker_key
       ]["logfoldchanges"][cluster]
   else:
       logfc = [np.nan] * len(genes)
   for rank, (gene, score, pval, fc) in enumerate(
       zip(genes, scores, pvals, logfc),
       start=1
   ):
       if rank > 10:
           break
       marker_rows.append({
           "cluster": cluster,
           "rank": rank,
           "gene": gene,
           "score": score,
           "logfoldchange": fc,
           "adjusted_pvalue": pval
       })
markers_top10_df = pd.DataFrame(marker_rows)
# ############################################################
# 17. SAVE TOP 10 MARKERS FOR EVERY CLUSTER
# Write one row per cluster-marker pair, including score and adjusted p-value.
# ############################################################
top10_csv = os.path.join(
   CSV_OUT,
   "balanced_cluster_top10_DEGs.csv"
)
markers_top10_df.to_csv(
   top10_csv,
   index=False
)
print("")
print("Top 10 DEGs saved to:")
print(top10_csv)
# ############################################################
# 18. SAVE A COMPACT TOP 10 TABLE
# Put one cluster per row and its ranked marker genes in separate columns.
# ############################################################
top10_wide = (
   markers_top10_df
   .pivot(
       index="cluster",
       columns="rank",
       values="gene"
   )
   .reset_index()
)
top10_wide.columns = [
   "cluster"
   if col == "cluster"
   else f"DEG_{int(col)}"
   for col in top10_wide.columns
]
top10_wide_csv = os.path.join(
   CSV_OUT,
   "balanced_cluster_top10_DEGs_wide.csv"
)
top10_wide.to_csv(
   top10_wide_csv,
   index=False
)
print("")
print("Wide-format top 10 DEG table saved to:")
print(top10_wide_csv)
# ############################################################
# 19. PRINT TOP MARKERS FOR EVERY CLUSTER
# List the ten highest-ranked genes for each cluster in the run log.
# ############################################################
print("")
print("==============================================")
print("TOP 10 DEGs PER CLUSTER")
print("==============================================")
for cluster in groups:
   cluster_genes = markers_top10_df[
       markers_top10_df["cluster"].astype(str)
       == str(cluster)
   ]
   print("")
   print(f"Cluster {cluster}:")
   for _, row in cluster_genes.iterrows():
       print(
           f"  {int(row['rank'])}. "
           f"{row['gene']} "
           f"(logFC={row['logfoldchange']:.2f}, "
           f"adj.P={row['adjusted_pvalue']:.2e})"
       )
# ############################################################
# 20. VISUALISE TOP MARKERS FOR EVERY CLUSTER
# Plot the top ten ranked genes for each cluster in one heatmap.
# ############################################################
print("")
print("Creating top-marker heatmap for all clusters...")
sc.pl.rank_genes_groups_heatmap(
   adata_balanced,
   groupby=balanced_key,
   key=marker_key,
   n_genes=10,
   show=False
)
plt.savefig(
   os.path.join(
       FIG_OUT,
       "balanced_all_clusters_top10_marker_heatmap.png"
   ),
   dpi=300,
   bbox_inches="tight"
)
plt.close()
# ############################################################
# 21. UMAP OF ALL LEIDEN CLUSTERS
# Display all nuclei in the balanced set, coloured by their Leiden cluster.
# ############################################################
print("")
print("Creating balanced cluster UMAP...")
sc.pl.umap(
   adata_balanced,
   color=balanced_key,
   legend_loc="on data",
   title="Balanced Leiden clustering",
   frameon=False,
   show=False
)
plt.savefig(
   os.path.join(
       FIG_OUT,
       "balanced_Leiden_clusters.png"
   ),
   dpi=300,
   bbox_inches="tight"
)
plt.close()
# ############################################################
# 22A. UMAP COLOURED BY SUBTYPE
# Show subtype locations across the full balanced-cell embedding.
# ############################################################
sc.pl.umap(
   adata_balanced,
   color=SUBTYPE_KEY,
   title="GIST subtype",
   frameon=False,
   show=False
)
plt.savefig(
   os.path.join(
       FIG_OUT,
       "balanced_UMAP_subtype.png"
   ),
   dpi=300,
   bbox_inches="tight"
)
plt.close()
# ############################################################
# 22B. UMAP COLOURED BY SAMPLE
# Show how samples contribute to clusters in the full balanced-cell embedding.
# ############################################################
sc.pl.umap(
   adata_balanced,
   color="sample",
   title="Sample",
   frameon=False,
   show=False
)
plt.savefig(
   os.path.join(
       FIG_OUT,
       "balanced_UMAP_sample.png"
   ),
   dpi=300,
   bbox_inches="tight"
)
plt.close()
# ############################################################
# 22C. UMAPS OF CANDIDATE MARKER EXPRESSION
# Plot each available candidate gene across all selected nuclei; these are not cluster-specific subsets.
# ############################################################
candidate_genes = [
   "KIT",
   "ANO1",
   "MS4A1",
   "CD79A",
   "IGHM",
   "IGHD",
   "BANK1",
   "AIM2",
   "VCAN",
   "ITGAX",
   "SELL"
]
present_genes = [
   gene for gene in candidate_genes
   if gene in adata_balanced.var_names
]
if len(present_genes) > 0:
   print("")
   print(
       "Plotting candidate markers across the balanced dataset:"
   )
   sc.pl.umap(
       adata_balanced,
       color=present_genes,
       ncols=3,
       frameon=False,
       show=False
   )
   plt.savefig(
       os.path.join(
           FIG_OUT,
           "balanced_candidate_marker_expression.png"
       ),
       dpi=300,
       bbox_inches="tight"
   )
   plt.close()
# ############################################################
# 23. SAVE THE FULL BALANCED ANALYSIS OBJECT
# Save all selected nuclei, Leiden labels, UMAP coordinates, and marker results.
# ############################################################
output_h5ad = os.path.join(
   H5AD_OUT,
   "scVI30_balanced_downsampling.h5ad"
)
adata_balanced.write_h5ad(
   output_h5ad
)
# ############################################################
# 24. SUMMARISE THE COMPLETE ANALYSIS
# Report balance counts, number of clusters, and output locations.
# ############################################################
print("")
print("==============================================")
print("BALANCED DOWNSAMPLING ANALYSIS COMPLETE")
print("==============================================")
print("")
print("Original subtype counts:")
print(original_counts)
print("")
print("Balanced subtype counts:")
print(balanced_counts)
print("")
print(
   f"Number of Leiden clusters: "
   f"{adata_balanced.obs[balanced_key].nunique()}"
)
print(
   f"Total nuclei analysed: "
   f"{len(adata_balanced)}"
)
print("Top 10 marker table (all clusters):")
print(top10_csv)
print("Full ranked marker table (all clusters):")
print(all_markers_csv)
print("")
print("Outputs:")
print(FIG_OUT)
print(CSV_OUT)
print(output_h5ad)
print("")
print("Done.")
