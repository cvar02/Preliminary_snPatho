#!/usr/bin/env python
# ############################################################
# 1. SCRIPT SETUP
# Load the packages used by this analysis.
# ############################################################
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt
# ############################################################
# 2. PATHS
# Set the clustering input and annotation output folders.
# ############################################################
INPUT_H5AD = (
   "/home/cv1u24/GIST/5_Clustering/"
   "Outputs/Leiden/h5ad/"
   "scVI30_leiden_parameter_screen.h5ad"
)
BASE_OUT = (
   "/home/cv1u24/GIST/5_Clustering/"
   "Outputs/ClusterAnnotation"
)
FIG_OUT = os.path.join(BASE_OUT, "Figures")
CSV_OUT = os.path.join(BASE_OUT, "CSVs")
H5AD_OUT = os.path.join(BASE_OUT, "h5ad")
os.makedirs(FIG_OUT, exist_ok=True)
os.makedirs(CSV_OUT, exist_ok=True)
os.makedirs(H5AD_OUT, exist_ok=True)
# ############################################################
# 3. PARAMETERS
# Select the Leiden cluster key and marker-ranking key.
# ############################################################
CLUSTER_KEY = "leiden_scVI30_n30_r0_6"
MARKER_KEY = "markers_leiden_scVI30_n30_r0_6"
# Number of top markers to export
TOP_MARKERS = [10, 20, 50, 200]
# ############################################################
# 4. LOAD
# Read the clustering object and confirm the selected cluster labels exist.
# ############################################################
print("")
print("Loading clustering object...")
adata = sc.read_h5ad(INPUT_H5AD)
print(adata)
if CLUSTER_KEY not in adata.obs:
   raise ValueError(
       f"{CLUSTER_KEY} not found in adata.obs"
   )
# ############################################################
# 5. CURATED LITERATURE-INFORMED MARKER SETS
# Define marker panels used to inspect broad cell populations.
# ############################################################
marker_genes = {
   # --------------------------------------------------------
   # GIST tumour
   # --------------------------------------------------------
   "GIST_tumour": [
       "KIT",
       "PDGFRA",
       "ANO1",
       "PTP4A3",
       "MOXD1"
   ],
   # --------------------------------------------------------
   # T cells
   # --------------------------------------------------------
   "T_cells": [
       "CD3D",
       "CD3E",
       "CD3G",
       "TRBC1",
       "TRBC2",
       "IL7R",
       "LTB"
   ],
   # --------------------------------------------------------
   # Cytotoxic T cells
   # --------------------------------------------------------
   "Cytotoxic_T_cells": [
       "CD8A",
       "CD8B",
       "NKG7",
       "GNLY",
       "GZMK",
       "GZMB",
       "PRF1"
   ],
   # --------------------------------------------------------
   # NK cells
   # --------------------------------------------------------
   "NK_cells": [
       "NKG7",
       "GNLY",
       "KLRD1",
       "FCGR3A",
       "TYROBP"
   ],
   # --------------------------------------------------------
   # B cells
   # --------------------------------------------------------
   "B_cells": [
       "CD79A",
       "MS4A1",
       "CD37",
       "CD74",
       "HLA-DRA"
   ],
   # --------------------------------------------------------
   # Plasma cells
   # --------------------------------------------------------
   "Plasma_cells": [
       "MZB1",
       "XBP1",
       "JCHAIN",
       "IGKC",
       "IGHG1",
       "CD38"
   ],
   # --------------------------------------------------------
   # Broad myeloid
   # --------------------------------------------------------
   "Myeloid": [
       "LYZ",
       "LST1",
       "CTSS",
       "FCER1G",
       "TYROBP",
       "SPI1"
   ],
   # --------------------------------------------------------
   # Monocytes
   # --------------------------------------------------------
   "Monocytes": [
       "LYZ",
       "S100A8",
       "S100A9",
       "FCN1",
       "CTSD",
       "LGALS3"
   ],
   # --------------------------------------------------------
   # Macrophages
   # --------------------------------------------------------
   "Macrophages": [
       "C1QA",
       "C1QB",
       "C1QC",
       "APOE",
       "CD68",
       "MSR1"
   ],
   # --------------------------------------------------------
   # Dendritic cells
   # --------------------------------------------------------
   "Dendritic_cells": [
       "FCER1A",
       "CD1C",
       "CLEC10A",
       "CST3"
   ],
   # --------------------------------------------------------
   # Neutrophils
   # --------------------------------------------------------
   "Neutrophils": [
       "S100A8",
       "S100A9",
       "FCGR3B",
       "CSF3R",
       "CXCR2",
       "MPO"
   ],
   # --------------------------------------------------------
   # Mast cells
   # --------------------------------------------------------
   "Mast_cells": [
       "TPSAB1",
       "TPSB2",
       "MS4A2",
       "HDC",
       "CPA3"
   ],
   # --------------------------------------------------------
   # Endothelial cells
   # --------------------------------------------------------
   "Endothelial": [
       "VWF",
       "PECAM1",
       "CDH5",
       "KDR",
       "FLT1",
       "EMCN",
       "ESAM",
       "EGFL7"
   ],
   # --------------------------------------------------------
   # Pericytes
   # --------------------------------------------------------
   "Pericytes": [
       "RGS5",
       "PDGFRB",
       "CSPG4",
       "MCAM",
       "COL4A1",
       "COL4A2"
   ],
   # --------------------------------------------------------
   # Fibroblasts
   # --------------------------------------------------------
   "Fibroblasts": [
       "COL1A1",
       "COL1A2",
       "COL3A1",
       "COL5A1",
       "COL6A1",
       "COL6A2",
       "DCN",
       "LUM"
   ],
   # --------------------------------------------------------
   # Smooth muscle
   # --------------------------------------------------------
   "Smooth_muscle": [
       "ACTA2",
       "TAGLN",
       "MYL9",
       "CNN1",
       "DES"
   ],
   # --------------------------------------------------------
   # Epithelial
   # --------------------------------------------------------
   "Epithelial": [
       "EPCAM",
       "KRT8",
       "KRT18",
       "KRT19",
       "KRT5",
       "KRT15",
       "TRIM29",
       "DSP"
   ]
}
# ############################################################
# 6. CHECK MARKER AVAILABILITY
# Retain and report only curated genes present in the object.
# ############################################################
available_genes = set(adata.var_names)
filtered_marker_genes = {}
print("")
print("Checking curated marker availability...")
for cell_type, genes in marker_genes.items():
   present = [
       gene for gene in genes
       if gene in available_genes
   ]
   filtered_marker_genes[cell_type] = present
   print(
       f"{cell_type}: "
       f"{len(present)}/{len(genes)} genes present"
   )
# ############################################################
# 7. FIND UNBIASED CLUSTER MARKERS
# Rank genes by cluster using the Wilcoxon method.
# ############################################################
print("")
print("Finding cluster marker genes...")
sc.tl.rank_genes_groups(
   adata,
   groupby=CLUSTER_KEY,
   method="wilcoxon",
   key_added=MARKER_KEY
)
# ############################################################
# 8. EXPORT ALL MARKERS
# Export the full ranked marker table.
# ############################################################
print("")
print("Exporting marker table...")
marker_rows = []
groups = adata.obs[CLUSTER_KEY].cat.categories
for cluster in groups:
   genes = adata.uns[MARKER_KEY]["names"][cluster]
   scores = adata.uns[MARKER_KEY]["scores"][cluster]
   pvals = adata.uns[MARKER_KEY]["pvals_adj"][cluster]
   for rank, (gene, score, pval) in enumerate(
       zip(genes, scores, pvals),
       start=1
   ):
       marker_rows.append({
           "cluster": cluster,
           "rank": rank,
           "gene": gene,
           "score": score,
           "adjusted_pvalue": pval
       })
markers_df = pd.DataFrame(marker_rows)
markers_csv = os.path.join(
   CSV_OUT,
   "scVI30_n30_r0_6_cluster_markers.csv"
)
markers_df.to_csv(
   markers_csv,
   index=False
)
# ############################################################
# 9. EXPORT TOP MARKERS
# Write separate ranked-marker tables at the selected cut-offs.
# ############################################################
for n in TOP_MARKERS:
   top_df = markers_df[
       markers_df["rank"] <= n
   ]
   output = os.path.join(
       CSV_OUT,
       f"scVI30_n30_r0_6_top{n}_markers.csv"
   )
   top_df.to_csv(
       output,
       index=False
   )
# ############################################################
# 10. TOP 10 MARKER SUMMARY
# Create a compact table of the top ten genes per cluster.
# ############################################################
top10 = (
   markers_df
   .query("rank <= 10")
   .groupby("cluster")["gene"]
   .apply(lambda x: "; ".join(x))
   .reset_index()
)
top10.columns = [
   "cluster",
   "top_10_markers"
]
top10_csv = os.path.join(
   CSV_OUT,
   "scVI30_n30_r0_6_top10_marker_summary.csv"
)
top10.to_csv(
   top10_csv,
   index=False
)
# ############################################################
# 11. CURATED MARKER DOTPLOT
# Plot expression of the available curated marker panels by cluster.
# ############################################################
print("")
print("Creating curated marker DotPlot...")
dot_genes = []
for genes in filtered_marker_genes.values():
   for gene in genes:
       if gene not in dot_genes:
           dot_genes.append(gene)
if len(dot_genes) > 0:
   sc.pl.dotplot(
       adata,
       var_names=filtered_marker_genes,
       groupby=CLUSTER_KEY,
       standard_scale="var",
       dendrogram=False,
       show=False
   )
   plt.savefig(
       os.path.join(
           FIG_OUT,
           "scVI30_n30_r0_6_curated_marker_dotplot.png"
       ),
       dpi=300,
       bbox_inches="tight"
   )
   plt.close()
# ############################################################
# 12. TOP 10 MARKER HEATMAP
# Plot the top ranked genes for each cluster.
# ############################################################
print("Creating top marker heatmap...")
sc.pl.rank_genes_groups_heatmap(
   adata,
   groupby=CLUSTER_KEY,
   key=MARKER_KEY,
   n_genes=10,
   show=False
)
plt.savefig(
   os.path.join(
       FIG_OUT,
       "scVI30_n30_r0_6_top10_marker_heatmap.png"
   ),
   dpi=300,
   bbox_inches="tight"
)
plt.close()
# ============================================================
# UMAP: CLUSTERS
# ============================================================
print("Creating cluster UMAP...")
sc.pl.umap(
   adata,
   color=CLUSTER_KEY,
   legend_loc="on data",
   title="scVI30 n30 r0.6 Leiden clusters",
   show=False
)
plt.savefig(
   os.path.join(
       FIG_OUT,
       "scVI30_n30_r0_6_UMAP_clusters.png"
   ),
   dpi=300,
   bbox_inches="tight"
)
plt.close()
# ============================================================
# UMAP: SAMPLE
# ============================================================
if "sample" in adata.obs:
   print("Creating sample UMAP...")
   sc.pl.umap(
       adata,
       color="sample",
       title="Sample",
       show=False
   )
   plt.savefig(
       os.path.join(
           FIG_OUT,
           "scVI30_n30_r0_6_UMAP_sample.png"
       ),
       dpi=300,
       bbox_inches="tight"
   )
   plt.close()
# ============================================================
# UMAP: GIST SUBTYPE
# ============================================================
subtype_column = None
for col in [
   "gist_subtype",
   "genotype_group",
   "subtype"
]:
   if col in adata.obs:
       subtype_column = col
       break
if subtype_column:
   print(
       f"Creating subtype UMAP using {subtype_column}..."
   )
   sc.pl.umap(
       adata,
       color=subtype_column,
       title="GIST subtype",
       show=False
   )
   plt.savefig(
       os.path.join(
           FIG_OUT,
           "scVI30_n30_r0_6_UMAP_subtype.png"
       ),
       dpi=300,
       bbox_inches="tight"
   )
   plt.close()
# ############################################################
# 13. SAVE OBJECT
# Write the AnnData object containing the marker-ranking results.
# ############################################################
output_h5ad = os.path.join(
   H5AD_OUT,
   "scVI30_n30_r0_6_cluster_annotation.h5ad"
)
adata.write_h5ad(output_h5ad)
# ############################################################
# 14. FINISH
# Print the output folders and key tables.
# ############################################################
print("")
print("========================================")
print("CLUSTER ANNOTATION COMPLETE")
print("========================================")
print("")
print("Figures:")
print(FIG_OUT)
print("")
print("CSVs:")
print(CSV_OUT)
print("")
print("H5AD:")
print(output_h5ad)
print("")
print("Key files:")
print(markers_csv)
print(top10_csv)
print("")
print("Done.")