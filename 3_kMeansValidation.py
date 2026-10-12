#!/usr/bin/env python
# ############################################################
# 1. SCRIPT SETUP
# Load the packages used by this analysis.
# ############################################################
import os
import scanpy as sc
import pandas as pd
from sklearn.cluster import KMeans
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
# ############################################################
# 2. PATHS
# Set the screened scVI30 input and K-means comparison CSV path.
# ############################################################
INPUT_H5AD = (
   "/home/cv1u24/GIST/5_Clustering/"
   "Outputs/Leiden/h5ad/"
   "scVI30_leiden_parameter_screen.h5ad"
)
OUTPUT_CSV = (
   "/home/cv1u24/GIST/5_Clustering/"
   "Outputs/kMeans/CSVs/"
   "KMeans_vs_scVI30_n30_r0_6_ARI_NMI.csv"
)
# ############################################################
# 3. LOAD
# Read the scVI30 Leiden object.
# ############################################################
print("\nLoading scVI30 Leiden object...")
adata = sc.read_h5ad(INPUT_H5AD)
print("Loaded successfully:")
print(adata)
# ############################################################
# 4. CHECK DATA
# Select the scVI latent representation and reference Leiden labels.
# ############################################################
X = adata.obsm["X_scVI"]
REFERENCE = "leiden_scVI30_n30_r0_6"
reference_clusters = adata.obs[REFERENCE].astype(str)
print("\n----------------------------------------")
print("Validation data")
print("----------------------------------------")
print(f"Cells: {X.shape[0]}")
print(f"scVI dimensions: {X.shape[1]}")
print(f"Reference: {REFERENCE}")
print(f"Reference clusters: {reference_clusters.nunique()}")
# ############################################################
# 5. K-MEANS
# Fit K-means across the specified cluster counts and compare labels.
# ############################################################
K_VALUES = [8, 10, 12, 14, 16, 18]
results = []
print("\n----------------------------------------")
print("Running K-means")
print("----------------------------------------")
for k in K_VALUES:
   print(f"K = {k}")
   kmeans = KMeans(
       n_clusters=k,
       random_state=0,
       n_init=20
   )
   kmeans_labels = kmeans.fit_predict(X)
   ari = adjusted_rand_score(
       reference_clusters,
       kmeans_labels
   )
   nmi = normalized_mutual_info_score(
       reference_clusters,
       kmeans_labels
   )
   results.append({
       "k": k,
       "ARI": ari,
       "NMI": nmi
   })
   print(
       f"  ARI = {ari:.4f}"
   )
   print(
       f"  NMI = {nmi:.4f}"
   )
# ############################################################
# 6. SAVE RESULTS
# Create the output folder and write ARI/NMI scores.
# ############################################################
results_df = pd.DataFrame(results)
os.makedirs(
   os.path.dirname(OUTPUT_CSV),
   exist_ok=True
)
results_df.to_csv(
   OUTPUT_CSV,
   index=False
)
# ############################################################
# 7. FINISH
# Print the validation results and output location.
# ############################################################
print("\n========================================")
print("K-MEANS VALIDATION COMPLETE")
print("========================================")
print("\nResults:")
print(results_df.to_string(index=False))
print("\nSaved:")
print(OUTPUT_CSV)
