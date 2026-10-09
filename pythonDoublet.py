# ############################################################
# 1. IMPORTS
# Load the packages used to read, process, and save the data.
# ############################################################
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
import scrublet as scr
import doubletdetection

# ############################################################
# 2. PATHS AND SAMPLE CONFIGURATION
# Set the sample-specific input and output paths.
# ############################################################
SAMPLE_DIR = os.getcwd()
SAMPLE = os.path.basename(SAMPLE_DIR)
INPUT_H5 = os.path.join(
    "/home/cv1u24/GIST/1_Raw_Processing/CellBender",
    SAMPLE,
    "CellBender_Outputs",
    f"{SAMPLE}_cellbender_filtered.h5"
)
OUTPUT_DIR = os.path.join(SAMPLE_DIR, "Output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
if not os.path.exists(INPUT_H5):
    raise FileNotFoundError(f"Input file not found: {INPUT_H5}")
print(f"Processing sample: {SAMPLE}")
print(f"Input: {INPUT_H5}")
print(f"Output directory: {OUTPUT_DIR}")

# ############################################################
# 3. LOAD THE CELLBENDER-FILTERED MATRIX
# Read the input matrix and ensure it is in sparse CSR format.
# ############################################################
adata = sc.read_10x_h5(INPUT_H5)
adata.var_names_make_unique()
X = adata.X
if not sparse.issparse(X):
    X = sparse.csr_matrix(X)
else:
    X = X.tocsr()
print(f"Cells: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")

# ############################################################
# 4. RUN SCRUBLET
# Detect doublets using Scrublet's default settings.
# ############################################################
print("\nRunning Scrublet...")
scrub = scr.Scrublet(X)
scrub_scores, scrub_calls = scrub.scrub_doublets()
scrub_scores = np.asarray(scrub_scores)
scrub_calls = np.asarray(scrub_calls).astype(bool)
print(f"Scrublet doublets detected: {scrub_calls.sum()}")

# ############################################################
# 5. RUN DOUBLETDETECTION
# Detect doublets using DoubletDetection's default settings.
# ############################################################
print("\nRunning DoubletDetection...")
clf = doubletdetection.BoostClassifier()
dd_calls = clf.fit(X).predict()
dd_calls = np.asarray(dd_calls).astype(bool)
try:
    dd_scores = np.asarray(clf.doublet_score())
except Exception:
    print("Warning: could not extract DoubletDetection scores. Saving NaN values.")
    dd_scores = np.full(adata.n_obs, np.nan)
print(f"DoubletDetection doublets detected: {dd_calls.sum()}")

# ############################################################
# 6. SAVE PER-BARCODE CALLS
# Write both methods' scores and calls for each barcode to CSV.
# ############################################################
barcode_df = pd.DataFrame({
    "sample": SAMPLE,
    "barcode": adata.obs_names,
    "scrublet_score": scrub_scores,
    "scrublet_doublet": scrub_calls,
    "doubletdetection_score": dd_scores,
    "doubletdetection_doublet": dd_calls
})
barcode_output = os.path.join(
    OUTPUT_DIR,
    f"{SAMPLE}_doublet_barcode_calls.csv"
)
barcode_df.to_csv(barcode_output, index=False)

# ############################################################
# 7. SAVE SUMMARY STATISTICS
# Summarise the number and scores of detected doublets.
# ############################################################
summary_df = pd.DataFrame([{
    "sample": SAMPLE,
    "n_cells": adata.n_obs,
    "n_genes": adata.n_vars,
    "scrublet_n_doublets": int(scrub_calls.sum()),
    "scrublet_pct_doublets": float(scrub_calls.mean() * 100),
    "scrublet_mean_score": float(np.nanmean(scrub_scores)),
    "scrublet_median_score": float(np.nanmedian(scrub_scores)),
    "scrublet_max_score": float(np.nanmax(scrub_scores)),
    "doubletdetection_n_doublets": int(dd_calls.sum()),
    "doubletdetection_pct_doublets": float(dd_calls.mean() * 100),
    "doubletdetection_mean_score": float(np.nanmean(dd_scores)),
    "doubletdetection_median_score": float(np.nanmedian(dd_scores)),
    "doubletdetection_max_score": float(np.nanmax(dd_scores))
}])
summary_output = os.path.join(
    OUTPUT_DIR,
    f"{SAMPLE}_doublet_summary.csv"
)
summary_df.to_csv(summary_output, index=False)

# ############################################################
# 8. SAVE THE ANNOTATED H5AD OBJECT
# Add the doublet calls and scores to the AnnData object and save it.
# ############################################################
adata.obs["scrublet_score"] = scrub_scores
adata.obs["scrublet_doublet"] = scrub_calls
adata.obs["doubletdetection_score"] = dd_scores
adata.obs["doubletdetection_doublet"] = dd_calls
h5ad_output = os.path.join(
    OUTPUT_DIR,
    f"{SAMPLE}_doublet_calls.h5ad"
)
adata.write_h5ad(h5ad_output)

print("\nSaved files:")
print(barcode_output)
print(summary_output)
print(h5ad_output)
print("\nDone.")
