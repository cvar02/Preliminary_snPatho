#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
# ############################################################
# 1. IMPORTS AND CONFIGURATION
# Load analysis packages and set input, output, and preprocessing parameters.
# ############################################################
BASE_DIR = "/home/cv1u24/GIST/4_Preprocessing_and_Integration"
H5AD_DIR = os.path.join(BASE_DIR, "Outputs", "h5ad")
FIG_DIR = os.path.join(BASE_DIR, "Outputs", "Figures", "PCA")
CSV_DIR = os.path.join(BASE_DIR, "Outputs", "CSVs", "PCA")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(CSV_DIR, exist_ok=True)
INPUT = os.path.join(
    H5AD_DIR,
    "merged_singlets_unintegrated.h5ad"
)
OUTPUT = os.path.join(
    H5AD_DIR,
    "merged_singlets_normalised_log_hvg_scaled_pca.h5ad"
)
TARGET_SUM = 1e4
N_HVGS = 3000
N_PCS = 50
SCALE_MAX = 10

# ############################################################
# 2. LOAD INPUT DATA
# Read the merged singlet AnnData object and report its dimensions.
# ############################################################
print("\n" + "=" * 70)
print("PREPROCESSING AND PCA")
print("=" * 70)
adata = sc.read_h5ad(INPUT)
print(f"\nLoaded: {INPUT}")
print(adata)
print(f"Nuclei: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")

# ############################################################
# 3. CHECK REQUIRED METADATA
# Confirm that sample and GIST subtype annotations are present.
# ############################################################
for column in ["sample", "gist_subtype"]:
    if column not in adata.obs:
        raise ValueError(f"Missing required metadata: {column}")
print("\nSamples:")
print(adata.obs["sample"].value_counts().sort_index())
print("\nGIST subtypes:")
print(adata.obs["gist_subtype"].value_counts())

# ############################################################
# 4. PRESERVE RAW COUNTS
# Store or restore the raw count matrix in the counts layer.
# ############################################################
if "counts" not in adata.layers:
    print("\nCreating counts layer from X...")
    adata.layers["counts"] = adata.X.copy()
else:
    print("\nRaw counts layer found.")
adata.X = adata.layers["counts"].copy()

# ############################################################
# 5. NORMALISE AND LOG-TRANSFORM
# Apply total-count normalisation and log1p transformation.
# ############################################################
print("\nNormalising counts...")
sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
print("Log transforming...")
sc.pp.log1p(adata)
adata.layers["lognorm"] = adata.X.copy()

# ############################################################
# 6. SELECT HIGHLY VARIABLE GENES
# Identify the requested number of highly variable genes and save their summary.
# ############################################################
print(f"\nSelecting {N_HVGS} highly variable genes...")
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=N_HVGS,
    subset=False
)
n_hvgs = int(adata.var["highly_variable"].sum())
print(f"HVGs selected: {n_hvgs}")
hvg_cols = [
    c for c in [
        "highly_variable",
        "means",
        "dispersions",
        "dispersions_norm"
    ]
    if c in adata.var
]
adata.var[hvg_cols].sort_values(
    "highly_variable",
    ascending=False
).to_csv(
    os.path.join(CSV_DIR, "highly_variable_genes_summary.csv")
)

# ############################################################
# 7. SCALE AND RUN PCA
# Scale the expression matrix and calculate principal components.
# ############################################################
print(f"\nScaling data (max_value={SCALE_MAX})...")
sc.pp.scale(adata, max_value=SCALE_MAX)
print(f"Running PCA using {N_HVGS} HVGs...")
sc.tl.pca(
    adata,
    n_comps=N_PCS,
    use_highly_variable=True,
    svd_solver="arpack"
)

# ############################################################
# 8. SUMMARISE PCA VARIANCE
# Calculate and save per-component and cumulative variance.
# ############################################################
variance = adata.uns["pca"]["variance_ratio"]
cumulative = np.cumsum(variance)
variance_df = pd.DataFrame({
    "PC": np.arange(1, len(variance) + 1),
    "variance_ratio": variance,
    "variance_percent": variance * 100,
    "cumulative_variance": cumulative,
    "cumulative_percent": cumulative * 100
})
variance_df.to_csv(
    os.path.join(CSV_DIR, "pca_variance_summary.csv"),
    index=False
)
print("\nPCA variance:")
for i in range(len(variance)):
    print(
        f"PC{i + 1}: "
        f"{variance[i]:.6f} "
        f"(cumulative {cumulative[i]:.6f})"
    )

# ############################################################
# 9. CREATE PCA FIGURES
# Define the plotting helper and generate PCA plots coloured by sample, subtype, genotype, and QC metrics.
# ############################################################
def save_pca_plot(color, filename):
    if color not in adata.obs:
        print(f"Skipping {color}: not found.")
        return
    sc.pl.pca(
        adata,
        color=color,
        components=["1,2"],
        show=False
    )
    plt.savefig(
        os.path.join(FIG_DIR, filename),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()
print("\nCreating PCA figures...")
# Always use sample
save_pca_plot(
    "sample",
    "PCA_PC1_PC2_by_sample.png"
)
# Always use GIST subtype
save_pca_plot(
    "gist_subtype",
    "PCA_PC1_PC2_by_gist_subtype.png"
)
# Genotype only if it actually exists
if "genotype_group" in adata.obs.columns:
    save_pca_plot(
        "genotype_group",
        "PCA_PC1_PC2_by_genotype.png"
    )
else:
    print("genotype_group not present - using gist_subtype instead.")
# QC metrics
for qc in [
    "total_counts",
    "n_genes_by_counts",
    "total_counts_mt",
    "pct_counts_mt"
]:
    save_pca_plot(
        qc,
        f"PCA_PC1_PC2_by_{qc}.png"
    )
# ############################################################
# 10. PLOT PCA VARIANCE
# Save the PCA variance visualisation.
# ############################################################
sc.pl.pca_variance_ratio(
    adata,
    n_pcs=N_PCS,
    log=True,
    show=False
)
plt.savefig(
    os.path.join(FIG_DIR, "PCA_variance_elbow.png"),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# ############################################################
# 11. EXPORT TOP PCA LOADINGS
# Write the highest-ranked gene loadings for the principal components.
# ############################################################
print("\nSaving top PCA genes...")
loadings = adata.varm["PCs"]
rows = []
for pc in range(min(20, loadings.shape[1])):
    values = loadings[:, pc]
    positive = np.argsort(values)[-20:][::-1]
    negative = np.argsort(values)[:20]
    for i in positive:
        rows.append({
           "PC": pc + 1,
           "direction": "positive",
           "gene": adata.var_names[i],
           "loading": values[i]
        })
    for i in negative:
        rows.append({
           "PC": pc + 1,
           "direction": "negative",
           "gene": adata.var_names[i],
           "loading": values[i]
        })
pd.DataFrame(rows).to_csv(
    os.path.join(CSV_DIR, "top_PC_loadings.csv"),
    index=False
)

# ############################################################
# 12. SAVE PCA OBJECT
# Save the processed AnnData object containing PCA results.
# ############################################################
adata.write_h5ad(OUTPUT)

# ############################################################
# 13. PRINT SUMMARY
# Report the completed preprocessing outputs.
# ############################################################
print("\n" + "=" * 70)
print("PREPROCESSING COMPLETE")
print("=" * 70)
print(f"\nNuclei: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")
print(f"HVGs: {n_hvgs}")
print(f"PCs: {N_PCS}")
print(f"\nProcessed object:")
print(OUTPUT)
print(f"\nFigures:")
print(FIG_DIR)
print(f"\nCSVs:")
print(CSV_DIR)
print("\nDone.")
