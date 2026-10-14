#!/usr/bin/env python3
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import silhouette_score, adjusted_rand_score
# ############################################################
# 1. IMPORTS AND CONFIGURATION
# Load analysis packages and define the input, output, and diagnostic settings.
# ############################################################
BASE_DIR = "/home/cv1u24/GIST/4_Preprocessing_and_Integration"
H5AD_DIR = os.path.join(BASE_DIR, "Outputs", "h5ad")
FIG_DIR = os.path.join(BASE_DIR, "Outputs", "Figures", "PCA_Diagnostics")
CSV_DIR = os.path.join(BASE_DIR, "Outputs", "CSVs", "PCA_Diagnostics")
os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(CSV_DIR, exist_ok=True)
INPUT = os.path.join(
    H5AD_DIR,
    "merged_singlets_normalised_log_hvg_scaled_pca.h5ad"
)
PC_OPTIONS = [10, 15, 20, 25, 30]
N_NEIGHBORS = 15
LEIDEN_RESOLUTION = 0.5
RANDOM_STATE = 42

# ############################################################
# 2. LOAD INPUT DATA
# Read the preprocessed AnnData object.
# ############################################################
print("\n" + "=" * 70)
print("PCA DIAGNOSTICS AND PC SELECTION")
print("=" * 70)
adata = sc.read_h5ad(INPUT)
print(f"\nLoaded: {INPUT}")
print(adata)
print(f"Nuclei: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")

# ############################################################
# 3. CHECK PCA RESULTS
# Confirm PCA fields and determine the available tested component counts.
# ############################################################
if "X_pca" not in adata.obsm:
    raise ValueError("X_pca not found in AnnData.")
if "pca" not in adata.uns:
    raise ValueError("PCA information not found in adata.uns.")
X_pca = adata.obsm["X_pca"]
variance = np.asarray(
    adata.uns["pca"]["variance_ratio"]
)
n_pcs_available = X_pca.shape[1]
PC_OPTIONS = [
    pc for pc in PC_OPTIONS
    if pc <= n_pcs_available
]
print(f"\nPCs available: {n_pcs_available}")
print(f"PC numbers being tested: {PC_OPTIONS}")

# ############################################################
# 4. CHECK SAMPLE METADATA
# Inspect sample and subtype annotations used in diagnostics.
# ############################################################
if "sample" not in adata.obs:
    raise ValueError("sample metadata not found.")
print("\nSample counts:")
print(
    adata.obs["sample"]
    .value_counts()
    .sort_index()
)
if "gist_subtype" in adata.obs:
    print("\nGIST subtype counts:")
    print(
        adata.obs["gist_subtype"]
        .value_counts()
    )
else:
    print("\nWARNING: gist_subtype not found.")

# ############################################################
# 5. SUMMARISE PCA VARIANCE
# Calculate per-component and cumulative variance and export the table.
# ############################################################
print("\n" + "=" * 70)
print("PCA VARIANCE")
print("=" * 70)
cumulative = np.cumsum(variance)
variance_df = pd.DataFrame({
    "PC": np.arange(1, len(variance) + 1),
    "variance_ratio": variance,
    "variance_percent": variance * 100,
    "cumulative_variance": cumulative,
    "cumulative_percent": cumulative * 100
})
variance_df.to_csv(
    os.path.join(
        CSV_DIR,
        "PCA_variance_table.csv"
    ),
    index=False
)
print(
    variance_df.head(20).to_string(index=False)
)

# ############################################################
# 6. PLOT PCA VARIANCE
# Save scree and cumulative-variance figures.
# ############################################################
print("\nCreating PCA variance figures...")
# Scree plot
plt.figure(figsize=(9, 6))
plt.plot(
    variance_df["PC"],
    variance_df["variance_percent"],
    marker="o",
    markersize=3
)
plt.xlabel("Principal Component")
plt.ylabel("Variance explained (%)")
plt.title("PCA Scree Plot")
plt.tight_layout()
plt.savefig(
    os.path.join(
        FIG_DIR,
        "PCA_scree_plot.png"
    ),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Cumulative variance
plt.figure(figsize=(9, 6))
plt.plot(
    variance_df["PC"],
    variance_df["cumulative_percent"],
    marker="o",
    markersize=3
)
plt.xlabel("Number of PCs")
plt.ylabel("Cumulative variance explained (%)")
plt.title("Cumulative PCA Variance")
plt.tight_layout()
plt.savefig(
    os.path.join(
        FIG_DIR,
        "PCA_cumulative_variance.png"
    ),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# ############################################################
# 7. SUMMARISE TESTED PC COUNTS
# Report variance captured at each tested number of PCs.
# ############################################################
pc_summary = []
for pc in PC_OPTIONS:
    pc_summary.append({
        "PCs": pc,
        "variance_explained_percent":
            variance_df.loc[
                variance_df["PC"] == pc,
                "variance_percent"
            ].iloc[0],
        "cumulative_variance_percent":
            variance_df.loc[
                variance_df["PC"] == pc,
                "cumulative_percent"
            ].iloc[0]
    })
pc_summary_df = pd.DataFrame(pc_summary)
pc_summary_df.to_csv(
    os.path.join(
        CSV_DIR,
        "PC_number_summary.csv"
    ),
    index=False
)
print("\nVariance at tested PC numbers:")
print(
    pc_summary_df.to_string(index=False)
)

# ############################################################
# 8. CREATE PCA SCATTER PLOTS
# Plot selected PC pairs coloured by sample and other metadata.
# ############################################################
print("\nCreating PCA scatter plots...")
scatter_pairs = [
    (1, 2),
    (1, 3),
    (2, 3)
]
for pc_x, pc_y in scatter_pairs:
    if pc_x > n_pcs_available or pc_y > n_pcs_available:
        continue
    df = pd.DataFrame({
        "PCx": X_pca[:, pc_x - 1],
        "PCy": X_pca[:, pc_y - 1],
        "sample": adata.obs["sample"].astype(str).values
    })
    plt.figure(figsize=(9, 7))
    sns.scatterplot(
        data=df,
        x="PCx",
        y="PCy",
        hue="sample",
        s=10,
        alpha=0.65,
        linewidth=0
    )
    plt.xlabel(
        f"PC{pc_x} "
        f"({variance[pc_x - 1] * 100:.2f}% variance)"
    )
    plt.ylabel(
        f"PC{pc_y} "
        f"({variance[pc_y - 1] * 100:.2f}% variance)"
    )
    plt.title(
        f"PCA: PC{pc_x} vs PC{pc_y} by Sample"
    )
    plt.legend(
        title="Sample",
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )
    plt.tight_layout()
    plt.savefig(
        os.path.join(
            FIG_DIR,
            f"PCA_PC{pc_x}_PC{pc_y}_by_sample.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

# PCA by subtype
if "gist_subtype" in adata.obs:
    df = pd.DataFrame({
        "PC1": X_pca[:, 0],
        "PC2": X_pca[:, 1],
        "gist_subtype":
            adata.obs["gist_subtype"].astype(str).values
    })
    plt.figure(figsize=(9, 7))
    sns.scatterplot(
        data=df,
        x="PC1",
        y="PC2",
        hue="gist_subtype",
        s=12,
        alpha=0.65,
        linewidth=0
    )
    plt.xlabel(
        f"PC1 ({variance[0] * 100:.2f}% variance)"
    )
    plt.ylabel(
        f"PC2 ({variance[1] * 100:.2f}% variance)"
    )
    plt.title("PCA: PC1 vs PC2 by GIST Subtype")
    plt.legend(
        title="GIST subtype",
        bbox_to_anchor=(1.02, 1),
        loc="upper left"
    )
    plt.tight_layout()
    plt.savefig(
        os.path.join(
            FIG_DIR,
            "PCA_PC1_PC2_by_gist_subtype.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

# ############################################################
# 9. TEST PC–QC CORRELATIONS
# Evaluate associations between principal-component scores and QC metrics.
# ############################################################
print("\n" + "=" * 70)
print("PC–QC CORRELATIONS")
print("=" * 70)
QC_VARIABLES = [
    "n_genes_by_counts",
    "total_counts",
    "total_counts_mt",
    "pct_counts_mt"
]
QC_VARIABLES = [
    x for x in QC_VARIABLES
    if x in adata.obs
]
correlation_rows = []
for pc_number in range(
    1,
    min(20, n_pcs_available) + 1
):
    print(f"\nPC{pc_number}")
    pc_values = X_pca[:, pc_number - 1]
    for qc in QC_VARIABLES:
        qc_values = pd.to_numeric(
            adata.obs[qc],
            errors="coerce"
        ).to_numpy()
        valid = (
            np.isfinite(pc_values)
&
            np.isfinite(qc_values)
        )
        if valid.sum() >= 3:
            correlation = np.corrcoef(
                pc_values[valid],
                qc_values[valid]
            )[0, 1]
        else:
            correlation = np.nan
        print(
            f"  {qc}: {correlation:.4f}"
        )
        correlation_rows.append({
            "PC": pc_number,
            "QC_metric": qc,
            "correlation": correlation,
            "n_valid": valid.sum()
        })
pd.DataFrame(
    correlation_rows
).to_csv(
    os.path.join(
        CSV_DIR,
        "PC_QC_correlations.csv"
    ),
    index=False
)

# ############################################################
# 10. PLOT PC SCORES BY SAMPLE
# Compare principal-component scores across samples.
# ############################################################
print("\nCreating PC score summaries...")
n_summary_pcs = min(20, n_pcs_available)
pc_df = pd.DataFrame(
    X_pca[:, :n_summary_pcs],
    columns=[
        f"PC{i}"
        for i in range(1, n_summary_pcs + 1)
    ],
    index=adata.obs_names
)
pc_df["sample"] = (
    adata.obs["sample"].astype(str).values
)
mean_pc = (
    pc_df
    .groupby("sample")
    .mean()
)
median_pc = (
    pc_df
    .groupby("sample")
    .median()
)
mean_pc.to_csv(
    os.path.join(
        CSV_DIR,
        "mean_PC_scores_by_sample.csv"
    )
)
median_pc.to_csv(
    os.path.join(
        CSV_DIR,
        "median_PC_scores_by_sample.csv"
    )
)

# ############################################################
# 11. SUMMARISE SAMPLE COUNTS
# Export the number of nuclei contributed by each sample.
# ############################################################
adata.obs["sample"].value_counts().sort_index().to_csv(
    os.path.join(
        CSV_DIR,
        "nuclei_per_sample.csv"
    )
)

# ############################################################
# 12. CLUSTER USING TESTED PC COUNTS
# Build neighbour graphs and Leiden clusters for each PC setting; calculate composition and diagnostics.
# ############################################################
print("\n" + "=" * 70)
print("PC SELECTION ANALYSIS")
print("=" * 70)
cluster_results = {}
cluster_summary = []
silhouette_results = []
purity_results = []
for n_pcs in PC_OPTIONS:
    print(
        f"\nRunning Leiden using {n_pcs} PCs..."
    )
    neighbors_key = f"neighbors_{n_pcs}PCs"
    cluster_key = f"leiden_{n_pcs}PCs"
    sc.pp.neighbors(
        adata,
        n_neighbors=N_NEIGHBORS,
        n_pcs=n_pcs,
        use_rep="X_pca",
        key_added=neighbors_key
    )
    sc.tl.leiden(
        adata,
        resolution=LEIDEN_RESOLUTION,
        neighbors_key=neighbors_key,
        key_added=cluster_key,
        flavor="igraph",
        directed=False,
        random_state=RANDOM_STATE
    )
    labels = adata.obs[
        cluster_key
    ].astype(str)
    cluster_results[n_pcs] = labels.copy()
    counts = labels.value_counts()
    cluster_summary.append({
        "PCs": n_pcs,
        "number_of_clusters": len(counts),
        "smallest_cluster": counts.min(),
        "largest_cluster": counts.max()
    })
    print(
        f"{n_pcs} PCs -> "
        f"{len(counts)} clusters | "
        f"smallest = {counts.min()} | "
        f"largest = {counts.max()}"
    )

    # --------------------------------------------------------
    # CLUSTER × SAMPLE
    # --------------------------------------------------------
    sample_table = pd.crosstab(
        labels,
        adata.obs["sample"]
    )
    sample_table.to_csv(
        os.path.join(
            CSV_DIR,
            f"cluster_sample_{n_pcs}PCs.csv"
        )
    )
    sample_percent = (
        sample_table
        .div(
            sample_table.sum(axis=1),
            axis=0
        )
        * 100
    )
    sample_percent.to_csv(
        os.path.join(
            CSV_DIR,
            f"cluster_sample_percentage_{n_pcs}PCs.csv"
        )
    )

    # --------------------------------------------------------
    # CLUSTER × GIST SUBTYPE
    # --------------------------------------------------------
    if "gist_subtype" in adata.obs:
        subtype_table = pd.crosstab(
            labels,
            adata.obs["gist_subtype"]
        )
        subtype_table.to_csv(
            os.path.join(
                CSV_DIR,
                f"cluster_gist_subtype_{n_pcs}PCs.csv"
            )
        )

    # --------------------------------------------------------
    # SAMPLE PURITY
    # --------------------------------------------------------
    cluster_purity = (
        sample_table.max(axis=1)
        /
        sample_table.sum(axis=1)
    )
    purity_results.append({
        "PCs": n_pcs,
        "mean_cluster_sample_purity":
            cluster_purity.mean(),
        "median_cluster_sample_purity":
            cluster_purity.median(),
        "max_cluster_sample_purity":
            cluster_purity.max()
    })

    # --------------------------------------------------------
    # SILHOUETTE
    # --------------------------------------------------------
    try:
        X = X_pca[:, :n_pcs]
        labels_array = labels.values
        if len(np.unique(labels_array)) > 1:
            max_cells = 10000
            if X.shape[0] > max_cells:
                rng = np.random.default_rng(
                   RANDOM_STATE
                )
                indices = rng.choice(
                   X.shape[0],
                   size=max_cells,
                   replace=False
                )
                X_sil = X[indices]
                labels_sil = labels_array[indices]
            else:
                X_sil = X
                labels_sil = labels_array
            silhouette = silhouette_score(
                X_sil,
                labels_sil
            )
        else:
            silhouette = np.nan
    except Exception as e:
        print(
            f"Silhouette failed for {n_pcs} PCs: {e}"
        )
        silhouette = np.nan
    silhouette_results.append({
        "PCs": n_pcs,
        "silhouette_score": silhouette
    })

# ############################################################
# 13. SAVE PC-SELECTION RESULTS
# Write cluster counts, silhouette scores, and sample association summaries.
# ############################################################
cluster_summary_df = pd.DataFrame(
    cluster_summary
)
silhouette_df = pd.DataFrame(
    silhouette_results
)
purity_df = pd.DataFrame(
    purity_results
)
cluster_summary_df.to_csv(
    os.path.join(
        CSV_DIR,
        "cluster_counts_by_PC_number.csv"
    ),
    index=False
)
silhouette_df.to_csv(
    os.path.join(
        CSV_DIR,
        "silhouette_scores_by_PC_number.csv"
    ),
    index=False
)
purity_df.to_csv(
    os.path.join(
        CSV_DIR,
        "sample_association_summary.csv"
    ),
    index=False
)

# ############################################################
# 14. CALCULATE ARI STABILITY
# Compare adjacent PC clustering solutions using adjusted Rand index.
# ############################################################
print("\nCalculating clustering stability...")
ari_results = []
for i in range(len(PC_OPTIONS) - 1):
    pc_a = PC_OPTIONS[i]
    pc_b = PC_OPTIONS[i + 1]
    ari = adjusted_rand_score(
        cluster_results[pc_a],
        cluster_results[pc_b]
    )
    ari_results.append({
        "comparison":
            f"{pc_a}_vs_{pc_b}",
        "PCs_A": pc_a,
        "PCs_B": pc_b,
        "ARI": ari
    })
    print(
        f"{pc_a} vs {pc_b} PCs: "
        f"ARI = {ari:.4f}"
    )
ari_df = pd.DataFrame(ari_results)
ari_df.to_csv(
    os.path.join(
        CSV_DIR,
        "clustering_stability_ARI.csv"
    ),
    index=False
)

# ############################################################
# 15. PLOT PC-SELECTION METRICS
# Save cluster-count, silhouette, and ARI comparison figures.
# ############################################################
print("\nCreating PC selection figures...")
# Cluster count
plt.figure(figsize=(8, 6))
plt.plot(
    cluster_summary_df["PCs"],
    cluster_summary_df["number_of_clusters"],
    marker="o"
)
plt.xlabel("Number of PCs")
plt.ylabel("Number of Leiden clusters")
plt.title("Number of Leiden Clusters vs Number of PCs")
plt.xticks(PC_OPTIONS)
plt.tight_layout()
plt.savefig(
    os.path.join(
        FIG_DIR,
        "cluster_count_vs_PCs.png"
    ),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# Silhouette
plt.figure(figsize=(8, 6))
plt.plot(
    silhouette_df["PCs"],
    silhouette_df["silhouette_score"],
    marker="o"
)
plt.xlabel("Number of PCs")
plt.ylabel("Silhouette score")
plt.title("Silhouette Score vs Number of PCs")
plt.xticks(PC_OPTIONS)
plt.tight_layout()
plt.savefig(
    os.path.join(
        FIG_DIR,
        "silhouette_vs_PCs.png"
    ),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# ARI
plt.figure(figsize=(8, 6))
plt.plot(
    ari_df["comparison"],
    ari_df["ARI"],
    marker="o"
)
plt.xlabel("PC comparison")
plt.ylabel("Adjusted Rand Index")
plt.title("Clustering Stability Between PC Solutions")
plt.ylim(0, 1.05)
plt.tight_layout()
plt.savefig(
    os.path.join(
        FIG_DIR,
        "clustering_stability_ARI.png"
    ),
    dpi=300,
    bbox_inches="tight"
)
plt.close()

# ############################################################
# 16. EXPORT TOP PCA LOADINGS
# Save positive and negative gene loadings for the first principal components.
# ############################################################
print("\nSaving top PC loadings...")
loadings = adata.varm["PCs"]
rows = []
for pc in range(
    min(20, loadings.shape[1])
):
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
    os.path.join(
        CSV_DIR,
        "top_PC_loadings.csv"
    ),
    index=False
)

# ############################################################
# 17. PRINT FINAL SUMMARY
# Report diagnostic results and output locations.
# ############################################################
print("\n" + "=" * 70)
print("PC DIAGNOSTIC ANALYSIS COMPLETE")
print("=" * 70)
print("\nPCs tested:")
print(PC_OPTIONS)
print("\nCluster counts:")
print(
    cluster_summary_df.to_string(index=False)
)
print("\nSilhouette scores:")
print(
    silhouette_df.to_string(index=False)
)
print("\nARI stability:")
print(
    ari_df.to_string(index=False)
)
print("\nSample association:")
print(
    purity_df.to_string(index=False)
)
print("\nFigures saved to:")
print(FIG_DIR)
print("\nCSVs saved to:")
print(CSV_DIR)
print("\nIMPORTANT:")
print(
    "Do not automatically select the final PC number "
    "from a single metric."
)
print(
    "Consider PCA variance, PC loadings, QC associations, "
    "sample association, clustering stability and "
    "silhouette scores together."
)
print("\nDone.")
