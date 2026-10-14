#!/usr/bin/env python
import os
import scanpy as sc
import scvi
import pandas as pd

# ############################################################
# 1. IMPORTS AND CONFIGURATION
# Load packages and define the input object and scVI comparison settings.
# ############################################################
INPUT = "/iridisfs/home/cv1u24/GIST/4_Preprocessing_and_Integration/Outputs/h5ad/merged_singlets_normalised_log_hvg_scaled_pca.h5ad"
OUT = "/iridisfs/home/cv1u24/GIST/4_Preprocessing_and_Integration/Outputs/BatchToolComparisons/scVI"
FIG_DIR = os.path.join(OUT, "Figures")
H5AD_DIR = os.path.join(OUT, "h5ad")
CSV_DIR = os.path.join(OUT, "CSV")
for d in [FIG_DIR, H5AD_DIR, CSV_DIR]:
    os.makedirs(d, exist_ok=True)
PC_VALUES = [20, 25, 30]

# ############################################################
# 2. LOAD INPUT DATA
# Read the preprocessed AnnData object and report its dimensions.
# ############################################################
adata = sc.read_h5ad(INPUT)
print("=" * 70)
print("scVI INTEGRATION")
print("=" * 70)
print(f"Loaded: {INPUT}")
print(f"Nuclei: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")

# ############################################################
# 3. PREPARE RAW COUNTS
# Confirm a counts layer is available and set the scVI random seed.
# ############################################################
if "counts" not in adata.layers:
    raise ValueError("counts layer not found.")

# scVI needs raw counts
scvi.settings.seed = 0

# ############################################################
# 4. TRAIN AND SUMMARISE scVI MODELS
# Train models for each requested latent dimension, export embeddings and UMAPs, and save the resulting AnnData objects.
# ############################################################
for n_pcs in PC_VALUES:
    print("\n" + "=" * 70)
    print(f"RUNNING scVI - {n_pcs} PCs")
    print("=" * 70)
    result = adata.copy()

    # --------------------------------------------------------
    # scVI uses the count matrix.
    # PCs are used only to determine the number of dimensions
    # retained for downstream comparison.
    # --------------------------------------------------------
    scvi.model.SCVI.setup_anndata(
        result,
        layer="counts",
        batch_key="sample"
    )
    model = scvi.model.SCVI(
        result,
        n_latent=n_pcs
    )
    print(f"Training scVI with n_latent={n_pcs}...")
    model.train()

    # --------------------------------------------------------
    # GET LATENT REPRESENTATION
    # --------------------------------------------------------
    result.obsm["X_scVI"] = model.get_latent_representation()
    print(
        f"scVI latent representation: "
        f"{result.obsm['X_scVI'].shape}"
    )

    # --------------------------------------------------------
    # SAVE LATENT SCORES
    # --------------------------------------------------------
    latent_df = pd.DataFrame(
        result.obsm["X_scVI"],
        index=result.obs_names,
        columns=[f"scVI_{i+1}" for i in range(n_pcs)]
    )
    latent_file = os.path.join(
        CSV_DIR,
        f"scVI_{n_pcs}PCs_latent.csv"
    )
    latent_df.to_csv(latent_file)

    # --------------------------------------------------------
    # NEIGHBOURS + UMAP
    # --------------------------------------------------------
    sc.pp.neighbors(
        result,
        use_rep="X_scVI",
        n_neighbors=15
    )
    sc.tl.umap(result)

    # --------------------------------------------------------
    # UMAP - SAMPLE
    # --------------------------------------------------------
    sc.pl.umap(
        result,
        color="sample",
        show=False,
        title=f"scVI - {n_pcs} PCs - Sample"
    )
    fig = sc.pl.umap(
        result,
        color="sample",
        show=False,
        title=f"scVI - {n_pcs} PCs - Sample",
        return_fig=True
    )
    fig.savefig(
        os.path.join(
            FIG_DIR,
            f"scVI_{n_pcs}PCs_sample.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )

    # --------------------------------------------------------
    # UMAP - GIST SUBTYPE
    # --------------------------------------------------------
    if "gist_subtype" in result.obs:
        fig = sc.pl.umap(
            result,
            color="gist_subtype",
            show=False,
            title=f"scVI - {n_pcs} PCs - GIST subtype",
            return_fig=True
        )
        fig.savefig(
            os.path.join(
                FIG_DIR,
                f"scVI_{n_pcs}PCs_subtype.png"
            ),
            dpi=300,
            bbox_inches="tight"
        )

    # --------------------------------------------------------
    # SAVE H5AD
    # --------------------------------------------------------
    h5ad_file = os.path.join(
        H5AD_DIR,
        f"merged_scVI_{n_pcs}PCs.h5ad"
    )
    result.write(h5ad_file)
    print(f"Saved H5AD: {h5ad_file}")
    print(f"Saved latent CSV: {latent_file}")
print("\n" + "=" * 70)
print("scVI COMPLETE")
print("=" * 70)
print(f"Outputs: {OUT}")
print(f"Figures: {FIG_DIR}")
print(f"H5AD:    {H5AD_DIR}")
print(f"CSV:     {CSV_DIR}")
