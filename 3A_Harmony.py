#!/usr/bin/env python3
import os
import scanpy as sc
import pandas as pd
import numpy as np
import harmonypy as hm

# ############################################################
# 1. IMPORTS AND CONFIGURATION
# Load packages and set the input, output, and PC comparison settings.
# ############################################################
INPUT = "/iridisfs/home/cv1u24/GIST/4_Preprocessing_and_Integration/Outputs/h5ad/merged_singlets_normalised_log_hvg_scaled_pca.h5ad"
OUTDIR = "/iridisfs/home/cv1u24/GIST/4_Preprocessing_and_Integration/Outputs/BatchToolComparisons/Harmony"
PC_VALUES = [20, 25, 30]
os.makedirs(OUTDIR, exist_ok=True)
print("=" * 70)
print("HARMONY INTEGRATION - PC COMPARISON")
print("=" * 70)

# ############################################################
# 2. LOAD AND CHECK INPUT DATA
# Read the PCA object and confirm its required representation and sample metadata.
# ############################################################
print("\nLoading:")
print(INPUT)
adata = sc.read_h5ad(INPUT)
print(adata)
print(f"Nuclei: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")
if "X_pca" not in adata.obsm:
    raise ValueError("X_pca not found.")
if "sample" not in adata.obs:
    raise ValueError("'sample' column not found.")
print(f"\nAvailable PCs: {adata.obsm['X_pca'].shape[1]}")
print(f"Testing PCs: {PC_VALUES}")

# ############################################################
# 3. RUN HARMONY FOR EACH PC SET
# Correct the selected PCA embedding for sample and store each result.
# ############################################################
for n_pcs in PC_VALUES:
    print("\n" + "=" * 70)
    print(f"HARMONY USING {n_pcs} PCs")
    print("=" * 70)

    # --------------------------------------------------------
    # Select ONLY requested PCs
    # --------------------------------------------------------
    X_pca = adata.obsm["X_pca"][:, :n_pcs]
    print(f"Input PCA shape: {X_pca.shape}")

    # --------------------------------------------------------
    # Run Harmony
    # --------------------------------------------------------
    harmony_output = hm.run_harmony(
        X_pca,
        adata.obs,
        "sample",
        max_iter_harmony=20,
        verbose=True
    )

    # --------------------------------------------------------
    # Get corrected embedding
    # --------------------------------------------------------
    X_harmony = harmony_output.Z_corr
    print(f"Raw Harmony output shape: {X_harmony.shape}")

    # --------------------------------------------------------
    # Safety check for orientation
    # --------------------------------------------------------
    if X_harmony.shape[0] != adata.n_obs:
        if X_harmony.shape[1] == adata.n_obs:
            X_harmony = X_harmony.T
            print("Transposed Harmony embedding.")
        else:
            raise ValueError(
                f"Unexpected Harmony shape: {X_harmony.shape}"
            )
    print(f"Final Harmony embedding shape: {X_harmony.shape}")

    # --------------------------------------------------------
    # Create result object
    # --------------------------------------------------------
    result = adata.copy()
    result.obsm["X_harmony"] = X_harmony
    result.uns["harmony"] = {
        "n_pcs": n_pcs,
        "batch_key": "sample",
        "max_iter_harmony": 20
    }

    # ############################################################
    # 4. BUILD NEIGHBOUR GRAPHS AND UMAPS
    # Calculate neighbours and a UMAP from each corrected Harmony embedding.
    # ############################################################
    print("\nCalculating neighbours...")
    sc.pp.neighbors(
        result,
        use_rep="X_harmony",
        n_neighbors=15
    )
    print("Calculating UMAP...")
    sc.tl.umap(
        result,
        random_state=42
    )

    # ############################################################
    # 5. SAVE HARMONY OBJECT
    # Write each Harmony-integrated AnnData object.
    # ############################################################
    outfile = os.path.join(
        OUTDIR,
        f"merged_Harmony_PC{n_pcs}.h5ad"
    )
    result.write(outfile)
    print(f"\nSaved:")
    print(outfile)

    # ############################################################
    # 6. PLOT UMAP BY SAMPLE
    # Save a sample-coloured UMAP for each tested PC setting.
    # ############################################################
    fig = sc.pl.umap(
        result,
        color="sample",
        show=False,
        return_fig=True,
        title=f"Harmony - PC{n_pcs} - Sample"
    )
    fig.savefig(
        os.path.join(
            OUTDIR,
            f"Harmony_PC{n_pcs}_Sample_UMAP.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )

    # ############################################################
    # 7. PLOT UMAP BY SUBTYPE
    # Save a subtype-coloured UMAP when subtype metadata is available.
    # ############################################################
    if "gist_subtype" in result.obs:
        fig = sc.pl.umap(
            result,
            color="gist_subtype",
            show=False,
            return_fig=True,
            title=f"Harmony - PC{n_pcs} - GIST Subtype"
        )
        fig.savefig(
            os.path.join(
                OUTDIR,
                f"Harmony_PC{n_pcs}_Subtype_UMAP.png"
            ),
            dpi=300,
            bbox_inches="tight"
        )

    # ############################################################
    # 8. EXPORT UMAP COORDINATES
    # Write UMAP coordinates and sample/subtype labels to CSV.
    # ############################################################
    umap_df = pd.DataFrame(
        result.obsm["X_umap"],
        index=result.obs_names,
        columns=["UMAP1", "UMAP2"]
    )
    umap_df["sample"] = result.obs["sample"].values
    if "gist_subtype" in result.obs:
        umap_df["gist_subtype"] = result.obs["gist_subtype"].values
    umap_df.to_csv(
        os.path.join(
            OUTDIR,
            f"Harmony_PC{n_pcs}_UMAP_coordinates.csv"
        )
    )
    print(f"Finished Harmony PC{n_pcs}")

# ############################################################
# 9. PRINT COMPLETION SUMMARY
# Report the PC settings tested and the output location.
# ############################################################
print("\n" + "=" * 70)
print("HARMONY PC COMPARISON COMPLETE")
print("=" * 70)
print("\nTested:")
for n in PC_VALUES:
    print(f"  Harmony PC{n}")
print("\nOutputs saved to:")
print(OUTDIR)

#!/usr/bin/env python3
import os
import scanpy as sc
import pandas as pd
import harmonypy as hm

# ############################################################
# 1. IMPORTS AND CONFIGURATION
# Load packages and define the input, output folders, and PC comparison settings.
# ############################################################
INPUT = "/iridisfs/home/cv1u24/GIST/4_Preprocessing_and_Integration/Outputs/h5ad/merged_singlets_normalised_log_hvg_scaled_pca.h5ad"
OUTDIR = "/iridisfs/home/cv1u24/GIST/4_Preprocessing_and_Integration/Outputs/BatchToolComparisons/Harmony"
PC_VALUES = [20, 25, 30]

# ############################################################
# 2. PREPARE OUTPUT DIRECTORIES
# Create figure, h5ad, and CSV output folders.
# ############################################################
FIGDIR = os.path.join(OUTDIR, "Figures")
H5ADDIR = os.path.join(OUTDIR, "h5ad")
CSVDIR = os.path.join(OUTDIR, "CSV")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(H5ADDIR, exist_ok=True)
os.makedirs(CSVDIR, exist_ok=True)

# ############################################################
# 3. LOAD AND CHECK INPUT DATA
# Read the PCA object and confirm its required representation and sample metadata.
# ############################################################
print("=" * 70)
print("HARMONY INTEGRATION - PC COMPARISON")
print("=" * 70)
print("\nLoading:")
print(INPUT)
adata = sc.read_h5ad(INPUT)
print(adata)
print(f"Nuclei: {adata.n_obs}")
print(f"Genes: {adata.n_vars}")
if "X_pca" not in adata.obsm:
    raise ValueError("X_pca not found.")
if "sample" not in adata.obs:
    raise ValueError("'sample' column not found.")
print(f"\nAvailable PCs: {adata.obsm['X_pca'].shape[1]}")
print(f"Testing PCs: {PC_VALUES}")

# ############################################################
# 4. RUN HARMONY FOR EACH PC SET
# Correct the selected PCA embedding for sample and store each result.
# ############################################################
for n_pcs in PC_VALUES:
    print("\n" + "=" * 70)
    print(f"HARMONY USING {n_pcs} PCs")
    print("=" * 70)

    # --------------------------------------------------------
    # Select ONLY requested PCs
    # --------------------------------------------------------
    X_pca = adata.obsm["X_pca"][:, :n_pcs]
    print(f"Input PCA shape: {X_pca.shape}")

    # --------------------------------------------------------
    # Run Harmony
    # --------------------------------------------------------
    harmony_output = hm.run_harmony(
        X_pca,
        adata.obs,
        "sample",
        max_iter_harmony=20,
        verbose=True
    )

    # --------------------------------------------------------
    # Get corrected embedding
    # --------------------------------------------------------
    X_harmony = harmony_output.Z_corr
    print(f"Raw Harmony output shape: {X_harmony.shape}")

    # Safety check
    if X_harmony.shape[0] != adata.n_obs:
        if X_harmony.shape[1] == adata.n_obs:
            X_harmony = X_harmony.T
            print("Transposed Harmony embedding.")
        else:
            raise ValueError(
                f"Unexpected Harmony shape: {X_harmony.shape}"
            )
    print(f"Final Harmony embedding shape: {X_harmony.shape}")

    # --------------------------------------------------------
    # Create result object
    # --------------------------------------------------------
    result = adata.copy()
    result.obsm["X_harmony"] = X_harmony
    result.uns["harmony"] = {
        "n_pcs": n_pcs,
        "batch_key": "sample",
        "max_iter_harmony": 20
    }

    # ############################################################
    # 5. BUILD NEIGHBOUR GRAPHS AND UMAPS
    # Calculate neighbours and a UMAP from each corrected Harmony embedding.
    # ############################################################
    print("\nCalculating neighbours...")
    sc.pp.neighbors(
        result,
        use_rep="X_harmony",
        n_neighbors=15
    )
    print("Calculating UMAP...")
    sc.tl.umap(
        result,
        random_state=42
    )

    # ############################################################
    # 6. SAVE HARMONY OBJECT
    # Write each Harmony-integrated AnnData object.
    # ############################################################
    h5ad_file = os.path.join(
        H5ADDIR,
        f"merged_Harmony_PC{n_pcs}.h5ad"
    )
    result.write(h5ad_file)
    print(f"\nSaved h5ad:")
    print(h5ad_file)

    # ############################################################
    # 7. PLOT UMAP BY SAMPLE
    # Save a sample-coloured UMAP for each tested PC setting.
    # ############################################################
    fig = sc.pl.umap(
        result,
        color="sample",
        show=False,
        return_fig=True,
        title=f"Harmony - PC{n_pcs} - Sample"
    )
    fig.savefig(
        os.path.join(
            FIGDIR,
            f"Harmony_PC{n_pcs}_Sample_UMAP.png"
        ),
        dpi=300,
        bbox_inches="tight"
    )

    # ############################################################
    # 8. PLOT UMAP BY SUBTYPE
    # Save a subtype-coloured UMAP when subtype metadata is available.
    # ############################################################
    if "gist_subtype" in result.obs:
        fig = sc.pl.umap(
            result,
            color="gist_subtype",
            show=False,
            return_fig=True,
            title=f"Harmony - PC{n_pcs} - GIST Subtype"
        )
        fig.savefig(
            os.path.join(
                FIGDIR,
                f"Harmony_PC{n_pcs}_Subtype_UMAP.png"
            ),
            dpi=300,
            bbox_inches="tight"
        )

    # ############################################################
    # 9. EXPORT UMAP COORDINATES
    # Write UMAP coordinates and sample/subtype labels to CSV.
    # ############################################################
    umap_df = pd.DataFrame(
        result.obsm["X_umap"],
        index=result.obs_names,
        columns=["UMAP1", "UMAP2"]
    )
    umap_df["sample"] = result.obs["sample"].values
    if "gist_subtype" in result.obs:
        umap_df["gist_subtype"] = result.obs["gist_subtype"].values
    csv_file = os.path.join(
        CSVDIR,
        f"Harmony_PC{n_pcs}_UMAP_coordinates.csv"
    )
    umap_df.to_csv(csv_file)
    print(f"Saved CSV:")
    print(csv_file)
    print(f"\nFinished Harmony PC{n_pcs}")
