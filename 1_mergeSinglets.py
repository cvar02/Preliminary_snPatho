#!/usr/bin/env python
import os
import glob
import pandas as pd
import scanpy as sc
import anndata as ad

# ############################################################
# 1. PATHS
# Define the source sample directory and output locations.
# ############################################################
INPUT_BASE = (
    "/home/cv1u24/GIST/3_Doublet_Removal"
)
OUTPUT_BASE = (
    "/home/cv1u24/GIST/4_Preprocessing_and_Integration"
)
H5AD_OUT = os.path.join(
    OUTPUT_BASE,
    "Outputs",
    "h5ad"
)
CSV_OUT = os.path.join(
    OUTPUT_BASE,
    "Outputs",
    "CSVs"
)
os.makedirs(
    H5AD_OUT,
    exist_ok=True
)
os.makedirs(
    CSV_OUT,
    exist_ok=True
)

# ############################################################
# 2. SAMPLE METADATA
# Map sample identifiers to GIST subtype labels and set the inclusion threshold.
# ############################################################
GIST_SUBTYPE_MAP = {
    "PATH01": "WT_undefined",
    "PATH02": "KIT",
    "PATH03": "NF1_associated",
    "PATH04": "SDH_deficient"
}

# Minimum number of nuclei required for inclusion
MIN_CELLS = 500

# ############################################################
# 3. FIND SAMPLE DIRECTORIES
# Locate the PATH0* sample folders to process.
# ############################################################
sample_paths = sorted(
    glob.glob(
        os.path.join(
            INPUT_BASE,
            "PATH0*"
        )
    )
)
if len(sample_paths) == 0:
    raise ValueError(
        f"No PATH0* sample directories found in:\n"
        f"{INPUT_BASE}"
    )
print("\n")
print("=" * 70)
print("MERGING FINAL SINGLET OBJECTS")
print("=" * 70)

# ############################################################
# 4. LOAD SAMPLES
# Load each final doublet-filtered object and record sample-level counts.
# ############################################################
adata_list = []
summary_rows = []
for sample_dir in sample_paths:
    sample = os.path.basename(
        os.path.normpath(sample_dir)
    )

    # --------------------------------------------------------
    # Expected final DoubletDetection-filtered object
    # --------------------------------------------------------
    input_h5ad = os.path.join(
        sample_dir,
        "Final",
        f"{sample}_final_doublet_filtered.h5ad"
    )
    if not os.path.exists(input_h5ad):
        print(
            f"\nWARNING: Missing final object for "
            f"{sample}"
        )
        print(
            f"Expected:\n{input_h5ad}"
        )
        continue

# --------------------------------------------------------
# Load object
# --------------------------------------------------------
    print("\n")
    print("-" * 70)
    print(f"Loading {sample}")
    print("-" * 70)
    adata = sc.read_h5ad(
        input_h5ad
    )
    n_cells = adata.n_obs
    n_genes = adata.n_vars

    # --------------------------------------------------------
    # Assign subtype
    # --------------------------------------------------------
    gist_subtype = GIST_SUBTYPE_MAP.get(
        sample,
        "Unknown"
    )

    # --------------------------------------------------------
    # Determine whether sample passes threshold
    # --------------------------------------------------------
    keep_sample = (
        n_cells >= MIN_CELLS
    )
    print(
        f"Cells: {n_cells}"
    )
    print(
        f"Genes: {n_genes}"
    )
    print(
        f"GIST subtype: {gist_subtype}"
    )
    print(
        f"Keep for merged analysis: "
        f"{keep_sample}"
    )

    # --------------------------------------------------------
    # Record sample summary
    # --------------------------------------------------------
    summary_rows.append({
        "sample":
            sample,
        "gist_subtype":
            gist_subtype,
        "n_cells_after_doublet_removal":
            n_cells,
        "n_genes":
            n_genes,
        "minimum_cell_threshold":
            MIN_CELLS,
        "passed_minimum_cell_threshold":
            keep_sample
    })

    # --------------------------------------------------------
    # Exclude low-cell samples
    # --------------------------------------------------------
    if not keep_sample:
        print(
            f"\nExcluding {sample}: "
            f"{n_cells} nuclei "
            f"(< {MIN_CELLS})"
        )
        continue

# ############################################################
# 5. ADD METADATA
# Attach sample, subtype, and doublet-filtering metadata.
# ############################################################
    adata.obs["sample"] = sample
    adata.obs["gist_subtype"] = gist_subtype

    # --------------------------------------------------------
    # Record that this object passed doublet filtering
    # --------------------------------------------------------
    adata.obs[
        "doublet_filter"
    ] = "DoubletDetection"

    # ############################################################
    # 6. PRESERVE RAW COUNTS
    # Ensure each sample object has a counts layer before merging.
    # ############################################################
    if "counts" not in adata.layers:
        adata.layers[
            "counts"
        ] = adata.X.copy()

# ############################################################
# 7. ADD TO MERGE LIST
# Collect samples that pass the minimum-nuclei threshold.
# ############################################################
    adata_list.append(
        adata
    )

# ############################################################
# 8. SAVE SAMPLE SUMMARY
# Write the per-sample inclusion summary to CSV.
# ############################################################
summary_df = pd.DataFrame(
    summary_rows
)
summary_csv = os.path.join(
    CSV_OUT,
    "sample_merge_summary.csv"
)
summary_df.to_csv(
    summary_csv,
    index=False
)
print("\n")
print("=" * 70)
print("SAMPLE SUMMARY")
print("=" * 70)
print(
    summary_df.to_string(
        index=False
    )
)
print(
    f"\nSaved sample summary:\n"
    f"{summary_csv}"
)

# ############################################################
# 9. CHECK THAT SAMPLES REMAIN
# Stop if no sample passed the inclusion threshold.
# ############################################################
if len(adata_list) == 0:
    raise ValueError(
        "No samples passed the minimum cell threshold."
    )

# ############################################################
# 10. MERGE SAMPLES
# Concatenate the retained sample objects using shared genes.
# ############################################################
print("\n")
print("=" * 70)
print("MERGING SAMPLES")
print("=" * 70)
sample_keys = [
    x.obs["sample"].iloc[0]
    for x in adata_list
]
adata_merged = ad.concat(
    adata_list,
    join="inner",
    label="sample_merge_id",
    keys=sample_keys,
    index_unique="-"
)
print(
    "\nMerged object created."
)
print(
    f"Total nuclei: "
    f"{adata_merged.n_obs}"
)
print(
    f"Genes: "
    f"{adata_merged.n_vars}"
)

# ############################################################
# 11. SAVE MERGED UNINTEGRATED OBJECT
# Write the merged, unintegrated AnnData object.
# ############################################################
merged_output = os.path.join(
    H5AD_OUT,
    "merged_singlets_unintegrated.h5ad"
)
adata_merged.write_h5ad(
    merged_output
)
print(
    f"\nSaved merged object:\n"
    f"{merged_output}"
)

# ############################################################
# 12. FINAL SUMMARY
# Print the samples and nuclei included in the merged object.
# ############################################################
print("\n")
print("=" * 70)
print("MERGE COMPLETE")
print("=" * 70)
print(
    f"\nSamples included: "
    f"{len(adata_list)}"
)
print(
    f"Total nuclei: "
    f"{adata_merged.n_obs}"
)
print(
    f"Genes: "
    f"{adata_merged.n_vars}"
)
print(
    f"\nOutput:\n"
    f"{merged_output}"
)
print("\nDone.")
