#!/usr/bin/env python3
import os
import scanpy as sc
import pandas as pd
import matplotlib.pyplot as plt

# ============================================================
# 1. PATHWAYS
# ============================================================
SAMPLE_DIR = os.getcwd()
SAMPLE = os.path.basename(SAMPLE_DIR)
RAW_MATRIX = os.path.join(
   "/scratch/cv1u24/GIST/snPatho",
   SAMPLE,
   "outs",
   "raw_feature_bc_matrix.h5"
)
CELLBENDER_MATRIX = os.path.join(
   "/home/cv1u24/GIST/1_Raw_Processing/CellBender",
   SAMPLE,
   "CellBender_Outputs",
   f"{SAMPLE}_cellbender_filtered.h5"
)
OUTPUT_DIR = os.path.join(SAMPLE_DIR, "Outputs")
FIG_DIR = os.path.join(OUTPUT_DIR, "Figures")
FILTERED_DIR = os.path.join(OUTPUT_DIR, "Filtered_h5ad")
SUMMARY_DIR = os.path.join(OUTPUT_DIR, "Summary")
for directory in [FIG_DIR, FILTERED_DIR, SUMMARY_DIR]:
   os.makedirs(directory, exist_ok=True)

# ============================================================
# 2. QC THRESHOLDS
#    Based on published snPATHO-seq methodology
# ============================================================
MIN_COUNTS = 200
MAX_COUNTS = 8000
MAX_MT = 10

# ============================================================
# 3. QC FUNCTIONS
# ============================================================
def add_qc_metrics(adata):
   adata.var_names_make_unique()
   # Identify mitochondrial genes
   adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
   sc.pp.calculate_qc_metrics(
       adata,
       qc_vars=["mt"],
       percent_top=None,
       log1p=False,
       inplace=True
   )
   return adata

def apply_qc_filter(adata):
   keep = (
       (adata.obs["total_counts"] >= MIN_COUNTS) &
       (adata.obs["total_counts"] <= MAX_COUNTS) &
       (adata.obs["pct_counts_mt"] < MAX_MT)
   )
   return adata[keep].copy()

def plot_qc(adata, title, output_path):
   fig, axes = plt.subplots(1, 3, figsize=(16, 5))
   # Total counts
   axes[0].violinplot(
       adata.obs["total_counts"],
       showmeans=True,
       showmedians=True
   )
   axes[0].axhline(MIN_COUNTS, linestyle="--", linewidth=1)
   axes[0].axhline(MAX_COUNTS, linestyle="--", linewidth=1)
   axes[0].set_title("Total counts")
   axes[0].set_ylabel("Counts")
   axes[0].set_xticks([])
   # Mitochondrial percentage
   axes[1].violinplot(
       adata.obs["pct_counts_mt"],
       showmeans=True,
       showmedians=True
   )
   axes[1].axhline(MAX_MT, linestyle="--", linewidth=1)
   axes[1].set_title("Mitochondrial content")
   axes[1].set_ylabel("Percent mitochondrial counts")
   axes[1].set_xticks([])
   # Counts versus genes detected
   axes[2].scatter(
       adata.obs["total_counts"],
       adata.obs["n_genes_by_counts"],
       s=3,
       alpha=0.4
   )
   axes[2].axvline(MIN_COUNTS, linestyle="--", linewidth=1)
   axes[2].axvline(MAX_COUNTS, linestyle="--", linewidth=1)
   axes[2].set_title("Counts vs genes detected")
   axes[2].set_xlabel("Total counts")
   axes[2].set_ylabel("Genes detected")
   fig.suptitle(title)
   plt.tight_layout()
   plt.savefig(output_path, dpi=300, bbox_inches="tight")
   plt.close()

# ============================================================
# 4. CHECK INPUT FILES
# ============================================================
if not os.path.exists(RAW_MATRIX):
   raise FileNotFoundError(
       f"Raw matrix not found: {RAW_MATRIX}"
   )
if not os.path.exists(CELLBENDER_MATRIX):
   raise FileNotFoundError(
       f"CellBender matrix not found: {CELLBENDER_MATRIX}"
   )
print(f"\nProcessing {SAMPLE}")
print(f"Raw matrix: {RAW_MATRIX}")
print(f"CellBender matrix: {CELLBENDER_MATRIX}")

# ============================================================
# 5. RAW MATRIX
# ============================================================
raw = sc.read_10x_h5(RAW_MATRIX)
raw = add_qc_metrics(raw)
raw_fig = os.path.join(
   FIG_DIR,
   f"{SAMPLE}_01_raw_qc.png"
)
plot_qc(
   raw,
   f"{SAMPLE} Raw 10X",
   raw_fig
)

# ============================================================
# 6. CELLBENDER MATRIX BEFORE QC
# ============================================================
cb = sc.read_10x_h5(CELLBENDER_MATRIX)
cb = add_qc_metrics(cb)
cb.obs["sample"] = SAMPLE
cellbender_fig = os.path.join(
   FIG_DIR,
   f"{SAMPLE}_02_cellbender_qc.png"
)
plot_qc(
   cb,
   f"{SAMPLE} CellBender",
   cellbender_fig
)

# ============================================================
# 7. APPLY QC
# ============================================================
qc = apply_qc_filter(cb)
qc.obs["sample"] = SAMPLE
qc_fig = os.path.join(
   FIG_DIR,
   f"{SAMPLE}_03_qc_filtered.png"
)
plot_qc(
   qc,
   f"{SAMPLE} QC filtered",
   qc_fig
)

# ============================================================
# 8. SAVE FILTERED OBJECT
# ============================================================
filtered_output = os.path.join(
   FILTERED_DIR,
   f"{SAMPLE}_cellbender_qc_filtered.h5ad"
)
qc.write_h5ad(filtered_output)

# ============================================================
# 9. SAVE SUMMARY
# ============================================================
summary = pd.DataFrame([{
   "sample": SAMPLE,
   "raw_barcodes": raw.n_obs,
   "cellbender_barcodes": cb.n_obs,
   "qc_barcodes": qc.n_obs,
   "removed_by_qc": cb.n_obs - qc.n_obs,
   "min_counts": MIN_COUNTS,
   "max_counts": MAX_COUNTS,
   "max_mt_percent": MAX_MT,
   "genes_filtered": False
}])
summary_output = os.path.join(
   SUMMARY_DIR,
   f"{SAMPLE}_cellbender_qc_summary.csv"
)
summary.to_csv(summary_output, index=False)

# ============================================================
# 10. REPORT
# ============================================================
print("\nQC summary:")
print(f"CellBender nuclei: {cb.n_obs}")
print(f"QC-passed nuclei: {qc.n_obs}")
print(f"Removed by QC: {cb.n_obs - qc.n_obs}")
print("\nQC thresholds:")
print(f"Total counts: {MIN_COUNTS}-{MAX_COUNTS}")
print(f"Mitochondrial: <{MAX_MT} percent")
print("Gene filtering: None")
print("\nSaved outputs:")
print(raw_fig)
print(cellbender_fig)
print(qc_fig)
print(filtered_output)
print(summary_output)