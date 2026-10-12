#!/usr/bin/env python3
"""Create the cluster annotation figure, marker dotplot and top-50 DEG table for PATH02."""

# ============================================================
# 1. IMPORTS AND GENERAL SETTINGS
# Load the analysis and plotting libraries and set shared parameters.
# ============================================================
from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

N_HVG = 3000
N_PCS = 50
RANDOM_STATE = 0
TARGET_SUM = 10_000
TOP_N = 50
FIGURE_MARKER_LIMIT = 8
ANNOTATION_EVIDENCE_N = 100

# ============================================================
# 2. SAMPLE CONFIGURATION
# Define the input, selected clustering and output-root paths for this sample.
# ============================================================
SAMPLE = "PATH02"
INPUT_H5AD = (
    "/iridisfs/home/cv1u24/GIST/3_Doublet_Removal/PATH02/Final/PATH02_final_doublet_filtered.h5ad"
)
BASE_DIR = (
    "/iridisfs/home/cv1u24/GIST/5B_Indiv_Clustering_Investigation/PATH02"
)
N_NEIGHBORS = 10
LEIDEN_RESOLUTION = 0.6
CLUSTER = "leiden_PCA_n10_r0_6"
DEG_MODE = "ranked"

# ============================================================
# 3. SAMPLE-SPECIFIC MARKER PANELS
# Keep annotation evidence and dotplot genes explicit for this sample.
# ============================================================
ANNOTATION_MARKERS = {
   "GIST-associated": [
       "KIT", "PDGFRA", "ANO1", "CD34",
       "GPR20", "PTP4A3", "ETV1", "PRKCQ"
   ],
   "Fibroblast / stromal": [
       "DCN", "COL1A1", "COL1A2", "COL3A1",
       "COL5A1", "COL6A1", "COL6A2", "COL6A3",
       "COL14A1", "LUM", "SFRP2", "FBLN1",
       "SFRP4", "CFD", "OGN", "DPT",
       "PDGFRA", "VCAN", "C7"
   ],
   "Perivascular / contractile": [
       "ACTA2", "TAGLN", "MYH11", "CNN1",
       "MYLK", "ACTG2", "RGS5", "MCAM",
       "DES", "LMOD1", "NOTCH3", "RBPMS",
       "PRKG1", "SMTN", "CACNA1C", "TNS2",
       "FHL2", "ITGA10"
   ],
   "Macrophage": [
       "CD163", "C1QA", "C1QB", "C1QC",
       "CSF1R", "TYROBP", "CD68", "MS4A6A",
       "STAB1", "MRC1", "FCGR2A"
   ],
   "Inflammatory myeloid": [
       "CSF3R", "S100A8", "S100A9", "TREM1",
       "FPR1", "FCAR", "CXCL8", "AQP9",
       "IL1R2", "TLR2", "FPR2"
   ],
   "T cell": [
       "CD3D", "CD3E", "TRAC", "CD2",
       "CD7", "ZAP70", "CD6", "CD28",
       "CD247", "LCK", "IL7R"
   ],
   "Cytotoxic T / NK": [
       "NKG7", "NCAM1", "NCR1", "GZMA",
       "GZMB", "PRF1", "KLRB1", "CX3CR1"
   ],
   "B cell": [
       "MS4A1", "CD79A", "CD79B", "CD37",
       "CD22", "PAX5", "BANK1", "BLK",
       "TNFRSF13C", "IGHM", "IGKC"
   ],
   "Plasma cell": [
       "MZB1", "JCHAIN", "XBP1", "SDC1",
       "CD38", "TENT5C", "IRF4", "IGHG1",
       "IGHG2", "IGHA1"
   ],
   "Mast cell": [
       "CPA3", "TPSAB1", "TPSB2", "HDC",
       "MS4A2", "CMA1", "GATA2", "HPGDS",
       "TPSG1", "IL1RL1", "SIGLEC6"
   ],
   "Dendritic / APC": [
       "LAMP3", "CCR7", "FSCN1", "CD1C",
       "CLEC10A", "CLEC9A", "XCR1", "IL3RA",
       "CSF2RA", "CSF2RB", "CD84"
   ],
   "Endothelial": [
       "PECAM1", "VWF", "CDH5", "PLVAP",
       "EGFL7", "CLEC14A", "CAVIN2", "PTPRB",
       "CLDN5", "RAMP3", "MMRN2"
   ],
   "Epithelial": [
       "KRT4", "KRT5", "KRT6A", "KRT13",
       "KRT19", "KRT18", "KRT8", "EPCAM",
       "KRT14", "SPRR3", "KRT78"
   ],
   "Skeletal muscle": [
       "MYH7", "ACTA1", "TNNT1", "TNNC1",
       "CKM", "MYL2", "MYL1", "TCAP",
       "TNNI1", "MYBPC1", "CRYAB"
   ]
}

DOTPLOT_MARKERS = {
   "GIST": [
       "KIT", "PDGFRA", "ANO1", "CD34", "GPR20", "PTP4A3"
   ],
   "Fibroblast": [
       "COL1A1", "COL1A2", "COL3A1", "COL6A1",
       "COL6A2", "DCN", "LUM", "FAP", "CFD",
       "C1S", "C1R"
   ],
   "Smooth_muscle": [
       "ACTA2", "TAGLN", "MYH11", "CNN1",
       "DES", "CALD1", "MYL9", "TPM2"
   ],
   "Endothelial": [
       "VWF", "PECAM1", "CDH5", "EMCN", "KDR",
       "ENG", "ESAM", "RAMP2", "EPAS1", "PLVAP"
   ],
   "Macrophage_myeloid": [
       "C1QA", "C1QB", "C1QC", "APOE", "CD68",
       "LGMN", "CTSD", "FCER1G", "TYROBP",
       "AIF1", "LYZ", "MS4A7", "LGALS3"
   ],
   "T_cell": [
       "CD3D", "CD3E", "CD3G", "TRBC1", "TRBC2",
       "IL7R", "LTB"
   ],
   "Cytotoxic_T": [
       "CD8A", "CD8B", "NKG7", "GNLY", "GZMB",
       "GZMH", "GZMK", "PRF1", "CTSW"
   ],
   "NK": [
       "NKG7", "GNLY", "KLRD1", "KLRB1",
       "FCER1G", "TYROBP", "XCL1", "XCL2"
   ],
   "B_cell": [
       "CD79A", "MS4A1", "CD37", "CD74",
       "HLA-DRA", "HLA-DPA1", "CD22", "CD19"
   ],
   "Plasma": [
       "JCHAIN", "MZB1", "SDC1", "XBP1",
       "IGHG1", "IGKC"
   ],
   "Mast": [
       "KIT", "TPSAB1", "TPSB2", "MS4A2",
       "HDC", "CPA3", "GATA2"
   ],
   "Dendritic": [
       "FCER1A", "CD1C", "CLEC10A", "CST3",
       "HLA-DRA", "HLA-DPA1", "HLA-DPB1", "CLEC4A"
   ],
   "Neutrophil": [
       "S100A8", "S100A9", "FCGR3B", "CSF3R",
       "CXCR2", "FPR1", "CTSG", "MMP9"
   ],
   "Pericyte": [
       "RGS5", "PDGFRB", "CSPG4", "MCAM",
       "DES", "ACTA2", "COL4A1", "COL4A2"
   ]
}

CONFIG = {
    "input": INPUT_H5AD,
    "base": BASE_DIR,
    "neighbors": N_NEIGHBORS,
    "resolution": LEIDEN_RESOLUTION,
    "cluster": CLUSTER,
    "deg_mode": DEG_MODE,
    "annotation_markers": ANNOTATION_MARKERS,
    "dotplot_markers": DOTPLOT_MARKERS,
}


# ============================================================
# 4. SMALL HELPER FUNCTIONS
# Provide sorting and count-normalisation helpers.
# ============================================================

def _numeric_key(value):
    """Sort numeric cluster IDs numerically, with a safe string fallback."""
    text = str(value)
    try:
        return (0, int(text))
    except ValueError:
        return (1, text)


def _prepare_counts_and_log_normalize(adata):
    """Use the counts layer when present, then apply the established normalization."""
    if "counts" in adata.layers:
        adata.X = adata.layers["counts"].copy()
    sc.pp.normalize_total(adata, target_sum=TARGET_SUM)
    sc.pp.log1p(adata)
    return adata


# ============================================================
# 5. SELECTED CLUSTERING
# Recreate the configured PCA, neighbour graph, UMAP and Leiden assignment.
# ============================================================

def _run_selected_clustering(sample, cfg):
    """Recreate the selected PCA/UMAP/Leiden workflow from each sample script."""
    adata = sc.read_h5ad(cfg["input"])
    print(f"{sample}: loaded {adata.n_obs:,} nuclei and {adata.n_vars:,} genes")
    _prepare_counts_and_log_normalize(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=N_HVG, flavor="seurat")
    adata = adata[:, adata.var["highly_variable"]].copy()
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, n_comps=N_PCS, svd_solver="arpack", random_state=RANDOM_STATE)

    neighbors_key = f"neighbors_PCA_n{cfg['neighbors']}"
    sc.pp.neighbors(
        adata, n_neighbors=cfg["neighbors"], n_pcs=N_PCS,
        use_rep="X_pca", key_added=neighbors_key
    )
    sc.tl.umap(adata, neighbors_key=neighbors_key, random_state=RANDOM_STATE)
    umap_key = f"X_umap_PCA_n{cfg['neighbors']}"
    adata.obsm[umap_key] = adata.obsm["X_umap"].copy()
    sc.tl.leiden(
        adata, resolution=cfg["resolution"], neighbors_key=neighbors_key,
        key_added=cfg["cluster"], random_state=RANDOM_STATE
    )
    return adata, umap_key


# ============================================================
# 6. WILCOXON DIFFERENTIAL EXPRESSION
# Rank genes across the full gene set using the selected clusters.
# ============================================================

def _rank_full_gene_expression(sample, cfg, cluster_labels):
    """Rank genes across all genes using the sample's saved cluster assignments."""
    expr = sc.read_h5ad(cfg["input"])
    if not expr.obs_names.equals(cluster_labels.index):
        if not expr.obs_names.isin(cluster_labels.index).all():
            raise ValueError(f"{sample}: expression and clustering barcodes do not match.")
        cluster_labels = cluster_labels.reindex(expr.obs_names)
    if cluster_labels.isna().any():
        raise ValueError(f"{sample}: missing cluster labels after barcode alignment.")
    expr.obs[cfg["cluster"]] = cluster_labels.astype(str).to_numpy()
    _prepare_counts_and_log_normalize(expr)

    # Scanpy ranks the full tested gene set; outputs below take either the top 50
    # as originally used for PATH02, or the positive significant subset for 03/04.
    sc.tl.rank_genes_groups(
        expr, groupby=cfg["cluster"], method="wilcoxon", reference="rest",
        n_genes=None, pts=True, key_added="DEG"
    )
    ranked = []
    groups = expr.uns["DEG"]["names"].dtype.names
    for cluster in sorted(groups, key=_numeric_key):
        table = sc.get.rank_genes_groups_df(expr, group=cluster, key="DEG")
        table = table.rename(columns={
            "names": "gene", "scores": "score", "logfoldchanges": "logFC",
            "pvals": "p_value", "pvals_adj": "adjusted_p_value",
            "pct_nz_group": "pct_expressing_in_cluster",
            "pct_nz_reference": "pct_expressing_in_rest",
        })
        table.insert(0, "cluster", str(cluster))
        table.insert(1, "rank", np.arange(1, len(table) + 1))
        ranked.append(table)
    return expr, pd.concat(ranked, ignore_index=True)


# ============================================================
# 7. TOP 50 GENES AND MARKER EVIDENCE
# Apply the sample-specific export criteria and identify annotation evidence.
# ============================================================

def _top50_table(ranked, mode):
    """Retain each source workflow's existing DEG inclusion/ranking criteria."""
    tables = []
    for cluster, table in ranked.groupby("cluster", sort=False):
        if mode == "positive_significant":
            table = table[(table["logFC"] > 0) & (table["adjusted_p_value"] < 0.05)]
        table = table.head(TOP_N).copy()
        table["rank"] = np.arange(1, len(table) + 1)
        tables.append(table)
    result = pd.concat(tables, ignore_index=True)
    columns = [
        "cluster", "rank", "gene", "score", "logFC", "p_value",
        "adjusted_p_value", "pct_expressing_in_cluster", "pct_expressing_in_rest",
    ]
    return result[[column for column in columns if column in result.columns]]


def _marker_matches(ranked, markers):
    """Match the annotation panel to the top 100 positive significant genes."""
    evidence = ranked[(ranked["logFC"] > 0) & (ranked["adjusted_p_value"] < 0.05)]
    matches = []
    for cluster, cluster_table in evidence.groupby("cluster", sort=False):
        cluster_table = cluster_table.head(ANNOTATION_EVIDENCE_N)
        for cell_type, genes in markers.items():
            selected = cluster_table[cluster_table["gene"].isin(genes)].sort_values("rank")
            if not selected.empty:
                matches.append({
                    "cluster": cluster,
                    "candidate_cell_type": cell_type,
                    "n_marker_matches": len(selected),
                    "matched_markers": selected["gene"].tolist(),
                })
    return matches


# ============================================================
# 8. FIGURES
# Create and save the UMAP annotation figure and marker dotplot.
# ============================================================

def _save_annotation_figure(sample, coords, cluster_labels, matches, output_path):
    """Save the cluster UMAP with marker evidence used as annotation support."""
    cluster_text = cluster_labels.astype(str).to_numpy()
    clusters = sorted(pd.unique(cluster_text), key=_numeric_key)
    by_cluster = {}
    for entry in matches:
        by_cluster.setdefault(entry["cluster"], []).append(entry)
    row_units = [max(1, len(by_cluster.get(cluster, []))) for cluster in clusters]
    total_units = max(1, sum(row_units))
    figure_height = max(8.0, 0.38 * total_units + 2.0)
    fig = plt.figure(figsize=(19, figure_height))
    grid = GridSpec(1, 2, width_ratios=[1.0, 1.8], wspace=0.08, figure=fig)
    ax_umap = fig.add_subplot(grid[0, 0])
    palette = plt.get_cmap("tab20")
    for index, cluster in enumerate(clusters):
        mask = cluster_text == cluster
        ax_umap.scatter(coords[mask, 0], coords[mask, 1], s=20, alpha=0.8,
                        color=palette(index % 20))
        ax_umap.text(coords[mask, 0].mean(), coords[mask, 1].mean(), cluster,
                     fontsize=10, fontweight="bold", ha="center", va="center")
    ax_umap.set_title("A. Leiden clusters", fontsize=14, fontweight="bold")
    ax_umap.set_xlabel("UMAP1")
    ax_umap.set_ylabel("UMAP2")
    ax_umap.set_xticks([])
    ax_umap.set_yticks([])

    ax_text = fig.add_subplot(grid[0, 1])
    ax_text.axis("off")
    ax_text.text(0.0, 1.0, "B. Canonical marker evidence among top 100 DEGs",
                 fontsize=14, fontweight="bold", va="top")
    ax_text.text(0.00, 0.94, "Cluster", fontsize=9, fontweight="bold")
    ax_text.text(0.12, 0.94, "Candidate annotations and marker genes",
                 fontsize=9, fontweight="bold")
    y = 0.91
    for index, cluster in enumerate(clusters):
        ax_text.text(0.00, y, cluster, fontsize=8, fontweight="bold", va="top")
        candidates = by_cluster.get(cluster, [])
        if candidates:
            lines = []
            for entry in candidates:
                genes = entry["matched_markers"][:FIGURE_MARKER_LIMIT]
                suffix = "…" if len(entry["matched_markers"]) > FIGURE_MARKER_LIMIT else ""
                lines.append(f"{entry['candidate_cell_type']} ({entry['n_marker_matches']}): "
                             f"{', '.join(genes)}{suffix}")
            text = "\n".join(lines)
        else:
            text = "No canonical markers in top 100"
        ax_text.text(0.12, y, text, fontsize=7.5, va="top")
        y -= 0.86 * row_units[index] / total_units
    fig.suptitle(f"{sample}: Leiden clustering and canonical marker evidence",
                 fontsize=17, fontweight="bold")
    fig.subplots_adjust(top=0.91, bottom=0.05, left=0.04, right=0.98)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _save_dotplot(expr, cluster_key, markers, output_path, sample):
    """Plot the sample's own curated annotation genes, retaining source styling."""
    available = set(expr.var_names)
    genes = list(dict.fromkeys(
        gene for panel in markers.values() for gene in panel if gene in available
    ))
    if not genes:
        raise ValueError(f"{sample}: none of the configured dotplot markers are in the data.")
    missing = [gene for panel in markers.values() for gene in panel if gene not in available]
    if missing:
        print(f"{sample}: {len(set(missing))} configured dotplot markers absent; omitted.")
    plot = sc.pl.dotplot(
        expr, var_names=genes, groupby=cluster_key, standard_scale="var",
        dendrogram=False, show=False, return_fig=True
    )
    plot.savefig(output_path, dpi=300, bbox_inches="tight")


# ============================================================
# 9. RUN WORKFLOW AND SAVE OUTPUTS
# Run this sample and report the three output paths.
# ============================================================

def main():
    """Run the three requested outputs for one configured sample."""
    sample = SAMPLE
    cfg = CONFIG
    if not Path(cfg["input"]).is_file():
        raise FileNotFoundError(f"{sample} input file not found: {cfg['input']}")
    outdir = Path(cfg["base"]) / "annotation_figure"
    outdir.mkdir(parents=True, exist_ok=True)
    figure_path = outdir / f"{sample}_Leiden_annotation.png"
    dotplot_path = outdir / f"{sample}_marker_dotplot.png"
    deg_path = outdir / f"{sample}_top50_DEGs_by_cluster.csv"

    clustered, umap_key = _run_selected_clustering(sample, cfg)
    labels = clustered.obs[cfg["cluster"]].astype(str)
    expr, ranked = _rank_full_gene_expression(sample, cfg, labels)
    top50 = _top50_table(ranked, cfg["deg_mode"])
    top50.to_csv(deg_path, index=False)
    matches = _marker_matches(ranked, cfg["annotation_markers"])
    _save_annotation_figure(sample, clustered.obsm[umap_key], labels, matches, figure_path)
    _save_dotplot(expr, cfg["cluster"], cfg["dotplot_markers"], dotplot_path, sample)

    print(f"\n{sample} complete. Outputs:")
    for path in (figure_path, dotplot_path, deg_path):
        print(f"  {path}")
    print(f"Clusters: {', '.join(sorted(labels.unique(), key=_numeric_key))}")


if __name__ == "__main__":
    main()
