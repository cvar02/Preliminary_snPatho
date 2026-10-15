# Preliminary_snPatho

This repository contains a staged workflow for processing and exploring single-nucleus RNA-seq data from GIST samples. The analyses cover Cell Ranger and CellBender processing, QC, doublet assessment, sample integration, clustering, and sample-level annotation.

## Branches

- `1_Raw_Processing` — Cell Ranger Multi and CellBender scripts.
- `2_QC_Figures` — QC filtering and diagnostic figures.
- `3_Doublet_Removal` — Python and R doublet analyses and singlet outputs.
- `4_Sample_Integration` — merging, preprocessing, batch correction, and preliminary clustering.
- `5_Clustering` — clustering, validation, annotation, and downsampling analyses.
- `5B_Indiv_Clustering_Investigation` — sample-specific cluster annotations for PATH02, PATH03, and PATH04.

Scripts use project-specific paths and, where specified, SLURM job settings. Review and update these settings for your environment before running.
