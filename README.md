# QC and figures

Applies per-nucleus count and mitochondrial-content filters to CellBender output, then saves QC figures, summaries, and filtered AnnData objects.

## Files

- `QCandFigures.py` — calculates QC metrics and applies thresholds of 200–8,000 total counts and <10% mitochondrial counts; comments identify these cut-offs as based on published snPATHO-seq methodology. It writes diagnostic plots and filtered data without gene filtering.
- `QC_Figures_All.sh` — submits the Python workflow for each `PATH0*` directory as a separate SLURM job

The thresholds are explicit analysis settings and should be reviewed if the assay, tissue quality, or sample set changes.
