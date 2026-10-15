# Doublet removal

Compares doublet calls from Python and R methods and provides batch launchers for sample-level analyses.

## Files

- `pythonDoublet.py` — runs Scrublet and DoubletDetection on CellBender-filtered counts using each tool’s defaults, then exports per-barcode calls and scores for comparison.
- `scDBLFinder.R` — runs `scDblFinder`, setting seed 1 for reproducible calls, and saves barcode-level results, summary statistics, and an R object.
- `1_pythonDBL_All.sh` — runs the Python doublet script across sample directories.
- `1_scDBLfinder_All.sh` — runs the R doublet script across sample directories.
- `Run_DBL_Investigation.sh` — launches `Scripts/3_DBL_Investigation.py` for a broader method comparison.
- `Run_DBLremove_and_QC.sh` — launches `Scripts/4_DBL_Remove_and_QC.py` to create doublet-filtered, QC outputs.
- `Run_scatterVisualization.sh` — launches `Scripts/2_scatterVisualization.py` for visual comparison of calls.
