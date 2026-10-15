# Raw processing

Processes PATH01–PATH04 with Cell Ranger Multi and CellBender.

## Files

- `CellRangerMulti.sh` — checks for each sample’s Cell Ranger CSV, runs cellranger multi, then aggregates molecule files with normalize=mapped to make sample-level outputs comparable.
- `RunCellBender.sh` — runs `cellbender remove-background` on one sample’s raw feature-barcode matrix using CellBender defaults, so background correction is performed without custom parameter overrides.
- `RunCellBenderAll.sh` — submits each sample’s `RunCellBender.sh` as an independent SLURM job.
