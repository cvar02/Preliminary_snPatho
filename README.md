# Individual sample clustering

Runs the same sample-level annotation workflow separately for PATH02, PATH03, and PATH04, preserving a marker panel and clustering configuration for each sample.

## Files

- `PATH02.py` — processes PATH02 and produces an annotation figure, marker dotplot, and top-50 DEG table; uses 3,000 HVGs, 50 PCs, 10 neighbours, Leiden resolution 0.6, and seed 0.
- `PATH03.py` — produces the same outputs for PATH03; uses 3,000 HVGs, 30 PCs, 20 neighbours, resolution 0.8, and seed 0.
- `PATH04.py` — produces the same outputs for PATH04; uses 3,000 HVGs, 50 PCs, 20 neighbours, resolution 0.8, and seed 0.

All three normalise to 10,000 counts and cap scaled values at 10; these shared settings make sample-level preprocessing comparable. Each script uses sample-specific marker panels, displays up to 8 genes per marker category, and uses the top 100 ranked genes as annotation evidence. PATH02 requests ranked DEG output; PATH03 and PATH04 request positive significant DEGs.
