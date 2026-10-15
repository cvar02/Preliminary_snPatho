# Sample integration

Combines filtered sample objects, preprocesses counts, evaluates dimensionality and batch correction, and screens preliminary clusters.

## Files

- `1_mergeSinglets.py` — labels samples by subtype, records nuclei counts, and merges final doublet-filtered objects; the 500-nucleus minimum is an inclusion threshold to avoid carrying very low-yield samples into joint analysis.
- `2A_Preprocessing.py` — normalises to 10,000 counts, selects 3,000 highly variable genes, scales values to ±10, and computes 50 PCs for downstream comparisons.
- `2B_PCAdiagnostics .py` — compares 10, 15, 20, 25, and 30 PCs while holding neighbours at 15, Leiden resolution at 0.5, and seed at 42, so dimensionality can be assessed under consistent settings.
- `3A_Harmony.py` — corrects sample-associated variation with Harmony using 20, 25, or 30 PCs; it uses sample as the batch key, 20 Harmony iterations, and 15 neighbours for comparable embeddings.
- `3B_scVI.py` — trains count-based scVI models with 20, 25, or 30 latent dimensions and `sample` as the batch key; this tests latent-space dimensionality while preserving raw counts.
- `4_PrelimClustering.py` — compares Harmony and scVI clustering across 15/25/50 neighbours and Leiden resolutions 0.2–1.0, using 30 dimensions to explore cluster granularity and representation choice.

The component, neighbour, and resolution grids are comparative screens, not claims that one setting is universally optimal.
