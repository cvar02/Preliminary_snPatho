# Clustering

Screens integrated scVI objects, checks cluster consistency, examines marker genes, and tests whether subtype imbalance affects clustering.

## Files

- `1_LeidenClustering.py` — screens scVI 25- and 30-dimensional objects with 10/15/20/30 neighbours and Leiden resolutions 0.2–1.0; seed 0 supports reproducible comparisons.
- `2_ClusteringValidation.py` — compares selected scVI25/scVI30 clusterings and exports the top 50 genes per cluster to assess marker and sample/subtype patterns.
- `3_kMeansValidation.py` — compares K-means solutions with 8, 10, 12, 14, 16, and 18 clusters to the scVI30 Leiden reference using ARI/NMI; seed 0 and 20 initialisations stabilise the comparison.
- `4_ClusteringAnnotation.py` — evaluates the selected scVI30 clustering (`n=30`, resolution 0.6) against curated marker panels and exports figures, tables, and an annotated object; marker lists cover GIST and major immune/stromal populations.
- `5_DownsamplingAnalysis.py` — draws 839 nuclei per subtype with seed 1, then reclusters the balanced set using 30 neighbours and resolution 0.6 in the 30-dimensional scVI space; equal sampling tests sensitivity to subtype imbalance.

These settings support parameter and robustness comparisons; marker-based labels remain provisional pending biological interpretation.
