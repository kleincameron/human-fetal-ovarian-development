# Xenium cortical centrality boundaries

`metadata/xenium_cortex_outer_boundary.csv` and `metadata/xenium_cortex_inner_boundary.csv` contain manually defined cortical boundary polylines used for Xenium cortical centrality analyses.

For each analyzed Xenium cell, the workflow computes:

- distance to the outer cortical boundary
- distance to the inner cortical boundary
- normalized cortical depth, where 0 is closest to the outer boundary and 1 is closest to the inner boundary
- cortical centrality, where 0 is near either cortical edge and 1 is near the midpoint between the two boundaries

The coordinate system is expected to match the cell centroid coordinates in the deposited Xenium `cells.csv.gz` file.
