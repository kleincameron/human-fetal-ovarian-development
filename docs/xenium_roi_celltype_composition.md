# Xenium ROI cell-type composition

This workflow quantifies broad cell-type composition across manually defined Xenium regions of interest.

Tracked ROI metadata:

- `metadata/xenium_medulla_roi_coordinates.csv`
- `metadata/xenium_primordial_follicle_roi_coordinates.csv`
- `metadata/xenium_degenerating_follicle_rois/degen*_coordinates.csv`

ROI assignment priority:

1. Degenerating follicle ROI
2. Inner cortex / primordial follicle ROI
3. Medulla
4. Outer cortex / cortex outside the primordial follicle ROI

The coordinate system is expected to match the cell centroid coordinates in the deposited Xenium `cells.csv.gz` file.
