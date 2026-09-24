# Xenium ROI expression heatmaps by annotation

This workflow compares Xenium expression across manually defined cortical regions for selected annotated cell states.

Inputs are generated or tracked within the public GitHub workflow structure:

- `config/xenium_samples.csv`, local ignored manifest pointing to the downloaded Xenium bundle
- `metadata/xenium_medulla_roi_coordinates.csv`
- `metadata/xenium_primordial_follicle_roi_coordinates.csv`
- `metadata/xenium_degenerating_follicle_rois/degen*_coordinates.csv`
- `<FETAL_OVARY_RESULTS_ROOT>/xenium_annotated_object/objects/fetal_ovary_xenium_annotated.rds`

Regions used for heatmaps:

1. Outer cortex
2. Inner cortex / primordial follicle ROI
3. Degenerating follicles

Medulla-assigned cells are excluded from the heatmap comparisons. Degenerating follicle ROIs take priority over the inner cortex / primordial follicle ROI when ROIs overlap.

Differential-expression contrasts:

1. Outer cortex vs. inner cortex
2. Degenerating follicle ROI vs. whole cortex, where whole cortex is outer cortex plus inner cortex

The figure keeps the row ordering by three enrichment categories but does not display the enrichment-category strip labels on the heatmaps.
