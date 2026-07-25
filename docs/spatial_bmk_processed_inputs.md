# BMK-ST spatial processed inputs

The BMK-ST spatial analyses start from sequencing-provider processed spatial matrices, not directly from FASTQ files. For each fetal ovary spatial sample, the required processed inputs are:

- `matrix.mtx.gz`
- `features.tsv.gz` or `genes.tsv.gz`
- `barcodes.tsv.gz`
- `barcodes_pos.tsv.gz`

The files should be organized as 10x-style sparse matrix directories, with `barcodes_pos.tsv.gz` containing barcode-level spatial coordinates. These processed files are expected to be deposited with the public data release and referenced in the data availability statement.

The repository code reads these processed inputs, validates raw count structure and coordinate alignment, and writes canonical AnnData H5AD files for downstream cell2location analysis. Raw FASTQ files alone are not sufficient to reproduce these spatial objects unless the sequencing-provider spatial processing pipeline, image inputs, coordinate-generation workflow, software versions, and parameters are also available.
