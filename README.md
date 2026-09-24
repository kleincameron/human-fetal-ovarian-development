# Human Fetal Ovarian Development

Code and analysis pipeline for the human fetal ovary development study. This repository contains the analysis workflows used for fetal snRNA-seq, fetal–adult integration, pseudotime and differentiation-potential analyses, cell–cell communication, differential expression and pathway analysis, BMK-ST spatial transcriptomics, Xenium analysis, and manuscript figure generation.

---

## 1. System Requirements

### Operating System

| Item | Version / Detail |
|---|---|
| OS | Linux / HPC environment |
| R | 4.4.3 used for the public R workflows |
| Python | Python 3; workflow-specific environments are recommended |
| Shell | Bash |

### Main R Dependencies

| Package / group | Usage |
|---|---|
| Seurat, SeuratObject, Matrix | Single-cell and Xenium object handling, normalization, clustering, label transfer, differential expression |
| SingleCellExperiment, scDblFinder | Single-cell interoperability and doublet detection |
| presto | Fast marker analysis |
| monocle3 | Pseudotime trajectory inference |
| CellChat, igraph | Cell–cell communication analysis and network visualization |
| edgeR | Fetal-versus-adult pseudobulk differential expression |
| fgsea | Hallmark and Reactome enrichment analysis |
| AnnotationDbi, org.Hs.eg.db, GO.db | Gene Ontology enrichment analysis |
| dplyr, tidyr, readr, tibble, stringr | Data manipulation and input/output |
| ggplot2, patchwork, cowplot, scales, ggtext, ggrepel, pheatmap | Figure generation |
| reticulate, cellxgene.census | Adult CELLxGENE reference handling |
| future, ragg | Workflow execution and high-resolution plotting |

`msigdbr` can optionally be used as a Hallmark gene-set source when a local Hallmark GMT file is not supplied.

### Main Python Dependencies

| Package / group | Usage |
|---|---|
| numpy, pandas, scipy | Numerical, tabular, and sparse-matrix operations |
| scanpy, anndata | BMK-ST spatial object handling |
| matplotlib | Spatial QC and visualization |
| torch, lightning, scvi-tools, cell2location | Cell2location spatial mapping |

The CytoTRACE2 workflow additionally requires a working `cytotrace2` command-line installation.

### Hardware

- A workstation or HPC node with substantial memory is recommended for full-dataset analyses; ≥64 GB RAM is a practical starting point.
- A CUDA-capable GPU can accelerate cell2location but is not required for all workflows.

---

## 2. Installation Guide

Clone the repository:

```bash
git clone https://github.com/kleincameron/human-fetal-ovarian-development.git
cd human-fetal-ovarian-development
```

### R Packages

Install the core CRAN and Bioconductor dependencies required for the workflows you plan to run. For example:

```r
install.packages(c(
  "Seurat", "SeuratObject", "Matrix",
  "dplyr", "tidyr", "readr", "tibble", "stringr",
  "ggplot2", "patchwork", "cowplot", "scales",
  "ggtext", "ggrepel", "pheatmap", "future", "igraph",
  "reticulate", "ragg", "msigdbr"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c(
  "SingleCellExperiment", "scDblFinder", "edgeR", "fgsea",
  "AnnotationDbi", "org.Hs.eg.db", "GO.db"
))
```

Specialized packages including `presto`, `monocle3`, `CellChat`, and `cellxgene.census` should be installed according to their current project documentation.

### Python Packages

A minimal environment for BMK-ST object construction and QC requires:

```bash
pip install numpy pandas scipy scanpy anndata matplotlib
```

Cell2location additionally requires a compatible PyTorch/scvi-tools environment:

```bash
pip install torch lightning scvi-tools cell2location
```

---

## 3. Instructions for Use

### Input Data

The repository expects the following external data resources:

- **Fetal snRNA-seq:** data associated with GSA-Human accession **HRA019091**.
- **Adult ovary reference:** Tabula Sapiens adult ovary from CELLxGENE  
  - collection: `e5f58829-1a66-40b5-a624-9046778e74f5`  
  - dataset: `584027d5-32d7-4696-9424-f61134ff2aa7`
- **BMK-ST spatial data:** provider-processed 10x-style matrices plus `barcodes_pos.tsv.gz`.
- **Xenium:** downloaded Xenium output bundle for the profiled fetal ovary sample.
- **Gene sets:** MSigDB Hallmark and Reactome GMT files for the fetal–adult enrichment workflow.

The fetal snRNA-seq workflow begins from filtered 10x-style feature-barcode matrices. FASTQ-to-matrix processing is an upstream prerequisite and is not implemented in this repository.

### Configuration

Run analyses from the repository root or set:

```bash
export PROJECT_ROOT=/path/to/human-fetal-ovarian-development
export FETAL_OVARY_DATA_ROOT=/path/to/fetal_ovary_data
export FETAL_OVARY_RESULTS_ROOT=/path/to/fetal_ovary_results
```

For BMK-ST data stored under a separate processed-data directory, the object builder also supports:

```bash
export FETAL_OVARY_SPATIAL_BMK_ROOT=/path/to/processed_spatial
```

Local machine-specific configuration can be supplied through an ignored `config/paths_local.R`.

Example manifests are provided for snRNA-seq, spatial, and Xenium inputs:

```text
config/snrna_samples_example.csv
config/xenium_samples_example.csv
metadata/spatial_bmk_samples_example.csv
```

Copy and edit the corresponding local manifest before running workflows that require it.

### Analysis Overview

| Directory | Analysis |
|---|---|
| `00_build_objects/` | Build canonical fetal snRNA-seq, adult reference, fetal–adult integrated, BMK-ST, and Xenium objects |
| `01_snRNAseq/` | Fetal snRNA-seq QC, annotation figures, proportions, pseudotime, CytoTRACE2, CellChat, differential expression, and enrichment |
| `02_fetal_adult/` | Fetal-versus-adult integration, composition, pseudobulk edgeR, volcano plots, and fgsea |
| `03_spatial/` | BMK-ST reference signatures, cell2location mapping, spatial QC, and cell-type highlight plots |
| `04_xenium/` | Xenium QC, medulla/cortex analyses, cortical centrality, ROI composition and expression, and differential expression |
| `config/` | Shared labels, colors, plotting settings, path configuration, and example manifests |
| `metadata/` | Public cell annotations, ROI coordinates, dataset manifests, and spatial/Xenium metadata |
| `docs/` | Additional documentation for spatial and Xenium inputs and ROI definitions |

### Reproduction of Results

A typical analysis order is:

```text
00_build_objects/
    ↓
01_snRNAseq/
    ↓
02_fetal_adult/
03_spatial/
04_xenium/
```

Build the canonical fetal snRNA-seq object and add public annotations first:

```bash
Rscript 00_build_objects/00_build_snRNAseq_canonical_object.R
Rscript 00_build_objects/01_add_snRNAseq_cell_annotations.R
```

Run the fetal snRNA-seq analyses in `01_snRNAseq/` as needed. Some downstream scripts depend on earlier outputs; for example:

```text
14_generate_deg_tables_celltype_subcluster.R
    ├── 15_gobp_ora_heatmaps_by_celltype.R
    ├── 16_subcluster_marker_expression_heatmaps.R
    ├── 17_two_subcluster_defining_deg_lollipop.R
    └── 18_hallmark_ora_heatmaps_by_celltype.R
```

For fetal–adult analyses, first prepare the adult reference and integrated object:

```bash
Rscript 00_build_objects/02_prepare_adult_cellxgene_reference.R
Rscript 00_build_objects/03_integrate_fetal_adult_snRNAseq.R
```

For BMK-ST spatial analysis:

```bash
python 00_build_objects/04_build_spatial_bmk_canonical_h5ad.py
Rscript 03_spatial/00_build_cell2location_reference_signatures.R
python 03_spatial/01_run_cell2location_bmk.py
python 03_spatial/02_qc_spatial_bmk.py
python 03_spatial/03_plot_spatial_bmk_celltype_highlights.py
```

For Xenium analysis, first build the annotated object and then run the workflows in `04_xenium/`:

```bash
Rscript 00_build_objects/05_build_xenium_annotated_object.R
```

Results are written outside the repository under `FETAL_OVARY_RESULTS_ROOT`. Workflow-specific directories contain generated objects, tables, figures, and logs.

---

## License

This project is licensed under the [MIT License](LICENSE).

Third-party datasets, software, and external resources remain subject to their respective access and reuse terms.
