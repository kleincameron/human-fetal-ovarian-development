#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(tibble)
})

options(bitmapType = "cairo")

project_root <- normalizePath(
  Sys.getenv("FETAL_OVARY_PROJECT_ROOT", unset = getwd()),
  mustWork = TRUE
)

results_root <- normalizePath(
  Sys.getenv(
    "FETAL_OVARY_RESULTS_ROOT",
    unset = file.path(dirname(project_root), "github_code_for_publication_results")
  ),
  mustWork = FALSE
)

input_rds <- file.path(
  results_root,
  "snRNAseq_annotated_object",
  "objects",
  "fetal_ovary_snRNAseq_canonical_annotated.rds"
)

out_root <- file.path(results_root, "cell2location_reference")
out_table_dir <- file.path(out_root, "tables")
out_log_dir <- file.path(out_root, "logs")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

n_hvg <- as.integer(Sys.getenv("CELL2LOC_N_HVG", unset = "5000"))
markers_per_state <- as.integer(Sys.getenv("CELL2LOC_MARKERS_PER_STATE", unset = "200"))
min_cells_per_state <- as.integer(Sys.getenv("CELL2LOC_MIN_CELLS_PER_STATE", unset = "20"))

message("Building cell2location reference signatures")
message("Input object: ", input_rds)
message("Output directory: ", out_root)

if (!file.exists(input_rds)) {
  stop("Annotated snRNA-seq object not found: ", input_rds)
}

if (!requireNamespace("presto", quietly = TRUE)) {
  stop("The presto package is required for fast marker selection. Install/load presto before running this script.")
}

obj <- readRDS(input_rds)
required_meta <- c("fine_subcluster", "major_cell_type")
missing_meta <- setdiff(required_meta, colnames(obj@meta.data))
if (length(missing_meta) > 0) {
  stop("Missing required metadata columns: ", paste(missing_meta, collapse = ", "))
}

DefaultAssay(obj) <- "RNA"

if ("Assay5" %in% class(obj[["RNA"]])) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")

if (!inherits(counts, "dgCMatrix")) {
  counts <- as(counts, "dgCMatrix")
}

cell_meta <- obj@meta.data[colnames(counts), , drop = FALSE] %>%
  rownames_to_column("cell_barcode") %>%
  mutate(
    fine_subcluster = as.character(fine_subcluster),
    major_cell_type = as.character(major_cell_type)
  )

if (any(is.na(cell_meta$fine_subcluster))) {
  stop("Some count-matrix cells are missing fine_subcluster metadata after alignment.")
}

if (length(unique(cell_meta$fine_subcluster)) < 2) {
  stop(
    "Fewer than two fine_subcluster levels found after metadata alignment. ",
    "Check that colnames(counts) match rownames(obj@meta.data)."
  )
}

state_counts <- cell_meta %>%
  count(fine_subcluster, major_cell_type, name = "n_cells") %>%
  arrange(fine_subcluster, desc(n_cells))

major_map <- state_counts %>%
  group_by(fine_subcluster) %>%
  slice_max(order_by = n_cells, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(fine_subcluster) %>%
  select(fine_subcluster, major_cell_type, n_cells)

small_states <- major_map %>%
  filter(n_cells < min_cells_per_state)

if (nrow(small_states) > 0) {
  stop(
    "Some fine_subcluster states have fewer than ",
    min_cells_per_state,
    " cells: ",
    paste(small_states$fine_subcluster, collapse = ", ")
  )
}
write_csv(
  major_map,
  file.path(out_table_dir, "fine_subcluster_to_major_cell_type.csv")
)

message("Selecting variable features.")
obj <- NormalizeData(obj, assay = "RNA", verbose = FALSE)
obj <- FindVariableFeatures(
  obj,
  assay = "RNA",
  selection.method = "vst",
  nfeatures = n_hvg,
  verbose = FALSE
)

hvg <- VariableFeatures(obj)

message("Running presto marker selection by fine_subcluster.")
marker_tbl <- presto::wilcoxauc(obj, group_by = "fine_subcluster")

gene_col <- if ("feature" %in% colnames(marker_tbl)) "feature" else "gene"
group_col <- if ("group" %in% colnames(marker_tbl)) "group" else "fine_subcluster"

if (!all(c(gene_col, group_col) %in% colnames(marker_tbl))) {
  stop("Could not identify gene/group columns in presto output.")
}

marker_tbl <- marker_tbl %>%
  rename(
    gene = all_of(gene_col),
    fine_subcluster = all_of(group_col)
  )

if (!"padj" %in% colnames(marker_tbl)) {
  marker_tbl$padj <- NA_real_
}
if (!"logFC" %in% colnames(marker_tbl)) {
  marker_tbl$logFC <- 0
}
if (!"pct_in" %in% colnames(marker_tbl)) {
  marker_tbl$pct_in <- 0
}

marker_tbl <- marker_tbl %>%
  mutate(
    gene = as.character(gene),
    fine_subcluster = as.character(fine_subcluster)
  )

selected_markers <- marker_tbl %>%
  filter(!is.na(gene)) %>%
  filter(is.na(padj) | padj <= 0.05) %>%
  filter(logFC > 0) %>%
  filter(pct_in >= 0.05) %>%
  group_by(fine_subcluster) %>%
  arrange(desc(logFC), desc(pct_in), .by_group = TRUE) %>%
  slice_head(n = markers_per_state) %>%
  ungroup()

write_csv(
  selected_markers,
  file.path(out_table_dir, "selected_reference_markers_by_fine_subcluster.csv")
)
excluded_gene <- grepl("^MT-", rownames(counts), ignore.case = FALSE) |
  grepl("^RPL", rownames(counts), ignore.case = FALSE) |
  grepl("^RPS", rownames(counts), ignore.case = FALSE) |
  grepl("^MRPL", rownames(counts), ignore.case = FALSE) |
  grepl("^MRPS", rownames(counts), ignore.case = FALSE)

genes_use <- union(hvg, selected_markers$gene)
genes_use <- intersect(genes_use, rownames(counts))
genes_use <- setdiff(genes_use, rownames(counts)[excluded_gene])
genes_use <- sort(unique(genes_use))

if (length(genes_use) < 500) {
  stop("Too few genes selected for cell2location reference: ", length(genes_use))
}

message("Selected genes: ", length(genes_use))

groups <- factor(cell_meta$fine_subcluster)

message("Fine subcluster levels after metadata alignment: ", nlevels(groups))
message("Cells with fine_subcluster metadata: ", sum(!is.na(groups)))

if (nlevels(groups) < 2) {
  stop("Cannot build reference signatures with fewer than two fine_subcluster levels.")
}

group_mat <- sparse.model.matrix(~ 0 + groups)
colnames(group_mat) <- sub("^groups", "", colnames(group_mat))

cell_counts_by_state <- Matrix::colSums(group_mat)
state_order <- names(cell_counts_by_state)

message("Computing average raw expression per fine_subcluster.")
counts_use <- counts[genes_use, , drop = FALSE]
gene_by_state_sum <- counts_use %*% group_mat
gene_by_state_mean <- sweep(
  as.matrix(gene_by_state_sum),
  2,
  as.numeric(cell_counts_by_state),
  FUN = "/"
)

colnames(gene_by_state_mean) <- state_order
rownames(gene_by_state_mean) <- genes_use

gene_by_state_mean <- gene_by_state_mean[, sort(colnames(gene_by_state_mean)), drop = FALSE]

signature_path <- file.path(out_root, "reference_signatures_fine_subcluster.csv")
write_csv(
  as.data.frame(gene_by_state_mean) %>% rownames_to_column("gene"),
  signature_path
)

gene_source <- tibble(
  gene = genes_use,
  in_hvg = genes_use %in% hvg,
  in_selected_marker = genes_use %in% selected_markers$gene,
  excluded_mito_or_ribo = FALSE
)

write_csv(
  gene_source,
  file.path(out_table_dir, "genes_used_HVG_union_markers.csv")
)

summary_lines <- c(
  "cell2location reference signature build",
  paste("Input object name:", basename(input_rds)),
  paste("Number of cells:", ncol(obj)),
  paste("Number of genes in object:", nrow(obj)),
  paste("Fine subcluster states:", paste(colnames(gene_by_state_mean), collapse = ", ")),
  paste("Number of selected genes:", length(genes_use)),
  paste("Number of HVGs requested:", n_hvg),
  paste("Markers per state requested:", markers_per_state),
  paste("Signature file:", basename(signature_path))
)

writeLines(summary_lines, file.path(out_log_dir, "cell2location_reference_summary.txt"))

sink(file.path(out_log_dir, "cell2location_reference_sessionInfo.txt"))
print(sessionInfo())
sink()

message("Done.")
message("Reference signatures: ", signature_path)
