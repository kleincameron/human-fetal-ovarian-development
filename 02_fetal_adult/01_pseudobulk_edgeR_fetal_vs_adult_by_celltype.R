# Pseudobulk fetal-vs-adult differential expression by shared cell type.
#
# This downstream analysis uses the fetal-adult integrated Seurat object for
# harmonized metadata, but it performs edgeR on RNA/counts, not on the
# integrated assay.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(edgeR)
  library(readr)
  library(tibble)
})

set.seed(42)
options(bitmapType = "cairo")

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()),
  mustWork = TRUE
)

if (!file.exists(file.path(project_root, "config", "paths_example.R"))) {
  stop("PROJECT_ROOT does not point to the repository root. Run from the repo root or set PROJECT_ROOT.")
}

paths_local <- file.path(project_root, "config", "paths_local.R")
if (file.exists(paths_local)) {
  source(paths_local)
}

results_base <- if (exists("results_root", inherits = FALSE)) {
  results_root
} else {
  Sys.getenv(
    "FETAL_OVARY_RESULTS_ROOT",
    unset = file.path(dirname(project_root), paste0(basename(project_root), "_results"))
  )
}
results_base <- normalizePath(results_base, mustWork = FALSE)


default_integrated_rds <- file.path(
  results_base,
  "fetal_adult_integration",
  "objects",
  "fetal_adult_integrated_snRNAseq.rds"
)

integrated_rds <- if (exists("fetal_adult_integrated_rds", inherits = FALSE)) {
  fetal_adult_integrated_rds
} else {
  default_integrated_rds
}

results_root <- file.path(results_base, "fetal_adult_pseudobulk_edgeR")

out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

min_cells_per_donor_default <- 20
min_cells_per_donor_germ <- 5
min_donors_per_group <- 2
min_cpm_filter <- 1
min_samples_filter <- 2

get_layers_safe <- function(object, assay = "RNA") {
  as.character(
    tryCatch(
      Layers(object[[assay]]),
      error = function(e) character(0)
    )
  )
}
get_layer_matrix <- function(object, assay = "RNA", layer) {
  matrix <- LayerData(object, assay = assay, layer = layer)

  if (!inherits(matrix, "dgCMatrix")) {
    matrix <- as(matrix, "dgCMatrix")
  }

  matrix
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

summarize_sparse_matrix <- function(matrix, dataset, layer) {
  values <- matrix@x

  fraction_noninteger <- if (length(values) > 0) {
    mean(abs(values - round(values)) > 1e-6)
  } else {
    0
  }

  value_range <- if (length(values) > 0) {
    range(values)
  } else {
    c(0, 0)
  }

  tibble(
    dataset = dataset,
    layer = layer,
    genes = nrow(matrix),
    cells = ncol(matrix),
    nonzero_values = length(values),
    fraction_noninteger = fraction_noninteger,
    min_value = value_range[[1]],
    max_value = value_range[[2]]
  )
}

pick_best_counts_layer_for_dataset <- function(
  object,
  assay = "RNA",
  dataset_value,
  dataset_col = "dataset",
  pattern = "^counts"
) {
  stopifnot(dataset_col %in% colnames(object@meta.data))

  layers <- get_layers_safe(object, assay = assay)
  candidate_layers <- layers[grepl(pattern, layers)]

  if (length(candidate_layers) == 0) {
    stop("No layers matching ", pattern, " in assay ", assay)
  }

  dataset_cells <- rownames(object@meta.data)[
    as.character(object@meta.data[[dataset_col]]) == dataset_value
  ]

  if (length(dataset_cells) == 0) {
    stop("No cells found for dataset_value = ", dataset_value)
  }

  coverage <- vapply(
    candidate_layers,
    function(layer) {
      layer_cells <- colnames(get_layer_matrix(object, assay = assay, layer = layer))
      length(intersect(dataset_cells, layer_cells)) / length(dataset_cells)
    },
    numeric(1)
  )
  best_coverage <- max(coverage)
  tied_layers <- candidate_layers[coverage == best_coverage]

  name_rank <- function(layer) {
    if (layer == "counts") return(0L)
    if (grepl("^counts\\.SeuratProject$", layer)) return(1L)
    if (grepl("^counts\\.adult_cellxgene", layer)) return(1L)
    if (grepl("^counts\\.[^.]+$", layer)) return(2L)
    if (grepl("^counts\\.[0-9]+\\.", layer)) return(50L)
    100L
  }

  chosen_layer <- tied_layers[
    which.min(vapply(tied_layers, name_rank, integer(1)))
  ]

  tibble(
    dataset = dataset_value,
    chosen_layer = chosen_layer,
    coverage = coverage[[chosen_layer]],
    candidate_layers = paste(candidate_layers, collapse = "; ")
  )
}

pseudobulk_sum <- function(counts_gene_cell, group_ids) {
  stopifnot(ncol(counts_gene_cell) == length(group_ids))

  group_ids <- as.factor(group_ids)
  model_matrix <- Matrix::sparse.model.matrix(~ 0 + group_ids)

  colnames(model_matrix) <- sub("^group_ids", "", colnames(model_matrix))

  counts_gene_cell %*% model_matrix
}

write_empty_result <- function(file, celltype, status) {
  empty <- tibble(
    celltype = celltype,
    comparison = "fetal_vs_adult",
    gene = NA_character_,
    logFC_direction = "positive_logFC_indicates_higher_in_fetal",
    logFC = NA_real_,
    logCPM = NA_real_,
    F = NA_real_,
    PValue = NA_real_,
    FDR = NA_real_,
    status = status
  )

  write_csv(empty, file)
}

message("Loading fetal-adult integrated object:\n  ", integrated_rds)
stopifnot(file.exists(integrated_rds))

obj <- readRDS(integrated_rds)
stopifnot(inherits(obj, "Seurat"))

required_cols <- c(
  "dataset",
  "cell_type",
  "adult_cell_type_harmonized",
  "donor_id",
  "sampleID"
)

missing_cols <- setdiff(required_cols, colnames(obj@meta.data))

if (length(missing_cols) > 0) {
  stop("Integrated object is missing required metadata: ", paste(missing_cols, collapse = ", "))
}
DefaultAssay(obj) <- "RNA"

obj$celltype_unified <- ifelse(
  as.character(obj$dataset) == "fetal",
  as.character(obj$cell_type),
  as.character(obj$adult_cell_type_harmonized)
)

obj$donor_unified <- ifelse(
  as.character(obj$dataset) == "fetal",
  as.character(obj$sampleID),
  as.character(obj$donor_id)
)

keep_cells <- which(
  !is.na(obj$dataset) &
    obj$dataset %in% c("fetal", "adult") &
    !is.na(obj$celltype_unified) &
    obj$celltype_unified != "" &
    !is.na(obj$donor_unified) &
    obj$donor_unified != ""
)

obj <- subset(obj, cells = colnames(obj)[keep_cells])
metadata <- obj@meta.data

message("Cells after metadata filtering: ", ncol(obj))
print(table(metadata$dataset, useNA = "ifany"))

message("RNA assay layers:")
print(get_layers_safe(obj, assay = "RNA"))

layer_selection <- bind_rows(
  pick_best_counts_layer_for_dataset(obj, assay = "RNA", dataset_value = "fetal"),
  pick_best_counts_layer_for_dataset(obj, assay = "RNA", dataset_value = "adult")
)

fetal_layer <- layer_selection$chosen_layer[layer_selection$dataset == "fetal"]
adult_layer <- layer_selection$chosen_layer[layer_selection$dataset == "adult"]

message("Chosen fetal counts layer: ", fetal_layer)
message("Chosen adult counts layer: ", adult_layer)

if (!grepl("^counts", fetal_layer)) {
  stop("Chosen fetal layer is not a counts layer: ", fetal_layer)
}

if (!grepl("^counts", adult_layer)) {
  stop("Chosen adult layer is not a counts layer: ", adult_layer)
}

if (grepl("^data|^scale|integrated", fetal_layer, ignore.case = TRUE)) {
  stop("Chosen fetal layer is data/scale/integrated-like, not counts: ", fetal_layer)
}

if (grepl("^data|^scale|integrated", adult_layer, ignore.case = TRUE)) {
  stop("Chosen adult layer is data/scale/integrated-like, not counts: ", adult_layer)
}

counts_fetal <- get_layer_matrix(obj, assay = "RNA", layer = fetal_layer)
counts_adult <- get_layer_matrix(obj, assay = "RNA", layer = adult_layer)

layer_qc <- bind_rows(
  summarize_sparse_matrix(counts_fetal, "fetal", fetal_layer),
  summarize_sparse_matrix(counts_adult, "adult", adult_layer)
)

write_csv(
  layer_selection,
  file.path(out_table_dir, "PSEUDOBULK_edgeR_count_layer_selection.csv")
)

write_csv(
  layer_qc,
  file.path(out_table_dir, "PSEUDOBULK_edgeR_selected_count_layer_qc.csv")
)

if (any(layer_qc$fraction_noninteger > 0.001)) {
  stop("Selected count layer contains non-integer values. Check count-layer selection.")
}
cells_fetal <- colnames(obj)[as.character(obj$dataset) == "fetal"]
cells_adult <- colnames(obj)[as.character(obj$dataset) == "adult"]

cells_fetal_in_layer <- intersect(cells_fetal, colnames(counts_fetal))
cells_adult_in_layer <- intersect(cells_adult, colnames(counts_adult))

if (length(cells_fetal_in_layer) < length(cells_fetal)) {
  message(
    "Fetal layer is missing ",
    length(setdiff(cells_fetal, cells_fetal_in_layer)),
    " fetal cells; dropping them from analysis."
  )
}

if (length(cells_adult_in_layer) < length(cells_adult)) {
  message(
    "Adult layer is missing ",
    length(setdiff(cells_adult, cells_adult_in_layer)),
    " adult cells; dropping them from analysis."
  )
}

keep_counts_cells <- c(cells_fetal_in_layer, cells_adult_in_layer)
obj <- subset(obj, cells = keep_counts_cells)
metadata <- obj@meta.data

counts_fetal <- counts_fetal[
  ,
  intersect(colnames(counts_fetal), colnames(obj)[metadata$dataset == "fetal"]),
  drop = FALSE
]

counts_adult <- counts_adult[
  ,
  intersect(colnames(counts_adult), colnames(obj)[metadata$dataset == "adult"]),
  drop = FALSE
]

genes_shared_for_de <- sort(intersect(rownames(counts_fetal), rownames(counts_adult)))

if (length(genes_shared_for_de) <= 2000) {
  stop("Too few shared genes for pseudobulk DE: ", length(genes_shared_for_de))
}

counts_fetal <- counts_fetal[genes_shared_for_de, , drop = FALSE]
counts_adult <- counts_adult[genes_shared_for_de, , drop = FALSE]

counts_gene_cell <- cbind(counts_fetal, counts_adult)
counts_gene_cell <- counts_gene_cell[, colnames(obj), drop = FALSE]

stopifnot(ncol(counts_gene_cell) == ncol(obj))
stopifnot(identical(colnames(counts_gene_cell), rownames(metadata)))

write_csv(
  tibble(gene = genes_shared_for_de),
  file.path(out_table_dir, "PSEUDOBULK_edgeR_shared_genes_used_for_DE.csv")
)

message("Combined count matrix:")
message("  genes: ", nrow(counts_gene_cell))
message("  cells: ", ncol(counts_gene_cell))

dataset_count_table <- metadata |>
  count(dataset, name = "n_cells") |>
  arrange(dataset)

write_csv(
  dataset_count_table,
  file.path(out_table_dir, "PSEUDOBULK_edgeR_cells_used_by_dataset.csv")
)

celltype_dataset_table <- metadata |>
  count(celltype_unified, dataset, name = "n_cells") |>
  tidyr::pivot_wider(
    names_from = dataset,
    values_from = n_cells,
    values_fill = 0
  ) |>
  arrange(celltype_unified)
write_csv(
  celltype_dataset_table,
  file.path(out_table_dir, "PSEUDOBULK_edgeR_celltype_counts_by_dataset.csv")
)

celltype_matrix <- table(metadata$celltype_unified, metadata$dataset)
shared_celltypes <- rownames(celltype_matrix)[
  celltype_matrix[, "fetal"] > 0 & celltype_matrix[, "adult"] > 0
]

shared_celltypes <- sort(shared_celltypes)

write_csv(
  tibble(celltype = shared_celltypes),
  file.path(out_table_dir, "PSEUDOBULK_edgeR_shared_celltypes.csv")
)

message("Cell types present in both fetal and adult:")
print(shared_celltypes)

qc_rows <- list()
pseudobulk_sample_rows <- list()

for (celltype in shared_celltypes) {
  message("\n============================================")
  message("Pseudobulk edgeR DE for cell type: ", celltype)

  min_cells_per_donor <- if (celltype == "germ") {
    min_cells_per_donor_germ
  } else {
    min_cells_per_donor_default
  }

  cells_celltype <- rownames(metadata)[metadata$celltype_unified == celltype]

  output_csv <- file.path(
    out_table_dir,
    paste0(
      "PSEUDOBULK_edgeR_DEG_fetal_vs_adult__",
      sanitize_filename(celltype),
      ".csv"
    )
  )

  if (length(cells_celltype) == 0) {
    write_empty_result(output_csv, celltype, "skip_no_cells")

    qc_rows[[length(qc_rows) + 1]] <- tibble(
      celltype = celltype,
      status = "skip_no_cells",
      min_cells_per_donor = min_cells_per_donor,
      adult_pseudobulks = NA_integer_,
      fetal_pseudobulks = NA_integer_,
      genes_tested = NA_integer_,
      output_file = output_csv
    )

    next
  }

  pseudobulk_id <- paste0(
    metadata[cells_celltype, "dataset"],
    "__",
    metadata[cells_celltype, "donor_unified"]
  )

  contribution <- table(pseudobulk_id)
  keep_pseudobulks <- names(contribution[contribution >= min_cells_per_donor])
  keep_cells_celltype <- cells_celltype[pseudobulk_id %in% keep_pseudobulks]

  if (length(keep_cells_celltype) == 0) {
    write_empty_result(output_csv, celltype, "skip_no_pseudobulk_meeting_min_cells")

    qc_rows[[length(qc_rows) + 1]] <- tibble(
      celltype = celltype,
      status = "skip_no_pseudobulk_meeting_min_cells",
      min_cells_per_donor = min_cells_per_donor,
      adult_pseudobulks = 0L,
      fetal_pseudobulks = 0L,
      genes_tested = NA_integer_,
      output_file = output_csv
    )

    next
  }
  pseudobulk_id_filtered <- paste0(
    metadata[keep_cells_celltype, "dataset"],
    "__",
    metadata[keep_cells_celltype, "donor_unified"]
  )

  sub_counts <- counts_gene_cell[, keep_cells_celltype, drop = FALSE]
  pseudobulk_counts <- pseudobulk_sum(sub_counts, pseudobulk_id_filtered)

  pseudobulk_samples <- colnames(pseudobulk_counts)
  pseudobulk_group <- ifelse(grepl("^fetal__", pseudobulk_samples), "fetal", "adult")
  pseudobulk_group <- factor(pseudobulk_group, levels = c("adult", "fetal"))

  group_n <- table(pseudobulk_group)

  adult_pseudobulks <- if ("adult" %in% names(group_n)) {
    as.integer(group_n[["adult"]])
  } else {
    0L
  }

  fetal_pseudobulks <- if ("fetal" %in% names(group_n)) {
    as.integer(group_n[["fetal"]])
  } else {
    0L
  }

  pseudobulk_sample_rows[[length(pseudobulk_sample_rows) + 1]] <- tibble(
    celltype = celltype,
    pseudobulk_sample = pseudobulk_samples,
    group = as.character(pseudobulk_group),
    donor_unified = sub("^[^_]+__", "", pseudobulk_samples),
    n_cells = as.integer(contribution[pseudobulk_samples]),
    min_cells_per_donor = min_cells_per_donor
  )

  message("Pseudobulk samples per group:")
  print(group_n)

  if (
    adult_pseudobulks < min_donors_per_group ||
      fetal_pseudobulks < min_donors_per_group
  ) {
    write_empty_result(output_csv, celltype, "skip_insufficient_pseudobulks_per_group")

    qc_rows[[length(qc_rows) + 1]] <- tibble(
      celltype = celltype,
      status = "skip_insufficient_pseudobulks_per_group",
      min_cells_per_donor = min_cells_per_donor,
      adult_pseudobulks = adult_pseudobulks,
      fetal_pseudobulks = fetal_pseudobulks,
      genes_tested = NA_integer_,
      output_file = output_csv
    )

    next
  }

  y <- DGEList(
    counts = as.matrix(pseudobulk_counts),
    group = pseudobulk_group
  )

  y <- calcNormFactors(y)

  keep_genes <- rowSums(cpm(y) >= min_cpm_filter) >= min_samples_filter

  if (sum(keep_genes) <= 10) {
    write_empty_result(output_csv, celltype, "skip_too_few_genes_after_cpm_filter")

    qc_rows[[length(qc_rows) + 1]] <- tibble(
      celltype = celltype,
      status = "skip_too_few_genes_after_cpm_filter",
      min_cells_per_donor = min_cells_per_donor,
      adult_pseudobulks = adult_pseudobulks,
      fetal_pseudobulks = fetal_pseudobulks,
      genes_tested = sum(keep_genes),
      output_file = output_csv
    )
    next
  }

  y <- y[keep_genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)

  design <- model.matrix(~ pseudobulk_group)

  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design)
  qlf <- glmQLFTest(fit, coef = "pseudobulk_groupfetal")

  result <- topTags(qlf, n = Inf)$table
  result$gene <- rownames(result)
  result$celltype <- celltype
  result$comparison <- "fetal_vs_adult"
  result$logFC_direction <- "positive_logFC_indicates_higher_in_fetal"

  result <- result |>
    relocate(celltype, comparison, gene, logFC_direction)

  write_csv(result, output_csv)

  message("Saved: ", output_csv)

  qc_rows[[length(qc_rows) + 1]] <- tibble(
    celltype = celltype,
    status = "ok_saved",
    min_cells_per_donor = min_cells_per_donor,
    adult_pseudobulks = adult_pseudobulks,
    fetal_pseudobulks = fetal_pseudobulks,
    genes_tested = nrow(result),
    output_file = output_csv
  )
}

qc_summary <- bind_rows(qc_rows)

write_csv(
  qc_summary,
  file.path(out_table_dir, "PSEUDOBULK_edgeR_QC_summary.csv")
)

if (length(pseudobulk_sample_rows) > 0) {
  write_csv(
    bind_rows(pseudobulk_sample_rows),
    file.path(out_table_dir, "PSEUDOBULK_edgeR_pseudobulk_sample_metadata.csv")
  )
}

summary_lines <- c(
  paste("Input integrated object:", integrated_rds),
  paste("Output directory:", results_root),
  paste("Cells after metadata/count-layer filtering:", ncol(obj)),
  paste("Genes used for DE:", length(genes_shared_for_de)),
  paste("Fetal count layer:", fetal_layer),
  paste("Adult count layer:", adult_layer),
  paste("Fetal count-layer coverage:", layer_selection$coverage[layer_selection$dataset == "fetal"]),
  paste("Adult count-layer coverage:", layer_selection$coverage[layer_selection$dataset == "adult"]),
  paste("Fetal count-layer fraction non-integer:", layer_qc$fraction_noninteger[layer_qc$dataset == "fetal"]),
  paste("Adult count-layer fraction non-integer:", layer_qc$fraction_noninteger[layer_qc$dataset == "adult"]),
  paste("Minimum cells per donor default:", min_cells_per_donor_default),
  paste("Minimum cells per donor germ:", min_cells_per_donor_germ),
  paste("Minimum donor pseudobulks per group:", min_donors_per_group),
  paste("Minimum CPM filter:", min_cpm_filter),
  paste("Minimum samples passing CPM filter:", min_samples_filter),
  paste("Shared cell types tested/skipped:", paste(shared_celltypes, collapse = ", ")),
  "",
  "QC status counts:",
  paste(capture.output(print(table(qc_summary$status, useNA = "ifany"))), collapse = "\n")
)

writeLines(
  summary_lines,
  file.path(out_log_dir, "PSEUDOBULK_edgeR_fetal_vs_adult_summary.txt")
)

sink(file.path(out_log_dir, "PSEUDOBULK_edgeR_fetal_vs_adult_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nPseudobulk edgeR complete.\n")
cat("Output directory:\n  ", results_root, "\n", sep = "")
cat("\nQC summary:\n")
print(qc_summary)
