# Prepare the adult ovary CELLxGENE reference as a Seurat object.
#
# The selected adult H5AD contains:
#   - X: normalized/log-like expression
#   - layers["decontXcounts"]: integer, count-like expression
#
# For downstream compatibility:
#   - RNA/counts stores decontXcounts
#   - RNA/data stores X

options(bitmapType = "cairo")

# Use the active conda/environment Python rather than reticulate's managed
# Python environment. This avoids network-dependent uv setup on clusters.
if (!nzchar(Sys.getenv("RETICULATE_PYTHON"))) {
  active_python <- Sys.which("python")
  if (nzchar(active_python)) {
    Sys.setenv(RETICULATE_PYTHON = active_python)
  }
}
Sys.setenv(RETICULATE_USE_MANAGED_VENV = "no")

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(reticulate)
  library(dplyr)
  library(readr)
  library(tibble)
})

project_root <- "/home/liyan/liyan/Final/github_code_for_publication"
adult_results_root <- "/home/liyan/liyan/Final/github_code_for_publication_results/adult_cellxgene_reference"

paths_local <- file.path(project_root, "config", "paths_local.R")
if (file.exists(paths_local)) {
  source(paths_local)
}

manifest_csv <- file.path(project_root, "metadata", "adult_reference_dataset_manifest.csv")
stopifnot(file.exists(manifest_csv))

manifest <- read_csv(manifest_csv, show_col_types = FALSE)

adult_dataset_id <- "584027d5-32d7-4696-9424-f61134ff2aa7"
adult_dataset <- manifest |>
  filter(dataset_id == adult_dataset_id)

if (nrow(adult_dataset) != 1) {
  stop("Adult dataset manifest must contain exactly one row for dataset_id: ", adult_dataset_id)
}
raw_dir <- file.path(adult_results_root, "raw")
object_dir <- file.path(adult_results_root, "objects")
table_dir <- file.path(adult_results_root, "tables")
log_dir <- file.path(adult_results_root, "logs")

dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

default_adult_h5ad <- file.path(raw_dir, adult_dataset$expected_h5ad_filename)

adult_h5ad <- if (exists("adult_cellxgene_h5ad", inherits = FALSE)) {
  adult_cellxgene_h5ad
} else {
  default_adult_h5ad
}

adult_rds <- file.path(
  object_dir,
  "adult_cellxgene_584027d5_seurat.rds"
)

feature_map_csv <- file.path(
  table_dir,
  "adult_cellxgene_584027d5_feature_map.csv"
)

cell_metadata_columns_csv <- file.path(
  table_dir,
  "adult_cellxgene_584027d5_cell_metadata_columns.csv"
)

matrix_summary_csv <- file.path(
  table_dir,
  "adult_cellxgene_584027d5_matrix_summary.csv"
)

summary_txt <- file.path(
  log_dir,
  "adult_cellxgene_reference_summary.txt"
)

session_txt <- file.path(
  log_dir,
  "adult_cellxgene_reference_sessionInfo.txt"
)

pick_first_existing <- function(df, candidates) {
  available <- colnames(df)
  if (is.null(available)) return(NA_character_)

  hit <- candidates[candidates %in% available]
  if (length(hit) == 0) return(NA_character_)

  hit[[1]]
}
as_dgC_from_scipy_sparse <- function(x) {
  scipy_sparse <- reticulate::import("scipy.sparse", convert = FALSE)

  if (!reticulate::py_to_r(scipy_sparse$issparse(x))) {
    dense <- reticulate::py_to_r(x)
    dense <- as.matrix(dense)
    return(as(dense, "dgCMatrix"))
  }

  x_csc <- x$tocsc()

  values <- as.numeric(reticulate::py_to_r(x_csc$data))
  row_indices <- as.integer(reticulate::py_to_r(x_csc$indices))
  column_pointers <- as.integer(reticulate::py_to_r(x_csc$indptr))
  matrix_shape <- as.integer(reticulate::py_to_r(x_csc$shape))

  storage.mode(values) <- "double"

  new(
    "dgCMatrix",
    x = values,
    i = row_indices,
    p = column_pointers,
    Dim = matrix_shape
  )
}

summarize_sparse_matrix <- function(matrix, source_name) {
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
    source = source_name,
    genes = nrow(matrix),
    cells = ncol(matrix),
    nonzero_values = length(values),
    fraction_noninteger = fraction_noninteger,
    min_value = value_range[[1]],
    max_value = value_range[[2]]
  )
}

download_adult_h5ad_if_missing <- function(dataset_id, output_file) {
  if (file.exists(output_file)) {
    cat("Using existing adult H5AD:\n  ", output_file, "\n", sep = "")
    return(invisible(output_file))
  }
  cat("Adult H5AD not found. Attempting CELLxGENE Census download.\n")

  if (!requireNamespace("cellxgene.census", quietly = TRUE)) {
    stop(
      "The adult H5AD is missing and the cellxgene.census R package is not installed.\n",
      "Either install cellxgene.census in the R environment or manually download the H5AD to:\n",
      output_file,
      "\nAlternatively, define adult_cellxgene_h5ad in config/paths_local.R."
    )
  }

  dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

  cellxgene.census::download_source_h5ad(
    dataset_id = dataset_id,
    file = output_file,
    overwrite = FALSE
  )

  if (!file.exists(output_file)) {
    stop("CELLxGENE H5AD download did not create expected file: ", output_file)
  }

  invisible(output_file)
}

download_adult_h5ad_if_missing(
  dataset_id = adult_dataset_id,
  output_file = adult_h5ad
)

stopifnot(file.exists(adult_h5ad))

cat("Python configuration:\n")
print(reticulate::py_config())

anndata <- reticulate::import("anndata", convert = FALSE)

cat("\nReading adult CELLxGENE H5AD:\n  ", adult_h5ad, "\n", sep = "")
adata <- anndata$read_h5ad(adult_h5ad)

py_builtins <- reticulate::import_builtins(convert = FALSE)

layer_names <- as.character(
  reticulate::py_to_r(
    py_builtins$list(adata$layers$keys())
  )
)

data_source <- "X"
counts_source <- if ("decontXcounts" %in% layer_names) {
  "layers/decontXcounts"
} else {
  NA_character_
}

if (is.na(counts_source)) {
  stop(
    "No count-like adult matrix was found. Expected layers['decontXcounts']; ",
    "available layers are: ",
    paste(layer_names, collapse = ", ")
  )
}
genes_raw <- as.character(
  reticulate::py_to_r(adata$var_names$to_list())
)

cells <- as.character(
  reticulate::py_to_r(adata$obs_names$to_list())
)

cat(
  "AnnData dimensions: ",
  length(cells),
  " cells x ",
  length(genes_raw),
  " genes\n",
  sep = ""
)

cat("AnnData layers: ", paste(layer_names, collapse = ", "), "\n", sep = "")
cat("Using RNA/counts source: ", counts_source, "\n", sep = "")
cat("Using RNA/data source: ", data_source, "\n", sep = "")

data_cells_by_genes <- as_dgC_from_scipy_sparse(adata$X)
counts_cells_by_genes <- as_dgC_from_scipy_sparse(adata$layers[["decontXcounts"]])

stopifnot(nrow(data_cells_by_genes) == length(cells))
stopifnot(ncol(data_cells_by_genes) == length(genes_raw))
stopifnot(nrow(counts_cells_by_genes) == length(cells))
stopifnot(ncol(counts_cells_by_genes) == length(genes_raw))

dimnames(data_cells_by_genes) <- list(cells, genes_raw)
dimnames(counts_cells_by_genes) <- list(cells, genes_raw)

data_mat <- Matrix::t(data_cells_by_genes)
counts_mat <- Matrix::t(counts_cells_by_genes)

stopifnot(identical(rownames(data_mat), rownames(counts_mat)))
stopifnot(identical(colnames(data_mat), colnames(counts_mat)))

var <- reticulate::py_to_r(adata$var)
if (!is.data.frame(var)) {
  var <- as.data.frame(var)
}
rownames(var) <- genes_raw

obs <- reticulate::py_to_r(adata$obs)
if (!is.data.frame(obs)) {
  obs <- as.data.frame(obs)
}
rownames(obs) <- cells

gene_symbol_candidates <- c(
  "feature_name",
  "feature_symbol",
  "gene_symbol",
  "symbol"
)

symbol_column <- pick_first_existing(
  var,
  gene_symbol_candidates
)
gene_symbol <- rep(NA_character_, nrow(var))

if (!is.na(symbol_column)) {
  gene_symbol <- as.character(var[[symbol_column]])
}

features_used <- genes_raw

if (!all(is.na(gene_symbol))) {
  features_used <- ifelse(
    is.na(gene_symbol) | gene_symbol == "",
    genes_raw,
    gene_symbol
  )
}

features_used <- gsub("_", "-", features_used)
features_used <- make.unique(features_used)

rownames(data_mat) <- features_used
rownames(counts_mat) <- features_used

matrix_summary <- bind_rows(
  summarize_sparse_matrix(counts_mat, counts_source),
  summarize_sparse_matrix(data_mat, data_source)
)

write_csv(matrix_summary, matrix_summary_csv)

counts_fraction_noninteger <- matrix_summary$fraction_noninteger[
  matrix_summary$source == counts_source
]

data_fraction_noninteger <- matrix_summary$fraction_noninteger[
  matrix_summary$source == data_source
]

counts_value_range <- c(
  matrix_summary$min_value[matrix_summary$source == counts_source],
  matrix_summary$max_value[matrix_summary$source == counts_source]
)

data_value_range <- c(
  matrix_summary$min_value[matrix_summary$source == data_source],
  matrix_summary$max_value[matrix_summary$source == data_source]
)

if (counts_fraction_noninteger > 0.001) {
  stop("Selected adult counts matrix is not integer/count-like.")
}

if (data_fraction_noninteger < 0.001) {
  warning("Selected adult data matrix appears integer-like; expected normalized/log-like X.")
}

cat(
  "\nAdult matrix summary:\n",
  "  counts source = ", counts_source,
  " | fraction non-integer = ", signif(counts_fraction_noninteger, 4),
  " | range = [", signif(counts_value_range[[1]], 4), ", ", signif(counts_value_range[[2]], 4), "]\n",
  "  data source = ", data_source,
  " | fraction non-integer = ", signif(data_fraction_noninteger, 4),
  " | range = [", signif(data_value_range[[1]], 4), ", ", signif(data_value_range[[2]], 4), "]\n",
  sep = ""
)

adult <- CreateSeuratObject(
  counts = counts_mat,
  project = "adult_cellxgene_584027d5"
)

data_mat <- data_mat[rownames(adult), colnames(adult), drop = FALSE]

adult <- SetAssayData(
  adult,
  assay = "RNA",
  layer = "data",
  new.data = data_mat
)
obs <- obs[colnames(adult), , drop = FALSE]
adult <- AddMetaData(adult, metadata = obs)

feature_map <- data.frame(
  feature_id = genes_raw,
  gene_symbol = gene_symbol,
  rowname_used = features_used,
  stringsAsFactors = FALSE
)

adult@misc$feature_map <- feature_map
adult@misc$cellxgene_manifest <- as.data.frame(adult_dataset)
adult@misc$matrix_source_counts <- counts_source
adult@misc$matrix_source_data <- data_source
adult@misc$note <- "RNA/counts stores layers['decontXcounts']; RNA/data stores adata$X normalized/log-like expression."

adult$dataset <- "adult"
adult$developmental_stage <- "adult"
adult$adult_reference_dataset <- "cellxgene_584027d5"
adult$adult_dataset_id <- adult_dataset_id
adult$source_database <- "CZ CELLxGENE Discover / CELLXGENE Census"
adult$source_counts <- counts_source
adult$source_data <- data_source
adult$counts_frac_noninteger <- counts_fraction_noninteger
adult$data_frac_noninteger <- data_fraction_noninteger
adult$X_frac_nonint <- data_fraction_noninteger

write_csv(feature_map, feature_map_csv)

write_csv(
  tibble(
    metadata_column = colnames(obs),
    class = vapply(obs, function(x) paste(class(x), collapse = ";"), character(1)),
    n_unique = vapply(obs, function(x) length(unique(x)), integer(1))
  ),
  cell_metadata_columns_csv
)

saveRDS(adult, adult_rds)

donor_summary <- if ("donor_id" %in% colnames(adult@meta.data)) {
  paste(capture.output(print(table(adult$donor_id, useNA = "ifany"))), collapse = "\n")
} else {
  "donor_id column not present"
}

cell_type_summary <- if ("cell_type" %in% colnames(adult@meta.data)) {
  paste(capture.output(print(sort(table(adult$cell_type, useNA = "ifany"), decreasing = TRUE))), collapse = "\n")
} else {
  "cell_type column not present"
}

rna_layers <- paste(Layers(adult[["RNA"]]), collapse = ", ")

summary_lines <- c(
  paste("Adult dataset name:", adult_dataset$dataset_name),
  paste("Source database:", adult_dataset$source_database),
  paste("Collection:", adult_dataset$collection_name),
  paste("Collection ID:", adult_dataset$collection_id),
  paste("Dataset ID:", adult_dataset$dataset_id),
  paste("Publication DOI:", adult_dataset$publication_doi),
  paste("Input H5AD:", adult_h5ad),
  paste("Output Seurat object:", adult_rds),
  paste("Cells:", ncol(adult)),
  paste("Genes:", nrow(adult)),
  paste("RNA/counts source:", counts_source),
  paste("RNA/data source:", data_source),
  paste("Counts fraction non-integer matrix values:", signif(counts_fraction_noninteger, 4)),
  paste("Data fraction non-integer matrix values:", signif(data_fraction_noninteger, 4)),
  paste("Counts value range:", paste(signif(counts_value_range, 4), collapse = " to ")),
  paste("Data value range:", paste(signif(data_value_range, 4), collapse = " to ")),
  paste("Gene symbol column:", ifelse(is.na(symbol_column), "none", symbol_column)),
  paste("RNA layers:", rna_layers),
  "",
  "Adult donors:",
  donor_summary,
  "",
  "Adult cell_type labels:",
  cell_type_summary
)
writeLines(summary_lines, summary_txt)

sink(session_txt)
print(sessionInfo())
sink()

cat("\nAdult cells:", ncol(adult), "\n")
cat("Adult genes:", nrow(adult), "\n")
cat("RNA layers:", rna_layers, "\n")
cat("\nSaved adult Seurat object:\n  ", adult_rds, "\n", sep = "")
cat("Summary written to:\n  ", summary_txt, "\n", sep = "")
