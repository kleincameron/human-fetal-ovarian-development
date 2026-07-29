#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)
options(bitmapType = "cairo")

project_root <- "/home/liyan/liyan/Final/github_code_for_publication"
results_base <- "/home/liyan/liyan/Final/github_code_for_publication_results"

source(file.path(project_root, "config", "plotting.R"))
source(file.path(project_root, "config", "labels_colors.R"))

paths_local <- file.path(project_root, "config", "paths_local.R")
if (file.exists(paths_local)) {
  source(paths_local)
}

manifest_file <- if (exists("xenium_manifest", inherits = FALSE)) {
  xenium_manifest
} else {
  file.path(project_root, "config", "xenium_samples.csv")
}

xenium_annotated_rds <- if (exists("xenium_annotated_rds", inherits = FALSE)) {
  xenium_annotated_rds
} else {
  file.path(
    results_base,
    "xenium_annotated_object",
    "objects",
    "fetal_ovary_xenium_annotated.rds"
  )
}

results_root <- file.path(results_base, "xenium_qc")
out_figure_dir <- file.path(results_root, "figures")
out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

# -------------------------------------------------------------------------
# Helpers.
# -------------------------------------------------------------------------
normalize_include <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

pick_column <- function(df, candidates, required = TRUE) {
  hit <- candidates[candidates %in% colnames(df)]

  if (length(hit) > 0) {
    return(hit[[1]])
  }

  if (isTRUE(required)) {
    stop("None of these columns were found: ", paste(candidates, collapse = ", "))
  }

  NA_character_
}

numeric_or_na <- function(df, candidates) {
  col <- pick_column(df, candidates, required = FALSE)

  if (is.na(col)) {
    return(rep(NA_real_, nrow(df)))
  }

  suppressWarnings(as.numeric(df[[col]]))
}
save_plot_pdf_png <- function(plot, prefix, width, height, dpi = 600) {
  pdf_file <- paste0(prefix, ".pdf")
  png_file <- paste0(prefix, ".png")

  ggsave(
    pdf_file,
    plot,
    width = width,
    height = height,
    useDingbats = FALSE
  )

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      png_file,
      plot,
      width = width,
      height = height,
      dpi = dpi,
      device = ragg::agg_png,
      background = "white"
    )
  } else {
    ggsave(
      png_file,
      plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white"
    )
  }

  invisible(c(pdf = pdf_file, png = png_file))
}

get_counts_matrix <- function(obj, assay) {
  mat <- tryCatch(
    SeuratObject::LayerData(obj, assay = assay, layer = "counts"),
    error = function(e) NULL
  )

  if (is.null(mat)) {
    mat <- tryCatch(
      GetAssayData(obj, assay = assay, slot = "counts"),
      error = function(e) NULL
    )
  }

  if (is.null(mat)) {
    stop("Could not retrieve counts from assay: ", assay)
  }

  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }

  mat
}

summarize_numeric_metric <- function(data, metric) {
  values <- suppressWarnings(as.numeric(data[[metric]]))
  values <- values[is.finite(values)]

  if (length(values) == 0) {
    return(tibble(
      metric = metric,
      n = 0L,
      min = NA_real_,
      q05 = NA_real_,
      median = NA_real_,
      mean = NA_real_,
      q95 = NA_real_,
      max = NA_real_
    ))
  }

  tibble(
    metric = metric,
    n = length(values),
    min = min(values),
    q05 = unname(quantile(values, 0.05)),
    median = median(values),
    mean = mean(values),
    q95 = unname(quantile(values, 0.95)),
    max = max(values)
  )
}

histogram_plot <- function(data, x, x_label, bins = 60) {
  ggplot(data, aes(x = .data[[x]])) +
    geom_histogram(bins = bins) +
    labs(x = x_label, y = "Cells") +
    theme_publication()
}

safe_annotation_label <- function(final_annotation, fallback_label) {
  out <- as.character(final_annotation)
  fallback_label <- as.character(fallback_label)

  bad <- is.na(out) | out == "" | fallback_label %in% c("Ambiguous_subcluster", "Ambiguous")
  out[bad] <- fallback_label[bad]
  out[is.na(out) | out == ""] <- "Ambiguous_subcluster"

  out
}
# -------------------------------------------------------------------------
# Inputs.
# -------------------------------------------------------------------------
if (!file.exists(manifest_file)) {
  stop(
    "Missing local Xenium manifest:\n  ",
    manifest_file,
    "\nCreate it from config/xenium_samples_example.csv and set xenium_dir."
  )
}

manifest <- read_csv(manifest_file, show_col_types = FALSE)

required_manifest_cols <- c("sample_id", "gestational_week", "xenium_dir", "include")
missing_manifest_cols <- setdiff(required_manifest_cols, colnames(manifest))

if (length(missing_manifest_cols) > 0) {
  stop("Manifest missing columns: ", paste(missing_manifest_cols, collapse = ", "))
}

manifest <- manifest |>
  mutate(
    include = normalize_include(include),
    sample_id = as.character(sample_id),
    gestational_week = as.character(gestational_week),
    xenium_dir = as.character(xenium_dir)
  ) |>
  filter(include)

if (nrow(manifest) != 1) {
  stop("Expected exactly one included Xenium sample; found: ", nrow(manifest))
}

sample_id <- manifest$sample_id[[1]]
gestational_week <- manifest$gestational_week[[1]]
xenium_dir <- manifest$xenium_dir[[1]]
cells_csv <- file.path(xenium_dir, "cells.csv.gz")

stopifnot(file.exists(xenium_annotated_rds))
stopifnot(file.exists(cells_csv))

cat("Loading annotated Xenium object:\n  ", xenium_annotated_rds, "\n", sep = "")
xenium <- readRDS(xenium_annotated_rds)
stopifnot(inherits(xenium, "Seurat"))

cat("Loading Xenium cell summary:\n  ", cells_csv, "\n", sep = "")
cells_raw <- read_csv(cells_csv, show_col_types = FALSE)

cell_id_col <- pick_column(cells_raw, c("cell_id", "CellID", "barcode"))

cells_qc <- cells_raw |>
  transmute(
    cell_id = as.character(.data[[cell_id_col]]),
    transcript_counts_cell_summary = numeric_or_na(cells_raw, c("transcript_counts")),
    total_counts_cell_summary = numeric_or_na(cells_raw, c("total_counts")),
    control_probe_counts = numeric_or_na(cells_raw, c("control_probe_counts")),
    control_codeword_counts = numeric_or_na(cells_raw, c("control_codeword_counts")),
    genomic_control_counts = numeric_or_na(cells_raw, c("genomic_control_counts")),
    unassigned_codeword_counts = numeric_or_na(cells_raw, c("unassigned_codeword_counts")),
    deprecated_codeword_counts = numeric_or_na(cells_raw, c("deprecated_codeword_counts")),
    cell_area = numeric_or_na(cells_raw, c("cell_area", "cell_area_um2")),
    nucleus_area = numeric_or_na(cells_raw, c("nucleus_area", "nucleus_area_um2")),
    nucleus_count = numeric_or_na(cells_raw, c("nucleus_count"))
  )

xen_assay <- if ("Xenium" %in% Assays(xenium)) {
  "Xenium"
} else {
  DefaultAssay(xenium)
}

DefaultAssay(xenium) <- xen_assay

counts <- get_counts_matrix(xenium, assay = xen_assay)

matrix_qc <- tibble(
  cell_id = colnames(counts),
  xenium_matrix_transcripts = as.numeric(Matrix::colSums(counts)),
  xenium_matrix_genes = as.numeric(Matrix::colSums(counts > 0))
)

cell_summary_join_rate <- mean(matrix_qc$cell_id %in% cells_qc$cell_id)

if (!is.finite(cell_summary_join_rate) || cell_summary_join_rate < 0.90) {
  stop(
    "Cell summary join rate is too low: ",
    round(100 * cell_summary_join_rate, 2),
    "%. Check whether cells.csv.gz cell IDs match the Seurat object cell IDs."
  )
}
metadata_qc <- xenium@meta.data |>
  rownames_to_column("cell_id") |>
  transmute(
    cell_id = cell_id,
    sample_id = as.character(sample_id),
    gestational_week = as.character(gestational_week),
    dataset = as.character(dataset),
    pred_major_cell_type_filtered = as.character(pred_major_cell_type_filtered),
    pred_major_cell_type_score = suppressWarnings(as.numeric(pred_major_cell_type_score)),
    pred_major_cell_type_gap = suppressWarnings(as.numeric(pred_major_cell_type_gap)),
    pred_fine_subcluster_2stage_filtered = as.character(pred_fine_subcluster_2stage_filtered),
    pred_final_annotation_2stage_filtered = as.character(pred_final_annotation_2stage_filtered),
    pred_fine_subcluster_2stage_score = suppressWarnings(as.numeric(pred_fine_subcluster_2stage_score)),
    pred_fine_subcluster_hybrid = as.character(pred_fine_subcluster_hybrid),
    pred_final_annotation_hybrid = as.character(pred_final_annotation_hybrid)
  )

qc <- matrix_qc |>
  left_join(cells_qc, by = "cell_id") |>
  left_join(metadata_qc, by = "cell_id") |>
  mutate(
    transcripts_per_cell = ifelse(
      is.finite(transcript_counts_cell_summary),
      transcript_counts_cell_summary,
      xenium_matrix_transcripts
    ),
    genes_per_cell = xenium_matrix_genes,
    negative_control_counts = rowSums(
      cbind(
        control_probe_counts,
        control_codeword_counts,
        genomic_control_counts,
        unassigned_codeword_counts,
        deprecated_codeword_counts
      ),
      na.rm = TRUE
    ),
    total_counts_for_fraction = ifelse(
      is.finite(total_counts_cell_summary),
      total_counts_cell_summary,
      transcripts_per_cell + negative_control_counts
    ),
    negative_control_fraction = ifelse(
      total_counts_for_fraction > 0,
      negative_control_counts / total_counts_for_fraction,
      NA_real_
    ),
    nucleus_to_cell_area_ratio = ifelse(
      cell_area > 0,
      nucleus_area / cell_area,
      NA_real_
    ),
    empty_cell = transcripts_per_cell == 0,
    log10_transcripts_per_cell = log10(transcripts_per_cell + 1),
    log10_genes_per_cell = log10(genes_per_cell + 1),
    log10_cell_area = log10(cell_area + 1),
    log10_nucleus_area = log10(nucleus_area + 1),
    percent_negative_control = 100 * negative_control_fraction,
    is_major_ambiguous = pred_major_cell_type_filtered == "Ambiguous",
    is_stage2_ambiguous = pred_fine_subcluster_2stage_filtered == "Ambiguous_subcluster",
    has_hybrid_annotation = !is.na(pred_fine_subcluster_hybrid) &
      pred_fine_subcluster_hybrid != "" &
      pred_fine_subcluster_hybrid != "Ambiguous_subcluster",
    stage2_annotation_label = safe_annotation_label(
      pred_final_annotation_2stage_filtered,
      pred_fine_subcluster_2stage_filtered
    ),
    hybrid_annotation_label = safe_annotation_label(
      pred_final_annotation_hybrid,
      pred_fine_subcluster_hybrid
    )
  )
cat("Cells in annotated object: ", ncol(xenium), "\n", sep = "")
cat("Cells joined to cells.csv.gz: ", sum(matrix_qc$cell_id %in% cells_qc$cell_id), "\n", sep = "")

# -------------------------------------------------------------------------
# Tables.
# -------------------------------------------------------------------------
raw_metrics <- c(
  "transcripts_per_cell",
  "genes_per_cell",
  "negative_control_counts",
  "negative_control_fraction",
  "cell_area",
  "nucleus_area",
  "nucleus_to_cell_area_ratio",
  "nucleus_count"
)

raw_qc_summary <- bind_rows(
  lapply(raw_metrics, function(metric) summarize_numeric_metric(qc, metric))
)

raw_qc_summary <- bind_rows(
  raw_qc_summary,
  tibble(
    metric = "empty_cell_fraction",
    n = nrow(qc),
    min = NA_real_,
    q05 = NA_real_,
    median = mean(qc$empty_cell, na.rm = TRUE),
    mean = mean(qc$empty_cell, na.rm = TRUE),
    q95 = NA_real_,
    max = NA_real_
  )
)

write_csv(
  raw_qc_summary,
  file.path(out_table_dir, "xenium_raw_qc_summary.csv")
)

annotation_qc_summary <- tibble(
  metric = c(
    "cells_total",
    "major_ambiguous_cells",
    "major_ambiguous_fraction",
    "stage2_ambiguous_cells",
    "stage2_ambiguous_fraction",
    "hybrid_annotated_cells",
    "hybrid_annotated_fraction",
    "median_major_prediction_score",
    "median_major_prediction_score_gap",
    "median_stage2_prediction_score"
  ),
  value = c(
    nrow(qc),
    sum(qc$is_major_ambiguous, na.rm = TRUE),
    mean(qc$is_major_ambiguous, na.rm = TRUE),
    sum(qc$is_stage2_ambiguous, na.rm = TRUE),
    mean(qc$is_stage2_ambiguous, na.rm = TRUE),
    sum(qc$has_hybrid_annotation, na.rm = TRUE),
    mean(qc$has_hybrid_annotation, na.rm = TRUE),
    median(qc$pred_major_cell_type_score, na.rm = TRUE),
    median(qc$pred_major_cell_type_gap, na.rm = TRUE),
    median(qc$pred_fine_subcluster_2stage_score, na.rm = TRUE)
  )
)

write_csv(
  annotation_qc_summary,
  file.path(out_table_dir, "xenium_annotation_qc_summary.csv")
)

annotation_counts <- bind_rows(
  qc |>
    count(pred_major_cell_type_filtered, sort = TRUE, name = "cells") |>
    transmute(
      annotation_level = "major_cell_type_filtered",
      label = pred_major_cell_type_filtered,
      cells = cells
    ),
  qc |>
    count(stage2_annotation_label, sort = TRUE, name = "cells") |>
    transmute(
      annotation_level = "fine_subcluster_2stage_final_annotation",
      label = stage2_annotation_label,
      cells = cells
    ),
  qc |>
    count(hybrid_annotation_label, sort = TRUE, name = "cells") |>
    transmute(
      annotation_level = "fine_subcluster_hybrid_final_annotation",
      label = hybrid_annotation_label,
      cells = cells
    )
)

write_csv(
  annotation_counts,
  file.path(out_table_dir, "xenium_annotation_counts.csv")
)

# -------------------------------------------------------------------------
# Figures: raw-data QC.
# -------------------------------------------------------------------------
p_raw_distributions <- (
  histogram_plot(qc, "log10_transcripts_per_cell", "log10(transcripts per cell + 1)") |
    histogram_plot(qc, "log10_genes_per_cell", "log10(genes detected per cell + 1)")
) / (
  histogram_plot(qc, "log10_cell_area", "log10(cell area + 1)") |
    histogram_plot(qc, "log10_nucleus_area", "log10(nucleus area + 1)")
) / (
  histogram_plot(qc, "percent_negative_control", "Negative/control counts (%)") |
    histogram_plot(qc, "nucleus_to_cell_area_ratio", "Nucleus/cell area ratio")
)

save_plot_pdf_png(
  p_raw_distributions,
  file.path(out_figure_dir, "xenium_raw_cell_qc_distributions"),
  width = 9,
  height = 8
)

# -------------------------------------------------------------------------
# Figures: annotation QC.
# -------------------------------------------------------------------------
p_major_score <- histogram_plot(
  qc,
  "pred_major_cell_type_score",
  "Major cell-type prediction score"
)

p_major_gap <- histogram_plot(
  qc,
  "pred_major_cell_type_gap",
  "Major cell-type score gap"
)

p_stage2_score <- histogram_plot(
  qc,
  "pred_fine_subcluster_2stage_score",
  "Stage-2 subcluster prediction score"
)
annotation_bar_data <- qc |>
  count(hybrid_annotation_label, sort = TRUE, name = "cells") |>
  slice_head(n = 30) |>
  mutate(hybrid_annotation_label = reorder(hybrid_annotation_label, cells))

p_annotation_counts <- ggplot(
  annotation_bar_data,
  aes(x = hybrid_annotation_label, y = cells)
) +
  geom_col() +
  coord_flip() +
  labs(
    x = NULL,
    y = "Cells"
  ) +
  theme_publication()

p_annotation_distributions <- (
  p_major_score | p_major_gap
) / (
  p_stage2_score | p_annotation_counts
)

save_plot_pdf_png(
  p_annotation_distributions,
  file.path(out_figure_dir, "xenium_annotation_qc_distributions"),
  width = 10,
  height = 8
)

# -------------------------------------------------------------------------
# Logs.
# -------------------------------------------------------------------------
summary_lines <- c(
  paste("Input Xenium annotated object:", xenium_annotated_rds),
  paste("Input Xenium manifest:", manifest_file),
  paste("Input cells.csv.gz:", cells_csv),
  paste("Sample ID:", sample_id),
  paste("Gestational week:", gestational_week),
  paste("Xenium assay:", xen_assay),
  paste("Cells:", nrow(qc)),
  paste("Genes/features in Xenium assay:", nrow(counts)),
  paste("Cell summary join rate:", round(100 * cell_summary_join_rate, 3), "%"),
  "",
  "Figures:",
  paste("  ", file.path(out_figure_dir, "xenium_raw_cell_qc_distributions.pdf")),
  paste("  ", file.path(out_figure_dir, "xenium_annotation_qc_distributions.pdf")),
  "",
  "Tables:",
  paste("  ", file.path(out_table_dir, "xenium_raw_qc_summary.csv")),
  paste("  ", file.path(out_table_dir, "xenium_annotation_qc_summary.csv")),
  paste("  ", file.path(out_table_dir, "xenium_annotation_counts.csv")),
  "",
  "Annotation QC uses final annotation names for fine-subcluster-derived labels rather than machine labels."
)

writeLines(
  summary_lines,
  file.path(out_log_dir, "xenium_qc_summary.txt")
)

sink(file.path(out_log_dir, "xenium_qc_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone.\n")
cat("Xenium QC outputs written to:\n  ", results_root, "\n", sep = "")
