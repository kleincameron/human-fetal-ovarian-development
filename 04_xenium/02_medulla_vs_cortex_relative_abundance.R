#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(ggplot2)
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


source(file.path(project_root, "config", "plotting.R"))
source(file.path(project_root, "config", "labels_colors.R"))


manifest_file <- if (exists("xenium_manifest", inherits = FALSE)) {
  xenium_manifest
} else {
  file.path(project_root, "config", "xenium_samples.csv")
}

annotation_csv <- if (exists("xenium_annotation_csv", inherits = FALSE)) {
  xenium_annotation_csv
} else {
  file.path(
    results_base,
    "xenium_annotated_object",
    "tables",
    "xenium_cell_annotations.csv"
  )
}

snrna_annotation_metadata <- file.path(
  project_root,
  "metadata",
  "snRNAseq_cell_annotations.csv"
)

medulla_roi_csv <- file.path(
  project_root,
  "metadata",
  "xenium_medulla_roi_coordinates.csv"
)

results_root <- file.path(results_base, "xenium_medulla_vs_cortex_relative_abundance")

if (dir.exists(results_root)) {
  unlink(results_root, recursive = TRUE, force = TRUE)
}

out_figure_dir <- file.path(results_root, "figures")
out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

min_cells_per_annotation <- 50
pseudocount_percent <- 0.01
max_spatial_plot_cells <- 100000

fallback_celltype_order <- c(
  "germ",
  "degenerated",
  "granulosa",
  "stroma",
  "endothelial",
  "mural",
  "immune",
  "erythroid"
)

fallback_celltype_colors <- c(
  germ = "#1B9E77",
  degenerated = "#CFCFCF",
  granulosa = "#D95F02",
  stroma = "#7570B3",
  endothelial = "#E7298A",
  mural = "#66A61E",
  immune = "#E6AB02",
  erythroid = "#A6761D"
)
if (!exists("celltype_order", inherits = FALSE)) {
  celltype_order <- fallback_celltype_order
}

if (!exists("celltype_colors", inherits = FALSE)) {
  celltype_colors <- fallback_celltype_colors
} else {
  missing_colors <- setdiff(names(fallback_celltype_colors), names(celltype_colors))
  if (length(missing_colors) > 0) {
    celltype_colors <- c(celltype_colors, fallback_celltype_colors[missing_colors])
  }
}

celltype_order <- intersect(celltype_order, names(celltype_colors))

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

read_roi_coordinates <- function(path) {
  stopifnot(file.exists(path))

  candidate_skips <- c(0, 1, 2, 3, 4)

  for (skip_n in candidate_skips) {
    dat <- tryCatch(
      read_csv(path, skip = skip_n, show_col_types = FALSE),
      error = function(e) NULL
    )

    if (is.null(dat) || ncol(dat) == 0) {
      next
    }

    x_col <- pick_column(
      dat,
      c("x", "X", "x_um", "X_um", "x_centroid", "x_centroid_um"),
      required = FALSE
    )

    y_col <- pick_column(
      dat,
      c("y", "Y", "y_um", "Y_um", "y_centroid", "y_centroid_um"),
      required = FALSE
    )

    if (!is.na(x_col) && !is.na(y_col)) {
      out <- tibble(
        x = suppressWarnings(as.numeric(dat[[x_col]])),
        y = suppressWarnings(as.numeric(dat[[y_col]]))
      ) |>
        filter(is.finite(x), is.finite(y))

      if (nrow(out) >= 3) {
        return(out)
      }
    }
  }
  stop("Could not read medulla ROI coordinates from: ", path)
}

point_in_polygon <- function(px, py, poly_x, poly_y) {
  n <- length(poly_x)
  inside <- rep(FALSE, length(px))
  j <- n

  for (i in seq_len(n)) {
    xi <- poly_x[i]
    yi <- poly_y[i]
    xj <- poly_x[j]
    yj <- poly_y[j]

    crosses <- ((yi > py) != (yj > py)) &
      (px < (xj - xi) * (py - yi) / (yj - yi + .Machine$double.eps) + xi)

    inside <- xor(inside, crosses)
    j <- i
  }

  inside
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

reorder_within <- function(x, by, within, fun = mean, sep = "___") {
  stats::reorder(paste(x, within, sep = sep), by, FUN = fun)
}

scale_y_reordered <- function(..., sep = "___") {
  scale_y_discrete(labels = function(x) gsub(paste0(sep, ".*$"), "", x), ...)
}

# -------------------------------------------------------------------------
# Inputs.
# -------------------------------------------------------------------------
stopifnot(file.exists(manifest_file))
stopifnot(file.exists(annotation_csv))
stopifnot(file.exists(snrna_annotation_metadata))
stopifnot(file.exists(medulla_roi_csv))

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

stopifnot(file.exists(cells_csv))

cat("Loading Xenium annotation table:\n  ", annotation_csv, "\n", sep = "")
anno <- read_csv(annotation_csv, show_col_types = FALSE)

required_annotation_cols <- c(
  "cell_id",
  "pred_fine_subcluster_hybrid",
  "pred_final_annotation_hybrid"
)

missing_annotation_cols <- setdiff(required_annotation_cols, colnames(anno))
if (length(missing_annotation_cols) > 0) {
  stop("Annotation table missing columns: ", paste(missing_annotation_cols, collapse = ", "))
}

cat("Loading snRNA-seq annotation metadata:\n  ", snrna_annotation_metadata, "\n", sep = "")
snrna_md <- read_csv(snrna_annotation_metadata, show_col_types = FALSE)

required_snrna_cols <- c("fine_subcluster", "major_cell_type", "final_annotation")
missing_snrna_cols <- setdiff(required_snrna_cols, colnames(snrna_md))
if (length(missing_snrna_cols) > 0) {
  stop("snRNA-seq metadata missing columns: ", paste(missing_snrna_cols, collapse = ", "))
}

subcluster_map <- snrna_md |>
  transmute(
    fine_subcluster = as.character(fine_subcluster),
    major_cell_type = recode(as.character(major_cell_type), quiescent = "degenerated"),
    final_annotation_reference = as.character(final_annotation)
  ) |>
  distinct(fine_subcluster, .keep_all = TRUE)

fine_to_celltype <- setNames(subcluster_map$major_cell_type, subcluster_map$fine_subcluster)
fine_to_annotation <- setNames(
  subcluster_map$final_annotation_reference,
  subcluster_map$fine_subcluster
)

anno <- anno |>
  transmute(
    cell_id = as.character(cell_id),
    fine_subcluster = as.character(pred_fine_subcluster_hybrid),
    final_annotation = as.character(pred_final_annotation_hybrid)
  ) |>
  mutate(
    final_annotation_from_reference = unname(fine_to_annotation[fine_subcluster]),
    final_annotation = ifelse(
      is.na(final_annotation) | final_annotation == "" | final_annotation == fine_subcluster,
      final_annotation_from_reference,
      final_annotation
    ),
    cell_type = unname(fine_to_celltype[fine_subcluster]),
    cell_type = recode(as.character(cell_type), quiescent = "degenerated"),
    is_ambiguous = fine_subcluster %in% c("Ambiguous", "Ambiguous_subcluster") |
      final_annotation %in% c("Ambiguous", "Ambiguous_subcluster") |
      is.na(fine_subcluster) |
      is.na(final_annotation)
  ) |>
  filter(!is_ambiguous, !is.na(cell_type), cell_type %in% celltype_order) |>
  select(cell_id, fine_subcluster, final_annotation, cell_type)
cat("Loading Xenium cell coordinates:\n  ", cells_csv, "\n", sep = "")
cells_raw <- read_csv(cells_csv, show_col_types = FALSE)

cell_id_col <- pick_column(cells_raw, c("cell_id", "CellID", "cell", "barcode"))
x_col <- pick_column(
  cells_raw,
  c("x_centroid", "x_center", "x_location", "xglobal_px", "X", "x_centroid_px", "x_centroid_um")
)
y_col <- pick_column(
  cells_raw,
  c("y_centroid", "y_center", "y_location", "yglobal_px", "Y", "y_centroid_px", "y_centroid_um")
)

cells <- cells_raw |>
  transmute(
    cell_id = as.character(.data[[cell_id_col]]),
    x = suppressWarnings(as.numeric(.data[[x_col]])),
    y = suppressWarnings(as.numeric(.data[[y_col]]))
  ) |>
  filter(!is.na(cell_id), is.finite(x), is.finite(y))

cell_join_rate <- mean(anno$cell_id %in% cells$cell_id)

if (!is.finite(cell_join_rate) || cell_join_rate < 0.90) {
  stop(
    "Cell coordinate join rate is too low: ",
    round(100 * cell_join_rate, 2),
    "%. Check that the annotation CSV and Xenium cells.csv.gz are from the same dataset."
  )
}

cat("Loading medulla ROI coordinates:\n  ", medulla_roi_csv, "\n", sep = "")
medulla <- read_roi_coordinates(medulla_roi_csv)

dat_all <- anno |>
  inner_join(cells, by = "cell_id") |>
  mutate(
    compartment = ifelse(
      point_in_polygon(x, y, medulla$x, medulla$y),
      "Medulla",
      "Cortex"
    ),
    compartment = factor(compartment, levels = c("Cortex", "Medulla")),
    cell_type = factor(cell_type, levels = celltype_order)
  )

annotation_totals <- dat_all |>
  count(cell_type, fine_subcluster, final_annotation, name = "total_xenium_cells") |>
  arrange(cell_type, desc(total_xenium_cells)) |>
  mutate(
    included_in_annotation_analysis = total_xenium_cells >= min_cells_per_annotation
  )

included_subclusters <- annotation_totals |>
  filter(included_in_annotation_analysis) |>
  pull(fine_subcluster)

dat_annotation <- dat_all |>
  filter(fine_subcluster %in% included_subclusters)

# -------------------------------------------------------------------------
# Tables.
# -------------------------------------------------------------------------
cell_assignments <- dat_all |>
  left_join(
    annotation_totals |>
      select(fine_subcluster, total_xenium_cells, included_in_annotation_analysis),
    by = "fine_subcluster"
  ) |>
  arrange(cell_type, final_annotation, cell_id)

write_csv(
  cell_assignments,
  file.path(out_table_dir, "xenium_medulla_cortex_cell_assignments.csv")
)

write_csv(
  annotation_totals,
  file.path(out_table_dir, "xenium_medulla_cortex_annotation_total_counts.csv")
)

celltype_counts <- dat_all |>
  count(compartment, cell_type, name = "n")

celltype_grid <- tidyr::expand_grid(
  compartment = factor(c("Cortex", "Medulla"), levels = c("Cortex", "Medulla")),
  cell_type = factor(celltype_order, levels = celltype_order)
)

celltype_abund <- celltype_grid |>
  left_join(celltype_counts, by = c("compartment", "cell_type")) |>
  mutate(n = ifelse(is.na(n), 0L, n)) |>
  group_by(compartment) |>
  mutate(
    total_cells = sum(n),
    percent = ifelse(total_cells > 0, 100 * n / total_cells, 0)
  ) |>
  ungroup()

write_csv(
  celltype_abund,
  file.path(out_table_dir, "xenium_medulla_cortex_abundance_by_cell_type.csv")
)
compartment_totals <- dat_all |>
  count(compartment, name = "total_nonambiguous_cells")

annotation_counts <- dat_annotation |>
  count(compartment, cell_type, fine_subcluster, final_annotation, name = "n")

annotation_grid <- tidyr::expand_grid(
  compartment = factor(c("Cortex", "Medulla"), levels = c("Cortex", "Medulla")),
  fine_subcluster = included_subclusters
) |>
  left_join(
    annotation_totals |>
      select(cell_type, fine_subcluster, final_annotation),
    by = "fine_subcluster"
  )

annotation_abund <- annotation_grid |>
  left_join(
    annotation_counts,
    by = c("compartment", "cell_type", "fine_subcluster", "final_annotation")
  ) |>
  left_join(compartment_totals, by = "compartment") |>
  mutate(
    n = ifelse(is.na(n), 0L, n),
    total_nonambiguous_cells = ifelse(is.na(total_nonambiguous_cells), 0L, total_nonambiguous_cells),
    percent = ifelse(total_nonambiguous_cells > 0, 100 * n / total_nonambiguous_cells, 0),
    cell_type = factor(as.character(cell_type), levels = celltype_order)
  ) |>
  arrange(cell_type, final_annotation, compartment)

write_csv(
  annotation_abund,
  file.path(out_table_dir, "xenium_medulla_cortex_abundance_by_annotation.csv")
)

annotation_enrichment <- annotation_abund |>
  select(cell_type, fine_subcluster, final_annotation, compartment, n, percent) |>
  pivot_wider(
    names_from = compartment,
    values_from = c(n, percent),
    values_fill = 0
  ) |>
  mutate(
    log2_medulla_vs_cortex_percent = log2(
      (percent_Medulla + pseudocount_percent) /
        (percent_Cortex + pseudocount_percent)
    ),
    delta_percent_medulla_minus_cortex = percent_Medulla - percent_Cortex,
    abs_delta_percent = abs(delta_percent_medulla_minus_cortex),
    total_xenium_cells = n_Cortex + n_Medulla,
    cell_type = factor(as.character(cell_type), levels = celltype_order)
  ) |>
  arrange(cell_type, log2_medulla_vs_cortex_percent)

write_csv(
  annotation_enrichment,
  file.path(out_table_dir, "xenium_medulla_cortex_annotation_enrichment.csv")
)

roi_summary <- tibble(
  sample_id = sample_id,
  gestational_week = gestational_week,
  xenium_cells_with_nonambiguous_annotation = nrow(dat_all),
  cells_in_cortex = sum(dat_all$compartment == "Cortex"),
  cells_in_medulla = sum(dat_all$compartment == "Medulla"),
  percent_cortex = 100 * cells_in_cortex / xenium_cells_with_nonambiguous_annotation,
  percent_medulla = 100 * cells_in_medulla / xenium_cells_with_nonambiguous_annotation,
  medulla_roi_vertices = nrow(medulla),
  min_cells_per_annotation = min_cells_per_annotation,
  annotations_included = length(included_subclusters),
  coordinate_join_rate = cell_join_rate
)

write_csv(
  roi_summary,
  file.path(out_table_dir, "xenium_medulla_cortex_roi_summary.csv")
)

# -------------------------------------------------------------------------
# Figures.
# -------------------------------------------------------------------------
p_celltype_stacked <- ggplot(
  celltype_abund,
  aes(x = compartment, y = percent, fill = cell_type)
) +
  geom_col(width = 0.75) +
  scale_fill_manual(values = celltype_colors[celltype_order], drop = FALSE) +
  labs(
    x = NULL,
    y = "Relative abundance (%)",
    fill = "Cell type"
  ) +
  theme_publication()

save_plot_pdf_png(
  p_celltype_stacked,
  file.path(out_figure_dir, "xenium_medulla_cortex_cell_type_relative_abundance"),
  width = 6,
  height = 4.5
)
top_n_each <- 10

top_bias <- bind_rows(
  annotation_enrichment |>
    arrange(desc(log2_medulla_vs_cortex_percent)) |>
    slice_head(n = top_n_each) |>
    mutate(direction = "Medulla-biased"),
  annotation_enrichment |>
    arrange(log2_medulla_vs_cortex_percent) |>
    slice_head(n = top_n_each) |>
    mutate(direction = "Cortex-biased")
) |>
  distinct(final_annotation, .keep_all = TRUE) |>
  arrange(log2_medulla_vs_cortex_percent) |>
  mutate(
    cell_type = factor(as.character(cell_type), levels = celltype_order),
    final_annotation = factor(final_annotation, levels = final_annotation)
  )

p_lollipop <- ggplot(
  top_bias,
  aes(
    x = final_annotation,
    y = log2_medulla_vs_cortex_percent,
    fill = cell_type
  )
) +
  geom_hline(yintercept = 0, linewidth = 0.35) +
  geom_col(width = 0.78) +
  coord_flip() +
  scale_fill_manual(values = celltype_colors[celltype_order], drop = FALSE) +
  labs(
    x = NULL,
    y = "log2((medulla % + 0.01) / (cortex % + 0.01))",
    fill = "Cell type"
  ) +
  theme_publication() +
  theme(
    legend.position = "right",
    axis.text.y = element_text(size = 7),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

save_plot_pdf_png(
  p_lollipop,
  file.path(out_figure_dir, "xenium_medulla_cortex_annotation_log2_enrichment"),
  width = 8.2,
  height = 6.3
)

plot_cells <- dat_all

if (nrow(plot_cells) > max_spatial_plot_cells) {
  plot_cells <- plot_cells |>
    slice_sample(n = max_spatial_plot_cells)
}

medulla_closed <- bind_rows(medulla, medulla[1, , drop = FALSE])

p_roi <- ggplot() +
  geom_point(
    data = plot_cells,
    aes(x = x, y = y, color = compartment),
    size = 0.08,
    alpha = 0.6
  ) +
  geom_path(
    data = medulla_closed,
    aes(x = x, y = y),
    linewidth = 0.45,
    color = "black"
  ) +
  coord_equal() +
  labs(
    x = "x",
    y = "y",
    color = "Compartment"
  ) +
  theme_publication()

save_plot_pdf_png(
  p_roi,
  file.path(out_figure_dir, "xenium_medulla_cortex_roi_assignment_sanity_check"),
  width = 7,
  height = 6
)
# -------------------------------------------------------------------------
# Logs.
# -------------------------------------------------------------------------
summary_lines <- c(
  paste("Input Xenium annotation table:", annotation_csv),
  paste("Input Xenium manifest:", manifest_file),
  paste("Input cells.csv.gz:", cells_csv),
  paste("Input medulla ROI:", medulla_roi_csv),
  paste("Sample ID:", sample_id),
  paste("Gestational week:", gestational_week),
  paste("Cell coordinate join rate:", round(100 * cell_join_rate, 3), "%"),
  paste("Nonambiguous annotated cells:", nrow(dat_all)),
  paste("Cortex cells:", sum(dat_all$compartment == "Cortex")),
  paste("Medulla cells:", sum(dat_all$compartment == "Medulla")),
  paste("Minimum cells per annotation:", min_cells_per_annotation),
  paste("Annotations included:", length(included_subclusters)),
  "",
  "Generated figures:",
  paste("  ", file.path(out_figure_dir, "xenium_medulla_cortex_cell_type_relative_abundance.pdf")),
  paste("  ", file.path(out_figure_dir, "xenium_medulla_cortex_annotation_log2_enrichment.pdf")),
  paste("  ", file.path(out_figure_dir, "xenium_medulla_cortex_roi_assignment_sanity_check.pdf")),
  "",
  "Generated tables:",
  paste("  ", file.path(out_table_dir, "xenium_medulla_cortex_cell_assignments.csv")),
  paste("  ", file.path(out_table_dir, "xenium_medulla_cortex_annotation_total_counts.csv")),
  paste("  ", file.path(out_table_dir, "xenium_medulla_cortex_abundance_by_cell_type.csv")),
  paste("  ", file.path(out_table_dir, "xenium_medulla_cortex_abundance_by_annotation.csv")),
  paste("  ", file.path(out_table_dir, "xenium_medulla_cortex_annotation_enrichment.csv")),
  paste("  ", file.path(out_table_dir, "xenium_medulla_cortex_roi_summary.csv"))
)

writeLines(
  summary_lines,
  file.path(out_log_dir, "xenium_medulla_cortex_relative_abundance_summary.txt")
)

sink(file.path(out_log_dir, "xenium_medulla_cortex_relative_abundance_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone.\n")
cat("Xenium medulla/cortex outputs written to:\n  ", results_root, "\n", sep = "")
