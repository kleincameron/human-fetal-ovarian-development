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

pf_roi_csv <- file.path(
  project_root,
  "metadata",
  "xenium_primordial_follicle_roi_coordinates.csv"
)

degen_roi_dir <- file.path(
  project_root,
  "metadata",
  "xenium_degenerating_follicle_rois"
)

results_root <- file.path(results_base, "xenium_roi_celltype_composition")

if (dir.exists(results_root)) {
  unlink(results_root, recursive = TRUE, force = TRUE)
}

out_figure_dir <- file.path(results_root, "figures")
out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

# ROI assignment priority:
# 1. Degenerating follicle ROI
# 2. Inner cortex / primordial follicle ROI
# 3. Medulla
# 4. Outer cortex / cortex outside the primordial follicle ROI
roi_order <- c(
  "Degen_ROI",
  "Inner_cortex_PF_ROI",
  "Medulla",
  "Outer_cortex_not_PF"
)

roi_labels <- c(
  Degen_ROI = "Degenerating follicle ROI",
  Inner_cortex_PF_ROI = "Inner cortex / PF ROI",
  Medulla = "Medulla",
  Outer_cortex_not_PF = "Outer cortex / not PF"
)

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
  mural = "#FF0000",
  immune = "#E6AB02",
  erythroid = "#FF0000"
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

celltype_order <- intersect(
  unique(c(celltype_order, fallback_celltype_order)),
  names(celltype_colors)
)

# Figure-specific color overrides requested for this ROI composition plot.
celltype_colors_plot <- celltype_colors
celltype_colors_plot["degenerated"] <- "#CFCFCF"
celltype_colors_plot[c("mural", "erythroid")] <- "#FF0000"

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

read_polygon_coordinates <- function(path, label) {
  stopifnot(file.exists(path))

  for (skip_n in c(0, 1, 2, 3, 4)) {
    dat <- suppressWarnings(
      tryCatch(
        read_csv(path, skip = skip_n, show_col_types = FALSE, progress = FALSE),
        error = function(e) NULL
      )
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

  stop("Could not read ", label, " polygon coordinates from: ", path)
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

point_in_any_polygon <- function(px, py, polygons) {
  inside <- rep(FALSE, length(px))

  for (i in seq_along(polygons)) {
    poly <- polygons[[i]]
    inside <- inside | point_in_polygon(px, py, poly$x, poly$y)
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

# -------------------------------------------------------------------------
# Inputs.
# -------------------------------------------------------------------------
stopifnot(file.exists(manifest_file))
stopifnot(file.exists(annotation_csv))
stopifnot(file.exists(snrna_annotation_metadata))
stopifnot(file.exists(medulla_roi_csv))
stopifnot(file.exists(pf_roi_csv))
stopifnot(dir.exists(degen_roi_dir))

degen_files <- list.files(
  degen_roi_dir,
  pattern = "^degen[0-9]+_coordinates\\.csv$",
  full.names = TRUE
)

degen_file_numbers <- as.integer(
  sub("^degen([0-9]+)_coordinates\\.csv$", "\\1", basename(degen_files))
)

degen_files <- degen_files[order(degen_file_numbers)]
degen_file_numbers <- sort(degen_file_numbers)

if (length(degen_files) == 0) {
  stop("No degenerating follicle ROI coordinate files found in: ", degen_roi_dir)
}

if (!identical(degen_file_numbers, seq_len(length(degen_file_numbers)))) {
  stop("Degenerating follicle ROI files must be consecutively numbered from degen1_coordinates.csv.")
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
cells_csv <- file.path(manifest$xenium_dir[[1]], "cells.csv.gz")

stopifnot(file.exists(cells_csv))

cat("Loading Xenium annotation table:\n  ", annotation_csv, "\n", sep = "")
anno_raw <- read_csv(annotation_csv, show_col_types = FALSE)

metadata_cols <- colnames(anno_raw)

fine_col <- if ("pred_fine_subcluster_hybrid_filtered" %in% metadata_cols) {
  "pred_fine_subcluster_hybrid_filtered"
} else if ("pred_fine_subcluster_hybrid" %in% metadata_cols) {
  "pred_fine_subcluster_hybrid"
} else if ("pred_fine_subcluster_2stage_filtered" %in% metadata_cols) {
  "pred_fine_subcluster_2stage_filtered"
} else {
  stop("Could not find a fine-subcluster annotation column in xenium_cell_annotations.csv.")
}

annotation_col <- if ("pred_final_annotation_hybrid_filtered" %in% metadata_cols) {
  "pred_final_annotation_hybrid_filtered"
} else if ("pred_final_annotation_hybrid" %in% metadata_cols) {
  "pred_final_annotation_hybrid"
} else if ("pred_final_annotation_2stage_filtered" %in% metadata_cols) {
  "pred_final_annotation_2stage_filtered"
} else {
  stop("Could not find a final-annotation column in xenium_cell_annotations.csv.")
}

cat("Using fine-subcluster column: ", fine_col, "\n", sep = "")
cat("Using final-annotation column: ", annotation_col, "\n", sep = "")

cat("Loading snRNA-seq annotation metadata:\n  ", snrna_annotation_metadata, "\n", sep = "")
snrna_md <- read_csv(snrna_annotation_metadata, show_col_types = FALSE)

celltype_col <- if ("major_cell_type" %in% colnames(snrna_md)) {
  "major_cell_type"
} else if ("cell_type" %in% colnames(snrna_md)) {
  "cell_type"
} else {
  stop("snRNA-seq metadata must contain major_cell_type or cell_type.")
}

required_snrna_cols <- c("fine_subcluster", celltype_col)
missing_snrna_cols <- setdiff(required_snrna_cols, colnames(snrna_md))

if (length(missing_snrna_cols) > 0) {
  stop("snRNA-seq metadata missing columns: ", paste(missing_snrna_cols, collapse = ", "))
}

subcluster_map <- snrna_md |>
  transmute(
    fine_subcluster = as.character(fine_subcluster),
    cell_type = recode(as.character(.data[[celltype_col]]), quiescent = "degenerated")
  ) |>
  distinct(fine_subcluster, .keep_all = TRUE)

fine_to_celltype <- setNames(subcluster_map$cell_type, subcluster_map$fine_subcluster)

anno <- anno_raw |>
  transmute(
    cell_id = as.character(cell_id),
    fine_subcluster = as.character(.data[[fine_col]]),
    final_annotation = as.character(.data[[annotation_col]])
  ) |>
  mutate(
    cell_type = unname(fine_to_celltype[fine_subcluster]),
    cell_type = recode(as.character(cell_type), quiescent = "degenerated"),
    is_ambiguous = fine_subcluster %in% c("Ambiguous", "Ambiguous_subcluster") |
      final_annotation %in% c("Ambiguous", "Ambiguous_subcluster") |
      is.na(fine_subcluster) |
      is.na(final_annotation) |
      is.na(cell_type)
  ) |>
  filter(!is_ambiguous) |>
  mutate(
    cell_type = factor(cell_type, levels = celltype_order)
  ) |>
  filter(!is.na(cell_type))
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

join_rate <- mean(anno$cell_id %in% cells$cell_id)

if (!is.finite(join_rate) || join_rate < 0.90) {
  stop(
    "Annotation-coordinate join rate is too low: ",
    round(100 * join_rate, 2),
    "%"
  )
}

cat("Loading ROI polygons.\n")
medulla_roi <- read_polygon_coordinates(medulla_roi_csv, "medulla ROI")
pf_roi <- read_polygon_coordinates(pf_roi_csv, "primordial follicle ROI")

degen_rois <- lapply(degen_files, function(f) {
  read_polygon_coordinates(f, basename(f))
})
names(degen_rois) <- basename(degen_files)

# -------------------------------------------------------------------------
# ROI assignment.
# -------------------------------------------------------------------------
cat("Assigning ROI membership.\n")

dat <- anno |>
  inner_join(cells, by = "cell_id") |>
  mutate(
    inside_medulla = point_in_polygon(x, y, medulla_roi$x, medulla_roi$y),
    inside_pf_roi = point_in_polygon(x, y, pf_roi$x, pf_roi$y),
    inside_degen_roi = point_in_any_polygon(x, y, degen_rois),
    roi_region = case_when(
      inside_degen_roi ~ "Degen_ROI",
      inside_pf_roi ~ "Inner_cortex_PF_ROI",
      inside_medulla ~ "Medulla",
      TRUE ~ "Outer_cortex_not_PF"
    ),
    roi_region = factor(roi_region, levels = roi_order),
    roi_label = factor(
      unname(roi_labels[as.character(roi_region)]),
      levels = unname(roi_labels[roi_order])
    ),
    cell_type = factor(as.character(cell_type), levels = celltype_order)
  )

roi_assignment_summary <- tibble(
  sample_id = sample_id,
  gestational_week = gestational_week,
  total_nonambiguous_cells_with_coordinates = nrow(dat),
  annotation_coordinate_join_rate = join_rate,
  cells_inside_medulla_polygon = sum(dat$inside_medulla),
  cells_inside_pf_roi = sum(dat$inside_pf_roi),
  cells_inside_degen_roi = sum(dat$inside_degen_roi),
  cells_inside_pf_and_degen_roi = sum(dat$inside_pf_roi & dat$inside_degen_roi),
  cells_inside_medulla_and_degen_roi = sum(dat$inside_medulla & dat$inside_degen_roi),
  cells_inside_medulla_and_pf_roi = sum(dat$inside_medulla & dat$inside_pf_roi),
  assigned_degen_roi = sum(dat$roi_region == "Degen_ROI"),
  assigned_inner_cortex_pf_roi = sum(dat$roi_region == "Inner_cortex_PF_ROI"),
  assigned_medulla = sum(dat$roi_region == "Medulla"),
  assigned_outer_cortex_not_pf = sum(dat$roi_region == "Outer_cortex_not_PF"),
  degen_roi_files = length(degen_files)
)

write_csv(
  roi_assignment_summary,
  file.path(out_table_dir, "xenium_roi_celltype_composition_assignment_summary.csv")
)

# -------------------------------------------------------------------------
# Relative abundance table.
# -------------------------------------------------------------------------
roi_totals <- dat |>
  count(roi_region, roi_label, name = "total_cells_in_roi")

celltype_grid <- expand_grid(
  roi_region = factor(roi_order, levels = roi_order),
  cell_type = factor(celltype_order, levels = celltype_order)
) |>
  mutate(
    roi_label = factor(
      unname(roi_labels[as.character(roi_region)]),
      levels = unname(roi_labels[roi_order])
    )
  )

celltype_composition <- dat |>
  count(roi_region, roi_label, cell_type, name = "n") |>
  right_join(
    celltype_grid,
    by = c("roi_region", "roi_label", "cell_type")
  ) |>
  mutate(n = ifelse(is.na(n), 0L, n)) |>
  left_join(
    roi_totals |>
      select(roi_region, total_cells_in_roi),
    by = "roi_region"
  ) |>
  mutate(
    total_cells_in_roi = ifelse(is.na(total_cells_in_roi), 0L, total_cells_in_roi),
    relative_abundance_percent = ifelse(
      total_cells_in_roi > 0,
      100 * n / total_cells_in_roi,
      0
    ),
    color = unname(celltype_colors_plot[as.character(cell_type)])
  ) |>
  arrange(roi_region, cell_type)

write_csv(
  celltype_composition,
  file.path(out_table_dir, "xenium_roi_celltype_composition.csv")
)
# -------------------------------------------------------------------------
# Pie chart figure.
# -------------------------------------------------------------------------
plot_dat <- celltype_composition |>
  filter(total_cells_in_roi > 0) |>
  mutate(
    roi_label = factor(as.character(roi_label), levels = unname(roi_labels[roi_order])),
    cell_type = factor(as.character(cell_type), levels = celltype_order)
  ) |>
  group_by(roi_label) |>
  arrange(cell_type, .by_group = TRUE) |>
  mutate(
    ymax = cumsum(relative_abundance_percent),
    ymin = lag(ymax, default = 0)
  ) |>
  ungroup()

p_pies <- ggplot(
  plot_dat,
  aes(
    ymax = ymax,
    ymin = ymin,
    xmax = 4,
    xmin = 2,
    fill = cell_type
  )
) +
  geom_rect(color = "white", linewidth = 0.25) +
  coord_polar(theta = "y") +
  xlim(c(0, 4)) +
  facet_wrap(~ roi_label, nrow = 1) +
  scale_fill_manual(values = celltype_colors_plot[celltype_order], drop = FALSE) +
  labs(
    fill = "Cell type"
  ) +
  theme_void(base_family = "Helvetica", base_size = 8) +
  theme(
    strip.text = element_text(size = 8, family = "Helvetica"),
    legend.position = "right",
    legend.title = element_text(size = 8, family = "Helvetica"),
    legend.text = element_text(size = 8, family = "Helvetica"),
    legend.key.size = unit(0.30, "cm"),
    panel.spacing = unit(0.10, "lines"),
    plot.margin = margin(4, 4, 4, 4, unit = "pt")
  )

save_plot_pdf_png(
  p_pies,
  file.path(out_figure_dir, "xenium_roi_celltype_composition_piecharts"),
  width = 8.5,
  height = 2.6
)

# -------------------------------------------------------------------------
# Log.
# -------------------------------------------------------------------------
summary_lines <- c(
  paste("Input Xenium annotation table:", annotation_csv),
  paste("Input cells.csv.gz:", cells_csv),
  paste("Input medulla ROI:", medulla_roi_csv),
  paste("Input primordial follicle ROI:", pf_roi_csv),
  paste("Input degenerating follicle ROI directory:", degen_roi_dir),
  paste("Degenerating follicle ROI files:", length(degen_files)),
  paste("Sample ID:", sample_id),
  paste("Gestational week:", gestational_week),
  paste("Fine-subcluster column:", fine_col),
  paste("Final annotation column:", annotation_col),
  paste("Annotation-coordinate join rate:", round(100 * join_rate, 3), "%"),
  "",
  "ROI assignment priority:",
  "1. Degenerating follicle ROI",
  "2. Inner cortex / primordial follicle ROI",
  "3. Medulla",
  "4. Outer cortex / cortex outside the primordial follicle ROI",
  "",
  "Generated figures:",
  paste("  ", file.path(out_figure_dir, "xenium_roi_celltype_composition_piecharts.pdf")),
  "",
  "Generated tables:",
  paste("  ", file.path(out_table_dir, "xenium_roi_celltype_composition_assignment_summary.csv")),
  paste("  ", file.path(out_table_dir, "xenium_roi_celltype_composition.csv")),
  "",
  "ROI assignment summary:",
  paste(capture.output(print(as.data.frame(roi_assignment_summary))), collapse = "\n")
)

writeLines(
  summary_lines,
  file.path(out_log_dir, "xenium_roi_celltype_composition_summary.txt")
)

sink(file.path(out_log_dir, "xenium_roi_celltype_composition_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone.\n")
cat("Xenium ROI cell-type composition outputs written to:\n  ", results_root, "\n", sep = "")
