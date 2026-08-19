#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
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

xenium_rds <- file.path(
  results_base,
  "xenium_annotated_object",
  "objects",
  "fetal_ovary_xenium_annotated.rds"
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

results_root <- file.path(
  results_base,
  "xenium_roi_expression_heatmaps_by_annotation"
)

if (dir.exists(results_root)) {
  unlink(results_root, recursive = TRUE, force = TRUE)
}

out_figure_dir <- file.path(results_root, "figures")
out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

target_annotations <- c(
  pf_oocytes = "Primordial follicle oocytes",
  pachytene_diplotene_germ = "Pachytene/diplotene germ cells",
  pf_granulosa = "Primordial follicle granulosa",
  cortical_stroma = "Cortical stroma"
)

roi_levels <- c(
  "outer_cortex",
  "inner_cortex",
  "degenerating"
)

roi_labels <- c(
  outer_cortex = "Outer cortex",
  inner_cortex = "Inner cortex",
  degenerating = "Degenerating follicles"
)

selection_group_levels <- c(
  "increased_outer_vs_inner",
  "increased_inner_vs_outer",
  "increased_degen_vs_cortex"
)

min_cells_per_de_group <- 25L
min_pct_detected <- 0.10
min_pct_delta <- 0.05
logfc_cut <- 0.25
padj_cut <- 0.05
top_genes_per_direction <- 15L

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

get_counts_matrix <- function(obj, assay) {
  mat <- tryCatch(
    LayerData(obj, assay = assay, layer = "counts"),
    error = function(e) NULL
  )

  if (is.null(mat)) {
    mat <- GetAssayData(obj, assay = assay, slot = "counts")
  }

  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }

  mat
}

get_data_matrix <- function(obj, assay) {
  mat <- tryCatch(
    LayerData(obj, assay = assay, layer = "data"),
    error = function(e) NULL
  )

  if (is.null(mat)) {
    mat <- GetAssayData(obj, assay = assay, slot = "data")
  }

  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }

  mat
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
empty_de_table <- function(label, ident_1, ident_2, n1 = 0L, n2 = 0L) {
  tibble(
    gene = character(),
    contrast = character(),
    avg_log2FC = numeric(),
    pct.1 = numeric(),
    pct.2 = numeric(),
    p_val = numeric(),
    p_val_adj = numeric(),
    ident_1 = character(),
    ident_2 = character(),
    cells_ident_1 = integer(),
    cells_ident_2 = integer()
  )
}

run_marker_test <- function(obj, group_col, ident_1, ident_2, label) {
  keep_cells <- colnames(obj)[obj@meta.data[[group_col]] %in% c(ident_1, ident_2)]

  if (length(keep_cells) == 0) {
    return(empty_de_table(label, ident_1, ident_2))
  }

  sub <- subset(obj, cells = keep_cells)
  Idents(sub) <- factor(sub@meta.data[[group_col]], levels = c(ident_1, ident_2))

  n1 <- sum(Idents(sub) == ident_1)
  n2 <- sum(Idents(sub) == ident_2)

  if (n1 < min_cells_per_de_group || n2 < min_cells_per_de_group) {
    warning(
      "Skipping ",
      label,
      " because one group has too few cells: ",
      ident_1,
      "=",
      n1,
      ", ",
      ident_2,
      "=",
      n2
    )

    return(empty_de_table(label, ident_1, ident_2, n1, n2))
  }

  res <- tryCatch(
    FindMarkers(
      sub,
      assay = "Xenium",
      ident.1 = ident_1,
      ident.2 = ident_2,
      test.use = "wilcox",
      min.pct = 0.01,
      logfc.threshold = 0,
      only.pos = FALSE,
      verbose = FALSE
    ),
    error = function(e) {
      warning("FindMarkers failed for ", label, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(res) || nrow(res) == 0) {
    return(empty_de_table(label, ident_1, ident_2, n1, n2))
  }

  res <- res |>
    rownames_to_column("gene") |>
    as_tibble()

  if (!"avg_log2FC" %in% colnames(res) && "avg_logFC" %in% colnames(res)) {
    res <- res |>
      rename(avg_log2FC = avg_logFC)
  }

  res |>
    mutate(
      contrast = label,
      ident_1 = ident_1,
      ident_2 = ident_2,
      cells_ident_1 = n1,
      cells_ident_2 = n2,
      .before = 2
    )
}

summarize_expression_by_roi <- function(counts_mat, data_mat, groups) {
  groups <- factor(groups, levels = roi_levels)

  bind_rows(lapply(levels(groups), function(g) {
    cells_g <- which(groups == g)

    if (length(cells_g) == 0) {
      return(tibble(
        gene = rownames(counts_mat),
        roi_region = g,
        cells = 0L,
        pct_detected = 0,
        mean_counts = 0,
        mean_log_normalized = 0
      ))
    }

    counts_g <- counts_mat[, cells_g, drop = FALSE]
    data_g <- data_mat[, cells_g, drop = FALSE]

    tibble(
      gene = rownames(counts_mat),
      roi_region = g,
      cells = length(cells_g),
      pct_detected = 100 * Matrix::rowMeans(counts_g > 0),
      mean_counts = Matrix::rowMeans(counts_g),
      mean_log_normalized = Matrix::rowMeans(data_g)
    )
  }))
}
select_heatmap_genes <- function(de_outer_vs_inner, de_degen_vs_cortex) {
  outer_up <- de_outer_vs_inner |>
    mutate(
      pct_delta_outer_minus_inner = pct.1 - pct.2,
      selection_group = "increased_outer_vs_inner",
      selection_score = avg_log2FC *
        pmax(pct_delta_outer_minus_inner, 0) *
        (-log10(pmax(p_val_adj, 1e-300)))
    ) |>
    filter(
      p_val_adj <= padj_cut,
      avg_log2FC >= logfc_cut,
      pct.1 >= min_pct_detected,
      pct_delta_outer_minus_inner >= min_pct_delta
    ) |>
    arrange(desc(selection_score)) |>
    mutate(rank_in_selection_group = row_number()) |>
    slice_head(n = top_genes_per_direction)

  inner_up <- de_outer_vs_inner |>
    mutate(
      pct_delta_inner_minus_outer = pct.2 - pct.1,
      selection_group = "increased_inner_vs_outer",
      selection_score = abs(avg_log2FC) *
        pmax(pct_delta_inner_minus_outer, 0) *
        (-log10(pmax(p_val_adj, 1e-300)))
    ) |>
    filter(
      p_val_adj <= padj_cut,
      avg_log2FC <= -logfc_cut,
      pct.2 >= min_pct_detected,
      pct_delta_inner_minus_outer >= min_pct_delta
    ) |>
    arrange(desc(selection_score)) |>
    mutate(rank_in_selection_group = row_number()) |>
    slice_head(n = top_genes_per_direction)

  degen_up <- de_degen_vs_cortex |>
    mutate(
      pct_delta_degen_minus_cortex = pct.1 - pct.2,
      selection_group = "increased_degen_vs_cortex",
      selection_score = avg_log2FC *
        pmax(pct_delta_degen_minus_cortex, 0) *
        (-log10(pmax(p_val_adj, 1e-300)))
    ) |>
    filter(
      p_val_adj <= padj_cut,
      avg_log2FC >= logfc_cut,
      pct.1 >= min_pct_detected,
      pct_delta_degen_minus_cortex >= min_pct_delta
    ) |>
    arrange(desc(selection_score)) |>
    mutate(rank_in_selection_group = row_number()) |>
    slice_head(n = top_genes_per_direction)

  bind_rows(outer_up, inner_up, degen_up) |>
    mutate(
      selection_group = factor(selection_group, levels = selection_group_levels)
    ) |>
    arrange(selection_group, rank_in_selection_group) |>
    distinct(gene, .keep_all = TRUE) |>
    mutate(
      heatmap_row_order = row_number()
    )
}

make_expression_heatmap <- function(expr_by_roi, heatmap_genes, target_label, target_id) {
  if (nrow(heatmap_genes) == 0) {
    warning("No significant heatmap genes for ", target_label, ". Skipping heatmap.")
    return(tibble())
  }

  gene_order <- heatmap_genes |>
    arrange(heatmap_row_order) |>
    mutate(gene_panel = paste(gene, selection_group, sep = "___"))

  plot_dat <- expr_by_roi |>
    inner_join(
      heatmap_genes |>
        select(gene, selection_group, heatmap_row_order),
      by = "gene"
    ) |>
    mutate(
      target_id = target_id,
      target_label = target_label,
      roi_region = factor(roi_region, levels = roi_levels),
      roi_label = factor(
        unname(roi_labels[as.character(roi_region)]),
        levels = unname(roi_labels[roi_levels])
      ),
      selection_group = factor(as.character(selection_group), levels = selection_group_levels),
      gene_panel = paste(gene, selection_group, sep = "___")
    ) |>
    group_by(gene_panel) |>
    mutate(
      row_mean = mean(mean_log_normalized, na.rm = TRUE),
      row_sd = sd(mean_log_normalized, na.rm = TRUE),
      row_z = ifelse(
        is.finite(row_sd) & row_sd > 0,
        (mean_log_normalized - row_mean) / row_sd,
        0
      ),
      row_z_clipped = pmax(pmin(row_z, 2), -2)
    ) |>
    ungroup() |>
    mutate(
      gene_panel = factor(gene_panel, levels = rev(gene_order$gene_panel))
    )
  p <- ggplot(
    plot_dat,
    aes(
      x = roi_label,
      y = gene_panel,
      fill = row_z_clipped
    )
  ) +
    geom_tile(color = "white", linewidth = 0.25) +
    facet_grid(selection_group ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) gsub("___.*$", "", x)) +
    scale_fill_gradient2(
      low = "#2166AC",
      mid = "white",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2, 2)
    ) +
    labs(
      title = target_label,
      x = NULL,
      y = NULL,
      fill = "Row-scaled\nmean expression"
    ) +
    theme_publication() +
    theme(
      plot.title = element_text(size = 8, family = "Helvetica", hjust = 0.5),
      axis.text.x = element_text(size = 8, family = "Helvetica", angle = 25, hjust = 1),
      axis.text.y = element_text(size = 7, family = "Helvetica"),
      strip.text.y = element_blank(),
      strip.background = element_blank(),
      legend.title = element_text(size = 7, family = "Helvetica"),
      legend.text = element_text(size = 7, family = "Helvetica"),
      panel.spacing.y = unit(0.25, "lines")
    )

  n_genes <- nrow(gene_order)
  heatmap_height <- max(3.2, min(11.5, 1.6 + 0.18 * n_genes))

  save_plot_pdf_png(
    p,
    file.path(out_figure_dir, paste0("xenium_", target_id, "_roi_expression_heatmap")),
    width = 4.8,
    height = heatmap_height
  )

  plot_dat
}

run_target_analysis <- function(target_id, target_label, analysis_md, counts) {
  cat("\n==============================\n")
  cat("Target: ", target_label, "\n", sep = "")
  cat("==============================\n")

  target_md <- analysis_md |>
    filter(final_annotation == target_label, !is.na(roi_region)) |>
    mutate(
      roi_region = factor(roi_region, levels = roi_levels),
      cortex_degen_group = case_when(
        roi_region == "degenerating" ~ "degenerating",
        roi_region %in% c("outer_cortex", "inner_cortex") ~ "whole_cortex",
        TRUE ~ NA_character_
      ),
      cortex_degen_group = factor(cortex_degen_group, levels = c("degenerating", "whole_cortex"))
    )

  if (nrow(target_md) == 0) {
    warning("No cells found for target: ", target_label)

    return(list(
      run_summary = tibble(
        target_id = target_id,
        target_label = target_label,
        cells_total = 0L,
        cells_outer_cortex = 0L,
        cells_inner_cortex = 0L,
        cells_degenerating = 0L,
        de_outer_vs_inner_genes = 0L,
        de_degen_vs_cortex_genes = 0L,
        genes_selected_for_heatmap = 0L
      ),
      de = tibble(),
      selected_genes = tibble(),
      heatmap_source = tibble()
    ))
  }

  target_cells <- intersect(target_md$cell_id, colnames(counts))

  target_md <- target_md |>
    filter(cell_id %in% target_cells) |>
    arrange(match(cell_id, target_cells))

  target_counts <- counts[, target_md$cell_id, drop = FALSE]

  target_obj <- CreateSeuratObject(
    counts = target_counts,
    assay = "Xenium",
    meta.data = target_md |>
      as.data.frame() |>
      column_to_rownames("cell_id")
  )

  target_obj <- NormalizeData(
    target_obj,
    assay = "Xenium",
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )

  counts_target <- get_counts_matrix(target_obj, assay = "Xenium")
  data_target <- get_data_matrix(target_obj, assay = "Xenium")

  expr_by_roi <- summarize_expression_by_roi(
    counts_mat = counts_target,
    data_mat = data_target,
    groups = target_obj$roi_region
  )

  roi_counts <- target_md |>
    count(roi_region, name = "cells") |>
    complete(
      roi_region = factor(roi_levels, levels = roi_levels),
      fill = list(cells = 0L)
    ) |>
    mutate(
      roi_label = unname(roi_labels[as.character(roi_region)]),
      percent = ifelse(sum(cells) > 0, 100 * cells / sum(cells), 0)
    )

  target_obj$roi_region <- factor(target_obj$roi_region, levels = roi_levels)
  target_obj$cortex_degen_group <- factor(
    target_obj$cortex_degen_group,
    levels = c("degenerating", "whole_cortex")
  )

  de_outer_vs_inner <- run_marker_test(
    obj = target_obj,
    group_col = "roi_region",
    ident_1 = "outer_cortex",
    ident_2 = "inner_cortex",
    label = "outer_cortex_vs_inner_cortex"
  )
  de_degen_vs_cortex <- run_marker_test(
    obj = target_obj,
    group_col = "cortex_degen_group",
    ident_1 = "degenerating",
    ident_2 = "whole_cortex",
    label = "degenerating_vs_whole_cortex"
  )

  heatmap_genes <- select_heatmap_genes(
    de_outer_vs_inner = de_outer_vs_inner,
    de_degen_vs_cortex = de_degen_vs_cortex
  )

  de_all <- bind_rows(de_outer_vs_inner, de_degen_vs_cortex) |>
    mutate(
      target_id = target_id,
      target_label = target_label,
      .before = 1
    )

  heatmap_genes <- heatmap_genes |>
    mutate(
      target_id = target_id,
      target_label = target_label,
      .before = 1
    )

  heatmap_source <- make_expression_heatmap(
    expr_by_roi = expr_by_roi,
    heatmap_genes = heatmap_genes,
    target_label = target_label,
    target_id = target_id
  )

  run_summary <- tibble(
    target_id = target_id,
    target_label = target_label,
    cells_total = nrow(target_md),
    cells_outer_cortex = roi_counts$cells[roi_counts$roi_region == "outer_cortex"],
    cells_inner_cortex = roi_counts$cells[roi_counts$roi_region == "inner_cortex"],
    cells_degenerating = roi_counts$cells[roi_counts$roi_region == "degenerating"],
    de_outer_vs_inner_genes = nrow(de_outer_vs_inner),
    de_degen_vs_cortex_genes = nrow(de_degen_vs_cortex),
    genes_selected_for_heatmap = nrow(heatmap_genes)
  )

  list(
    run_summary = run_summary,
    de = de_all,
    selected_genes = heatmap_genes,
    heatmap_source = heatmap_source
  )
}

# -------------------------------------------------------------------------
# Inputs.
# -------------------------------------------------------------------------
stopifnot(file.exists(manifest_file))
stopifnot(file.exists(xenium_rds))
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

cat("Loading Xenium object:\n  ", xenium_rds, "\n", sep = "")
xenium <- readRDS(xenium_rds)
stopifnot(inherits(xenium, "Seurat"))

if (length(xenium@images) > 0) {
  xenium@images <- list()
}

xen_assay <- if ("Xenium" %in% Assays(xenium)) "Xenium" else DefaultAssay(xenium)
DefaultAssay(xenium) <- xen_assay

if (xen_assay != "Xenium") {
  stop("Expected Xenium assay, found: ", xen_assay)
}

metadata_cols <- colnames(xenium@meta.data)
fine_col <- if ("pred_fine_subcluster_hybrid_filtered" %in% metadata_cols) {
  "pred_fine_subcluster_hybrid_filtered"
} else if ("pred_fine_subcluster_hybrid" %in% metadata_cols) {
  "pred_fine_subcluster_hybrid"
} else {
  stop("Could not find fine-subcluster metadata column.")
}

annotation_col <- if ("pred_final_annotation_hybrid_filtered" %in% metadata_cols) {
  "pred_final_annotation_hybrid_filtered"
} else if ("pred_final_annotation_hybrid" %in% metadata_cols) {
  "pred_final_annotation_hybrid"
} else {
  stop("Could not find final-annotation metadata column.")
}

cat("Using fine-subcluster column: ", fine_col, "\n", sep = "")
cat("Using final-annotation column: ", annotation_col, "\n", sep = "")

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

cat("Loading medulla ROI:\n  ", medulla_roi_csv, "\n", sep = "")
medulla_roi <- read_polygon_coordinates(medulla_roi_csv, "medulla ROI")

cat("Loading primordial follicle ROI:\n  ", pf_roi_csv, "\n", sep = "")
pf_roi <- read_polygon_coordinates(pf_roi_csv, "primordial follicle ROI")

cat("Loading degenerating follicle ROIs from:\n  ", degen_roi_dir, "\n", sep = "")
degen_rois <- lapply(degen_files, function(f) {
  read_polygon_coordinates(f, basename(f))
})
names(degen_rois) <- basename(degen_files)

metadata <- xenium@meta.data |>
  rownames_to_column("cell_id") |>
  transmute(
    cell_id = cell_id,
    fine_subcluster = as.character(.data[[fine_col]]),
    final_annotation = as.character(.data[[annotation_col]])
  ) |>
  mutate(
    is_ambiguous = fine_subcluster %in% c("Ambiguous", "Ambiguous_subcluster") |
      final_annotation %in% c("Ambiguous", "Ambiguous_subcluster") |
      is.na(fine_subcluster) |
      is.na(final_annotation)
  ) |>
  filter(!is_ambiguous)

join_rate <- mean(metadata$cell_id %in% cells$cell_id)

if (!is.finite(join_rate) || join_rate < 0.90) {
  stop("Coordinate join rate is too low: ", round(100 * join_rate, 2), "%")
}

cat("Assigning ROI membership.\n")

analysis_md <- metadata |>
  inner_join(cells, by = "cell_id") |>
  mutate(
    inside_degen_roi = point_in_any_polygon(x, y, degen_rois),
    inside_pf_roi = point_in_polygon(x, y, pf_roi$x, pf_roi$y),
    inside_medulla = point_in_polygon(x, y, medulla_roi$x, medulla_roi$y),
    roi_region = case_when(
      inside_degen_roi ~ "degenerating",
      inside_pf_roi ~ "inner_cortex",
      inside_medulla ~ NA_character_,
      TRUE ~ "outer_cortex"
    ),
    roi_region = factor(roi_region, levels = roi_levels)
  )

assignment_summary <- tibble(
  sample_id = sample_id,
  gestational_week = gestational_week,
  total_nonambiguous_cells_with_coordinates = nrow(analysis_md),
  coordinate_join_rate = join_rate,
  cells_inside_degen_roi = sum(analysis_md$inside_degen_roi),
  cells_inside_pf_roi = sum(analysis_md$inside_pf_roi),
  cells_inside_medulla = sum(analysis_md$inside_medulla),
  cells_assigned_outer_cortex = sum(analysis_md$roi_region == "outer_cortex", na.rm = TRUE),
  cells_assigned_inner_cortex = sum(analysis_md$roi_region == "inner_cortex", na.rm = TRUE),
  cells_assigned_degenerating = sum(analysis_md$roi_region == "degenerating", na.rm = TRUE),
  cells_excluded_medulla = sum(is.na(analysis_md$roi_region)),
  degen_roi_files = length(degen_files)
)
write_csv(
  assignment_summary,
  file.path(out_table_dir, "xenium_roi_expression_heatmap_assignment_summary.csv")
)

counts <- get_counts_matrix(xenium, assay = "Xenium")

target_results <- lapply(names(target_annotations), function(target_id) {
  run_target_analysis(
    target_id = target_id,
    target_label = target_annotations[[target_id]],
    analysis_md = analysis_md,
    counts = counts
  )
})

names(target_results) <- names(target_annotations)

run_summary <- bind_rows(lapply(target_results, `[[`, "run_summary"))
de_all <- bind_rows(lapply(target_results, `[[`, "de"))
selected_genes <- bind_rows(lapply(target_results, `[[`, "selected_genes"))
heatmap_source <- bind_rows(lapply(target_results, `[[`, "heatmap_source"))

write_csv(
  run_summary,
  file.path(out_table_dir, "xenium_roi_expression_heatmap_run_summary.csv")
)

write_csv(
  de_all,
  file.path(out_table_dir, "xenium_roi_expression_heatmap_differential_expression.csv")
)

write_csv(
  selected_genes,
  file.path(out_table_dir, "xenium_roi_expression_heatmap_selected_genes.csv")
)

write_csv(
  heatmap_source,
  file.path(out_table_dir, "xenium_roi_expression_heatmap_source_data.csv")
)

summary_lines <- c(
  paste("Input Xenium object:", xenium_rds),
  paste("Input Xenium manifest:", manifest_file),
  paste("Input cells.csv.gz:", cells_csv),
  paste("Input medulla ROI:", medulla_roi_csv),
  paste("Input primordial follicle ROI:", pf_roi_csv),
  paste("Input degenerating follicle ROI directory:", degen_roi_dir),
  paste("Degenerating follicle ROI files:", length(degen_files)),
  paste("Sample ID:", sample_id),
  paste("Gestational week:", gestational_week),
  paste("Fine-subcluster column:", fine_col),
  paste("Final-annotation column:", annotation_col),
  paste("Coordinate join rate:", round(100 * join_rate, 3), "%"),
  "",
  "ROI column order in heatmaps:",
  paste(unname(roi_labels[roi_levels]), collapse = " | "),
  "",
  "Significance contrasts:",
  "1. outer_cortex vs inner_cortex",
  "2. degenerating vs whole_cortex, where whole_cortex = outer_cortex + inner_cortex",
  "",
  "Gene ordering in heatmaps:",
  "1. genes increased in outer cortex vs inner cortex",
  "2. genes increased in inner cortex vs outer cortex",
  "3. genes increased in degenerating ROI vs whole cortex",
  "",
  paste("Minimum cells per DE group:", min_cells_per_de_group),
  paste("Candidate filter: adjusted P <=", padj_cut),
  paste("Candidate filter: absolute log2FC >=", logfc_cut),
  paste("Candidate filter: detection in enriched group >=", min_pct_detected),
  paste("Candidate filter: detection delta >=", min_pct_delta),
  paste("Top genes per direction:", top_genes_per_direction),
  "",
  "Generated figures:",
  paste(
    "  ",
    file.path(out_figure_dir, paste0("xenium_", names(target_annotations), "_roi_expression_heatmap.pdf")),
    collapse = "\n"
  ),
  "",
  "Generated tables:",
  paste("  ", file.path(out_table_dir, "xenium_roi_expression_heatmap_assignment_summary.csv")),
  paste("  ", file.path(out_table_dir, "xenium_roi_expression_heatmap_run_summary.csv")),
  paste("  ", file.path(out_table_dir, "xenium_roi_expression_heatmap_differential_expression.csv")),
  paste("  ", file.path(out_table_dir, "xenium_roi_expression_heatmap_selected_genes.csv")),
  paste("  ", file.path(out_table_dir, "xenium_roi_expression_heatmap_source_data.csv")),
  "",
  "Assignment summary:",
  paste(capture.output(print(as.data.frame(assignment_summary))), collapse = "\n"),
  "",
  "Run summary:",
  paste(capture.output(print(as.data.frame(run_summary))), collapse = "\n")
)
writeLines(
  summary_lines,
  file.path(out_log_dir, "xenium_roi_expression_heatmaps_by_annotation_summary.txt")
)

sink(file.path(out_log_dir, "xenium_roi_expression_heatmaps_by_annotation_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone.\n")
cat("Xenium ROI expression heatmap outputs written to:\n  ", results_root, "\n", sep = "")
