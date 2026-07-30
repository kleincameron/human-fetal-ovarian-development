#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
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

outer_boundary_csv <- file.path(
  project_root,
  "metadata",
  "xenium_cortex_outer_boundary.csv"
)

inner_boundary_csv <- file.path(
  project_root,
  "metadata",
  "xenium_cortex_inner_boundary.csv"
)

results_root <- file.path(results_base, "xenium_cortical_centrality")

if (dir.exists(results_root)) {
  unlink(results_root, recursive = TRUE, force = TRUE)
}

out_figure_dir <- file.path(results_root, "figures")
out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

max_spatial_plot_cells <- 100000
# -------------------------------------------------------------------------
# Exact Xenium Explorer group colors used for manuscript Xenium visualization.
# -------------------------------------------------------------------------
xenium_group_colors <- c(
  "Stroma: cortical" = "#9C27B0",
  "Granulosa: supportive pre-granulosa" = "#E69F00",
  "Granulosa: primordial follicle" = "#FDB462",
  "Immune: tissue macrophages" = "#FFD92F",
  "Granulosa: morphogenetic RELN+" = "#FFA54F",
  "Granulosa: signaling" = "#D55E00",
  "Germ: primordial follicle oocytes" = "#2C7FB8",
  "Ambiguous_subcluster" = "#999999",
  "Germ: stalled meiotic" = "#4DBBD5",
  "Endothelial: angiogenic" = "#FF4FA3",
  "Immune: NK T cells" = "#E6AB02",
  "Germ: pachytene/diplotene" = "#00C1A2",
  "Germ: leptotene/zygotene" = "#00A087",
  "Stroma: proliferative progenitors" = "#8E63CE",
  "Mural: pericytes" = "#E41A1C",
  "Granulosa: epithelial-like" = "#F28E2B",
  "Degenerated: clearance-associated follicle cells" = "#B0B0B0",
  "Stroma: medullary" = "#7B61FF",
  "Granulosa: proliferative progenitors" = "#FF7F00",
  "Mural: contractile VSMC" = "#B22222",
  "Stroma: signaling" = "#6A3D9A",
  "Degenerated: atresia-stressed follicle cells" = "#6B6B6B",
  "Erythroid: late" = "#6F3F1F",
  "Stroma: perineural" = "#C77CFF",
  "Granulosa: matrix-remodeling" = "#C77C2E",
  "Erythroid: early" = "#A66A3F",
  "Germ: Degenerating germ cells" = "#1B9E77"
)

# Map GitHub final annotation labels to Xenium Explorer group labels.
final_annotation_to_xenium_group <- c(
  "Cortical stroma" = "Stroma: cortical",
  "Supportive pre-granulosa" = "Granulosa: supportive pre-granulosa",
  "Primordial follicle granulosa" = "Granulosa: primordial follicle",
  "Tissue macrophages" = "Immune: tissue macrophages",
  "Morphogenetic granulosa RELN+" = "Granulosa: morphogenetic RELN+",
  "Signaling granulosa" = "Granulosa: signaling",
  "Primordial follicle oocytes" = "Germ: primordial follicle oocytes",
  "Ambiguous_subcluster" = "Ambiguous_subcluster",
  "Stalled meiotic germ cells" = "Germ: stalled meiotic",
  "Angiogenic endothelial" = "Endothelial: angiogenic",
  "Angiogenic endothelia" = "Endothelial: angiogenic",
  "NK T cells" = "Immune: NK T cells",
  "Pachytene/diplotene germ cells" = "Germ: pachytene/diplotene",
  "Leptotene/zygotene germ cells" = "Germ: leptotene/zygotene",
  "Proliferative stromal progenitors" = "Stroma: proliferative progenitors",
  "Pericytes" = "Mural: pericytes",
  "Epithelial-like granulosa" = "Granulosa: epithelial-like",
  "Clearance-associated degenerating follicle cells" = "Degenerated: clearance-associated follicle cells",
  "Medullary stroma" = "Stroma: medullary",
  "Proliferative granulosa progenitors" = "Granulosa: proliferative progenitors",
  "Contractile VSMC" = "Mural: contractile VSMC",
  "Signaling stroma" = "Stroma: signaling",
  "Atresia-stressed degenerating follicle cells" = "Degenerated: atresia-stressed follicle cells",
  "Late erythroid" = "Erythroid: late",
  "Perineural stroma" = "Stroma: perineural",
  "Matrix-remodeling granulosa" = "Granulosa: matrix-remodeling",
  "Early erythroid" = "Erythroid: early",
  "Degenerating germ cells" = "Germ: Degenerating germ cells"
)

germ_annotation_order <- c(
  "Stalled meiotic germ cells",
  "Leptotene/zygotene germ cells",
  "Pachytene/diplotene germ cells",
  "Primordial follicle oocytes"
)
lineage_specs <- list(
  germ = list(
    cell_type = "germ",
    min_cells = 0L,
    target_annotations = germ_annotation_order,
    order_mode = "manual"
  ),
  granulosa = list(
    cell_type = "granulosa",
    min_cells = 100L,
    target_annotations = NULL,
    order_mode = "mean_centrality"
  )
)

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

read_coordinate_file <- function(path, label) {
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

      if (nrow(out) >= 2) {
        return(out)
      }
    }
  }

  stop("Could not read ", label, " coordinates from: ", path)
}

point_to_segment_distance <- function(px, py, ax, ay, bx, by) {
  abx <- bx - ax
  aby <- by - ay
  apx <- px - ax
  apy <- py - ay

  ab2 <- abx * abx + aby * aby

  if (ab2 == 0) {
    dx <- px - ax
    dy <- py - ay
    return(sqrt(dx * dx + dy * dy))
  }
  t <- (apx * abx + apy * aby) / ab2
  t <- max(0, min(1, t))

  qx <- ax + t * abx
  qy <- ay + t * aby

  dx <- px - qx
  dy <- py - qy

  sqrt(dx * dx + dy * dy)
}

point_to_polyline_distance <- function(px, py, line_df) {
  n <- nrow(line_df)

  if (n < 2) {
    stop("Boundary polyline must have at least 2 points.")
  }

  dists <- numeric(n - 1)

  for (i in seq_len(n - 1)) {
    dists[i] <- point_to_segment_distance(
      px,
      py,
      line_df$x[i],
      line_df$y[i],
      line_df$x[i + 1],
      line_df$y[i + 1]
    )
  }

  min(dists)
}

compute_cortical_metrics <- function(data, outer, inner) {
  n_cells <- nrow(data)
  d_outer <- numeric(n_cells)
  d_inner <- numeric(n_cells)

  for (i in seq_len(n_cells)) {
    d_outer[i] <- point_to_polyline_distance(data$x[i], data$y[i], outer)
    d_inner[i] <- point_to_polyline_distance(data$x[i], data$y[i], inner)

    if (i %% 1000 == 0 || i == n_cells) {
      cat("  processed ", i, " of ", n_cells, " cells\n", sep = "")
    }
  }

  data |>
    mutate(
      d_outer = d_outer,
      d_inner = d_inner,
      depth01 = d_outer / (d_outer + d_inner),
      centrality = 1 - 2 * abs(depth01 - 0.5),
      centrality = pmin(pmax(centrality, 0), 1)
    )
}

p_to_stars <- function(p) {
  if (is.na(p)) return("ns")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

valid_hex_color <- function(x) {
  !is.na(x) & x != "" & grepl("^#[0-9A-Fa-f]{6}$", x)
}

make_annotation_color_vector <- function(final_annotations) {
  final_annotations <- as.character(final_annotations)
  xenium_groups <- unname(final_annotation_to_xenium_group[final_annotations])
  colors <- unname(xenium_group_colors[xenium_groups])

  missing <- !valid_hex_color(colors)
  if (any(missing)) {
    stop(
      "Missing Xenium colors for final annotations:\n  ",
      paste(unique(final_annotations[missing]), collapse = "\n  ")
    )
  }

  names(colors) <- final_annotations
  colors
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

safe_kruskal <- function(data, metric) {
  groups <- unique(as.character(data$final_annotation))

  if (length(groups) < 2) {
    return(tibble(
      metric = metric,
      test = "Kruskal-Wallis",
      statistic = NA_real_,
      df = NA_real_,
      p_value = NA_real_
    ))
  }

  form <- as.formula(paste(metric, "~ final_annotation"))
  fit <- kruskal.test(form, data = data)

  tibble(
    metric = metric,
    test = "Kruskal-Wallis",
    statistic = unname(fit$statistic),
    df = unname(fit$parameter),
    p_value = fit$p.value
  )
}

make_adjacent_stats <- function(data, ordered_groups) {
  if (length(ordered_groups) < 2) {
    return(tibble(
      group_1 = character(),
      group_2 = character(),
      p_value = numeric(),
      p_adj_BH = numeric(),
      label = character(),
      x1 = numeric(),
      x2 = numeric()
    ))
  }

  pairs <- lapply(seq_len(length(ordered_groups) - 1), function(i) {
    c(ordered_groups[[i]], ordered_groups[[i + 1]])
  })

  out <- bind_rows(lapply(seq_along(pairs), function(i) {
    pair <- pairs[[i]]
    x <- data$centrality[data$final_annotation == pair[[1]]]
    y <- data$centrality[data$final_annotation == pair[[2]]]

    p <- if (length(x) >= 2 && length(y) >= 2) {
      suppressWarnings(wilcox.test(x, y)$p.value)
    } else {
      NA_real_
    }

    tibble(
      group_1 = pair[[1]],
      group_2 = pair[[2]],
      p_value = p,
      x1 = i,
      x2 = i + 1
    )
  }))

  out |>
    mutate(
      p_adj_BH = p.adjust(p_value, method = "BH"),
      label = vapply(p_adj_BH, p_to_stars, character(1))
    )
}

plot_centrality_boxplot <- function(data, adjacent_stats, lineage, colors) {
  y_max <- max(data$centrality, na.rm = TRUE)
  y_min <- min(data$centrality, na.rm = TRUE)
  y_range <- y_max - y_min

  if (!is.finite(y_range) || y_range <= 0) {
    y_range <- 1
  }
  sig_df <- adjacent_stats |>
    mutate(
      y = y_max + seq(0.08, by = 0.08, length.out = n()) * y_range
    )

  p <- ggplot(data, aes(x = final_annotation, y = centrality, fill = final_annotation)) +
    geom_boxplot(
      color = "black",
      outlier.size = 0.18,
      width = 0.7,
      linewidth = 0.35
    ) +
    scale_fill_manual(values = colors, drop = FALSE) +
    labs(
      x = NULL,
      y = "Cortical centrality"
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 35, hjust = 1),
      legend.position = "none"
    )

  if (nrow(sig_df) > 0) {
    p <- p +
      coord_cartesian(
        ylim = c(
          min(data$centrality, na.rm = TRUE),
          max(sig_df$y) + 0.08 * y_range
        )
      )

    for (i in seq_len(nrow(sig_df))) {
      p <- p +
        geom_segment(
          data = sig_df[i, ],
          aes(x = x1, xend = x2, y = y, yend = y),
          inherit.aes = FALSE,
          linewidth = 0.35
        ) +
        geom_segment(
          data = sig_df[i, ],
          aes(x = x1, xend = x1, y = y - 0.025 * y_range, yend = y),
          inherit.aes = FALSE,
          linewidth = 0.35
        ) +
        geom_segment(
          data = sig_df[i, ],
          aes(x = x2, xend = x2, y = y - 0.025 * y_range, yend = y),
          inherit.aes = FALSE,
          linewidth = 0.35
        ) +
        annotate(
          "text",
          x = mean(c(sig_df$x1[i], sig_df$x2[i])),
          y = sig_df$y[i] + 0.02 * y_range,
          label = sig_df$label[i],
          size = 3
        )
    }
  }

  if (lineage == "granulosa") {
    p <- p + theme(axis.text.x = element_text(angle = 40, hjust = 1))
  }

  p
}

plot_spatial_sanity <- function(data, outer, inner, colors) {
  plot_data <- data

  if (nrow(plot_data) > max_spatial_plot_cells) {
    plot_data <- plot_data |>
      slice_sample(n = max_spatial_plot_cells)
  }

  ggplot() +
    geom_path(
      data = outer,
      aes(x = x, y = y),
      linewidth = 0.45,
      color = "black"
    ) +
    geom_path(
      data = inner,
      aes(x = x, y = y),
      linewidth = 0.45,
      color = "gray45"
    ) +
    geom_point(
      data = plot_data,
      aes(x = x, y = y, color = final_annotation),
      size = 0.12,
      alpha = 0.7
    ) +
    scale_color_manual(values = colors, drop = FALSE) +
    coord_equal() +
    labs(
      x = "x",
      y = "y",
      color = "Annotation"
    ) +
    theme_publication()
}
run_lineage <- function(lineage, spec, all_data, outer, inner) {
  cat("\n==============================\n")
  cat("Lineage: ", lineage, "\n", sep = "")
  cat("==============================\n")

  lineage_data <- all_data |>
    filter(cell_type == spec$cell_type)

  if (!is.null(spec$target_annotations)) {
    lineage_data <- lineage_data |>
      filter(final_annotation %in% spec$target_annotations)
  }

  counts_before <- lineage_data |>
    count(final_annotation, name = "total_cells") |>
    arrange(desc(total_cells))

  if (spec$min_cells > 0) {
    keep_annotations <- counts_before |>
      filter(total_cells >= spec$min_cells) |>
      pull(final_annotation)

    lineage_data <- lineage_data |>
      filter(final_annotation %in% keep_annotations)
  }

  if (nrow(lineage_data) == 0) {
    stop("No cells remained for lineage: ", lineage)
  }

  cat("Cells per annotation before centrality calculation:\n")
  print(as.data.frame(counts_before))

  cat("\nComputing cortical metrics for ", nrow(lineage_data), " cells.\n", sep = "")
  lineage_data <- compute_cortical_metrics(lineage_data, outer, inner)

  summary_df <- lineage_data |>
    group_by(final_annotation) |>
    summarise(
      n = n(),
      mean_depth01 = mean(depth01, na.rm = TRUE),
      median_depth01 = median(depth01, na.rm = TRUE),
      sd_depth01 = sd(depth01, na.rm = TRUE),
      mean_centrality = mean(centrality, na.rm = TRUE),
      median_centrality = median(centrality, na.rm = TRUE),
      sd_centrality = sd(centrality, na.rm = TRUE),
      .groups = "drop"
    )

  if (spec$order_mode == "manual") {
    annotation_order <- intersect(spec$target_annotations, summary_df$final_annotation)
  } else if (spec$order_mode == "mean_centrality") {
    annotation_order <- summary_df |>
      arrange(mean_centrality) |>
      pull(final_annotation)
  } else {
    annotation_order <- summary_df |>
      arrange(final_annotation) |>
      pull(final_annotation)
  }

  lineage_data <- lineage_data |>
    mutate(
      lineage = lineage,
      final_annotation = factor(final_annotation, levels = annotation_order)
    )

  summary_df <- summary_df |>
    mutate(
      lineage = lineage,
      cell_type = spec$cell_type,
      final_annotation = factor(final_annotation, levels = annotation_order)
    ) |>
    arrange(final_annotation) |>
    select(lineage, cell_type, final_annotation, everything())

  adjacent_stats <- make_adjacent_stats(lineage_data, annotation_order) |>
    mutate(
      lineage = lineage,
      cell_type = spec$cell_type,
      .before = 1
    )

  global_stats <- safe_kruskal(lineage_data, "centrality") |>
    mutate(
      lineage = lineage,
      cell_type = spec$cell_type,
      .before = 1
    )

  colors <- make_annotation_color_vector(annotation_order)
  p_centrality <- plot_centrality_boxplot(
    lineage_data,
    adjacent_stats,
    lineage,
    colors
  )

  p_spatial <- plot_spatial_sanity(
    lineage_data,
    outer,
    inner,
    colors
  )

  centrality_width <- if (lineage == "granulosa") 9.5 else 6.8

  save_plot_pdf_png(
    p_centrality,
    file.path(out_figure_dir, paste0("xenium_", lineage, "_cortical_centrality_boxplot")),
    width = centrality_width,
    height = 4.8
  )

  save_plot_pdf_png(
    p_spatial,
    file.path(out_figure_dir, paste0("xenium_", lineage, "_cortical_boundary_spatial_sanity_check")),
    width = 7,
    height = 6
  )

  list(
    per_cell = lineage_data |>
      arrange(lineage, final_annotation, cell_id),
    summary = summary_df,
    global_stats = global_stats,
    adjacent_stats = adjacent_stats,
    color_map = tibble(
      lineage = lineage,
      final_annotation = annotation_order,
      xenium_group = unname(final_annotation_to_xenium_group[annotation_order]),
      color = unname(colors)
    ),
    run_summary = tibble(
      lineage = lineage,
      cell_type = spec$cell_type,
      cells_analyzed = nrow(lineage_data),
      annotations_analyzed = length(annotation_order),
      min_cells_filter = spec$min_cells,
      order_mode = spec$order_mode,
      color_source = "hardcoded Xenium Explorer manuscript group colors"
    )
  )
}

# -------------------------------------------------------------------------
# Inputs.
# -------------------------------------------------------------------------
stopifnot(file.exists(manifest_file))
stopifnot(file.exists(annotation_csv))
stopifnot(file.exists(snrna_annotation_metadata))
stopifnot(file.exists(outer_boundary_csv))
stopifnot(file.exists(inner_boundary_csv))

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
  filter(!is_ambiguous, !is.na(cell_type)) |>
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
    "%. Check that the annotation table and Xenium cells.csv.gz are from the same dataset."
  )
}

cat("Loading cortical boundaries.\n")
outer_boundary <- read_coordinate_file(outer_boundary_csv, "outer cortical boundary")
inner_boundary <- read_coordinate_file(inner_boundary_csv, "inner cortical boundary")

all_data <- anno |>
  inner_join(cells, by = "cell_id")

cat("Cells available after annotation-coordinate join: ", nrow(all_data), "\n", sep = "")

lineage_results <- lapply(names(lineage_specs), function(lineage) {
  run_lineage(
    lineage = lineage,
    spec = lineage_specs[[lineage]],
    all_data = all_data,
    outer = outer_boundary,
    inner = inner_boundary
  )
})

names(lineage_results) <- names(lineage_specs)

per_cell <- bind_rows(lapply(lineage_results, `[[`, "per_cell"))
summary_table <- bind_rows(lapply(lineage_results, `[[`, "summary"))
global_stats <- bind_rows(lapply(lineage_results, `[[`, "global_stats"))
adjacent_stats <- bind_rows(lapply(lineage_results, `[[`, "adjacent_stats"))
color_map <- bind_rows(lapply(lineage_results, `[[`, "color_map"))
run_summary <- bind_rows(lapply(lineage_results, `[[`, "run_summary"))

write_csv(
  per_cell,
  file.path(out_table_dir, "xenium_cortical_centrality_cell_metrics.csv")
)

write_csv(
  summary_table,
  file.path(out_table_dir, "xenium_cortical_centrality_summary_by_annotation.csv")
)

write_csv(
  global_stats,
  file.path(out_table_dir, "xenium_cortical_centrality_global_stats.csv")
)

write_csv(
  adjacent_stats,
  file.path(out_table_dir, "xenium_cortical_centrality_adjacent_pairwise_stats.csv")
)

write_csv(
  color_map,
  file.path(out_table_dir, "xenium_cortical_centrality_color_map.csv")
)

write_csv(
  run_summary,
  file.path(out_table_dir, "xenium_cortical_centrality_run_summary.csv")
)

summary_lines <- c(
  paste("Input Xenium annotation table:", annotation_csv),
  paste("Input Xenium manifest:", manifest_file),
  paste("Input cells.csv.gz:", cells_csv),
  paste("Input outer boundary:", outer_boundary_csv),
  paste("Input inner boundary:", inner_boundary_csv),
  paste("Sample ID:", sample_id),
  paste("Gestational week:", gestational_week),
  paste("Cell coordinate join rate:", round(100 * cell_join_rate, 3), "%"),
  paste("Cells available after annotation-coordinate join:", nrow(all_data)),
  "",
  "Generated figures:",
  paste("  ", file.path(out_figure_dir, "xenium_germ_cortical_centrality_boxplot.pdf")),
  paste("  ", file.path(out_figure_dir, "xenium_germ_cortical_boundary_spatial_sanity_check.pdf")),
  paste("  ", file.path(out_figure_dir, "xenium_granulosa_cortical_centrality_boxplot.pdf")),
  paste("  ", file.path(out_figure_dir, "xenium_granulosa_cortical_boundary_spatial_sanity_check.pdf")),
  "",
  "Generated tables:",
  paste("  ", file.path(out_table_dir, "xenium_cortical_centrality_cell_metrics.csv")),
  paste("  ", file.path(out_table_dir, "xenium_cortical_centrality_summary_by_annotation.csv")),
  paste("  ", file.path(out_table_dir, "xenium_cortical_centrality_global_stats.csv")),
  paste("  ", file.path(out_table_dir, "xenium_cortical_centrality_adjacent_pairwise_stats.csv")),
  paste("  ", file.path(out_table_dir, "xenium_cortical_centrality_color_map.csv")),
  paste("  ", file.path(out_table_dir, "xenium_cortical_centrality_run_summary.csv")),
  "",
  "Run summaries:",
  paste(capture.output(print(as.data.frame(run_summary))), collapse = "\n")
)

writeLines(
  summary_lines,
  file.path(out_log_dir, "xenium_cortical_centrality_summary.txt")
)

sink(file.path(out_log_dir, "xenium_cortical_centrality_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone.\n")
cat("Xenium cortical centrality outputs written to:\n  ", results_root, "\n", sep = "")
