#!/usr/bin/env Rscript

required_packages <- c(
  "dplyr",
  "readr",
  "tibble",
  "stringr",
  "ggplot2"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tibble)
  library(stringr)
  library(ggplot2)
})

set.seed(42)
options(bitmapType = "cairo")

project_root <- normalizePath(
  Sys.getenv("PROJECT_ROOT", unset = getwd()),
  mustWork = TRUE
)

if (!file.exists(file.path(project_root, "config", "plotting.R"))) {
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

deg_results_root <- if (exists("snrna_deg_results_root", inherits = FALSE)) {
  snrna_deg_results_root
} else {
  file.path(results_base, "snRNAseq_DEGs_celltype_subcluster")
}

out_root <- if (exists("snrna_two_subcluster_deg_results_root", inherits = FALSE)) {
  snrna_two_subcluster_deg_results_root
} else {
  file.path(results_base, "snRNAseq_two_subcluster_DEG_lollipop")
}

subcluster_top_table <- file.path(
  deg_results_root,
  "tables",
  "subcluster_top50",
  "DEG_subcluster__top50_adaptive_thresholds_combined.csv"
)

if (!file.exists(subcluster_top_table)) {
  stop(
    "Adaptive subcluster top DEG table not found: ", subcluster_top_table,
    "\nRun 01_snRNAseq/14_generate_deg_tables_celltype_subcluster.R first."
  )
}

out_figure_dir <- file.path(out_root, "figures")
out_table_dir <- file.path(out_root, "tables")
out_log_dir <- file.path(out_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

target_cell_types <- trimws(unlist(strsplit(
  Sys.getenv("TWO_SUBCLUSTER_DEG_CELL_TYPES", unset = "degenerated,mural,immune,erythroid"),
  ","
)))
target_cell_types <- target_cell_types[target_cell_types != ""]

top_genes_per_side <- as.integer(Sys.getenv("TWO_SUBCLUSTER_DEG_TOP_GENES_PER_SIDE", unset = "10"))

figure_width <- as.numeric(Sys.getenv("TWO_SUBCLUSTER_DEG_WIDTH", unset = "5.2"))
figure_height <- as.numeric(Sys.getenv("TWO_SUBCLUSTER_DEG_HEIGHT", unset = "3.6"))
figure_dpi <- as.integer(Sys.getenv("TWO_SUBCLUSTER_DEG_DPI", unset = "600"))

safe_file_component <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

natural_subcluster_order <- function(x) {
  x <- unique(as.character(x))
  prefix <- sub("_[0-9]+$", "", x)
  number <- suppressWarnings(as.integer(sub("^.*_([0-9]+)$", "\\1", x)))

  tibble(x = x, prefix = prefix, number = number) |>
    arrange(prefix, number, x) |>
    pull(x)
}

order_two_subclusters <- function(cell_type, subclusters) {
  subclusters <- unique(as.character(subclusters))

  if (cell_type == "degenerated") {
    preferred <- c("germ_0", "germ_6")
    return(c(intersect(preferred, subclusters), setdiff(natural_subcluster_order(subclusters), preferred)))
  }

  natural_subcluster_order(subclusters)
}

mode_string <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return(NA_character_)
  }

  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

shorten_label <- function(x) {
  dplyr::recode(
    as.character(x),
    "Atresia-stressed degenerating follicle cells" = "Atresia-stressed",
    "Clearance-associated degenerating follicle cells" = "Clearance-associated",
    .default = as.character(x)
  )
}

subcluster_plot_colors <- c(
  "Atresia-stressed degenerating follicle cells" = "#6B6B6B",
  "Clearance-associated degenerating follicle cells" = "#B0B0B0",
  "Pericytes" = "#E41A1C",
  "Contractile VSMC" = "#B22222",
  "Tissue macrophages" = "#FFD92F",
  "NK T cells" = "#E6AB02",
  "Late erythroid" = "#6F3F1F",
  "Early erythroid" = "#A66A3F"
)
get_subcluster_plot_color <- function(final_annotation, cell_type) {
  final_annotation <- as.character(final_annotation)

  if (length(final_annotation) == 1 && final_annotation %in% names(subcluster_plot_colors)) {
    return(unname(subcluster_plot_colors[[final_annotation]]))
  }

  if (exists("celltype_colors", inherits = TRUE) && cell_type %in% names(celltype_colors)) {
    return(unname(celltype_colors[[cell_type]]))
  }

  "#999999"
}

read_subcluster_top_table <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE)

  required <- c(
    "cell_type",
    "fine_subcluster",
    "final_annotation",
    "gene",
    "avg_log2FC",
    "pct.1",
    "pct.2",
    "p_val_adj",
    "score"
  )

  missing_required <- setdiff(required, colnames(df))
  if (length(missing_required) > 0) {
    stop(
      "Adaptive subcluster top DEG table is missing required columns: ",
      paste(missing_required, collapse = ", ")
    )
  }

  df |>
    mutate(
      top_table_row = row_number(),
      gene = as.character(gene),
      cell_type = as.character(cell_type),
      fine_subcluster = as.character(fine_subcluster),
      final_annotation = as.character(final_annotation),
      avg_log2FC = as.numeric(avg_log2FC),
      pct.1 = as.numeric(pct.1),
      pct.2 = as.numeric(pct.2),
      p_val_adj = as.numeric(p_val_adj),
      score = as.numeric(score),
      pct_difference = pct.1 - pct.2
    )
}

read_celltype_deg_tables <- function(cell_type_i) {
  out <- subcluster_top_combined |>
    filter(.data$cell_type == .env$cell_type_i)

  if (nrow(out) == 0) {
    warning("Skipping ", cell_type_i, ": no rows found in adaptive subcluster top DEG table.")
  }

  out
}

select_two_sided_markers <- function(deg_tbl, cell_type_i) {
  subclusters <- order_two_subclusters(cell_type_i, unique(deg_tbl$fine_subcluster))

  if (length(subclusters) != 2) {
    warning(
      "Skipping ", cell_type_i,
      ": expected exactly two subclusters, found ", length(subclusters), "."
    )
    return(tibble())
  }

  selected <- deg_tbl |>
    filter(
      !is.na(gene),
      !is.na(avg_log2FC),
      !is.na(pct.1),
      !is.na(pct.2),
      !is.na(p_val_adj),
      !is.na(score)
    ) |>
    mutate(selection_subcluster = fine_subcluster) |>
    group_by(fine_subcluster) |>
    arrange(top_table_row, .by_group = TRUE) |>
    slice_head(n = top_genes_per_side) |>
    ungroup()

  if (nrow(selected) == 0) {
    return(tibble())
  }

  selected |>
    mutate(
      side = ifelse(selection_subcluster == subclusters[[1]], "left", "right"),
      plot_log2FC = ifelse(side == "left", -abs(avg_log2FC), abs(avg_log2FC)),
      side_order = ifelse(side == "left", 1L, 2L)
    ) |>
    arrange(side_order, desc(abs(plot_log2FC)), top_table_row, gene)
}

save_plot_dual <- function(plot_obj, filename_prefix, width, height, dpi) {
  pdf_file <- paste0(filename_prefix, ".pdf")
  png_file <- paste0(filename_prefix, ".png")

  ggsave(
    filename = pdf_file,
    plot = plot_obj,
    width = width,
    height = height,
    device = "pdf",
    family = publication_font_family,
    useDingbats = FALSE,
    bg = "white",
    limitsize = FALSE
  )

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      filename = png_file,
      plot = plot_obj,
      width = width,
      height = height,
      dpi = dpi,
      device = ragg::agg_png,
      background = "white",
      limitsize = FALSE
    )
  } else {
    ggsave(
      filename = png_file,
      plot = plot_obj,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }

  invisible(c(pdf_file, png_file))
}

make_lollipop_plot <- function(plot_tbl, cell_type) {
  subclusters <- order_two_subclusters(cell_type, unique(plot_tbl$fine_subcluster))

  labels_tbl <- plot_tbl |>
    group_by(fine_subcluster) |>
    summarise(
      final_annotation = mode_string(final_annotation),
      label = shorten_label(final_annotation),
      .groups = "drop"
    )

  label_map <- labels_tbl$label
  names(label_map) <- labels_tbl$fine_subcluster

  if (any(is.na(label_map[subclusters]))) {
    label_map[subclusters] <- subclusters
  }

  color_map <- vapply(subclusters, function(sc) {
    annotation_i <- labels_tbl$final_annotation[match(sc, labels_tbl$fine_subcluster)]
    get_subcluster_plot_color(annotation_i[[1]], cell_type)
  }, character(1))
  names(color_map) <- label_map[subclusters]

  plot_tbl <- plot_tbl |>
    mutate(
      gene = factor(gene, levels = rev(unique(gene))),
      direction_label = factor(label_map[selection_subcluster], levels = label_map[subclusters])
    )

  max_abs <- max(abs(plot_tbl$plot_log2FC), na.rm = TRUE)
  x_limit <- max(0.5, ceiling(max_abs * 10) / 10)

  x_label <- paste0(
    "log2FC (",
    label_map[[subclusters[[1]]]],
    " \u2190  \u2192 ",
    label_map[[subclusters[[2]]]],
    ")"
  )
  ggplot(plot_tbl, aes(x = plot_log2FC, y = gene, color = direction_label)) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey70") +
    geom_segment(
      aes(x = 0, xend = plot_log2FC, yend = gene),
      linewidth = 0.45,
      alpha = 0.85
    ) +
    geom_point(size = 1.8) +
    scale_x_continuous(limits = c(-x_limit, x_limit)) +
    scale_color_manual(values = color_map, name = NULL) +
    labs(
      title = paste0(cell_type, " \u2014 defining DEGs"),
      x = x_label,
      y = NULL
    ) +
    theme_publication() +
    theme(
      plot.title = element_text(size = publication_base_size, hjust = 0.5),
      axis.text.y = element_text(size = max(4, publication_base_size - 2)),
      axis.text.x = element_text(size = max(5, publication_base_size - 1)),
      axis.title.x = element_text(size = publication_base_size),
      legend.position = "none",
      plot.margin = margin(t = 6, r = 10, b = 6, l = 6)
    )
}

subcluster_top_combined <- read_subcluster_top_table(subcluster_top_table)

write_csv(
  subcluster_top_combined |>
    count(cell_type, fine_subcluster, final_annotation, name = "n_top_table_rows") |>
    arrange(cell_type, fine_subcluster),
  file.path(out_table_dir, "QC_adaptive_subcluster_top_table_rows.csv")
)

summary_rows <- list()
combined_selected <- list()

for (cell_type in target_cell_types) {
  cat("\n============================================================\n")
  cat("Generating two-subcluster DEG lollipop: ", cell_type, "\n", sep = "")
  cat("============================================================\n")

  deg_tbl <- read_celltype_deg_tables(cell_type)

  if (nrow(deg_tbl) == 0) {
    next
  }

  selected <- select_two_sided_markers(deg_tbl, cell_type)

  if (nrow(selected) == 0) {
    next
  }

  write_csv(
    selected,
    file.path(out_table_dir, paste0("two_subcluster_selected_DEGs_", safe_file_component(cell_type), ".csv"))
  )

  p <- make_lollipop_plot(selected, cell_type)

  save_plot_dual(
    plot_obj = p,
    filename_prefix = file.path(
      out_figure_dir,
      paste0("Two_subcluster_defining_DEGs_", safe_file_component(cell_type))
    ),
    width = figure_width,
    height = figure_height,
    dpi = figure_dpi
  )

  combined_selected[[cell_type]] <- selected |>
    mutate(cell_type_plot = cell_type) |>
    relocate(cell_type_plot)

  summary_rows[[cell_type]] <- tibble(
    cell_type = cell_type,
    n_subclusters = length(unique(deg_tbl$fine_subcluster)),
    n_adaptive_top_table_rows_available = nrow(deg_tbl),
    n_selected_genes = nrow(selected),
    top_genes_per_side = top_genes_per_side
  )
}

if (length(combined_selected) > 0) {
  write_csv(
    bind_rows(combined_selected),
    file.path(out_table_dir, "COMBINED_two_subcluster_selected_DEGs.csv")
  )
}

if (length(summary_rows) > 0) {
  write_csv(
    bind_rows(summary_rows),
    file.path(out_table_dir, "two_subcluster_DEG_lollipop_summary.csv")
  )
}

writeLines(
  c(
    paste("Input adaptive subcluster top DEG table:", subcluster_top_table),
    paste("Output root:", out_root),
    paste("Target cell types:", paste(target_cell_types, collapse = ", ")),
    paste("Top genes per side:", top_genes_per_side),
    "Selected genes are taken directly from the adaptive subcluster top DEG table generated by script 14.",
    "This preserves the script-14 marker-selection logic, including score ranking, adaptive thresholds, priority high-prevalence/high-score genes, and top-table exclusion rules.",
    "For each two-subcluster cell type, the first ordered subcluster is plotted to the left and the second ordered subcluster is plotted to the right."
  ),
  file.path(out_log_dir, "two_subcluster_DEG_lollipop_notes.txt")
)

sink(file.path(out_log_dir, "sessionInfo_two_subcluster_DEG_lollipop.txt"))
print(sessionInfo())
sink()

cat("\nDone. Two-subcluster DEG lollipop figures written to:\n", out_root, "\n")
