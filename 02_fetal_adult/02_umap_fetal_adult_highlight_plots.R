# Fetal-adult integrated UMAP highlight plots.
#
# This script uses the integrated fetal-adult object and generates two matched
# UMAPs:
#   1. adult harmonized cell types colored, fetal cells greyed
#   2. fetal major cell types colored, adult cells greyed

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(cowplot)
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


source(file.path(project_root, "config", "plotting.R"))
source(file.path(project_root, "config", "labels_colors.R"))


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

results_root <- file.path(results_base, "fetal_adult_umap_figures")

out_figure_dir <- file.path(results_root, "figures")
out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

min_cells_per_group <- 10

point_size <- 0.16
point_alpha <- 1
legend_point_size <- 2.4
grey_background <- "grey85"

umap_panel_width <- 4.8
legend_width <- 2.0
figure_height <- 3.6
dpi_out <- 600
adult_extra_colors <- c(
  epithelial = "#8DD3C7",
  theca = "#FB8072",
  glial = "#BEBADA",
  "smooth muscle cell" = "#80B1D3",
  unknown = "#B3B3B3"
)

make_fallback_palette <- function(labels) {
  if (length(labels) == 0) {
    return(character(0))
  }

  colors <- grDevices::hcl.colors(length(labels), palette = "Dynamic")
  names(colors) <- labels

  colors
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

save_umap_with_fixed_panel <- function(
  plot,
  filename_prefix,
  legend_side = "right",
  umap_width = umap_panel_width,
  legend_width_in = legend_width,
  height = figure_height,
  dpi = dpi_out
) {
  legend_side <- match.arg(legend_side, c("right", "left"))

  legend <- cowplot::get_legend(
    plot +
      theme(
        legend.position = legend_side,
        legend.box.margin = margin(0, 0, 0, 0)
      )
  )

  umap <- plot + theme(legend.position = "none")

  combined <- if (legend_side == "right") {
    cowplot::plot_grid(
      umap,
      legend,
      nrow = 1,
      rel_widths = c(umap_width, legend_width_in)
    )
  } else {
    cowplot::plot_grid(
      legend,
      umap,
      nrow = 1,
      rel_widths = c(legend_width_in, umap_width)
    )
  }

  png_file <- file.path(out_figure_dir, paste0(filename_prefix, ".png"))
  pdf_file <- file.path(out_figure_dir, paste0(filename_prefix, ".pdf"))

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      png_file,
      combined,
      width = umap_width + legend_width_in,
      height = height,
      dpi = dpi,
      device = ragg::agg_png,
      background = "white",
      limitsize = FALSE
    )
  } else {
    ggsave(
      png_file,
      combined,
      width = umap_width + legend_width_in,
      height = height,
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
  ggsave(
    pdf_file,
    combined,
    width = umap_width + legend_width_in,
    height = height,
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )

  invisible(c(png = png_file, pdf = pdf_file))
}

make_umap_plot <- function(
  background_df,
  foreground_df,
  color_column,
  color_values
) {
  ggplot() +
    geom_point(
      data = background_df,
      aes(UMAP_1, UMAP_2),
      shape = 16,
      size = point_size,
      alpha = point_alpha,
      stroke = 0,
      color = grey_background
    ) +
    geom_point(
      data = foreground_df,
      aes(
        UMAP_1,
        UMAP_2,
        color = .data[[color_column]]
      ),
      shape = 16,
      size = point_size,
      alpha = point_alpha,
      stroke = 0
    ) +
    scale_color_manual(values = color_values, drop = TRUE) +
    guides(
      color = guide_legend(
        override.aes = list(
          size = legend_point_size,
          alpha = 1
        )
      )
    ) +
    labs(
      x = "UMAP 1",
      y = "UMAP 2",
      color = NULL
    ) +
    theme_publication() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(
        family = publication_font_family,
        size = publication_base_size
      ),
      legend.key.height = grid::unit(0.30, "cm"),
      legend.key.width = grid::unit(0.30, "cm")
    )
}
message("Loading fetal-adult integrated object:\n  ", integrated_rds)
stopifnot(file.exists(integrated_rds))

integrated <- readRDS(integrated_rds)
stopifnot(inherits(integrated, "Seurat"))

required_cols <- c(
  "dataset",
  "cell_type",
  "adult_cell_type_harmonized"
)

missing_cols <- setdiff(required_cols, colnames(integrated@meta.data))

if (length(missing_cols) > 0) {
  stop("Integrated object is missing required metadata: ", paste(missing_cols, collapse = ", "))
}

if (!"umap" %in% names(integrated@reductions)) {
  stop("Integrated object does not contain a UMAP reduction.")
}

embedding <- Embeddings(integrated, reduction = "umap")

if (ncol(embedding) < 2) {
  stop("UMAP reduction has fewer than two dimensions.")
}

metadata <- integrated@meta.data

plot_df <- tibble(
  cell_barcode = rownames(embedding),
  UMAP_1 = embedding[, 1],
  UMAP_2 = embedding[, 2],
  dataset = as.character(metadata[rownames(embedding), "dataset"]),
  fetal_cell_type = NA_character_,
  adult_cell_type_harmonized = NA_character_
)

is_fetal <- plot_df$dataset == "fetal"
is_adult <- plot_df$dataset == "adult"

plot_df$fetal_cell_type[is_fetal] <- as.character(
  metadata[plot_df$cell_barcode[is_fetal], "cell_type"]
)

plot_df$adult_cell_type_harmonized[is_adult] <- as.character(
  metadata[plot_df$cell_barcode[is_adult], "adult_cell_type_harmonized"]
)

if (any(is.na(plot_df$dataset))) {
  stop("Some UMAP cells are missing dataset labels.")
}

adult_counts <- plot_df |>
  filter(dataset == "adult", !is.na(adult_cell_type_harmonized), adult_cell_type_harmonized != "") |>
  count(adult_cell_type_harmonized, name = "n_cells") |>
  arrange(desc(n_cells))

fetal_counts <- plot_df |>
  filter(dataset == "fetal", !is.na(fetal_cell_type), fetal_cell_type != "") |>
  count(fetal_cell_type, name = "n_cells") |>
  arrange(desc(n_cells))

adult_keep <- adult_counts |>
  filter(n_cells >= min_cells_per_group) |>
  pull(adult_cell_type_harmonized)

fetal_keep <- fetal_counts |>
  filter(n_cells >= min_cells_per_group) |>
  pull(fetal_cell_type)

plot_df <- plot_df |>
  mutate(
    adult_cell_type_harmonized = if_else(
      dataset == "adult" & adult_cell_type_harmonized %in% adult_keep,
      adult_cell_type_harmonized,
      NA_character_
    ),
    fetal_cell_type = if_else(
      dataset == "fetal" & fetal_cell_type %in% fetal_keep,
      fetal_cell_type,
      NA_character_
    )
  )
adult_levels <- sort(unique(na.omit(plot_df$adult_cell_type_harmonized)))
fetal_levels <- celltype_order[celltype_order %in% unique(na.omit(plot_df$fetal_cell_type))]

remaining_fetal_levels <- setdiff(
  sort(unique(na.omit(plot_df$fetal_cell_type))),
  fetal_levels
)

fetal_levels <- c(fetal_levels, remaining_fetal_levels)

adult_base_colors <- c(
  celltype_colors,
  adult_extra_colors
)

missing_adult_colors <- setdiff(adult_levels, names(adult_base_colors))
adult_colors <- c(
  adult_base_colors[intersect(names(adult_base_colors), adult_levels)],
  make_fallback_palette(missing_adult_colors)
)
adult_colors <- adult_colors[adult_levels]

missing_fetal_colors <- setdiff(fetal_levels, names(celltype_colors))
fetal_colors <- c(
  celltype_colors[intersect(names(celltype_colors), fetal_levels)],
  make_fallback_palette(missing_fetal_colors)
)
fetal_colors <- fetal_colors[fetal_levels]

adult_foreground <- plot_df |>
  filter(dataset == "adult", !is.na(adult_cell_type_harmonized)) |>
  mutate(adult_cell_type_harmonized = factor(adult_cell_type_harmonized, levels = adult_levels))

fetal_foreground <- plot_df |>
  filter(dataset == "fetal", !is.na(fetal_cell_type)) |>
  mutate(fetal_cell_type = factor(fetal_cell_type, levels = fetal_levels))

adult_background <- plot_df |>
  filter(dataset == "adult")

fetal_background <- plot_df |>
  filter(dataset == "fetal")

p_adult <- make_umap_plot(
  background_df = fetal_background,
  foreground_df = adult_foreground,
  color_column = "adult_cell_type_harmonized",
  color_values = adult_colors
)

p_fetal <- make_umap_plot(
  background_df = adult_background,
  foreground_df = fetal_foreground,
  color_column = "fetal_cell_type",
  color_values = fetal_colors
)

adult_files <- save_umap_with_fixed_panel(
  p_adult,
  filename_prefix = "fetal_adult_umap_adult_colored_fetal_grey",
  legend_side = "right"
)

fetal_files <- save_umap_with_fixed_panel(
  p_fetal,
  filename_prefix = "fetal_adult_umap_fetal_colored_adult_grey",
  legend_side = "right"
)
write_csv(
  adult_counts |>
    mutate(kept_for_plot = adult_cell_type_harmonized %in% adult_keep),
  file.path(out_table_dir, "fetal_adult_umap_adult_harmonized_counts.csv")
)

write_csv(
  fetal_counts |>
    mutate(kept_for_plot = fetal_cell_type %in% fetal_keep),
  file.path(out_table_dir, "fetal_adult_umap_fetal_cell_type_counts.csv")
)

write_csv(
  bind_rows(
    tibble(
      plot = "adult_colored_fetal_grey",
      label = names(adult_colors),
      color = unname(adult_colors)
    ),
    tibble(
      plot = "fetal_colored_adult_grey",
      label = names(fetal_colors),
      color = unname(fetal_colors)
    )
  ),
  file.path(out_table_dir, "fetal_adult_umap_colors.csv")
)

summary_lines <- c(
  paste("Input integrated object:", integrated_rds),
  paste("Output directory:", results_root),
  paste("Cells:", nrow(plot_df)),
  paste("Adult cells:", sum(plot_df$dataset == "adult")),
  paste("Fetal cells:", sum(plot_df$dataset == "fetal")),
  paste("Minimum cells per plotted group:", min_cells_per_group),
  paste("Adult groups kept:", paste(adult_levels, collapse = ", ")),
  paste("Fetal groups kept:", paste(fetal_levels, collapse = ", ")),
  paste("Adult colored PNG:", adult_files[["png"]]),
  paste("Adult colored PDF:", adult_files[["pdf"]]),
  paste("Fetal colored PNG:", fetal_files[["png"]]),
  paste("Fetal colored PDF:", fetal_files[["pdf"]])
)

writeLines(
  summary_lines,
  file.path(out_log_dir, "fetal_adult_umap_figures_summary.txt")
)

sink(file.path(out_log_dir, "fetal_adult_umap_figures_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Wrote fetal-adult UMAP figures to:\n  ", out_figure_dir, "\n", sep = "")
