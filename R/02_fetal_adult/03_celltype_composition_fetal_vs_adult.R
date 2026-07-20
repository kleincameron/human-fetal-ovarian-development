# Fetal-versus-adult cell-type composition figure.
#
# This script uses the fetal-adult integrated Seurat object and generates a
# manuscript-style stacked composition plot with simple connector lines between
# the fetal and adult bars.

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(tibble)
  library(scales)
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

results_root <- file.path(results_base, "fetal_adult_celltype_composition")

out_figure_dir <- file.path(results_root, "figures")
out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")
dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

png_file <- file.path(out_figure_dir, "fetal_adult_celltype_composition.png")
pdf_file <- file.path(out_figure_dir, "fetal_adult_celltype_composition.pdf")

exclude_categories <- c("unknown")
label_min_prop <- 0.02

dataset_levels <- c("fetal", "adult")

# Legend order, top to bottom
celltype_levels <- c(
  "germ",
  "granulosa",
  "stroma",
  "endothelial",
  "mural",
  "immune",
  "erythroid",
  "degenerated",
  "epithelial",
  "theca",
  "glial cell",
  "plasma cell",
  "smooth muscle cell"
)

# Shared manuscript colors + extra adult-only colors matched to the provided legend
celltype_colors_extended <- c(
  germ = "#1B9E77",
  granulosa = "#D95F02",
  stroma = "#7570B3",
  endothelial = "#E7298A",
  mural = "#66A61E",
  immune = "#E6AB02",
  erythroid = "#A6761D",
  degenerated = "#CFCFCF",
  epithelial = "#BFE6E2",
  theca = "#FB8072",
  "glial cell" = "#F04E37",
  "plasma cell" = "#00E41A",
  "smooth muscle cell" = "#1F2BFF"
)

make_fallback_palette <- function(labels) {
  if (length(labels) == 0) {
    return(character(0))
  }
  cols <- grDevices::hcl.colors(length(labels), palette = "Dynamic")
  names(cols) <- labels
  cols
}

save_publication_plot_pdf_png <- function(plot, pdf_file, png_file, width, height, dpi = 600) {
  ggsave(
    pdf_file,
    plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      png_file,
      plot,
      width = width,
      height = height,
      dpi = dpi,
      device = ragg::agg_png,
      background = "white",
      limitsize = FALSE
    )
  } else {
    ggsave(
      png_file,
      plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
}
build_positions <- function(comp_df, stack_levels, dataset_levels) {
  parts <- lapply(dataset_levels, function(ds) {
    d <- comp_df |>
      filter(dataset == ds) |>
      mutate(stack_order = factor(display_cell_type, levels = stack_levels)) |>
      arrange(stack_order)

    d$ymin <- c(0, head(cumsum(d$prop), -1))
    d$ymax <- cumsum(d$prop)
    d$ymid <- (d$ymin + d$ymax) / 2
    d
  })

  bind_rows(parts)
}

message("Loading fetal-adult integrated object:\n  ", integrated_rds)
stopifnot(file.exists(integrated_rds))

integrated <- readRDS(integrated_rds)
stopifnot(inherits(integrated, "Seurat"))

required_cols <- c("dataset", "adult_cell_type_harmonized")
missing_cols <- setdiff(required_cols, colnames(integrated@meta.data))
if (length(missing_cols) > 0) {
  stop("Integrated object is missing required metadata: ", paste(missing_cols, collapse = ", "))
}

meta_cols <- colnames(integrated@meta.data)

# Prefer a fetal metadata column that retains the degenerating group if available.
fetal_label_column <- if ("major_cell_type" %in% meta_cols) {
  "major_cell_type"
} else if ("cell_type_quiescent" %in% meta_cols) {
  "cell_type_quiescent"
} else if ("cell_type" %in% meta_cols) {
  "cell_type"
} else {
  stop("No suitable fetal cell-type column found. Expected one of: major_cell_type, cell_type_quiescent, cell_type")
}

metadata <- integrated@meta.data |>
  as_tibble(rownames = "cell_barcode") |>
  mutate(
    dataset = as.character(dataset),
    fetal_label_raw = as.character(.data[[fetal_label_column]]),
    adult_label_raw = as.character(adult_cell_type_harmonized),
    display_cell_type = case_when(
      dataset == "fetal" ~ fetal_label_raw,
      dataset == "adult" ~ adult_label_raw,
      TRUE ~ NA_character_
    ),
    display_cell_type = dplyr::recode(
      display_cell_type,
      quiescent = "degenerated"
    )
  ) |>
  filter(
    dataset %in% dataset_levels,
    !is.na(display_cell_type),
    display_cell_type != "",
    !display_cell_type %in% exclude_categories
  )

counts_df <- metadata |>
  count(dataset, display_cell_type, name = "n_cells")

extra_levels <- setdiff(unique(counts_df$display_cell_type), celltype_levels)
if (length(extra_levels) > 0) {
  celltype_levels <- c(celltype_levels, sort(extra_levels))
}

counts_df <- metadata |>
  count(dataset, display_cell_type, name = "n_cells") |>
  complete(
    dataset = dataset_levels,
    display_cell_type = celltype_levels,
    fill = list(n_cells = 0)
  ) |>
  group_by(dataset) |>
  mutate(
    dataset_total = sum(n_cells),
    prop = ifelse(dataset_total > 0, n_cells / dataset_total, 0)
  ) |>
  ungroup()

counts_df$display_cell_type <- factor(
  counts_df$display_cell_type,
  levels = celltype_levels
)
write_csv(
  counts_df,
  file.path(out_table_dir, "fetal_adult_celltype_composition_counts_and_proportions.csv")
)

observed_levels <- counts_df |>
  group_by(display_cell_type) |>
  summarise(total_cells = sum(n_cells), .groups = "drop") |>
  filter(total_cells > 0) |>
  pull(display_cell_type) |>
  as.character()

celltype_levels_present <- celltype_levels[celltype_levels %in% observed_levels]

missing_colors <- setdiff(celltype_levels_present, names(celltype_colors_extended))
plot_colors <- c(
  celltype_colors_extended[intersect(names(celltype_colors_extended), celltype_levels_present)],
  make_fallback_palette(missing_colors)
)
plot_colors <- plot_colors[celltype_levels_present]

write_csv(
  tibble(
    display_cell_type = names(plot_colors),
    color = unname(plot_colors)
  ),
  file.path(out_table_dir, "fetal_adult_celltype_composition_colors.csv")
)

# Bottom-to-top stack is reverse of legend order
stack_levels <- rev(celltype_levels_present)

plot_df <- counts_df |>
  filter(as.character(display_cell_type) %in% celltype_levels_present)

plot_df <- build_positions(
  comp_df = plot_df,
  stack_levels = stack_levels,
  dataset_levels = dataset_levels
)

x_positions <- c(fetal = 1, adult = 2)
bar_width <- 0.70
gap_left <- x_positions[["fetal"]] + bar_width / 2
gap_right <- x_positions[["adult"]] - bar_width / 2

plot_df <- plot_df |>
  mutate(
    x = x_positions[dataset],
    xmin = x - bar_width / 2,
    xmax = x + bar_width / 2,
    label = ifelse(prop >= label_min_prop, percent(prop, accuracy = 0.1), NA_character_)
  )

# Plain connector lines between corresponding cumulative boundaries
boundary_df <- plot_df |>
  select(dataset, display_cell_type, ymin, ymax) |>
  pivot_wider(
    names_from = dataset,
    values_from = c(ymin, ymax)
  ) |>
  filter(
    !is.na(ymin_fetal),
    !is.na(ymax_fetal),
    !is.na(ymin_adult),
    !is.na(ymax_adult)
  )
line_df <- bind_rows(
  boundary_df |>
    transmute(y_left = ymin_fetal, y_right = ymin_adult),
  boundary_df |>
    transmute(y_left = ymax_fetal, y_right = ymax_adult)
) |>
  distinct() |>
  arrange(y_left, y_right)

p <- ggplot() +
  geom_segment(
    data = line_df,
    aes(
      x = gap_left,
      xend = gap_right,
      y = y_left,
      yend = y_right
    ),
    inherit.aes = FALSE,
    color = "grey80",
    linewidth = 0.25
  ) +
  geom_rect(
    data = plot_df,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax,
      fill = display_cell_type
    ),
    color = NA
  ) +
  geom_rect(
    data = plot_df |>
      group_by(dataset) |>
      summarise(
        xmin = min(xmin),
        xmax = max(xmax),
        ymin = 0,
        ymax = 1,
        .groups = "drop"
      ),
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    fill = NA,
    color = "grey35",
    linewidth = 0.25
  ) +
  geom_text(
    data = plot_df |>
      filter(!is.na(label)),
    aes(x = x, y = ymid, label = label),
    family = publication_font_family,
    size = 2.2,
    color = "white",
    fontface = "bold"
  ) +
  scale_fill_manual(
    values = plot_colors,
    breaks = celltype_levels_present,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = c(1, 2),
    labels = c("fetal", "adult"),
    limits = c(0.45, 2.55),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  labs(
    x = NULL,
    y = "Percent of cells",
    fill = NULL
  ) +
  theme_publication() +
  theme(
    legend.position = "right",
    legend.text = element_text(
      family = publication_font_family,
      size = publication_base_size
    ),
    legend.key.height = grid::unit(0.38, "cm"),
    legend.key.width = grid::unit(0.38, "cm"),
    axis.text.x = element_text(
      family = publication_font_family,
      size = publication_base_size
    ),
    axis.text.y = element_text(
      family = publication_font_family,
      size = publication_base_size
    )
  )
save_publication_plot_pdf_png(
  plot = p,
  pdf_file = pdf_file,
  png_file = png_file,
  width = 4.8,
  height = 3.3,
  dpi = 600
)

summary_lines <- c(
  paste("Input integrated object:", integrated_rds),
  paste("Output directory:", results_root),
  paste("Fetal label column used:", fetal_label_column),
  paste("Total cells used:", nrow(metadata)),
  paste("Fetal cells used:", sum(metadata$dataset == "fetal")),
  paste("Adult cells used:", sum(metadata$dataset == "adult")),
  paste("Excluded categories:", paste(exclude_categories, collapse = ", ")),
  paste("Displayed cell types:", paste(celltype_levels_present, collapse = ", ")),
  paste("PDF:", pdf_file),
  paste("PNG:", png_file)
)

writeLines(
  summary_lines,
  file.path(out_log_dir, "fetal_adult_celltype_composition_summary.txt")
)

sink(file.path(out_log_dir, "fetal_adult_celltype_composition_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Wrote fetal-adult cell-type composition figure to:\n  ", out_figure_dir, "\n", sep = "")
