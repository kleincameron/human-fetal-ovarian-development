suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(cowplot)
})

set.seed(42)

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

input_root <- file.path(results_base, "snRNAseq_annotated_object")
results_root <- file.path(results_base, "snRNAseq_subcluster_proportions")

source(file.path(project_root, "config", "plotting.R"))
source(file.path(project_root, "config", "labels_colors.R"))

annotated_rds <- file.path(
  input_root,
  "objects",
  "fetal_ovary_snRNAseq_canonical_annotated.rds"
)

out_figure_dir <- file.path(results_root, "figures")
out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

fine_subcluster_colors <- c(
  germ_1 = "#F8766D",
  germ_2 = "#C49A00",
  germ_3 = "#53B400",
  germ_4 = "#00C094",
  germ_5 = "#00B6EB",
  germ_7 = "#A58AFF",
  germ_8 = "#FB61D7",
  germ_0 = "#6B6B6B",
  germ_6 = "#B0B0B0",

  granulosa_0 = "#F8766D",
  granulosa_1 = "#D39200",
  granulosa_2 = "#93AA00",
  granulosa_3 = "#00BA38",
  granulosa_4 = "#00C19F",
  granulosa_5 = "#00B9E3",
  granulosa_6 = "#619CFF",
  granulosa_7 = "#DB72FB",
  granulosa_8 = "#FF61C3",

  stroma_0 = "#F8766D",
  stroma_1 = "#A3A500",
  stroma_2 = "#00BF7D",
  stroma_3 = "#00B0F6",
  stroma_4 = "#E76BF3",

  endothelial_0 = "#E7298A",

  mural_0 = "#66A61E",
  mural_1 = "#B2DF8A",

  immune_0 = "#E6AB02",
  immune_1 = "#FDBF6F",

  erythroid_0 = "#A6761D",
  erythroid_1 = "#D95F0E"
)
plot_definitions <- list(
  list(
    plot_id = "germ_including_degenerated",
    output_prefix = "snRNAseq_subcluster_proportions_germ_including_degenerated",
    include_cell_types = c("germ", "degenerated"),
    title = "Germ and degenerating follicle-cell states"
  ),
  list(
    plot_id = "germ_only",
    output_prefix = "snRNAseq_subcluster_proportions_germ_only",
    include_cell_types = c("germ"),
    title = "Germ-cell states"
  ),
  list(
    plot_id = "degenerated",
    output_prefix = "snRNAseq_subcluster_proportions_degenerated",
    include_cell_types = c("degenerated"),
    title = "Degenerating follicle-cell states"
  ),
  list(
    plot_id = "granulosa",
    output_prefix = "snRNAseq_subcluster_proportions_granulosa",
    include_cell_types = c("granulosa"),
    title = "Granulosa states"
  ),
  list(
    plot_id = "stroma",
    output_prefix = "snRNAseq_subcluster_proportions_stroma",
    include_cell_types = c("stroma"),
    title = "Stromal states"
  ),
  list(
    plot_id = "endothelial",
    output_prefix = "snRNAseq_subcluster_proportions_endothelial",
    include_cell_types = c("endothelial"),
    title = "Endothelial states"
  ),
  list(
    plot_id = "mural",
    output_prefix = "snRNAseq_subcluster_proportions_mural",
    include_cell_types = c("mural"),
    title = "Mural states"
  ),
  list(
    plot_id = "immune",
    output_prefix = "snRNAseq_subcluster_proportions_immune",
    include_cell_types = c("immune"),
    title = "Immune states"
  ),
  list(
    plot_id = "erythroid",
    output_prefix = "snRNAseq_subcluster_proportions_erythroid",
    include_cell_types = c("erythroid"),
    title = "Erythroid states"
  )
)

display_annotation <- function(fine_subcluster, final_annotation) {
  case_when(
    fine_subcluster == "germ_0" ~ "Atresia-stressed",
    fine_subcluster == "germ_6" ~ "Clearance-associated",
    TRUE ~ final_annotation
  )
}

order_subclusters <- function(metadata) {
  metadata |>
    distinct(major_cell_type, fine_subcluster, final_annotation) |>
    mutate(
      major_cell_type = factor(as.character(major_cell_type), levels = celltype_order),
      subcluster_number = suppressWarnings(as.integer(sub("^.*_", "", fine_subcluster)))
    ) |>
    arrange(major_cell_type, subcluster_number, fine_subcluster) |>
    pull(fine_subcluster) |>
    unique()
}

make_fallback_colors <- function(missing_subclusters) {
  if (length(missing_subclusters) == 0) {
    return(character(0))
  }

  fallback <- grDevices::hcl.colors(length(missing_subclusters), palette = "Dynamic")
  names(fallback) <- missing_subclusters
  fallback
}
save_plot_with_fixed_panel_width <- function(
  plot,
  filename,
  panel_width = 4.0,
  legend_width = 2.0,
  height = 3.0,
  dpi = 600
) {
  legend <- cowplot::get_legend(
    plot +
      theme(
        legend.position = "right",
        legend.background = element_rect(fill = "white", color = NA),
        legend.box.background = element_rect(fill = "white", color = NA)
      )
  )

  plot_panel <- plot +
    theme(
      legend.position = "none",
      plot.margin = margin(4, 10, 10, 4),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  combined_plot <- cowplot::plot_grid(
    plot_panel,
    legend,
    nrow = 1,
    rel_widths = c(panel_width, legend_width),
    align = "h"
  )

  combined_plot <- cowplot::ggdraw(combined_plot) +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  save_publication_plot(
    combined_plot,
    filename,
    width = panel_width + legend_width,
    height = height,
    dpi = dpi
  )
}

make_proportion_plot <- function(metadata, plot_definition, weeks) {
  plot_id <- plot_definition$plot_id
  include_cell_types <- plot_definition$include_cell_types

  plot_metadata <- metadata |>
    filter(major_cell_type %in% include_cell_types)

  if (nrow(plot_metadata) == 0) {
    message("Skipping ", plot_id, ": no cells found.")
    return(NULL)
  }

  subcluster_levels <- order_subclusters(plot_metadata)

  annotation_lookup <- plot_metadata |>
    distinct(fine_subcluster, final_annotation, major_cell_type) |>
    mutate(
      display_annotation = display_annotation(fine_subcluster, final_annotation),
      fine_subcluster = factor(fine_subcluster, levels = subcluster_levels)
    ) |>
    arrange(fine_subcluster)

  display_annotation_levels <- annotation_lookup$display_annotation
  names(display_annotation_levels) <- annotation_lookup$fine_subcluster

  count_table <- plot_metadata |>
    count(gestational_week, fine_subcluster, final_annotation, major_cell_type, name = "n_cells") |>
    complete(
      gestational_week = weeks,
      fine_subcluster = subcluster_levels,
      fill = list(n_cells = 0L)
    ) |>
    left_join(annotation_lookup, by = "fine_subcluster", suffix = c("", "_lookup")) |>
    mutate(
      final_annotation = coalesce(final_annotation, final_annotation_lookup),
      major_cell_type = coalesce(as.character(major_cell_type), as.character(major_cell_type_lookup))
    ) |>
    select(-any_of(c("final_annotation_lookup", "major_cell_type_lookup"))) |>
    group_by(gestational_week) |>
    mutate(
      total_cells = sum(n_cells),
      proportion = if_else(total_cells > 0, n_cells / total_cells, 0)
    ) |>
    ungroup() |>
    mutate(
      plot_id = plot_id,
      fine_subcluster = factor(fine_subcluster, levels = subcluster_levels),
      final_annotation = factor(final_annotation),
      display_annotation = factor(display_annotation, levels = display_annotation_levels)
    ) |>
    arrange(gestational_week, fine_subcluster)

  subcluster_colors <- fine_subcluster_colors[subcluster_levels]
  missing_colors <- setdiff(subcluster_levels, names(subcluster_colors[!is.na(subcluster_colors)]))

  if (length(missing_colors) > 0) {
    subcluster_colors <- c(subcluster_colors, make_fallback_colors(missing_colors))
  }

  subcluster_colors <- subcluster_colors[subcluster_levels]
  display_annotation_colors <- subcluster_colors
  names(display_annotation_colors) <- display_annotation_levels

  figure_width <- 5.0
  figure_height <- 3.0
  panel_width <- 4.0
  legend_width <- 2.0
  p <- ggplot(
    count_table,
    aes(
      x = gestational_week,
      y = proportion,
      fill = display_annotation,
      group = display_annotation
    )
  ) +
    geom_area(color = NA, alpha = 1) +
    scale_fill_manual(values = display_annotation_colors, drop = FALSE) +
    scale_y_continuous(
      labels = label_percent(accuracy = 1),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_x_continuous(
      breaks = weeks,
      labels = paste0(weeks, "w"),
      expand = expansion(mult = c(0.02, 0.05))
    ) +
    labs(
      x = "Gestational age",
      y = "Proportion within group",
      fill = NULL
    ) +
    theme_publication() +
    theme(
      legend.position = "right",
      legend.title = element_blank(),
      legend.text = element_text(
        family = publication_font_family,
        size = publication_base_size
      ),
      legend.key.height = unit(0.30, "cm"),
      legend.key.width = unit(0.30, "cm")
    )

  save_plot_with_fixed_panel_width(
    p,
    file.path(out_figure_dir, paste0(plot_definition$output_prefix, ".png")),
    panel_width = panel_width,
    legend_width = legend_width,
    height = figure_height,
    dpi = 600
  )

  write_csv(
    count_table,
    file.path(out_table_dir, paste0(plot_definition$output_prefix, "_counts_and_proportions.csv"))
  )

  count_table
}

cat("Loading annotated snRNA-seq object.\n")
stopifnot(file.exists(annotated_rds))

seu <- readRDS(annotated_rds)
stopifnot(inherits(seu, "Seurat"))

required_cols <- c("gestational_week", "major_cell_type", "fine_subcluster", "final_annotation")
missing_cols <- setdiff(required_cols, colnames(seu@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

metadata <- seu@meta.data |>
  mutate(
    gestational_week = as.numeric(gestational_week),
    major_cell_type = as.character(major_cell_type),
    fine_subcluster = as.character(fine_subcluster),
    final_annotation = as.character(final_annotation)
  ) |>
  filter(
    !is.na(gestational_week),
    !is.na(major_cell_type),
    !is.na(fine_subcluster),
    !is.na(final_annotation),
    major_cell_type != "",
    fine_subcluster != "",
    final_annotation != ""
  )

extra_cell_types <- setdiff(unique(metadata$major_cell_type), celltype_order)
if (length(extra_cell_types) > 0) {
  stop("Unexpected major cell types: ", paste(extra_cell_types, collapse = ", "))
}
weeks <- sort(unique(metadata$gestational_week))

plot_tables <- lapply(
  plot_definitions,
  function(plot_definition) make_proportion_plot(metadata, plot_definition, weeks)
)

plot_tables <- plot_tables[!vapply(plot_tables, is.null, logical(1))]

all_plot_table <- bind_rows(plot_tables)

write_csv(
  all_plot_table,
  file.path(out_table_dir, "snRNAseq_subcluster_proportions_by_week_all_plots.csv")
)

color_table <- tibble(
  fine_subcluster = names(fine_subcluster_colors),
  color = unname(fine_subcluster_colors)
) |>
  left_join(
    metadata |>
      distinct(fine_subcluster, final_annotation, major_cell_type) |>
      mutate(display_annotation = display_annotation(fine_subcluster, final_annotation)),
    by = "fine_subcluster"
  ) |>
  arrange(major_cell_type, fine_subcluster)

write_csv(
  color_table,
  file.path(out_table_dir, "snRNAseq_subcluster_proportion_colors.csv")
)

summary_lines <- c(
  paste("Input annotated object:", annotated_rds),
  paste("Cells:", ncol(seu)),
  paste("Genes:", nrow(seu)),
  paste("Gestational weeks:", paste(paste0(weeks, "w"), collapse = ", ")),
  paste("Plots generated:", length(plot_tables)),
  paste("Plot IDs:", paste(vapply(plot_definitions, `[[`, character(1), "plot_id"), collapse = ", ")),
  "Figure-only abbreviated labels: germ_0 = Atresia-stressed; germ_6 = Clearance-associated",
  "Fixed panel layout: graphical panel width is held constant across plots.",
  paste("Output directory:", results_root)
)

writeLines(summary_lines, file.path(out_log_dir, "snRNAseq_subcluster_proportions_summary.txt"))

sink(file.path(out_log_dir, "snRNAseq_subcluster_proportions_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Done.\n")
cat("Subcluster proportion figures written to: ", out_figure_dir, "\n", sep = "")
