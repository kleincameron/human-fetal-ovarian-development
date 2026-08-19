suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(scales)
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
results_root <- file.path(results_base, "snRNAseq_cell_type_proportions")

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

cat("Loading annotated snRNA-seq object.\n")
stopifnot(file.exists(annotated_rds))

seu <- readRDS(annotated_rds)
stopifnot(inherits(seu, "Seurat"))

required_cols <- c("gestational_week", "major_cell_type")
missing_cols <- setdiff(required_cols, colnames(seu@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}
metadata <- seu@meta.data |>
  mutate(
    gestational_week = as.numeric(gestational_week),
    major_cell_type = as.character(major_cell_type)
  ) |>
  filter(!is.na(gestational_week), !is.na(major_cell_type), major_cell_type != "")

extra_cell_types <- setdiff(unique(metadata$major_cell_type), celltype_order)
if (length(extra_cell_types) > 0) {
  stop("Unexpected major cell types: ", paste(extra_cell_types, collapse = ", "))
}

weeks <- sort(unique(metadata$gestational_week))
cell_types_present <- intersect(celltype_order, unique(metadata$major_cell_type))

proportion_table <- metadata |>
  count(gestational_week, major_cell_type, name = "n_cells") |>
  complete(
    gestational_week = weeks,
    major_cell_type = celltype_order,
    fill = list(n_cells = 0L)
  ) |>
  group_by(gestational_week) |>
  mutate(
    total_cells = sum(n_cells),
    proportion = if_else(total_cells > 0, n_cells / total_cells, 0)
  ) |>
  ungroup() |>
  mutate(
    major_cell_type = factor(major_cell_type, levels = celltype_order),
    week_label = paste0(gestational_week, "w")
  ) |>
  arrange(gestational_week, major_cell_type)

write_csv(
  proportion_table,
  file.path(out_table_dir, "snRNAseq_cell_type_counts_and_proportions_by_week.csv")
)

plot_df <- proportion_table |>
  filter(major_cell_type %in% cell_types_present)

plot_colors <- celltype_colors[celltype_order]
plot_colors <- plot_colors[!is.na(plot_colors)]

p_area <- ggplot(
  plot_df,
  aes(
    x = gestational_week,
    y = proportion,
    fill = major_cell_type,
    group = major_cell_type
  )
) +
  geom_area(color = NA, alpha = 1) +
  scale_fill_manual(values = plot_colors, drop = FALSE) +
  scale_y_continuous(
    labels = label_percent(accuracy = 1),
    limits = c(0, 1),
    expand = c(0, 0)
  ) +
  scale_x_continuous(
    breaks = weeks,
    labels = paste0(weeks, "w"),
    expand = expansion(mult = c(0.01, 0.01))
  ) +
  labs(
    x = "Gestational age",
    y = "Cell proportion",
    fill = NULL
  ) +
  theme_publication() +
  theme(
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(
      family = publication_font_family,
      size = publication_base_size
    )
  )

save_publication_plot(
  p_area,
  file.path(out_figure_dir, "snRNAseq_cell_type_proportions_by_week_area.png"),
  width = 4.8,
  height = 2.8,
  dpi = 600
)

summary_lines <- c(
  paste("Input annotated object:", annotated_rds),
  paste("Cells:", ncol(seu)),
  paste("Genes:", nrow(seu)),
  paste("Gestational weeks:", paste(paste0(weeks, "w"), collapse = ", ")),
  paste("Major cell types:", paste(cell_types_present, collapse = ", ")),
  paste("Output directory:", results_root)
)

writeLines(summary_lines, file.path(out_log_dir, "snRNAseq_cell_type_proportions_summary.txt"))

sink(file.path(out_log_dir, "snRNAseq_cell_type_proportions_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Done.\n")
cat("Cell-type proportion figure written to: ", out_figure_dir, "\n", sep = "")
