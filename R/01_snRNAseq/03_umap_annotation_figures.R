suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(grid)
})

set.seed(42)

project_root <- "/home/liyan/liyan/Final/github_code_for_publication"
input_root <- "/home/liyan/liyan/Final/github_code_for_publication_results/snRNAseq_annotated_object"
results_root <- "/home/liyan/liyan/Final/github_code_for_publication_results/snRNAseq_umap_annotation_figures"

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

make_distinct_palette <- function(n) {
  if (n <= 0) return(character(0))

  h <- seq(0, 1, length.out = n + 1)[seq_len(n)]
  perm <- (seq_len(n) * 7) %% n
  perm[perm == 0] <- n
  h <- h[perm]

  s <- rep(c(0.95, 0.85, 1.00, 0.90), length.out = n)
  v <- rep(c(0.95, 1.00, 0.85, 0.90), length.out = n)

  grDevices::hsv(h = h, s = s, v = v)
}
make_fine_subcluster_colors <- function(fine_subclusters, major_cell_type_colors) {
  fine_subclusters <- as.character(fine_subclusters)
  fine_subclusters <- fine_subclusters[!is.na(fine_subclusters)]

  zero_labels <- fine_subclusters[grepl("_0$", fine_subclusters)]

  zero_colors <- character(0)
  for (label in zero_labels) {
    prefix <- sub("_.*$", "", label)

    if (prefix %in% names(major_cell_type_colors)) {
      zero_colors[label] <- major_cell_type_colors[[prefix]]
    } else {
      zero_colors[label] <- "black"
    }
  }

  other_labels <- setdiff(fine_subclusters, names(zero_colors))
  other_colors <- make_distinct_palette(length(other_labels))
  names(other_colors) <- other_labels

  colors <- c(zero_colors, other_colors)
  colors[fine_subclusters]
}

cat("Loading annotated snRNA-seq object.\n")
stopifnot(file.exists(annotated_rds))

seu <- readRDS(annotated_rds)
stopifnot(inherits(seu, "Seurat"))

required_cols <- c(
  "sample_id",
  "gestational_week",
  "major_cell_type",
  "fine_subcluster",
  "final_annotation",
  "UMAP_1_publication",
  "UMAP_2_publication"
)

missing_cols <- setdiff(required_cols, colnames(seu@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

extra_cell_types <- setdiff(unique(as.character(seu$major_cell_type)), celltype_order)
if (length(extra_cell_types) > 0) {
  stop("Unexpected major cell types: ", paste(extra_cell_types, collapse = ", "))
}

metadata <- seu@meta.data |>
  mutate(
    major_cell_type = factor(as.character(major_cell_type), levels = celltype_order),
    gestational_week = as.numeric(gestational_week),
    week_label = paste0(gestational_week, "w"),
    fine_subcluster = as.character(fine_subcluster),
    final_annotation = as.character(final_annotation),
    UMAP_1_publication = as.numeric(UMAP_1_publication),
    UMAP_2_publication = as.numeric(UMAP_2_publication)
  )
if (anyNA(metadata$UMAP_1_publication) || anyNA(metadata$UMAP_2_publication)) {
  stop("Publication UMAP coordinate columns contain missing values.")
}

week_levels <- metadata |>
  distinct(gestational_week, week_label) |>
  arrange(gestational_week) |>
  pull(week_label)

fine_order <- metadata |>
  distinct(fine_subcluster, final_annotation, major_cell_type) |>
  mutate(
    major_cell_type = factor(as.character(major_cell_type), levels = celltype_order),
    subcluster_number = suppressWarnings(as.integer(sub("^.*_", "", fine_subcluster)))
  ) |>
  arrange(major_cell_type, subcluster_number, fine_subcluster) |>
  pull(fine_subcluster)

final_annotation_levels <- metadata |>
  distinct(fine_subcluster, final_annotation, major_cell_type) |>
  mutate(fine_subcluster = factor(fine_subcluster, levels = fine_order)) |>
  arrange(fine_subcluster) |>
  pull(final_annotation)

metadata <- metadata |>
  mutate(
    week_label = factor(week_label, levels = week_levels),
    fine_subcluster = factor(fine_subcluster, levels = fine_order),
    final_annotation = factor(final_annotation, levels = final_annotation_levels)
  )

plot_df <- data.frame(
  UMAP_1 = metadata$UMAP_1_publication,
  UMAP_2 = metadata$UMAP_2_publication,
  major_cell_type = metadata$major_cell_type,
  gestational_week = metadata$week_label,
  fine_subcluster = metadata$fine_subcluster,
  final_annotation = metadata$final_annotation,
  stringsAsFactors = FALSE
)

major_cell_type_colors <- celltype_colors[celltype_order]
major_cell_type_colors <- major_cell_type_colors[!is.na(major_cell_type_colors)]

fine_subcluster_colors <- make_fine_subcluster_colors(
  fine_subclusters = fine_order,
  major_cell_type_colors = major_cell_type_colors
)
final_annotation_colors <- fine_subcluster_colors[fine_order]
names(final_annotation_colors) <- final_annotation_levels

write_csv(
  metadata |>
    count(major_cell_type, name = "n_cells") |>
    arrange(major_cell_type),
  file.path(out_table_dir, "snRNAseq_umap_counts_by_major_cell_type.csv")
)

write_csv(
  metadata |>
    count(fine_subcluster, final_annotation, major_cell_type, name = "n_cells") |>
    arrange(major_cell_type, fine_subcluster),
  file.path(out_table_dir, "snRNAseq_umap_counts_by_final_annotation.csv")
)

write_csv(
  metadata |>
    count(gestational_week, week_label, name = "n_cells") |>
    arrange(gestational_week),
  file.path(out_table_dir, "snRNAseq_umap_counts_by_gestational_week.csv")
)

point_size <- 0.12
point_alpha <- 1

legend_theme <- theme(
  legend.title = element_blank(),
  legend.text = element_text(
    family = publication_font_family,
    size = publication_base_size
  ),
  legend.key.height = unit(0.36, "cm"),
  legend.key.width = unit(0.36, "cm")
)

base_umap_theme <- theme_publication() +
  legend_theme +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

p_major_cell_type <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = major_cell_type)) +
  geom_point(shape = 16, size = point_size, alpha = point_alpha, stroke = 0) +
  scale_color_manual(values = major_cell_type_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 2.0, alpha = 1))) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  base_umap_theme

p_gestational_week <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = gestational_week)) +
  geom_point(shape = 16, size = point_size, alpha = point_alpha, stroke = 0) +
  guides(color = guide_legend(override.aes = list(size = 2.0, alpha = 1))) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  base_umap_theme
p_final_annotation <- ggplot(plot_df, aes(UMAP_1, UMAP_2, color = final_annotation)) +
  geom_point(shape = 16, size = point_size, alpha = point_alpha, stroke = 0) +
  scale_color_manual(values = final_annotation_colors, drop = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 1.8, alpha = 1), ncol = 1)) +
  labs(x = "UMAP 1", y = "UMAP 2") +
  base_umap_theme +
  theme(
    legend.key.height = unit(0.28, "cm"),
    legend.key.width = unit(0.32, "cm")
  )

save_publication_plot(
  p_major_cell_type,
  file.path(out_figure_dir, "snRNAseq_umap_major_cell_type.png"),
  width = 4.0,
  height = 3.6,
  dpi = 600
)

save_publication_plot(
  p_gestational_week,
  file.path(out_figure_dir, "snRNAseq_umap_gestational_week.png"),
  width = 4.0,
  height = 3.6,
  dpi = 600
)

save_publication_plot(
  p_final_annotation,
  file.path(out_figure_dir, "snRNAseq_umap_final_annotation.png"),
  width = 6.4,
  height = 5.8,
  dpi = 600
)

summary_lines <- c(
  paste("Input annotated object:", annotated_rds),
  paste("Cells:", ncol(seu)),
  paste("Genes:", nrow(seu)),
  paste("UMAP source:", "metadata columns UMAP_1_publication and UMAP_2_publication"),
  paste("Major cell types:", paste(intersect(celltype_order, unique(as.character(metadata$major_cell_type))), collapse = ", ")),
  paste("Gestational weeks:", paste(week_levels, collapse = ", ")),
  paste("Fine subclusters:", length(fine_order)),
  paste("Final annotations:", length(final_annotation_levels)),
  paste("Output directory:", results_root)
)

writeLines(summary_lines, file.path(out_log_dir, "snRNAseq_umap_annotation_figures_summary.txt"))

sink(file.path(out_log_dir, "snRNAseq_umap_annotation_figures_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Done.\n")
cat("UMAP annotation figures written to: ", out_figure_dir, "\n", sep = "")
