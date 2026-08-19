suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
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

input_root <- file.path(results_base, "snRNAseq_initial_processing")
results_root <- file.path(results_base, "snRNAseq_qc_figure")

source(file.path(project_root, "config", "plotting.R"))
source(file.path(project_root, "config", "labels_colors.R"))

canonical_rds <- file.path(
  input_root,
  "objects",
  "fetal_ovary_snRNAseq_canonical.rds"
)

out_figure_dir <- file.path(results_root, "figures")
out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading canonical snRNA-seq object.\n")
stopifnot(file.exists(canonical_rds))

seu <- readRDS(canonical_rds)
stopifnot(inherits(seu, "Seurat"))

DefaultAssay(seu) <- "RNA"
required_cols <- c("nFeature_RNA", "nCount_RNA", "sample_id", "gestational_week")
missing_cols <- setdiff(required_cols, colnames(seu@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

if (!"percent.mt" %in% colnames(seu@meta.data)) {
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
}

metadata <- seu@meta.data |>
  mutate(
    gestational_week = as.numeric(gestational_week),
    week_label = paste0(gestational_week, "w"),
    nCount_RNA = as.numeric(nCount_RNA),
    nFeature_RNA = as.numeric(nFeature_RNA),
    percent.mt = as.numeric(percent.mt)
  ) |>
  filter(!is.na(gestational_week), !is.na(week_label))

week_levels <- metadata |>
  distinct(gestational_week, week_label) |>
  arrange(gestational_week) |>
  pull(week_label)

metadata <- metadata |>
  mutate(week_label = factor(week_label, levels = week_levels))

set.seed(42)
violin_metadata <- metadata |>
  group_by(week_label) |>
  group_modify(~ {
    n_take <- min(5000, nrow(.x))
    slice_sample(.x, n = n_take)
  }) |>
  ungroup()

set.seed(42)
scatter_n <- min(120000, nrow(metadata))
scatter_metadata <- metadata |>
  slice_sample(n = scatter_n)

qc_summary <- metadata |>
  group_by(gestational_week, week_label) |>
  summarise(
    n_cells = n(),
    median_nCount_RNA = median(nCount_RNA),
    median_nFeature_RNA = median(nFeature_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  ) |>
  arrange(gestational_week)

write_csv(qc_summary, file.path(out_table_dir, "snRNAseq_QC_summary_by_week.csv"))
axis_x_week <- theme(
  axis.text.x = element_text(
    angle = 45,
    hjust = 1,
    family = publication_font_family,
    size = publication_base_size
  )
)

p_umi <- ggplot(violin_metadata, aes(x = week_label, y = nCount_RNA)) +
  geom_violin(scale = "width", trim = TRUE, fill = "grey85", color = "black", linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.size = 0.15, linewidth = 0.2, fill = "white") +
  scale_y_continuous(labels = label_comma()) +
  labs(x = NULL, y = "UMIs per nucleus") +
  theme_publication() +
  axis_x_week

p_genes <- ggplot(violin_metadata, aes(x = week_label, y = nFeature_RNA)) +
  geom_violin(scale = "width", trim = TRUE, fill = "grey85", color = "black", linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.size = 0.15, linewidth = 0.2, fill = "white") +
  scale_y_continuous(labels = label_comma()) +
  labs(x = NULL, y = "Genes per nucleus") +
  theme_publication() +
  axis_x_week

p_mt <- ggplot(violin_metadata, aes(x = week_label, y = percent.mt)) +
  geom_violin(scale = "width", trim = TRUE, fill = "grey85", color = "black", linewidth = 0.2) +
  geom_boxplot(width = 0.12, outlier.size = 0.15, linewidth = 0.2, fill = "white") +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = NULL, y = "Mitochondrial reads") +
  theme_publication() +
  axis_x_week

p_gene_umi <- ggplot(scatter_metadata, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(size = 0.12, alpha = 0.18, color = "black") +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = "UMIs per nucleus", y = "Genes per nucleus") +
  theme_publication()

p_mt_umi <- ggplot(scatter_metadata, aes(x = nCount_RNA, y = percent.mt)) +
  geom_point(size = 0.12, alpha = 0.18, color = "black") +
  scale_x_continuous(labels = label_comma()) +
  scale_y_continuous(labels = function(x) paste0(x, "%")) +
  labs(x = "UMIs per nucleus", y = "Mitochondrial reads") +
  theme_publication()

p_counts <- metadata |>
  count(week_label, name = "n_cells") |>
  ggplot(aes(x = week_label, y = n_cells)) +
  geom_col(fill = "grey70", color = "black", linewidth = 0.2) +
  scale_y_continuous(labels = label_comma()) +
  labs(x = NULL, y = "Nuclei") +
  theme_publication() +
  axis_x_week
qc_figure <- (p_umi | p_genes | p_mt) /
  (p_gene_umi | p_mt_umi | p_counts) +
  plot_layout(guides = "collect") &
  theme_publication()

save_publication_plot(
  qc_figure,
  file.path(out_figure_dir, "snRNAseq_QC_publication_figure.png"),
  width = 7.2,
  height = 4.2,
  dpi = 600
)

summary_lines <- c(
  paste("Input object:", canonical_rds),
  paste("Cells:", ncol(seu)),
  paste("Genes:", nrow(seu)),
  paste("Gestational weeks:", paste(week_levels, collapse = ", ")),
  paste("Violin cells per week cap:", 5000),
  paste("Scatter cell cap:", 120000),
  paste("Figure PDF:", file.path(out_figure_dir, "snRNAseq_QC_publication_figure.pdf")),
  paste("Figure PNG:", file.path(out_figure_dir, "snRNAseq_QC_publication_figure.png"))
)

writeLines(summary_lines, file.path(out_log_dir, "snRNAseq_QC_figure_summary.txt"))

sink(file.path(out_log_dir, "snRNAseq_QC_figure_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Done.\n")
cat("QC figure written to: ", out_figure_dir, "\n", sep = "")
