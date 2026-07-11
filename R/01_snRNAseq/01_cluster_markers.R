suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
  library(presto)
})

set.seed(42)

project_root <- "/home/liyan/liyan/Final/github_code_for_publication"
input_root <- "/home/liyan/liyan/Final/github_code_for_publication_results/snRNAseq_initial_processing"
results_root <- "/home/liyan/liyan/Final/github_code_for_publication_results/snRNAseq_cluster_markers"

canonical_rds <- file.path(
  input_root,
  "objects",
  "fetal_ovary_snRNAseq_canonical.rds"
)

out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

is_mt_or_ribo <- function(genes) {
  grepl("^MT-", genes) | grepl("^(RP[SL]|MRP[SL])", genes)
}

cat("Loading canonical snRNA-seq object.\n")
stopifnot(file.exists(canonical_rds))
seu <- readRDS(canonical_rds)
stopifnot(inherits(seu, "Seurat"))
stopifnot("seurat_clusters" %in% colnames(seu@meta.data))

DefaultAssay(seu) <- "RNA"

if (!requireNamespace("presto", quietly = TRUE)) {
  stop(
    "The presto package is required for this script. ",
    "Install it in the R environment before running marker analysis."
  )
}

cat("Running presto Wilcoxon marker analysis.\n")
markers <- presto::wilcoxauc(
  X = seu,
  group_by = "seurat_clusters",
  seurat_assay = "RNA",
  assay = "data"
)
markers <- markers %>%
  rename(
    gene = feature,
    cluster = group,
    p_val = pval,
    p_val_adj = padj,
    pct.1 = pct_in,
    pct.2 = pct_out,
    avg_log2FC = logFC
  ) %>%
  mutate(
    cluster = as.character(cluster),
    is_mt_or_ribo = is_mt_or_ribo(gene)
  ) %>%
  arrange(as.numeric(cluster), p_val_adj, desc(avg_log2FC))

markers_filtered <- markers %>%
  filter(
    avg_log2FC > 0.10,
    pct.1 >= 0.10,
    p_val_adj <= 0.05,
    !is_mt_or_ribo
  ) %>%
  arrange(as.numeric(cluster), p_val_adj, desc(avg_log2FC))

write_csv(
  markers,
  file.path(out_table_dir, "snRNAseq_cluster_markers_all_presto.csv")
)

write_csv(
  markers_filtered,
  file.path(out_table_dir, "snRNAseq_cluster_markers_up_filtered_presto.csv")
)

summary_lines <- c(
  paste("Input object:", canonical_rds),
  paste("Cells:", ncol(seu)),
  paste("Genes:", nrow(seu)),
  paste("Clusters:", paste(sort(unique(seu$seurat_clusters)), collapse = ", ")),
  paste("Total marker rows:", nrow(markers)),
  paste("Filtered marker rows:", nrow(markers_filtered)),
  paste("Output directory:", results_root),
  paste("presto version:", as.character(packageVersion("presto")))
)

writeLines(summary_lines, file.path(out_log_dir, "snRNAseq_cluster_markers_summary.txt"))

sink(file.path(out_log_dir, "snRNAseq_cluster_markers_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Done.\n")
cat("Marker tables written to: ", out_table_dir, "\n", sep = "")
