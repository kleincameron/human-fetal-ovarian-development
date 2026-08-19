suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

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
results_root <- file.path(results_base, "snRNAseq_annotated_object")

canonical_rds <- file.path(
  input_root,
  "objects",
  "fetal_ovary_snRNAseq_canonical.rds"
)

annotation_csv <- file.path(
  project_root,
  "metadata",
  "snRNAseq_cell_annotations.csv"
)

out_object_dir <- file.path(results_root, "objects")
out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

annotated_rds <- file.path(
  out_object_dir,
  "fetal_ovary_snRNAseq_canonical_annotated.rds"
)

required_annotation_cols <- c(
  "cell_barcode",
  "sample_id",
  "gestational_week",
  "seurat_cluster",
  "major_cell_type",
  "fine_subcluster",
  "final_annotation",
  "annotation_version",
  "annotation_source",
  "UMAP_1_publication",
  "UMAP_2_publication"
)
cat("Loading canonical snRNA-seq object.\n")
stopifnot(file.exists(canonical_rds))
seu <- readRDS(canonical_rds)
stopifnot(inherits(seu, "Seurat"))

cat("Loading cell annotation metadata.\n")
stopifnot(file.exists(annotation_csv))
annotations <- read_csv(annotation_csv, show_col_types = FALSE)

missing_annotation_cols <- setdiff(required_annotation_cols, colnames(annotations))
if (length(missing_annotation_cols) > 0) {
  stop(
    "Annotation metadata is missing required columns: ",
    paste(missing_annotation_cols, collapse = ", ")
  )
}

annotations <- annotations |>
  mutate(
    cell_barcode = as.character(cell_barcode),
    sample_id = as.character(sample_id),
    gestational_week = as.numeric(gestational_week),
    seurat_cluster = as.character(seurat_cluster),
    major_cell_type = as.character(major_cell_type),
    fine_subcluster = as.character(fine_subcluster),
    final_annotation = as.character(final_annotation),
    annotation_version = as.character(annotation_version),
    annotation_source = as.character(annotation_source),
    UMAP_1_publication = as.numeric(UMAP_1_publication),
    UMAP_2_publication = as.numeric(UMAP_2_publication)
  )

if (anyDuplicated(annotations$cell_barcode) > 0) {
  stop("Annotation metadata contains duplicated cell_barcode values.")
}

canonical_cells <- colnames(seu)

missing_from_annotations <- setdiff(canonical_cells, annotations$cell_barcode)
extra_in_annotations <- setdiff(annotations$cell_barcode, canonical_cells)

if (length(missing_from_annotations) > 0 || length(extra_in_annotations) > 0) {
  cat("Example canonical cells missing from annotation metadata:\n")
  print(head(missing_from_annotations, 10))

  cat("\nExample annotation cells absent from canonical object:\n")
  print(head(extra_in_annotations, 10))

  stop("Canonical object cells and annotation metadata cells do not match exactly.")
}
annotations_ordered <- annotations[match(canonical_cells, annotations$cell_barcode), ]

if (!identical(annotations_ordered$cell_barcode, canonical_cells)) {
  stop("Annotation metadata could not be ordered to match canonical object cells.")
}

if (anyNA(annotations_ordered$major_cell_type)) {
  stop("major_cell_type contains missing values.")
}

if (anyNA(annotations_ordered$fine_subcluster)) {
  stop("fine_subcluster contains missing values.")
}

if (anyNA(annotations_ordered$final_annotation)) {
  stop("final_annotation contains missing values.")
}

if (anyNA(annotations_ordered$UMAP_1_publication) || anyNA(annotations_ordered$UMAP_2_publication)) {
  stop("Publication UMAP coordinates contain missing values.")
}

annotation_metadata <- annotations_ordered |>
  select(
    cell_barcode,
    major_cell_type,
    fine_subcluster,
    final_annotation,
    annotation_version,
    annotation_source,
    UMAP_1_publication,
    UMAP_2_publication
  )

annotation_metadata <- as.data.frame(annotation_metadata)
rownames(annotation_metadata) <- annotation_metadata$cell_barcode
annotation_metadata$cell_barcode <- NULL

seu <- AddMetaData(seu, metadata = annotation_metadata)

if (!all(c("UMAP_1_publication", "UMAP_2_publication") %in% colnames(seu@meta.data))) {
  stop("Publication UMAP coordinates were not added to Seurat metadata.")
}

saveRDS(seu, annotated_rds)

write_csv(
  seu@meta.data |>
    count(major_cell_type, name = "n_cells") |>
    arrange(major_cell_type),
  file.path(out_table_dir, "snRNAseq_annotated_counts_by_major_cell_type.csv")
)
write_csv(
  seu@meta.data |>
    count(fine_subcluster, final_annotation, major_cell_type, name = "n_cells") |>
    arrange(major_cell_type, fine_subcluster),
  file.path(out_table_dir, "snRNAseq_annotated_counts_by_subcluster.csv")
)

write_csv(
  seu@meta.data |>
    count(gestational_week, major_cell_type, name = "n_cells") |>
    arrange(gestational_week, major_cell_type),
  file.path(out_table_dir, "snRNAseq_annotated_counts_by_week_and_major_cell_type.csv")
)

summary_lines <- c(
  paste("Input canonical object:", canonical_rds),
  paste("Input annotation metadata:", annotation_csv),
  paste("Output annotated object:", annotated_rds),
  paste("Cells:", ncol(seu)),
  paste("Genes:", nrow(seu)),
  paste("Major cell types:", paste(unique(annotations_ordered$major_cell_type), collapse = ", ")),
  paste("Fine subclusters:", length(unique(annotations_ordered$fine_subcluster))),
  paste("Annotation version:", paste(unique(annotations_ordered$annotation_version), collapse = ", ")),
  paste("Annotation source:", paste(unique(annotations_ordered$annotation_source), collapse = ", ")),
  paste("Publication UMAP coordinates included:", all(c("UMAP_1_publication", "UMAP_2_publication") %in% colnames(seu@meta.data)))
)

writeLines(summary_lines, file.path(out_log_dir, "snRNAseq_annotated_object_summary.txt"))

sink(file.path(out_log_dir, "snRNAseq_annotated_object_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Done.\n")
cat("Annotated object written to: ", annotated_rds, "\n", sep = "")
