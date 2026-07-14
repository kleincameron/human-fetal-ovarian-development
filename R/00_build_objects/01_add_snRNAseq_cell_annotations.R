suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

set.seed(42)

project_root <- "/home/liyan/liyan/Final/github_code_for_publication"
input_root <- "/home/liyan/liyan/Final/github_code_for_publication_results/snRNAseq_initial_processing"
results_root <- "/home/liyan/liyan/Final/github_code_for_publication_results/snRNAseq_annotated_object"

source(file.path(project_root, "config", "labels_colors.R"))

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

out_rds <- file.path(
  out_object_dir,
  "fetal_ovary_snRNAseq_canonical_annotated.rds"
)

cat("Loading canonical snRNA-seq object.\n")
stopifnot(file.exists(canonical_rds))
seu <- readRDS(canonical_rds)
stopifnot(inherits(seu, "Seurat"))
cat("Loading cell annotation metadata.\n")
stopifnot(file.exists(annotation_csv))
annotations <- read_csv(annotation_csv, show_col_types = FALSE)

required_annotation_cols <- c(
  "cell_barcode",
  "sample_id",
  "gestational_week",
  "seurat_cluster",
  "major_cell_type",
  "fine_subcluster",
  "final_annotation",
  "annotation_version",
  "annotation_source"
)

missing_annotation_cols <- setdiff(required_annotation_cols, colnames(annotations))
if (length(missing_annotation_cols) > 0) {
  stop(
    "Missing required annotation columns: ",
    paste(missing_annotation_cols, collapse = ", ")
  )
}

required_object_cols <- c("sample_id", "gestational_week", "seurat_clusters")
missing_object_cols <- setdiff(required_object_cols, colnames(seu@meta.data))
if (length(missing_object_cols) > 0) {
  stop(
    "Missing required object metadata columns: ",
    paste(missing_object_cols, collapse = ", ")
  )
}

if (anyDuplicated(annotations$cell_barcode) > 0) {
  duplicated_barcodes <- annotations$cell_barcode[duplicated(annotations$cell_barcode)]
  stop(
    "Duplicated cell barcodes in annotation file. First duplicated barcode: ",
    duplicated_barcodes[1]
  )
}

if (nrow(annotations) != ncol(seu)) {
  stop(
    "Annotation row count does not match canonical object cell count. ",
    "Annotation rows: ", nrow(annotations), "; object cells: ", ncol(seu)
  )
}
missing_from_annotations <- setdiff(colnames(seu), annotations$cell_barcode)
extra_in_annotations <- setdiff(annotations$cell_barcode, colnames(seu))

if (length(missing_from_annotations) > 0) {
  stop(
    "Some canonical object cells are missing from the annotation file. Example: ",
    missing_from_annotations[1]
  )
}

if (length(extra_in_annotations) > 0) {
  stop(
    "Some annotation-file cells are not present in the canonical object. Example: ",
    extra_in_annotations[1]
  )
}

annotations <- annotations[match(colnames(seu), annotations$cell_barcode), ]

if (!identical(annotations$cell_barcode, colnames(seu))) {
  stop("Failed to reorder annotations to match canonical object cell order.")
}

sample_match <- as.character(annotations$sample_id) == as.character(seu$sample_id)
if (!all(sample_match)) {
  stop("Sample IDs do not match between annotation file and canonical object.")
}

week_match <- as.numeric(annotations$gestational_week) == as.numeric(seu$gestational_week)
if (!all(week_match)) {
  stop("Gestational weeks do not match between annotation file and canonical object.")
}

cluster_match <- as.character(annotations$seurat_cluster) == as.character(seu$seurat_clusters)
if (!all(cluster_match)) {
  stop("Seurat clusters do not match between annotation file and canonical object.")
}

unknown_major_cell_types <- sort(setdiff(unique(annotations$major_cell_type), celltype_order))
if (length(unknown_major_cell_types) > 0) {
  stop(
    "Unknown major_cell_type values not present in celltype_order:\n",
    paste(unknown_major_cell_types, collapse = "\n")
  )
}
unknown_subclusters <- sort(setdiff(unique(annotations$fine_subcluster), names(subcluster_label_map)))
if (length(unknown_subclusters) > 0) {
  stop(
    "Unknown fine_subcluster values not present in subcluster_label_map:\n",
    paste(unknown_subclusters, collapse = "\n")
  )
}

expected_final_annotations <- unname(subcluster_label_map[annotations$fine_subcluster])
label_match <- as.character(annotations$final_annotation) == expected_final_annotations
if (!all(label_match)) {
  mismatch <- annotations[which(!label_match)[1], ]
  stop(
    "final_annotation does not match config/labels_colors.R for fine_subcluster ",
    mismatch$fine_subcluster,
    ". Annotation file has: ",
    mismatch$final_annotation,
    "; expected: ",
    expected_final_annotations[which(!label_match)[1]]
  )
}

annotation_metadata <- annotations |>
  transmute(
    major_cell_type = factor(major_cell_type, levels = celltype_order),
    fine_subcluster = factor(fine_subcluster, levels = names(subcluster_label_map)),
    final_annotation = factor(final_annotation, levels = unname(subcluster_label_map)),
    annotation_version = as.character(annotation_version),
    annotation_source = as.character(annotation_source)
  ) |>
  as.data.frame()

rownames(annotation_metadata) <- annotations$cell_barcode

if (!identical(rownames(annotation_metadata), colnames(seu))) {
  stop("Annotation metadata row names do not match Seurat object cell names.")
}

seu <- AddMetaData(seu, metadata = annotation_metadata)

major_counts <- seu@meta.data |>
  count(major_cell_type, name = "n_cells") |>
  arrange(factor(major_cell_type, levels = celltype_order))

subcluster_counts <- seu@meta.data |>
  count(fine_subcluster, final_annotation, major_cell_type, name = "n_cells") |>
  arrange(factor(fine_subcluster, levels = names(subcluster_label_map)))

week_major_counts <- seu@meta.data |>
  count(gestational_week, sample_id, major_cell_type, name = "n_cells") |>
  arrange(gestational_week, factor(major_cell_type, levels = celltype_order))

write_csv(
  major_counts,
  file.path(out_table_dir, "snRNAseq_annotated_counts_by_major_cell_type.csv")
)
write_csv(
  subcluster_counts,
  file.path(out_table_dir, "snRNAseq_annotated_counts_by_subcluster.csv")
)

write_csv(
  week_major_counts,
  file.path(out_table_dir, "snRNAseq_annotated_counts_by_week_and_major_cell_type.csv")
)

saveRDS(seu, out_rds)

summary_lines <- c(
  paste("Input canonical object:", canonical_rds),
  paste("Input annotation metadata:", annotation_csv),
  paste("Output annotated object:", out_rds),
  paste("Cells:", ncol(seu)),
  paste("Genes:", nrow(seu)),
  paste("Major cell types:", paste(celltype_order, collapse = ", ")),
  paste("Fine subclusters:", length(unique(seu$fine_subcluster))),
  paste("Annotation version:", paste(unique(seu$annotation_version), collapse = ", ")),
  paste("Annotation source:", paste(unique(seu$annotation_source), collapse = ", "))
)

writeLines(summary_lines, file.path(out_log_dir, "snRNAseq_annotated_object_summary.txt"))

sink(file.path(out_log_dir, "snRNAseq_annotated_object_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Done.\n")
cat("Annotated object written to: ", out_rds, "\n", sep = "")
