#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(dplyr)
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


source(file.path(project_root, "config", "labels_colors.R"))


manifest_file <- if (exists("xenium_manifest", inherits = FALSE)) {
  xenium_manifest
} else {
  file.path(project_root, "config", "xenium_samples.csv")
}

snrna_annotated_rds <- if (exists("fetal_annotated_rds", inherits = FALSE)) {
  fetal_annotated_rds
} else {
  file.path(
    results_base,
    "snRNAseq_annotated_object",
    "objects",
    "fetal_ovary_snRNAseq_canonical_annotated.rds"
  )
}

results_root <- file.path(results_base, "xenium_annotated_object")
out_object_dir <- file.path(results_root, "objects")
out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

out_rds <- file.path(out_object_dir, "fetal_ovary_xenium_annotated.rds")

# -------------------------------------------------------------------------
# Label-transfer settings.
# -------------------------------------------------------------------------
celltype_score_thresh <- 0.60
use_score_gap_rejection <- TRUE
celltype_min_gap <- 0.10

subcluster_score_thresh <- 0.50

min_query_cells_per_type <- 200
min_ref_cells_per_type <- 200

markers_per_type_start <- 200
min_stage1_features <- 600
max_stage1_features <- 3000

stage2_features_per_type <- 2500
stage2_min_features <- 300

transfer_dims <- 1:30
major_cell_type_levels <- c(
  "germ",
  "degenerated",
  "granulosa",
  "stroma",
  "endothelial",
  "mural",
  "immune",
  "erythroid"
)

required_xenium_files <- c(
  "cell_feature_matrix.h5",
  "cells.csv.gz",
  "cell_boundaries.csv.gz"
)

recommended_xenium_files <- c(
  "nucleus_boundaries.csv.gz",
  "transcripts.csv.gz",
  "transcripts.parquet",
  "morphology_focus.ome.tif",
  "experiment.xenium",
  "gene_panel.json",
  "analysis_summary.html"
)

# -------------------------------------------------------------------------
# Helpers.
# -------------------------------------------------------------------------
normalize_include <- function(x) {
  tolower(trimws(as.character(x))) %in% c("true", "t", "1", "yes", "y")
}

get_layer_matrix <- function(obj, assay, layer, slot_fallback = NULL) {
  mat <- tryCatch(
    SeuratObject::LayerData(obj, assay = assay, layer = layer),
    error = function(e) NULL
  )

  if (is.null(mat) && !is.null(slot_fallback)) {
    mat <- tryCatch(
      GetAssayData(obj, assay = assay, slot = slot_fallback),
      error = function(e) NULL
    )
  }

  if (is.null(mat)) {
    stop("Could not retrieve assay layer: ", assay, "/", layer)
  }

  if (!inherits(mat, "dgCMatrix")) {
    mat <- as(mat, "dgCMatrix")
  }

  mat
}

is_bad_feature <- function(gene) {
  grepl("^RPL", gene) |
    grepl("^RPS", gene) |
    grepl("^MT-", gene)
}

make_prediction_gap <- function(metadata, label_levels) {
  score_cols <- intersect(label_levels, colnames(metadata))

  if (length(score_cols) < 2) {
    prefixed <- paste0("prediction.score.", label_levels)
    score_cols <- intersect(prefixed, colnames(metadata))
  }

  if (length(score_cols) < 2) {
    return(rep(NA_real_, nrow(metadata)))
  }

  score_mat <- as.matrix(metadata[, score_cols, drop = FALSE])
  top1 <- apply(score_mat, 1, max, na.rm = TRUE)
  top2 <- apply(score_mat, 1, function(v) sort(v, decreasing = TRUE)[2])

  top1 - top2
}
top_variable_genes_sparse <- function(mat, features, n = 2000) {
  features <- intersect(features, rownames(mat))
  if (length(features) == 0) return(character(0))

  m <- mat[features, , drop = FALSE]

  mu <- Matrix::rowMeans(m)
  mu2 <- Matrix::rowMeans(m^2)
  variance <- mu2 - mu^2
  variance[is.na(variance)] <- 0

  features[order(variance, decreasing = TRUE)][seq_len(min(n, length(features)))]
}

map_final_annotation <- function(fine_subcluster, annotation_map) {
  out <- unname(annotation_map[as.character(fine_subcluster)])

  unresolved <- is.na(out) & !is.na(fine_subcluster)
  out[unresolved] <- as.character(fine_subcluster)[unresolved]
  out[is.na(out)] <- NA_character_

  out
}

map_to_full <- function(values, source_cells, target_cells) {
  unname(values[match(target_cells, source_cells)])
}

check_xenium_directory <- function(xenium_dir) {
  stopifnot(dir.exists(xenium_dir))

  required <- tibble(
    file_or_directory = required_xenium_files,
    required_for_script = TRUE,
    recommended_for_public_archive = TRUE,
    exists = file.exists(file.path(xenium_dir, required_xenium_files)),
    path = file.path(xenium_dir, required_xenium_files)
  )

  missing_required <- required |> filter(!exists)

  if (nrow(missing_required) > 0) {
    stop(
      "Xenium directory is missing required files:\n  ",
      paste(missing_required$file_or_directory, collapse = "\n  "),
      "\nDirectory checked:\n  ",
      xenium_dir
    )
  }

  optional <- tibble(
    file_or_directory = recommended_xenium_files,
    required_for_script = FALSE,
    recommended_for_public_archive = TRUE,
    exists = file.exists(file.path(xenium_dir, recommended_xenium_files)),
    path = file.path(xenium_dir, recommended_xenium_files)
  )

  morphology_focus_dir <- tibble(
    file_or_directory = "morphology_focus/",
    required_for_script = FALSE,
    recommended_for_public_archive = TRUE,
    exists = dir.exists(file.path(xenium_dir, "morphology_focus")),
    path = file.path(xenium_dir, "morphology_focus")
  )

  bind_rows(required, optional, morphology_focus_dir)
}

load_xenium_object <- function(xenium_dir, sample_id) {
  load_args <- list(
    data.dir = xenium_dir,
    fov = sample_id,
    assay = "Xenium"
  )
  lx_formals <- names(formals(LoadXenium))

  if ("molecule.coordinates" %in% lx_formals) {
    load_args$molecule.coordinates <- FALSE
  }

  if ("cell.centroids" %in% lx_formals) {
    load_args$cell.centroids <- TRUE
  }

  if ("segmentations" %in% lx_formals) {
    load_args$segmentations <- "cell"
  }

  message("Loading Xenium directory:\n  ", xenium_dir)
  do.call(LoadXenium, load_args)
}

# -------------------------------------------------------------------------
# Inputs.
# -------------------------------------------------------------------------
if (!file.exists(manifest_file)) {
  stop(
    "Missing local Xenium manifest:\n  ",
    manifest_file,
    "\nCreate it from config/xenium_samples_example.csv and set xenium_dir."
  )
}

manifest <- read_csv(manifest_file, show_col_types = FALSE)

required_manifest_cols <- c("sample_id", "gestational_week", "xenium_dir", "include")
missing_manifest_cols <- setdiff(required_manifest_cols, colnames(manifest))

if (length(missing_manifest_cols) > 0) {
  stop("Manifest missing columns: ", paste(missing_manifest_cols, collapse = ", "))
}

manifest <- manifest |>
  mutate(
    include = normalize_include(include),
    sample_id = as.character(sample_id),
    gestational_week = as.character(gestational_week),
    xenium_dir = as.character(xenium_dir)
  ) |>
  filter(include)

if (nrow(manifest) != 1) {
  stop("Expected exactly one included Xenium sample; found: ", nrow(manifest))
}

sample_id <- manifest$sample_id[[1]]
gestational_week <- manifest$gestational_week[[1]]
xenium_dir <- manifest$xenium_dir[[1]]
source_accession <- if ("source_accession" %in% colnames(manifest)) {
  as.character(manifest$source_accession[[1]])
} else {
  NA_character_
}

if (!dir.exists(xenium_dir)) {
  stop("Xenium directory does not exist:\n  ", xenium_dir)
}

file_audit <- check_xenium_directory(xenium_dir)

write_csv(
  file_audit,
  file.path(out_table_dir, "xenium_input_file_audit.csv")
)

message("Loading annotated snRNA-seq reference:\n  ", snrna_annotated_rds)
stopifnot(file.exists(snrna_annotated_rds))

ref <- readRDS(snrna_annotated_rds)
stopifnot(inherits(ref, "Seurat"))

required_ref_cols <- c("major_cell_type", "fine_subcluster", "final_annotation")
missing_ref_cols <- setdiff(required_ref_cols, colnames(ref@meta.data))

if (length(missing_ref_cols) > 0) {
  stop("Reference object missing metadata columns: ", paste(missing_ref_cols, collapse = ", "))
}

ref@meta.data <- ref@meta.data |>
  mutate(
    major_cell_type = recode(as.character(major_cell_type), quiescent = "degenerated"),
    fine_subcluster = as.character(fine_subcluster),
    final_annotation = as.character(final_annotation)
  )

ref$label_transfer_major_cell_type <- factor(
  ref$major_cell_type,
  levels = major_cell_type_levels
)
ref$label_transfer_fine_subcluster <- factor(ref$fine_subcluster)

annotation_map_df <- ref@meta.data |>
  distinct(fine_subcluster, final_annotation) |>
  arrange(fine_subcluster)

annotation_map <- setNames(
  annotation_map_df$final_annotation,
  annotation_map_df$fine_subcluster
)

fine_to_major_df <- ref@meta.data |>
  distinct(fine_subcluster, major_cell_type) |>
  arrange(fine_subcluster)

fine_to_major <- setNames(
  fine_to_major_df$major_cell_type,
  fine_to_major_df$fine_subcluster
)

message("Reference major cell-type counts:")
print(as.data.frame(sort(table(ref$label_transfer_major_cell_type), decreasing = TRUE)))

xenium_full <- load_xenium_object(xenium_dir, sample_id = sample_id)
stopifnot(inherits(xenium_full, "Seurat"))

xenium_full$sample_id <- sample_id
xenium_full$gestational_week <- gestational_week
xenium_full$dataset <- "xenium"
xenium_full$source_accession <- source_accession

xenium <- xenium_full

if (length(xenium@images) > 0) {
  message("Dropping images from mapping copy only; full Xenium object keeps images/FOV.")
  xenium@images <- list()
}

ref_assay <- if ("RNA" %in% Assays(ref)) "RNA" else DefaultAssay(ref)
xen_assay <- if ("Xenium" %in% Assays(xenium)) "Xenium" else DefaultAssay(xenium)

DefaultAssay(ref) <- ref_assay
DefaultAssay(xenium) <- xen_assay
DefaultAssay(xenium_full) <- xen_assay

message("Reference assay: ", ref_assay)
message("Xenium assay: ", xen_assay)
message("Reference dimensions: ", nrow(ref), " genes x ", ncol(ref), " cells")
message("Xenium dimensions: ", nrow(xenium), " genes x ", ncol(xenium), " cells")

# -------------------------------------------------------------------------
# Stage 1: major cell-type transfer.
# -------------------------------------------------------------------------
message("Normalizing reference and Xenium objects.")

ref <- NormalizeData(ref, verbose = FALSE)
xenium <- NormalizeData(xenium, verbose = FALSE)

common_genes <- intersect(rownames(ref), rownames(xenium))

if (length(common_genes) < 200) {
  stop("Too few shared genes between reference and Xenium: ", length(common_genes))
}

message("Shared genes: ", length(common_genes))
message("Stage 1: finding reference major cell-type markers.")

Idents(ref) <- ref$label_transfer_major_cell_type

markers <- FindAllMarkers(
  ref,
  only.pos = TRUE,
  min.pct = 0.10,
  logfc.threshold = 0.25,
  return.thresh = 0.05
)
markers <- markers |>
  filter(!is.na(gene), gene %in% common_genes, !is_bad_feature(gene))

get_stage1_features <- function(markers_per_type_local) {
  markers |>
    group_by(cluster) |>
    arrange(desc(avg_log2FC), desc(pct.1 - pct.2), .by_group = TRUE) |>
    slice_head(n = markers_per_type_local) |>
    ungroup() |>
    pull(gene) |>
    unique() |>
    intersect(common_genes) |>
    (\(x) x[!is_bad_feature(x)])()
}

markers_per_type <- markers_per_type_start
features_stage1 <- get_stage1_features(markers_per_type)

while (length(features_stage1) < min_stage1_features && markers_per_type < 600) {
  markers_per_type <- markers_per_type + 100
  message("Stage 1 features too few; expanding markers_per_type to ", markers_per_type)
  features_stage1 <- get_stage1_features(markers_per_type)
}

if (length(features_stage1) > max_stage1_features) {
  features_stage1 <- features_stage1[seq_len(max_stage1_features)]
}

if (length(features_stage1) < 200) {
  stop("Too few stage 1 transfer features after filtering: ", length(features_stage1))
}

message("Stage 1 transfer features: ", length(features_stage1))

message("Stage 1: transferring major cell type.")
anchors_major <- FindTransferAnchors(
  reference = ref,
  query = xenium,
  reference.assay = ref_assay,
  query.assay = xen_assay,
  features = features_stage1,
  reduction = "pcaproject",
  dims = transfer_dims,
  verbose = TRUE
)

pred_major <- TransferData(
  anchorset = anchors_major,
  refdata = ref$label_transfer_major_cell_type,
  dims = transfer_dims,
  verbose = TRUE
)

xenium <- AddMetaData(xenium, pred_major)

xenium$pred_major_cell_type <- as.character(xenium$predicted.id)
xenium$pred_major_cell_type_score <- xenium$prediction.score.max
xenium$pred_major_cell_type_gap <- make_prediction_gap(
  xenium@meta.data,
  label_levels = levels(ref$label_transfer_major_cell_type)
)

xenium$pred_major_cell_type_filtered <- xenium$pred_major_cell_type
xenium$pred_major_cell_type_filtered[
  is.na(xenium$pred_major_cell_type_score) |
    xenium$pred_major_cell_type_score < celltype_score_thresh
] <- "Ambiguous"

if (use_score_gap_rejection && !all(is.na(xenium$pred_major_cell_type_gap))) {
  xenium$pred_major_cell_type_filtered[
    xenium$pred_major_cell_type_gap < celltype_min_gap
  ] <- "Ambiguous"
}

stage1_counts <- as.data.frame(
  sort(table(xenium$pred_major_cell_type_filtered), decreasing = TRUE)
)
colnames(stage1_counts) <- c("label", "cells")
stage1_counts$annotation_level <- "major_cell_type_filtered"

# -------------------------------------------------------------------------
# One-step fine-subcluster transfer.
# -------------------------------------------------------------------------
message("One-step transfer: transferring fine_subcluster.")

anchors_fine_1step <- FindTransferAnchors(
  reference = ref,
  query = xenium,
  reference.assay = ref_assay,
  query.assay = xen_assay,
  features = features_stage1,
  reduction = "pcaproject",
  dims = transfer_dims,
  verbose = TRUE
)
pred_fine_1step <- TransferData(
  anchorset = anchors_fine_1step,
  refdata = ref$label_transfer_fine_subcluster,
  dims = transfer_dims,
  verbose = TRUE
)

xenium <- AddMetaData(xenium, pred_fine_1step)

xenium$pred_fine_subcluster_1step <- as.character(xenium$predicted.id)
xenium$pred_fine_subcluster_1step_score <- xenium$prediction.score.max
xenium$pred_final_annotation_1step <- map_final_annotation(
  xenium$pred_fine_subcluster_1step,
  annotation_map
)

# -------------------------------------------------------------------------
# Stage 2: within-major-cell-type fine-subcluster transfer.
# -------------------------------------------------------------------------
message("Creating slim mapping objects for stage 2.")

ref_map <- DietSeurat(
  ref,
  assays = ref_assay,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  dimreducs = NULL,
  graphs = NULL
)

xenium_map <- DietSeurat(
  xenium,
  assays = xen_assay,
  counts = TRUE,
  data = TRUE,
  scale.data = FALSE,
  dimreducs = NULL,
  graphs = NULL
)

DefaultAssay(ref_map) <- ref_assay
DefaultAssay(xenium_map) <- xen_assay

validObject(ref_map)
validObject(xenium_map)

ref_counts_all <- get_layer_matrix(ref_map, assay = ref_assay, layer = "counts", slot_fallback = "counts")
ref_data_all <- get_layer_matrix(ref_map, assay = ref_assay, layer = "data", slot_fallback = "data")
xen_counts_all <- get_layer_matrix(xenium_map, assay = xen_assay, layer = "counts", slot_fallback = "counts")
xen_data_all <- get_layer_matrix(xenium_map, assay = xen_assay, layer = "data", slot_fallback = "data")

stopifnot(identical(colnames(ref_counts_all), colnames(ref_map)))
stopifnot(identical(colnames(ref_data_all), colnames(ref_map)))
stopifnot(identical(colnames(xen_counts_all), colnames(xenium_map)))
stopifnot(identical(colnames(xen_data_all), colnames(xenium_map)))

ref_fine_vec <- setNames(as.character(ref_map$fine_subcluster), colnames(ref_map))

xenium_map$pred_fine_subcluster_2stage <- NA_character_
xenium_map$pred_fine_subcluster_2stage_score <- NA_real_

stage2_summary_rows <- list()

cell_types_to_run <- intersect(
  major_cell_type_levels,
  unique(as.character(ref_map$label_transfer_major_cell_type))
)
message("Stage 2: within-major-cell-type fine_subcluster transfer.")

for (ct in cell_types_to_run) {
  message("---- Major cell type: ", ct)

  q_cells <- colnames(xenium_map)[xenium_map$pred_major_cell_type_filtered == ct]
  r_cells <- colnames(ref_map)[ref_map$label_transfer_major_cell_type == ct]

  n_q <- length(q_cells)
  n_r <- length(r_cells)

  message("     reference cells: ", n_r, "; query cells: ", n_q)

  status <- "ok"

  if (n_q < min_query_cells_per_type) {
    message("     skipping: too few query cells.")
    status <- "skip_too_few_query_cells"
  }

  if (n_r < min_ref_cells_per_type) {
    message("     skipping: too few reference cells.")
    status <- "skip_too_few_reference_cells"
  }

  if (status != "ok") {
    stage2_summary_rows[[length(stage2_summary_rows) + 1]] <- tibble(
      major_cell_type = ct,
      reference_cells = n_r,
      query_cells = n_q,
      features_used = NA_integer_,
      status = status
    )
    next
  }

  ref_counts <- ref_counts_all[, r_cells, drop = FALSE]
  ref_data <- ref_data_all[, r_cells, drop = FALSE]
  xen_counts <- xen_counts_all[, q_cells, drop = FALSE]
  xen_data <- xen_data_all[, q_cells, drop = FALSE]

  refdata_ct <- factor(ref_fine_vec[r_cells])

  ref_ct_obj <- CreateSeuratObject(
    counts = ref_counts,
    assay = ref_assay,
    meta.data = ref_map@meta.data[r_cells, , drop = FALSE]
  )

  ref_ct_obj <- SetAssayData(
    ref_ct_obj,
    assay = ref_assay,
    layer = "data",
    new.data = ref_data
  )

  xen_ct_obj <- CreateSeuratObject(
    counts = xen_counts,
    assay = xen_assay,
    meta.data = xenium_map@meta.data[q_cells, , drop = FALSE]
  )

  xen_ct_obj <- SetAssayData(
    xen_ct_obj,
    assay = xen_assay,
    layer = "data",
    new.data = xen_data
  )

  DefaultAssay(ref_ct_obj) <- ref_assay
  DefaultAssay(xen_ct_obj) <- xen_assay

  common_ct <- intersect(rownames(ref_ct_obj), rownames(xen_ct_obj))

  feats_ct <- top_variable_genes_sparse(
    ref_data[common_ct, , drop = FALSE],
    common_ct,
    n = stage2_features_per_type
  )

  feats_ct <- feats_ct[!is_bad_feature(feats_ct)]

  if (length(feats_ct) < stage2_min_features) {
    feats_ct <- common_ct[!is_bad_feature(common_ct)]
  }

  if (length(feats_ct) > stage2_features_per_type) {
    feats_ct <- feats_ct[seq_len(stage2_features_per_type)]
  }
  message("     features: ", length(feats_ct))

  if (length(feats_ct) < 100) {
    warning("Too few features for ", ct, "; skipping stage 2.")
    stage2_summary_rows[[length(stage2_summary_rows) + 1]] <- tibble(
      major_cell_type = ct,
      reference_cells = n_r,
      query_cells = n_q,
      features_used = length(feats_ct),
      status = "skip_too_few_features"
    )
    next
  }

  anchors_ct <- FindTransferAnchors(
    reference = ref_ct_obj,
    query = xen_ct_obj,
    reference.assay = ref_assay,
    query.assay = xen_assay,
    features = feats_ct,
    reduction = "pcaproject",
    dims = transfer_dims,
    verbose = FALSE
  )

  pred_ct <- TransferData(
    anchorset = anchors_ct,
    refdata = refdata_ct,
    dims = transfer_dims,
    verbose = FALSE
  )

  xenium_map$pred_fine_subcluster_2stage[q_cells] <- as.character(pred_ct$predicted.id)
  xenium_map$pred_fine_subcluster_2stage_score[q_cells] <- pred_ct$prediction.score.max

  stage2_summary_rows[[length(stage2_summary_rows) + 1]] <- tibble(
    major_cell_type = ct,
    reference_cells = n_r,
    query_cells = n_q,
    features_used = length(feats_ct),
    status = "ok"
  )
}

stage2_summary <- bind_rows(stage2_summary_rows)

xenium_map$pred_fine_subcluster_2stage_filtered <- ifelse(
  !is.na(xenium_map$pred_fine_subcluster_2stage) &
    xenium_map$pred_fine_subcluster_2stage_score >= subcluster_score_thresh,
  xenium_map$pred_fine_subcluster_2stage,
  "Ambiguous_subcluster"
)

stage2_counts <- as.data.frame(
  sort(table(xenium_map$pred_fine_subcluster_2stage_filtered), decreasing = TRUE)
)
colnames(stage2_counts) <- c("label", "cells")
stage2_counts$annotation_level <- "fine_subcluster_2stage_filtered"

# -------------------------------------------------------------------------
# Final hybrid annotation.
# -------------------------------------------------------------------------
xenium$pred_fine_subcluster_2stage <- map_to_full(
  xenium_map$pred_fine_subcluster_2stage,
  colnames(xenium_map),
  colnames(xenium)
)
xenium$pred_fine_subcluster_2stage_score <- map_to_full(
  xenium_map$pred_fine_subcluster_2stage_score,
  colnames(xenium_map),
  colnames(xenium)
)

xenium$pred_fine_subcluster_2stage_filtered <- map_to_full(
  xenium_map$pred_fine_subcluster_2stage_filtered,
  colnames(xenium_map),
  colnames(xenium)
)

xenium$pred_final_annotation_2stage_filtered <- map_final_annotation(
  xenium$pred_fine_subcluster_2stage_filtered,
  annotation_map
)

one_step_major <- unname(fine_to_major[xenium$pred_fine_subcluster_1step])

replace_with_1step <- xenium$pred_fine_subcluster_2stage_filtered == "Ambiguous_subcluster" &
  !is.na(xenium$pred_fine_subcluster_1step) &
  !(one_step_major %in% c("germ", "degenerated"))

xenium$pred_fine_subcluster_hybrid <- xenium$pred_fine_subcluster_2stage_filtered
xenium$pred_fine_subcluster_hybrid[replace_with_1step] <- xenium$pred_fine_subcluster_1step[
  replace_with_1step
]

xenium$pred_final_annotation_hybrid <- map_final_annotation(
  xenium$pred_fine_subcluster_hybrid,
  annotation_map
)

xenium$pred_final_annotation_hybrid[
  xenium$pred_fine_subcluster_hybrid == "Ambiguous_subcluster"
] <- "Ambiguous_subcluster"

metadata_cols_to_copy <- c(
  "sample_id",
  "gestational_week",
  "dataset",
  "source_accession",
  "pred_major_cell_type",
  "pred_major_cell_type_score",
  "pred_major_cell_type_gap",
  "pred_major_cell_type_filtered",
  "pred_fine_subcluster_1step",
  "pred_fine_subcluster_1step_score",
  "pred_final_annotation_1step",
  "pred_fine_subcluster_2stage",
  "pred_fine_subcluster_2stage_score",
  "pred_fine_subcluster_2stage_filtered",
  "pred_final_annotation_2stage_filtered",
  "pred_fine_subcluster_hybrid",
  "pred_final_annotation_hybrid"
)

for (col in metadata_cols_to_copy) {
  xenium_full[[col]] <- map_to_full(
    xenium@meta.data[[col]],
    colnames(xenium),
    colnames(xenium_full)
  )
}

xenium_full@misc$xenium_annotation <- list(
  script = "00_build_objects/05_build_xenium_annotated_object.R",
  sample_id = sample_id,
  gestational_week = gestational_week,
  source_accession = source_accession,
  reference_object = "snRNAseq_annotated_object/objects/fetal_ovary_snRNAseq_canonical_annotated.rds",
  celltype_score_thresh = celltype_score_thresh,
  use_score_gap_rejection = use_score_gap_rejection,
  celltype_min_gap = celltype_min_gap,
  subcluster_score_thresh = subcluster_score_thresh,
  stage1_features = length(features_stage1)
)

saveRDS(xenium_full, out_rds)
# -------------------------------------------------------------------------
# Minimal outputs.
# -------------------------------------------------------------------------
annotation_table <- xenium_full@meta.data |>
  rownames_to_column("cell_id") |>
  transmute(
    cell_id = cell_id,
    sample_id = as.character(sample_id),
    gestational_week = as.character(gestational_week),
    dataset = as.character(dataset),
    source_accession = as.character(source_accession),
    pred_major_cell_type_filtered = as.character(pred_major_cell_type_filtered),
    pred_major_cell_type_score = as.numeric(pred_major_cell_type_score),
    pred_major_cell_type_gap = as.numeric(pred_major_cell_type_gap),
    pred_fine_subcluster_1step = as.character(pred_fine_subcluster_1step),
    pred_fine_subcluster_1step_score = as.numeric(pred_fine_subcluster_1step_score),
    pred_fine_subcluster_2stage_filtered = as.character(pred_fine_subcluster_2stage_filtered),
    pred_fine_subcluster_2stage_score = as.numeric(pred_fine_subcluster_2stage_score),
    pred_fine_subcluster_hybrid = as.character(pred_fine_subcluster_hybrid),
    pred_final_annotation_hybrid = as.character(pred_final_annotation_hybrid)
  )

write_csv(
  annotation_table,
  file.path(out_table_dir, "xenium_cell_annotations.csv")
)

hybrid_counts <- annotation_table |>
  count(pred_fine_subcluster_hybrid, pred_final_annotation_hybrid, sort = TRUE, name = "cells") |>
  mutate(annotation_level = "fine_subcluster_hybrid") |>
  rename(label = pred_fine_subcluster_hybrid, final_annotation = pred_final_annotation_hybrid) |>
  select(annotation_level, label, final_annotation, cells)

stage1_counts_out <- stage1_counts |>
  mutate(final_annotation = NA_character_) |>
  select(annotation_level, label, final_annotation, cells)

stage2_counts_out <- stage2_counts |>
  mutate(final_annotation = map_final_annotation(label, annotation_map)) |>
  select(annotation_level, label, final_annotation, cells)

annotation_counts <- bind_rows(
  stage1_counts_out,
  stage2_counts_out,
  hybrid_counts
)

write_csv(
  annotation_counts,
  file.path(out_table_dir, "xenium_annotation_counts.csv")
)

write_csv(
  stage2_summary,
  file.path(out_table_dir, "xenium_stage2_transfer_summary.csv")
)
summary_table <- tibble(
  sample_id = sample_id,
  gestational_week = gestational_week,
  source_accession = source_accession,
  xenium_dir = xenium_dir,
  snrna_annotated_rds = snrna_annotated_rds,
  output_rds = out_rds,
  reference_cells = ncol(ref),
  reference_genes = nrow(ref),
  xenium_cells = ncol(xenium_full),
  xenium_genes = nrow(xenium_full),
  shared_genes = length(common_genes),
  stage1_transfer_features = length(features_stage1),
  celltype_score_thresh = celltype_score_thresh,
  use_score_gap_rejection = use_score_gap_rejection,
  celltype_min_gap = celltype_min_gap,
  subcluster_score_thresh = subcluster_score_thresh,
  min_query_cells_per_type = min_query_cells_per_type,
  min_ref_cells_per_type = min_ref_cells_per_type
)

write_csv(
  summary_table,
  file.path(out_table_dir, "xenium_build_summary.csv")
)

summary_lines <- c(
  paste("Input snRNA-seq annotated object:", snrna_annotated_rds),
  paste("Input Xenium manifest:", manifest_file),
  paste("Sample ID:", sample_id),
  paste("Gestational week:", gestational_week),
  paste("Source accession:", source_accession),
  paste("Xenium directory:", xenium_dir),
  paste("Output object:", out_rds),
  paste("Reference assay:", ref_assay),
  paste("Xenium assay:", xen_assay),
  paste("Reference cells:", ncol(ref)),
  paste("Reference genes:", nrow(ref)),
  paste("Xenium cells:", ncol(xenium_full)),
  paste("Xenium genes:", nrow(xenium_full)),
  paste("Shared genes:", length(common_genes)),
  paste("Stage 1 transfer features:", length(features_stage1)),
  paste("Stage 1 cell-type score threshold:", celltype_score_thresh),
  paste("Stage 1 score-gap rejection:", use_score_gap_rejection),
  paste("Stage 1 minimum score gap:", celltype_min_gap),
  paste("Stage 2 subcluster score threshold:", subcluster_score_thresh),
  "",
  "Generated outputs:",
  paste("  ", out_rds),
  paste("  ", file.path(out_table_dir, "xenium_input_file_audit.csv")),
  paste("  ", file.path(out_table_dir, "xenium_cell_annotations.csv")),
  paste("  ", file.path(out_table_dir, "xenium_annotation_counts.csv")),
  paste("  ", file.path(out_table_dir, "xenium_stage2_transfer_summary.csv")),
  paste("  ", file.path(out_table_dir, "xenium_build_summary.csv")),
  "",
  "Stage 2 transfer summary:",
  paste(capture.output(print(as.data.frame(stage2_summary))), collapse = "\n")
)

writeLines(
  summary_lines,
  file.path(out_log_dir, "xenium_annotated_object_summary.txt")
)
sink(file.path(out_log_dir, "xenium_annotated_object_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone.\n")
cat("Annotated Xenium object written to:\n  ", out_rds, "\n", sep = "")
cat("Tables written to:\n  ", out_table_dir, "\n", sep = "")
cat("Logs written to:\n  ", out_log_dir, "\n", sep = "")
