suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(SingleCellExperiment)
  library(scDblFinder)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

project_root <- "/home/liyan/liyan/Final/github_code_for_publication"
results_root <- "/home/liyan/liyan/Final/github_code_for_publication_results/snRNAseq_initial_processing"

source(file.path(project_root, "config", "plotting.R"))
source(file.path(project_root, "config", "labels_colors.R"))

sample_manifest <- file.path(project_root, "config", "snrna_samples.csv")

out_object_dir <- file.path(results_root, "objects")
out_table_dir <- file.path(results_root, "tables")
out_figure_dir <- file.path(results_root, "figures")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_object_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

canonical_rds <- file.path(out_object_dir, "fetal_ovary_snRNAseq_canonical.rds")

read_sample_manifest <- function(path) {
  stopifnot(file.exists(path))

  manifest <- readr::read_csv(path, show_col_types = FALSE)

  required <- c("sample_id", "gestational_week", "matrix_dir", "include")
  missing <- setdiff(required, colnames(manifest))
  if (length(missing) > 0) {
    stop("Missing columns in sample manifest: ", paste(missing, collapse = ", "))
  }

  manifest <- manifest |>
    mutate(
      sample_id = as.character(sample_id),
      gestational_week = as.numeric(gestational_week),
      matrix_dir = as.character(matrix_dir),
      include = as.logical(include)
    ) |>
    filter(include)

  if (nrow(manifest) == 0) {
    stop("No samples marked include=TRUE in sample manifest.")
  }
  duplicated_samples <- manifest$sample_id[duplicated(manifest$sample_id)]
  if (length(duplicated_samples) > 0) {
    stop("Duplicate sample_id values: ", paste(unique(duplicated_samples), collapse = ", "))
  }

  missing_dirs <- manifest$matrix_dir[!dir.exists(manifest$matrix_dir)]
  if (length(missing_dirs) > 0) {
    stop("These matrix directories do not exist:\n", paste(missing_dirs, collapse = "\n"))
  }

  manifest
}

read_one_sample <- function(sample_id, gestational_week, matrix_dir) {
  cat("Loading ", sample_id, "\n", sep = "")

  counts <- Read10X(data.dir = matrix_dir)

  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      counts <- counts[[1]]
    }
  }

  if (!inherits(counts, "dgCMatrix")) {
    counts <- as(counts, "dgCMatrix")
  }

  colnames(counts) <- paste(sample_id, colnames(counts), sep = "_")

  obj <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.cells = 1,
    min.features = 100
  )

  obj$sample_id <- sample_id
  obj$sampleID <- sample_id
  obj$gestational_week <- gestational_week
  obj$Week <- gestational_week

  obj
}
join_layers_if_needed <- function(obj) {
  DefaultAssay(obj) <- "RNA"

  layers <- tryCatch(
    Layers(obj[["RNA"]]),
    error = function(e) character()
  )

  if (length(layers) > 0 && any(grepl("^counts\\.", layers))) {
    obj <- JoinLayers(obj, assay = "RNA")
  }

  obj
}

count_cells <- function(stage, obj) {
  obj@meta.data |>
    count(sample_id, name = "n_cells") |>
    mutate(stage = stage, .before = 1)
}

run_scdblfinder_by_sample <- function(obj, dbr = 0.06) {
  cat("Running scDblFinder per sample.\n")

  all_cells <- colnames(obj)
  dbl_class <- setNames(rep(NA_character_, length(all_cells)), all_cells)
  dbl_score <- setNames(rep(NA_real_, length(all_cells)), all_cells)
  singlet_cells <- character()

  for (sid in sort(unique(obj$sample_id))) {
    cells <- colnames(obj)[obj$sample_id == sid]
    sub <- subset(obj, cells = cells)
    sce <- as.SingleCellExperiment(sub)

    sce <- scDblFinder(sce, dbr = dbr, verbose = FALSE)

    cls <- as.character(SingleCellExperiment::colData(sce)$scDblFinder.class)
    scr <- as.numeric(SingleCellExperiment::colData(sce)$scDblFinder.score)

    names(cls) <- colnames(sce)
    names(scr) <- colnames(sce)

    keep <- cls == "singlet"

    dbl_class[colnames(sce)] <- cls
    dbl_score[colnames(sce)] <- scr
    singlet_cells <- c(singlet_cells, colnames(sce)[keep])

    cat("  ", sid, ": kept ", sum(keep), " singlets of ", length(keep), " cells.\n", sep = "")
  }

  obj$scDblFinder_class <- unname(dbl_class[colnames(obj)])
  obj$scDblFinder_score <- unname(dbl_score[colnames(obj)])

  subset(obj, cells = singlet_cells)
}
add_lineage_scores <- function(obj) {
  gene_sets <- list(
    germ = c("DDX4", "DPPA3", "GDF9", "ZP2", "ZP3", "FIGLA", "SYCP1", "SYCP3"),
    granulosa = c("FOXL2", "AMHR2", "INHBA", "INHBB", "HSD17B1", "WT1"),
    stroma = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "FN1", "COL6A1", "COL6A3", "THY1"),
    endothelial = c("PECAM1", "VWF", "KDR", "KLF2", "CLDN5", "ESAM"),
    mural = c("RGS5", "PDGFRB", "ACTA2", "TAGLN", "MYH11", "NOTCH3"),
    immune = c("PTPRC", "LYZ", "TYROBP", "HLA-DRA", "HLA-DRB1", "CD3D", "MS4A1"),
    ose = c("KRT8", "KRT18", "KRT19", "EPCAM", "WT1", "CLDN3", "CLDN4", "MSLN"),
    neuronal = c("MAP2", "RBFOX3", "NEFL", "NEFM", "STMN2", "TUBB3", "SNAP25")
  )

  present <- rownames(obj)
  gene_sets <- lapply(gene_sets, intersect, present)
  gene_sets <- gene_sets[lengths(gene_sets) > 0]

  for (nm in names(gene_sets)) {
    obj <- AddModuleScore(
      object = obj,
      features = list(gene_sets[[nm]]),
      name = paste0("score_", nm),
      assay = "RNA"
    )
  }

  obj
}

write_lineage_score_summary <- function(obj) {
  score_cols <- grep("^score_.*1$", colnames(obj@meta.data), value = TRUE)

  if (length(score_cols) == 0) {
    return(invisible(NULL))
  }

  thresholds <- sapply(score_cols, function(cn) {
    v <- obj[[cn]][, 1]
    median(v, na.rm = TRUE) + 0.25 * IQR(v, na.rm = TRUE)
  })

  out <- lapply(score_cols, function(cn) {
    lineage <- sub("^score_", "", sub("1$", "", cn))
    is_pos <- obj[[cn]][, 1] > thresholds[[cn]]

    by_sample <- data.frame(
      lineage = lineage,
      sample_id = obj$sample_id,
      positive = is_pos
    ) |>
      group_by(lineage, sample_id) |>
      summarise(pct_pos = mean(positive) * 100, .groups = "drop")

    global <- data.frame(
      lineage = lineage,
      sample_id = "ALL",
      pct_pos = mean(is_pos) * 100
    )

    bind_rows(global, by_sample)
  }) |>
    bind_rows() |>
    arrange(lineage, sample_id)

  write_csv(out, file.path(out_table_dir, "snRNAseq_lineage_score_summary_by_sample.csv"))
}
plot_outputs <- function(obj) {
  qc_plot <- VlnPlot(
    obj,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "sample_id",
    ncol = 3,
    pt.size = 0
  ) &
    theme_publication()

  save_publication_plot(
    qc_plot,
    file.path(out_figure_dir, "snRNAseq_qc_violin_by_sample.png"),
    width = 7.5,
    height = 3.2
  )

  cluster_plot <- DimPlot(
    obj,
    group.by = "seurat_clusters",
    label = TRUE,
    repel = TRUE,
    reduction = "umap"
  ) +
    ggtitle("Clusters") +
    theme_publication()

  save_publication_plot(
    cluster_plot,
    file.path(out_figure_dir, "snRNAseq_umap_by_cluster.png"),
    width = 4.2,
    height = 3.8
  )

  sample_plot <- DimPlot(
    obj,
    group.by = "sample_id",
    reduction = "umap"
  ) +
    ggtitle("Samples") +
    theme_publication()

  save_publication_plot(
    sample_plot,
    file.path(out_figure_dir, "snRNAseq_umap_by_sample.png"),
    width = 4.2,
    height = 3.8
  )
  score_cols <- grep("^score_.*1$", colnames(obj@meta.data), value = TRUE)

  if (length(score_cols) > 0) {
    score_plots <- lapply(score_cols, function(cn) {
      title <- sub("^score_", "", sub("1$", "", cn))
      FeaturePlot(
        obj,
        features = cn,
        min.cutoff = "q05",
        max.cutoff = "q95",
        reduction = "umap"
      ) +
        ggtitle(title) +
        theme_publication()
    })

    score_plot <- wrap_plots(score_plots, ncol = 2)

    save_publication_plot(
      score_plot,
      file.path(out_figure_dir, "snRNAseq_umap_lineage_scores.png"),
      width = 7.5,
      height = 9.5
    )
  }
}

cat("Reading sample manifest.\n")
samples <- read_sample_manifest(sample_manifest)

cat("Loading samples.\n")
objs <- lapply(seq_len(nrow(samples)), function(i) {
  read_one_sample(
    sample_id = samples$sample_id[[i]],
    gestational_week = samples$gestational_week[[i]],
    matrix_dir = samples$matrix_dir[[i]]
  )
})
names(objs) <- samples$sample_id

cell_counts <- bind_rows(lapply(names(objs), function(nm) {
  count_cells("loaded", objs[[nm]])
}))
cat("Merging samples.\n")
seu <- Reduce(function(x, y) merge(x, y), objs)
seu <- join_layers_if_needed(seu)

cell_counts <- bind_rows(cell_counts, count_cells("merged", seu))

seu <- run_scdblfinder_by_sample(seu, dbr = 0.06)
cell_counts <- bind_rows(cell_counts, count_cells("after_scDblFinder_singlet_filter", seu))

cat("Calculating QC metrics and filtering cells.\n")
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

seu <- subset(
  seu,
  subset =
    nFeature_RNA >= 300 &
    nFeature_RNA <= 20000 &
    nCount_RNA >= 500 &
    percent.mt <= 15
)

cell_counts <- bind_rows(cell_counts, count_cells("after_qc_filter", seu))

cat("Normalizing and clustering.\n")
DefaultAssay(seu) <- "RNA"
seu <- NormalizeData(seu, verbose = FALSE)
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
seu <- ScaleData(seu, features = VariableFeatures(seu), verbose = FALSE)
seu <- RunPCA(seu, features = VariableFeatures(seu), npcs = 50, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30, verbose = FALSE)
seu <- FindNeighbors(seu, dims = 1:30, verbose = FALSE)
seu <- FindClusters(seu, resolution = 0.6, verbose = FALSE)

seu$seurat_clusters <- as.character(Idents(seu))

cat("Adding lineage scores.\n")
seu <- add_lineage_scores(seu)

cat("Writing tables.\n")
write_csv(cell_counts, file.path(out_table_dir, "snRNAseq_cell_counts_by_sample_and_stage.csv"))

qc_summary <- seu@meta.data |>
  group_by(sample_id, gestational_week) |>
  summarise(
    n_cells = n(),
    median_nFeature_RNA = median(nFeature_RNA),
    median_nCount_RNA = median(nCount_RNA),
    median_percent_mt = median(percent.mt),
    .groups = "drop"
  ) |>
  arrange(gestational_week, sample_id)
write_csv(qc_summary, file.path(out_table_dir, "snRNAseq_qc_summary_by_sample.csv"))

cluster_counts <- seu@meta.data |>
  count(seurat_clusters, sample_id, name = "n_cells") |>
  arrange(as.numeric(seurat_clusters), sample_id)

write_csv(cluster_counts, file.path(out_table_dir, "snRNAseq_cluster_counts.csv"))

write_lineage_score_summary(seu)

cat("Saving canonical object.\n")
saveRDS(seu, canonical_rds)

cat("Writing figures.\n")
plot_outputs(seu)


cat("Writing logs.\n")
summary_lines <- c(
  paste("Canonical object:", canonical_rds),
  paste("Cells:", ncol(seu)),
  paste("Genes:", nrow(seu)),
  paste("Samples:", paste(sort(unique(seu$sample_id)), collapse = ", ")),
  paste("Gestational weeks:", paste(sort(unique(seu$gestational_week)), collapse = ", ")),
  paste("Clusters:", paste(sort(unique(seu$seurat_clusters)), collapse = ", ")),
  paste("Output directory:", results_root)
)

writeLines(summary_lines, file.path(out_log_dir, "snRNAseq_initial_processing_summary.txt"))

sink(file.path(out_log_dir, "snRNAseq_initial_processing_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Done.\n")
cat("Canonical object: ", canonical_rds, "\n", sep = "")
cat("Output directory: ", results_root, "\n", sep = "")
