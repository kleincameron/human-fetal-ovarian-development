suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(readr)
  library(tibble)
  library(pheatmap)
  library(grid)
})

options(bitmapType = "cairo")

# ==========================================================
# PATHS
# ==========================================================
project_root <- "/home/liyan/liyan/Final/github_code_for_publication"
results_root <- "/home/liyan/liyan/Final/github_code_for_publication_results/snRNAseq_selective_cell_death_pathways"

input_rds <- file.path(
  "/home/liyan/liyan/Final/github_code_for_publication_results",
  "snRNAseq_annotated_object",
  "objects",
  "fetal_ovary_snRNAseq_canonical_annotated.rds"
)

out_figure_dir <- file.path(results_root, "figures")
out_table_dir  <- file.path(results_root, "tables")
out_log_dir    <- file.path(results_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir,    recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_rds))

# ==========================================================
# FIGURE STYLE
# ==========================================================
publication_font_family <- "Helvetica"
publication_base_size   <- 8

heatmap_width_in  <- 12.0
heatmap_height_in <- 6.8
png_dpi           <- 600

heatmap_title <- "Selective cell death pathways — snRNA-seq fine subclusters"

# ==========================================================
# ESTABLISHED CELL-TYPE COLORS
# ==========================================================
major_lineage_colors <- c(
  germ        = "#1B9E77",
  degenerated = "#B0B0B0",
  granulosa   = "#D95F02",
  stroma      = "#7570B3",
  endothelial = "#E7298A",
  mural       = "#66A61E",
  immune      = "#E6AB02",
  erythroid   = "#A6761D"
)

# ==========================================================
# CURATED PATHWAY GENE SETS
# Adjust here later if you want to refine the pathway definitions.
# ==========================================================
pathway_gene_sets <- list(
  Lysosomal = c(
    "LAMP1", "LAMP2", "CTSB", "CTSD", "CTSL", "GAA", "ATP6V1A", "ATP6V0D1",
    "MCOLN1", "NPC1"
  ),
  Pyroptosis = c(
    "CASP1", "CASP4", "CASP5", "GSDMD", "NLRP3", "PYCARD", "IL1B", "IL18"
  ),
  Apoptosis = c(
    "BAX", "BAK1", "BCL2L11", "CASP3", "CASP7", "CASP8", "CASP9", "APAF1",
    "FAS", "PMAIP1"
  ),
  Necroptosis = c(
    "RIPK1", "RIPK3", "MLKL", "FADD", "CASP8", "TNFRSF1A", "TNF", "TLR3",
    "ZBP1", "CYLD"
  ),
  Ferroptosis = c(
    "GPX4", "SLC7A11", "ACSL4", "ALOX15", "TFRC", "NCOA4", "LPCAT3",
    "FTH1", "FTL", "SAT1"
  ),
  Autophagy = c(
    "BECN1", "ATG5", "ATG7", "ATG12", "MAP1LC3B", "SQSTM1", "ULK1",
    "WIPI1", "GABARAPL1", "ATG3"
  ),
  Cuproptosis = c(
    "FDX1", "LIAS", "LIPT1", "DLD", "DLAT", "PDHA1", "PDHB", "MTF1",
    "GLS", "CDKN2A"
  ),
  Parthanatos = c(
    "PARP1", "PARP2", "AIFM1", "XRCC1", "POLB", "MIF", "HMGB1", "SIRT6"
  )
)

pathway_order <- names(pathway_gene_sets)

# ==========================================================
# HELPERS
# ==========================================================
mode_string <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) {
    return(NA_character_)
  }
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

shorten_subcluster_label <- function(x) {
  dplyr::recode(
    x,
    "Atresia-stressed degenerating follicle cells" = "Atresia-stressed",
    "Clearance-associated degenerating follicle cells" = "Clearance-associated",
    .default = x
  )
}
infer_major_lineage <- function(fine_subcluster, final_annotation) {
  lbl <- tolower(ifelse(is.na(final_annotation), "", final_annotation))
  fs  <- tolower(ifelse(is.na(fine_subcluster), "", fine_subcluster))

  if (grepl("atresia|clearance-associated degenerating follicle|degenerating follicle", lbl)) {
    return("degenerated")
  }
  if (grepl("germ|oogonia|oocyte|meiotic|pachytene|leptotene|zygote", lbl)) {
    return("germ")
  }
  if (grepl("granulosa", lbl)) {
    return("granulosa")
  }
  if (grepl("stroma", lbl)) {
    return("stroma")
  }
  if (grepl("endothelial", lbl)) {
    return("endothelial")
  }
  if (grepl("mural|pericyte|vsmc|smooth muscle", lbl)) {
    return("mural")
  }
  if (grepl("immune|macrophage|nk t|tissue macrophages|nk", lbl)) {
    return("immune")
  }
  if (grepl("erythroid", lbl)) {
    return("erythroid")
  }

  prefix <- sub("_.*$", "", fs)

  dplyr::case_when(
    prefix == "germ"        ~ "germ",
    prefix == "granulosa"   ~ "granulosa",
    prefix == "stroma"      ~ "stroma",
    prefix == "endothelial" ~ "endothelial",
    prefix == "mural"       ~ "mural",
    prefix == "immune"      ~ "immune",
    prefix == "erythroid"   ~ "erythroid",
    TRUE                    ~ "germ"
  )
}

save_pheatmap_dual <- function(ph, pdf_file, png_file, width, height, dpi = 600) {
  pdf(
    pdf_file,
    width = width,
    height = height,
    useDingbats = FALSE,
    family = publication_font_family
  )
  grid.newpage()
  grid.draw(ph$gtable)
  dev.off()

  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(
      filename = png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      background = "white"
    )
  } else {
    png(
      filename = png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white"
    )
  }

  grid.newpage()
  grid.draw(ph$gtable)
  dev.off()
}

# ==========================================================
# LOAD OBJECT
# ==========================================================
cat("Loading annotated snRNA-seq object:\n ", input_rds, "\n", sep = "")
obj <- readRDS(input_rds)

stopifnot(inherits(obj, "Seurat"))
stopifnot("RNA" %in% names(obj@assays))
stopifnot("fine_subcluster" %in% colnames(obj@meta.data))

DefaultAssay(obj) <- "RNA"

# Join layers if needed, then ensure a normalized data layer exists
if ("JoinLayers" %in% getNamespaceExports("SeuratObject")) {
  obj <- tryCatch(
    JoinLayers(obj, assay = "RNA"),
    error = function(e) obj
  )
}

rna_layers <- Layers(obj[["RNA"]])

if (!("data" %in% rna_layers)) {
  cat("No RNA data layer found. Running NormalizeData().\n")
  obj <- NormalizeData(obj, verbose = FALSE)
  rna_layers <- Layers(obj[["RNA"]])
}
stopifnot("data" %in% rna_layers)

meta <- obj@meta.data |>
  tibble::rownames_to_column("cell_id")

if ("final_annotation" %in% colnames(meta)) {
  label_column <- "final_annotation"
} else {
  label_column <- "fine_subcluster"
}

# ==========================================================
# DEFINE GROUP LABELS
# ==========================================================
group_map <- meta |>
  group_by(fine_subcluster) |>
  summarise(
    final_annotation = mode_string(.data[[label_column]]),
    n_cells = dplyr::n(),
    .groups = "drop"
  ) |>
  mutate(
    display_label = shorten_subcluster_label(final_annotation),
    major_lineage = vapply(
      seq_len(n()),
      function(i) infer_major_lineage(
        fine_subcluster = fine_subcluster[[i]],
        final_annotation = final_annotation[[i]]
      ),
      character(1)
    )
  )

stopifnot(nrow(group_map) > 0)

# ==========================================================
# AVERAGE NORMALIZED EXPRESSION BY FINE SUBCLUSTER
# ==========================================================
cat("Calculating average normalized expression by fine_subcluster.\n")

# Do not use AverageExpression() here. In Seurat v5, it rewrites identity
# names containing underscores into dashes, which can break downstream
# matching to fine_subcluster metadata. Instead, compute group means directly
# from the normalized RNA data layer while preserving original labels.

data_mat <- tryCatch(
  GetAssayData(obj, assay = "RNA", layer = "data"),
  error = function(e) GetAssayData(obj, assay = "RNA", slot = "data")
)

stopifnot(is.matrix(data_mat) || inherits(data_mat, "Matrix"))

valid_cells <- intersect(meta$cell_id, colnames(data_mat))
if (length(valid_cells) == 0) {
  stop("No metadata cell IDs overlap the RNA data matrix column names.")
}

meta <- meta |>
  filter(cell_id %in% valid_cells)

group_map <- group_map |>
  filter(fine_subcluster %in% unique(meta$fine_subcluster))

if (nrow(group_map) == 0) {
  stop("No fine_subcluster groups remain after matching metadata to RNA data matrix.")
}

fine_subcluster_levels <- group_map$fine_subcluster

avg_expr_list <- lapply(fine_subcluster_levels, function(group_name) {
  cells_group <- meta$cell_id[meta$fine_subcluster == group_name]
  cells_group <- intersect(cells_group, colnames(data_mat))

  if (length(cells_group) == 0) {
    stop("No cells found for fine_subcluster: ", group_name)
  }

  Matrix::rowMeans(data_mat[, cells_group, drop = FALSE])
})

avg_expr <- do.call(cbind, avg_expr_list)
rownames(avg_expr) <- rownames(data_mat)
colnames(avg_expr) <- fine_subcluster_levels

stopifnot(identical(colnames(avg_expr), group_map$fine_subcluster))

# ==========================================================
# PATHWAY SCORING
# Mean normalized expression of pathway genes, then row z-score
# ==========================================================
score_matrix <- matrix(
  NA_real_,
  nrow = length(pathway_gene_sets),
  ncol = ncol(avg_expr),
  dimnames = list(names(pathway_gene_sets), colnames(avg_expr))
)

gene_set_table <- lapply(names(pathway_gene_sets), function(pathway_name) {
  genes_requested <- unique(pathway_gene_sets[[pathway_name]])
  genes_present   <- intersect(genes_requested, rownames(avg_expr))

  if (length(genes_present) == 0) {
    stop("No genes from pathway '", pathway_name, "' are present in the RNA assay.")
  }

  score_matrix[pathway_name, ] <<- colMeans(avg_expr[genes_present, , drop = FALSE])

  tibble(
    pathway = pathway_name,
    gene_requested = genes_requested,
    gene_present = genes_requested %in% genes_present
  )
})

gene_set_table <- bind_rows(gene_set_table)
write_csv(gene_set_table, file.path(out_table_dir, "selective_cell_death_pathway_gene_sets.csv"))

score_matrix <- score_matrix[pathway_order, , drop = FALSE]

scale_row_z <- function(x) {
  s <- stats::sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

score_matrix_z <- t(apply(score_matrix, 1, scale_row_z))
score_matrix_z <- score_matrix_z[pathway_order, , drop = FALSE]

# Clip to match the visual dynamic range of the original figure
clip_low  <- -2
clip_high <- 4
score_matrix_z_clipped <- pmin(pmax(score_matrix_z, clip_low), clip_high)
# ==========================================================
# SAVE TABLES
# ==========================================================
score_long <- as.data.frame(score_matrix) |>
  rownames_to_column("pathway") |>
  tidyr::pivot_longer(
    cols = -pathway,
    names_to = "fine_subcluster",
    values_to = "mean_normalized_expression"
  ) |>
  left_join(group_map, by = "fine_subcluster")

score_z_long <- as.data.frame(score_matrix_z) |>
  rownames_to_column("pathway") |>
  tidyr::pivot_longer(
    cols = -pathway,
    names_to = "fine_subcluster",
    values_to = "row_zscore"
  ) |>
  left_join(group_map, by = "fine_subcluster")

write_csv(score_long, file.path(out_table_dir, "selective_cell_death_pathway_scores_mean_expression.csv"))
write_csv(score_z_long, file.path(out_table_dir, "selective_cell_death_pathway_scores_row_zscore.csv"))
write_csv(group_map, file.path(out_table_dir, "selective_cell_death_pathway_subcluster_labels.csv"))

# ==========================================================
# HEATMAP ANNOTATION
# ==========================================================
annotation_col <- data.frame(
  major_lineage = factor(
    group_map$major_lineage,
    levels = c("germ", "degenerated", "granulosa", "stroma", "endothelial", "mural", "immune", "erythroid")
  ),
  row.names = group_map$fine_subcluster,
  stringsAsFactors = FALSE
)

annotation_colors <- list(
  major_lineage = major_lineage_colors
)

display_labels <- group_map$display_label
names(display_labels) <- group_map$fine_subcluster

heat_colors <- colorRampPalette(c("#00008B", "#F7F7F7", "#D73027"))(100)
heat_breaks <- seq(clip_low, clip_high, length.out = length(heat_colors) + 1)

# ==========================================================
# DRAW HEATMAP
# ==========================================================
cat("Generating heatmap.\n")

ph <- pheatmap(
  mat = score_matrix_z_clipped,
  color = heat_colors,
  breaks = heat_breaks,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_col = annotation_col,
  annotation_colors = annotation_colors,
  labels_col = display_labels[colnames(score_matrix_z_clipped)],
  angle_col = 90,
  border_color = NA,
  fontsize = publication_base_size,
  fontsize_row = 8,
  fontsize_col = 6.5,
  annotation_names_col = TRUE,
  annotation_legend = TRUE,
  legend = TRUE,
  treeheight_row = 50,
  treeheight_col = 35,
  main = heatmap_title,
  silent = TRUE
)

figure_pdf <- file.path(
  out_figure_dir,
  "snRNAseq_selective_cell_death_pathways_fine_subcluster_heatmap.pdf"
)
figure_png <- file.path(
  out_figure_dir,
  "snRNAseq_selective_cell_death_pathways_fine_subcluster_heatmap.png"
)

save_pheatmap_dual(
  ph = ph,
  pdf_file = figure_pdf,
  png_file = figure_png,
  width = heatmap_width_in,
  height = heatmap_height_in,
  dpi = png_dpi
)

# ==========================================================
# SAVE ORDERED MATRICES (helpful for review/manuscript)
# ==========================================================
row_order <- rownames(score_matrix_z_clipped)[ph$tree_row$order]
col_order <- colnames(score_matrix_z_clipped)[ph$tree_col$order]

ordered_z <- score_z_long |>
  mutate(
    pathway = factor(pathway, levels = row_order),
    fine_subcluster = factor(fine_subcluster, levels = col_order)
  ) |>
  arrange(pathway, fine_subcluster)

write_csv(
  ordered_z,
  file.path(out_table_dir, "selective_cell_death_pathway_scores_row_zscore_clustered_order.csv")
)
# ==========================================================
# SUMMARY + SESSION INFO
# ==========================================================
summary_lines <- c(
  paste("Input object:", input_rds),
  paste("Label column used:", label_column),
  paste("Number of fine subclusters plotted:", nrow(group_map)),
  paste("Number of pathways plotted:", length(pathway_gene_sets)),
  paste("Figure PDF:", figure_pdf),
  paste("Figure PNG:", figure_png),
  paste("Pathway-score method: mean normalized RNA expression across genes in each curated pathway, then row z-score across fine subclusters."),
  paste("Z-score clipping range:", paste0("[", clip_low, ", ", clip_high, "]"))
)

writeLines(
  summary_lines,
  file.path(out_log_dir, "selective_cell_death_pathway_heatmap_summary.txt")
)

sink(file.path(out_log_dir, "selective_cell_death_pathway_heatmap_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone.\n")
cat("Figure written to:\n  ", figure_pdf, "\n  ", figure_png, "\n", sep = "")
cat("Results root:\n  ", results_root, "\n", sep = "")
