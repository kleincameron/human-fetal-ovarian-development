suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(scales)
  library(grid)
  library(ggtext)
})

set.seed(42)

project_root <- "/home/liyan/liyan/Final/github_code_for_publication"
input_root <- "/home/liyan/liyan/Final/github_code_for_publication_results/snRNAseq_annotated_object"
results_root <- "/home/liyan/liyan/Final/github_code_for_publication_results/snRNAseq_marker_dotplots"

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
markers_by_group <- list(
  "Germ: mitotic oogonia" = c("POU5F1", "NANOG", "LIN28A", "DDX4", "DPPA3", "UTF1"),
  "Germ: pre-meiotic" = c("STRA8", "IRF1", "ZGLP1"),
  "Germ: meiotic prophase" = c("SYCP3", "SYCP1", "SMC1B", "REC8", "IL13RA2", "MEIOC"),
  "Germ: oocytes and primordial follicles" = c("DDX4", "FIGLA", "ZP3", "ZP1", "GDF9", "NOBOX", "LHX8", "PECAM1", "NPM2"),

  "Granulosa: general" = c("FOXL2", "WNT4", "RSPO1", "GATA4", "WT1", "WNT6"),
  "Granulosa: first-wave" = c("FOXL2", "RSPO1"),
  "Granulosa: second-wave" = c("FOXL2", "AMHR2", "EMX2", "BMP2"),
  "Granulosa: primordial follicle" = c("FOXL2", "NOTCH3", "HEYL", "NR1H4", "PBX3"),

  "Stroma" = c("TCF21", "PDGFRA", "NR2F2", "NR2F1", "POSTN"),
  "Endothelial" = c("PECAM1", "VWF", "CDH5", "ESM1", "ROBO4"),
  "Mural" = c("ACTA2", "PDGFRB", "RGS5", "MYH11", "TAGLN"),

  "Immune: macrophage" = c("CD68", "CD14", "CSF1R", "HLA-DRA", "CD163", "SIGLEC1"),
  "Immune: lymphoid" = c("CD3D", "KLRB1", "NKG7", "TRDC"),

  "Erythroid" = c("HBA1", "HBA2", "HBB", "HBG1", "HBG2", "ALAS2", "AHSP", "GYPA")
)

genes_to_remove <- unique(c(
  "AMH", "PAX8", "OSR1", "LGR5", "AXIN2", "GJA1", "FST",
  "GATA2", "COL3A1", "DES", "NCAM1", "DCN"
))

deduplicate_markers <- function(markers_list) {
  seen <- character(0)
  retained <- list()
  duplicated <- list()

  for (group_name in names(markers_list)) {
    genes <- unique(markers_list[[group_name]])
    genes <- setdiff(genes, genes_to_remove)

    duplicate_genes <- genes[genes %in% seen]
    retained_genes <- genes[!genes %in% seen]

    if (length(duplicate_genes) > 0) {
      duplicated[[group_name]] <- duplicate_genes
    }

    if (length(retained_genes) > 0) {
      retained[[group_name]] <- retained_genes
    }

    seen <- c(seen, retained_genes)
  }

  list(retained = retained, duplicated = duplicated)
}
safe_zscore <- function(x) {
  if (all(is.na(x))) {
    return(rep(NA_real_, length(x)))
  }

  s <- stats::sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)

  if (is.na(s) || s == 0) {
    return(rep(0, length(x)))
  }

  (x - m) / s
}

html_escape <- function(x) {
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x
}

make_colored_axis_labels <- function(levels, colors) {
  colors <- colors[levels]
  colors[is.na(colors)] <- "black"

  labels <- sprintf(
    "<span style='color:%s'>%s</span>",
    colors,
    html_escape(levels)
  )
  names(labels) <- levels
  labels
}

order_final_annotations <- function(metadata) {
  metadata |>
    distinct(major_cell_type, fine_subcluster, final_annotation) |>
    mutate(
      major_cell_type = factor(as.character(major_cell_type), levels = celltype_order),
      subcluster_number = suppressWarnings(as.integer(sub("^.*_", "", fine_subcluster)))
    ) |>
    arrange(major_cell_type, subcluster_number, fine_subcluster) |>
    pull(final_annotation) |>
    unique()
}

make_marker_dotplot <- function(
  seu,
  group_col,
  group_levels,
  group_label,
  output_prefix,
  marker_map,
  y_label_colors
) {
  stopifnot(group_col %in% colnames(seu@meta.data))

  metadata <- seu@meta.data |>
    mutate(
      plot_group = as.character(.data[[group_col]])
    ) |>
    filter(!is.na(plot_group), plot_group != "")

  group_levels <- intersect(group_levels, unique(metadata$plot_group))
  if (length(group_levels) == 0) {
    stop("No valid group levels for ", group_col)
  }

  genes_use <- intersect(as.character(marker_map$gene), rownames(seu))
  if (length(genes_use) == 0) {
    stop("No marker genes are present in the object.")
  }
  expr <- tryCatch(
    GetAssayData(seu, assay = "RNA", layer = "data"),
    error = function(e) GetAssayData(seu, assay = "RNA", slot = "data")
  )

  expr <- expr[genes_use, rownames(metadata), drop = FALSE]
  plot_group <- factor(metadata$plot_group, levels = group_levels)

  avg_mat <- sapply(group_levels, function(group_name) {
    cells_group <- which(plot_group == group_name)
    if (length(cells_group) == 0) {
      return(rep(NA_real_, length(genes_use)))
    }
    Matrix::rowMeans(expr[, cells_group, drop = FALSE])
  })

  pct_mat <- sapply(group_levels, function(group_name) {
    cells_group <- which(plot_group == group_name)
    if (length(cells_group) == 0) {
      return(rep(NA_real_, length(genes_use)))
    }
    Matrix::rowMeans(expr[, cells_group, drop = FALSE] > 0) * 100
  })

  rownames(avg_mat) <- genes_use
  rownames(pct_mat) <- genes_use
  colnames(avg_mat) <- group_levels
  colnames(pct_mat) <- group_levels

  dotplot_df <- as.data.frame(as.table(avg_mat), stringsAsFactors = FALSE) |>
    rename(gene = Var1, group_level = Var2, average_expression = Freq) |>
    inner_join(
      as.data.frame(as.table(pct_mat), stringsAsFactors = FALSE) |>
        rename(gene = Var1, group_level = Var2, percent_expressing = Freq),
      by = c("gene", "group_level")
    ) |>
    left_join(as.data.frame(marker_map), by = "gene") |>
    group_by(gene) |>
    mutate(scaled_average_expression = safe_zscore(average_expression)) |>
    ungroup() |>
    mutate(
      scaled_average_expression = pmax(pmin(scaled_average_expression, 2.5), -2.5),
      gene = factor(gene, levels = levels(marker_map$gene)),
      group_level = factor(group_level, levels = rev(group_levels))
    )

  group_boundaries <- marker_map |>
    as.data.frame() |>
    mutate(gene_index = as.integer(gene)) |>
    group_by(gene_group) |>
    summarise(
      start = min(gene_index),
      end = max(gene_index),
      midpoint = (start + end) / 2,
      .groups = "drop"
    )

  write_csv(
    dotplot_df,
    file.path(out_table_dir, paste0(output_prefix, "_dotplot_values.csv"))
  )
  plot_width <- max(7.2, 0.16 * length(genes_use) + 2.0)
  plot_height <- max(2.8, 0.22 * length(group_levels) + 1.6)

  y_axis_labels <- make_colored_axis_labels(group_levels, y_label_colors)

  p <- ggplot(dotplot_df, aes(x = gene, y = group_level)) +
    geom_point(aes(size = percent_expressing, color = scaled_average_expression)) +
    scale_size(
      range = c(0.2, 2.8),
      limits = c(0, 100),
      breaks = c(25, 50, 75, 100),
      name = "% expressing"
    ) +
    scale_color_gradient2(
      low = "#2166AC",
      mid = "#F7F7F7",
      high = "#B2182B",
      midpoint = 0,
      limits = c(-2.5, 2.5),
      oob = scales::squish,
      name = "Mean expression\n(z score)"
    ) +
    scale_y_discrete(labels = y_axis_labels) +
    labs(x = NULL, y = group_label) +
    theme_publication() +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1,
        vjust = 1,
        family = publication_font_family,
        size = publication_base_size
      ),
      axis.text.y = ggtext::element_markdown(
        family = publication_font_family,
        size = publication_base_size
      ),
      legend.title = element_text(
        family = publication_font_family,
        size = publication_base_size
      ),
      legend.text = element_text(
        family = publication_font_family,
        size = publication_base_size
      ),
      legend.key.height = unit(0.32, "cm"),
      legend.key.width = unit(0.32, "cm"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )

  if (nrow(group_boundaries) > 1) {
    p <- p +
      geom_vline(
        xintercept = group_boundaries$end[-nrow(group_boundaries)] + 0.5,
        linewidth = 0.2,
        alpha = 0.5
      )
  }
  save_publication_plot(
    p,
    file.path(out_figure_dir, paste0(output_prefix, ".png")),
    width = plot_width,
    height = plot_height,
    dpi = 600
  )

  invisible(p)
}

cat("Loading annotated snRNA-seq object.\n")
stopifnot(file.exists(annotated_rds))
seu <- readRDS(annotated_rds)
stopifnot(inherits(seu, "Seurat"))

DefaultAssay(seu) <- "RNA"

required_cols <- c("major_cell_type", "fine_subcluster", "final_annotation")
missing_cols <- setdiff(required_cols, colnames(seu@meta.data))
if (length(missing_cols) > 0) {
  stop("Missing required metadata columns: ", paste(missing_cols, collapse = ", "))
}

extra_cell_types <- setdiff(unique(as.character(seu$major_cell_type)), celltype_order)
if (length(extra_cell_types) > 0) {
  stop("Unexpected major cell types: ", paste(extra_cell_types, collapse = ", "))
}

deduplicated <- deduplicate_markers(markers_by_group)

markers_present <- lapply(
  deduplicated$retained,
  function(genes) intersect(genes, rownames(seu))
)

markers_present <- markers_present[lengths(markers_present) > 0]

if (length(markers_present) == 0) {
  stop("No marker genes are present in the object.")
}

marker_map <- bind_rows(lapply(names(markers_present), function(group_name) {
  data.frame(
    gene_group = group_name,
    gene = markers_present[[group_name]],
    stringsAsFactors = FALSE
  )
}))
marker_map$gene_group <- factor(marker_map$gene_group, levels = unique(marker_map$gene_group))
marker_map$gene <- factor(marker_map$gene, levels = unique(marker_map$gene))

missing_genes <- setdiff(
  unique(unlist(deduplicated$retained, use.names = FALSE)),
  rownames(seu)
)

removed_duplicate_table <- bind_rows(lapply(names(deduplicated$duplicated), function(group_name) {
  tibble(
    gene_group = group_name,
    removed_duplicate_gene = deduplicated$duplicated[[group_name]]
  )
}))

if (nrow(removed_duplicate_table) == 0) {
  removed_duplicate_table <- tibble(
    gene_group = character(),
    removed_duplicate_gene = character()
  )
}

write_csv(marker_map, file.path(out_table_dir, "snRNAseq_marker_dotplot_gene_map.csv"))

write_csv(
  tibble(missing_gene = missing_genes),
  file.path(out_table_dir, "snRNAseq_marker_dotplot_missing_genes.csv")
)

write_csv(
  removed_duplicate_table,
  file.path(out_table_dir, "snRNAseq_marker_dotplot_removed_duplicate_genes.csv")
)

major_cell_type_levels <- intersect(celltype_order, unique(as.character(seu$major_cell_type)))

final_annotation_levels <- order_final_annotations(seu@meta.data)

major_cell_type_y_colors <- celltype_colors[major_cell_type_levels]

final_annotation_to_celltype <- seu@meta.data |>
  distinct(final_annotation, major_cell_type) |>
  mutate(
    final_annotation = as.character(final_annotation),
    major_cell_type = as.character(major_cell_type)
  )

if (anyDuplicated(final_annotation_to_celltype$final_annotation) > 0) {
  stop("One or more final annotations map to multiple major cell types.")
}

final_annotation_y_colors <- celltype_colors[final_annotation_to_celltype$major_cell_type]
names(final_annotation_y_colors) <- final_annotation_to_celltype$final_annotation
final_annotation_y_colors <- final_annotation_y_colors[final_annotation_levels]

make_marker_dotplot(
  seu = seu,
  group_col = "major_cell_type",
  group_levels = major_cell_type_levels,
  group_label = "Major cell type",
  output_prefix = "snRNAseq_marker_dotplot_major_cell_type",
  marker_map = marker_map,
  y_label_colors = major_cell_type_y_colors
)

make_marker_dotplot(
  seu = seu,
  group_col = "final_annotation",
  group_levels = final_annotation_levels,
  group_label = "Final annotation",
  output_prefix = "snRNAseq_marker_dotplot_final_annotation",
  marker_map = marker_map,
  y_label_colors = final_annotation_y_colors
)
summary_lines <- c(
  paste("Input annotated object:", annotated_rds),
  paste("Cells:", ncol(seu)),
  paste("Genes:", nrow(seu)),
  paste("Marker groups:", length(markers_present)),
  paste("Marker genes plotted:", length(unique(as.character(marker_map$gene)))),
  paste("Missing requested genes:", length(missing_genes)),
  paste("Duplicate requested genes removed:", nrow(removed_duplicate_table)),
  paste("Major cell type groups plotted:", length(major_cell_type_levels)),
  paste("Final annotation groups plotted:", length(final_annotation_levels)),
  paste("Output directory:", results_root)
)

writeLines(summary_lines, file.path(out_log_dir, "snRNAseq_marker_dotplots_summary.txt"))

sink(file.path(out_log_dir, "snRNAseq_marker_dotplots_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Done.\n")
cat("Marker dotplots written to: ", out_figure_dir, "\n", sep = "")
