suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(stringr)
  library(tibble)
})

options(bitmapType = "cairo")

# ==========================================================
# PATHS
# ==========================================================
input_rds <- "/home/liyan/liyan/Final/github_code_for_publication_results/snRNAseq_annotated_object/objects/fetal_ovary_snRNAseq_canonical_annotated.rds"
results_root <- "/home/liyan/liyan/Final/github_code_for_publication_results/snRNAseq_germ_stage_validation_panel"

out_figure_dir <- file.path(results_root, "figures")
out_table_dir  <- file.path(results_root, "tables")
out_log_dir    <- file.path(results_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir,    recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_rds))

# ==========================================================
# SETTINGS
# ==========================================================
# Top-to-bottom order to match the earlier validation panel.
subcluster_order_top_to_bottom <- c(
  "germ_7",
  "germ_5",
  "germ_4",
  "germ_0",
  "germ_1",
  "germ_6",
  "germ_3",
  "germ_8",
  "germ_2"
)

gene_modules <- list(
  "Mitotic / PGC-like" = c(
    "POU5F1", "TFAP2C", "SOX17", "NANOG", "DPPA3", "PRDM1",
    "KIT", "UTF1", "IFITM3", "MKI67", "TOP2A"
  ),
  "Meiotic entry / prophase" = c(
    "STRA8", "REC8", "SYCP3", "SMC1B", "DMC1", "MEIOC",
    "HORMAD1", "HORMAD2", "TEX12", "SPO11", "SYCE2"
  ),
  "Late meiotic / oocyte" = c(
    "FIGLA", "LHX8", "NOBOX", "GDF9", "ZP3", "NPM2",
    "OOEP", "PATL2", "TUBB8", "KHDC3L", "ZAR1"
  ),
  "Stress / degeneration" = c(
    "FOS", "JUN", "HSPA1A", "HSPA1B", "DNAJB1", "ATF3"
  )
)

figure_width  <- 12.0
figure_height <- 4.8
figure_dpi    <- 600
base_family   <- "Helvetica"
base_size     <- 8

# ==========================================================
# HELPERS
# ==========================================================
save_pdf_png <- function(plot_obj, pdf_file, png_file,
                         width, height, dpi = 600, bg = "white") {
  grDevices::pdf(
    file = pdf_file,
    width = width,
    height = height,
    useDingbats = FALSE,
    family = base_family
  )
  print(plot_obj)
  dev.off()

  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(
      filename = png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      background = bg
    )
    print(plot_obj)
    dev.off()
  } else {
    grDevices::png(
      filename = png_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = bg
    )
    print(plot_obj)
    dev.off()
  }
}

mode_label <- function(x) {
  x <- x[!is.na(x) & x != ""]
  if (length(x) == 0) return(NA_character_)
  tb <- sort(table(x), decreasing = TRUE)
  names(tb)[1]
}

shorten_label <- function(x) {
  x <- as.character(x)

  x <- dplyr::recode(
    x,
    "Atresia-stressed degenerating follicle cells" = "Atresia-stressed",
    "Clearance-associated degenerating follicle cells" = "Clearance-associated",
    .default = x
  )

  x
}
# fallback only if final_annotation is not available
fallback_label_map <- c(
  germ_0 = "Atresia-stressed",
  germ_1 = "Pachytene/diplotene germ cells",
  germ_2 = "Meiotic-entry germ cells",
  germ_3 = "Mitotic oogonia",
  germ_4 = "Stalled meiotic germ cells",
  germ_5 = "Leptotene/zygotene germ cells",
  germ_6 = "Clearance-associated",
  germ_7 = "Primordial follicle oocytes",
  germ_8 = "Degenerating germ cells"
)

# ==========================================================
# LOAD OBJECT
# ==========================================================
cat("Loading annotated snRNA-seq object:\n ", input_rds, "\n", sep = "")
obj <- readRDS(input_rds)

stopifnot("fine_subcluster" %in% colnames(obj@meta.data))

meta <- obj@meta.data |>
  rownames_to_column("cell_id")

germ_cells <- meta |>
  filter(grepl("^germ_", fine_subcluster)) |>
  pull(cell_id)

if (length(germ_cells) == 0) {
  stop("No germ_* cells found in fine_subcluster.")
}

obj_germ <- subset(obj, cells = germ_cells)
meta_germ <- obj_germ@meta.data |>
  rownames_to_column("cell_id")

# ==========================================================
# SUBCLUSTER LABEL MAP
# ==========================================================
present_subclusters <- unique(as.character(meta_germ$fine_subcluster))
present_order_top_to_bottom <- subcluster_order_top_to_bottom[
  subcluster_order_top_to_bottom %in% present_subclusters
]

if (length(present_order_top_to_bottom) == 0) {
  stop("None of the expected germ subclusters were found.")
}

if ("final_annotation" %in% colnames(meta_germ)) {
  label_map <- meta_germ |>
    group_by(fine_subcluster) |>
    summarise(
      final_annotation = mode_label(final_annotation),
      .groups = "drop"
    ) |>
    mutate(display_label = shorten_label(final_annotation))
} else {
  label_map <- tibble(
    fine_subcluster = present_order_top_to_bottom,
    final_annotation = fallback_label_map[present_order_top_to_bottom],
    display_label = shorten_label(final_annotation)
  )
}

# if labels are duplicated, append raw subcluster for uniqueness
dup_labels <- label_map$display_label[duplicated(label_map$display_label)]
if (length(dup_labels) > 0) {
  label_map <- label_map |>
    mutate(
      display_label = ifelse(
        display_label %in% dup_labels,
        paste0(display_label, " (", fine_subcluster, ")"),
        display_label
      )
    )
}

label_map <- label_map |>
  filter(fine_subcluster %in% present_order_top_to_bottom)

write_csv(
  label_map,
  file.path(out_table_dir, "germ_subcluster_label_map.csv")
)

# ==========================================================
# GENES
# ==========================================================
feature_order <- unlist(gene_modules, use.names = FALSE)
feature_order <- unique(feature_order)

genes_present <- feature_order[feature_order %in% rownames(obj_germ)]
genes_missing <- setdiff(feature_order, genes_present)

if (length(genes_present) == 0) {
  stop("None of the requested marker genes are present in the object.")
}

write_lines(
  genes_missing,
  file.path(out_table_dir, "germ_stage_validation_missing_genes.txt")
)
gene_module_map <- tibble(
  gene = unlist(gene_modules, use.names = FALSE),
  module = rep(names(gene_modules), lengths(gene_modules))
) |>
  distinct(gene, .keep_all = TRUE) |>
  filter(gene %in% genes_present)

# ==========================================================
# DOTPLOT INPUT
# ==========================================================
# Use a synthetic plotting group to avoid any underscore/dash identity rewriting.
plot_group_ids <- paste0("group_", seq_along(present_order_top_to_bottom))
names(plot_group_ids) <- present_order_top_to_bottom

obj_germ$plot_group <- unname(plot_group_ids[as.character(obj_germ$fine_subcluster)])
obj_germ$plot_group <- factor(
  obj_germ$plot_group,
  levels = plot_group_ids[present_order_top_to_bottom]
)

group_key <- tibble(
  plot_group = plot_group_ids[present_order_top_to_bottom],
  fine_subcluster = present_order_top_to_bottom
) |>
  left_join(label_map, by = "fine_subcluster")

dot_data <- DotPlot(
  object = obj_germ,
  assay = "RNA",
  features = genes_present,
  group.by = "plot_group",
  scale = TRUE,
  dot.scale = 6
)$data |>
  as_tibble() |>
  left_join(group_key, by = c("id" = "plot_group")) |>
  left_join(gene_module_map, by = c("features.plot" = "gene")) |>
  mutate(
    avg.exp.scaled = pmax(pmin(avg.exp.scaled, 2), -2),
    features.plot = factor(features.plot, levels = genes_present),
    module = factor(module, levels = names(gene_modules)),
    display_label = factor(
      display_label,
      levels = rev(group_key$display_label)  # bottom-to-top factor so top matches requested order
    )
  )

write_csv(
  dot_data,
  file.path(out_table_dir, "germ_stage_validation_dotplot_data.csv")
)

# ==========================================================
# PLOT
# ==========================================================
figure_title <- "Germ subcluster stage-validation panel"

p <- ggplot(
  dot_data,
  aes(x = features.plot, y = display_label)
) +
  geom_point(
    aes(size = pct.exp, color = avg.exp.scaled)
  ) +
  facet_grid(
    . ~ module,
    scales = "free_x",
    space = "free_x"
  ) +
  scale_color_gradient2(
    low = "#D9D9D9",
    mid = "#8E75E6",
    high = "#1F3CFF",
    midpoint = 0,
    limits = c(-2, 2),
    breaks = c(-2, -1, 0, 1, 2),
    name = "Average expression\n(scaled)"
  ) +
  scale_size(
    range = c(0.1, 6.0),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75),
    name = "Percent expressed"
  ) +
  labs(
    title = figure_title,
    x = NULL,
    y = NULL
  ) +
  theme_classic(base_size = base_size, base_family = base_family) +
  theme(
    plot.title = element_text(size = 10, face = "bold", hjust = 0),
    axis.text.x = element_text(
      angle = 60,
      hjust = 1,
      vjust = 1,
      size = 6
    ),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 8),
    strip.background = element_blank(),
    strip.text = element_text(size = 7, face = "bold"),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 7),
    legend.position = "right",
    legend.box = "vertical",
    legend.box.just = "left",
    panel.grid = element_blank(),
    panel.spacing.x = grid::unit(0.4, "lines"),
    plot.margin = margin(t = 8, r = 18, b = 6, l = 6)
  ) +
  guides(
    color = guide_colorbar(
      title.position = "top",
      barheight = grid::unit(2.8, "cm"),
      barwidth = grid::unit(0.45, "cm")
    ),
    size = guide_legend(
      title.position = "top",
      override.aes = list(alpha = 1)
    )
  )
# ==========================================================
# SAVE
# ==========================================================
figure_pdf <- file.path(
  out_figure_dir,
  "germ_subcluster_stage_validation_panel.pdf"
)
figure_png <- file.path(
  out_figure_dir,
  "germ_subcluster_stage_validation_panel.png"
)

save_pdf_png(
  plot_obj = p,
  pdf_file = figure_pdf,
  png_file = figure_png,
  width = figure_width,
  height = figure_height,
  dpi = figure_dpi
)

summary_lines <- c(
  paste("Input object:", input_rds),
  paste("Output PDF:", figure_pdf),
  paste("Output PNG:", figure_png),
  "",
  "Subcluster order (top to bottom):",
  paste(present_order_top_to_bottom, collapse = ", "),
  "",
  "Display labels (top to bottom):",
  paste(label_map$display_label[match(present_order_top_to_bottom, label_map$fine_subcluster)], collapse = " | "),
  "",
  paste("Genes requested:", length(feature_order)),
  paste("Genes present:", length(genes_present)),
  paste("Genes missing:", length(genes_missing))
)

write_lines(
  summary_lines,
  file.path(out_log_dir, "germ_stage_validation_panel_summary.txt")
)

sink(file.path(out_log_dir, "germ_stage_validation_panel_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("\nDone.\n")
cat("Figure written to:\n  ", figure_pdf, "\n  ", figure_png, "\n", sep = "")
cat("Results root:\n  ", results_root, "\n", sep = "")
