# Volcano plots for fetal-vs-adult pseudobulk edgeR results.
#
# This script reads the standardized pseudobulk edgeR output tables and creates
# one volcano plot per successfully tested shared cell type.

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggrepel)
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


source(file.path(project_root, "config", "plotting.R"))
source(file.path(project_root, "config", "labels_colors.R"))


deg_results_root <- file.path(results_base, "fetal_adult_pseudobulk_edgeR")
deg_table_dir <- file.path(deg_results_root, "tables")

results_root <- file.path(results_base, "fetal_adult_pseudobulk_volcano_plots")

out_figure_dir <- file.path(results_root, "figures")
out_table_dir <- file.path(results_root, "tables")
out_log_dir <- file.path(results_root, "logs")

dir.create(out_figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log_dir, recursive = TRUE, showWarnings = FALSE)

qc_summary_csv <- file.path(
  deg_table_dir,
  "PSEUDOBULK_edgeR_QC_summary.csv"
)

fdr_cut <- 0.05
lfc_cut <- 0.25
label_top_each_side <- 10

# Keep this explicit, even though the QC file should already skip germ.
skip_celltypes <- c("germ")

# Manuscript colors. Degenerated is intentionally light gray for this manuscript.
volcano_celltype_colors <- c(
  celltype_colors,
  degenerated = "#CFCFCF",
  epithelial = "#BFE6E2",
  theca = "#FB8072",
  "glial cell" = "#F04E37",
  "plasma cell" = "#00E41A",
  "smooth muscle cell" = "#1F2BFF"
)

lighten_hex <- function(hex, amount = 0.65) {
  hex <- gsub("^#", "", hex)

  if (nchar(hex) != 6) {
    stop("Invalid hex color: ", hex)
  }

  r <- strtoi(substr(hex, 1, 2), 16L)
  g <- strtoi(substr(hex, 3, 4), 16L)
  b <- strtoi(substr(hex, 5, 6), 16L)

  r2 <- round(r + (255 - r) * amount)
  g2 <- round(g + (255 - g) * amount)
  b2 <- round(b + (255 - b) * amount)

  sprintf("#%02X%02X%02X", r2, g2, b2)
}

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

save_publication_plot_pdf_png <- function(plot, pdf_file, png_file, width, height, dpi = 600) {
  ggsave(
    pdf_file,
    plot,
    width = width,
    height = height,
    device = cairo_pdf,
    bg = "white",
    limitsize = FALSE
  )

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(
      png_file,
      plot,
      width = width,
      height = height,
      dpi = dpi,
      device = ragg::agg_png,
      background = "white",
      limitsize = FALSE
    )
  } else {
    ggsave(
      png_file,
      plot,
      width = width,
      height = height,
      dpi = dpi,
      bg = "white",
      limitsize = FALSE
    )
  }
}
prepare_deg_table <- function(tab) {
  required_cols <- c("celltype", "gene", "logFC", "FDR")
  missing_cols <- setdiff(required_cols, colnames(tab))

  if (length(missing_cols) > 0) {
    stop("DE table missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  tab |>
    mutate(
      celltype = as.character(celltype),
      gene = as.character(gene),
      logFC = as.numeric(logFC),
      FDR = as.numeric(FDR)
    ) |>
    filter(!is.na(gene), gene != "", !is.na(logFC), !is.na(FDR)) |>
    mutate(
      neglog10FDR = -log10(pmax(FDR, 1e-300)),

      # The underlying edgeR table stores logFC as fetal vs. adult.
      # For manuscript volcano plots, flip the plotted x-axis so that
      # fetal-enriched genes appear on the left and adult-enriched genes
      # appear on the right.
      plot_logFC = -logFC,

      status = case_when(
        FDR < fdr_cut & plot_logFC <= -lfc_cut ~ "Higher in fetal",
        FDR < fdr_cut & plot_logFC >=  lfc_cut ~ "Higher in adult",
        TRUE ~ "Not significant"
      )
    )
}

make_volcano <- function(tab, celltype, xmax, ymax) {
  celltype_key <- as.character(celltype)

  base_color <- volcano_celltype_colors[[celltype_key]]
  if (is.null(base_color) || is.na(base_color)) {
    base_color <- "#4D4D4D"
  }

  light_color <- lighten_hex(base_color, amount = 0.70)

  plot_df <- tab |>
    mutate(
      point_color_group = if_else(status == "Not significant", "Not significant", "Significant"),
      status = factor(
        status,
        levels = c("Higher in fetal", "Higher in adult", "Not significant")
      )
    )

  label_df <- bind_rows(
    plot_df |>
      filter(status == "Higher in fetal") |>
      arrange(FDR, desc(abs(logFC))) |>
      slice_head(n = label_top_each_side),
    plot_df |>
      filter(status == "Higher in adult") |>
      arrange(FDR, desc(abs(logFC))) |>
      slice_head(n = label_top_each_side)
  ) |>
    distinct(gene, .keep_all = TRUE)
  ggplot(plot_df, aes(x = plot_logFC, y = neglog10FDR)) +
    geom_point(
      aes(color = point_color_group, shape = status),
      alpha = 0.85,
      size = 0.75,
      stroke = 0
    ) +
    geom_vline(
      xintercept = c(-lfc_cut, lfc_cut),
      linetype = "dashed",
      linewidth = 0.25,
      color = "grey45"
    ) +
    geom_hline(
      yintercept = -log10(fdr_cut),
      linetype = "dashed",
      linewidth = 0.25,
      color = "grey45"
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(label = gene),
      family = publication_font_family,
      size = 2.1,
      max.overlaps = Inf,
      box.padding = 0.25,
      point.padding = 0.15,
      min.segment.length = 0,
      linewidth = 0.15,
      seed = 42
    ) +
    coord_cartesian(
      xlim = c(-xmax, xmax),
      ylim = c(0, ymax)
    ) +
    scale_color_manual(
      values = c(
        "Not significant" = light_color,
        "Significant" = base_color
      ),
      breaks = c("Significant", "Not significant"),
      name = NULL
    ) +
    scale_shape_manual(
      values = c(
        "Higher in fetal" = 16,
        "Higher in adult" = 16,
        "Not significant" = 16
      ),
      drop = FALSE,
      name = NULL
    ) +
    labs(
      x = "logFC, adult vs. fetal",
      y = expression(-log[10]("FDR"))
    ) +
    theme_publication() +
    theme(
      legend.position = "none"
    )
}

stopifnot(dir.exists(deg_table_dir))
stopifnot(file.exists(qc_summary_csv))

qc_summary <- read_csv(qc_summary_csv, show_col_types = FALSE)

required_qc_cols <- c("celltype", "status", "output_file")
missing_qc_cols <- setdiff(required_qc_cols, colnames(qc_summary))

if (length(missing_qc_cols) > 0) {
  stop("QC summary missing required columns: ", paste(missing_qc_cols, collapse = ", "))
}
plot_qc <- qc_summary |>
  mutate(
    celltype = as.character(celltype),
    status = as.character(status),
    output_file = as.character(output_file)
  ) |>
  filter(
    status == "ok_saved",
    !tolower(celltype) %in% skip_celltypes
  )

if (nrow(plot_qc) == 0) {
  stop("No successfully tested cell types available for volcano plotting.")
}

missing_files <- plot_qc$output_file[!file.exists(plot_qc$output_file)]
if (length(missing_files) > 0) {
  stop("Missing DE result files:\n", paste(missing_files, collapse = "\n"))
}

deg_tables <- lapply(plot_qc$output_file, function(file) {
  prepare_deg_table(read_csv(file, show_col_types = FALSE))
})

names(deg_tables) <- plot_qc$celltype

all_deg <- bind_rows(deg_tables)

if (nrow(all_deg) == 0) {
  stop("No DE rows available after filtering.")
}

xmax <- max(abs(all_deg$plot_logFC), na.rm = TRUE)
xmax <- ceiling(xmax * 10) / 10

ymax <- max(all_deg$neglog10FDR, na.rm = TRUE)
ymax <- ceiling(ymax * 10) / 10

write_csv(
  tibble(
    fdr_cut = fdr_cut,
    lfc_cut = lfc_cut,
    label_top_each_side = label_top_each_side,
    xmax = xmax,
    ymax = ymax,
    x_axis_note = "plot_logFC = -edgeR_logFC; negative values indicate higher fetal expression; positive values indicate higher adult expression"
  ),
  file.path(out_table_dir, "PSEUDOBULK_edgeR_volcano_plot_settings.csv")
)

plot_summary_rows <- list()
label_rows <- list()

for (celltype in names(deg_tables)) {
  message("Creating volcano plot for: ", celltype)

  tab <- deg_tables[[celltype]]

  p <- make_volcano(
    tab = tab,
    celltype = celltype,
    xmax = xmax,
    ymax = ymax
  )

  file_prefix <- paste0(
    "PSEUDOBULK_edgeR_volcano_fetal_vs_adult__",
    sanitize_filename(celltype)
  )
  pdf_file <- file.path(out_figure_dir, paste0(file_prefix, ".pdf"))
  png_file <- file.path(out_figure_dir, paste0(file_prefix, ".png"))

  save_publication_plot_pdf_png(
    plot = p,
    pdf_file = pdf_file,
    png_file = png_file,
    width = 3.2,
    height = 3.0,
    dpi = 600
  )

  label_df <- bind_rows(
    tab |>
      filter(status == "Higher in fetal") |>
      arrange(FDR, desc(abs(logFC))) |>
      slice_head(n = label_top_each_side),
    tab |>
      filter(status == "Higher in adult") |>
      arrange(FDR, desc(abs(logFC))) |>
      slice_head(n = label_top_each_side)
  ) |>
    distinct(gene, .keep_all = TRUE) |>
    mutate(celltype = celltype) |>
    select(celltype, gene, logFC, plot_logFC, FDR, status)

  label_rows[[length(label_rows) + 1]] <- label_df

  plot_summary_rows[[length(plot_summary_rows) + 1]] <- tibble(
    celltype = celltype,
    n_genes = nrow(tab),
    n_higher_in_fetal = sum(tab$status == "Higher in fetal", na.rm = TRUE),
    n_higher_in_adult = sum(tab$status == "Higher in adult", na.rm = TRUE),
    n_not_significant = sum(tab$status == "Not significant", na.rm = TRUE),
    color = volcano_celltype_colors[[celltype]],
    pdf_file = pdf_file,
    png_file = png_file
  )
}

plot_summary <- bind_rows(plot_summary_rows)

write_csv(
  plot_summary,
  file.path(out_table_dir, "PSEUDOBULK_edgeR_volcano_plot_summary.csv")
)

if (length(label_rows) > 0) {
  write_csv(
    bind_rows(label_rows),
    file.path(out_table_dir, "PSEUDOBULK_edgeR_volcano_labeled_genes.csv")
  )
}

summary_lines <- c(
  paste("Input pseudobulk table directory:", deg_table_dir),
  paste("Input QC summary:", qc_summary_csv),
  paste("Output directory:", results_root),
  paste("FDR cutoff:", fdr_cut),
  paste("logFC cutoff:", lfc_cut),
  paste("Label top genes per side:", label_top_each_side),
  paste("Skipped cell types:", paste(skip_celltypes, collapse = ", ")),
  paste("Plotted cell types:", paste(plot_summary$celltype, collapse = ", ")),
  paste("Shared x-axis maximum:", xmax),
  paste("Shared y-axis maximum:", ymax),
  "Volcano x-axis note: plot_logFC = -edgeR logFC; negative values indicate higher fetal expression; positive values indicate higher adult expression."
)

writeLines(
  summary_lines,
  file.path(out_log_dir, "PSEUDOBULK_edgeR_volcano_plots_summary.txt")
)

sink(file.path(out_log_dir, "PSEUDOBULK_edgeR_volcano_plots_sessionInfo.txt"))
print(sessionInfo())
sink()

cat("Wrote volcano plots to:\n  ", out_figure_dir, "\n", sep = "")
