options(bitmapType = "cairo")

publication_font_family <- "Helvetica"
publication_base_size <- 8

theme_publication <- function(base_size = publication_base_size,
                              base_family = publication_font_family) {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      text = ggplot2::element_text(family = base_family, size = base_size),
      axis.text = ggplot2::element_text(color = "black", size = base_size),
      axis.title = ggplot2::element_text(color = "black", size = base_size),
      plot.title = ggplot2::element_text(color = "black", size = base_size, hjust = 0),
      legend.title = ggplot2::element_text(color = "black", size = base_size),
      legend.text = ggplot2::element_text(color = "black", size = base_size),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(color = "black", size = base_size),
      panel.grid = ggplot2::element_blank()
    )
}

save_publication_plot <- function(plot,
                                  filename,
                                  width,
                                  height,
                                  dpi = 600) {
  dir.create(dirname(filename), recursive = TRUE, showWarnings = FALSE)

  pdf_file <- sub("\\.[A-Za-z0-9]+$", ".pdf", filename)
  png_file <- sub("\\.[A-Za-z0-9]+$", ".png", filename)

  ggplot2::ggsave(
    filename = pdf_file,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE
  )

  if (requireNamespace("ragg", quietly = TRUE)) {
    ggplot2::ggsave(
      filename = png_file,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi,
      device = ragg::agg_png
    )
  } else {
    ggplot2::ggsave(
      filename = png_file,
      plot = plot,
      width = width,
      height = height,
      units = "in",
      dpi = dpi
    )
  }

  invisible(c(pdf = pdf_file, png = png_file))
}
