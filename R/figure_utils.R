# Canonical manuscript-figure outputs.

manuscript_figure_manifest <- function() {
  tibble::tribble(
    ~id, ~filename, ~width, ~height, ~dpi,
    "bird_sigma", "bird_sigma.png", 7.287, 4.5, 600,
    "allometry", "S2-allometric-calibration-relationships.png", 7.287, 4.5, 600,
    "bird_sigma_models", "S3-bird-sigma-model-comparison.png", 7.287, 4.5, 600,
    "gompertz_wolff", "gompertz_wolff.png", 7.287, 4.5, 600,
    "mammal_loess", "S3-mammal-loess-parameters-all-curves.png", 7.287, 4.5, 600,
    "bird_loess", "S3-bird-loess-parameters-by-diet.png", 7.287, 4.5, 600,
    "area_curve", "area_curve.png", 6.4, 3.3, 600,
    "spatial_process", "spatial_process.png", 8.5, 9.1, 300,
    "persistence_comparison", "results_comparison_4panel.png", 7.0, 6.0, 600
  )
}

save_manuscript_figure <- function(plot, id, figure_dir = "Figures") {
  assert(inherits(plot, c("gg", "ggplot", "grob", "gtable")), paste0("Figure '", id, "' is not a ggplot-compatible object."))
  manifest <- manuscript_figure_manifest()
  spec <- manifest[manifest$id == id, , drop = FALSE]
  assert(nrow(spec) == 1L, paste0("Unknown manuscript figure id: ", id))

  ensure_dir(figure_dir)
  path <- file.path(figure_dir, spec$filename)
  ggplot2::ggsave(
    filename = path,
    plot = plot,
    width = spec$width,
    height = spec$height,
    dpi = spec$dpi,
    units = "in",
    bg = "white"
  )
  need_file(path, paste0("saved figure '", id, "'"))
  assert(file.info(path)$size > 0, paste0("Saved figure is empty: ", path))
  message("Wrote figure: ", path)
  invisible(path)
}
