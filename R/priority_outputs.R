# Output writers for the priority-removal sequence.

if (!exists("assert", mode = "function")) {
  source(file.path("R", "project_utils.R"))
}

write_removal_order_raster <- function(removal_order_by_cell, template_raster, output_path) {
  assert(!is.null(template_raster), "template_raster must be supplied.")

  output_raster <- terra::rast(template_raster[[1L]])
  assert(
    terra::ncell(output_raster) == length(removal_order_by_cell),
    "Length of removal_order_by_cell does not match template_raster."
  )

  terra::values(output_raster) <- as.integer(removal_order_by_cell)
  terra::writeRaster(
    x = output_raster,
    filename = output_path,
    datatype = "INT4S",
    overwrite = TRUE
  )

  invisible(output_path)
}

write_rankmap_raster <- function(
  removal_order_by_cell,
  removal_events,
  template_raster,
  output_path,
  value_col = "cum_prop_cells_removed"
) {
  assert(!is.null(template_raster), "template_raster must be supplied.")
  assert(is.data.frame(removal_events), "removal_events must be a data.frame or data.table.")
  need_cols(removal_events, c("removal_step", value_col), "removal_events")

  output_raster <- terra::rast(template_raster[[1L]])
  assert(
    terra::ncell(output_raster) == length(removal_order_by_cell),
    "Length of removal_order_by_cell does not match template_raster."
  )

  rank_values <- rep(NA_real_, length(removal_order_by_cell))
  ranked_cells <- which(!is.na(removal_order_by_cell))

  if (length(ranked_cells)) {
    step_to_rank <- stats::setNames(
      as.numeric(removal_events[[value_col]]),
      as.character(as.integer(removal_events$removal_step))
    )
    rank_values[ranked_cells] <- step_to_rank[as.character(removal_order_by_cell[ranked_cells])]
    assert(
      all(is.finite(rank_values[ranked_cells])),
      "rankmap values could not be assigned for one or more ranked cells."
    )
  }

  finite_rank_values <- rank_values[is.finite(rank_values)]
  assert(
    !length(finite_rank_values) ||
      (min(finite_rank_values) >= 0 && max(finite_rank_values) <= 1),
    "rankmap values must be in [0, 1]."
  )

  terra::values(output_raster) <- rank_values
  terra::writeRaster(
    x = output_raster,
    filename = output_path,
    overwrite = TRUE,
    wopt = list(
      datatype = "FLT4S",
      gdal = c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6", "TILED=YES", "BIGTIFF=IF_SAFER")
    )
  )

  invisible(output_path)
}
