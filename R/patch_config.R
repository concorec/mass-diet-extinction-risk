source(file.path("R", "project_utils.R"))
source(file.path("R", "analysis_contract.R"))

required_patch_packages <- function(include_figures = FALSE) {
  pkgs <- c("terra", "sf", "igraph", "dplyr", "readr", "stringr", "tools")
  if (isTRUE(include_figures)) {
    pkgs <- c(pkgs, "ggplot2", "tidyterra", "cowplot", "rnaturalearth", "scales")
  }
  unique(pkgs)
}

load_patch_packages <- function(include_figures = FALSE) {
  pkgs <- required_patch_packages(include_figures)
  load_packages(pkgs)

  sf::sf_use_s2(TRUE)
  invisible(pkgs)
}

patch_paths <- function(params) {
  list(
    species_csv      = params$species_csv,
    landcover_tif    = params$landcover_tif,
    patch_dir        = params$patch_dir %||% "Data/Clean/Patches",
    patch_lookup_rds = params$patch_lookup_rds %||% "Data/Clean/all_patch_lookup.rds",
    connectivity_rds = params$connectivity_rds %||% "Data/Clean/all_connectivity.rds"
  )
}

patch_run_options <- function(params) {
  list(
    do_mammals              = isTRUE(params$do_mammals),
    do_birds                = isTRUE(params$do_birds),
    grass_dir               = resolve_optional_grass_dir(params$grass_dir),
    overwrite_patch_rasters = isTRUE(params$overwrite_patch_rasters),
    stop_on_species_error   = isTRUE(params$stop_on_species_error)
  )
}

patch_roi <- function(params) {
  terra::ext(params$roi_xmin, params$roi_xmax, params$roi_ymin, params$roi_ymax)
}

validate_patch_config <- function(paths, run_options, roi, production = TRUE) {
  assert(inherits(roi, "SpatExtent"), "ROI parameters did not create a valid SpatExtent.")
  need_file(paths$species_csv, "species table from script 4")
  need_file(paths$landcover_tif, "ESA CCI land-cover raster")
  if (!is.na(optional_path(run_options$grass_dir))) {
    need_dir(run_options$grass_dir, "GRASS GIS installation directory")
  }

  if (isTRUE(production)) {
    assert(
      run_options$do_mammals || run_options$do_birds,
      "At least one of params$do_mammals or params$do_birds must be TRUE."
    )
    dir.create(paths$patch_dir, recursive = TRUE, showWarnings = FALSE)
    need_dir(paths$patch_dir, "patch raster output directory")
    dir.create(dirname(paths$patch_lookup_rds), recursive = TRUE, showWarnings = FALSE)
    dir.create(dirname(paths$connectivity_rds), recursive = TRUE, showWarnings = FALSE)
  }

  init_faster_raster(run_options$grass_dir)
  invisible(TRUE)
}
