# Configuration and validation helpers for 3_gompertz_loess_models.Rmd.

source(file.path("R", "project_utils.R"))
source(file.path("R", "analysis_contract.R"))

gomp_abort <- function(...) stop(paste0(...), call. = FALSE)

gompertz_curves <- function() {
  persistence_curves()
}

required_gompertz_packages <- function() {
  c(
    "dplyr", "readr", "tibble", "ggplot2", "cowplot", "scales",
    "viridisLite", "compositions", "zCompositions"
  )
}

load_gompertz_packages <- function() {
  load_packages(required_gompertz_packages())
}

validate_gompertz_config <- function(params, curves, loess_span, loess_z) {
  if (!isTRUE(params$run_mammals) && !isTRUE(params$run_birds)) {
    gomp_abort("At least one of run_mammals or run_birds must be TRUE.")
  }
  if (!is.finite(loess_span) || loess_span <= 0 || loess_span > 1) {
    gomp_abort("loess_span must be in (0, 1].")
  }
  if (!is.finite(loess_z) || loess_z <= 0) {
    gomp_abort("loess_z must be positive.")
  }
  invisible(TRUE)
}

check_gompertz_inputs <- function(paths, run_mammals, run_birds, write_models) {
  points_cols <- c("combo", "mass_idx", "mass_g", "curve", "K", "p_eval")

  if (isTRUE(run_mammals)) {
    need_file(paths$mammals_points, "mammal persistence points")
    need_cols(readr::read_csv(paths$mammals_points, n_max = 1, show_col_types = FALSE),
              points_cols, "persist_points_mammals.csv")
  }
  if (isTRUE(run_birds)) {
    need_file(paths$birds_points, "bird persistence points")
    need_file(paths$bird_ilr_coords, "bird ILR coordinates")
    need_file(paths$bird_sigma_draws, "bird sigma posterior draws")
    need_cols(readr::read_csv(paths$birds_points, n_max = 1, show_col_types = FALSE),
              points_cols, "persist_points_birds.csv")
  }
  if (isTRUE(write_models)) ensure_dir(dirname(paths$loess_models_rds))

  invisible(TRUE)
}
