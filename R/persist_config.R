# Configuration and validation helpers for 2_persist_points_crn.Rmd.

source(file.path("R", "project_utils.R"))
source(file.path("R", "analysis_contract.R"))

required_persist_packages <- function() {
  c("data.table", "Rcpp")
}

load_persist_packages <- function() {
  load_packages(required_persist_packages())
}

validate_persist_config <- function(params, sim, grid) {
  run_mammals <- isTRUE(params$run_mammals)
  run_birds <- isTRUE(params$run_birds)

  assert(is.logical(run_mammals) && length(run_mammals) == 1, "run_mammals must be TRUE/FALSE.")
  assert(is.logical(run_birds) && length(run_birds) == 1, "run_birds must be TRUE/FALSE.")
  assert(run_mammals || run_birds, "At least one of run_mammals or run_birds must be TRUE.")

  assert(
    sim$quasi_extinction_abundance > 0 && sim$years > 0,
    "quasi_extinction_abundance and years must be positive."
  )
  assert(sim$n_draws > 0 && sim$reps > 0 && sim$chunk_size > 0, "n_draws, reps, and chunk_size must be positive.")
  assert(
    grid$k_search$start_k > sim$quasi_extinction_abundance,
    "start_k must exceed the quasi-extinction threshold."
  )
  assert(all(grid$p_anchor > 0 & grid$p_anchor < 1) && all(diff(grid$p_anchor) > 0),
         "p_anchor must be strictly increasing in (0, 1).")
  assert(all(grid$p_grid > 0 & grid$p_grid < 1) && all(diff(grid$p_grid) > 0),
         "p_grid must be strictly increasing in (0, 1).")

  invisible(TRUE)
}

describe_persist_params <- function(params, source_label = "effective") {
  paste0(
    source_label,
    " | run_mammals=", isTRUE(params$run_mammals),
    " run_birds=", isTRUE(params$run_birds),
    " force_rerun_bird_combos=", isTRUE(params$force_rerun_bird_combos),
    " rebuild_cpp=", isTRUE(params$rebuild_cpp),
    " verbose=", isTRUE(params$verbose),
    " run_crn_self_test=", isTRUE(params$run_crn_self_test)
  )
}
