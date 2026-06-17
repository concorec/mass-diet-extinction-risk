# Orchestrator for 7.3_persist_cmp.Rmd.
# The implementation is split into focused calculation, reporting, and figure
# files while the Rmd keeps one stable source() call.

stage7_compare_env <- environment()

source_stage7_compare_file <- function(path, quiet = TRUE) {
  run_source <- function() {
    source(path, local = stage7_compare_env, print.eval = !isTRUE(quiet))
  }

  if (isTRUE(quiet)) {
    invisible(capture.output(run_source(), type = "output"))
  } else {
    run_source()
  }

  invisible(NULL)
}

source_stage7_compare_file(file.path("R", "stage7_persistence_compare_core.R"), quiet = TRUE)
source_stage7_compare_file(file.path("R", "stage7_persistence_compare_stats.R"), quiet = TRUE)
source_stage7_compare_file(file.path("R", "stage7_persistence_compare_figures.R"), quiet = TRUE)
source_stage7_compare_file(file.path("R", "stage7_persistence_compare_report.R"), quiet = FALSE)

if (exists("fig_cmp_4panel", envir = stage7_compare_env, inherits = FALSE)) {
  source(file.path("R", "figure_utils.R"), local = stage7_compare_env)
  figure <- get("fig_cmp_4panel", envir = stage7_compare_env, inherits = FALSE)
  save_manuscript_figure(figure, "persistence_comparison")
  print(figure)
}
