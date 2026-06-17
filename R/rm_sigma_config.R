# Configuration and validation helpers for 1_rm_sigma_models.Rmd.

source(file.path("R", "project_utils.R"))

bayes_config <- function(params) {
  list(
    tau_beta  = params$tau_beta,
    sigma_min = params$sigma_min,
    sigma_max = params$sigma_max,
    n_chains  = params$n_chains,
    n_adapt   = params$n_adapt,
    n_iter    = params$n_iter,
    thin      = params$thin,
    seed      = params$seed
  )
}

required_rm_sigma_packages <- function(run_fits) {
  pkgs <- c(
    "dplyr", "readr", "tibble", "compositions", "zCompositions",
    "ggplot2", "cowplot", "ggtern", "scales", "viridisLite"
  )
  if (isTRUE(run_fits)) pkgs <- c(pkgs, "rjags", "coda")
  unique(pkgs)
}

assert_raw_inputs_exist <- function(paths) {
  missing_files <- unlist(paths$raw)[!file.exists(unlist(paths$raw))]
  if (length(missing_files) > 0) {
    stop("Missing required input file(s):\n", paste(missing_files, collapse = "\n"), call. = FALSE)
  }
  invisible(TRUE)
}

assert_posterior_inputs_exist <- function(paths) {
  needed_posts <- unlist(paths$clean[c(
    "post_mammal_rm", "post_mammal_sigma", "post_bird_rm", "post_bird_sigma_ilr3"
  )])
  missing_posts <- needed_posts[!file.exists(needed_posts)]

  if (length(missing_posts) > 0) {
    stop(
      "run_fits is FALSE, but required posterior CSV(s) are missing:\n",
      paste(missing_posts, collapse = "\n"),
      call. = FALSE
    )
  }

  invisible(TRUE)
}

assert_posterior_schema <- function(paths, expected_draws = 1500L) {
  specs <- list(
    mammal_rm = list(
      path = paths$clean$post_mammal_rm,
      cols = c("alpha", "beta_logM", "sigma")
    ),
    mammal_sigma = list(
      path = paths$clean$post_mammal_sigma,
      cols = c("alpha", "beta_logM", "sigma")
    ),
    bird_rm = list(
      path = paths$clean$post_bird_rm,
      cols = c("alpha", "beta_logM", "sigma")
    ),
    bird_sigma = list(
      path = paths$clean$post_bird_sigma_ilr3,
      cols = c("alpha", "beta_logM", "beta_ilr1", "beta_ilr2", "sigma")
    )
  )

  for (nm in names(specs)) {
    spec <- specs[[nm]]
    need_file(spec$path, nm)
    header <- readr::read_csv(spec$path, n_max = 0, show_col_types = FALSE)
    assert_has_cols(header, spec$cols, nm)

    n_rows <- nrow(readr::read_csv(spec$path, show_col_types = FALSE))
    if (!identical(as.integer(n_rows), as.integer(expected_draws))) {
      stop(
        nm, " should contain ", expected_draws, " posterior draws, but contains ", n_rows, ".",
        call. = FALSE
      )
    }
  }

  invisible(TRUE)
}
