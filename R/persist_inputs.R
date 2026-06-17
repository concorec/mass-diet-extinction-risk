# Input loading, C++ compilation, and CRN setup for persistence-point simulations.

check_persist_files <- function(paths, run_mammals, run_birds) {
  need_file(paths$cpp_file, "C++ simulator file")
  ensure_writable_dir(paths$results_dir, "results directory")
  ensure_writable_dir(paths$rcpp_cache, "Rcpp cache directory")

  if (isTRUE(run_mammals)) {
    lapply(c(paths$mammal_mass_grid, paths$post_m_rm, paths$post_m_sigma), need_file)
  }
  if (isTRUE(run_birds)) {
    lapply(c(paths$bird_mass_grid, paths$bird_ilr_coords, paths$post_b_rm, paths$post_b_sigma), need_file)
  }

  invisible(TRUE)
}

compile_persist_cpp <- function(paths, rebuild_cpp = FALSE) {
  cache_dir <- paths$rcpp_cache
  if (isTRUE(rebuild_cpp)) {
    cache_dir <- tempfile("sourceCpp_", tmpdir = paths$rcpp_cache)
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }

  Rcpp::sourceCpp(
    file = paths$cpp_file,
    cacheDir = cache_dir,
    rebuild = isTRUE(rebuild_cpp),
    verbose = FALSE
  )

  assert(exists("crn_context_create"), "C++ symbol missing: crn_context_create")
  assert(exists("simulate_persist_qvec_cpp"), "C++ symbol missing: simulate_persist_qvec_cpp")
  assert(exists("simulate_persist_probs_cpp"), "C++ symbol missing: simulate_persist_probs_cpp")
  assert(exists("persist_cpp_contract"), "C++ symbol missing: persist_cpp_contract")
  assert(
    identical(as.character(persist_cpp_contract()), "initial-threshold-v2"),
    "Compiled C++ simulator is stale. Rebuild src/simulate_persist_probs_cpp.cpp."
  )

  invisible(TRUE)
}

load_mammal_demographic_inputs <- function(paths) {
  hdr_mass <- data.table::fread(paths$mammal_mass_grid, nrows = 0)
  hdr_rm <- data.table::fread(paths$post_m_rm, nrows = 0)
  hdr_sg <- data.table::fread(paths$post_m_sigma, nrows = 0)

  assert_has_cols(hdr_mass, "Mass_g", "mammal_mass_grid")
  assert_has_cols(hdr_rm, c("alpha", "beta_logM"), "post_mammal_rm_coefs")
  assert_has_cols(hdr_sg, c("alpha", "beta_logM"), "post_mammal_sigma_coefs")

  list(
    mass = data.table::fread(paths$mammal_mass_grid, select = "Mass_g")[["Mass_g"]],
    post_rm = data.table::fread(paths$post_m_rm, select = c("alpha", "beta_logM")),
    post_sigma = data.table::fread(paths$post_m_sigma, select = c("alpha", "beta_logM"))
  )
}

load_bird_demographic_inputs <- function(paths) {
  hdr_mass <- data.table::fread(paths$bird_mass_grid, nrows = 0)
  hdr_ilr <- data.table::fread(paths$bird_ilr_coords, nrows = 0)
  hdr_rm <- data.table::fread(paths$post_b_rm, nrows = 0)
  hdr_sg <- data.table::fread(paths$post_b_sigma, nrows = 0)

  assert_has_cols(hdr_mass, "Mass_g", "bird_mass_grid")
  assert_has_cols(hdr_ilr, c("combo", "ilr1", "ilr2"), "bird_ilr_combo_coords")
  assert_has_cols(hdr_rm, c("alpha", "beta_logM"), "post_bird_r_coefs")
  assert_has_cols(hdr_sg, c("alpha", "beta_logM", "beta_ilr1", "beta_ilr2"), "post_bird_sigma_ilr_coefs")

  bird_ilr <- data.table::fread(paths$bird_ilr_coords, select = c("combo", "ilr1", "ilr2"))
  assert(
    bird_ilr[, data.table::uniqueN(combo)] == nrow(bird_ilr),
    paste0("Duplicate combo values found in ", paths$bird_ilr_coords, ".")
  )
  validate_bird_ilr_grid(bird_ilr, paths$bird_ilr_coords)

  list(
    mass = data.table::fread(paths$bird_mass_grid, select = "Mass_g")[["Mass_g"]],
    ilr = bird_ilr,
    post_rm = data.table::fread(paths$post_b_rm, select = c("alpha", "beta_logM")),
    post_sigma = data.table::fread(paths$post_b_sigma, select = c("alpha", "beta_logM", "beta_ilr1", "beta_ilr2"))
  )
}

validate_bird_ilr_grid <- function(bird_ilr, label) {
  assert(nrow(bird_ilr) == 66L, paste0(label, " must contain the 66-row 10% ternary diet grid."))

  parts <- data.table::tstrsplit(bird_ilr$combo, "_", fixed = TRUE, type.convert = TRUE)
  assert(length(parts) == 3L, paste0(label, " combo values must have format inv_plants_vertfishscav."))

  part_dt <- data.table::data.table(inv = parts[[1]], plants = parts[[2]], vertfishscav = parts[[3]])
  assert(
    all(vapply(part_dt, is.numeric, logical(1))) &&
      all(part_dt >= 0) &&
      all(part_dt %% 10 == 0) &&
      all(rowSums(part_dt) == 100),
    paste0(label, " combo values must be non-negative 10% increments summing to 100.")
  )

  invisible(TRUE)
}

sample_posterior_indices <- function(n_rm, n_sigma, n_draws, base_seed, rm_offset, sigma_offset) {
  set.seed(base_seed + rm_offset)
  idx_rm <- sample.int(n_rm, n_draws, replace = TRUE)
  set.seed(base_seed + sigma_offset)
  idx_sigma <- sample.int(n_sigma, n_draws, replace = TRUE)
  list(rm = idx_rm, sigma = idx_sigma)
}

create_crn_context <- function(sim) {
  crn_context_create(
    seed = as.integer(sim$base_seed),
    n_draws = as.integer(sim$n_draws),
    reps = as.integer(sim$reps),
    years = as.integer(sim$years),
    chunk_size = as.integer(sim$chunk_size)
  )
}

run_crn_self_test <- function(sim) {
  tiny_ctx <- crn_context_create(
    seed = as.integer(sim$base_seed),
    n_draws = 2L,
    reps = 3L,
    years = 2L,
    chunk_size = 2L
  )

  q1 <- simulate_persist_qvec_cpp(
    r = c(0.1, 0.2),
    sigma = c(0.2, 0.3),
    K = 1000,
    ext_thr = as.integer(sim$quasi_extinction_abundance),
    cap_factor = as.double(sim$cap_factor),
    crn_ctx = tiny_ctx
  )
  q2 <- simulate_persist_qvec_cpp(
    r = c(0.1, 0.2),
    sigma = c(0.2, 0.3),
    K = 1000,
    ext_thr = as.integer(sim$quasi_extinction_abundance),
    cap_factor = as.double(sim$cap_factor),
    crn_ctx = tiny_ctx
  )

  assert(isTRUE(all.equal(q1, q2, tolerance = 0)),
         "CRN self-test failed: repeated calls did not match exactly.")

  q_thr <- simulate_persist_qvec_cpp(
    r = c(0.1, 0.2),
    sigma = c(0.2, 0.3),
    K = as.double(sim$quasi_extinction_abundance),
    ext_thr = as.integer(sim$quasi_extinction_abundance),
    cap_factor = as.double(sim$cap_factor),
    crn_ctx = tiny_ctx
  )
  q_below <- simulate_persist_qvec_cpp(
    r = c(0.1, 0.2),
    sigma = c(0.2, 0.3),
    K = as.double(sim$quasi_extinction_abundance - 1L),
    ext_thr = as.integer(sim$quasi_extinction_abundance),
    cap_factor = as.double(sim$cap_factor),
    crn_ctx = tiny_ctx
  )

  assert(all(q_thr == 0) && all(q_below == 0),
         "CRN self-test failed: K <= ext_thr must have zero persistence.")
}

load_persist_inputs <- function(paths, params, sim) {
  run_mammals <- isTRUE(params$run_mammals)
  run_birds <- isTRUE(params$run_birds)

  check_persist_files(paths, run_mammals, run_birds)
  compile_persist_cpp(paths, rebuild_cpp = params$rebuild_cpp)

  out <- list(mammals = NULL, birds = NULL)

  if (run_mammals) {
    m <- load_mammal_demographic_inputs(paths)
    m$idx <- sample_posterior_indices(
      n_rm = nrow(m$post_rm),
      n_sigma = nrow(m$post_sigma),
      n_draws = sim$n_draws,
      base_seed = sim$base_seed,
      rm_offset = 7001L,
      sigma_offset = 8001L
    )
    out$mammals <- m
  }

  if (run_birds) {
    b <- load_bird_demographic_inputs(paths)
    b$idx <- sample_posterior_indices(
      n_rm = nrow(b$post_rm),
      n_sigma = nrow(b$post_sigma),
      n_draws = sim$n_draws,
      base_seed = sim$base_seed,
      rm_offset = 9001L,
      sigma_offset = 10001L
    )
    out$birds <- b
  }

  out
}

validate_points_schema <- function(path, label, min_k = NULL) {
  need_file(path, label)
  hdr <- data.table::fread(path, nrows = 0)
  assert_has_cols(hdr, c("combo", "mass_idx", "mass_g", "curve", "K", "p_eval"), label)

  if (!is.null(min_k)) {
    k_vals <- data.table::fread(path, select = "K")[["K"]]
    bad_k <- !is.finite(k_vals) | k_vals <= min_k
    assert(
      !any(bad_k),
      paste0(
        label, " contains ", sum(bad_k),
        " K value(s) at or below the quasi-extinction threshold of ", min_k,
        ". Persistence points must be generated only for K > Nq."
      )
    )
  }

  invisible(TRUE)
}
