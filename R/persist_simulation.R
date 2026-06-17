# Simulation helpers for persistence-point generation.

make_mammal_sampler <- function(inputs, sim) {
  function(mass_g) {
    logM <- log10(mass_g)
    r <- sim$r_buffer * 10^(
      inputs$post_rm$alpha[inputs$idx$rm] +
        inputs$post_rm$beta_logM[inputs$idx$rm] * logM
    )
    s <- 10^(
      inputs$post_sigma$alpha[inputs$idx$sigma] +
        inputs$post_sigma$beta_logM[inputs$idx$sigma] * logM
    )
    list(r = as.numeric(r), sigma = as.numeric(s))
  }
}

make_bird_sampler <- function(inputs, sim, ilr1, ilr2) {
  force(ilr1)
  force(ilr2)
  function(mass_g) {
    logM <- log10(mass_g)
    r <- sim$r_buffer * 10^(
      inputs$post_rm$alpha[inputs$idx$rm] +
        inputs$post_rm$beta_logM[inputs$idx$rm] * logM
    )
    s <- 10^(
      inputs$post_sigma$alpha[inputs$idx$sigma] +
        inputs$post_sigma$beta_logM[inputs$idx$sigma] * logM +
        inputs$post_sigma$beta_ilr1[inputs$idx$sigma] * ilr1 +
        inputs$post_sigma$beta_ilr2[inputs$idx$sigma] * ilr2
    )
    list(r = as.numeric(r), sigma = as.numeric(s))
  }
}

eval_qvec_crn <- function(r, sigma, K, crn_ctx, sim) {
  simulate_persist_qvec_cpp(
    r = r,
    sigma = sigma,
    K = as.double(K),
    ext_thr = as.integer(sim$quasi_extinction_abundance),
    cap_factor = as.double(sim$cap_factor),
    crn_ctx = crn_ctx
  )
}

find_K_for_target <- function(p_target, get_q_at_K, low0, high0, rel_tol) {
  low <- as.numeric(low0)
  high <- as.numeric(high0)
  if (high < low) {
    stop("Bisection search requires high0 >= low0.", call. = FALSE)
  }

  while (get_q_at_K(high) < p_target) {
    low <- high
    high <- high * 2
  }
  hi_bracket <- high

  while ((high - low) > rel_tol * high) {
    mid <- (low + high) / 2
    if (get_q_at_K(mid) >= p_target) high <- mid else low <- mid
  }

  list(K = (low + high) / 2, hi_bracket = hi_bracket)
}

fit_gompertz_ab <- function(K_anchor, p_anchor, sim) {
  x <- pmax(as.numeric(K_anchor) - sim$quasi_extinction_abundance, 1e-12)
  p <- as.numeric(p_anchor)

  y <- log(-log(p))
  lx <- log(x)
  lm0 <- stats::lm(y ~ lx)
  b0 <- -as.numeric(stats::coef(lm0)[["lx"]])
  a0 <- exp(as.numeric(stats::coef(lm0)[["(Intercept)"]]))

  if (!is.finite(b0) || b0 <= 0) b0 <- 1
  if (!is.finite(a0) || a0 <= 0) a0 <- 1
  b0 <- max(1e-4, min(80, b0))
  a0 <- max(1e-12, a0)

  fml <- p ~ exp(-exp(loga) * x^(-exp(logb)))
  ctrl <- stats::nls.control(maxiter = 500, tol = 1e-7, minFactor = 1 / 2048, warnOnly = TRUE)

  mult <- c(0.5, 0.8, 1, 1.25, 1.6, 2, 3, 5)
  b_grid <- unique(pmax(1e-4, pmin(80, b0 * mult)))
  i_mid <- match(0.50, p_anchor)

  fit <- NULL
  for (b_try in b_grid) {
    a_try <- -log(p[i_mid]) * x[i_mid]^b_try
    if (!is.finite(a_try) || a_try <= 0) a_try <- a0

    fit <- tryCatch(
      stats::nls(
        fml,
        data = data.frame(p = p, x = x),
        start = list(loga = log(a_try), logb = log(b_try)),
        algorithm = "port",
        control = ctrl
      ),
      error = function(e) NULL
    )
    if (!is.null(fit)) break
  }

  if (is.null(fit)) stop("Gompertz pre-fit failed while building the K grid.", call. = FALSE)

  co <- stats::coef(fit)
  list(a = exp(as.numeric(co[["loga"]])), b = exp(as.numeric(co[["logb"]])))
}

invert_gompertz <- function(a, b, p_grid, sim, k_search) {
  p <- pmin(pmax(as.numeric(p_grid), 1e-6), 1 - 1e-8)
  K <- sim$quasi_extinction_abundance + (a / -log(p))^(1 / b)
  if (isTRUE(k_search$round_k)) round(K) else ceiling(K)
}

run_one_scenario <- function(masses_g, sampler, out_file, crn_ctx, combo, sim, grid,
                             append = FALSE, verbose = FALSE) {
  if (!isTRUE(append) && file.exists(out_file)) file.remove(out_file)
  wrote_any <- file.exists(out_file)

  write_points <- function(dt) {
    data.table::fwrite(dt, out_file, append = wrote_any, col.names = !wrote_any)
    wrote_any <<- TRUE
  }

  for (mi in seq_along(masses_g)) {
    mass <- masses_g[mi]
    pars <- sampler(mass)
    cache <- new.env(parent = emptyenv())

    evalK <- function(K, who = "") {
      key <- sprintf("%.17g", K)
      v <- cache[[key]]
      if (!is.null(v)) return(v)

      qv <- eval_qvec_crn(pars$r, pars$sigma, K, crn_ctx, sim)
      if (isTRUE(verbose)) {
        log_msg(
          "CRN", who,
          "| mass_idx=", mi,
          "mass_g=", formatC(mass, format = "f", digits = 0),
          "K=", formatC(K, format = "f", digits = 0),
          "q50=", sprintf("%.4f", qv["q50"]),
          "q16=", sprintf("%.4f", qv["q16"]),
          "q025=", sprintf("%.4f", qv["q025"]),
          flush = verbose
        )
      }

      cache[[key]] <- qv
      qv
    }

    for (qtag in names(grid$q_levels)) {
      K_anchor <- numeric(length(grid$p_anchor))
      hi_prev <- as.numeric(grid$k_search$start_k)

      for (ai in seq_along(grid$p_anchor)) {
        pt <- grid$p_anchor[ai]
        who <- paste0("| stage=anchor q=", qtag, " p_target=", sprintf("%.3f", pt))
        get_q <- function(K) evalK(K, who = who)[[qtag]]

        if (ai == 1L) {
          low0 <- hi_prev
          high0 <- hi_prev
        } else {
          low0 <- hi_prev / 2
          high0 <- hi_prev
        }

        res <- find_K_for_target(pt, get_q, low0 = low0, high0 = high0,
                                 rel_tol = grid$k_search$rel_tol)
        K_anchor[ai] <- res$K
        hi_prev <- res$hi_bracket
      }

      ab <- fit_gompertz_ab(K_anchor, grid$p_anchor, sim)
      K_guess <- invert_gompertz(ab$a, ab$b, grid$p_grid, sim, grid$k_search)
      if (any(!is.finite(K_guess) | K_guess <= sim$quasi_extinction_abundance)) {
        stop(
          "Gompertz prefit produced K values at or below the quasi-extinction threshold. ",
          "This indicates that the five-anchor prefit is outside the valid persistence domain.",
          call. = FALSE
        )
      }

      p_hat <- numeric(length(K_guess))
      for (ki in seq_along(K_guess)) {
        who <- paste0("| stage=grid q=", qtag, " p_target=", sprintf("%.3f", grid$p_grid[ki]))
        p_hat[ki] <- evalK(K_guess[ki], who = who)[[qtag]]
      }

      write_points(data.table::data.table(
        combo = combo,
        mass_idx = as.integer(mi),
        mass_g = as.numeric(mass),
        curve = qtag,
        K = as.numeric(K_guess),
        p_eval = as.numeric(p_hat)
      ))

      if (isTRUE(verbose)) {
        x_anchor <- as.numeric(K_anchor) - sim$quasi_extinction_abundance
        p_fit <- exp(-ab$a * x_anchor^(-ab$b))
        delta <- p_fit - as.numeric(grid$p_anchor)
        log_msg(
          "GOMP | stage=pregrid q=", qtag,
          "| mass_idx=", mi,
          "mass_g=", formatC(mass, format = "f", digits = 0),
          "| deltas(p_fit - p_target)=",
          paste(sprintf("p=%.3f:%+.5f", grid$p_anchor, delta), collapse = " "),
          flush = verbose
        )
      }
    }

    log_msg("DONE | mass_idx=", mi, "/", length(masses_g), " mass_g=",
            formatC(mass, format = "f", digits = 0), flush = verbose)
  }

  invisible(out_file)
}

expected_rows_per_combo <- function(masses_g, grid) {
  length(masses_g) * length(grid$q_levels) * length(grid$p_grid)
}

remove_partial_combos <- function(path, expected_rows) {
  if (!file.exists(path)) return(character(0))

  counts <- data.table::fread(path, select = "combo")[, .N, by = combo]
  partial <- counts[N != expected_rows, as.character(combo)]

  if (length(partial)) {
    log_msg("Removing partial combo rows from ", path, ": ", paste(partial, collapse = ", "))
    dt <- data.table::fread(path)
    dt <- dt[!combo %chin% partial]
    data.table::fwrite(dt, path)
  }

  partial
}

completed_combos <- function(path, expected_rows) {
  if (!file.exists(path)) return(character(0))
  counts <- data.table::fread(path, select = "combo")[, .N, by = combo]
  counts[N == expected_rows, as.character(combo)]
}

run_bird_combo_set <- function(combo_ids, bird_inputs, paths, sim, grid, crn_ctx,
                               force = FALSE, verbose = FALSE) {
  expected_rows <- expected_rows_per_combo(bird_inputs$mass, grid)
  remove_partial_combos(paths$birds_points, expected_rows)

  done <- completed_combos(paths$birds_points, expected_rows)
  combos_to_run <- combo_ids
  if (!isTRUE(force)) combos_to_run <- setdiff(combos_to_run, done)

  if (!length(combos_to_run)) {
    log_msg("BIRDS | all requested combos already complete in: ", paths$birds_points)
    return(invisible(paths$birds_points))
  }

  if (file.exists(paths$birds_points)) {
    existing_dt <- data.table::fread(paths$birds_points)
    existing_dt <- existing_dt[!combo %chin% combos_to_run]
    data.table::fwrite(existing_dt, paths$birds_points)
  }

  log_msg("BIRDS | expected rows per complete combo = ", expected_rows)
  log_msg("BIRDS | combos to simulate now: ", paste(combos_to_run, collapse = ", "))

  tmp_dir <- file.path(paths$results_dir, "_tmp_bird_persistence_points")
  dir.create(tmp_dir, recursive = TRUE, showWarnings = FALSE)

  for (combo_i in combos_to_run) {
    row_i <- bird_inputs$ilr[combo == combo_i]
    if (nrow(row_i) != 1L) stop("Expected one ILR row for combo ", combo_i, call. = FALSE)

    log_msg("BIRDS | combo start = ", combo_i)

    tmp_combo_file <- file.path(tmp_dir, paste0("persist_points_birds_", combo_i, ".csv"))
    if (file.exists(tmp_combo_file)) file.remove(tmp_combo_file)

    run_one_scenario(
      masses_g = bird_inputs$mass,
      sampler = make_bird_sampler(bird_inputs, sim, as.numeric(row_i$ilr1), as.numeric(row_i$ilr2)),
      out_file = tmp_combo_file,
      crn_ctx = crn_ctx,
      combo = combo_i,
      sim = sim,
      grid = grid,
      append = FALSE,
      verbose = verbose
    )

    tmp_dt <- data.table::fread(tmp_combo_file)
    if (nrow(tmp_dt) != expected_rows) {
      stop(
        "Temporary output for combo ", combo_i, " has ", nrow(tmp_dt),
        " rows, but expected ", expected_rows, ". Not appending to target file.",
        call. = FALSE
      )
    }

    data.table::fwrite(tmp_dt, paths$birds_points,
                       append = file.exists(paths$birds_points),
                       col.names = !file.exists(paths$birds_points))
    file.remove(tmp_combo_file)

    log_msg("BIRDS | combo complete = ", combo_i)
  }

  invisible(paths$birds_points)
}

finalize_bird_output <- function(bird_inputs, paths, grid) {
  bird_combos_all <- bird_inputs$ilr[order(combo), as.character(combo)]
  expected_rows <- expected_rows_per_combo(bird_inputs$mass, grid)

  bird_dt <- data.table::fread(paths$birds_points)
  bird_dt <- bird_dt[combo %chin% bird_combos_all]
  data.table::setorder(bird_dt, combo, mass_idx, curve, K)
  data.table::fwrite(bird_dt, paths$birds_points)

  combo_counts <- bird_dt[, .N, by = combo][order(combo)]
  incomplete <- combo_counts[N != expected_rows, as.character(combo)]
  missing <- setdiff(bird_combos_all, combo_counts$combo)

  if (length(missing) || length(incomplete)) {
    stop(
      "Bird persistence output is incomplete. Missing combos: ",
      paste(missing, collapse = ", "),
      "; incomplete combos: ",
      paste(incomplete, collapse = ", "),
      call. = FALSE
    )
  }

  combo_counts
}
