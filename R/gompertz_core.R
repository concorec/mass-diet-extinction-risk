# Core shifted-Gompertz and LOESS helpers for script 3.

clamp01 <- function(p) pmin(pmax(p, 1e-8), 1 - 1e-8)

predict_p <- function(alpha, beta, K, K0) {
  x <- pmax(K - K0, 1e-12)
  exp(-alpha * x^(-beta))
}

fit_gompertz_ab <- function(K, p_eval, K0) {
  df <- tibble::tibble(K = as.numeric(K), p = clamp01(as.numeric(p_eval))) |>
    dplyr::filter(is.finite(K), is.finite(p), K > K0)

  x <- pmax(df$K - K0, 1e-12)
  y <- log(-log(df$p))
  lx <- log(x)

  lm0 <- stats::lm(y ~ lx)
  beta0 <- -as.numeric(stats::coef(lm0)[["lx"]])
  alpha0 <- exp(as.numeric(stats::coef(lm0)[["(Intercept)"]]))
  if (!is.finite(beta0) || beta0 <= 0) beta0 <- 1
  if (!is.finite(alpha0) || alpha0 <= 0) alpha0 <- 1

  fml <- p ~ exp(-exp(log_alpha) * pmax(K - K0, 1e-12)^(-exp(log_beta)))
  ctrl <- stats::nls.control(maxiter = 400, tol = 1e-7, minFactor = 1 / 2048, warnOnly = TRUE)

  i_mid <- which.min(abs(df$p - 0.5))
  beta_grid <- unique(pmax(1e-4, pmin(80, beta0 * c(0.6, 0.85, 1, 1.25, 1.7, 2.5, 4))))
  fit <- NULL
  for (b_try in beta_grid) {
    a_try <- -log(df$p[i_mid]) * pmax(df$K[i_mid] - K0, 1e-12)^b_try
    if (!is.finite(a_try) || a_try <= 0) a_try <- alpha0
    fit <- tryCatch(
      stats::nls(
        fml, data = df,
        start = list(log_alpha = log(a_try), log_beta = log(b_try)),
        algorithm = "port", control = ctrl
      ),
      error = function(e) NULL
    )
    if (!is.null(fit)) break
  }

  if (is.null(fit)) {
    return(tibble::tibble(alpha = NA_real_, beta = NA_real_, fit_ok = FALSE,
                          n_pts = nrow(df), rmse_lin = NA_real_, fit_msg = "NLS failed"))
  }

  co <- stats::coef(fit)
  alpha_hat <- exp(as.numeric(co[["log_alpha"]]))
  beta_hat <- exp(as.numeric(co[["log_beta"]]))
  p_hat <- predict_p(alpha_hat, beta_hat, df$K, K0 = K0)
  rmse_lin <- sqrt(mean((log(-log(df$p)) - log(-log(clamp01(p_hat))))^2))

  tibble::tibble(alpha = alpha_hat, beta = beta_hat, fit_ok = TRUE,
                 n_pts = nrow(df), rmse_lin = rmse_lin, fit_msg = "")
}

read_persistence_points <- function(path, group_name, curves, k0) {
  readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::transmute(
      group = group_name,
      combo = as.character(combo),
      mass_idx = as.integer(mass_idx),
      mass_g = as.numeric(mass_g),
      curve = as.character(curve),
      K = as.numeric(K),
      p_eval = as.numeric(p_eval)
    ) |>
    dplyr::filter(curve %in% curves, is.finite(mass_g), mass_g > 0, is.finite(K), K > k0)
}

load_all_persistence_points <- function(paths, run_mammals, run_birds, curves, k0) {
  dplyr::bind_rows(
    if (isTRUE(run_mammals)) read_persistence_points(paths$mammals_points, "Mammals", curves, k0) else tibble::tibble(),
    if (isTRUE(run_birds)) read_persistence_points(paths$birds_points, "Birds", curves, k0) else tibble::tibble()
  )
}

validate_bird_combo_coverage <- function(curves_tbl, paths, run_birds) {
  if (!isTRUE(run_birds)) return(invisible(TRUE))

  ilr <- readr::read_csv(paths$bird_ilr_coords, show_col_types = FALSE, col_select = "combo") |>
    dplyr::mutate(combo = as.character(combo))
  bird_combos <- unique(curves_tbl$combo[curves_tbl$group == "Birds"])

  missing <- setdiff(ilr$combo, bird_combos)
  extra <- setdiff(bird_combos, ilr$combo)

  if (length(missing) || length(extra)) {
    gomp_abort(
      "Bird persistence points do not match the canonical bird ILR combo grid. ",
      "Missing combos: ", paste(missing, collapse = ", "),
      "; extra combos: ", paste(extra, collapse = ", "),
      ". Rerun 2_persist_points_crn.Rmd after regenerating the 66-combo ILR grid."
    )
  }

  invisible(TRUE)
}

summarise_persistence_points <- function(curves_tbl) {
  curves_tbl |>
    dplyr::group_by(group, combo, curve) |>
    dplyr::summarise(
      n_mass = dplyr::n_distinct(mass_idx),
      n_points = dplyr::n(),
      K_min = min(K),
      K_max = max(K),
      .groups = "drop"
    )
}

fit_all_gompertz_parameters <- function(curves_tbl, k0) {
  curves_tbl |>
    dplyr::group_by(group, combo, mass_idx, mass_g, curve) |>
    dplyr::group_modify(~fit_gompertz_ab(.x$K, .x$p_eval, K0 = k0)) |>
    dplyr::ungroup()
}

summarise_gompertz_fits <- function(gomp_params) {
  gomp_params |>
    dplyr::group_by(group, curve) |>
    dplyr::summarise(
      n_fits = dplyr::n(),
      n_failed = sum(!fit_ok),
      median_rmse_lin = stats::median(rmse_lin, na.rm = TRUE),
      max_rmse_lin = max(rmse_lin, na.rm = TRUE),
      .groups = "drop"
    )
}

fit_loess_logparam <- function(mass_g, param_pos, span = 0.75) {
  df <- tibble::tibble(logM = log10(as.numeric(mass_g)), lp = log(as.numeric(param_pos))) |>
    dplyr::filter(is.finite(logM), is.finite(lp))
  stats::loess(lp ~ logM, data = df, span = span, degree = 2, family = "gaussian")
}

predict_loess_exp <- function(fit, mass_g, z = 1.96) {
  nd <- data.frame(logM = log10(as.numeric(mass_g)))
  pr <- stats::predict(fit, newdata = nd, se = TRUE)
  tibble::tibble(
    mass_g = as.numeric(mass_g),
    mid = exp(pr$fit),
    lo = exp(pr$fit - z * pr$se.fit),
    hi = exp(pr$fit + z * pr$se.fit)
  )
}

fit_curve_loess <- function(d, span) {
  d <- d |>
    dplyr::filter(
      is.finite(mass_g), mass_g > 0,
      is.finite(alpha), alpha > 0,
      is.finite(beta), beta > 0
    )
  list(
    alpha = tryCatch(fit_loess_logparam(d$mass_g, d$alpha, span = span), error = function(e) NULL),
    beta = tryCatch(fit_loess_logparam(d$mass_g, d$beta, span = span), error = function(e) NULL),
    n = nrow(d)
  )
}

fit_loess_set <- function(d, curves, span) {
  stats::setNames(
    lapply(curves, function(cc) fit_curve_loess(dplyr::filter(d, curve == cc), span = span)),
    curves
  )
}

predict_one_curve <- function(fit_obj, mass_grid, curve_name, z) {
  if (is.null(fit_obj$alpha) || is.null(fit_obj$beta)) return(NULL)
  alpha <- predict_loess_exp(fit_obj$alpha, mass_grid, z = z) |>
    dplyr::transmute(mass_g, curve = curve_name, alpha_mid = mid, alpha_lo = lo, alpha_hi = hi)
  beta <- predict_loess_exp(fit_obj$beta, mass_grid, z = z) |>
    dplyr::transmute(mass_g, curve = curve_name, beta_mid = mid, beta_lo = lo, beta_hi = hi)
  dplyr::left_join(alpha, beta, by = c("mass_g", "curve"))
}

predict_loess_set <- function(d, fits, curves, z, n_grid = 700) {
  d <- d |>
    dplyr::filter(is.finite(mass_g), mass_g > 0, is.finite(alpha), alpha > 0, is.finite(beta), beta > 0)
  if (nrow(d) == 0) return(tibble::tibble())
  mass_grid <- log_space(min(d$mass_g), max(d$mass_g), n_grid)
  dplyr::bind_rows(lapply(curves, function(cc) predict_one_curve(fits[[cc]], mass_grid, cc, z = z)))
}

build_gompertz_loess_models <- function(gomp_params, paths, run_mammals, run_birds, curves, span, z, k0) {
  pm_mammals <- gomp_params |> dplyr::filter(group == "Mammals", fit_ok)
  pm_birds <- gomp_params |> dplyr::filter(group == "Birds", fit_ok)

  loess_models <- list(
    meta = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      K0 = k0,
      curves = curves,
      span = span,
      x_transform = "log10(mass_g)",
      y_transform = "log(parameter)",
      source_files = paths[c("mammals_points", "birds_points")]
    ),
    mammals = if (isTRUE(run_mammals)) fit_loess_set(pm_mammals, curves, span) else NULL,
    birds = NULL
  )

  if (isTRUE(run_birds)) {
    bird_by_combo <- split(pm_birds, pm_birds$combo)
    loess_models$birds <- stats::setNames(
      lapply(bird_by_combo, fit_loess_set, curves = curves, span = span),
      names(bird_by_combo)
    )
  }

  loess_pred <- dplyr::bind_rows(
    if (isTRUE(run_mammals)) {
      predict_loess_set(pm_mammals, loess_models$mammals, curves = curves, z = z) |>
        dplyr::mutate(group = "Mammals", combo = NA_character_)
    } else tibble::tibble(),
    if (isTRUE(run_birds) && length(loess_models$birds) > 0) {
      dplyr::bind_rows(lapply(names(loess_models$birds), function(cb) {
        predict_loess_set(dplyr::filter(pm_birds, combo == cb), loess_models$birds[[cb]], curves = curves, z = z) |>
          dplyr::mutate(group = "Birds", combo = cb)
      }))
    } else tibble::tibble()
  )

  list(models = loess_models, prediction_grid = loess_pred)
}
