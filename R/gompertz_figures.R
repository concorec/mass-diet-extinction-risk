# Figure helpers for 3_gompertz_loess_models.Rmd.

fmt_mass_label <- function(mass_g) {
  if (!is.finite(mass_g)) return("Mass: NA")
  if (mass_g < 1000) return(paste0("Mass: ", scales::label_number(accuracy = 1)(mass_g), " g"))
  kg <- mass_g / 1000
  acc <- if (kg >= 100) 1 else 0.1
  paste0("Mass: ", scales::label_number(accuracy = acc, big.mark = ",")(kg), " kg")
}

fmt_mass_axis <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_character_, length(x))
  ok <- is.finite(x)
  if (!any(ok)) return(out)
  x_ok <- x[ok]
  is_g <- x_ok < 1000
  out_ok <- character(length(x_ok))
  if (any(is_g)) out_ok[is_g] <- paste0(scales::label_number(accuracy = 1, big.mark = ",")(x_ok[is_g]), " g")
  if (any(!is_g)) {
    kg <- x_ok[!is_g] / 1000
    acc <- ifelse(kg >= 100, 1, 0.1)
    out_ok[!is_g] <- paste0(scales::label_number(accuracy = acc, big.mark = ",")(kg), " kg")
  }
  out[ok] <- out_ok
  out
}

mammal_K_range <- function(curves_tbl, curve_use) {
  rng <- curves_tbl |>
    dplyr::filter(group == "Mammals", curve == curve_use) |>
    dplyr::summarise(K_min = min(K, na.rm = TRUE), K_max = max(K, na.rm = TRUE))
  list(K_min = as.numeric(rng$K_min), K_max = as.numeric(rng$K_max))
}

safe_K_left <- function(K_left, K0) {
  K_left_safe <- max(as.numeric(K_left), 1e-6)
  if (is.finite(K0) && K_left_safe <= K0) K_left_safe <- max(K_left_safe, K0 * 1.0005)
  K_left_safe
}

build_fitted_curves <- function(pm, K_grid, K0) {
  pm |>
    dplyr::mutate(mass_g = as.numeric(mass_g), logM = log10(as.numeric(mass_g))) |>
    dplyr::arrange(logM) |>
    dplyr::group_by(mass_idx, mass_g, logM) |>
    dplyr::group_modify(~{
      tibble::tibble(K = K_grid, p = predict_p(.x$alpha[1], .x$beta[1], K_grid, K0 = K0))
    }) |>
    dplyr::ungroup()
}

mass_legend_breaks <- function(pm, idx = c(1L, 16L, 31L), labels = c("2 g", "3.1 kg", "4,750 kg")) {
  logM_all <- if ("logM" %in% names(pm)) as.numeric(pm$logM) else log10(as.numeric(pm$mass_g))
  g <- pm |> dplyr::filter(mass_idx %in% idx) |> dplyr::arrange(mass_idx) |> dplyr::pull(mass_g)
  list(limits = range(logM_all, na.rm = TRUE), breaks = log10(as.numeric(g)), labels = labels)
}

tick_K <- function(K_max, K_left_safe) {
  breaksK <- c(500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000)
  breaksK <- breaksK[breaksK >= K_left_safe & breaksK <= K_max * 1.001]
  if (length(breaksK) < 3) {
    breaksK <- scales::breaks_log(n = 4)(c(K_left_safe, K_max))
    breaksK <- breaksK[is.finite(breaksK) & breaksK >= K_left_safe & breaksK <= K_max * 1.001]
  }
  breaksK
}

r2_key_layout <- function(K_left, K_max,
                          x_box = c(0.05, 0.48), x_seg = c(0.08, 0.17), x_num = 0.20,
                          y_box = c(0.79, 0.965), y_title = 0.945, y_vals = c(0.890, 0.835)) {
  logL <- log10(K_left)
  logR <- log10(K_max)
  rng <- logR - logL
  list(
    box_x0 = 10^(logL + x_box[1] * rng), box_x1 = 10^(logL + x_box[2] * rng),
    seg_x0 = 10^(logL + x_seg[1] * rng), seg_x1 = 10^(logL + x_seg[2] * rng),
    num_x = 10^(logL + x_num * rng),
    box_y0 = y_box[1], box_y1 = y_box[2], y_title = y_title, y_vals = y_vals
  )
}

p_wolff <- function(K, c, A) exp(-A * exp(-c * K))

r2_prob <- function(y, yhat) {
  rss <- sum((y - yhat)^2, na.rm = TRUE)
  tss <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  if (is.finite(tss) && tss > 0) 1 - rss / tss else NA_real_
}

fit_wolff_c_start <- function(K, p, A) {
  p <- clamp01(p)
  y <- log(-log(p)) - log(A)
  fit <- stats::lm(y ~ 0 + K)
  c0 <- -as.numeric(stats::coef(fit)[["K"]])
  if (!is.finite(c0) || c0 <= 0) c0 <- 1e-8
  c0
}

fit_wolff_c_nls <- function(K, p, A = 10) {
  p <- clamp01(p)
  c0 <- fit_wolff_c_start(K, p, A = A)
  c_grid <- c0 * c(0.25, 0.5, 1, 2, 4)

  attempt <- function(c_start, iters = 250) {
    suppressWarnings(
      try(
        stats::nls(
          p ~ exp(-A * exp(-c * K)),
          start = list(c = c_start), algorithm = "port", lower = c(c = 1e-12),
          control = stats::nls.control(warnOnly = TRUE, maxiter = iters)
        ),
        silent = TRUE
      )
    )
  }

  for (cc in c_grid) {
    fit <- attempt(cc)
    if (!inherits(fit, "try-error") && !is.null(fit$convInfo) && isTRUE(fit$convInfo$isConv)) {
      return(as.numeric(stats::coef(fit)[["c"]]))
    }
  }

  old_seed <- if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
    get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  } else NULL
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    } else {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    }
  }, add = TRUE)

  set.seed(123)
  for (j in 1:20) {
    cc <- c0 * exp(stats::runif(1, log(0.2), log(5)))
    fit <- attempt(cc, iters = 350)
    if (!inherits(fit, "try-error") && !is.null(fit$convInfo) && isTRUE(fit$convInfo$isConv)) {
      return(as.numeric(stats::coef(fit)[["c"]]))
    }
  }

  stop("Wolff NLS fit failed for this panel.")
}

