# JAGS model fitting helpers for log10-scale calibration regressions.

jags_lm_model <- "
model {
  for (i in 1:N) {
    y[i] ~ dnorm(mu[i], tau)
    mu[i] <- alpha + inprod(beta[1:P], X[i,1:P])
  }
  alpha ~ dnorm(0.0, tau_beta)
  for (j in 1:P) { beta[j] ~ dnorm(0.0, tau_beta) }
  sigma ~ dunif(sigma_min, sigma_max)
  tau <- pow(sigma, -2)
}
"

assert_jags_available <- function(run_fits) {
  if (isTRUE(run_fits) &&
      !tryCatch(length(rjags::jags.version()) > 0, error = function(e) FALSE)) {
    stop(
      "run_fits is TRUE, but rjags cannot find a working JAGS installation. ",
      "Install JAGS or set params$run_fits to false and provide posterior CSVs.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

fit_jags_lm <- function(y, X, coef_names, out_path, bayes) {
  X <- as.matrix(X)
  data_jags <- list(
    N = length(y),
    P = ncol(X),
    y = as.numeric(y),
    X = X,
    tau_beta = bayes$tau_beta,
    sigma_min = bayes$sigma_min,
    sigma_max = bayes$sigma_max
  )

  inits <- lapply(seq_len(bayes$n_chains), function(ch) {
    list(
      alpha = rnorm(1, 0, 1),
      beta = rnorm(data_jags$P, 0, 1),
      sigma = runif(1, 0.05, 0.5),
      .RNG.name = "base::Wichmann-Hill",
      .RNG.seed = bayes$seed + ch
    )
  })

  jm <- rjags::jags.model(
    textConnection(jags_lm_model),
    data = data_jags,
    inits = inits,
    n.chains = bayes$n_chains,
    n.adapt = bayes$n_adapt
  )

  beta_mon <- paste0("beta[", seq_len(data_jags$P), "]")
  mcmc <- rjags::coda.samples(
    jm,
    variable.names = c("alpha", beta_mon, "sigma"),
    n.iter = bayes$n_iter,
    thin = bayes$thin
  )

  draws <- tibble::as_tibble(as.data.frame(as.matrix(mcmc), check.names = FALSE))
  out <- draws[, c("alpha", beta_mon, "sigma")]
  names(out) <- c("alpha", paste0("beta_", coef_names), "sigma")

  readr::write_csv(out, out_path)
  out
}

fit_loglog_allometry <- function(mass_g, y_pos, out_path, bayes) {
  fit_jags_lm(
    y = log10(y_pos),
    X = matrix(log10(mass_g), ncol = 1),
    coef_names = "logM",
    out_path = out_path,
    bayes = bayes
  )
}

fit_all_demographic_models <- function(cal, paths, bayes) {
  fit_loglog_allometry(cal$mammal_rmax$Mass_g, cal$mammal_rmax$rm,
                       paths$clean$post_mammal_rm, bayes)
  fit_loglog_allometry(cal$mammal_sigma$Mass_g, cal$mammal_sigma$sigma,
                       paths$clean$post_mammal_sigma, bayes)
  fit_loglog_allometry(cal$bird_rmax$Mass_g, cal$bird_rmax$rm,
                       paths$clean$post_bird_rm, bayes)

  diet_parts <- c("Diet_Inv", "Diet_AllPlants", "Diet_VertFishScav")
  ilr <- to_ilr3(cal$bird_sigma_data[, diet_parts] / 100)
  X_bird_sigma <- cbind(log10(cal$bird_sigma_data$Mass_g), ilr$ilr1, ilr$ilr2)

  fit_jags_lm(
    y = log10(cal$bird_sigma_data$sigma),
    X = X_bird_sigma,
    coef_names = c("logM", "ilr1", "ilr2"),
    out_path = paths$clean$post_bird_sigma_ilr3,
    bayes = bayes
  )
}
