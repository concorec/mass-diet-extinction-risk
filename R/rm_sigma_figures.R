# Figure helpers for 1_rm_sigma_models.Rmd.

mass_grid_from_data <- function(mass_g, n = 200) {
  tibble::tibble(Mass_g = log_space(min(mass_g, na.rm = TRUE), max(mass_g, na.rm = TRUE), n))
}

simplex_grid_10pct_prop <- function() {
  ternary_diet_grid_10pct() |>
    dplyr::mutate(
      Diet_Inv = Diet_Inv / 100,
      Diet_AllPlants = Diet_AllPlants / 100,
      Diet_VertFishScav = Diet_VertFishScav / 100
    )
}

simplex_grid_10pct <- function() {
  k <- 10
  grid <- expand.grid(i = 0:k, j = 0:k)
  grid$k <- k - grid$i - grid$j
  grid <- grid[grid$k >= 0, , drop = FALSE]
  tibble::tibble(
    Diet_Inv = grid$i / k,
    Diet_AllPlants = grid$j / k,
    Diet_VertFishScav = grid$k / k
  )
}

bird_diet_palette <- function(diet_values) {
  diet_order <- c("FruiNect", "Invertebrate", "Omnivore", "PlantSeed", "VertFishScav", "Unknown")
  diet_levels <- c(
    intersect(diet_order, unique(as.character(diet_values))),
    setdiff(sort(unique(as.character(diet_values))), diet_order)
  )
  diet_labels <- c(
    FruiNect = "Fruit + nectar",
    Invertebrate = "Invertebrate",
    Omnivore = "Omnivore",
    PlantSeed = "Plant + seed",
    VertFishScav = "Vertebrate + fish + scavenging",
    Unknown = "Unknown"
  )[diet_levels]
  diet_pal <- stats::setNames(
    viridisLite::viridis(length(diet_levels), option = "D", end = 0.9),
    diet_levels
  )
  if ("Unknown" %in% diet_levels) diet_pal["Unknown"] <- "grey70"

  list(levels = diet_levels, labels = diet_labels, palette = diet_pal)
}

band_from_posterior <- function(mass_grid, draws_path) {
  draws <- readr::read_csv(
    draws_path,
    show_col_types = FALSE,
    col_select = dplyr::all_of(c("alpha", "beta_logM"))
  )

  logM <- log10(mass_grid)
  pred_log10 <- draws$alpha + tcrossprod(draws$beta_logM, logM)
  pred <- 10^pred_log10

  qs <- apply(pred, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
  tibble::tibble(Mass = mass_grid, lo = qs[1, ], mid = qs[2, ], hi = qs[3, ])
}

band_from_ols <- function(df, mass_grid) {
  fit <- stats::lm(log10(y) ~ log10(Mass), data = df)
  pr <- stats::predict(
    fit,
    newdata = data.frame(Mass = mass_grid),
    interval = "confidence",
    level = 0.95
  )

  tibble::tibble(
    Mass = mass_grid,
    lo = 10^pr[, "lwr"],
    mid = 10^pr[, "fit"],
    hi = 10^pr[, "upr"]
  )
}

make_allometry_panel <- function(df, y_lab, x_limits, x_breaks, draws_path = NULL) {
  df <- df[df$Mass > 0 & df$y > 0 & is.finite(df$Mass) & is.finite(df$y), , drop = FALSE]

  mass_grid <- log_space(x_limits[1], x_limits[2], n = 300)
  band <- if (is.null(draws_path)) band_from_ols(df, mass_grid) else band_from_posterior(mass_grid, draws_path)

  r2 <- summary(stats::lm(log10(y) ~ log10(Mass), data = df))$r.squared
  ann <- sprintf("N = %d\nR^2 = %.2f", nrow(df), r2)

  ggplot2::ggplot(df, ggplot2::aes(Mass, y)) +
    ggplot2::geom_ribbon(
      data = band,
      ggplot2::aes(x = Mass, ymin = lo, ymax = hi),
      inherit.aes = FALSE,
      fill = "grey60",
      alpha = 0.25
    ) +
    ggplot2::geom_line(
      data = band,
      ggplot2::aes(x = Mass, y = mid),
      inherit.aes = FALSE,
      linewidth = 0.75,
      color = "grey10"
    ) +
    ggplot2::geom_point(size = 1.2, alpha = 0.65, color = "grey10") +
    ggplot2::annotate(
      "label",
      x = Inf, y = Inf,
      label = ann,
      hjust = 1.02,
      vjust = 1.05,
      label.size = 0,
      fill = scales::alpha("white", 0.85),
      size = 3.6,
      fontface = "bold",
      color = "grey10"
    ) +
    ggplot2::scale_x_log10(
      limits = x_limits,
      breaks = x_breaks,
      labels = scales::label_number(big.mark = ",")
    ) +
    ggplot2::scale_y_log10(labels = scales::label_number(big.mark = ",")) +
    ggplot2::labs(x = "Body mass (g)", y = y_lab) +
    theme_pub(base_size = 12) +
    ggplot2::theme(aspect.ratio = 0.72)
}

