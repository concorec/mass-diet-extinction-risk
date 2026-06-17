# Compositional-data helpers for avian diet predictors.

zero_replace <- function(prop_df) {
  zCompositions::cmultRepl(as.data.frame(prop_df), method = "CZM", label = 0) |>
    as.data.frame()
}

percent_to_prop <- function(d, parts) {
  X <- as.data.frame(d[, parts, drop = FALSE]) / 100
  X / rowSums(X)
}

to_ilr3 <- function(prop_df) {
  X <- as.data.frame(prop_df)
  X <- X / rowSums(X)
  X <- zero_replace(X)
  X <- X / rowSums(X)

  Z <- compositions::ilr(compositions::acomp(X)) |>
    as.data.frame()
  names(Z) <- c("ilr1", "ilr2")
  Z
}

add_ilr_from_percent <- function(d, parts, prefix = "ilr") {
  X <- percent_to_prop(d, parts)
  Z <- compositions::ilr(compositions::acomp(zero_replace(X))) |>
    as.data.frame()
  names(Z) <- paste0(prefix, seq_len(ncol(Z)))
  dplyr::bind_cols(d, Z)
}

ternary_diet_grid_10pct <- function() {
  vals <- seq(0, 100, by = 10)
  grid <- expand.grid(
    Diet_Inv = vals,
    Diet_AllPlants = vals,
    KEEP.OUT.ATTRS = FALSE
  )
  grid$Diet_VertFishScav <- 100 - grid$Diet_Inv - grid$Diet_AllPlants
  grid <- grid[grid$Diet_VertFishScav >= 0, , drop = FALSE]

  grid <- tibble::as_tibble(grid) |>
    dplyr::arrange(Diet_Inv, Diet_AllPlants, Diet_VertFishScav) |>
    dplyr::mutate(
      combo = paste(Diet_Inv, Diet_AllPlants, Diet_VertFishScav, sep = "_")
    ) |>
    dplyr::select(combo, Diet_Inv, Diet_AllPlants, Diet_VertFishScav)

  if (nrow(grid) != 66L || dplyr::n_distinct(grid$combo) != 66L) {
    stop("The 10% ternary bird diet grid must contain exactly 66 unique combinations.", call. = FALSE)
  }

  grid
}

bird_ilr_combo_design <- function() {
  grid <- ternary_diet_grid_10pct()
  ilr <- to_ilr3(grid[, c("Diet_Inv", "Diet_AllPlants", "Diet_VertFishScav")] / 100)
  out <- dplyr::bind_cols(grid["combo"], ilr)

  if (any(!is.finite(as.matrix(out[, c("ilr1", "ilr2")])))) {
    stop("Generated bird ILR combo coordinates contain non-finite values.", call. = FALSE)
  }

  out
}
