# Summary statistics for 7.3_persist_cmp.Rmd.
# Requires the core 7.3 persistence comparison objects.

# ============================================================
# 18_results_text_statistics
#
# Main-text and Supplementary Methods statistics.
#
# Computes:
#   1. late-stage mean persistence loss reduction
#   2. species-level signed AUC
#   3. area-normalized species persistence advantage
#   4. interpolated threshold-crossing asymmetry
#   5. checkpoint tail-risk counts
#   6. PU redundancy-loss counts
# ============================================================

suppressPackageStartupMessages({
  library(data.table)
})

# -----------------------------
# Reload standard 7.3 outputs if needed
# -----------------------------

if (!exists("stage_sum")) {
  if (exists("out_stage") && file.exists(out_stage)) {
    stage_sum <- fread(out_stage)
  } else {
    stop("`stage_sum` not found and `out_stage` is unavailable.", call. = FALSE)
  }
}

if (!exists("sp_long")) {
  if (exists("out_sp") && file.exists(out_sp)) {
    sp_long <- fread(out_sp)
  } else {
    stop("`sp_long` not found and `out_sp` is unavailable.", call. = FALSE)
  }
}

if (!exists("pu_long")) {
  if (exists("out_pu") && file.exists(out_pu)) {
    pu_long <- fread(out_pu)
  } else {
    stop("`pu_long` not found and `out_pu` is unavailable.", call. = FALSE)
  }
}

if (!exists("rank_method")) rank_method <- "abf"

if (!exists("rank_label")) {
  rank_label <- benchmark_rank_label(rank_method)
}

metric_eps <- 1e-9

metric_median <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  median(x, na.rm = TRUE)
}

metric_quantile <- function(x, p) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  as.numeric(quantile(x, probs = p, na.rm = TRUE, names = FALSE))
}

trapz_xy <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]
  y <- y[ok]

  if (length(x) < 2L) return(NA_real_)

  ord <- order(x)
  x <- x[ord]
  y <- y[ord]

  dd <- data.table(x = x, y = y)
  dd <- dd[, .(y = mean(y, na.rm = TRUE)), by = x]
  setorder(dd, x)

  if (nrow(dd) < 2L) return(NA_real_)

  sum(diff(dd$x) * (head(dd$y, -1L) + tail(dd$y, -1L)) / 2)
}

x_to_pct <- function(x) {
  if (max(stage_sum$x, na.rm = TRUE) <= 1.0001) 100 * x else x
}

target_to_x <- function(pct) {
  if (max(stage_sum$x, na.rm = TRUE) <= 1.0001) pct / 100 else pct
}

sid_fallback <- function(scientific_name) {
  if (exists("sid", mode = "function")) {
    sid(scientific_name)
  } else {
    x <- trimws(as.character(scientific_name))
    x <- gsub("\\s+", "_", x)
    x <- gsub("[/\\\\:<>\"|?*]+", "_", x)
    x <- gsub("[^A-Za-z0-9_]+", "_", x)
    x <- gsub("_+", "_", x)
    gsub("^_|_$", "", x)
  }
}

resolve_species_id <- function(scientific_name, dt = sp_long) {
  candidate <- sid_fallback(scientific_name)

  if (candidate %in% unique(dt$species)) {
    return(candidate)
  }

  hit <- unique(dt[
    tolower(trimws(scientificName)) == tolower(trimws(scientific_name)),
    species
  ])

  if (length(hit)) return(hit[1L])

  NA_character_
}

focal_name_results <- if (
  exists("params") &&
    !is.null(params$focus_species) &&
    length(params$focus_species) >= 1L &&
    nzchar(trimws(as.character(params$focus_species[[1]])))
) {
  as.character(params$focus_species[[1]])
} else {
  "Pterocles personatus"
}

focal_id_results <- resolve_species_id(focal_name_results, sp_long)

# ============================================================
# 1. Late-stage mean persistence and loss reduction
# ============================================================

target_pcts <- c(90, 95, 99)

stage_key <- unique(stage_sum[, .(stage, stage_order, x)])
setorder(stage_key, stage_order)

checkpoint_key <- rbindlist(lapply(target_pcts, function(pct_i) {
  target_x <- target_to_x(pct_i)
  ii <- which.min(abs(stage_key$x - target_x))

  out <- copy(stage_key[ii])
  out[, target_pct_removed := pct_i]
  out[, actual_pct_removed := x_to_pct(x)]
  out
}), use.names = TRUE, fill = TRUE)

checkpoint_key <- checkpoint_key[!duplicated(stage)]

checkpoint_stage_sum <- merge(
  stage_sum,
  checkpoint_key[
    ,
    .(
      stage,
      stage_order,
      x,
      target_pct_removed,
      actual_pct_removed
    )
  ],
  by = c("stage", "stage_order", "x"),
  all.x = FALSE,
  sort = FALSE
)

checkpoint_stage_sum[
  ,
  persistence_loss := 1 - mean_persist
]

checkpoint_stage_wide <- dcast(
  checkpoint_stage_sum,
  target_pct_removed + actual_pct_removed + stage + stage_order + x ~ method,
  value.var = c(
    "mean_persist",
    "persistence_loss",
    "median_persist",
    "frac_gt_05",
    "total_pu_area_km2"
  )
)

checkpoint_stage_wide[
  ,
  `:=`(
    mean_persist_diff_pipe_minus_rank =
      mean_persist_pipe - mean_persist_rank,
    persistence_loss_reduction_abs =
      persistence_loss_rank - persistence_loss_pipe,
    persistence_loss_reduction_pct_vs_rank =
      fifelse(
        persistence_loss_rank > metric_eps,
        100 * (persistence_loss_rank - persistence_loss_pipe) /
          persistence_loss_rank,
        NA_real_
      )
  )
]

setorder(checkpoint_stage_wide, target_pct_removed)

# ============================================================
# 2. Species-level signed AUC
# ============================================================

auc_dt <- dcast(
  sp_long[
    ,
    .(
      species,
      scientificName,
      className,
      stage,
      stage_order,
      x,
      method,
      sp_persist
    )
  ],
  species + scientificName + className + stage + stage_order + x ~ method,
  value.var = "sp_persist"
)

auc_dt <- auc_dt[is.finite(pipe) & is.finite(rank)]
auc_dt[, gap_pipe_minus_rank := pipe - rank]
setorder(auc_dt, species, stage_order)

species_auc <- auc_dt[
  ,
  .(
    n_stages = .N,
    auc_diff_pipe_minus_rank =
      trapz_xy(x, gap_pipe_minus_rank),
    mean_gap_over_sequence =
      trapz_xy(x, gap_pipe_minus_rank) /
      (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)),
    max_gap_pipe_minus_rank =
      max(gap_pipe_minus_rank, na.rm = TRUE),
    min_gap_pipe_minus_rank =
      min(gap_pipe_minus_rank, na.rm = TRUE),
    x_at_max_gap =
      x[which.max(gap_pipe_minus_rank)]
  ),
  by = .(species, scientificName, className)
]

species_auc <- species_auc[is.finite(auc_diff_pipe_minus_rank)]
setorder(species_auc, -auc_diff_pipe_minus_rank)
species_auc[, positive_auc_rank := seq_len(.N)]

n_species_auc <- nrow(species_auc)
n_positive_auc <- species_auc[
  ,
  sum(auc_diff_pipe_minus_rank > metric_eps, na.rm = TRUE)
]
n_negative_auc <- species_auc[
  ,
  sum(auc_diff_pipe_minus_rank < -metric_eps, na.rm = TRUE)
]

focal_auc_row <- species_auc[species == focal_id_results]

# ============================================================
# 3. Area-normalized persistence advantage
# ============================================================

area_norm_one_species <- function(dd, n_grid = 101L) {
  sp_id <- dd$species[1L]
  sp_name <- dd$scientificName[1L]
  class_name <- dd$className[1L]

  prep_curve <- function(method_i) {
    z <- copy(dd[
      method == method_i &
        is.finite(total_pu_area_km2) &
        is.finite(sp_persist),
      .(
        area = as.numeric(total_pu_area_km2),
        persist = as.numeric(sp_persist)
      )
    ])

    # Same convention used in earlier area-normalized metrics:
    # anchor at origin for the end of the area-persistence curve.
    z <- rbindlist(
      list(z, data.table(area = 0, persist = 0)),
      use.names = TRUE,
      fill = TRUE
    )

    z <- z[
      is.finite(area) & is.finite(persist) & area >= 0,
      .(persist = mean(persist, na.rm = TRUE)),
      by = area
    ]

    setorder(z, area)
    z
  }

  pipe_curve <- prep_curve("pipe")
  rank_curve <- prep_curve("rank")

  if (nrow(pipe_curve) < 2L || nrow(rank_curve) < 2L) {
    return(data.table(
      species = sp_id,
      scientificName = sp_name,
      className = class_name,
      mean_gap_area_normalized = NA_real_,
      frac_area_grid_pipe_gt_rank = NA_real_
    ))
  }

  area_min_common <- max(
    min(pipe_curve$area, na.rm = TRUE),
    min(rank_curve$area, na.rm = TRUE)
  )

  area_max_common <- min(
    max(pipe_curve$area, na.rm = TRUE),
    max(rank_curve$area, na.rm = TRUE)
  )

  if (
    !is.finite(area_min_common) ||
      !is.finite(area_max_common) ||
      area_max_common <= area_min_common
  ) {
    return(data.table(
      species = sp_id,
      scientificName = sp_name,
      className = class_name,
      mean_gap_area_normalized = NA_real_,
      frac_area_grid_pipe_gt_rank = NA_real_
    ))
  }

  area_grid <- seq(area_min_common, area_max_common, length.out = n_grid)

  pipe_interp <- approx(
    x = pipe_curve$area,
    y = pipe_curve$persist,
    xout = area_grid,
    rule = 1,
    ties = "ordered"
  )$y

  rank_interp <- approx(
    x = rank_curve$area,
    y = rank_curve$persist,
    xout = area_grid,
    rule = 1,
    ties = "ordered"
  )$y

  gap <- pipe_interp - rank_interp
  area_width <- area_max_common - area_min_common

  data.table(
    species = sp_id,
    scientificName = sp_name,
    className = class_name,
    mean_gap_area_normalized =
      trapz_xy(area_grid, gap) / area_width,
    median_gap_area_normalized =
      metric_median(gap),
    max_gap_area_normalized =
      max(gap, na.rm = TRUE),
    min_gap_area_normalized =
      min(gap, na.rm = TRUE),
    frac_area_grid_pipe_gt_rank =
      mean(gap > metric_eps, na.rm = TRUE)
  )
}

area_norm_dt <- rbindlist(
  lapply(
    split(
      sp_long[
        ,
        .(
          method,
          species,
          scientificName,
          className,
          stage,
          stage_order,
          x,
          total_pu_area_km2,
          sp_persist
        )
      ],
      by = "species",
      keep.by = TRUE
    ),
    area_norm_one_species
  ),
  use.names = TRUE,
  fill = TRUE
)

area_norm_dt <- area_norm_dt[is.finite(mean_gap_area_normalized)]

n_area_species <- nrow(area_norm_dt)
n_area_positive <- area_norm_dt[
  ,
  sum(mean_gap_area_normalized > metric_eps, na.rm = TRUE)
]
n_area_negative <- area_norm_dt[
  ,
  sum(mean_gap_area_normalized < -metric_eps, na.rm = TRUE)
]

# ============================================================
# 4. Continuous/interpolated threshold crossing
# ============================================================

interp_first_crossing <- function(x, y, stage_order, threshold) {
  ok <- is.finite(x) & is.finite(y) & is.finite(stage_order)
  x <- x[ok]
  y <- y[ok]
  stage_order <- stage_order[ok]

  if (!length(x)) {
    return(data.table(
      crossed = FALSE,
      x_cross_continuous = NA_real_,
      x_cross_or_end = NA_real_
    ))
  }

  ord <- order(stage_order)
  x <- x[ord]
  y <- y[ord]

  if (y[1L] < threshold) {
    return(data.table(
      crossed = TRUE,
      x_cross_continuous = x[1L],
      x_cross_or_end = x[1L]
    ))
  }

  trans <- which(head(y, -1L) >= threshold & tail(y, -1L) < threshold)

  if (!length(trans)) {
    return(data.table(
      crossed = FALSE,
      x_cross_continuous = NA_real_,
      x_cross_or_end = max(x, na.rm = TRUE)
    ))
  }

  ii <- trans[1L]
  x1 <- x[ii]
  x2 <- x[ii + 1L]
  y1 <- y[ii]
  y2 <- y[ii + 1L]

  x_cross <- if (!is.finite(y2 - y1) || abs(y2 - y1) < .Machine$double.eps) {
    x2
  } else {
    x1 + (threshold - y1) * (x2 - x1) / (y2 - y1)
  }

  data.table(
    crossed = TRUE,
    x_cross_continuous = x_cross,
    x_cross_or_end = x_cross
  )
}

threshold_values <- c(0.9, 0.5, 0.1)

threshold_long <- rbindlist(lapply(threshold_values, function(thr) {
  sp_long[
    ,
    {
      out <- interp_first_crossing(
        x = x,
        y = sp_persist,
        stage_order = stage_order,
        threshold = thr
      )
      out[, threshold := thr]
      out
    },
    by = .(method, species, scientificName, className)
  ]
}), use.names = TRUE, fill = TRUE)

threshold_wide <- dcast(
  threshold_long,
  threshold + species + scientificName + className ~ method,
  value.var = c(
    "crossed",
    "x_cross_continuous",
    "x_cross_or_end"
  )
)

threshold_wide[
  ,
  crossing_state := fcase(
    crossed_pipe & crossed_rank, "both_crossed",
    !crossed_pipe & crossed_rank, "rank_only",
    crossed_pipe & !crossed_rank, "pipe_only",
    !crossed_pipe & !crossed_rank, "neither_crossed",
    default = NA_character_
  )
]

threshold_wide[
  ,
  delay_x_pipe_minus_rank_both_crossed :=
    fifelse(
      crossed_pipe & crossed_rank,
      x_cross_continuous_pipe - x_cross_continuous_rank,
      NA_real_
    )
]

threshold_summary <- threshold_wide[
  ,
  .(
    n_species = .N,
    n_both_crossed = sum(crossing_state == "both_crossed", na.rm = TRUE),
    n_rank_only = sum(crossing_state == "rank_only", na.rm = TRUE),
    n_pipe_only = sum(crossing_state == "pipe_only", na.rm = TRUE),
    n_neither_crossed = sum(crossing_state == "neither_crossed", na.rm = TRUE),
    mean_delay_both_crossed =
      mean(delay_x_pipe_minus_rank_both_crossed, na.rm = TRUE),
    median_delay_both_crossed =
      metric_median(delay_x_pipe_minus_rank_both_crossed)
  ),
  by = threshold
]

setorder(threshold_summary, -threshold)

# ============================================================
# 5. Checkpoint tail-risk counts
# ============================================================

tail_targets <- c(95, 99)
tail_thresholds <- c(0.9, 0.5, 0.1)

tail_checkpoint_key <- rbindlist(lapply(tail_targets, function(pct_i) {
  target_x <- target_to_x(pct_i)
  ii <- which.min(abs(stage_key$x - target_x))
  out <- copy(stage_key[ii])
  out[, target_pct_removed := pct_i]
  out[, actual_pct_removed := x_to_pct(x)]
  out
}), use.names = TRUE, fill = TRUE)

tail_checkpoint_key <- tail_checkpoint_key[!duplicated(stage)]

tail_sp <- merge(
  sp_long[stage %in% tail_checkpoint_key$stage],
  tail_checkpoint_key[
    ,
    .(
      stage,
      stage_order,
      x,
      target_pct_removed,
      actual_pct_removed
    )
  ],
  by = c("stage", "stage_order", "x"),
  all.x = FALSE,
  sort = FALSE
)

tail_pair <- dcast(
  tail_sp[
    ,
    .(
      target_pct_removed,
      actual_pct_removed,
      stage,
      stage_order,
      x,
      species,
      scientificName,
      className,
      method,
      sp_persist
    )
  ],
  target_pct_removed + actual_pct_removed + stage + stage_order + x +
    species + scientificName + className ~ method,
  value.var = "sp_persist"
)

tail_summary <- rbindlist(lapply(tail_thresholds, function(thr) {
  tail_pair[
    ,
    .(
      threshold = thr,
      n_species = .N,
      n_pipe_kept_above_rank_below =
        sum(pipe >= thr & rank < thr, na.rm = TRUE),
      n_rank_kept_above_pipe_below =
        sum(rank >= thr & pipe < thr, na.rm = TRUE),
      n_both_above =
        sum(pipe >= thr & rank >= thr, na.rm = TRUE),
      n_both_below =
        sum(pipe < thr & rank < thr, na.rm = TRUE)
    ),
    by = .(
      target_pct_removed,
      actual_pct_removed,
      stage,
      stage_order,
      x
    )
  ]
}), use.names = TRUE, fill = TRUE)

setorder(tail_summary, target_pct_removed, -threshold)

# ============================================================
# 6. PU redundancy-loss timing
# ============================================================

method_stage_grid <- unique(sp_long[, .(method, stage, stage_order, x)])
species_grid <- unique(sp_long[, .(species, scientificName, className)])

method_stage_grid[, tmp_join := 1L]
species_grid[, tmp_join := 1L]

full_species_stage_grid <- merge(
  method_stage_grid,
  species_grid,
  by = "tmp_join",
  allow.cartesian = TRUE,
  sort = FALSE
)[, tmp_join := NULL]

pu_redundancy_obs <- pu_long[
  ,
  .(
    n_pu_current = .N,
    n_pu_gt_09 = sum(P_pu > 0.9, na.rm = TRUE),
    n_pu_gt_05 = sum(P_pu > 0.5, na.rm = TRUE),
    expected_persisting_pu = sum(P_pu, na.rm = TRUE),
    max_P_pu = max(P_pu, na.rm = TRUE)
  ),
  by = .(
    method,
    stage,
    stage_order,
    x,
    species,
    scientificName
  )
]

pu_redundancy_full <- merge(
  full_species_stage_grid,
  pu_redundancy_obs,
  by = c(
    "method",
    "stage",
    "stage_order",
    "x",
    "species",
    "scientificName"
  ),
  all.x = TRUE,
  sort = FALSE
)

for (cc in c(
  "n_pu_current",
  "n_pu_gt_09",
  "n_pu_gt_05",
  "expected_persisting_pu",
  "max_P_pu"
)) {
  pu_redundancy_full[is.na(get(cc)), (cc) := 0]
}

first_condition_time <- function(x, y, stage_order, condition_fun) {
  ok <- is.finite(x) & is.finite(y) & is.finite(stage_order)
  x <- x[ok]
  y <- y[ok]
  stage_order <- stage_order[ok]

  if (!length(x)) {
    return(data.table(
      condition_reached = FALSE,
      x_condition = NA_real_,
      x_condition_or_end = NA_real_
    ))
  }

  ord <- order(stage_order)
  x <- x[ord]
  y <- y[ord]

  cond <- condition_fun(y)
  reached <- any(cond, na.rm = TRUE)

  if (reached) {
    ii <- which(cond)[1L]
    x_condition <- x[ii]
  } else {
    x_condition <- NA_real_
  }

  data.table(
    condition_reached = reached,
    x_condition = x_condition,
    x_condition_or_end = if (reached) x_condition else max(x, na.rm = TRUE)
  )
}

redundancy_specs <- data.table(
  target_id = c(
    "no_PU_gt_0.9",
    "no_PU_gt_0.5",
    "expected_persisting_PU_lt_1"
  ),
  metric = c(
    "n_pu_gt_09",
    "n_pu_gt_05",
    "expected_persisting_pu"
  ),
  threshold = c(0, 0, 1),
  condition_type = c("le", "le", "lt"),
  manuscript_label = c(
    "lost all PUs with P_PU > 0.9",
    "lost all PUs with P_PU > 0.5",
    "fell below one expected persisting PU"
  )
)

redundancy_long <- rbindlist(lapply(seq_len(nrow(redundancy_specs)), function(ii) {
  metric_i <- redundancy_specs$metric[ii]
  threshold_i <- redundancy_specs$threshold[ii]
  condition_type_i <- redundancy_specs$condition_type[ii]

  condition_fun <- switch(
    condition_type_i,
    le = function(z) z <= threshold_i,
    lt = function(z) z < threshold_i,
    stop("Unknown condition type.")
  )

  out <- pu_redundancy_full[
    ,
    first_condition_time(
      x = x,
      y = get(metric_i),
      stage_order = stage_order,
      condition_fun = condition_fun
    ),
    by = .(method, species, scientificName, className)
  ]

  out[, target_id := redundancy_specs$target_id[ii]]
  out[, metric := metric_i]
  out[, threshold := threshold_i]
  out[, condition_type := condition_type_i]
  out[, manuscript_label := redundancy_specs$manuscript_label[ii]]
  out
}), use.names = TRUE, fill = TRUE)

redundancy_wide <- dcast(
  redundancy_long,
  target_id + metric + threshold + condition_type + manuscript_label +
    species + scientificName + className ~ method,
  value.var = c(
    "condition_reached",
    "x_condition",
    "x_condition_or_end"
  )
)

redundancy_wide[
  ,
  loss_state := fcase(
    condition_reached_pipe & condition_reached_rank, "both_lost",
    !condition_reached_pipe & condition_reached_rank, "rank_only_lost",
    condition_reached_pipe & !condition_reached_rank, "pipe_only_lost",
    !condition_reached_pipe & !condition_reached_rank, "neither_lost",
    default = NA_character_
  )
]

redundancy_summary <- redundancy_wide[
  ,
  .(
    n_species = .N,
    n_both_lost = sum(loss_state == "both_lost", na.rm = TRUE),
    n_rank_only_lost = sum(loss_state == "rank_only_lost", na.rm = TRUE),
    n_pipe_only_lost = sum(loss_state == "pipe_only_lost", na.rm = TRUE),
    n_neither_lost = sum(loss_state == "neither_lost", na.rm = TRUE),
    mean_censored_delay =
      mean(x_condition_or_end_pipe - x_condition_or_end_rank, na.rm = TRUE),
    median_censored_delay =
      metric_median(x_condition_or_end_pipe - x_condition_or_end_rank)
  ),
  by = .(target_id, manuscript_label)
]


