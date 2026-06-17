# Console report for 7.3_persist_cmp.Rmd.
# Core trajectory construction, summary statistics, and figure construction are
# sourced before this file runs.

suppressPackageStartupMessages({
  library(data.table)
})

fmt_num <- function(x, digits = 3) {
  x <- as.numeric(x)[1L]
  if (!is.finite(x)) return("NA")
  formatC(x, format = "f", digits = digits, big.mark = ",")
}

fmt_int <- function(x) {
  x <- as.numeric(x)[1L]
  if (!is.finite(x)) return("NA")
  formatC(as.integer(round(x)), format = "d", big.mark = ",")
}

fmt_pct <- function(x, digits = 1) {
  x <- as.numeric(x)[1L]
  if (!is.finite(x)) return("NA")
  paste0(formatC(x, format = "f", digits = digits, big.mark = ","), "%")
}

fmt_x_pct <- function(x, digits = 1) {
  fmt_pct(x_to_pct(as.numeric(x)[1L]), digits = digits)
}

fmt_x_delta <- function(x, digits = 1) {
  x <- as.numeric(x)[1L]
  if (!is.finite(x)) return("NA")
  if (max(stage_sum$x, na.rm = TRUE) <= 1.0001) x <- 100 * x
  paste0(formatC(x, format = "f", digits = digits, big.mark = ","), " percentage points")
}

say <- function(...) {
  cat(..., "\n", sep = "")
}

section <- function(title) {
  cat("\n", title, "\n", sep = "")
  cat(paste(rep("-", nchar(title)), collapse = ""), "\n", sep = "")
}

safe_top <- function(dt, n = 3L) {
  if (!nrow(dt)) return(dt[0])
  dt[seq_len(min(n, nrow(dt)))]
}

print_auc_examples <- function(title, dt, n = 3L) {
  say(title)
  ex <- safe_top(dt, n)
  if (!nrow(ex)) {
    say("  None.")
    return(invisible(NULL))
  }

  for (ii in seq_len(nrow(ex))) {
    row <- ex[ii]
    say(
      "  ", ii, ". ", row$scientificName, " [", row$className, "]",
      ": signed AUC = ", fmt_num(row$auc_diff_pipe_minus_rank),
      "; mean gap = ", fmt_num(row$mean_gap_over_sequence),
      "; max gap = ", fmt_num(row$max_gap_pipe_minus_rank),
      " at ", fmt_x_pct(row$x_at_max_gap), " removed"
    )
  }

  invisible(NULL)
}

print_area_examples <- function(title, dt, n = 3L) {
  say(title)
  ex <- safe_top(dt, n)
  if (!nrow(ex)) {
    say("  None.")
    return(invisible(NULL))
  }

  for (ii in seq_len(nrow(ex))) {
    row <- ex[ii]
    say(
      "  ", ii, ". ", row$scientificName, " [", row$className, "]",
      ": mean area-normalized gap = ", fmt_num(row$mean_gap_area_normalized),
      "; median gap = ", fmt_num(row$median_gap_area_normalized),
      "; max gap = ", fmt_num(row$max_gap_area_normalized)
    )
  }

  invisible(NULL)
}

print_cross_examples <- function(threshold, state, title, n = 3L) {
  say(title)
  threshold_i <- as.numeric(threshold)[1L]
  ex <- threshold_wide[
    abs(threshold - threshold_i) < metric_eps &
      crossing_state == state
  ]

  if (state == "rank_only") {
    setorder(ex, x_cross_continuous_rank)
  } else if (state == "pipe_only") {
    setorder(ex, x_cross_continuous_pipe)
  }

  ex <- safe_top(ex, n)
  if (!nrow(ex)) {
    say("  None.")
    return(invisible(NULL))
  }

  for (ii in seq_len(nrow(ex))) {
    row <- ex[ii]
    if (state == "rank_only") {
      say(
        "  ", ii, ". ", row$scientificName, " [", row$className, "]",
        ": benchmark crossed at ", fmt_x_pct(row$x_cross_continuous_rank),
        "; persistence-based run did not cross"
      )
    } else {
      say(
        "  ", ii, ". ", row$scientificName, " [", row$className, "]",
        ": persistence-based run crossed at ", fmt_x_pct(row$x_cross_continuous_pipe),
        "; benchmark did not cross"
      )
    }
  }

  invisible(NULL)
}

print_tail_examples <- function(target_pct, threshold, direction, title, n = 3L) {
  say(title)
  ex <- tail_pair[target_pct_removed == target_pct]

  if (direction == "pipe_above_rank_below") {
    ex <- ex[pipe >= threshold & rank < threshold]
    ex[, gap := pipe - rank]
    setorder(ex, -gap)
  } else {
    ex <- ex[rank >= threshold & pipe < threshold]
    ex[, gap := rank - pipe]
    setorder(ex, -gap)
  }

  ex <- safe_top(ex, n)
  if (!nrow(ex)) {
    say("  None.")
    return(invisible(NULL))
  }

  for (ii in seq_len(nrow(ex))) {
    row <- ex[ii]
    gap_value <- if (direction == "pipe_above_rank_below") {
      row$pipe - row$rank
    } else {
      row$rank - row$pipe
    }
    gap_label <- if (direction == "pipe_above_rank_below") {
      "persistence-minus-benchmark gap"
    } else {
      "benchmark-minus-persistence gap"
    }
    say(
      "  ", ii, ". ", row$scientificName, " [", row$className, "]",
      ": persistence = ", fmt_num(row$pipe),
      "; benchmark = ", fmt_num(row$rank),
      "; ", gap_label, " = ", fmt_num(gap_value),
      " at ", fmt_pct(row$actual_pct_removed), " removal"
    )
  }

  invisible(NULL)
}

print_redundancy_examples <- function(target_id_i, state, title, n = 3L) {
  say(title)
  ex <- redundancy_wide[target_id == target_id_i & loss_state == state]

  if (state == "rank_only_lost") {
    setorder(ex, x_condition_rank)
  } else if (state == "pipe_only_lost") {
    setorder(ex, x_condition_pipe)
  }

  ex <- safe_top(ex, n)
  if (!nrow(ex)) {
    say("  None.")
    return(invisible(NULL))
  }

  for (ii in seq_len(nrow(ex))) {
    row <- ex[ii]
    if (state == "rank_only_lost") {
      say(
        "  ", ii, ". ", row$scientificName, " [", row$className, "]",
        ": benchmark reached the loss condition at ",
        fmt_x_pct(row$x_condition_rank),
        "; persistence-based run did not"
      )
    } else {
      say(
        "  ", ii, ". ", row$scientificName, " [", row$className, "]",
        ": persistence-based run reached the loss condition at ",
        fmt_x_pct(row$x_condition_pipe),
        "; benchmark did not"
      )
    }
  }

  invisible(NULL)
}

cat("\n============================================================\n")
say("PERSISTENCE-COMPARISON CONSOLE REPORT")
say("Benchmark: ", rank_label)
say("Focal species: ", focal_name_results)
cat("============================================================\n")

section("Run Scope")
say("Pipeline run ID: ", run_id)
say("Comparison directory: ", out_root)
say(
  "Compared ", fmt_int(uniqueN(sp_long$species)), " species across ",
  fmt_int(uniqueN(stage_sum$stage)), " matched stages."
)
say(
  "Population-unit trajectory rows: ", fmt_int(nrow(pu_long)),
  "; species trajectory rows: ", fmt_int(nrow(sp_long)), "."
)
say(
  "Removal axis span in comparison outputs: ",
  fmt_pct(min(x_to_pct(stage_sum$x), na.rm = TRUE)), " to ",
  fmt_pct(max(x_to_pct(stage_sum$x), na.rm = TRUE)), "."
)

section("Late-Stage Mean Persistence")
for (ii in seq_len(nrow(checkpoint_stage_wide))) {
  row <- checkpoint_stage_wide[ii]
  say(
    "Target ", fmt_pct(row$target_pct_removed), " removal",
    " (nearest stage ", row$stage, "; actual ", fmt_pct(row$actual_pct_removed), "):",
    " mean persistence ", fmt_num(row$mean_persist_pipe), " vs ",
    fmt_num(row$mean_persist_rank),
    " for ", rank_label,
    "; difference = ", fmt_num(row$mean_persist_diff_pipe_minus_rank),
    "; loss reduction = ",
    fmt_pct(row$persistence_loss_reduction_pct_vs_rank), "."
  )
  say(
    "  Median persistence: ", fmt_num(row$median_persist_pipe),
    " vs ", fmt_num(row$median_persist_rank),
    "; species with S > 0.5: ",
    fmt_pct(100 * row$frac_gt_05_pipe), " vs ",
    fmt_pct(100 * row$frac_gt_05_rank), "."
  )
}

section("Sequence-Integrated Species Advantage")
n_near_zero_auc <- n_species_auc - n_positive_auc - n_negative_auc
say(
  "Positive signed AUC species: ", fmt_int(n_positive_auc), " / ",
  fmt_int(n_species_auc), " (", fmt_pct(100 * n_positive_auc / n_species_auc), ")."
)
say(
  "Negative signed AUC species: ", fmt_int(n_negative_auc), " / ",
  fmt_int(n_species_auc), " (", fmt_pct(100 * n_negative_auc / n_species_auc), "); ",
  "near-zero species: ", fmt_int(n_near_zero_auc), "."
)
say(
  "Signed AUC distribution: mean = ",
  fmt_num(mean(species_auc$auc_diff_pipe_minus_rank, na.rm = TRUE)),
  "; median = ", fmt_num(metric_median(species_auc$auc_diff_pipe_minus_rank)),
  "; 10th to 90th percentile = ",
  fmt_num(metric_quantile(species_auc$auc_diff_pipe_minus_rank, 0.10)),
  " to ",
  fmt_num(metric_quantile(species_auc$auc_diff_pipe_minus_rank, 0.90)), "."
)
if (nrow(focal_auc_row)) {
  say(
    focal_name_results, " signed AUC rank: ",
    fmt_int(focal_auc_row$positive_auc_rank), " / ", fmt_int(n_species_auc),
    "; signed AUC = ", fmt_num(focal_auc_row$auc_diff_pipe_minus_rank),
    "; mean sequence gap = ", fmt_num(focal_auc_row$mean_gap_over_sequence), "."
  )
} else {
  say("Focal species was not found in the signed-AUC table.")
}

print_auc_examples(
  "Largest persistence-based advantages by signed AUC:",
  species_auc[order(-auc_diff_pipe_minus_rank)],
  n = 5L
)
print_auc_examples(
  "Largest benchmark advantages by signed AUC:",
  species_auc[order(auc_diff_pipe_minus_rank)],
  n = 5L
)

section("Area-Normalized Persistence Advantage")
say(
  "Positive area-normalized mean-gap species: ",
  fmt_int(n_area_positive), " / ", fmt_int(n_area_species),
  " (", fmt_pct(100 * n_area_positive / n_area_species), ")."
)
say(
  "Negative area-normalized mean-gap species: ",
  fmt_int(n_area_negative), " / ", fmt_int(n_area_species),
  " (", fmt_pct(100 * n_area_negative / n_area_species), ")."
)
say(
  "Area-normalized mean-gap distribution: mean = ",
  fmt_num(mean(area_norm_dt$mean_gap_area_normalized, na.rm = TRUE)),
  "; median = ", fmt_num(metric_median(area_norm_dt$mean_gap_area_normalized)),
  "; 10th to 90th percentile = ",
  fmt_num(metric_quantile(area_norm_dt$mean_gap_area_normalized, 0.10)),
  " to ",
  fmt_num(metric_quantile(area_norm_dt$mean_gap_area_normalized, 0.90)), "."
)
print_area_examples(
  "Largest persistence-based advantages after matching retained PU area:",
  area_norm_dt[order(-mean_gap_area_normalized)],
  n = 3L
)
print_area_examples(
  "Largest benchmark advantages after matching retained PU area:",
  area_norm_dt[order(mean_gap_area_normalized)],
  n = 3L
)

section("Threshold-Crossing Asymmetry")
for (thr in c(0.5, 0.1)) {
  row <- threshold_summary[abs(threshold - thr) < metric_eps][1L]
  say(
    "S < ", fmt_num(thr, digits = 1), ": benchmark-only crossings = ",
    fmt_int(row$n_rank_only),
    "; persistence-only crossings = ", fmt_int(row$n_pipe_only),
    "; both crossed = ", fmt_int(row$n_both_crossed),
    "; neither crossed = ", fmt_int(row$n_neither_crossed), "."
  )
  say(
    "  Among species crossing under both methods, median persistence-minus-benchmark crossing delay = ",
    fmt_x_delta(row$median_delay_both_crossed), "."
  )
}
print_cross_examples(
  threshold = 0.5,
  state = "rank_only",
  title = "Examples where benchmark crossed below S = 0.5 but persistence-based ranking did not:",
  n = 3L
)
print_cross_examples(
  threshold = 0.5,
  state = "pipe_only",
  title = "Examples where persistence-based ranking crossed below S = 0.5 but benchmark did not:",
  n = 3L
)
print_cross_examples(
  threshold = 0.1,
  state = "rank_only",
  title = "Examples where benchmark crossed below S = 0.1 but persistence-based ranking did not:",
  n = 3L
)
print_cross_examples(
  threshold = 0.1,
  state = "pipe_only",
  title = "Examples where persistence-based ranking crossed below S = 0.1 but benchmark did not:",
  n = 3L
)

section("Late-Checkpoint Tail Risk")
for (target_i in c(95, 99)) {
  for (thr in c(0.5, 0.1)) {
    row <- tail_summary[
      target_pct_removed == target_i &
        abs(threshold - thr) < metric_eps
    ][1L]
    say(
      "At target ", fmt_pct(target_i), " removal",
      " (actual ", fmt_pct(row$actual_pct_removed), "), threshold S = ",
      fmt_num(thr, digits = 1), ": persistence above / benchmark below = ",
      fmt_int(row$n_pipe_kept_above_rank_below),
      "; benchmark above / persistence below = ",
      fmt_int(row$n_rank_kept_above_pipe_below),
      "; both above = ", fmt_int(row$n_both_above),
      "; both below = ", fmt_int(row$n_both_below), "."
    )
  }
}
print_tail_examples(
  target_pct = 99,
  threshold = 0.5,
  direction = "pipe_above_rank_below",
  title = "Largest 99% checkpoint gaps where persistence stayed above S = 0.5 and benchmark fell below:",
  n = 3L
)
print_tail_examples(
  target_pct = 99,
  threshold = 0.5,
  direction = "rank_above_pipe_below",
  title = "Largest 99% checkpoint gaps where benchmark stayed above S = 0.5 and persistence fell below:",
  n = 3L
)

section("Population-Unit Redundancy")
for (target_i in c("no_PU_gt_0.9", "no_PU_gt_0.5", "expected_persisting_PU_lt_1")) {
  row <- redundancy_summary[target_id == target_i][1L]
  say(
    row$manuscript_label,
    ": benchmark-only losses = ", fmt_int(row$n_rank_only_lost),
    "; persistence-only losses = ", fmt_int(row$n_pipe_only_lost),
    "; both lost = ", fmt_int(row$n_both_lost),
    "; neither lost = ", fmt_int(row$n_neither_lost),
    "; median persistence-minus-benchmark delay = ",
    fmt_x_delta(row$median_censored_delay), "."
  )
}
print_redundancy_examples(
  target_id_i = "no_PU_gt_0.5",
  state = "rank_only_lost",
  title = "Examples where benchmark lost all PUs with P_PU > 0.5 but persistence-based ranking did not:",
  n = 3L
)
print_redundancy_examples(
  target_id_i = "no_PU_gt_0.5",
  state = "pipe_only_lost",
  title = "Examples where persistence-based ranking lost all PUs with P_PU > 0.5 but benchmark did not:",
  n = 3L
)

section("Files Checked Or Updated")
say("PU trajectories: ", out_pu)
say("Species trajectories: ", out_sp)
say("Stage summaries: ", out_stage)
say("Stage comparisons: ", out_compare)
say("The 4-panel persistence metrics comparison figure is printed below this console report.")

cat("\n============================================================\n")
say("END PERSISTENCE-COMPARISON CONSOLE REPORT")
cat("============================================================\n")
