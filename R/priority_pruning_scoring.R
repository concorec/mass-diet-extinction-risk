# =====================================================================
# HELPER: compute current PU-level log scores
# =====================================================================

compute_pu_log_scores <- function(pu_area_table, species_params) {
  # Join species-level parameters onto current PU-area rows.
  pu_data <- species_params[pu_area_table, on = .(species)]

  # Read current PU areas in km^2.
  pu_area_km2 <- pu_data$pu_area_km2

  # Read species density values.
  species_density <- pu_data$density

  # Read fitted Gompertz alpha parameters.
  persistence_a <- pu_data$a_pred

  # Read fitted Gompertz beta parameters.
  persistence_b <- pu_data$b_pred

  # Read the canonical species-specific minimum population-area thresholds.
  pu_area_threshold_km2 <- pu_data$min_population_area_km2

  # Read species names aligned with PU rows.
  species_name <- pu_data$species

  # Read PU IDs aligned with PU rows.
  pu_id <- pu_data$pu_id

  # Compute each PU's area above the quasi-extinction area threshold.
  area_above_threshold <- pu_area_km2 - pu_area_threshold_km2

  # Take log area-above-threshold once for reuse.
  log_area_above_threshold <- log(area_above_threshold)

  # Compute log(a * density^(-b)).
  log_a_density_term <-
    log(persistence_a) - persistence_b * log(species_density)

  # Compute log(a * density^(-b) * area_above_threshold^(-b)).
  log_inner_term <-
    log_a_density_term -
    persistence_b * log_area_above_threshold

  # Convert the inner term back to ordinary scale.
  inner_term <- exp(log_inner_term)

  # Compute current PU persistence.
  pu_persistence <- exp(-inner_term)

  # Compute log(1 - PU persistence) stably.
  log_one_minus_persistence <- log1p(-pu_persistence)

  # Sum log failure probabilities across PUs within each species.
  species_sum_log_failure <- tapply(log_one_minus_persistence, species_name, sum)

  # Remove each PU's own failure term to get the product over all other PUs.
  log_other_pu_term <-
    species_sum_log_failure[species_name] - log_one_minus_persistence

  # Compute log derivative of PU persistence with respect to PU area.
  log_persistence_derivative <-
    log_a_density_term +
    log(persistence_b) -
    (persistence_b + 1) * log_area_above_threshold -
    inner_term

  # Combine derivative and redundancy terms into the final PU log marginal score.
  log_score <- log_persistence_derivative + log_other_pu_term

  # Stop if any PU score is non-finite.
  if (any(!is.finite(log_score))) {
    stop("compute_pu_log_scores() produced non-finite PU scores.")
  }

  # Return one log-score row per PU.
  data.table(
    species = species_name,
    pu_id = pu_id,
    log_score = log_score
  )
}



# =====================================================================
# HELPER: update patch-score cache for dirty species only
# =====================================================================

update_patch_score_cache <- function(
  patch_table,
  species_params,
  patch_score_cache = NULL,
  dirty_species = NULL
) {
  # Identify species that currently have at least one surviving patch.
  current_species <- sort(unique(as.character(patch_table$species)))

  # Stop if the patch table is unexpectedly empty.
  if (!length(current_species)) {
    stop("update_patch_score_cache() received an empty patch_table.")
  }

  # Initialize an empty cache if this is the first scoring pass.
  if (is.null(patch_score_cache)) {
    # Create an empty list for per-species score vectors.
    patch_score_cache <- list()

    # Force all current species to be recomputed on the first pass.
    dirty_species <- current_species
  }

  # If dirty species were not supplied, conservatively recompute all current species.
  if (is.null(dirty_species)) {
    dirty_species <- current_species
  }

  # Keep only dirty species that still exist in the current patch table.
  dirty_species <- intersect(
    sort(unique(as.character(dirty_species))),
    current_species
  )

  # Identify cached species that no longer have surviving patches.
  stale_species <- setdiff(names(patch_score_cache), current_species)

  # Remove stale species from the cache.
  if (length(stale_species)) {
    patch_score_cache[stale_species] <- NULL
  }

  # Recompute score vectors only for dirty species.
  if (length(dirty_species)) {
    # Restrict the patch table to species whose state changed.
    dirty_patch_table <- patch_table[
      species %in% dirty_species
    ]

    # Stop if dirty species were requested but no rows exist for them.
    if (!nrow(dirty_patch_table)) {
      stop("Dirty species were requested, but dirty_patch_table is empty.")
    }

    # Aggregate dirty patch rows to current PU areas.
    dirty_pu_area_table <- dirty_patch_table[
      ,
      list(pu_area_km2 = sum(patch_area_km2)),
      by = .(species, pu_id)
    ]

    # Compute current PU-level log scores for dirty species.
    dirty_pu_score_table <- compute_pu_log_scores(
      pu_area_table = dirty_pu_area_table,
      species_params = species_params
    )

    # Attach each dirty patch to its current PU log score.
    dirty_patch_score_table <- dirty_patch_table[
      ,
      list(species, patch_id, pu_id)
    ][
      dirty_pu_score_table,
      on = .(species, pu_id),
      nomatch = 0L
    ][
      ,
      list(
        species = as.character(species),
        patch_id = as.integer(patch_id),
        score = log_score
      )
    ]

    # Stop if no patch scores were produced for dirty species.
    if (!nrow(dirty_patch_score_table)) {
      stop("update_patch_score_cache() produced no scores for dirty species.")
    }

    # Identify dirty species that actually received score rows.
    scored_dirty_species <- sort(unique(dirty_patch_score_table$species))

    # Identify dirty species missing from the score table.
    missing_dirty_species <- setdiff(dirty_species, scored_dirty_species)

    # Stop if any dirty surviving species did not receive scores.
    if (length(missing_dirty_species)) {
      stop(
        "No patch scores were produced for dirty species: ",
        paste(missing_dirty_species, collapse = ", ")
      )
    }

    # Stop if patch IDs cannot be used as positive integer vector indices.
    if (any(is.na(dirty_patch_score_table$patch_id) |
            dirty_patch_score_table$patch_id < 1L)) {
      stop("Patch IDs must be positive integers for direct patch-score indexing.")
    }

    # Stop if any patch score is non-finite.
    if (any(!is.finite(dirty_patch_score_table$score))) {
      stop("update_patch_score_cache() produced non-finite patch scores.")
    }

    # Split row numbers by species for per-species score-vector construction.
    row_index_by_species <- split(
      seq_len(nrow(dirty_patch_score_table)),
      dirty_patch_score_table$species
    )

    # Rebuild one direct patch_id -> score vector per dirty species.
    for (species_name in names(row_index_by_species)) {
      # Read row positions for this species.
      row_index <- row_index_by_species[[species_name]]

      # Read patch IDs for this species.
      patch_ids <- dirty_patch_score_table$patch_id[row_index]

      # Read log scores for this species.
      scores <- dirty_patch_score_table$score[row_index]

      # Stop if the score table contains duplicate patch IDs for one species.
      if (anyDuplicated(patch_ids)) {
        stop("Duplicate patch IDs found while building score cache for: ", species_name)
      }

      # Allocate direct lookup vector, filling absent patch IDs with -Inf.
      score_by_patch_id <- rep(-Inf, max(patch_ids))

      # Store finite patch scores at their patch-ID positions.
      score_by_patch_id[patch_ids] <- scores

      # Store this species' direct lookup vector in the cache.
      patch_score_cache[[species_name]] <- list(
        score_by_patch_id = score_by_patch_id
      )
    }
  }

  # Stop if the score cache is unexpectedly empty.
  if (!length(patch_score_cache)) {
    stop("update_patch_score_cache() produced an empty patch-score cache.")
  }

  # Return the updated score cache.
  patch_score_cache
}



# =====================================================================
# HELPER: score frontier cells
# =====================================================================

score_frontier_cells <- function(
  frontier_cells,
  patch_scores_by_species,
  patch_id_by_species_env
) {
  # Count frontier cells to score.
  n_frontier <- length(frontier_cells)

  # Return an empty score vector if there are no frontier cells.
  if (!n_frontier) {
    return(numeric(0L))
  }

  # Allocate per-frontier-cell maximum log contribution for log-sum-exp.
  max_log_term <- rep(-Inf, n_frontier)

  # Allocate per-frontier-cell shifted exponential sum for log-sum-exp.
  sum_exp_terms <- numeric(n_frontier)

  # Track whether each frontier position has at least one species contribution.
  has_score <- logical(n_frontier)

  # Loop over species with current cached patch-score vectors.
  for (species_name in names(patch_scores_by_species)) {
    # Read this species' current full-grid cell -> patch_id vector.
    species_patch_ids <- get(
      species_name,
      envir = patch_id_by_species_env,
      inherits = FALSE
    )

    # Subset the species patch vector to current frontier cells.
    patch_ids_in_frontier_cells <- species_patch_ids[frontier_cells]

    # Identify frontier positions where this species is present.
    valid_position <- which(
      !is.na(patch_ids_in_frontier_cells) &
        patch_ids_in_frontier_cells >= 1L
    )

    # Skip this species if it is absent from the current frontier.
    if (!length(valid_position)) {
      next
    }

    # Read this species' direct patch_id -> log-score lookup vector.
    score_by_patch_id <-
      patch_scores_by_species[[species_name]]$score_by_patch_id

    # Read patch IDs at frontier positions where this species is present.
    patch_ids_valid <- as.integer(patch_ids_in_frontier_cells[valid_position])

    # Keep only patch IDs that can index the score vector.
    in_range <- !is.na(patch_ids_valid) &
      patch_ids_valid >= 1L &
      patch_ids_valid <= length(score_by_patch_id)

    # Skip this species if all candidate patch IDs are stale or out of range.
    if (!any(in_range)) {
      next
    }

    # Keep frontier positions with in-range patch IDs.
    valid_position <- valid_position[in_range]

    # Keep corresponding patch IDs.
    patch_ids_valid <- patch_ids_valid[in_range]

    # Look up log scores for the valid patch IDs.
    log_score <- score_by_patch_id[patch_ids_valid]

    # Identify finite score lookups.
    finite_score <- is.finite(log_score)

    # Skip this species if no finite score lookups remain.
    if (!any(finite_score)) {
      next
    }

    # Keep only frontier positions with finite log scores.
    valid_position <- valid_position[finite_score]

    # Keep only finite log scores.
    log_score <- log_score[finite_score]

    # Identify positions receiving their first species contribution.
    first_score <- !has_score[valid_position]

    # Process first contributions.
    if (any(first_score)) {
      # Extract frontier positions receiving first contributions.
      first_pos <- valid_position[first_score]

      # Store first contribution as current max log term.
      max_log_term[first_pos] <- log_score[first_score]

      # Initialize shifted sum to one for first contribution.
      sum_exp_terms[first_pos] <- 1

      # Mark these frontier positions as scored.
      has_score[first_pos] <- TRUE
    }

    # Identify positions receiving additional species contributions.
    repeat_score <- !first_score

    # Process repeated contributions.
    if (any(repeat_score)) {
      # Extract frontier positions receiving repeated contributions.
      repeat_pos <- valid_position[repeat_score]

      # Extract new log scores for those positions.
      repeat_log_score <- log_score[repeat_score]

      # Read old max log terms for those positions.
      old_max <- max_log_term[repeat_pos]

      # Identify where the new score becomes the new maximum.
      new_is_larger <- repeat_log_score > old_max

      # Handle cases where old maximum remains larger.
      if (any(!new_is_larger)) {
        # Extract positions where old max remains the max.
        pos <- repeat_pos[!new_is_larger]

        # Add new shifted contribution to the running sum.
        sum_exp_terms[pos] <-
          sum_exp_terms[pos] +
          exp(repeat_log_score[!new_is_larger] - max_log_term[pos])
      }

      # Handle cases where new contribution becomes the maximum.
      if (any(new_is_larger)) {
        # Extract positions where new log score is largest.
        pos <- repeat_pos[new_is_larger]

        # Extract the new maximum values.
        new_max <- repeat_log_score[new_is_larger]

        # Store old maximum values before replacement.
        old_max_for_pos <- max_log_term[pos]

        # Rescale existing shifted sum to the new maximum and add the new term.
        sum_exp_terms[pos] <-
          sum_exp_terms[pos] * exp(old_max_for_pos - new_max) + 1

        # Replace the max log term with the new maximum.
        max_log_term[pos] <- new_max
      }
    }
  }

  # Allocate the final score vector.
  frontier_scores <- rep(-Inf, n_frontier)

  # Identify positions that received at least one valid species contribution.
  scored_position <- has_score & sum_exp_terms > 0 & is.finite(max_log_term)

  # Convert log-sum-exp marginal loss to the script's negative-loss score convention.
  frontier_scores[scored_position] <-
    -(max_log_term[scored_position] + log(sum_exp_terms[scored_position]))

  # Stop if any frontier cell has no valid finite score.
  if (any(!is.finite(frontier_scores))) {
    # Identify bad frontier positions.
    bad <- which(!is.finite(frontier_scores))

    # Report example global raster-cell IDs.
    stop(
      "score_frontier_cells() produced non-finite frontier scores for ",
      length(bad),
      " frontier cells. This indicates missing/stale patch-score lookups. ",
      "Example raster cell ids: ",
      paste(head(frontier_cells[bad], 10L), collapse = ", ")
    )
  }

  # Return one score per frontier cell, aligned with frontier_cells.
  frontier_scores
}



# =====================================================================
# HELPER: choose which frontier cells to remove
# =====================================================================

choose_frontier_cells_to_remove <- function(
  frontier_cells,
  frontier_scores,
  cells_to_remove_per_step
) {
  # Count available frontier cells.
  n_frontier <- length(frontier_cells)

  # Stop if the requested removal batch size is invalid.
  if (cells_to_remove_per_step < 1L) {
    stop("cells_to_remove_per_step must be >= 1.")
  }

  # If the request is at least the full frontier, return the full frontier deterministically.
  if (cells_to_remove_per_step >= n_frontier) {
    # Order higher score first, then lower cell index.
    ordered_positions <- order(-frontier_scores, frontier_cells)

    # Return ordered frontier cell IDs.
    return(frontier_cells[ordered_positions])
  }

  # Find the cutoff score for the top requested number of frontier cells.
  cutoff_score <- Rfast::nth(
    frontier_scores,
    k = cells_to_remove_per_step,
    descending = TRUE,
    index.return = FALSE
  )

  # Keep frontier positions whose score is at least the cutoff score.
  selected_positions <- which(frontier_scores >= cutoff_score)

  # If ties over-selected cells, break ties deterministically.
  if (length(selected_positions) > cells_to_remove_per_step) {
    # Order tied candidates by higher score, then lower cell index.
    tie_break_order <- order(
      -frontier_scores[selected_positions],
      frontier_cells[selected_positions]
    )

    # Keep exactly the requested number of selected positions.
    selected_positions <-
      selected_positions[tie_break_order][seq_len(cells_to_remove_per_step)]
  }

  # Convert selected positions back to global raster-cell IDs.
  frontier_cells[selected_positions]
}



# =====================================================================
# HELPER: summarize removed area by species-specific patch
# =====================================================================

summarize_removed_patch_area <- function(
  selected_cells,
  removed_cell_area,
  species_names,
  patch_id_by_species_env
) {
  # Combine one per-species summary table into a single data.table.
  rbindlist(
    # Loop over all species present in the current pruning stage.
    lapply(species_names, function(species_name) {
      # Read this species' patch IDs in the selected cells.
      patch_ids <- get(
        species_name,
        envir = patch_id_by_species_env,
        inherits = FALSE
      )[selected_cells]

      # Mark selected cells where this species is present.
      valid <- !is.na(patch_ids)

      # Return no rows if this species is absent from all selected cells.
      if (!any(valid)) {
        return(NULL)
      }

      # Sum removed cell area by patch ID for this species.
      area_sum <- rowsum(
        removed_cell_area[valid],
        patch_ids[valid],
        reorder = FALSE
      )

      # Return one row per affected species-specific patch.
      data.table(
        species = species_name,
        patch_id = as.integer(rownames(area_sum)),
        area_removed = as.numeric(area_sum[, 1])
      )
    }),
    use.names = TRUE,
    fill = TRUE
  )
}



# =====================================================================
# HELPER: find connected components among alive nodes in one PU graph
