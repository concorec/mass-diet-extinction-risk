# =====================================================================
# PRUNING-STAGE CHUNK
# =====================================================================
#
# Runs one full pruning stage: repeated frontier-cell pruning, live state
# updates, removed-cell tracking, changed-patch tracking, and early stopping
# when no removable frontier cells remain.
#
# =====================================================================


# =====================================================================
# HELPER: print one pruning-stage log message
# =====================================================================
#
# This helper prints exactly one summary message per completed pruning stage.
#
# It does not perform any expensive spatial reconstruction. All values passed
# to it are already available from the final stage state or from the stage
# accumulators.
#

log_pruning_stage <- function(
  stage_index,                      # global stage index
  completed_pruning_iterations,      # number of pruning iterations actually completed
  removed_cell_count,                # number of unique cells removed during this pruning stage
  changed_patch_count,               # number of unique patches changed during this pruning stage
  frontier_exhausted,                # whether pruning stopped because the frontier was exhausted
  remaining_patch_count,             # number of patch-table rows remaining after pruning
  remaining_pu_count,                # number of species-PU combinations remaining after pruning
  remaining_alive_cell_count         # number of globally alive cells remaining after pruning
) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  message(
    sprintf(
      paste(
        "[%s] pruning_stage |",
        "stage=%d completed_iterations=%d removed_cells=%d changed_patches=%d",
        "frontier_exhausted=%s remaining_patches=%d remaining_pus=%d remaining_alive_cells=%d"
      ),
      timestamp,
      as.integer(stage_index),
      as.integer(completed_pruning_iterations),
      as.integer(removed_cell_count),
      as.integer(changed_patch_count),
      if (isTRUE(frontier_exhausted)) "TRUE" else "FALSE",
      as.integer(remaining_patch_count),
      as.integer(remaining_pu_count),
      as.integer(remaining_alive_cell_count)
    )
  )

  invisible(NULL)
}


# =====================================================================
# MAIN FUNCTION: run one full pruning stage
# =====================================================================
#
# This function is the stage-level driver for the pruning portion of the
# spatial prioritization pipeline.
#
# It repeatedly calls run_pruning_iteration(), updating the live state after
# each call. It does not write rasters or masks. Instead, it returns the
# per-iteration removed-cell vectors so the outer pipeline can record removal
# order in memory and write one final removal-order raster later.
#
run_pruning_stage <- function(
  stage_index,                      # global stage index for this pruning stage
  pruning_iterations_per_stage,     # maximum pruning iterations to run in this stage
  cells_to_remove_per_iteration,    # number of frontier cells to remove per iteration
  patch_table,                      # current patch table
  pu_graphs_by_key,                 # current PU graph objects keyed by "species|pu_id"
  alive_species_count_by_cell,      # current number of surviving species in each cell
  patch_id_by_species_env,          # live environment: species -> cell-to-patch vector
  patch_cell_index_by_species_env,  # live environment: species -> patch-to-cells index
  cell_area_by_cell,                # area of each raster cell in km^2
  frontier_state,                   # current incremental frontier state
  rook_neighbor_index,              # compact rook-neighbor lookup
  species_params                    # species-level thresholds and persistence parameters
) {
  # -------------------------------------------------------------------
  # 1) Validate scalar stage-level arguments
  # -------------------------------------------------------------------

  if (!is.numeric(stage_index) || length(stage_index) != 1L || is.na(stage_index)) {
    stop("stage_index must be one non-missing numeric value.")
  }

  if (!is.numeric(pruning_iterations_per_stage) ||
      length(pruning_iterations_per_stage) != 1L ||
      is.na(pruning_iterations_per_stage) ||
      pruning_iterations_per_stage < 1) {
    stop("pruning_iterations_per_stage must be a single integer >= 1.")
  }

  if (!is.numeric(cells_to_remove_per_iteration) ||
      length(cells_to_remove_per_iteration) != 1L ||
      is.na(cells_to_remove_per_iteration) ||
      cells_to_remove_per_iteration < 1) {
    stop("cells_to_remove_per_iteration must be a single integer >= 1.")
  }

  pruning_iterations_per_stage <- as.integer(pruning_iterations_per_stage)

  cells_to_remove_per_iteration <- as.integer(cells_to_remove_per_iteration)


  # -------------------------------------------------------------------
  # 2) Prepare pruning-stage state shared across iterations
  # -------------------------------------------------------------------

  # Build the species-name vector once for this pruning stage.
  species_names <- sort(ls(envir = patch_id_by_species_env, all.names = TRUE))

  # Start this pruning stage with no patch-score cache.
  patch_score_cache <- NULL

  # Force the first pruning iteration to score all species currently in patch_table.
  dirty_species_for_scoring <- sort(unique(as.character(patch_table$species)))


  # -------------------------------------------------------------------
  # 3) Allocate per-iteration accumulators
  # -------------------------------------------------------------------

  # Allocate one removed-cell vector slot per possible pruning iteration.
  removed_cells_each_iteration <- vector("list", pruning_iterations_per_stage)

  # Allocate one changed-patch table slot per possible pruning iteration.
  changed_patches_each_iteration <- vector("list", pruning_iterations_per_stage)

  # Count how many pruning iterations actually completed.
  completed_pruning_iterations <- 0L

  # Track whether the stage ended because the frontier was exhausted.
  frontier_exhausted <- FALSE


  # -------------------------------------------------------------------
  # 4) Run repeated pruning iterations
  # -------------------------------------------------------------------

  # Loop over the maximum allowed pruning iterations for this stage.
  for (iteration_in_stage in seq_len(pruning_iterations_per_stage)) {
    # Run one pruning iteration using the current live state.
    pruning_result <- run_pruning_iteration(
      prune_iteration                 = iteration_in_stage,
      cells_to_remove_per_step        = cells_to_remove_per_iteration,
      patch_table                     = patch_table,
      pu_graphs_by_key                = pu_graphs_by_key,
      alive_species_count_by_cell     = alive_species_count_by_cell,
      cell_area_by_cell               = cell_area_by_cell,
      species_params                  = species_params,
      patch_id_by_species_env         = patch_id_by_species_env,
      patch_cell_index_by_species_env = patch_cell_index_by_species_env,
      frontier_state                  = frontier_state,
      rook_neighbor_index             = rook_neighbor_index,
      species_names                   = species_names,
      patch_score_cache               = patch_score_cache,
      dirty_species_for_scoring       = dirty_species_for_scoring
    )

    # Replace the stage patch table with the updated patch table.
    patch_table <- pruning_result$patch_table

    # Replace the stage PU graph list with the updated graph list.
    pu_graphs_by_key <- pruning_result$pu_graphs

    # Replace the stage alive-species count vector with the updated vector.
    alive_species_count_by_cell <- pruning_result$alive_species_count_by_cell

    # Replace the stage frontier state with the updated frontier state.
    frontier_state <- pruning_result$frontier_state

    # Carry forward the updated patch-score cache.
    patch_score_cache <- pruning_result$patch_score_cache

    # Carry forward the dirty-species set for the next pruning iteration.
    dirty_species_for_scoring <- pruning_result$dirty_species_for_scoring

    # Extract cells that became newly empty during this pruning iteration.
    removed_cells_this_iteration <- pruning_result$newly_empty_cells

    # Extract patches changed or removed during this pruning iteration.
    changed_patches_this_iteration <- pruning_result$changed_patches

    # Store this iteration's newly empty cells.
    removed_cells_each_iteration[[iteration_in_stage]] <- removed_cells_this_iteration

    # Store this iteration's changed patches.
    changed_patches_each_iteration[[iteration_in_stage]] <- changed_patches_this_iteration

    # Mark this pruning iteration as completed.
    completed_pruning_iterations <- as.integer(iteration_in_stage)

    # If this iteration removed no cells, the frontier is exhausted.
    if (length(removed_cells_this_iteration) == 0L) {
      # Record that the pruning stage ended by frontier exhaustion.
      frontier_exhausted <- TRUE

      # Stop the pruning-stage loop early.
      break
    }
  }


  # -------------------------------------------------------------------
  # 5) Trim per-iteration accumulators to completed iterations
  # -------------------------------------------------------------------

  # If at least one pruning iteration completed, keep only completed slots.
  if (completed_pruning_iterations > 0L) {
    # Keep removed-cell vectors for completed iterations only.
    removed_cells_each_iteration <-
      removed_cells_each_iteration[seq_len(completed_pruning_iterations)]

    # Keep changed-patch tables for completed iterations only.
    changed_patches_each_iteration <-
      changed_patches_each_iteration[seq_len(completed_pruning_iterations)]
  } else {
    # Use an empty list if no pruning iterations completed.
    removed_cells_each_iteration <- list()

    # Use an empty list if no pruning iterations completed.
    changed_patches_each_iteration <- list()
  }


  # -------------------------------------------------------------------
  # 6) Consolidate removed cells across the full pruning stage
  # -------------------------------------------------------------------

  # If no iterations completed, no cells were removed.
  if (!length(removed_cells_each_iteration)) {
    # Store an empty integer vector for stage-level removed cells.
    removed_cells_in_stage <- integer(0L)
  } else {
    # Combine, deduplicate, sort, and coerce removed cells to integer IDs.
    removed_cells_in_stage <- sort(
      unique(
        as.integer(
          unlist(
            removed_cells_each_iteration,
            use.names = FALSE
          )
        )
      )
    )
  }


  # -------------------------------------------------------------------
  # 7) Consolidate changed patches across the full pruning stage
  # -------------------------------------------------------------------

  # Keep only changed-patch objects that are data frames with at least one row.
  non_empty_changed_patch_tables <- Filter(
    f = function(x) {
      is.data.frame(x) && nrow(x) > 0L
    },
    x = changed_patches_each_iteration
  )

  # If no patches changed, return an empty changed-patch table.
  if (!length(non_empty_changed_patch_tables)) {
    # Create an empty two-column changed-patch table.
    changed_patches_in_stage <- data.table::data.table(
      species = character(),
      patch_id = integer()
    )
  } else {
    # Combine all changed-patch tables and deduplicate by species and patch ID.
    changed_patches_in_stage <- unique(
      data.table::rbindlist(
        non_empty_changed_patch_tables,
        use.names = TRUE,
        fill = TRUE
      ),
      by = c("species", "patch_id")
    )
  }


  # -------------------------------------------------------------------
  # 8) Compute final stage summaries for logging
  # -------------------------------------------------------------------

  # Count remaining patch-table rows.
  remaining_patch_count <- nrow(patch_table)

  # Count remaining unique species-PU combinations.
  remaining_pu_count <- data.table::uniqueN(
    patch_table,
    by = c("species", "pu_id")
  )

  # Count cells that still contain at least one surviving species.
  remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

  # Count unique cells removed during this pruning stage.
  removed_cell_count <- length(removed_cells_in_stage)

  # Count unique patches changed during this pruning stage.
  changed_patch_count <- nrow(changed_patches_in_stage)


  # -------------------------------------------------------------------
  # 9) Print exactly one pruning-stage summary log message
  # -------------------------------------------------------------------

  # Print the concise pruning-stage summary.
  log_pruning_stage(
    stage_index                  = stage_index,
    completed_pruning_iterations = completed_pruning_iterations,
    removed_cell_count           = removed_cell_count,
    changed_patch_count          = changed_patch_count,
    frontier_exhausted           = frontier_exhausted,
    remaining_patch_count        = remaining_patch_count,
    remaining_pu_count           = remaining_pu_count,
    remaining_alive_cell_count   = remaining_alive_cell_count
  )


  # -------------------------------------------------------------------
  # 10) Return updated state and pruning-stage bookkeeping
  # -------------------------------------------------------------------

  # Return the updated live state and compact stage bookkeeping.
  list(
    patch_table                   = patch_table,
    pu_graphs                     = pu_graphs_by_key,
    alive_species_count_by_cell   = alive_species_count_by_cell,
    frontier_state                = frontier_state,
    removed_cells_each_iteration  = removed_cells_each_iteration,
    removed_cells_in_stage        = removed_cells_in_stage,
    changed_patches_in_stage      = changed_patches_in_stage,
    completed_pruning_iterations  = completed_pruning_iterations,
    frontier_exhausted            = frontier_exhausted
  )
}