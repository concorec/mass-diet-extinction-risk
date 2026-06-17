# Pruning-iteration driver.
source(file.path("R", "project_utils.R"))
source(file.path("R", "analysis_contract.R"))
source(file.path("R", "priority_pruning_frontier.R"))
source(file.path("R", "priority_pruning_scoring.R"))
source(file.path("R", "priority_pruning_graph.R"))

# =====================================================================
# PRUNING-ITERATION CHUNK
# =====================================================================
#
# Helper functions for one frontier-cell pruning iteration: scoring,
# partial selection, patch loss, CSR graph repair, and threshold bookkeeping.
#
# =====================================================================


# =====================================================================
# HELPER: print one pruning-iteration log message
# =====================================================================

log_pruning_iteration <- function(
  prune_iteration,
  frontier_cell_count,
  selected_cell_count = NULL,
  newly_empty_cell_count = NULL,
  touched_patch_count = NULL,
  dropped_patch_threshold_count = NULL,
  dropped_component_threshold_count = NULL,
  dropped_pu_threshold_count = NULL,
  dropped_pu_count = NULL,
  new_pu_count = NULL,
  frontier_exhausted = FALSE
) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  # If the frontier is exhausted, print the short stopping message.
  if (isTRUE(frontier_exhausted)) {
    # Print exactly one frontier-exhausted log line.
    message(
      sprintf(
        paste(
          "[%s] prune_step |",
          "iter=%d frontier_exhausted frontier_cells=%d"
        ),
        timestamp,
        as.integer(prune_iteration),
        as.integer(frontier_cell_count)
      )
    )

    return(invisible(NULL))
  }

  # Print exactly one normal pruning-iteration summary line.
  message(
    sprintf(
      paste(
        "[%s] prune_step |",
        "iter=%d frontier_cells=%d selected_cells=%d newly_empty_cells=%d",
        "touched_patches=%d dropped_patch_threshold=%d",
        "dropped_component_threshold=%d dropped_pu_threshold=%d",
        "dropped_pus=%d new_pus_created=%d"
      ),
      timestamp,
      as.integer(prune_iteration),
      as.integer(frontier_cell_count),
      as.integer(selected_cell_count),
      as.integer(newly_empty_cell_count),
      as.integer(touched_patch_count),
      as.integer(dropped_patch_threshold_count),
      as.integer(dropped_component_threshold_count),
      as.integer(dropped_pu_threshold_count),
      as.integer(dropped_pu_count),
      as.integer(new_pu_count)
    )
  )

  invisible(NULL)
}

run_pruning_iteration <- function(
  prune_iteration,
  cells_to_remove_per_step,
  patch_table,
  pu_graphs_by_key,
  alive_species_count_by_cell,
  cell_area_by_cell,
  species_params,
  patch_id_by_species_env,
  patch_cell_index_by_species_env,
  frontier_state,
  rook_neighbor_index,
  species_names,
  patch_score_cache,
  dirty_species_for_scoring
) {
  # -------------------------------------------------------------------
  # STEP 1. Identify the current frontier
  # -------------------------------------------------------------------

  # Read all currently removable frontier cells.
  frontier_cells <- get_frontier_cells(frontier_state)

  # Stop this pruning stage if the frontier is too small for a full batch.
  if (length(frontier_cells) < cells_to_remove_per_step) {
    # Print the single frontier-exhausted log line.
    log_pruning_iteration(
      prune_iteration = prune_iteration,
      frontier_cell_count = length(frontier_cells),
      frontier_exhausted = TRUE
    )

    # Return unchanged state and empty change bookkeeping.
    return(list(
      patch_table = patch_table,
      pu_graphs = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      frontier_state = frontier_state,
      newly_empty_cells = integer(0L),
      changed_patches = data.table(species = character(), patch_id = integer()),
      patch_score_cache = patch_score_cache,
      dirty_species_for_scoring = dirty_species_for_scoring
    ))
  }

  # Snapshot current patch -> PU mapping before this iteration changes anything.
  patch_to_pu_before_step <- patch_table[, .(species, patch_id, pu_id)]


  # -------------------------------------------------------------------
  # STEP 2. Score the frontier and choose cells to remove
  # -------------------------------------------------------------------

  # Update cached patch scores only for species whose state changed.
  patch_score_cache <- update_patch_score_cache(
    patch_table = patch_table,
    species_params = species_params,
    patch_score_cache = patch_score_cache,
    dirty_species = dirty_species_for_scoring
  )

  # Score every current frontier cell.
  frontier_scores <- score_frontier_cells(
    frontier_cells = frontier_cells,
    patch_scores_by_species = patch_score_cache,
    patch_id_by_species_env = patch_id_by_species_env
  )

  selected_cells <- choose_frontier_cells_to_remove(
    frontier_cells = frontier_cells,
    frontier_scores = frontier_scores,
    cells_to_remove_per_step = cells_to_remove_per_step
  )


  # -------------------------------------------------------------------
  # STEP 3. Remove selected cells at the global cell level
  # -------------------------------------------------------------------

  # Store selected cells as integer cell IDs.
  newly_empty_cells <- as.integer(selected_cells)

  # Mark selected cells as globally removed.
  alive_species_count_by_cell[selected_cells] <- 0L

  # Read exact cell areas for patch-area decrementing.
  removed_cell_area <- cell_area_by_cell[selected_cells]


  # -------------------------------------------------------------------
  # STEP 4. Translate removed cells into removed area by patch
  # -------------------------------------------------------------------

  # Summarize selected-cell area by species-specific patch.
  removed_area_by_patch <- summarize_removed_patch_area(
    selected_cells = selected_cells,
    removed_cell_area = removed_cell_area,
    species_names = species_names,
    patch_id_by_species_env = patch_id_by_species_env
  )

  # Identify species directly present in at least one selected cell.
  directly_touched_species <- unique(removed_area_by_patch$species)

  # Clear selected cells from only species that occupied selected cells.
  for (species_name in directly_touched_species) {
    # Read this species' mutable cell -> patch_id vector.
    patch_ids <- get(
      species_name,
      envir = patch_id_by_species_env,
      inherits = FALSE
    )

    # Clear this species from selected cells.
    patch_ids[selected_cells] <- NA_integer_

    # Write updated species vector back to the environment.
    assign(species_name, patch_ids, envir = patch_id_by_species_env)
  }


  # -------------------------------------------------------------------
  # STEP 5. Decrement patch areas
  # -------------------------------------------------------------------

  # Join removed-area rows to the pre-step patch -> PU mapping.
  patch_area_decrements <- removed_area_by_patch[
    patch_to_pu_before_step,
    on = .(species, patch_id),
    nomatch = 0L
  ][
    ,
    .(species, patch_id, pu_id, area_removed)
  ]

  # Subtract removed area from affected patch-table rows.
  patch_table[
    patch_area_decrements,
    patch_area_km2 := patch_area_km2 - i.area_removed,
    on = .(species, patch_id, pu_id)
  ]

  # Record directly touched patches.
  touched_patches <- unique(
    patch_area_decrements[, .(species, patch_id)]
  )


  # -------------------------------------------------------------------
  # STEP 6. Drop patches below the patch-level threshold
  # -------------------------------------------------------------------

  # Identify touched patches that now fall below their species patch threshold.
  patches_dropped_by_patch_threshold <- patch_table[
    touched_patches,
    on = .(species, patch_id)
  ][
    species_params,
    on = .(species)
  ][
    !area_exceeds_threshold(patch_area_km2, min_patch_area_km2),
    .(species, patch_id)
  ]

  # Remove patch-threshold failures from the patch table.
  patch_table <- patch_table[
    !patches_dropped_by_patch_threshold,
    on = .(species, patch_id)
  ]


  # -------------------------------------------------------------------
  # STEP 7. Repair affected PUs after patch losses
  # -------------------------------------------------------------------

  # Create one pending patch-drop bucket per species.
  pending_patch_drops_by_species <-
    setNames(vector("list", length(species_names)), species_names)

  # Split patch-threshold failures by species.
  dropped_patch_ids_by_species <- split(
    patches_dropped_by_patch_threshold$patch_id,
    patches_dropped_by_patch_threshold$species
  )

  # Add patch-threshold failures to pending patch-drop buckets.
  for (species_name in names(dropped_patch_ids_by_species)) {
    # Append this species' dropped patch IDs.
    pending_patch_drops_by_species[[species_name]] <- c(
      pending_patch_drops_by_species[[species_name]],
      dropped_patch_ids_by_species[[species_name]]
    )
  }

  # Identify original species-PU pairs containing newly dropped patches.
  affected_species_pu_pairs <- merge(
    patches_dropped_by_patch_threshold,
    patch_to_pu_before_step,
    by = c("species", "patch_id")
  )

  # Allocate collector for component-threshold patch drops.
  component_threshold_drop_records <- vector("list", 0L)

  # Initialize count of new PUs created by splitting.
  new_pus_created <- 0L

  # Process each species with patch-threshold drops.
  for (species_name in unique(affected_species_pu_pairs$species)) {
    # Mark rows for this species.
    species_rows <- affected_species_pu_pairs$species == species_name

    # Split dropped patch IDs by original PU.
    dropped_patch_ids_by_pu <- split(
      affected_species_pu_pairs$patch_id[species_rows],
      affected_species_pu_pairs$pu_id[species_rows]
    )

    # Find next available PU ID for this species.
    next_available_pu_id <- if (any(patch_table$species == species_name)) {
      max(patch_table$pu_id[patch_table$species == species_name])
    } else {
      0L
    }

    # Read this species' PU-area threshold.
    pu_area_threshold <- species_params[.(species_name), min_population_area_km2]

    # Build current patch_id -> patch_area lookup for this species.
    patch_area_lookup <- patch_table[
      species == species_name,
      .(patch_id, patch_area_km2)
    ]

    # Process each affected original PU.
    for (pu_id_string in names(dropped_patch_ids_by_pu)) {
      # Build graph key for this species-PU pair.
      graph_key <- paste0(species_name, "|", pu_id_string)

      # Read current CSR graph for this PU.
      pu_graph <- pu_graphs_by_key[[graph_key]]

      # Skip missing graphs.
      if (is.null(pu_graph)) {
        next
      }

      # Read node-aligned patch IDs.
      node_patch_ids <- pu_graph$id2patch

      # Read patch IDs that died in this PU.
      dead_patch_ids <- dropped_patch_ids_by_pu[[pu_id_string]]

      # Initialize all graph nodes as alive.
      alive_nodes <- rep(TRUE, length(node_patch_ids))

      # Mark nodes belonging to dead patches as not alive.
      alive_nodes[node_patch_ids %in% dead_patch_ids] <- FALSE

      # Allocate node-aligned patch-area vector.
      node_patch_area <- numeric(length(node_patch_ids))

      # Match graph-node patch IDs to current surviving patch-area rows.
      area_match <- fastmatch::fmatch(
        node_patch_ids,
        patch_area_lookup$patch_id
      )

      # Identify graph nodes with current patch-area rows.
      valid_area <- !is.na(area_match)

      # Fill current patch areas for surviving graph nodes.
      node_patch_area[valid_area] <-
        patch_area_lookup$patch_area_km2[area_match[valid_area]]

      # Rebuild this PU after patch loss.
      pu_update <- rebuild_pu_after_patch_loss(
        pu_graph = pu_graph,
        alive_nodes = alive_nodes,
        patch_area = node_patch_area,
        pu_area_threshold = pu_area_threshold,
        next_available_pu_id = next_available_pu_id
      )

      # Carry forward updated next available PU ID.
      next_available_pu_id <- pu_update$next_available_pu_id

      # Read patches dropped because surviving components were too small.
      component_dropped_patch_ids <- pu_update$dropped_patch_ids

      # Queue and remove component-threshold failures.
      if (length(component_dropped_patch_ids)) {
        # Record component-threshold patch drops.
        component_threshold_drop_records[[length(component_threshold_drop_records) + 1L]] <-
          data.table(species = species_name, patch_id = component_dropped_patch_ids)

        # Add component-threshold patch drops to pending cell-level removals.
        pending_patch_drops_by_species[[species_name]] <- c(
          pending_patch_drops_by_species[[species_name]],
          component_dropped_patch_ids
        )

        # Remove component-threshold failures from the patch table.
        patch_table <- patch_table[
          !(species == species_name & patch_id %in% component_dropped_patch_ids)
        ]
      }

      # Read surviving component graphs.
      surviving_pu_graphs <- pu_update$surviving_pu_graphs

      # Remove original graph if no component survives.
      if (!length(surviving_pu_graphs)) {
        # Delete the old graph entry.
        pu_graphs_by_key[[graph_key]] <- NULL
      } else {
        # Count new PUs if this PU split into multiple surviving components.
        if (length(surviving_pu_graphs) > 1L) {
          new_pus_created <- new_pus_created + (length(surviving_pu_graphs) - 1L)
        }

        # Store all surviving component graphs.
        for (j in seq_along(surviving_pu_graphs)) {
          # Read one surviving component graph.
          component_graph <- surviving_pu_graphs[[j]]

          # Store first component under original graph key.
          if (j == 1L) {
            # Replace original graph with first surviving component.
            pu_graphs_by_key[[graph_key]] <- list(
              species = species_name,
              pu_id = component_graph$pu_id,
              id2patch = component_graph$id2patch,
              row_ptr = component_graph$row_ptr,
              col_idx = component_graph$col_idx
            )
          } else {
            # Build graph key for split-off component.
            new_key <- paste0(species_name, "|", component_graph$pu_id)

            # Store split-off component under new graph key.
            pu_graphs_by_key[[new_key]] <- list(
              species = species_name,
              pu_id = component_graph$pu_id,
              id2patch = component_graph$id2patch,
              row_ptr = component_graph$row_ptr,
              col_idx = component_graph$col_idx
            )
          }

          # Update patch-table PU IDs for patches in this component.
          patch_table[
            species == species_name & patch_id %in% component_graph$id2patch,
            pu_id := component_graph$pu_id
          ]
        }
      }
    }
  }

  # Collapse component-threshold patch drops into one unique table.
  patches_dropped_by_component_threshold <- if (length(component_threshold_drop_records)) {
    unique(
      rbindlist(component_threshold_drop_records, use.names = TRUE, fill = TRUE),
      by = c("species", "patch_id")
    )
  } else {
    data.table(species = character(), patch_id = integer())
  }


  # -------------------------------------------------------------------
  # STEP 8. Drop whole PUs that now fail the PU-level threshold
  # -------------------------------------------------------------------

  # Build candidate patches whose area loss could make a PU subthreshold.
  candidate_patches <- unique(
    rbindlist(
      list(
        removed_area_by_patch[, .(species, patch_id)],
        patches_dropped_by_patch_threshold
      ),
      use.names = TRUE,
      fill = TRUE
    ),
    by = c("species", "patch_id")
  )

  # Identify original PUs containing those candidate patches.
  affected_pus_before <- unique(
    patch_to_pu_before_step[
      candidate_patches,
      on = .(species, patch_id),
      nomatch = 0L
    ][
      ,
      .(species, pu_id)
    ],
    by = c("species", "pu_id")
  )

  # Recover original patches from affected original PUs, keeping only current survivors.
  affected_patches_that_still_exist <- patch_to_pu_before_step[
    affected_pus_before,
    on = .(species, pu_id),
    nomatch = 0L
  ][
    ,
    .(species, patch_id)
  ][
    patch_table,
    on = .(species, patch_id),
    nomatch = 0L
  ]

  candidate_pus_after <- unique(
    affected_patches_that_still_exist[, .(species, pu_id)],
    by = c("species", "pu_id")
  )

  # Recompute total area for only candidate PUs.
  pu_area_after_update <- patch_table[
    candidate_pus_after,
    .(pu_area_km2 = sum(patch_area_km2)),
    on = .(species, pu_id),
    by = .EACHI
  ]

  # Identify candidate PUs now below the PU-area threshold.
  pus_dropped_by_pu_threshold <- species_params[
    pu_area_after_update,
    on = .(species)
  ][
    !area_exceeds_threshold(pu_area_km2, min_population_area_km2),
    .(species, pu_id)
  ]

  # Allocate collector for patches dropped with whole PUs.
  pu_threshold_drop_records <- vector("list", 0L)

  # Process each PU that fails the PU threshold.
  for (i in seq_len(nrow(pus_dropped_by_pu_threshold))) {
    # Read species name for this dropped PU.
    species_name <- pus_dropped_by_pu_threshold$species[i]

    # Read PU ID for this dropped PU.
    pu_id_value <- pus_dropped_by_pu_threshold$pu_id[i]

    # Build graph key for this dropped PU.
    graph_key <- paste0(species_name, "|", pu_id_value)

    # Read current patch IDs in the dropped PU.
    patch_ids_in_dropped_pu <- patch_table[
      species == species_name & pu_id == pu_id_value,
      patch_id
    ]

    # If patches remain in the dropped PU, queue them for cell-level removal.
    if (length(patch_ids_in_dropped_pu)) {
      # Record patches lost because this whole PU was dropped.
      pu_threshold_drop_records[[length(pu_threshold_drop_records) + 1L]] <-
        data.table(species = species_name, patch_id = patch_ids_in_dropped_pu)

      # Add these patches to the species' pending-drop bucket.
      pending_patch_drops_by_species[[species_name]] <- c(
        pending_patch_drops_by_species[[species_name]],
        patch_ids_in_dropped_pu
      )
    }

    # Remove this PU's CSR graph.
    pu_graphs_by_key[[graph_key]] <- NULL

    # Remove this PU's patch rows from the patch table.
    patch_table <- patch_table[
      !(species == species_name & pu_id == pu_id_value)
    ]
  }

  # Collapse PU-threshold patch drops into one unique table.
  patches_dropped_by_pu_threshold <- if (length(pu_threshold_drop_records)) {
    unique(
      rbindlist(pu_threshold_drop_records, use.names = TRUE, fill = TRUE),
      by = c("species", "patch_id")
    )
  } else {
    data.table(species = character(), patch_id = integer())
  }


  # -------------------------------------------------------------------
  # STEP 9. Apply accumulated patch drops to cell-level state
  # -------------------------------------------------------------------

  # Loop over species with possible pending dropped patches.
  for (species_name in names(pending_patch_drops_by_species)) {
    # Deduplicate dropped patch IDs for this species.
    patch_ids_to_drop <- unique(pending_patch_drops_by_species[[species_name]])

    # Skip species with no pending dropped patches.
    if (!length(patch_ids_to_drop)) {
      next
    }

    patch_ids_to_drop <- sort(unique(as.integer(patch_ids_to_drop)))

    # Keep only valid positive patch IDs.
    patch_ids_to_drop <- patch_ids_to_drop[
      !is.na(patch_ids_to_drop) &
        patch_ids_to_drop >= 1L
    ]

    # Skip if no valid patch IDs remain.
    if (!length(patch_ids_to_drop)) {
      next
    }

    # Read species-specific patch_id -> cells index.
    patch_index <- get(
      species_name,
      envir = patch_cell_index_by_species_env,
      inherits = FALSE
    )

    # Skip if no patch-cell index exists for this species.
    if (is.null(patch_index)) {
      next
    }

    # Read species-specific full-grid cell -> patch_id vector.
    patch_ids <- get(
      species_name,
      envir = patch_id_by_species_env,
      inherits = FALSE
    )

    # Look up all cells belonging to dropped patches for this species.
    dropped_patch_cells <- get_patch_cells_batch(
      patch_index = patch_index,
      patch_ids = patch_ids_to_drop
    )

    # Skip if dropped patch IDs have no indexed cells.
    if (!nrow(dropped_patch_cells)) {
      next
    }

    # Deduplicate cells to clear for this species.
    cells_to_clear <- sort(unique(as.integer(dropped_patch_cells$cell)))

    # Keep only valid raster-cell IDs.
    cells_to_clear <- cells_to_clear[
      !is.na(cells_to_clear) &
        cells_to_clear >= 1L &
        cells_to_clear <= length(alive_species_count_by_cell)
    ]

    # Skip if no valid cells remain.
    if (!length(cells_to_clear)) {
      next
    }

    # Read global alive-species counts before removing this species.
    count_before_species_loss <- alive_species_count_by_cell[cells_to_clear]

    # Identify cells still globally alive before this species is removed.
    still_alive_cells <- cells_to_clear[count_before_species_loss > 0L]

    # Identify cells that become globally empty after this species is removed.
    cells_becoming_empty <- cells_to_clear[count_before_species_loss == 1L]

    # Decrement alive-species counts for cells that were still globally alive.
    if (length(still_alive_cells)) {
      alive_species_count_by_cell[still_alive_cells] <-
        alive_species_count_by_cell[still_alive_cells] - 1L
    }

    # Add cascading newly empty cells to the iteration-level list.
    if (length(cells_becoming_empty)) {
      newly_empty_cells <- c(newly_empty_cells, cells_becoming_empty)
    }

    # Clear this species from all cells belonging to dropped patches.
    patch_ids[cells_to_clear] <- NA_integer_

    # Write updated species vector back to the environment.
    assign(species_name, patch_ids, envir = patch_id_by_species_env)
  }


  # -------------------------------------------------------------------
  # STEP 10. Final bookkeeping for the stage-level orchestrator
  # -------------------------------------------------------------------

  # Deduplicate and sort all cells that became globally empty in this iteration.
  newly_empty_cells <- sort(unique(as.integer(newly_empty_cells)))

  # Update incremental frontier state using only newly empty cells.
  frontier_state <- update_frontier_state_after_empty_cells(
    frontier_state = frontier_state,
    newly_empty_cells = newly_empty_cells,
    rook_neighbor_index = rook_neighbor_index
  )

  # Collect every patch changed or removed during this pruning iteration.
  changed_patches <- unique(
    rbindlist(
      list(
        removed_area_by_patch[, .(species, patch_id)],
        patches_dropped_by_patch_threshold,
        patches_dropped_by_component_threshold,
        patches_dropped_by_pu_threshold
      ),
      use.names = TRUE,
      fill = TRUE
    ),
    by = c("species", "patch_id")
  )

  dirty_species_for_next_iteration <-
    sort(unique(as.character(changed_patches$species)))


  # -------------------------------------------------------------------
  # STEP 11. Print exactly one pruning-iteration log line
  # -------------------------------------------------------------------

  # Print the single concise production log message for this iteration.
  log_pruning_iteration(
    prune_iteration = prune_iteration,
    frontier_cell_count = length(frontier_cells),
    selected_cell_count = length(selected_cells),
    newly_empty_cell_count = length(newly_empty_cells),
    touched_patch_count = nrow(touched_patches),
    dropped_patch_threshold_count = nrow(patches_dropped_by_patch_threshold),
    dropped_component_threshold_count = nrow(patches_dropped_by_component_threshold),
    dropped_pu_threshold_count = nrow(patches_dropped_by_pu_threshold),
    dropped_pu_count = nrow(pus_dropped_by_pu_threshold),
    new_pu_count = as.integer(new_pus_created),
    frontier_exhausted = FALSE
  )


  # -------------------------------------------------------------------
  # STEP 12. Return updated state
  # -------------------------------------------------------------------

  # Return exactly the objects needed by the pruning-stage orchestrator.
  list(
    patch_table = patch_table,
    pu_graphs = pu_graphs_by_key,
    alive_species_count_by_cell = alive_species_count_by_cell,
    frontier_state = frontier_state,
    newly_empty_cells = newly_empty_cells,
    changed_patches = changed_patches,
    patch_score_cache = patch_score_cache,
    dirty_species_for_scoring = dirty_species_for_next_iteration
  )
}
