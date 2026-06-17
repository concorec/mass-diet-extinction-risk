# Distance-connectivity stage driver.
source(file.path("R", "project_utils.R"))
source(file.path("R", "analysis_contract.R"))
source(file.path("R", "priority_distance_geometry.R"))
source(file.path("R", "priority_distance_graph.R"))

# The fragmentation stage identifies patches whose dispersal links may have
# changed. This driver rechecks only those local candidate patches, removes
# links beyond the species dispersal threshold, and rebuilds affected
# persistence units.

log_distance_connectivity_stage <- function(
  stage_index,
  touched_species_count,
  input_recheck_patch_count,
  removed_cell_count,
  remaining_patch_count,
  remaining_pu_count,
  remaining_alive_cell_count
) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  message(
    sprintf(
      paste(
        "[%s] distance_stage |",
        "stage=%d touched_species=%d input_recheck_patches=%d",
        "removed_cells=%d remaining_patches=%d remaining_pus=%d remaining_alive_cells=%d"
      ),
      timestamp,
      as.integer(stage_index),
      as.integer(touched_species_count),
      as.integer(input_recheck_patch_count),
      as.integer(removed_cell_count),
      as.integer(remaining_patch_count),
      as.integer(remaining_pu_count),
      as.integer(remaining_alive_cell_count)
    )
  )

  invisible(NULL)
}


# =====================================================================
# MAIN FUNCTION: run the distance-based connectivity stage
# =====================================================================

run_distance_connectivity_stage <- function(
  stage_index,
  patches_requiring_distance_recheck,
  patch_table,
  pu_graphs_by_key,
  alive_species_count_by_cell,
  patch_id_by_species_env,
  patch_cell_index_by_species_env,
  species_params,
  species_raster_template_stack
) {
  # -------------------------------------------------------------------
  # 1. Standardize table inputs
  # -------------------------------------------------------------------

  patch_table <- data.table::as.data.table(patch_table)

  species_params <- data.table::as.data.table(species_params)

  patches_requiring_distance_recheck <-
    data.table::as.data.table(patches_requiring_distance_recheck)

  input_recheck_patch_count <- nrow(patches_requiring_distance_recheck)


  # -------------------------------------------------------------------
  # 2. Validate package dependency
  # -------------------------------------------------------------------

  # Stop early if sf is unavailable, because local distance predicates require sf.
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("The distance stage requires the 'sf' package.")
  }


  # -------------------------------------------------------------------
  # 3. Initialize distance-stage accumulator
  # -------------------------------------------------------------------

  # Store cells that become globally empty during distance-stage updates.
  removed_cells_in_distance <- integer(0L)


  # -------------------------------------------------------------------
  # 4. Fast return if no patches were flagged for distance recheck
  # -------------------------------------------------------------------

  # If there are no recheck patches, there is no distance-stage work.
  if (nrow(patches_requiring_distance_recheck) == 0L) {
    # Count current patch-table rows.
    remaining_patch_count <- nrow(patch_table)

    # Count current unique species-PU combinations.
    remaining_pu_count <- data.table::uniqueN(patch_table, by = c("species", "pu_id"))

    # Count globally alive cells.
    remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

    log_distance_connectivity_stage(
      stage_index = stage_index,
      touched_species_count = 0L,
      input_recheck_patch_count = 0L,
      removed_cell_count = 0L,
      remaining_patch_count = remaining_patch_count,
      remaining_pu_count = remaining_pu_count,
      remaining_alive_cell_count = remaining_alive_cell_count
    )

    # Return unchanged state.
    return(list(
      patch_table = patch_table,
      pu_graphs = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      removed_cells_in_distance = integer(0L),
      touched_species = character()
    ))
  }


  # -------------------------------------------------------------------
  # 5. Restrict recheck patches to species available in the live env
  # -------------------------------------------------------------------

  # Read species currently represented by live cell -> patch vectors.
  species_available_in_env <- sort(ls(envir = patch_id_by_species_env, all.names = TRUE))

  # Keep only recheck rows whose species has a live species vector.
  patches_requiring_distance_recheck <- patches_requiring_distance_recheck[
    species %in% species_available_in_env
  ]

  # Identify touched species after filtering.
  touched_species <- sort(unique(as.character(
    patches_requiring_distance_recheck$species
  )))

  # If no touched species remain, exit cleanly.
  if (!length(touched_species)) {
    # Count current patch-table rows.
    remaining_patch_count <- nrow(patch_table)

    # Count current unique species-PU combinations.
    remaining_pu_count <- data.table::uniqueN(patch_table, by = c("species", "pu_id"))

    # Count globally alive cells.
    remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

    log_distance_connectivity_stage(
      stage_index = stage_index,
      touched_species_count = 0L,
      input_recheck_patch_count = input_recheck_patch_count,
      removed_cell_count = 0L,
      remaining_patch_count = remaining_patch_count,
      remaining_pu_count = remaining_pu_count,
      remaining_alive_cell_count = remaining_alive_cell_count
    )

    # Return unchanged state.
    return(list(
      patch_table = patch_table,
      pu_graphs = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      removed_cells_in_distance = integer(0L),
      touched_species = character()
    ))
  }


  # ===================================================================
  # 6. Process each touched species
  # ===================================================================

  # Loop over species that have at least one surviving recheck patch.
  for (species_name in touched_species) {
    # Read recheck patch IDs for this species.
    recheck_patch_ids_for_species <- sort(unique(as.integer(
      patches_requiring_distance_recheck$patch_id[
        patches_requiring_distance_recheck$species == species_name
      ]
    )))

    # Read this species' live cell -> patch_id vector before distance updates.
    species_patch_ids_before <- get(
      species_name,
      envir = patch_id_by_species_env,
      inherits = FALSE
    )

    # Create an editable copy for this distance stage.
    species_patch_ids_after <- species_patch_ids_before

    # Read current patch rows for this species.
    species_patch_rows_before <- patch_table[species == species_name]

    # Skip species with no surviving patch rows.
    if (nrow(species_patch_rows_before) == 0L) {
      next
    }

    # Keep only recheck patches still present in the current patch table.
    recheck_patch_ids_for_species <- intersect(
      recheck_patch_ids_for_species,
      species_patch_rows_before$patch_id
    )

    # Skip this species if no recheck patches survive.
    if (!length(recheck_patch_ids_for_species)) {
      next
    }

    # Identify current PUs containing at least one recheck patch.
    affected_pu_ids <- sort(unique(
      species_patch_rows_before$pu_id[
        species_patch_rows_before$patch_id %in% recheck_patch_ids_for_species
      ]
    ))

    # Skip this species if no affected PUs exist.
    if (!length(affected_pu_ids)) {
      next
    }

    # Read this species' dispersal distance in kilometers.
    dispersal_threshold_km <- get_species_dispersal_threshold_km(
      species_name = species_name,
      species_params = species_params
    )

    # Materialize a temporary current species raster from the live env state.
    species_patch_raster <- materialize_species_patch_raster(
      species_name = species_name,
      species_raster_template_stack = species_raster_template_stack,
      patch_id_by_species_env = patch_id_by_species_env
    )

    # Build the minimal candidate patch set needed for local geometry.
    candidate_patch_ids <- extract_candidate_patch_ids_for_distance(
      species_name = species_name,
      recheck_patch_ids_for_species = recheck_patch_ids_for_species,
      affected_pu_ids = affected_pu_ids,
      pu_graphs_by_key = pu_graphs_by_key
    )

    # Keep only candidate patches still present in the patch table.
    candidate_patch_ids <- intersect(
      candidate_patch_ids,
      species_patch_rows_before$patch_id
    )

    # Skip if no candidate patches remain.
    if (!length(candidate_patch_ids)) {
      next
    }

    # Build local polygons only for candidate patches.
    candidate_patch_polygons <- build_local_patch_polygons(
      species_patch_raster = species_patch_raster,
      candidate_patch_ids = candidate_patch_ids
    )

    # Skip if no candidate geometry could be reconstructed.
    if (nrow(candidate_patch_polygons) == 0L) {
      next
    }

    # Build sparse Boolean within-distance lookup for candidate polygons.
    distance_predicate_lookup <- build_distance_predicate_lookup(
      candidate_patch_polygons = candidate_patch_polygons,
      dispersal_threshold_km = dispersal_threshold_km
    )


    # -----------------------------------------------------------------
    # 6A. Prepare affected/unaffected patch rows and graph collectors
    # -----------------------------------------------------------------

    # Keep patch rows from unaffected PUs unchanged.
    unaffected_species_patch_rows <- species_patch_rows_before[
      !(pu_id %in% affected_pu_ids),
      list(
        species = species,
        patch_id = patch_id,
        pu_id = pu_id,
        patch_area_km2 = patch_area_km2
      )
    ]

    # Allocate rebuilt patch-row collector for affected PUs.
    rebuilt_patch_rows_for_affected_pus <- vector("list", 0L)

    # Allocate rebuilt graph collector for affected PUs.
    rebuilt_graphs_for_affected_pus <- list()

    # Track patches dropped after distance-edge filtering.
    dropped_patch_ids_from_distance <- integer(0L)

    # Initialize next available PU ID for possible PU splits.
    next_available_pu_id <- if (nrow(species_patch_rows_before)) {
      max(species_patch_rows_before$pu_id)
    } else {
      0L
    }

    # Read this species' PU-area threshold.
    pu_area_threshold <- species_params[
      species == species_name
    ]$min_population_area_km2[1L]


    # -----------------------------------------------------------------
    # 6B. Filter and rebuild each affected PU
    # -----------------------------------------------------------------

    # Loop over affected PUs.
    for (affected_pu_id in affected_pu_ids) {
      # Read current patch rows for this PU.
      current_rows_in_this_pu <- species_patch_rows_before[
        pu_id == affected_pu_id
      ]

      # Build graph key for this species-PU.
      graph_key <- paste0(species_name, "|", affected_pu_id)

      # Read current graph object for this PU.
      old_pu_graph <- pu_graphs_by_key[[graph_key]]

      if (is.null(old_pu_graph)) {
        stop(
          "Missing PU graph during distance-stage update for species '",
          species_name,
          "' and pu_id ",
          affected_pu_id,
          ". This indicates inconsistent evolving state between patch_table ",
          "and pu_graphs_by_key."
        )
      }

      # Restrict recheck patches to this PU.
      recheck_patch_ids_in_pu <- intersect(
        recheck_patch_ids_for_species,
        current_rows_in_this_pu$patch_id
      )

      # If this affected PU has no recheck patches after restriction, keep it.
      if (!length(recheck_patch_ids_in_pu)) {
        # Materialize lazy graphs so the distance stage leaves standard CSR graphs.
        graph_to_keep <- if (is_lazy_added_node_overlay_graph(old_pu_graph)) {
          materialize_lazy_added_node_overlay_graph(old_pu_graph)
        } else {
          old_pu_graph
        }

        rebuilt_graphs_for_affected_pus[[graph_key]] <- graph_to_keep

        rebuilt_patch_rows_for_affected_pus[[
          length(rebuilt_patch_rows_for_affected_pus) + 1L
        ]] <- current_rows_in_this_pu

        # Continue to the next affected PU.
        next
      }

      # Filter distance-invalid edges from this PU graph.
      filtered_pu_graph <- filter_distance_invalid_edges_for_pu(
        pu_graph = old_pu_graph,
        recheck_patch_ids_in_pu = recheck_patch_ids_in_pu,
        distance_predicate_lookup = distance_predicate_lookup
      )

      # Build named patch_id -> patch_area vector for this PU.
      patch_area_by_patch <- setNames(
        object = current_rows_in_this_pu$patch_area_km2,
        nm = as.character(current_rows_in_this_pu$patch_id)
      )

      rebuilt_pu <- rebuild_pu_after_edge_filter(
        pu_graph = filtered_pu_graph,
        patch_area_by_patch = patch_area_by_patch,
        pu_area_threshold = pu_area_threshold,
        next_available_pu_id = next_available_pu_id
      )

      # Carry forward the next available PU ID.
      next_available_pu_id <- rebuilt_pu$next_available_pu_id

      # Accumulate dropped patch IDs from subthreshold split components.
      if (length(rebuilt_pu$dropped_patch_ids)) {
        dropped_patch_ids_from_distance <- c(
          dropped_patch_ids_from_distance,
          rebuilt_pu$dropped_patch_ids
        )
      }

      # Store surviving PU graph objects.
      for (surviving_graph in rebuilt_pu$surviving_pu_graphs) {
        # Build graph key for this surviving component.
        surviving_graph_key <- paste0(species_name, "|", surviving_graph$pu_id)

        # Store standard CSR graph for this surviving component.
        rebuilt_graphs_for_affected_pus[[surviving_graph_key]] <- list(
          species = species_name,
          pu_id = as.integer(surviving_graph$pu_id),
          id2patch = as.integer(surviving_graph$id2patch),
          row_ptr = as.integer(surviving_graph$row_ptr),
          col_idx = as.integer(surviving_graph$col_idx)
        )
      }

      for (surviving_graph in rebuilt_pu$surviving_pu_graphs) {
        # Read patch IDs in this surviving component.
        component_patch_ids <- as.integer(surviving_graph$id2patch)

        # Extract component patch rows in graph patch order.
        component_rows <- current_rows_in_this_pu[
          match(component_patch_ids, patch_id)
        ]

        # Assign this component's PU ID.
        component_rows[, pu_id := as.integer(surviving_graph$pu_id)]

        # Store component rows.
        rebuilt_patch_rows_for_affected_pus[[
          length(rebuilt_patch_rows_for_affected_pus) + 1L
        ]] <- component_rows
      }
    }


    # -----------------------------------------------------------------
    # 6C. Apply patch drops from distance-stage PU rebuilding
    # -----------------------------------------------------------------

    # Deduplicate dropped patch IDs for this species.
    dropped_patch_ids_from_distance <- sort(unique(as.integer(
      dropped_patch_ids_from_distance
    )))

    # Clear dropped distance-stage patches from this species vector.
    if (length(dropped_patch_ids_from_distance)) {
      species_patch_ids_after[
        species_patch_ids_after %in% dropped_patch_ids_from_distance
      ] <- NA_integer_
    }


    # -----------------------------------------------------------------
    # 6D. Rebuild this species' patch rows
    # -----------------------------------------------------------------

    # Combine rebuilt affected-PU rows, or create an empty table.
    rebuilt_species_patch_rows <- if (length(rebuilt_patch_rows_for_affected_pus)) {
      data.table::rbindlist(
        rebuilt_patch_rows_for_affected_pus,
        use.names = TRUE,
        fill = TRUE
      )
    } else {
      data.table::data.table(
        species = character(),
        patch_id = integer(),
        pu_id = integer(),
        patch_area_km2 = numeric()
      )
    }

    # Remove old rows for this species from the global patch table.
    patch_table <- patch_table[species != species_name]

    # Append unchanged rows and rebuilt affected-PU rows.
    patch_table <- data.table::rbindlist(
      list(
        patch_table,
        unaffected_species_patch_rows,
        rebuilt_species_patch_rows
      ),
      use.names = TRUE,
      fill = TRUE
    )


    # -----------------------------------------------------------------
    # 6E. Replace affected PU graphs
    # -----------------------------------------------------------------

    # Remove old graph entries for affected original PUs.
    for (affected_pu_id in affected_pu_ids) {
      pu_graphs_by_key[[paste0(species_name, "|", affected_pu_id)]] <- NULL
    }

    # Insert rebuilt graph entries.
    for (graph_key in names(rebuilt_graphs_for_affected_pus)) {
      pu_graphs_by_key[[graph_key]] <- rebuilt_graphs_for_affected_pus[[graph_key]]
    }


    # -----------------------------------------------------------------
    # 6F. Update alive-species counts for cells where this species vanished
    # -----------------------------------------------------------------

    # Identify cells where this species was present before but absent after.
    species_cells_removed_in_distance <- which(
      !is.na(species_patch_ids_before) &
        is.na(species_patch_ids_after)
    )

    # Read global alive-species counts before removing this species.
    count_before_species_loss <- alive_species_count_by_cell[
      species_cells_removed_in_distance
    ]

    # Identify cells that remain globally alive after losing this species.
    cells_still_alive_after_species_loss <- species_cells_removed_in_distance[
      count_before_species_loss > 0L
    ]

    # Identify cells that become globally empty after losing this species.
    cells_becoming_empty <- species_cells_removed_in_distance[
      count_before_species_loss == 1L
    ]

    # Decrement alive-species counts for cells that were still globally alive.
    if (length(cells_still_alive_after_species_loss)) {
      alive_species_count_by_cell[cells_still_alive_after_species_loss] <-
        alive_species_count_by_cell[cells_still_alive_after_species_loss] - 1L
    }

    # Add newly empty cells to the distance-stage accumulator.
    if (length(cells_becoming_empty)) {
      removed_cells_in_distance <- c(
        removed_cells_in_distance,
        cells_becoming_empty
      )
    }


    # -----------------------------------------------------------------
    # 6G. Write updated live species state back to environments
    # -----------------------------------------------------------------

    assign(
      species_name,
      species_patch_ids_after,
      envir = patch_id_by_species_env
    )

    # Rebuild this species' compact patch -> cells index.
    updated_species_patch_index <- rebuild_species_patch_index(
      species_patch_ids_after
    )

    assign(
      species_name,
      updated_species_patch_index,
      envir = patch_cell_index_by_species_env
    )
  }


  # -------------------------------------------------------------------
  # 7. Finalize distance-stage removed cells
  # -------------------------------------------------------------------

  # Deduplicate and sort globally empty cells created by the distance stage.
  removed_cells_in_distance <- sort(
    unique(as.integer(removed_cells_in_distance))
  )


  # -------------------------------------------------------------------
  # 8. Compute final stage summaries and log one line
  # -------------------------------------------------------------------

  # Count remaining patch-table rows.
  remaining_patch_count <- nrow(patch_table)

  # Count remaining unique species-PU combinations.
  remaining_pu_count <- data.table::uniqueN(
    patch_table,
    by = c("species", "pu_id")
  )

  # Count remaining globally alive cells.
  remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

  log_distance_connectivity_stage(
    stage_index = stage_index,
    touched_species_count = length(touched_species),
    input_recheck_patch_count = input_recheck_patch_count,
    removed_cell_count = length(removed_cells_in_distance),
    remaining_patch_count = remaining_patch_count,
    remaining_pu_count = remaining_pu_count,
    remaining_alive_cell_count = remaining_alive_cell_count
  )


  # -------------------------------------------------------------------
  # 9. Validate that no lazy overlay graphs remain
  # -------------------------------------------------------------------
  #
  # This is not a timing or diagnostic log. It protects the invariant that the
  # distance stage consumes lazy overlay graphs and returns standard CSR graphs.
  #

  # Count any lazy overlay graphs still present after the distance update.
  remaining_lazy_overlay_graphs <- sum(vapply(
    pu_graphs_by_key,
    is_lazy_added_node_overlay_graph,
    logical(1L)
  ))

  if (remaining_lazy_overlay_graphs > 0L) {
    stop(
      "Distance stage completed with ",
      remaining_lazy_overlay_graphs,
      " lazy overlay graph(s) still present. ",
      "This means at least one lazy overlay was not consumed into a standard CSR graph."
    )
  }


  # -------------------------------------------------------------------
  # 10. Return updated evolving state
  # -------------------------------------------------------------------

  # Return the updated state needed by the outer pipeline.
  list(
    patch_table = patch_table,
    pu_graphs = pu_graphs_by_key,
    alive_species_count_by_cell = alive_species_count_by_cell,
    removed_cells_in_distance = removed_cells_in_distance,
    touched_species = touched_species
  )
}
