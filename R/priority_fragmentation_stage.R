# Fragmentation-stage driver.
source(file.path("R", "project_utils.R"))
source(file.path("R", "analysis_contract.R"))
source(file.path("R", "priority_fragmentation_patches.R"))
source(file.path("R", "priority_fragmentation_graph.R"))

# After each pruning stage, only changed patches are relabeled. Newly split
# fragments are filtered by the minimum patch-size rule, affected persistence
# units are rebuilt, and patches that remain eligible for distance filtering are
# returned to the distance-connectivity stage.

log_fragmentation_stage <- function(
  stage_index,
  touched_species_count,
  input_changed_patch_count,
  removed_cell_count,
  recheck_patch_count,
  remaining_patch_count,
  remaining_pu_count,
  remaining_alive_cell_count
) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  message(
    sprintf(
      paste(
        "[%s] fragmentation_stage |",
        "stage=%d touched_species=%d input_changed_patches=%d",
        "removed_cells=%d recheck_patches=%d",
        "remaining_patches=%d remaining_pus=%d remaining_alive_cells=%d"
      ),
      timestamp,
      as.integer(stage_index),
      as.integer(touched_species_count),
      as.integer(input_changed_patch_count),
      as.integer(removed_cell_count),
      as.integer(recheck_patch_count),
      as.integer(remaining_patch_count),
      as.integer(remaining_pu_count),
      as.integer(remaining_alive_cell_count)
    )
  )

  invisible(NULL)
}


# =====================================================================
# MAIN FUNCTION: run the fragmentation stage
# =====================================================================
#
# This function updates only the species and patches changed during the
# preceding pruning stage.
#
# It:
#
#   1. checks changed patches for rook-contiguity fragmentation
#   2. relabels newly disconnected fragments
#   3. drops fragments below the patch-size threshold
#   4. rebuilds only affected PU graph structures
#   5. drops PU components below the PU-size threshold
#   6. updates live environment-based species state
#   7. returns cells removed during fragmentation
#   8. returns surviving patches needing distance-connectivity recheck
#
# It does not write rasters or masks.
#

run_fragmentation_stage <- function(
  stage_index,
  changed_patches_in_stage,
  patch_table,
  pu_graphs_by_key,
  alive_species_count_by_cell,
  patch_id_by_species_env,
  patch_cell_index_by_species_env,
  species_params,
  cell_area_by_cell,
  rook_neighbor_index
) {
  # -------------------------------------------------------------------
  # 1. Standardize inputs
  # -------------------------------------------------------------------

  patch_table <- data.table::as.data.table(patch_table)

  species_params <- data.table::as.data.table(species_params)

  changed_patches_in_stage <- data.table::as.data.table(changed_patches_in_stage)


  # -------------------------------------------------------------------
  # 2. Validate the required rook-neighbor index
  # -------------------------------------------------------------------

  if (
    is.null(rook_neighbor_index) ||
      is.null(rook_neighbor_index$start) ||
      is.null(rook_neighbor_index$end) ||
      is.null(rook_neighbor_index$to) ||
      is.null(rook_neighbor_index$n_cells)
  ) {
    stop(
      "run_fragmentation_stage() requires a valid rook_neighbor_index. ",
      "Build it once with build_rook_neighbor_index() and pass it into ",
      "the fragmentation stage."
    )
  }


  # -------------------------------------------------------------------
  # 3. Initialize fragmentation-stage accumulators
  # -------------------------------------------------------------------

  # Store cells that become globally empty during fragmentation.
  removed_cells_in_fragmentation <- integer(0L)

  input_changed_patch_count <- nrow(changed_patches_in_stage)


  # -------------------------------------------------------------------
  # 4. Fast return if no patches changed during pruning
  # -------------------------------------------------------------------

  # If no changed patches were supplied, fragmentation has nothing to update.
  if (nrow(changed_patches_in_stage) == 0L) {
    # Count current patch rows.
    remaining_patch_count <- nrow(patch_table)

    # Count current unique species-PU combinations.
    remaining_pu_count <- data.table::uniqueN(patch_table, by = c("species", "pu_id"))

    # Count current globally alive cells.
    remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

    log_fragmentation_stage(
      stage_index = stage_index,
      touched_species_count = 0L,
      input_changed_patch_count = 0L,
      removed_cell_count = 0L,
      recheck_patch_count = 0L,
      remaining_patch_count = remaining_patch_count,
      remaining_pu_count = remaining_pu_count,
      remaining_alive_cell_count = remaining_alive_cell_count
    )

    # Return unchanged state and empty fragmentation bookkeeping.
    return(list(
      patch_table = patch_table,
      pu_graphs = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      removed_cells_in_fragmentation = integer(0L),
      patches_requiring_distance_recheck = data.table::data.table(
        species = character(),
        patch_id = integer()
      ),
      touched_species = character()
    ))
  }


  # -------------------------------------------------------------------
  # 5. Restrict changed patches to species available in the live env
  # -------------------------------------------------------------------

  # Read species currently represented in the live cell -> patch environment.
  species_available_in_env <- sort(ls(envir = patch_id_by_species_env, all.names = TRUE))

  # Keep only changed patches whose species has a live environment vector.
  changed_patches_in_stage <- changed_patches_in_stage[
    species %in% species_available_in_env
  ]

  # Identify species touched by the changed patches that remain after filtering.
  touched_species <- sort(unique(as.character(changed_patches_in_stage$species)))

  # If filtering removed all touched species, no fragmentation update is possible.
  if (!length(touched_species)) {
    # Count current patch rows.
    remaining_patch_count <- nrow(patch_table)

    # Count current unique species-PU combinations.
    remaining_pu_count <- data.table::uniqueN(patch_table, by = c("species", "pu_id"))

    # Count current globally alive cells.
    remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

    log_fragmentation_stage(
      stage_index = stage_index,
      touched_species_count = 0L,
      input_changed_patch_count = input_changed_patch_count,
      removed_cell_count = 0L,
      recheck_patch_count = 0L,
      remaining_patch_count = remaining_patch_count,
      remaining_pu_count = remaining_pu_count,
      remaining_alive_cell_count = remaining_alive_cell_count
    )

    # Return unchanged state and empty fragmentation bookkeeping.
    return(list(
      patch_table = patch_table,
      pu_graphs = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      removed_cells_in_fragmentation = integer(0L),
      patches_requiring_distance_recheck = data.table::data.table(
        species = character(),
        patch_id = integer()
      ),
      touched_species = character()
    ))
  }


  # -------------------------------------------------------------------
  # 6. Allocate one distance-recheck table slot per touched species
  # -------------------------------------------------------------------

  # Create one output slot per touched species.
  recheck_patch_rows_by_species <- vector("list", length(touched_species))


  # ===================================================================
  # 7. Process each touched species
  # ===================================================================

  # Loop over species touched by pruning-stage patch changes.
  for (species_index in seq_along(touched_species)) {
    species_name <- touched_species[species_index]

    # Read changed patch IDs for this species.
    changed_patch_ids_for_species <- sort(unique(as.integer(
      changed_patches_in_stage$patch_id[changed_patches_in_stage$species == species_name]
    )))

    species_patch_ids_before <- get(
      species_name,
      envir = patch_id_by_species_env,
      inherits = FALSE
    )

    # Make an editable copy of the species cell -> patch_id vector.
    species_patch_ids_after <- species_patch_ids_before

    # Read current patch-table rows for this species.
    species_patch_rows_before <- patch_table[species == species_name]

    # Skip this species if it has no surviving patch rows.
    if (nrow(species_patch_rows_before) == 0L) {
      next
    }

    # Keep only changed patches that still survive in this species' patch rows.
    changed_patch_ids_for_species <- intersect(
      changed_patch_ids_for_species,
      species_patch_rows_before$patch_id
    )

    # Skip this species if none of its changed patches survive.
    if (!length(changed_patch_ids_for_species)) {
      next
    }


    # -----------------------------------------------------------------
    # 7A. Initialize patch origin map for this species
    # -----------------------------------------------------------------

    # Build a table mapping current patches to their origin patch and PU.
    patch_origin_map <- species_patch_rows_before[
      ,
      list(
        patch_id = as.integer(patch_id),
        origin_patch_id = as.integer(patch_id),
        pu_id = as.integer(pu_id)
      )
    ]

    # Initialize new patch IDs above the current species-specific maximum.
    next_available_patch_id <- max(species_patch_rows_before$patch_id)

    # Read this species' patch -> cells index from the live environment.
    species_patch_cell_index_before <- get(
      species_name,
      envir = patch_cell_index_by_species_env,
      inherits = FALSE
    )

    if (is.null(species_patch_cell_index_before)) {
      stop(
        "Missing patch-cell index for species '",
        species_name,
        "' during fragmentation repair."
      )
    }

    # Accumulate newly created fragment-origin rows here.
    new_patch_origin_rows <- vector("list", 0L)


    # -----------------------------------------------------------------
    # 7B. Pass 1: relabel changed patches by rook contiguity
    # -----------------------------------------------------------------

    # Loop over changed patches for this species.
    for (original_patch_id in changed_patch_ids_for_species) {
      # Retrieve cells belonging to this patch at the start of fragmentation.
      patch_cells <- get_patch_cells(
        patch_index = species_patch_cell_index_before,
        patch_id = original_patch_id
      )

      # Skip if the patch has no indexed cells.
      if (!length(patch_cells)) {
        next
      }

      # Keep only cells still assigned to this original patch ID.
      patch_cells <- patch_cells[
        species_patch_ids_after[patch_cells] == original_patch_id
      ]

      # Skip if no current cells still carry this patch ID.
      if (!length(patch_cells)) {
        next
      }

      # A one-cell patch cannot split into multiple rook components.
      if (length(patch_cells) == 1L) {
        next
      }

      # Label rook-connected components within the current cells of this patch.
      component_result <- label_patch_cell_components(
        patch_cells = patch_cells,
        rook_neighbor_index = rook_neighbor_index
      )

      # Read component-labeled cells.
      patch_cells <- component_result$cells

      # Read component labels aligned with patch_cells.
      fragment_labels_in_patch_cells <- component_result$component_id

      unique_fragment_labels <- seq_len(component_result$component_count)

      if (length(patch_cells) != length(fragment_labels_in_patch_cells)) {
        stop(
          "Fragmentation component labeling returned mismatched lengths for species '",
          species_name,
          "', patch_id ",
          original_patch_id,
          ": length(patch_cells)=",
          length(patch_cells),
          ", length(fragment_labels_in_patch_cells)=",
          length(fragment_labels_in_patch_cells),
          "."
        )
      }

      # If the patch remains one connected component, no relabeling is needed.
      if (length(unique_fragment_labels) <= 1L) {
        next
      }

      # Sum exact cell area by fragment label.
      fragment_area_sum <- rowsum(
        cell_area_by_cell[patch_cells],
        fragment_labels_in_patch_cells,
        reorder = FALSE
      )

      # Build a fragment-area table.
      fragment_area_table <- data.table::data.table(
        fragment_label = as.integer(rownames(fragment_area_sum)),
        fragment_area_km2 = as.numeric(fragment_area_sum[, 1L])
      )

      # Choose the largest fragment to keep the original patch ID.
      keep_fragment_label <- fragment_area_table[
        order(-fragment_area_km2, fragment_label)
      ]$fragment_label[1L]

      inherited_pu_id <- patch_origin_map[
        patch_id == original_patch_id
      ]$pu_id[1L]

      if (is.na(inherited_pu_id)) {
        stop(
          "Could not find inherited PU id for species '",
          species_name,
          "', patch_id ",
          original_patch_id,
          " during fragmentation repair."
        )
      }

      # Assign new patch IDs to all non-largest fragments.
      for (fragment_label in unique_fragment_labels) {
        # Skip the largest fragment because it keeps the original patch ID.
        if (fragment_label == keep_fragment_label) {
          next
        }

        # Read cells belonging to this non-largest fragment.
        fragment_cells <- patch_cells[
          fragment_labels_in_patch_cells == fragment_label
        ]

        # Skip empty fragment labels defensively.
        if (!length(fragment_cells)) {
          next
        }

        # Advance the species-specific patch ID counter.
        next_available_patch_id <- next_available_patch_id + 1L

        # Relabel this fragment in the editable species vector.
        species_patch_ids_after[fragment_cells] <- next_available_patch_id

        new_patch_origin_rows[[length(new_patch_origin_rows) + 1L]] <-
          data.table::data.table(
            patch_id = as.integer(next_available_patch_id),
            origin_patch_id = as.integer(original_patch_id),
            pu_id = as.integer(inherited_pu_id)
          )
      }
    }

    # Append all newly created fragment-origin rows at once.
    if (length(new_patch_origin_rows)) {
      patch_origin_map <- data.table::rbindlist(
        c(
          list(patch_origin_map),
          new_patch_origin_rows
        ),
        use.names = TRUE,
        fill = TRUE
      )
    }


    # -----------------------------------------------------------------
    # 7C. Summarize patch areas after fragmentation relabeling
    # -----------------------------------------------------------------

    # Recompute current patch areas for this species from its updated cell vector.
    species_patch_rows_current <- summarize_species_patch_areas(
      species_name = species_name,
      species_patch_ids = species_patch_ids_after,
      cell_area_by_cell = cell_area_by_cell
    )

    # Attach origin patch IDs and PU IDs to current patch rows.
    species_patch_rows_current <- species_patch_rows_current[
      patch_origin_map,
      on = "patch_id",
      nomatch = 0L
    ]


    # -----------------------------------------------------------------
    # 7D. Pass 2: enforce minimum patch-size threshold
    # -----------------------------------------------------------------

    # Read this species' minimum patch-size threshold.
    patch_area_threshold <- species_params[
      species == species_name
    ]$min_patch_area_km2[1L]

    # Identify current patches below the patch-size threshold.
    dropped_patch_ids_by_fragmentation <- species_patch_rows_current[
      !area_exceeds_threshold(patch_area_km2, patch_area_threshold)
    ]$patch_id

    # Remove below-threshold patches from the editable species vector.
    if (length(dropped_patch_ids_by_fragmentation)) {
      species_patch_ids_after[
        species_patch_ids_after %in% dropped_patch_ids_by_fragmentation
      ] <- NA_integer_
    }

    # Recompute current patch areas after dropping below-threshold patches.
    species_patch_rows_current <- summarize_species_patch_areas(
      species_name = species_name,
      species_patch_ids = species_patch_ids_after,
      cell_area_by_cell = cell_area_by_cell
    )

    # Reattach origin patch IDs and PU IDs to surviving current patch rows.
    species_patch_rows_current <- species_patch_rows_current[
      patch_origin_map,
      on = "patch_id",
      nomatch = 0L
    ]


    # -----------------------------------------------------------------
    # 7E. Identify affected PUs
    # -----------------------------------------------------------------

    # Affected PUs are those containing changed origin patches.
    affected_pu_ids <- sort(unique(
      patch_origin_map$pu_id[
        patch_origin_map$origin_patch_id %in% changed_patch_ids_for_species
      ]
    ))

    # Keep existing species patch rows from PUs unaffected by fragmentation.
    unaffected_species_patch_rows <- species_patch_rows_before[
      !(pu_id %in% affected_pu_ids),
      list(
        species = species,
        patch_id = patch_id,
        pu_id = pu_id,
        patch_area_km2 = patch_area_km2
      )
    ]

    # Keep current patch rows from affected PUs only.
    affected_species_patch_rows_current <- species_patch_rows_current[
      pu_id %in% affected_pu_ids
    ]

    # Allocate surviving patch rows rebuilt for affected PUs.
    rebuilt_patch_rows_for_affected_pus <- vector("list", 0L)

    # Allocate rebuilt graph objects for affected PUs.
    rebuilt_graphs_for_affected_pus <- list()

    next_available_pu_id <- if (nrow(species_patch_rows_before)) {
      max(species_patch_rows_before$pu_id)
    } else {
      0L
    }

    pu_area_threshold <- species_params[
      species == species_name
    ]$min_population_area_km2[1L]


    # -----------------------------------------------------------------
    # 7F. Pass 3: rebuild only affected PUs
    # -----------------------------------------------------------------

    # Loop over affected original PUs.
    for (affected_pu_id in affected_pu_ids) {
      # Read current patch rows for this affected PU.
      current_rows_in_this_pu <- affected_species_patch_rows_current[
        pu_id == affected_pu_id
      ]

      # Build the old graph key.
      old_graph_key <- paste0(species_name, "|", affected_pu_id)

      old_pu_graph <- pu_graphs_by_key[[old_graph_key]]

      # If all patches in this old PU disappeared, there is no graph to rebuild.
      if (nrow(current_rows_in_this_pu) == 0L) {
        next
      }

      # A surviving PU must have a corresponding old graph.
      if (is.null(old_pu_graph)) {
        stop(
          "Missing PU graph during fragmentation-stage update for species '",
          species_name,
          "' and pu_id ",
          affected_pu_id,
          ". This indicates inconsistent evolving state between patch_table ",
          "and pu_graphs_by_key."
        )
      }

      # Read old graph patch IDs.
      old_patch_ids <- as.integer(old_pu_graph$id2patch)

      # Read current patch IDs in this affected PU.
      current_patch_ids <- as.integer(current_rows_in_this_pu$patch_id)

      # Deduplicate and sort old patch IDs.
      old_patch_ids_unique <- sort(unique(old_patch_ids))

      # Deduplicate and sort current patch IDs.
      current_patch_ids_unique <- sort(unique(current_patch_ids))

      # Identify current patches newly added by fragmentation.
      new_patch_nodes <- setdiff(current_patch_ids_unique, old_patch_ids_unique)

      # Identify old graph patches removed by fragmentation or thresholding.
      removed_old_patch_nodes <- setdiff(old_patch_ids_unique, current_patch_ids_unique)

      # Determine whether the PU graph node set changed.
      pu_patch_node_set_changed <-
        length(new_patch_nodes) > 0L ||
        length(removed_old_patch_nodes) > 0L

      # Compute current total area for this PU.
      current_pu_area_km2 <- sum(current_rows_in_this_pu$patch_area_km2)


      # ---------------------------------------------------------------
      # Fast path 1: no patch nodes were added or removed
      # ---------------------------------------------------------------

      # If the node set did not change, graph topology can be reused.
      if (!pu_patch_node_set_changed) {
        # If the PU is now below threshold, remove all of its current patches.
        if (!area_exceeds_threshold(current_pu_area_km2, pu_area_threshold)) {
          # Clear this PU's patches from the editable species vector.
          species_patch_ids_after[
            species_patch_ids_after %in% current_patch_ids_unique
          ] <- NA_integer_
        } else {
          # Reuse the old CSR graph unchanged.
          rebuilt_graphs_for_affected_pus[[old_graph_key]] <- list(
            species = species_name,
            pu_id = as.integer(affected_pu_id),
            id2patch = as.integer(old_pu_graph$id2patch),
            row_ptr = as.integer(old_pu_graph$row_ptr),
            col_idx = as.integer(old_pu_graph$col_idx)
          )

          # Copy current patch rows for this PU.
          current_rows_fast <- data.table::copy(current_rows_in_this_pu)

          # Preserve the affected PU ID.
          current_rows_fast[, pu_id := as.integer(affected_pu_id)]

          # Store current rows for final patch-table rebuild.
          rebuilt_patch_rows_for_affected_pus[[
            length(rebuilt_patch_rows_for_affected_pus) + 1L
          ]] <- current_rows_fast
        }

        # Continue to the next affected PU.
        next
      }


      # ---------------------------------------------------------------
      # Fast path 2: added nodes only, no old nodes removed
      # ---------------------------------------------------------------

      # Check whether fragmentation only added patch nodes.
      only_added_nodes_no_old_nodes_removed <-
        length(new_patch_nodes) > 0L &&
        length(removed_old_patch_nodes) == 0L

      # If only new nodes were added, the old graph cannot have disconnected.
      if (only_added_nodes_no_old_nodes_removed) {
        # If the PU is below threshold, remove all current patches in this PU.
        if (!area_exceeds_threshold(current_pu_area_km2, pu_area_threshold)) {
          # Clear this PU's patches from the editable species vector.
          species_patch_ids_after[
            species_patch_ids_after %in% current_patch_ids_unique
          ] <- NA_integer_
        } else {
          # Build a lazy overlay graph containing only edges for added nodes.
          lazy_overlay_graph <- make_lazy_added_node_overlay_pu_graph(
            current_pu_patch_rows = current_rows_in_this_pu,
            old_pu_graph = old_pu_graph,
            new_patch_nodes = new_patch_nodes
          )

          rebuilt_graphs_for_affected_pus[[old_graph_key]] <- lazy_overlay_graph

          # Copy current patch rows for this PU.
          current_rows_fast <- data.table::copy(current_rows_in_this_pu)

          # Preserve the affected PU ID.
          current_rows_fast[, pu_id := as.integer(affected_pu_id)]

          # Store current rows for final patch-table rebuild.
          rebuilt_patch_rows_for_affected_pus[[
            length(rebuilt_patch_rows_for_affected_pus) + 1L
          ]] <- current_rows_fast
        }

        # Continue to the next affected PU.
        next
      }


      # ---------------------------------------------------------------
      # Full path: old patch nodes were removed
      # ---------------------------------------------------------------

      # Build the full provisional graph for current patches in this PU.
      provisional_pu_graph <- build_provisional_pu_graph(
        current_pu_patch_rows = current_rows_in_this_pu,
        old_pu_graph = old_pu_graph
      )

      # Rebuild this PU after patch loss and possible graph disconnection.
      rebuilt_pu <- rebuild_pu_after_patch_loss(
        pu_graph = provisional_pu_graph,
        alive_nodes = rep(TRUE, nrow(current_rows_in_this_pu)),
        patch_area = current_rows_in_this_pu$patch_area_km2,
        pu_area_threshold = pu_area_threshold,
        next_available_pu_id = next_available_pu_id
      )

      # Carry forward the next available PU ID after possible splitting.
      next_available_pu_id <- rebuilt_pu$next_available_pu_id

      # If any components are too small, clear their patches from the species vector.
      if (length(rebuilt_pu$dropped_patch_ids)) {
        # Remove dropped component patches from the editable species vector.
        species_patch_ids_after[
          species_patch_ids_after %in% rebuilt_pu$dropped_patch_ids
        ] <- NA_integer_
      }

      # Store surviving PU graph objects.
      for (surviving_graph in rebuilt_pu$surviving_pu_graphs) {
        # Build graph key for this surviving PU.
        surviving_graph_key <- paste0(species_name, "|", surviving_graph$pu_id)

        rebuilt_graphs_for_affected_pus[[surviving_graph_key]] <- list(
          species = species_name,
          pu_id = as.integer(surviving_graph$pu_id),
          id2patch = as.integer(surviving_graph$id2patch),
          row_ptr = as.integer(surviving_graph$row_ptr),
          col_idx = as.integer(surviving_graph$col_idx)
        )
      }

      # Store patch rows for surviving PU components with updated PU IDs.
      for (surviving_graph in rebuilt_pu$surviving_pu_graphs) {
        # Read patch IDs in this surviving component.
        component_patch_ids <- as.integer(surviving_graph$id2patch)

        # Extract current rows for these component patches.
        component_rows <- current_rows_in_this_pu[
          patch_id %in% component_patch_ids
        ]

        # Store component rows if any exist.
        if (nrow(component_rows)) {
          # Copy component rows before mutating pu_id.
          component_rows <- data.table::copy(component_rows)

          # Assign the surviving component PU ID.
          component_rows[, pu_id := as.integer(surviving_graph$pu_id)]

          # Store component rows for final patch-table rebuild.
          rebuilt_patch_rows_for_affected_pus[[
            length(rebuilt_patch_rows_for_affected_pus) + 1L
          ]] <- component_rows
        }
      }
    }


    # -----------------------------------------------------------------
    # 7G. Combine unaffected and rebuilt patch rows for this species
    # -----------------------------------------------------------------

    # Combine rebuilt affected-PU patch rows, or create an empty table.
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
        patch_area_km2 = numeric(),
        origin_patch_id = integer(),
        pu_id = integer()
      )
    }

    # Combine unaffected rows and rebuilt affected rows, retaining origin IDs.
    final_species_patch_rows_with_origin <- data.table::rbindlist(
      list(
        unaffected_species_patch_rows[
          ,
          list(
            species = species,
            patch_id = patch_id,
            pu_id = pu_id,
            patch_area_km2 = patch_area_km2,
            origin_patch_id = patch_id
          )
        ],
        rebuilt_species_patch_rows[
          ,
          list(
            species = species,
            patch_id = patch_id,
            pu_id = pu_id,
            patch_area_km2 = patch_area_km2,
            origin_patch_id = origin_patch_id
          )
        ]
      ),
      use.names = TRUE,
      fill = TRUE
    )


    # -----------------------------------------------------------------
    # 7H. Replace this species' patch rows in the global patch table
    # -----------------------------------------------------------------

    # Remove old patch rows for this species from the global patch table.
    patch_table <- patch_table[species != species_name]

    # Append updated patch rows for this species.
    patch_table <- data.table::rbindlist(
      list(
        patch_table,
        final_species_patch_rows_with_origin[
          ,
          list(
            species = species,
            patch_id = patch_id,
            pu_id = pu_id,
            patch_area_km2 = patch_area_km2
          )
        ]
      ),
      use.names = TRUE,
      fill = TRUE
    )


    # -----------------------------------------------------------------
    # 7I. Replace affected PU graphs in the global graph list
    # -----------------------------------------------------------------

    # Remove old graph entries for affected original PUs.
    for (affected_pu_id in affected_pu_ids) {
      # Delete the old graph entry for this affected PU.
      pu_graphs_by_key[[paste0(species_name, "|", affected_pu_id)]] <- NULL
    }

    # Add rebuilt graph entries for affected surviving PUs.
    for (graph_key in names(rebuilt_graphs_for_affected_pus)) {
      pu_graphs_by_key[[graph_key]] <- rebuilt_graphs_for_affected_pus[[graph_key]]
    }


    # -----------------------------------------------------------------
    # 7J. Update alive-species counts for cells where this species vanished
    # -----------------------------------------------------------------

    # Identify cells where this species was present before but absent after fragmentation.
    species_cells_removed_in_fragmentation <- which(
      !is.na(species_patch_ids_before) &
        is.na(species_patch_ids_after)
    )

    # Read global alive-species counts before removing this species from those cells.
    count_before_species_loss <- alive_species_count_by_cell[
      species_cells_removed_in_fragmentation
    ]

    # Identify cells that remain globally alive after losing this species.
    cells_still_alive_after_species_loss <- species_cells_removed_in_fragmentation[
      count_before_species_loss > 0L
    ]

    # Identify cells that become globally empty after losing this species.
    cells_becoming_empty <- species_cells_removed_in_fragmentation[
      count_before_species_loss == 1L
    ]

    # Decrement alive-species counts for cells that were still globally alive.
    if (length(cells_still_alive_after_species_loss)) {
      alive_species_count_by_cell[cells_still_alive_after_species_loss] <-
        alive_species_count_by_cell[cells_still_alive_after_species_loss] - 1L
    }

    # Add newly empty cells to the fragmentation-stage accumulator.
    if (length(cells_becoming_empty)) {
      removed_cells_in_fragmentation <- c(
        removed_cells_in_fragmentation,
        cells_becoming_empty
      )
    }


    # -----------------------------------------------------------------
    # 7K. Write updated species state back to live environments
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


    # -----------------------------------------------------------------
    # 7L. Record surviving patches requiring distance recheck
    # -----------------------------------------------------------------

    # Store patches whose geometry or identity changed and survived this stage.
    recheck_patch_rows_by_species[[species_index]] <- unique(
      final_species_patch_rows_with_origin[
        origin_patch_id %in% changed_patch_ids_for_species |
          patch_id != origin_patch_id,
        list(
          species = species,
          patch_id = patch_id
        )
      ],
      by = c("species", "patch_id")
    )
  }


  # -------------------------------------------------------------------
  # 8. Finalize removed cells and distance-recheck patch table
  # -------------------------------------------------------------------

  # Deduplicate and sort cells that became globally empty in fragmentation.
  removed_cells_in_fragmentation <- sort(
    unique(as.integer(removed_cells_in_fragmentation))
  )

  # Keep non-empty distance-recheck tables.
  non_empty_recheck_tables <- Filter(
    f = function(x) is.data.frame(x) && nrow(x) > 0L,
    x = recheck_patch_rows_by_species
  )

  # Combine distance-recheck tables, or return an empty table.
  patches_requiring_distance_recheck <- if (length(non_empty_recheck_tables)) {
    # Combine and deduplicate recheck rows.
    unique(
      data.table::rbindlist(non_empty_recheck_tables, use.names = TRUE, fill = TRUE),
      by = c("species", "patch_id")
    )
  } else {
    # Return an empty recheck table.
    data.table::data.table(
      species = character(),
      patch_id = integer()
    )
  }


  # -------------------------------------------------------------------
  # 9. Log one summary line for this fragmentation stage
  # -------------------------------------------------------------------

  # Count remaining patch-table rows.
  remaining_patch_count <- nrow(patch_table)

  # Count remaining unique species-PU combinations.
  remaining_pu_count <- data.table::uniqueN(patch_table, by = c("species", "pu_id"))

  # Count cells still containing at least one surviving species.
  remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

  log_fragmentation_stage(
    stage_index = stage_index,
    touched_species_count = length(touched_species),
    input_changed_patch_count = input_changed_patch_count,
    removed_cell_count = length(removed_cells_in_fragmentation),
    recheck_patch_count = nrow(patches_requiring_distance_recheck),
    remaining_patch_count = remaining_patch_count,
    remaining_pu_count = remaining_pu_count,
    remaining_alive_cell_count = remaining_alive_cell_count
  )


  # -------------------------------------------------------------------
  # 10. Return updated live state
  # -------------------------------------------------------------------

  # Return the updated state needed by the outer pipeline and distance stage.
  list(
    patch_table = patch_table,
    pu_graphs = pu_graphs_by_key,
    alive_species_count_by_cell = alive_species_count_by_cell,
    removed_cells_in_fragmentation = removed_cells_in_fragmentation,
    patches_requiring_distance_recheck = patches_requiring_distance_recheck,
    touched_species = touched_species
  )
}
