source(file.path("R", "project_utils.R"))
source(file.path("R", "priority_outputs.R"))

# Stage-wise reverse-removal priority pipeline.


# =====================================================================
# HELPER: write the patch lookup table for one completed global stage
# =====================================================================

write_stage_patch_lookup_table <- function(
  patch_table,                 # current patch table at the end of the stage
  stage_index,                 # completed global stage index
  patch_lookup_output_dir      # directory where stage patch tables are written
) {
  # Build the output path using zero-padded stage numbering.
  output_path <- file.path(
    patch_lookup_output_dir,
    sprintf("stage_patch_lookup_stage_%04d.csv", as.integer(stage_index))
  )

  data.table::fwrite(
    x = patch_table,
    file = output_path
  )

  # Return the output path invisibly for optional caller use.
  invisible(output_path)
}


# =====================================================================
# HELPER: print one end-of-stage pipeline summary message
# =====================================================================

log_priority_pipeline_stage <- function(
  stage_index,                 # completed global stage index
  frontier_exhausted,          # whether pruning exhausted the frontier
  removal_step_count,          # number of removal-order events recorded so far
  remaining_patch_count,       # number of patch-table rows remaining
  remaining_pu_count,          # number of species-PU combinations remaining
  remaining_alive_cell_count   # number of cells with at least one species remaining
) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  message(
    sprintf(
      paste(
        "[%s] priority_pipeline |",
        "stage=%d frontier_exhausted=%s removal_steps=%d",
        "remaining_patches=%d remaining_pus=%d remaining_alive_cells=%d"
      ),
      timestamp,
      as.integer(stage_index),
      if (isTRUE(frontier_exhausted)) "TRUE" else "FALSE",
      as.integer(removal_step_count),
      as.integer(remaining_patch_count),
      as.integer(remaining_pu_count),
      as.integer(remaining_alive_cell_count)
    )
  )

  invisible(NULL)
}


# =====================================================================
# HELPER: print one boundary-start pipeline message
# =====================================================================

log_priority_pipeline_boundary_start <- function(
  stage_index,                 # global stage index about to enter a boundary step
  boundary_step,               # boundary step label, e.g. fragmentation_start
  workload_label,              # name of the input workload count
  workload_count               # size of the input workload
) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  message(
    sprintf(
      "[%s] priority_pipeline_boundary | stage=%d step=%s %s=%d",
      timestamp,
      as.integer(stage_index),
      as.character(boundary_step),
      as.character(workload_label),
      as.integer(workload_count)
    )
  )

  invisible(NULL)
}


# =====================================================================
# MAIN FUNCTION: run the full stage-wise prioritization pipeline
# =====================================================================

run_priority_pipeline <- function(
  cells_to_remove_per_iteration,    # frontier cells removed per pruning iteration
  pruning_iterations_per_stage,     # maximum pruning iterations per global stage
  patch_table,                      # current patch table
  pu_graphs_by_key,                 # current PU graphs keyed by "species|pu_id"
  alive_species_count_by_cell,      # current surviving-species count by raster cell
  patch_id_by_species_env,          # live env: species -> cell-to-patch vector
  patch_cell_index_by_species_env,  # live env: species -> patch-to-cells index
  cell_area_by_cell,                # area of each raster cell in km^2
  rook_neighbor_pairs,              # rook-neighbor cell pairs for frontier tracking
  species_params,                   # species thresholds, persistence params, dispersal distance
  species_raster_template_stack,    # raster templates needed by the distance stage
  mask_template_raster,             # single-layer grid template for removal_order.tif
  output_dir,                       # root output directory
  max_stages = Inf                  # maximum number of full global stages
) {
  # -------------------------------------------------------------------
  # 1. Validate scalar control inputs
  # -------------------------------------------------------------------

  if (!is.numeric(cells_to_remove_per_iteration) ||
      length(cells_to_remove_per_iteration) != 1L ||
      is.na(cells_to_remove_per_iteration) ||
      cells_to_remove_per_iteration < 1) {
    stop("cells_to_remove_per_iteration must be a single integer >= 1.")
  }

  if (!is.numeric(pruning_iterations_per_stage) ||
      length(pruning_iterations_per_stage) != 1L ||
      is.na(pruning_iterations_per_stage) ||
      pruning_iterations_per_stage < 1) {
    stop("pruning_iterations_per_stage must be a single integer >= 1.")
  }

  if (!is.numeric(max_stages) ||
      length(max_stages) != 1L ||
      is.na(max_stages) ||
      max_stages <= 0) {
    stop("max_stages must be a single positive number or Inf.")
  }
  
  cells_to_remove_per_iteration <- as.integer(cells_to_remove_per_iteration)

  pruning_iterations_per_stage <- as.integer(pruning_iterations_per_stage)

  max_stages <- if (is.finite(max_stages)) as.integer(max_stages) else Inf

  # -------------------------------------------------------------------
  # 2. Standardize table inputs
  # -------------------------------------------------------------------

  patch_table <- data.table::as.data.table(patch_table)

  species_params <- data.table::as.data.table(species_params)


  # -------------------------------------------------------------------
  # 3. Validate required columns and template objects
  # -------------------------------------------------------------------

  if (!("dispersal_distance_km" %in% names(species_params))) {
    stop(
      "species_params must contain the canonical column 'dispersal_distance_km'."
    )
  }

  if (!is.character(output_dir) || length(output_dir) != 1L || !nzchar(output_dir)) {
    stop("output_dir must be one non-empty character path.")
  }

  if (is.null(mask_template_raster)) {
    stop("mask_template_raster must be supplied for writing removal_order.tif.")
  }

  if (is.null(species_raster_template_stack)) {
    stop("species_raster_template_stack must be supplied for the distance stage.")
  }


  # -------------------------------------------------------------------
  # 4. Create the minimal production output directories
  # -------------------------------------------------------------------

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  patch_lookup_output_dir <- file.path(output_dir, "patch_lookup_tables")

  dir.create(patch_lookup_output_dir, recursive = TRUE, showWarnings = FALSE)

  removal_order_output_path <- file.path(output_dir, "removal_order.tif")

  rankmap_output_path <- file.path(output_dir, "rankmap.tif")
  
  removal_events_output_path <- file.path(output_dir, "removal_events.csv")


  # -------------------------------------------------------------------
  # 5. Initialize global stage-loop state
  # -------------------------------------------------------------------

  n_cells <- length(alive_species_count_by_cell)
  
  # Record which cells were alive at the start of the prioritization run.
  initial_alive_by_cell <- alive_species_count_by_cell > 0L
  
  # Count initially alive cells once.
  initial_alive_cell_count <- sum(initial_alive_by_cell)
  
  # Sum initial alive habitat area once.
  initial_alive_area_km2 <- sum(
    cell_area_by_cell[initial_alive_by_cell],
    na.rm = TRUE
  )
  
  # Allocate the removal-order vector; NA means not removed or not initially alive.
  removal_order_by_cell <- rep(NA_integer_, n_cells)
  
  removal_step_counter <- 1L
  
  # Store one small metadata row per non-empty removal event.
  removal_event_rows <- vector("list", 0L)

  # Count completed full global stages.
  completed_stages <- 0L

  # Track whether pruning has exhausted the frontier.
  frontier_exhausted <- FALSE


  # -------------------------------------------------------------------
  # 6. Define local helper to record removal order in memory
  # -------------------------------------------------------------------
  
  record_removed_cells <- function(
    removed_cells,
    stage_index,
    event_type,
    pruning_iteration = NA_integer_
  ) {
    # Return immediately if this event removed no cells.
    if (!length(removed_cells)) {
      return(invisible(0L))
    }
  
    removed_cells <- sort(unique(as.integer(removed_cells)))
  
    # Keep only valid raster-cell IDs.
    removed_cells <- removed_cells[
      !is.na(removed_cells) &
        removed_cells >= 1L &
        removed_cells <= n_cells
    ]
  
    # Return if no valid cell IDs remain.
    if (!length(removed_cells)) {
      return(invisible(0L))
    }
  
    # Keep only cells that:
    #   1. were alive at the start of the run, and
    #   2. have not already been assigned a removal order.
    eligible_cells <- removed_cells[
      initial_alive_by_cell[removed_cells] &
        is.na(removal_order_by_cell[removed_cells])
    ]
  
    # Return if this event contains no newly recordable cells.
    if (!length(eligible_cells)) {
      return(invisible(0L))
    }
  
    this_step <- as.integer(removal_step_counter)
  
    # Assign this removal step to all eligible cells in the removal-order raster vector.
    removal_order_by_cell[eligible_cells] <<- this_step
  
    # Compute how many cells were actually recorded in this event.
    cells_removed_this_event <- length(eligible_cells)
  
    # Compute the area removed in this event.
    area_removed_this_event_km2 <- sum(
      cell_area_by_cell[eligible_cells],
      na.rm = TRUE
    )
  
    # Append one metadata row for this removal event.
    removal_event_rows[[length(removal_event_rows) + 1L]] <<-
      data.table::data.table(
        removal_step = this_step,
        stage = as.integer(stage_index),
        event_type = as.character(event_type),
        pruning_iteration = as.integer(pruning_iteration),
        cells_removed = as.integer(cells_removed_this_event),
        area_removed_km2 = as.numeric(area_removed_this_event_km2)
      )
  
    # Advance the removal-event counter once per non-empty recorded event.
    removal_step_counter <<- removal_step_counter + 1L
  
    # Return the number of cells recorded in this event.
    invisible(cells_removed_this_event)
  }
  
  
  
  record_final_retained_cells <- function(stage_index) {
    # Identify initially alive cells that never received a removal-order value.
    final_retained_cells <- which(
      initial_alive_by_cell &
        is.na(removal_order_by_cell)
    )
  
    # Assign one final removal-order step to all retained-to-end cells. If all
    # cells have already been removed, keep a zero-cell terminal event so the
    # event table still has an explicit endpoint.
    this_step <- as.integer(removal_step_counter)
  
    if (length(final_retained_cells)) {
      removal_order_by_cell[final_retained_cells] <<- this_step
    }
  
    # Compute retained-to-end area.
    final_retained_area_km2 <- sum(
      cell_area_by_cell[final_retained_cells],
      na.rm = TRUE
    )
  
    # Add a final event row so removal_order.tif and removal_events.csv stay consistent.
    removal_event_rows[[length(removal_event_rows) + 1L]] <<-
      data.table::data.table(
        removal_step = this_step,
        stage = as.integer(stage_index),
        event_type = "final_retained",
        pruning_iteration = NA_integer_,
        cells_removed = as.integer(length(final_retained_cells)),
        area_removed_km2 = as.numeric(final_retained_area_km2)
      )
  
    # Advance the event counter.
    removal_step_counter <<- removal_step_counter + 1L
  
    invisible(length(final_retained_cells))
  }


  # -------------------------------------------------------------------
  # 7. Build frontier-tracking objects once
  # -------------------------------------------------------------------

  # Build the compact rook-neighbor index once from the input neighbor pairs.
  rook_neighbor_index <- build_rook_neighbor_index(
    rook_neighbor_pairs = rook_neighbor_pairs,
    n_cells = n_cells,
    force_symmetric = TRUE
  )

  frontier_state <- initialize_frontier_state(
    alive_species_count_by_cell = alive_species_count_by_cell,
    rook_neighbor_index = rook_neighbor_index
  )


  # -------------------------------------------------------------------
  # 8. Run the global stage loop
  # -------------------------------------------------------------------

  stage_index <- 0L

  # Run until the frontier is exhausted or max_stages is reached.
  while (!isTRUE(frontier_exhausted) && stage_index < max_stages) {
    # Advance to the next global stage.
    stage_index <- stage_index + 1L


    # ---------------------------------------------------------------
    # 8A. Run the pruning stage
    # ---------------------------------------------------------------

    # Run repeated pruning iterations for this global stage.
    pruning_stage_result <- run_pruning_stage(
      stage_index = stage_index,
      pruning_iterations_per_stage = pruning_iterations_per_stage,
      cells_to_remove_per_iteration = cells_to_remove_per_iteration,
      patch_table = patch_table,
      pu_graphs_by_key = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      patch_id_by_species_env = patch_id_by_species_env,
      patch_cell_index_by_species_env = patch_cell_index_by_species_env,
      cell_area_by_cell = cell_area_by_cell,
      frontier_state = frontier_state,
      rook_neighbor_index = rook_neighbor_index,
      species_params = species_params
    )

    patch_table <- pruning_stage_result$patch_table

    pu_graphs_by_key <- pruning_stage_result$pu_graphs

    alive_species_count_by_cell <- pruning_stage_result$alive_species_count_by_cell

    frontier_state <- pruning_stage_result$frontier_state

    # Record whether pruning exhausted the frontier in this stage.
    frontier_exhausted <- isTRUE(pruning_stage_result$frontier_exhausted)
    
    # Record removal order for each pruning iteration in chronological order.
    if (length(pruning_stage_result$removed_cells_each_iteration)) {
      for (iteration_index in seq_along(pruning_stage_result$removed_cells_each_iteration)) {
        record_removed_cells(
          removed_cells = pruning_stage_result$removed_cells_each_iteration[[iteration_index]],
          stage_index = stage_index,
          event_type = "prune",
          pruning_iteration = iteration_index
        )
      }
    }


    # ---------------------------------------------------------------
    # 8B. Run the fragmentation stage
    # ---------------------------------------------------------------

    log_priority_pipeline_boundary_start(
      stage_index = stage_index,
      boundary_step = "fragmentation_start",
      workload_label = "input_changed_patches",
      workload_count = nrow(pruning_stage_result$changed_patches_in_stage)
    )

    # Run fragmentation repair on patches changed by pruning.
    fragmentation_stage_result <- run_fragmentation_stage(
      stage_index = stage_index,
      changed_patches_in_stage = pruning_stage_result$changed_patches_in_stage,
      patch_table = patch_table,
      pu_graphs_by_key = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      patch_id_by_species_env = patch_id_by_species_env,
      patch_cell_index_by_species_env = patch_cell_index_by_species_env,
      species_params = species_params,
      cell_area_by_cell = cell_area_by_cell,
      rook_neighbor_index = rook_neighbor_index
    )

    patch_table <- fragmentation_stage_result$patch_table

    pu_graphs_by_key <- fragmentation_stage_result$pu_graphs

    alive_species_count_by_cell <- fragmentation_stage_result$alive_species_count_by_cell

    frontier_state <- update_frontier_state_after_empty_cells(
      frontier_state = frontier_state,
      newly_empty_cells = fragmentation_stage_result$removed_cells_in_fragmentation,
      rook_neighbor_index = rook_neighbor_index
    )

    # Record fragmentation removals as one removal-order event.
    record_removed_cells(
      removed_cells = fragmentation_stage_result$removed_cells_in_fragmentation,
      stage_index = stage_index,
      event_type = "fragmentation",
      pruning_iteration = NA_integer_
    )


    # ---------------------------------------------------------------
    # 8C. Run the distance-connectivity stage
    # ---------------------------------------------------------------

    log_priority_pipeline_boundary_start(
      stage_index = stage_index,
      boundary_step = "distance_start",
      workload_label = "input_recheck_patches",
      workload_count = nrow(fragmentation_stage_result$patches_requiring_distance_recheck)
    )

    # Run distance-based connectivity repair on patches flagged by fragmentation.
    distance_stage_result <- run_distance_connectivity_stage(
      stage_index = stage_index,
      patches_requiring_distance_recheck =
        fragmentation_stage_result$patches_requiring_distance_recheck,
      patch_table = patch_table,
      pu_graphs_by_key = pu_graphs_by_key,
      alive_species_count_by_cell = alive_species_count_by_cell,
      patch_id_by_species_env = patch_id_by_species_env,
      patch_cell_index_by_species_env = patch_cell_index_by_species_env,
      species_params = species_params,
      species_raster_template_stack = species_raster_template_stack
    )

    patch_table <- distance_stage_result$patch_table

    pu_graphs_by_key <- distance_stage_result$pu_graphs

    alive_species_count_by_cell <- distance_stage_result$alive_species_count_by_cell

    frontier_state <- update_frontier_state_after_empty_cells(
      frontier_state = frontier_state,
      newly_empty_cells = distance_stage_result$removed_cells_in_distance,
      rook_neighbor_index = rook_neighbor_index
    )

    # Record distance-stage removals as one removal-order event.
    record_removed_cells(
      removed_cells = distance_stage_result$removed_cells_in_distance,
      stage_index = stage_index,
      event_type = "distance",
      pruning_iteration = NA_integer_
    )


    # ---------------------------------------------------------------
    # 8D. Write the patch lookup table for this completed full stage
    # ---------------------------------------------------------------

    write_stage_patch_lookup_table(
      patch_table = patch_table,
      stage_index = stage_index,
      patch_lookup_output_dir = patch_lookup_output_dir
    )


    # ---------------------------------------------------------------
    # 8E. Print one end-of-stage pipeline summary
    # ---------------------------------------------------------------

    # Count remaining patch-table rows.
    remaining_patch_count <- nrow(patch_table)

    # Count remaining unique species-PU combinations.
    remaining_pu_count <- data.table::uniqueN(
      patch_table,
      by = c("species", "pu_id")
    )

    # Count remaining globally alive cells.
    remaining_alive_cell_count <- sum(alive_species_count_by_cell > 0L)

    log_priority_pipeline_stage(
      stage_index = stage_index,
      frontier_exhausted = frontier_exhausted,
      removal_step_count = removal_step_counter - 1L,
      remaining_patch_count = remaining_patch_count,
      remaining_pu_count = remaining_pu_count,
      remaining_alive_cell_count = remaining_alive_cell_count
    )


    # ---------------------------------------------------------------
    # 8F. Record completed stage count
    # ---------------------------------------------------------------

    completed_stages <- as.integer(stage_index)
  }
  
  
  
  # -------------------------------------------------------------------
  # 9. Assign the terminal retained layer after the stopping condition
  # -------------------------------------------------------------------
  
  # Zonation rank maps assign every analysis-domain cell a priority rank.
  # Therefore, regardless of whether this run stopped because the frontier was
  # exhausted or because max_stages was reached, assign every still-unordered
  # initially alive cell to one terminal removal-order layer. Use a synthetic
  # stage after the last completed stage so stage_meta rows remain aligned with
  # actual patch_lookup_tables/stage_patch_lookup_stage_####.csv files.
  final_retained_cell_count <- 0L

  final_retained_cell_count <- record_final_retained_cells(
    stage_index = as.integer(completed_stages + 1L)
  )


  # -------------------------------------------------------------------
  # 10. Write the final removal-order raster once
  # -------------------------------------------------------------------

  write_removal_order_raster(
    removal_order_by_cell = removal_order_by_cell,
    template_raster = mask_template_raster,
    output_path = removal_order_output_path
  )
  
  
  
  # -------------------------------------------------------------------
  # -------------------------------------------------------------------
  
  # Combine one-row event tables, or create an empty event table if no cells were removed.
  removal_events <- if (length(removal_event_rows)) {
    data.table::rbindlist(
      removal_event_rows,
      use.names = TRUE,
      fill = TRUE
    )
  } else {
    data.table::data.table(
      removal_step = integer(),
      stage = integer(),
      event_type = character(),
      pruning_iteration = integer(),
      cells_removed = integer(),
      area_removed_km2 = numeric()
    )
  }
  
  # Add cumulative summaries if at least one removal event occurred.
  if (nrow(removal_events)) {
    data.table::setorder(removal_events, removal_step)
  
    # Cumulative number of initially alive cells removed.
    removal_events[, cum_cells_removed := cumsum(cells_removed)]
  
    # Cumulative area removed.
    removal_events[, cum_area_removed_km2 := cumsum(area_removed_km2)]
  
    # Cumulative proportion of initially alive cells removed.
    removal_events[, cum_prop_cells_removed :=
      cum_cells_removed / as.numeric(initial_alive_cell_count)
    ]
  
    # Cumulative proportion of initially alive area removed.
    removal_events[, cum_prop_area_removed :=
      cum_area_removed_km2 / as.numeric(initial_alive_area_km2)
    ]
  
    # Remaining initially alive cells after this event.
    removal_events[, cells_retained :=
      as.integer(initial_alive_cell_count - cum_cells_removed)
    ]
  
    # Remaining initially alive area after this event.
    removal_events[, area_retained_km2 :=
      as.numeric(initial_alive_area_km2 - cum_area_removed_km2)
    ]
  } else {
    # Add the same columns to the empty table so downstream code has a stable schema.
    removal_events[, `:=`(
      cum_cells_removed = integer(),
      cum_area_removed_km2 = numeric(),
      cum_prop_cells_removed = numeric(),
      cum_prop_area_removed = numeric(),
      cells_retained = integer(),
      area_retained_km2 = numeric()
    )]
  }
  
  data.table::fwrite(
    x = removal_events,
    file = removal_events_output_path
  )
  
  # proportions. The final retained layer maps to 1. Cells outside the initial
  # alive domain remain NA.
  write_rankmap_raster(
    removal_order_by_cell = removal_order_by_cell,
    removal_events = removal_events,
    template_raster = mask_template_raster,
    output_path = rankmap_output_path,
    value_col = "cum_prop_cells_removed"
  )


  # -------------------------------------------------------------------
  # 11. Compute compact final-run summaries
  # -------------------------------------------------------------------

  # Count cells that were alive at the start of the run.
  initial_alive_cell_count <- sum(initial_alive_by_cell)

  # Count initially alive cells that received a removal order.
  removed_initial_alive_cell_count <- sum(
    initial_alive_by_cell &
      !is.na(removal_order_by_cell)
  )

  # Count initially alive cells that still have no removal order.
  unremoved_initial_alive_cell_count <- sum(
    initial_alive_by_cell &
      is.na(removal_order_by_cell)
  )


  # -------------------------------------------------------------------
  # 12. Return the final in-memory state and output paths
  # -------------------------------------------------------------------

  # Return the final state and compact output metadata.
  list(
    patch_table = patch_table,
    pu_graphs = pu_graphs_by_key,
    alive_species_count_by_cell = alive_species_count_by_cell,
    frontier_state = frontier_state,
    completed_stages = completed_stages,
    frontier_exhausted = frontier_exhausted,
    removal_order_path = removal_order_output_path,
    rankmap_path = rankmap_output_path,
    removal_events_path = removal_events_output_path,
    patch_lookup_output_dir = patch_lookup_output_dir,
    removal_step_count = removal_step_counter - 1L,
    initial_alive_cell_count = initial_alive_cell_count,
    removed_initial_alive_cell_count = removed_initial_alive_cell_count,
    unremoved_initial_alive_cell_count = unremoved_initial_alive_cell_count
  )
}
