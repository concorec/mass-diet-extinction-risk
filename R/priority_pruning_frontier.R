# =====================================================================
# HELPER: get cells for many patches from one species patch index
# =====================================================================

get_patch_cells_batch <- function(patch_index, patch_ids) {
  # Convert requested patch IDs to sorted unique integers.
  patch_ids <- sort(unique(as.integer(patch_ids)))

  # Keep only valid positive patch IDs.
  patch_ids <- patch_ids[
    !is.na(patch_ids) &
      patch_ids >= 1L
  ]

  # Return an empty lookup table if no valid patch IDs were requested.
  if (!length(patch_ids)) {
    return(data.table::data.table(
      patch_id = integer(),
      cell = integer()
    ))
  }

  # Return an empty lookup table if the species patch-cell index is missing.
  if (is.null(patch_index) ||
      is.null(patch_index$pid) ||
      is.null(patch_index$cell) ||
      !length(patch_index$pid)) {
    return(data.table::data.table(
      patch_id = integer(),
      cell = integer()
    ))
  }

  # Find the first index-row position for each requested patch ID.
  lo <- findInterval(patch_ids - 0.5, patch_index$pid) + 1L

  # Find the last index-row position for each requested patch ID.
  hi <- findInterval(patch_ids + 0.5, patch_index$pid)

  # Mark requested patch IDs that are present in the sorted patch index.
  patch_present <- lo <= hi

  # Return an empty lookup table if none of the requested patch IDs are indexed.
  if (!any(patch_present)) {
    return(data.table::data.table(
      patch_id = integer(),
      cell = integer()
    ))
  }

  # Keep only requested patch IDs that are actually present.
  patch_ids <- patch_ids[patch_present]

  # Keep the corresponding first row positions.
  lo <- lo[patch_present]

  # Keep the corresponding last row positions.
  hi <- hi[patch_present]

  # Count how many cells belong to each requested patch.
  n_cells_by_patch <- as.integer(hi - lo + 1L)

  # Build all row positions in patch_index$cell for all requested patches at once.
  patch_index_rows <- sequence(n_cells_by_patch) +
    rep.int(lo - 1L, n_cells_by_patch)

  # Return one row per patch-cell membership.
  data.table::data.table(
    patch_id = rep.int(patch_ids, n_cells_by_patch),
    cell = as.integer(patch_index$cell[patch_index_rows])
  )
}



# =====================================================================
# HELPER: get all cells that currently belong to one patch
# =====================================================================

get_patch_cells <- function(patch_index, patch_id) {
  lo <- findInterval(patch_id - 0.5, patch_index$pid) + 1L
  hi <- findInterval(patch_id + 0.5, patch_index$pid)

  if (lo <= hi) {
    patch_index$cell[lo:hi]
  } else {
    integer(0L)
  }
}



# =====================================================================
# HELPER: connected components of one patch using cell indices
# =====================================================================

label_patch_cell_components <- function(patch_cells, rook_neighbor_index) {
  patch_cells <- sort(unique(as.integer(patch_cells)))

  patch_cells <- patch_cells[
    !is.na(patch_cells) &
      patch_cells >= 1L &
      patch_cells <= rook_neighbor_index$n_cells
  ]

  n_patch_cells <- length(patch_cells)

  if (!n_patch_cells) {
    return(list(
      cells = integer(0L),
      component_id = integer(0L),
      component_count = 0L
    ))
  }

  if (n_patch_cells == 1L) {
    return(list(
      cells = patch_cells,
      component_id = 1L,
      component_count = 1L
    ))
  }

  component_id_by_position <- integer(n_patch_cells)
  component_count <- 0L
  queue <- integer(n_patch_cells)

  for (seed_position in seq_len(n_patch_cells)) {
    if (component_id_by_position[seed_position] != 0L) {
      next
    }

    component_count <- component_count + 1L

    queue_head <- 1L
    queue_tail <- 1L
    queue[queue_tail] <- seed_position

    component_id_by_position[seed_position] <- component_count

    while (queue_head <= queue_tail) {
      current_position <- queue[queue_head]
      queue_head <- queue_head + 1L

      current_cell <- patch_cells[current_position]

      start_pos <- rook_neighbor_index$start[current_cell]

      if (is.na(start_pos)) {
        next
      }

      end_pos <- rook_neighbor_index$end[current_cell]
      neighbor_cells <- rook_neighbor_index$to[start_pos:end_pos]

      neighbor_positions <- fastmatch::fmatch(
        neighbor_cells,
        patch_cells,
        nomatch = 0L
      )

      neighbor_positions <- neighbor_positions[neighbor_positions > 0L]

      if (!length(neighbor_positions)) {
        next
      }

      neighbor_positions <- unique(neighbor_positions)

      neighbor_positions <- neighbor_positions[
        component_id_by_position[neighbor_positions] == 0L
      ]

      if (!length(neighbor_positions)) {
        next
      }

      component_id_by_position[neighbor_positions] <- component_count

      n_new <- length(neighbor_positions)

      queue[(queue_tail + 1L):(queue_tail + n_new)] <- neighbor_positions
      queue_tail <- queue_tail + n_new
    }
  }

  list(
    cells = patch_cells,
    component_id = as.integer(component_id_by_position),
    component_count = as.integer(component_count)
  )
}



# =====================================================================
# HELPER: build compact rook-neighbor index
# =====================================================================

build_rook_neighbor_index <- function(rook_neighbor_pairs, n_cells, force_symmetric = TRUE) {
  # Stop if the neighbor-pair object does not have at least two columns.
  if (is.null(rook_neighbor_pairs) || ncol(rook_neighbor_pairs) < 2L) {
    stop("rook_neighbor_pairs must have at least two columns.")
  }

  # Store the total cell count as an integer.
  n_cells <- as.integer(n_cells)

  # Read the first neighbor-pair column as source cell IDs.
  from <- as.integer(rook_neighbor_pairs[, 1L])

  # Read the second neighbor-pair column as target cell IDs.
  to <- as.integer(rook_neighbor_pairs[, 2L])

  # Keep only valid within-grid non-self neighbor pairs.
  valid <- !is.na(from) & !is.na(to) &
    from >= 1L & from <= n_cells &
    to >= 1L & to <= n_cells &
    from != to

  # Filter source cell IDs to valid pairs.
  from <- from[valid]

  # Filter target cell IDs to valid pairs.
  to <- to[valid]

  # If requested, make the neighbor relation symmetric.
  if (isTRUE(force_symmetric)) {
    # Add both directions for every pair.
    pairs <- data.table::data.table(
      from = c(from, to),
      to = c(to, from)
    )
  } else {
    # Keep the provided directionality unchanged.
    pairs <- data.table::data.table(
      from = from,
      to = to
    )
  }

  # Remove duplicate directed neighbor pairs.
  pairs <- unique(pairs, by = c("from", "to"))

  # Sort by source cell and then target cell.
  data.table::setorder(pairs, from, to)

  # Count how many neighbors each source cell has.
  neighbor_count <- tabulate(pairs$from, nbins = n_cells)

  # Compute the ending row position for each source cell in the compact edge list.
  end_index <- cumsum(neighbor_count)

  # Compute the starting row position for each source cell in the compact edge list.
  start_index <- end_index - neighbor_count + 1L

  # Mark cells with zero neighbors as having no start position.
  start_index[neighbor_count == 0L] <- NA_integer_

  # Mark cells with zero neighbors as having no end position.
  end_index[neighbor_count == 0L] <- NA_integer_

  # Return the compact neighbor index.
  list(
    from = as.integer(pairs$from),
    to = as.integer(pairs$to),
    start = as.integer(start_index),
    end = as.integer(end_index),
    n_cells = n_cells
  )
}



# =====================================================================
# HELPER: initialize incremental frontier state
# =====================================================================

initialize_frontier_state <- function(alive_species_count_by_cell, rook_neighbor_index) {
  # Mark cells as alive when at least one species remains in the cell.
  alive_by_cell <- alive_species_count_by_cell > 0L

  # Mark directed neighbor pairs where both cells are currently alive.
  alive_pair <- alive_by_cell[rook_neighbor_index$from] &
    alive_by_cell[rook_neighbor_index$to]

  # Count alive rook neighbors for each source cell.
  alive_neighbor_count_by_cell <- tabulate(
    rook_neighbor_index$from[alive_pair],
    nbins = rook_neighbor_index$n_cells
  )

  # A frontier cell is alive and has fewer than four alive rook neighbors.
  frontier_flag_by_cell <- alive_by_cell & alive_neighbor_count_by_cell < 4L

  # Return all vectors needed for incremental frontier updates.
  list(
    alive_by_cell = alive_by_cell,
    alive_neighbor_count_by_cell = as.integer(alive_neighbor_count_by_cell),
    frontier_flag_by_cell = frontier_flag_by_cell
  )
}



# =====================================================================
# HELPER: read current frontier cells
# =====================================================================

get_frontier_cells <- function(frontier_state) {
  # Return global raster-cell IDs currently marked as frontier cells.
  which(frontier_state$frontier_flag_by_cell)
}



# =====================================================================
# HELPER: update frontier state after cells become empty
# =====================================================================

update_frontier_state_after_empty_cells <- function(
  frontier_state,
  newly_empty_cells,
  rook_neighbor_index
) {
  # Convert newly empty cell IDs to unique integers.
  newly_empty_cells <- unique(as.integer(newly_empty_cells))

  # Keep only valid cell IDs inside the raster domain.
  newly_empty_cells <- newly_empty_cells[
    !is.na(newly_empty_cells) &
      newly_empty_cells >= 1L &
      newly_empty_cells <= rook_neighbor_index$n_cells
  ]

  # Keep only cells that the frontier state still considers alive.
  newly_empty_cells <- newly_empty_cells[
    frontier_state$alive_by_cell[newly_empty_cells]
  ]

  # Return unchanged state if there are no new empty cells to process.
  if (!length(newly_empty_cells)) {
    return(frontier_state)
  }

  # Mark newly empty cells as no longer alive.
  frontier_state$alive_by_cell[newly_empty_cells] <- FALSE

  # Mark newly empty cells as no longer frontier cells.
  frontier_state$frontier_flag_by_cell[newly_empty_cells] <- FALSE

  # Preallocate a small vector for neighbors affected by newly empty cells.
  affected_neighbors <- integer(length(newly_empty_cells) * 4L)

  # Initialize the write position in affected_neighbors.
  write_pos <- 0L

  # Loop over newly empty cells.
  for (cell in newly_empty_cells) {
    # Read this cell's first neighbor-edge position.
    start_pos <- rook_neighbor_index$start[cell]

    # Skip cells with no indexed rook neighbors.
    if (is.na(start_pos)) {
      next
    }

    # Read this cell's last neighbor-edge position.
    end_pos <- rook_neighbor_index$end[cell]

    # Extract this cell's rook-neighbor cell IDs.
    neighbors <- rook_neighbor_index$to[start_pos:end_pos]

    # Keep only neighbors that are still alive after the update.
    live_neighbors <- neighbors[frontier_state$alive_by_cell[neighbors]]

    # Skip if no live neighbors remain.
    if (!length(live_neighbors)) {
      next
    }

    # Decrement the alive-neighbor count of each live neighbor by one.
    frontier_state$alive_neighbor_count_by_cell[live_neighbors] <-
      frontier_state$alive_neighbor_count_by_cell[live_neighbors] - 1L

    # Count live neighbors affected by this newly empty cell.
    n_live <- length(live_neighbors)

    # Store these affected live-neighbor cell IDs.
    affected_neighbors[(write_pos + 1L):(write_pos + n_live)] <- live_neighbors

    # Advance the write position.
    write_pos <- write_pos + n_live
  }

  # If any live neighbors were affected, refresh only their frontier flags.
  if (write_pos > 0L) {
    # Deduplicate affected live-neighbor cell IDs.
    affected_neighbors <- unique(affected_neighbors[seq_len(write_pos)])

    # Recompute frontier status for affected live-neighbor cells only.
    frontier_state$frontier_flag_by_cell[affected_neighbors] <-
      frontier_state$alive_by_cell[affected_neighbors] &
      frontier_state$alive_neighbor_count_by_cell[affected_neighbors] < 4L
  }

  # Return the updated frontier state.
  frontier_state
}



# =====================================================================
# HELPER: compute current PU-level log scores
