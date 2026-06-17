source(file.path("R", "priority_csr_graph.R"))



# =====================================================================
# HELPER: rebuild one PU after distance-invalid edges are removed
# =====================================================================

rebuild_pu_after_edge_filter <- function(
  pu_graph,              # filtered CSR graph for one PU
  patch_area_by_patch,   # named numeric vector: patch_id -> area
  pu_area_threshold,     # species PU-area threshold
  next_available_pu_id   # next available PU ID if splitting occurs
) {
  pu_area_threshold <- validate_area_threshold(
    pu_area_threshold,
    "pu_area_threshold"
  )

  # Read graph node -> patch IDs.
  id2patch <- as.integer(pu_graph$id2patch)

  # Read CSR row pointer.
  row_ptr <- as.integer(pu_graph$row_ptr)

  # Read CSR column index.
  col_idx <- as.integer(pu_graph$col_idx)

  # Count graph nodes.
  n_nodes <- length(id2patch)

  # Return no surviving graphs if the PU has no nodes.
  if (n_nodes == 0L) {
    return(list(
      surviving_pu_graphs = list(),
      dropped_patch_ids = integer(0L),
      next_available_pu_id = next_available_pu_id
    ))
  }

  # Handle one-node graph without component search.
  if (n_nodes == 1L) {
    # Read the single patch area.
    patch_area <- as.numeric(patch_area_by_patch[as.character(id2patch[1L])])

    # Drop the single patch if its area is missing or below threshold.
    if (!area_exceeds_threshold(patch_area, pu_area_threshold)) {
      return(list(
        surviving_pu_graphs = list(),
        dropped_patch_ids = id2patch,
        next_available_pu_id = next_available_pu_id
      ))
    }

    # Otherwise return a one-node surviving PU graph.
    return(list(
      surviving_pu_graphs = list(list(
        pu_id = pu_graph$pu_id,
        id2patch = id2patch,
        row_ptr = c(0L, 0L),
        col_idx = integer(0L)
      )),
      dropped_patch_ids = integer(0L),
      next_available_pu_id = next_available_pu_id
    ))
  }

  # Label connected components in the filtered graph.
  component_id_by_node <- find_csr_components(row_ptr, col_idx)

  # Read patch area for each graph node.
  patch_area_by_node <- as.numeric(patch_area_by_patch[as.character(id2patch)])

  # Sum patch area by connected component.
  component_area <- as.numeric(
    rowsum(
      patch_area_by_node,
      component_id_by_node,
      reorder = FALSE
    )
  )

  # Components strictly above threshold survive.
  surviving_component_labels <- which(
    area_exceeds_threshold(component_area, pu_area_threshold)
  )

  # Components at or below threshold are dropped.
  dropped_component_labels <- which(
    !area_exceeds_threshold(component_area, pu_area_threshold)
  )

  # Allocate surviving graph list.
  surviving_pu_graphs <- vector("list", length(surviving_component_labels))

  # Assign PU IDs to surviving components.
  if (length(surviving_component_labels) > 0L) {
    # If only one component survives, keep the original PU ID.
    if (length(surviving_component_labels) == 1L) {
      surviving_pu_ids <- pu_graph$pu_id
    } else {
      # First surviving component keeps the old ID; others receive new IDs.
      surviving_pu_ids <- c(
        pu_graph$pu_id,
        next_available_pu_id + seq_len(length(surviving_component_labels) - 1L)
      )

      # Advance the next available PU ID.
      next_available_pu_id <- surviving_pu_ids[length(surviving_pu_ids)]
    }

    # Build one CSR subgraph per surviving component.
    for (j in seq_along(surviving_component_labels)) {
      # Read this surviving component label.
      component_label <- surviving_component_labels[j]

      # Identify original node indices in this component.
      component_node_indices <- which(component_id_by_node == component_label)

      # Build CSR subgraph for the component.
      component_subgraph <- build_csr_subgraph(
        row_ptr = row_ptr,
        col_idx = col_idx,
        kept_node_indices = component_node_indices
      )

      # Store this surviving PU graph.
      surviving_pu_graphs[[j]] <- list(
        pu_id = as.integer(surviving_pu_ids[j]),
        id2patch = as.integer(id2patch[component_node_indices]),
        row_ptr = as.integer(component_subgraph$row_ptr),
        col_idx = as.integer(component_subgraph$col_idx)
      )
    }
  }

  # Convert dropped component nodes to dropped patch IDs.
  dropped_patch_ids <- id2patch[
    component_id_by_node %in% dropped_component_labels
  ]

  # Return surviving graphs, dropped patches, and next available PU ID.
  list(
    surviving_pu_graphs = surviving_pu_graphs,
    dropped_patch_ids = as.integer(dropped_patch_ids),
    next_available_pu_id = as.integer(next_available_pu_id)
  )
}



# =====================================================================
# HELPER: materialize lazy overlay graph as a standard CSR graph
# =====================================================================
#
# This is used when a lazy graph reaches the distance stage but the current
# affected PU no longer contains any recheck patches after restriction.
# The graph is still converted to standard CSR so no lazy graphs remain after
# the distance stage.
#

materialize_lazy_added_node_overlay_graph <- function(
  pu_graph              # standard CSR graph or lazy overlay graph
) {
  # Return standard graphs unchanged.
  if (!is_lazy_added_node_overlay_graph(pu_graph)) {
    return(pu_graph)
  }

  # Read effective node -> patch IDs.
  id2patch <- as.integer(pu_graph$id2patch)

  # Count effective graph nodes.
  n_nodes <- length(id2patch)

  # Read wrapped base graph.
  base_graph <- pu_graph$base_graph

  # Read base graph node count.
  base_node_count <- as.integer(pu_graph$base_node_count)

  # Allocate adjacency list.
  adjacency_by_row <- vector("list", n_nodes)

  # Initialize each adjacency row as empty.
  for (row_index in seq_len(n_nodes)) {
    adjacency_by_row[[row_index]] <- integer(0L)
  }

  # Read base CSR row pointer.
  base_row_ptr <- as.integer(base_graph$row_ptr)

  # Read base CSR column index.
  base_col_idx <- as.integer(base_graph$col_idx)

  # Copy base graph edges.
  for (row_index in seq_len(base_node_count)) {
    # Compute first base edge position.
    first_edge_index <- base_row_ptr[row_index] + 1L

    # Compute last base edge position.
    last_edge_index <- base_row_ptr[row_index + 1L]

    # Copy base neighbors if this row has edges.
    if (last_edge_index >= first_edge_index) {
      adjacency_by_row[[row_index]] <- as.integer(
        base_col_idx[first_edge_index:last_edge_index]
      )
    }
  }

  # Read overlay edges.
  overlay_edges <- pu_graph$overlay_edges

  # Add overlay edges if any exist.
  if (!is.null(overlay_edges) && nrow(overlay_edges)) {
    # Loop over overlay-edge rows.
    for (edge_row in seq_len(nrow(overlay_edges))) {
      # Read first endpoint.
      u <- as.integer(overlay_edges$u[edge_row])

      # Read second endpoint.
      v <- as.integer(overlay_edges$v[edge_row])

      # Add v to u's adjacency.
      adjacency_by_row[[u]] <- c(adjacency_by_row[[u]], v)

      # Add u to v's adjacency.
      adjacency_by_row[[v]] <- c(adjacency_by_row[[v]], u)
    }
  }

  # Deduplicate and sort each adjacency row.
  adjacency_by_row <- lapply(adjacency_by_row, function(neighbor_rows) {
    sort(unique(as.integer(neighbor_rows)))
  })

  # Count neighbors in each row.
  row_lengths <- vapply(adjacency_by_row, length, integer(1L))

  # Build CSR row pointer.
  row_ptr <- as.integer(c(0L, cumsum(row_lengths)))

  # Build CSR column index.
  col_idx <- as.integer(unlist(adjacency_by_row, use.names = FALSE))

  # Return standard CSR graph.
  list(
    species = pu_graph$species,
    pu_id = as.integer(pu_graph$pu_id),
    id2patch = id2patch,
    row_ptr = row_ptr,
    col_idx = col_idx
  )
}



# =====================================================================
# HELPER: filter one PU graph by distance predicate validity
# =====================================================================
#
# Only edges incident to recheck patches are re-evaluated.
# Existing edges not incident to recheck patches are kept as-is.
#

filter_distance_invalid_edges_for_pu <- function(
  pu_graph,                       # standard CSR or lazy overlay graph
  recheck_patch_ids_in_pu,        # patch IDs whose incident edges need recheck
  distance_predicate_lookup       # sparse within-threshold lookup
) {
  # Read effective node -> patch IDs.
  id2patch <- get_pu_graph_id2patch(pu_graph)

  # Count graph nodes.
  n_nodes <- length(id2patch)

  # Return an empty-edge graph for zero- or one-node PUs.
  if (n_nodes <= 1L) {
    return(list(
      pu_id = pu_graph$pu_id,
      id2patch = id2patch,
      row_ptr = if (n_nodes == 0L) 0L else c(0L, 0L),
      col_idx = integer(0L)
    ))
  }

  # Coerce recheck patch IDs to integer.
  recheck_patch_ids_in_pu <- as.integer(recheck_patch_ids_in_pu)

  # Build a named lookup for O(1)-style membership checks.
  recheck_lookup <- stats::setNames(
    rep(TRUE, length(recheck_patch_ids_in_pu)),
    as.character(recheck_patch_ids_in_pu)
  )

  # Allocate adjacency list for the filtered graph.
  adjacency_by_row <- vector("list", n_nodes)

  # Initialize each adjacency row as empty.
  for (row_index in seq_len(n_nodes)) {
    adjacency_by_row[[row_index]] <- integer(0L)
  }

  # Define local helper for one undirected graph edge.
  process_undirected_edge <- function(u, v) {
    # Coerce first graph row to integer.
    u <- as.integer(u)

    # Coerce second graph row to integer.
    v <- as.integer(v)

    # Ignore self-edges.
    if (u == v) {
      return(invisible(NULL))
    }

    # Convert first graph row to patch ID.
    patch_u <- as.integer(id2patch[u])

    # Convert second graph row to patch ID.
    patch_v <- as.integer(id2patch[v])

    # Check whether patch_u needs incident-edge rechecking.
    patch_u_needs_recheck <- !is.na(recheck_lookup[as.character(patch_u)])

    # Check whether patch_v needs incident-edge rechecking.
    patch_v_needs_recheck <- !is.na(recheck_lookup[as.character(patch_v)])

    # Keep non-rechecked edges unchanged.
    if (!(patch_u_needs_recheck || patch_v_needs_recheck)) {
      keep_edge <- TRUE
    } else {
      # Convert patch_u to a character lookup key.
      patch_u_key <- as.character(patch_u)

      # Convert patch_v to a character lookup key.
      patch_v_key <- as.character(patch_v)

      # Confirm patch_u has geometry in the distance predicate lookup.
      row_u <- unname(distance_predicate_lookup$patch_row_by_id[patch_u_key])

      # Confirm patch_v has geometry in the distance predicate lookup.
      row_v <- unname(distance_predicate_lookup$patch_row_by_id[patch_v_key])

      # Stop if either endpoint lacks geometry.
      if (is.na(row_u) || is.na(row_v)) {
        stop(
          "Distance predicate geometry missing for patch pair (",
          patch_u, ", ", patch_v, ")."
        )
      }

      # Read patches within threshold of patch_u.
      within_patch_ids_for_u <-
        distance_predicate_lookup$within_patch_ids_by_patch_id[[patch_u_key]]

      # Keep this edge only if patch_v is within patch_u's threshold neighborhood.
      keep_edge <- patch_v %in% within_patch_ids_for_u
    }

    # If valid, add both directed adjacency entries.
    if (keep_edge) {
      adjacency_by_row[[u]] <<- c(adjacency_by_row[[u]], v)
      adjacency_by_row[[v]] <<- c(adjacency_by_row[[v]], u)
    }

    # Return invisibly because this helper mutates adjacency_by_row.
    invisible(NULL)
  }



  # -------------------------------------------------------------------
  # Process edges from lazy overlay graphs
  # -------------------------------------------------------------------

  if (is_lazy_added_node_overlay_graph(pu_graph)) {
    # Read base graph.
    base_graph <- pu_graph$base_graph

    # Read base node count.
    base_node_count <- as.integer(pu_graph$base_node_count)

    # Read base CSR row pointer.
    base_row_ptr <- as.integer(base_graph$row_ptr)

    # Read base CSR column index.
    base_col_idx <- as.integer(base_graph$col_idx)

    # Process base graph edges once as undirected edges.
    for (u in seq_len(base_node_count)) {
      # Compute first base edge position.
      first_edge_index <- base_row_ptr[u] + 1L

      # Compute last base edge position.
      last_edge_index <- base_row_ptr[u + 1L]

      # Skip base rows with no edges.
      if (last_edge_index < first_edge_index) {
        next
      }

      # Read neighbor rows.
      neighbor_rows <- as.integer(base_col_idx[first_edge_index:last_edge_index])

      # Process each undirected edge once.
      for (v in neighbor_rows) {
        if (v <= u) {
          next
        }

        process_undirected_edge(u, v)
      }
    }

    # Read overlay edges.
    overlay_edges <- pu_graph$overlay_edges

    # Process overlay edges if present.
    if (!is.null(overlay_edges) && nrow(overlay_edges)) {
      # Loop over overlay-edge rows.
      for (edge_row in seq_len(nrow(overlay_edges))) {
        # Process this overlay edge.
        process_undirected_edge(
          u = overlay_edges$u[edge_row],
          v = overlay_edges$v[edge_row]
        )
      }
    }



  # -------------------------------------------------------------------
  # Process edges from standard CSR graphs
  # -------------------------------------------------------------------

  } else {
    # Read CSR row pointer.
    row_ptr <- as.integer(pu_graph$row_ptr)

    # Read CSR column index.
    col_idx <- as.integer(pu_graph$col_idx)

    # Loop over graph rows.
    for (u in seq_len(n_nodes)) {
      # Compute first edge position.
      first_edge_index <- row_ptr[u] + 1L

      # Compute last edge position.
      last_edge_index <- row_ptr[u + 1L]

      # Skip rows with no edges.
      if (last_edge_index < first_edge_index) {
        next
      }

      # Loop over row edges.
      for (edge_index in first_edge_index:last_edge_index) {
        # Read neighbor row.
        v <- as.integer(col_idx[edge_index])

        # Process each undirected edge only once.
        if (v <= u) {
          next
        }

        process_undirected_edge(u, v)
      }
    }
  }

  # Deduplicate and sort each filtered adjacency row.
  adjacency_by_row <- lapply(adjacency_by_row, function(neighbor_rows) {
    sort(unique(as.integer(neighbor_rows)))
  })

  # Count neighbors by row.
  row_lengths <- vapply(adjacency_by_row, length, integer(1L))

  # Build CSR row pointer.
  new_row_ptr <- as.integer(c(0L, cumsum(row_lengths)))

  # Build CSR column index.
  new_col_idx <- as.integer(unlist(adjacency_by_row, use.names = FALSE))

  # Return standard CSR graph after distance filtering.
  list(
    pu_id = pu_graph$pu_id,
    id2patch = id2patch,
    row_ptr = new_row_ptr,
    col_idx = new_col_idx
  )
}



# =====================================================================
