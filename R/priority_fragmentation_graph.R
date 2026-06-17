# =====================================================================
# HELPER: build a provisional PU graph after fragmentation
# =====================================================================
#
# This is used only for the hard case where old patch nodes disappeared
# from a PU. In that case, the PU graph may have disconnected, so we need
# a provisional graph that can then be passed to rebuild_pu_after_patch_loss().
#
# Fragmentation-stage provisional rules:
#
#   - fragments from the same origin patch are provisionally linked
#   - fragments inherit the old origin patch's graph neighbors
#   - the later distance stage removes geometry-invalid links
#

build_provisional_pu_graph <- function(
  current_pu_patch_rows,
  old_pu_graph
) {
  # Count current patch rows in this PU.
  n_current_patches <- nrow(current_pu_patch_rows)

  # If no current patches remain, return an empty graph.
  if (n_current_patches == 0L) {
    return(list(
      pu_id = integer(0L),
      id2patch = integer(0L),
      row_ptr = 0L,
      col_idx = integer(0L)
    ))
  }

  # A provisional graph must be based on an old graph.
  if (is.null(old_pu_graph)) {
    stop("build_provisional_pu_graph() received NULL old_pu_graph.")
  }

  # Read current patch IDs.
  current_patch_ids <- as.integer(current_pu_patch_rows$patch_id)

  # Read the origin patch ID for each current patch.
  current_origin_patch_ids <- as.integer(current_pu_patch_rows$origin_patch_id)

  # Read the PU ID shared by this current PU row set.
  pu_id_value <- as.integer(current_pu_patch_rows$pu_id[1L])

  # Stop if patch IDs are duplicated inside this PU.
  if (anyDuplicated(current_patch_ids)) {
    duplicated_patch_ids <- unique(current_patch_ids[duplicated(current_patch_ids)])

    stop(
      "Duplicate patch_id values found inside one provisional PU graph: ",
      paste(utils::head(duplicated_patch_ids, 20L), collapse = ", "),
      if (length(duplicated_patch_ids) > 20L) " ..." else ""
    )
  }

  # Stop if required patch identity fields contain missing values.
  if (any(is.na(current_patch_ids)) || any(is.na(current_origin_patch_ids))) {
    stop("current_pu_patch_rows contains NA patch_id or origin_patch_id values.")
  }

  # Build a direct patch_id -> local graph row lookup.
  local_index_by_patch_id <- stats::setNames(
    seq_len(n_current_patches),
    as.character(current_patch_ids)
  )

  # Split current patch IDs by their origin patch.
  current_patches_by_origin <- split(
    current_patch_ids,
    current_origin_patch_ids
  )

  # Allocate a list of edge tables.
  edge_tables <- vector("list", 0L)

  # Read old graph patch IDs.
  old_patch_ids <- as.integer(old_pu_graph$id2patch)

  # Read old CSR row pointer.
  old_row_ptr <- as.integer(old_pu_graph$row_ptr)

  # Read old CSR column index.
  old_col_idx <- as.integer(old_pu_graph$col_idx)



  # -------------------------------------------------------------------
  # 1. Expand old graph edges through current origin groups
  # -------------------------------------------------------------------

  # Loop over old graph rows.
  for (old_u in seq_along(old_patch_ids)) {
    # Compute first edge position for old node old_u.
    first_edge_index <- old_row_ptr[old_u] + 1L

    # Compute last edge position for old node old_u.
    last_edge_index <- old_row_ptr[old_u + 1L]

    # Skip old nodes with no neighbors.
    if (last_edge_index < first_edge_index) {
      next
    }

    # Read old neighboring graph rows.
    old_neighbor_rows <- as.integer(old_col_idx[first_edge_index:last_edge_index])

    # Process each undirected old edge only once.
    old_neighbor_rows <- old_neighbor_rows[old_neighbor_rows > old_u]

    # Skip if no forward neighbor rows remain.
    if (!length(old_neighbor_rows)) {
      next
    }

    # Read the old origin patch ID for node old_u.
    origin_u <- old_patch_ids[old_u]

    # Read current patches descended from origin_u.
    current_u_patches <- current_patches_by_origin[[as.character(origin_u)]]

    # Skip if origin_u has no current descendants.
    if (is.null(current_u_patches) || !length(current_u_patches)) {
      next
    }

    # Loop over forward old-neighbor rows.
    for (old_v in old_neighbor_rows) {
      # Read the old origin patch ID for neighbor old_v.
      origin_v <- old_patch_ids[old_v]

      # Read current patches descended from origin_v.
      current_v_patches <- current_patches_by_origin[[as.character(origin_v)]]

      # Skip if origin_v has no current descendants.
      if (is.null(current_v_patches) || !length(current_v_patches)) {
        next
      }

      # Add all current descendant pairs as provisional edges.
      edge_tables[[length(edge_tables) + 1L]] <- data.table::CJ(
        patch_u = as.integer(current_u_patches),
        patch_v = as.integer(current_v_patches),
        sorted = FALSE
      )
    }
  }



  # -------------------------------------------------------------------
  # 2. Add sibling edges among fragments from the same origin patch
  # -------------------------------------------------------------------

  # Keep only origin groups with more than one current fragment.
  multi_fragment_origins <- current_patches_by_origin[
    lengths(current_patches_by_origin) > 1L
  ]

  # Add sibling-fragment edges.
  if (length(multi_fragment_origins)) {
    # Loop over origin groups that produced multiple fragments.
    for (fragment_patch_ids in multi_fragment_origins) {
      # Convert fragment patch IDs to integers.
      fragment_patch_ids <- as.integer(fragment_patch_ids)

      # Handle the common two-fragment case without combn().
      if (length(fragment_patch_ids) == 2L) {
        # Add the single sibling edge.
        edge_tables[[length(edge_tables) + 1L]] <- data.table::data.table(
          patch_u = fragment_patch_ids[1L],
          patch_v = fragment_patch_ids[2L]
        )
      } else {
        # Build all pairwise sibling edges for larger fragment groups.
        sibling_pairs <- utils::combn(fragment_patch_ids, 2L)

        # Add sibling edges to the edge table list.
        edge_tables[[length(edge_tables) + 1L]] <- data.table::data.table(
          patch_u = as.integer(sibling_pairs[1L, ]),
          patch_v = as.integer(sibling_pairs[2L, ])
        )
      }
    }
  }



  # -------------------------------------------------------------------
  # 3. Convert provisional patch-pair edges to CSR graph format
  # -------------------------------------------------------------------

  # If no edges were created, return an isolated-node graph.
  if (!length(edge_tables)) {
    return(list(
      pu_id = pu_id_value,
      id2patch = current_patch_ids,
      row_ptr = as.integer(rep(0L, n_current_patches + 1L)),
      col_idx = integer(0L)
    ))
  }

  # Combine all provisional edge tables.
  edge_table <- data.table::rbindlist(
    edge_tables,
    use.names = TRUE,
    fill = TRUE
  )

  # Drop self-edges.
  edge_table <- edge_table[patch_u != patch_v]

  # If all provisional edges were self-edges, return an isolated-node graph.
  if (!nrow(edge_table)) {
    return(list(
      pu_id = pu_id_value,
      id2patch = current_patch_ids,
      row_ptr = as.integer(rep(0L, n_current_patches + 1L)),
      col_idx = integer(0L)
    ))
  }

  # Canonicalize undirected patch edges.
  edge_table[, `:=`(
    patch_low = pmin(patch_u, patch_v),
    patch_high = pmax(patch_u, patch_v)
  )]

  # Deduplicate undirected patch edges.
  edge_table <- unique(
    edge_table[, .(patch_low, patch_high)],
    by = c("patch_low", "patch_high")
  )

  # Sort edges deterministically.
  data.table::setorder(edge_table, patch_low, patch_high)

  # Convert patch IDs to local graph row IDs.
  edge_table[, `:=`(
    u = as.integer(local_index_by_patch_id[as.character(patch_low)]),
    v = as.integer(local_index_by_patch_id[as.character(patch_high)])
  )]

  # Stop if any edge failed to map to a graph row.
  if (any(is.na(edge_table$u)) || any(is.na(edge_table$v))) {
    stop("Failed to map one or more provisional graph edges to local patch indices.")
  }

  # Create directed edges from undirected edges.
  directed_edges <- data.table::rbindlist(
    list(
      edge_table[, .(from = u, to = v)],
      edge_table[, .(from = v, to = u)]
    ),
    use.names = TRUE
  )

  # Sort directed edges by source row and then target row.
  data.table::setorder(directed_edges, from, to)

  # Count outgoing edges by source row.
  row_counts <- directed_edges[, .N, by = from]

  # Allocate CSR row pointer.
  row_ptr <- integer(n_current_patches + 1L)

  # Store row edge counts one position ahead for cumulative summing.
  row_ptr[row_counts$from + 1L] <- as.integer(row_counts$N)

  # Convert counts to CSR cumulative row pointer.
  row_ptr <- as.integer(cumsum(row_ptr))

  # Store CSR column index.
  col_idx <- as.integer(directed_edges$to)

  # Return the provisional CSR graph.
  list(
    pu_id = pu_id_value,
    id2patch = current_patch_ids,
    row_ptr = row_ptr,
    col_idx = col_idx
  )
}



# =====================================================================
# HELPER: make a lazy overlay graph for additions-only fragmentation
# =====================================================================
#
# This is the fast path for cases where fragmentation added new patch nodes
# but did not remove any old patch nodes.
#
# Adding nodes cannot disconnect an already connected graph, so this helper
# preserves the old CSR graph and stores only the provisional edges involving
# newly added fragment nodes. The distance stage can then consume this lazy
# graph and materialize or update it as needed.
#

make_lazy_added_node_overlay_pu_graph <- function(
  current_pu_patch_rows,
  old_pu_graph,
  new_patch_nodes
) {
  # A lazy overlay graph must wrap an existing old graph.
  if (is.null(old_pu_graph)) {
    stop("make_lazy_added_node_overlay_pu_graph() received NULL old_pu_graph.")
  }

  # Read current patch IDs.
  current_patch_ids <- as.integer(current_pu_patch_rows$patch_id)

  # Read current origin patch IDs.
  current_origin_patch_ids <- as.integer(current_pu_patch_rows$origin_patch_id)

  # Read old graph patch IDs.
  old_patch_ids <- as.integer(old_pu_graph$id2patch)

  # Read old CSR row pointer.
  old_row_ptr <- as.integer(old_pu_graph$row_ptr)

  # Read old CSR column index.
  old_col_idx <- as.integer(old_pu_graph$col_idx)

  # Stop if current patch IDs are duplicated.
  if (anyDuplicated(current_patch_ids)) {
    duplicated_patch_ids <- unique(current_patch_ids[duplicated(current_patch_ids)])

    stop(
      "Duplicate patch_id values found inside current_pu_patch_rows: ",
      paste(utils::head(duplicated_patch_ids, 20L), collapse = ", "),
      if (length(duplicated_patch_ids) > 20L) " ..." else ""
    )
  }

  # Stop if required current patch identity fields contain missing values.
  if (any(is.na(current_patch_ids)) || any(is.na(current_origin_patch_ids))) {
    stop("current_pu_patch_rows contains NA patch_id or origin_patch_id values.")
  }

  # Confirm that no old patch nodes are missing.
  missing_old_patch_ids <- setdiff(old_patch_ids, current_patch_ids)

  # Stop if this additions-only helper was called after old nodes disappeared.
  if (length(missing_old_patch_ids)) {
    stop(
      "Lazy overlay graph was requested even though old patch nodes are missing: ",
      paste(utils::head(missing_old_patch_ids, 20L), collapse = ", "),
      if (length(missing_old_patch_ids) > 20L) " ..." else ""
    )
  }

  # Preserve current patch-row order for added patch nodes.
  added_patch_ids <- unique(as.integer(
    current_patch_ids[current_patch_ids %in% as.integer(new_patch_nodes)]
  ))

  # Stop if no added nodes were supplied.
  if (!length(added_patch_ids)) {
    stop("make_lazy_added_node_overlay_pu_graph() was called with no added patch nodes.")
  }

  # Count original graph nodes.
  base_node_count <- length(old_patch_ids)

  # Build old patch_id -> old graph row lookup.
  old_row_by_patch_id <- stats::setNames(
    seq_along(old_patch_ids),
    as.character(old_patch_ids)
  )

  # Build added patch_id -> new graph row lookup.
  added_row_by_patch_id <- stats::setNames(
    base_node_count + seq_along(added_patch_ids),
    as.character(added_patch_ids)
  )

  # Combine old and added patch row lookups.
  full_row_by_patch_id <- c(
    old_row_by_patch_id,
    added_row_by_patch_id
  )

  # Build patch_id -> origin_patch_id lookup for current rows.
  current_origin_by_patch_id <- stats::setNames(
    current_origin_patch_ids,
    as.character(current_patch_ids)
  )

  # Split current patch IDs by origin patch.
  current_patches_by_origin <- split(
    current_patch_ids,
    current_origin_patch_ids
  )

  # Allocate overlay-edge table list.
  overlay_edge_tables <- vector("list", length(added_patch_ids))

  # Initialize overlay-edge table write index.
  overlay_edge_table_index <- 0L

  # Loop over newly added patch nodes.
  for (added_patch_id in added_patch_ids) {
    # Convert added patch ID to character key.
    added_patch_key <- as.character(added_patch_id)

    # Read graph row for this added patch.
    added_row <- as.integer(added_row_by_patch_id[added_patch_key])

    # Read the origin patch for this added patch.
    origin_patch_id <- as.integer(current_origin_by_patch_id[added_patch_key])

    # Stop if origin lookup failed.
    if (is.na(origin_patch_id)) {
      stop("Could not find origin_patch_id for added patch ", added_patch_id, ".")
    }

    # Read the old graph row corresponding to the origin patch.
    origin_row <- as.integer(old_row_by_patch_id[as.character(origin_patch_id)])

    # Stop if the origin patch is not present in the old graph.
    if (is.na(origin_row)) {
      stop(
        "Added patch ",
        added_patch_id,
        " has origin_patch_id ",
        origin_patch_id,
        ", but that origin is not present in the old PU graph."
      )
    }

    # Compute first old CSR edge position for the origin patch.
    first_edge_index <- old_row_ptr[origin_row] + 1L

    # Compute last old CSR edge position for the origin patch.
    last_edge_index <- old_row_ptr[origin_row + 1L]

    # Initialize neighbor origin patch IDs.
    neighbor_origin_patch_ids <- integer(0L)

    # If the origin node has old graph neighbors, read their origin patch IDs.
    if (last_edge_index >= first_edge_index) {
      # Read old neighboring graph rows.
      neighbor_rows <- as.integer(old_col_idx[first_edge_index:last_edge_index])

      # Convert old neighboring graph rows to old patch IDs.
      neighbor_origin_patch_ids <- as.integer(old_patch_ids[neighbor_rows])
    }

    # Added nodes provisionally connect to same-origin fragments and old-neighbor descendants.
    target_origin_ids <- unique(as.integer(c(
      origin_patch_id,
      neighbor_origin_patch_ids
    )))

    # Keep target origins that currently have descendants.
    target_origin_names <- intersect(
      as.character(target_origin_ids),
      names(current_patches_by_origin)
    )

    # Read all current patch IDs descended from those target origins.
    target_patch_ids <- unique(as.integer(unlist(
      current_patches_by_origin[target_origin_names],
      use.names = FALSE
    )))

    # Do not connect the added patch to itself.
    target_patch_ids <- target_patch_ids[target_patch_ids != added_patch_id]

    # Skip this added patch if it has no provisional targets.
    if (!length(target_patch_ids)) {
      next
    }

    # Map target patch IDs to graph rows.
    target_rows <- as.integer(
      full_row_by_patch_id[as.character(target_patch_ids)]
    )

    # Stop if any target row failed to map.
    if (any(is.na(target_rows))) {
      stop(
        "Failed to map one or more target patches to graph rows while ",
        "building lazy overlay edges for added patch ",
        added_patch_id,
        "."
      )
    }

    # Canonicalize overlay edges by graph-row ID.
    edge_low <- pmin(added_row, target_rows)

    # Canonicalize overlay edges by graph-row ID.
    edge_high <- pmax(added_row, target_rows)

    # Advance overlay-edge table write index.
    overlay_edge_table_index <- overlay_edge_table_index + 1L

    # Store overlay edges for this added patch.
    overlay_edge_tables[[overlay_edge_table_index]] <- data.table::data.table(
      u = as.integer(edge_low),
      v = as.integer(edge_high)
    )
  }

  # Combine overlay edge tables if any were created.
  overlay_edges <- if (overlay_edge_table_index > 0L) {
    # Bind only initialized edge-table slots.
    overlay_edges_raw <- data.table::rbindlist(
      overlay_edge_tables[seq_len(overlay_edge_table_index)],
      use.names = TRUE
    )

    # Drop self-edges.
    overlay_edges_raw <- overlay_edges_raw[u != v]

    # Deduplicate overlay edges.
    overlay_edges_raw <- unique(overlay_edges_raw, by = c("u", "v"))

    # Sort overlay edges deterministically.
    data.table::setorder(overlay_edges_raw, u, v)

    # Return the overlay edge table.
    overlay_edges_raw
  } else {
    # Return an empty overlay edge table.
    data.table::data.table(
      u = integer(),
      v = integer()
    )
  }

  # Return the lazy overlay graph object.
  list(
    graph_type = "lazy_added_node_overlay",
    species = old_pu_graph$species,
    pu_id = as.integer(old_pu_graph$pu_id),
    id2patch = as.integer(c(old_patch_ids, added_patch_ids)),
    base_graph = old_pu_graph,
    base_node_count = as.integer(base_node_count),
    added_patch_ids = as.integer(added_patch_ids),
    added_origin_patch_ids = as.integer(
      current_origin_by_patch_id[as.character(added_patch_ids)]
    ),
    overlay_edges = overlay_edges
  )
}
