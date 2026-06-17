# =====================================================================
# HELPER: read the species-specific dispersal threshold
# =====================================================================

get_species_dispersal_threshold_km <- function(
  species_name,        # character scalar naming the focal species
  species_params       # data.table containing species-level parameters
) {
  # Read this species' dispersal threshold in kilometers.
  dispersal_threshold_km <- species_params[
    species == species_name
  ]$dispersal_distance_km[1L]

  # Stop if the dispersal threshold is missing, non-finite, or non-positive.
  if (!is.finite(dispersal_threshold_km) || dispersal_threshold_km <= 0) {
    stop("Invalid dispersal_distance_km for species: ", species_name)
  }

  # Return the valid dispersal threshold.
  dispersal_threshold_km
}



# =====================================================================
# HELPER: identify lazy added-node overlay PU graphs
# =====================================================================

is_lazy_added_node_overlay_graph <- function(
  pu_graph              # one PU graph object
) {
  # Return TRUE only for list objects explicitly marked as lazy overlay graphs.
  is.list(pu_graph) &&
    !is.null(pu_graph$graph_type) &&
    identical(pu_graph$graph_type, "lazy_added_node_overlay")
}



# =====================================================================
# HELPER: get effective patch IDs from standard or lazy graph
# =====================================================================

get_pu_graph_id2patch <- function(
  pu_graph              # standard CSR graph or lazy overlay graph
) {
  # Lazy overlay graphs already expose the full effective id2patch vector.
  if (is_lazy_added_node_overlay_graph(pu_graph)) {
    return(as.integer(pu_graph$id2patch))
  }

  # Standard CSR graphs also store id2patch directly.
  as.integer(pu_graph$id2patch)
}



# =====================================================================
# HELPER: get graph-neighbor rows from standard or lazy graph
# =====================================================================

get_pu_graph_neighbor_rows <- function(
  pu_graph,             # standard CSR graph or lazy overlay graph
  row_index             # graph row whose neighbors should be returned
) {
  # Read the effective patch vector for this graph.
  id2patch <- get_pu_graph_id2patch(pu_graph)

  # Count graph nodes.
  n_nodes <- length(id2patch)

  # Coerce requested graph row to integer.
  row_index <- as.integer(row_index)

  # Stop if the requested row is outside the graph.
  if (row_index < 1L || row_index > n_nodes) {
    stop("Requested graph row is outside the graph node range.")
  }

  # Standard CSR graph path.
  if (!is_lazy_added_node_overlay_graph(pu_graph)) {
    # Read CSR row pointer.
    row_ptr <- as.integer(pu_graph$row_ptr)

    # Read CSR column index.
    col_idx <- as.integer(pu_graph$col_idx)

    # Compute first edge position for this row.
    first_edge_index <- row_ptr[row_index] + 1L

    # Compute last edge position for this row.
    last_edge_index <- row_ptr[row_index + 1L]

    # Return an empty vector if this graph row has no neighbors.
    if (last_edge_index < first_edge_index) {
      return(integer(0L))
    }

    # Return this row's CSR neighbor rows.
    return(as.integer(col_idx[first_edge_index:last_edge_index]))
  }

  # Lazy overlay graph path.
  base_graph <- pu_graph$base_graph

  # Read the number of nodes in the base graph.
  base_node_count <- as.integer(pu_graph$base_node_count)

  # Initialize base-neighbor vector.
  base_neighbors <- integer(0L)

  # If this row belongs to the base graph, read base CSR neighbors.
  if (row_index <= base_node_count) {
    # Read base CSR row pointer.
    base_row_ptr <- as.integer(base_graph$row_ptr)

    # Read base CSR column index.
    base_col_idx <- as.integer(base_graph$col_idx)

    # Compute first base edge position.
    first_edge_index <- base_row_ptr[row_index] + 1L

    # Compute last base edge position.
    last_edge_index <- base_row_ptr[row_index + 1L]

    # Read base neighbors if this base row has edges.
    if (last_edge_index >= first_edge_index) {
      base_neighbors <- as.integer(base_col_idx[first_edge_index:last_edge_index])
    }
  }

  # Initialize overlay-neighbor vector.
  overlay_neighbors <- integer(0L)

  # Read overlay-edge table.
  overlay_edges <- pu_graph$overlay_edges

  # If overlay edges exist, read edges incident to this row.
  if (!is.null(overlay_edges) && nrow(overlay_edges)) {
    # Collect neighbors from rows where this row is u or v.
    overlay_neighbors <- as.integer(c(
      overlay_edges[u == row_index, v],
      overlay_edges[v == row_index, u]
    ))
  }

  # Return unique sorted neighbors from base and overlay sources.
  sort(unique(as.integer(c(base_neighbors, overlay_neighbors))))
}



# =====================================================================
# HELPER: materialize one species' current patch raster from the live env
# =====================================================================
#
# The live state stores one integer vector per species:
#
#   cell -> patch_id, or NA if absent
#
# The distance stage needs a raster object only for local polygonization.
# This helper creates a temporary species raster from the live vector and a
# template raster layer.
#

materialize_species_patch_raster <- function(
  species_name,                    # character scalar naming the species
  species_raster_template_stack,   # multilayer SpatRaster template stack
  patch_id_by_species_env          # environment storing live species vectors
) {
  # Read the template layer names.
  template_layer_names <- names(species_raster_template_stack)

  # Stop if this species has no matching template layer.
  if (!(species_name %in% template_layer_names)) {
    stop("No template raster layer found for species: ", species_name)
  }

  # Extract the first matching template layer for this species.
  template_raster <- species_raster_template_stack[[
    which(template_layer_names == species_name)[1L]
  ]]

  # Create a new single-layer raster with the same geometry.
  species_raster <- terra::rast(template_raster)

  # Read this species' live cell -> patch_id vector.
  species_patch_ids <- get(
    species_name,
    envir = patch_id_by_species_env,
    inherits = FALSE
  )

  # Write the live patch IDs into the raster values.
  terra::values(species_raster) <- species_patch_ids

  # Return the temporary species patch raster.
  species_raster
}



# =====================================================================
# HELPER: extract candidate patch IDs for local distance rechecking
# =====================================================================
#
# The distance stage needs geometry for:
#
#   1. patches flagged by fragmentation, and
#   2. graph neighbors of those patches inside affected PUs.
#
# This avoids polygonizing every patch for a species.
#

extract_candidate_patch_ids_for_distance <- function(
  species_name,                         # focal species
  recheck_patch_ids_for_species,        # patches flagged by fragmentation
  affected_pu_ids,                      # PUs containing recheck patches
  pu_graphs_by_key                      # graph list keyed by "species|pu_id"
) {
  # Initialize candidate patches with the explicit recheck patches.
  candidate_patch_ids <- as.integer(recheck_patch_ids_for_species)

  # Loop over affected PUs.
  for (pu_id_value in affected_pu_ids) {
    # Build the graph key for this species-PU pair.
    graph_key <- paste0(species_name, "|", pu_id_value)

    # Read the PU graph.
    pu_graph <- pu_graphs_by_key[[graph_key]]

    # Skip missing graph entries.
    if (is.null(pu_graph)) {
      next
    }

    # Read the graph's effective node -> patch mapping.
    id2patch <- get_pu_graph_id2patch(pu_graph)

    # Identify graph rows corresponding to recheck patch IDs.
    recheck_rows <- which(id2patch %in% recheck_patch_ids_for_species)

    # Skip this PU if none of its rows are recheck patches.
    if (!length(recheck_rows)) {
      next
    }

    # Loop over recheck graph rows.
    for (row_index in recheck_rows) {
      # Read graph-neighbor rows for this recheck row.
      neighbor_rows <- get_pu_graph_neighbor_rows(
        pu_graph = pu_graph,
        row_index = row_index
      )

      # Add neighbor patch IDs to the candidate geometry set.
      if (length(neighbor_rows)) {
        candidate_patch_ids <- c(candidate_patch_ids, id2patch[neighbor_rows])
      }
    }
  }

  # Return sorted unique candidate patch IDs.
  sort(unique(as.integer(candidate_patch_ids)))
}



# =====================================================================
# HELPER: build sparse distance-predicate lookup for candidate patches
# =====================================================================
#
# This helper stores Boolean "within dispersal threshold" relationships.
# It does not store a dense numeric distance matrix.
#

build_distance_predicate_lookup <- function(
  candidate_patch_polygons,     # sf object with patch_id and geometry
  dispersal_threshold_km        # species-specific dispersal threshold in km
) {
  # Return an empty lookup if no polygons were supplied.
  if (nrow(candidate_patch_polygons) == 0L) {
    return(list(
      patch_row_by_id = integer(0L),
      within_patch_ids_by_patch_id = list(),
      directed_true_links = 0L
    ))
  }

  # Stop if the polygon object lacks patch IDs.
  if (!("patch_id" %in% names(candidate_patch_polygons))) {
    stop("candidate_patch_polygons must contain a patch_id column.")
  }

  # Read candidate patch IDs.
  patch_ids <- as.integer(candidate_patch_polygons$patch_id)

  # Stop if patch IDs are duplicated in the polygon table.
  if (anyDuplicated(patch_ids)) {
    duplicated_patch_ids <- unique(patch_ids[duplicated(patch_ids)])

    stop(
      "candidate_patch_polygons contains duplicated patch_id values: ",
      paste(utils::head(duplicated_patch_ids, 20L), collapse = ", "),
      if (length(duplicated_patch_ids) > 20L) " ..." else ""
    )
  }

  # Convert threshold from kilometers to meters for sf distance operations.
  threshold_m <- dispersal_threshold_km * 1000

  # Compute sparse within-distance row indices.
  within_rows <- sf::st_is_within_distance(
    x = candidate_patch_polygons,
    y = candidate_patch_polygons,
    dist = threshold_m,
    sparse = TRUE
  )

  # Convert row-index neighbors to patch-ID neighbors.
  within_patch_ids_by_patch_id <- stats::setNames(
    lapply(within_rows, function(row_indices) {
      as.integer(patch_ids[row_indices])
    }),
    as.character(patch_ids)
  )

  # Build patch_id -> polygon row lookup.
  patch_row_by_id <- stats::setNames(
    seq_along(patch_ids),
    as.character(patch_ids)
  )

  # Return the sparse predicate lookup.
  list(
    patch_row_by_id = patch_row_by_id,
    within_patch_ids_by_patch_id = within_patch_ids_by_patch_id,
    directed_true_links = sum(lengths(within_rows))
  )
}



# =====================================================================
# HELPER: build local patch polygons for one species
# =====================================================================
#
# This helper polygonizes only a local crop around candidate patches instead
# of polygonizing the full species raster.
#

build_local_patch_polygons <- function(
  species_patch_raster,       # temporary raster for one species
  candidate_patch_ids         # patch IDs needing local geometry
) {
  # Read current patch IDs from the raster.
  species_patch_ids <- terra::values(species_patch_raster, mat = FALSE)

  # Identify raster cells belonging to candidate patches.
  candidate_cells <- which(
    !is.na(species_patch_ids) &
      species_patch_ids %in% candidate_patch_ids
  )

  # Return an empty sf object if no candidate cells are present.
  if (!length(candidate_cells)) {
    return(sf::st_sf(
      patch_id = integer(),
      geometry = sf::st_sfc(
        crs = sf::st_crs(terra::crs(species_patch_raster, proj = TRUE))
      )
    ))
  }

  # Read x-y coordinates for candidate cells.
  candidate_xy <- terra::xyFromCell(species_patch_raster, candidate_cells)

  # Read raster resolution.
  raster_resolution <- terra::res(species_patch_raster)

  # Build a local crop extent with a one-cell margin.
  local_extent <- terra::ext(
    min(candidate_xy[, 1]) - raster_resolution[1],
    max(candidate_xy[, 1]) + raster_resolution[1],
    min(candidate_xy[, 2]) - raster_resolution[2],
    max(candidate_xy[, 2]) + raster_resolution[2]
  )

  # Crop the species raster to the local extent.
  local_raster <- terra::crop(
    species_patch_raster,
    local_extent,
    snap = "out"
  )

  # Read local raster patch IDs.
  local_patch_ids <- terra::values(local_raster, mat = FALSE)

  # Blank out all non-candidate patch IDs in the local raster.
  local_patch_ids[!(local_patch_ids %in% candidate_patch_ids)] <- NA_integer_

  # Write filtered local patch IDs back into the local raster.
  terra::values(local_raster) <- local_patch_ids

  # Return an empty sf object if the local raster now has no candidate cells.
  if (all(is.na(local_patch_ids))) {
    return(sf::st_sf(
      patch_id = integer(),
      geometry = sf::st_sfc(
        crs = sf::st_crs(terra::crs(species_patch_raster, proj = TRUE))
      )
    ))
  }

  # Polygonize the local candidate-patch raster.
  patch_polygons_vect <- terra::as.polygons(
    local_raster,
    dissolve = TRUE,
    na.rm = TRUE
  )

  # Convert terra vector polygons to sf.
  patch_polygons_sf <- sf::st_as_sf(patch_polygons_vect)

  # Read the active sf geometry column name.
  geometry_column_name <- attr(patch_polygons_sf, "sf_column")

  # Identify non-geometry attribute columns.
  non_geometry_columns <- setdiff(names(patch_polygons_sf), geometry_column_name)

  # Stop if no patch-ID attribute column was created.
  if (!length(non_geometry_columns)) {
    stop("Could not identify the patch-ID attribute column after polygonization.")
  }

  # Rename the first non-geometry column to patch_id.
  names(patch_polygons_sf)[
    names(patch_polygons_sf) == non_geometry_columns[1L]
  ] <- "patch_id"

  # Coerce polygon patch IDs to integer.
  patch_polygons_sf$patch_id <- as.integer(patch_polygons_sf$patch_id)

  # Repair invalid geometries before distance predicates are computed.
  patch_polygons_sf <- sf::st_make_valid(patch_polygons_sf)

  # Return the local candidate patch polygons.
  patch_polygons_sf
}



