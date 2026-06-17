# =====================================================================
# HELPER: rebuild the compact patch -> cells index for one species
# =====================================================================
#
# The live species state is a full cell -> patch_id vector. This helper
# rebuilds the compact sorted index used for fast patch-cell lookup.
#

rebuild_species_patch_index <- function(species_patch_ids) {
  # Find cells where this species is currently present.
  occupied_cells <- which(!is.na(species_patch_ids))

  # If the species has no remaining occupied cells, return no index.
  if (!length(occupied_cells)) {
    return(NULL)
  }

  # Read patch IDs for occupied cells.
  occupied_patch_ids <- species_patch_ids[occupied_cells]

  # Sort occupied cells by patch ID so each patch occupies a contiguous block.
  ordering <- order(occupied_patch_ids)

  # Return sorted patch IDs and their matching cell IDs.
  list(
    pid  = as.integer(occupied_patch_ids[ordering]),
    cell = as.integer(occupied_cells[ordering])
  )
}



# =====================================================================
# HELPER: summarize one species' current patch areas from its cell vector
# =====================================================================
#
# After fragmentation relabeling or patch dropping, this helper rebuilds
# patch areas directly from the live cell -> patch_id vector.
#

summarize_species_patch_areas <- function(
  species_name,
  species_patch_ids,
  cell_area_by_cell
) {
  # Find cells where this species is currently present.
  occupied_cells <- which(!is.na(species_patch_ids))

  # Return an empty patch-area table if the species has no occupied cells.
  if (!length(occupied_cells)) {
    return(data.table::data.table(
      species = character(),
      patch_id = integer(),
      patch_area_km2 = numeric()
    ))
  }

  # Read the patch ID assigned to each occupied cell.
  occupied_patch_ids <- species_patch_ids[occupied_cells]

  # Sum exact cell areas by patch ID.
  area_sum <- rowsum(
    cell_area_by_cell[occupied_cells],
    occupied_patch_ids,
    reorder = FALSE
  )

  # Return one row per current patch.
  data.table::data.table(
    species = species_name,
    patch_id = as.integer(rownames(area_sum)),
    patch_area_km2 = as.numeric(area_sum[, 1L])
  )
}



