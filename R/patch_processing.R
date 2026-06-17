log_patch_step <- function(species, status) {
  message(sprintf("[%s] [%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), species, status))
}

build_pu_connectivity <- function(nb, patches_sf, patch_to_pu_kept, final_patch_ids, species_name) {
  edge_i <- integer()
  edge_j <- integer()

  for (i in seq_along(nb)) {
    nbrs <- as.integer(nb[[i]])
    nbrs <- nbrs[nbrs > i]
    if (length(nbrs)) {
      edge_i <- c(edge_i, rep.int(i, length(nbrs)))
      edge_j <- c(edge_j, nbrs)
    }
  }

  edges <- data.frame(
    patch_i = patches_sf$patch_id[edge_i],
    patch_j = patches_sf$patch_id[edge_j]
  )

  keep_prefilter_ids <- as.integer(names(final_patch_ids))
  edges <- edges[
    edges$patch_i %in% keep_prefilter_ids & edges$patch_j %in% keep_prefilter_ids,
    ,
    drop = FALSE
  ]

  if (nrow(edges)) {
    edges$patch_i <- as.integer(final_patch_ids[as.character(edges$patch_i)])
    edges$patch_j <- as.integer(final_patch_ids[as.character(edges$patch_j)])
    edges <- edges[stats::complete.cases(edges), , drop = FALSE]
  }

  out <- list()

  for (pu in sort(unique(patch_to_pu_kept$pu_id))) {
    patch_ids <- sort(unique(patch_to_pu_kept$patch_id[patch_to_pu_kept$pu_id == pu]))
    patch_ids <- patch_ids[!is.na(patch_ids)]
    if (!length(patch_ids)) next

    idx_map <- seq_along(patch_ids)
    names(idx_map) <- as.character(patch_ids)

    sub_edges <- edges[
      edges$patch_i %in% patch_ids & edges$patch_j %in% patch_ids,
      ,
      drop = FALSE
    ]

    if (!nrow(sub_edges)) {
      row_ptr <- as.integer(c(0L, rep.int(0L, length(patch_ids))))
      col_idx <- integer(0)
    } else {
      symmetric_edges <- rbind(
        data.frame(i = sub_edges$patch_i, j = sub_edges$patch_j),
        data.frame(i = sub_edges$patch_j, j = sub_edges$patch_i)
      )
      symmetric_edges <- unique(symmetric_edges[symmetric_edges$i != symmetric_edges$j, , drop = FALSE])
      symmetric_edges <- symmetric_edges[order(symmetric_edges$i, symmetric_edges$j), , drop = FALSE]

      row_local <- as.integer(idx_map[as.character(symmetric_edges$i)])
      col_local <- as.integer(idx_map[as.character(symmetric_edges$j)])
      ok <- !is.na(row_local) & !is.na(col_local)
      row_local <- row_local[ok]
      col_local <- col_local[ok]

      split_neighbors <- split(col_local - 1L, row_local)
      adj_list <- vector("list", length(patch_ids))
      for (i in seq_along(adj_list)) {
        adj_list[[i]] <- as.integer(sort(unique(split_neighbors[[as.character(i)]] %||% integer(0))))
      }

      lens <- vapply(adj_list, length, integer(1))
      row_ptr <- as.integer(c(0L, cumsum(lens)))
      col_idx <- as.integer(unlist(adj_list, use.names = FALSE))
    }

    key <- paste0(species_name, "|", pu)
    out[[key]] <- list(
      species   = species_name,
      pu_id     = as.integer(pu),
      patch_ids = as.integer(patch_ids),
      row_ptr   = as.integer(row_ptr),
      col_idx   = as.integer(col_idx)
    )
  }

  out
}

patch_components_for_species <- function(row, mapped_habitat, species_name = row$scientificName) {
  mapped_habitat_area_km2 <- raster_area_km2(mapped_habitat)

  mapped_habitat_trim <- terra::trim(mapped_habitat)
  patch_raw <- clump_binary_raster(mapped_habitat_trim, diagonal = FALSE)
  names(patch_raw) <- "patch_id_raw"

  cell_area <- terra::cellSize(mapped_habitat_trim, unit = "km")
  patch_area <- as.data.frame(terra::zonal(cell_area, patch_raw, "sum", na.rm = TRUE))
  if (!nrow(patch_area)) {
    return(list(status = "skipped_no_patch_ids", species = species_name))
  }

  names(patch_area) <- c("patch_id_raw", "patch_area_km2")
  patch_area <- patch_area[!is.na(patch_area$patch_id_raw), , drop = FALSE]
  if (!nrow(patch_area)) {
    return(list(status = "skipped_no_valid_patch_ids", species = species_name))
  }

  keep_raw_patches <- patch_area$patch_id_raw[
    area_exceeds_threshold(patch_area$patch_area_km2, row$min_patch_km2)
  ]
  if (!length(keep_raw_patches)) {
    return(list(
      status = "skipped_no_patches_above_min_patch",
      species = species_name,
      mapped_habitat_area_km2 = mapped_habitat_area_km2,
      n_raw_patches = nrow(patch_area)
    ))
  }

  patch_prefilter_trim <- terra::classify(
    patch_raw,
    rcl    = cbind(as.integer(keep_raw_patches), seq_along(keep_raw_patches)),
    others = NA_integer_
  )
  names(patch_prefilter_trim) <- "patch_id"

  patch_prefilter <- terra::extend(patch_prefilter_trim, mapped_habitat)
  names(patch_prefilter) <- "patch_id"

  patch_table_prefilter <- data.frame(
    patch_id_prefilter = seq_along(keep_raw_patches),
    patch_area_km2 = patch_area$patch_area_km2[match(keep_raw_patches, patch_area$patch_id_raw)]
  )

  patch_polygons <- tryCatch(
    sf::st_as_sf(terra::as.polygons(patch_prefilter, values = TRUE, dissolve = TRUE, na.rm = TRUE)),
    error = function(e) NULL
  )
  if (is.null(patch_polygons) || !nrow(patch_polygons)) {
    return(list(status = "skipped_polygonization_failed", species = species_name))
  }

  value_col <- setdiff(names(patch_polygons), attr(patch_polygons, "sf_column"))[1]
  patch_polygons$patch_id <- as.integer(patch_polygons[[value_col]])
  patch_polygons <- patch_polygons[!is.na(patch_polygons$patch_id), ]
  if (!nrow(patch_polygons)) {
    return(list(status = "skipped_no_patch_polygons", species = species_name))
  }

  neighbors <- sf::st_is_within_distance(patch_polygons, patch_polygons, dist = row$disp_km * 1000)
  patch_graph <- igraph::graph_from_adj_list(neighbors, mode = "out") |>
    igraph::as_undirected("collapse") |>
    igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE)

  component_id <- igraph::components(patch_graph)$membership
  patch_polygons$pu_id_prefilter <- as.integer(match(component_id, sort(unique(component_id))))

  patch_to_pu_prefilter <- data.frame(
    patch_id_prefilter = patch_polygons$patch_id,
    pu_id_prefilter = patch_polygons$pu_id_prefilter
  )

  patch_area_with_pu <- merge(patch_table_prefilter, patch_to_pu_prefilter, by = "patch_id_prefilter")
  pu_area <- aggregate(patch_area_km2 ~ pu_id_prefilter, patch_area_with_pu, sum, na.rm = TRUE)
  names(pu_area) <- c("pu_id_prefilter", "pu_area_km2")

  keep_prefilter_pu <- sort(pu_area$pu_id_prefilter[
    area_exceeds_threshold(pu_area$pu_area_km2, row$min_pop_km2)
  ])
  if (!length(keep_prefilter_pu)) {
    return(list(
      status = "skipped_no_pus_above_min_pop",
      species = species_name,
      mapped_habitat_area_km2 = mapped_habitat_area_km2,
      n_raw_patches = nrow(patch_area),
      n_patches_after_patch_filter = length(keep_raw_patches),
      n_candidate_pus = nrow(pu_area)
    ))
  }

  patch_to_pu_prefilter$pu_id <- match(patch_to_pu_prefilter$pu_id_prefilter, keep_prefilter_pu)
  patch_to_pu_kept <- patch_to_pu_prefilter[!is.na(patch_to_pu_prefilter$pu_id), , drop = FALSE]

  keep_prefilter_patch_ids <- sort(unique(patch_to_pu_kept$patch_id_prefilter))
  final_patch_ids <- seq_along(keep_prefilter_patch_ids)
  names(final_patch_ids) <- as.character(keep_prefilter_patch_ids)

  patch_final <- terra::classify(
    patch_prefilter,
    rcl    = cbind(as.integer(keep_prefilter_patch_ids), final_patch_ids),
    others = NA_integer_
  )
  names(patch_final) <- "patch_id"

  patch_to_pu_kept$patch_id <- as.integer(final_patch_ids[as.character(patch_to_pu_kept$patch_id_prefilter)])

  patch_lookup <- data.frame(
    scientificName = species_name,
    patch_id = patch_to_pu_kept$patch_id,
    pu_id = patch_to_pu_kept$pu_id,
    patch_area_km2 = patch_table_prefilter$patch_area_km2[
      match(patch_to_pu_kept$patch_id_prefilter, patch_table_prefilter$patch_id_prefilter)
    ],
    stringsAsFactors = FALSE
  )

  connectivity <- build_pu_connectivity(
    nb = neighbors,
    patches_sf = patch_polygons,
    patch_to_pu_kept = patch_to_pu_kept[, c("patch_id", "pu_id")],
    final_patch_ids = final_patch_ids,
    species_name = species_name
  )

  list(
    status = "retained",
    species = species_name,
    mapped_habitat_area_km2 = mapped_habitat_area_km2,
    n_raw_patches = nrow(patch_area),
    n_patches_after_patch_filter = length(keep_raw_patches),
    n_candidate_pus = nrow(pu_area),
    n_final_pus = length(keep_prefilter_pu),
    n_final_patches = length(unique(patch_lookup$patch_id)),
    patch_lookup = patch_lookup,
    connectivity = connectivity,
    patch_final = patch_final
  )
}

process_patch_species <- function(row, paths, run_options, habitat_masks, template, label_map = habitat_label_map()) {
  species_name <- row$scientificName
  log_patch_step(species_name, "START")

  log_patch_step(species_name, "MAPPED_HABITAT")
  mapped <- mapped_habitat_for_species(row, habitat_masks, template, label_map)
  if (!identical(mapped$status, "ok")) {
    return(list(status = mapped$status, species = species_name))
  }
  if (!has_non_na(mapped$mapped_habitat)) {
    return(list(status = "skipped_empty_mapped_habitat", species = species_name))
  }

  log_patch_step(species_name, "PATCHES_AND_POPULATION_UNITS")
  components <- patch_components_for_species(row, mapped$mapped_habitat, species_name)
  if (!identical(components$status, "retained")) {
    return(components)
  }

  out_file <- file.path(paths$patch_dir, patch_filename_from_scientific(species_name))
  if (file.exists(out_file) && !run_options$overwrite_patch_rasters) {
    patch_abort(paste0("Output patch raster exists and overwrite_patch_rasters is FALSE: ", out_file))
  }

  log_patch_step(species_name, "WRITE_PATCH_RASTER")
  terra::writeRaster(
    components$patch_final,
    filename = out_file,
    overwrite = run_options$overwrite_patch_rasters,
    wopt = list(
      datatype = "INT4U",
      gdal = c("TILED=YES", "BLOCKXSIZE=512", "BLOCKYSIZE=512", "COMPRESS=LZW", "BIGTIFF=IF_SAFER")
    )
  )

  components$patch_final <- NULL
  components$class_lc <- row$class_lc
  components$sdm_method <- row$sdm_method
  components$patch_raster <- out_file
  components
}

patch_status_summary <- function(results) {
  dplyr::bind_rows(lapply(results, function(x) {
    data.frame(
      species = x$species %||% NA_character_,
      class_lc = x$class_lc %||% NA_character_,
      sdm_method = x$sdm_method %||% NA_character_,
      status = x$status %||% NA_character_,
      mapped_habitat_area_km2 = x$mapped_habitat_area_km2 %||% NA_real_,
      n_raw_patches = x$n_raw_patches %||% NA_integer_,
      n_patches_after_patch_filter = x$n_patches_after_patch_filter %||% NA_integer_,
      n_candidate_pus = x$n_candidate_pus %||% NA_integer_,
      n_final_pus = x$n_final_pus %||% NA_integer_,
      n_final_patches = x$n_final_patches %||% NA_integer_,
      patch_raster = x$patch_raster %||% NA_character_,
      stringsAsFactors = FALSE
    )
  }))
}

run_patch_pipeline <- function(selected_species, paths, run_options, habitat_masks) {
  template <- habitat_masks[[1]]
  label_map <- habitat_label_map()
  results <- vector("list", nrow(selected_species))

  for (i in seq_len(nrow(selected_species))) {
    row <- selected_species[i, ]
    results[[i]] <- tryCatch(
      process_patch_species(row, paths, run_options, habitat_masks, template, label_map),
      error = function(e) {
        msg <- paste0("ERROR: ", conditionMessage(e))
        log_patch_step(row$scientificName, msg)
        if (run_options$stop_on_species_error) stop(e)
        list(status = "error", species = row$scientificName, error = conditionMessage(e))
      }
    )
    gc(FALSE)
  }

  retained <- results[vapply(results, function(x) identical(x$status, "retained"), logical(1))]
  all_patch_lookup <- dplyr::bind_rows(lapply(retained, `[[`, "patch_lookup"))
  all_connectivity <- unlist(lapply(retained, `[[`, "connectivity"), recursive = FALSE)

  attr(all_connectivity, "csr_version") <- 2L
  attr(all_connectivity, "row_ptr_base") <- "0-based"
  attr(all_connectivity, "col_idx_space") <- "pu_local_index"
  attr(all_connectivity, "col_idx_base") <- "0-based"

  list(
    results = results,
    status_summary = patch_status_summary(results),
    all_patch_lookup = all_patch_lookup,
    all_connectivity = all_connectivity
  )
}
