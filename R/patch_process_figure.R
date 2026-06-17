madagascar_outline <- function(template_rast) {
  mg <- rnaturalearth::ne_countries(
    country = "Madagascar",
    scale = "medium",
    returnclass = "sf"
  )
  sf::st_transform(mg, crs = terra::crs(template_rast, proj = TRUE))
}

tight_extent_from_sf <- function(x, pad_x = 0.035, pad_y = 0.025) {
  bb <- sf::st_bbox(x)
  dx <- unname(bb["xmax"] - bb["xmin"])
  dy <- unname(bb["ymax"] - bb["ymin"])
  terra::ext(
    bb["xmin"] - pad_x * dx,
    bb["xmax"] + pad_x * dx,
    bb["ymin"] - pad_y * dy,
    bb["ymax"] + pad_y * dy
  )
}

fmt_int <- function(x) scales::label_comma(accuracy = 1)(x)
fmt_area_km2 <- function(x) paste0(scales::label_comma(accuracy = 1)(round(x)), " km\u00b2")

fmt_habitat_subtitle <- function(lbls) {
  prefix <- if (length(lbls) == 1L) "Suitable habitat: " else "Suitable habitats: "
  paste0(prefix, paste(unique(lbls), collapse = ", "))
}

fmt_species_subtitle <- function(x) bquote(italic(.(stringr::str_squish(as.character(x)[1]))))

geom_spatraster_silent <- function(...) suppressMessages(tidyterra::geom_spatraster(...))

candidate_pu_pal <- function(n) {
  if (n <= 0) return(character(0))
  grDevices::hcl(
    h = seq(15, 375, length.out = n + 1)[seq_len(n)],
    c = rep(c(38, 32, 35), length.out = n),
    l = rep(c(74, 69, 72), length.out = n),
    fixup = TRUE
  )
}

final_pu_pal <- function(n) {
  if (n <= 0) return(character(0))
  base_pal <- c("#4E79A7", "#59A14F", "#7B6FB2", "#3F7F93")
  if (n <= length(base_pal)) return(base_pal[seq_len(n)])
  c(base_pal, grDevices::hcl.colors(n - length(base_pal), palette = "Dark 3"))
}

theme_map <- ggplot2::theme_void(base_size = 10.5) +
  ggplot2::theme(
    text = ggplot2::element_text(colour = "grey10"),
    plot.title = ggplot2::element_text(hjust = 0.5, size = 11.6, face = "bold", margin = ggplot2::margin(b = 1.5)),
    plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 8.5, colour = "grey30", margin = ggplot2::margin(b = 2.3)),
    plot.title.position = "plot",
    plot.margin = ggplot2::margin(3, 3, 3, 3),
    panel.border = ggplot2::element_rect(color = "grey85", fill = NA, linewidth = 0.30),
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.background = ggplot2::element_rect(fill = "white", colour = NA)
  )

coord_ext <- function(r, e) {
  ggplot2::coord_sf(
    crs = terra::crs(r, proj = TRUE),
    xlim = c(e$xmin, e$xmax),
    ylim = c(e$ymin, e$ymax),
    expand = FALSE
  )
}

add_mg_fill <- function(mg_sf, fill = "white") {
  ggplot2::geom_sf(data = mg_sf, fill = fill, color = NA, inherit.aes = FALSE)
}

add_mg_outline <- function(mg_sf, col = "grey15", lwd = 0.28) {
  ggplot2::geom_sf(data = mg_sf, fill = NA, color = col, linewidth = lwd, inherit.aes = FALSE)
}

bin_plot <- function(r01, title, subtitle, fill_col, mg_sf, e, maxcell = 8e5, land_fill = "white") {
  r01 <- bin01(r01)
  ggplot2::ggplot() +
    add_mg_fill(mg_sf, fill = land_fill) +
    geom_spatraster_silent(data = r01, maxcell = maxcell, na.rm = TRUE) +
    ggplot2::scale_fill_gradient(low = fill_col, high = fill_col, na.value = NA, guide = "none") +
    add_mg_outline(mg_sf) +
    coord_ext(r01, e) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    theme_map
}

patch_filter_plot <- function(status_r, title, subtitle, mg_sf, e, discarded_col, retained_col, maxcell = 8e5, land_fill = "white") {
  status_plot <- terra::as.factor(status_r)
  levels(status_plot) <- data.frame(value = c(1, 2), label = c("Discarded", "Retained"))
  ggplot2::ggplot() +
    add_mg_fill(mg_sf, fill = land_fill) +
    geom_spatraster_silent(data = status_plot, maxcell = maxcell, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = c("Discarded" = discarded_col, "Retained" = retained_col), guide = "none", na.translate = FALSE) +
    add_mg_outline(mg_sf) +
    coord_ext(status_plot, e) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    theme_map
}

unit_plot <- function(r_id, title, subtitle, mg_sf, e, palette, outline_sf = NULL, maxcell = 8e5, land_fill = "white") {
  vals <- sort(unique(terra::values(r_id, mat = FALSE)))
  vals <- vals[!is.na(vals)]
  labs <- paste0("PU_", vals)
  r_plot <- terra::as.factor(r_id)
  if (length(vals)) levels(r_plot) <- data.frame(value = vals, label = labs)

  pal <- if (is.null(names(palette))) {
    stats::setNames(rep(palette, length.out = length(vals)), labs)
  } else {
    palette[labs]
  }
  missing_cols <- is.na(pal)
  if (any(missing_cols)) pal[missing_cols] <- candidate_pu_pal(sum(missing_cols))

  p <- ggplot2::ggplot() +
    add_mg_fill(mg_sf, fill = land_fill) +
    geom_spatraster_silent(data = r_plot, maxcell = maxcell, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = pal, guide = "none", na.translate = FALSE) +
    add_mg_outline(mg_sf) +
    coord_ext(r_plot, e) +
    ggplot2::labs(title = title, subtitle = subtitle) +
    theme_map

  if (!is.null(outline_sf)) {
    p <- p + ggplot2::geom_sf(
      data = outline_sf,
      fill = NA,
      color = scales::alpha("grey10", 0.12),
      linewidth = 0.11,
      inherit.aes = FALSE
    )
  }

  p
}

build_single_species_layers <- function(target_row, paths, roi) {
  habitat_masks <- load_habitat_masks(paths$landcover_tif, roi)
  template <- habitat_masks[[1]]
  mapped <- mapped_habitat_for_species(target_row, habitat_masks, template)
  assert(identical(mapped$status, "ok"), "Target species has no habitat labels represented in the land-cover crosswalk.")
  assert(has_non_na(mapped$mapped_habitat), "Target species has empty mapped habitat after intersecting habitat and distribution.")

  mapped_trim <- terra::trim(mapped$mapped_habitat)
  patch_all <- clump_binary_raster(mapped_trim, diagonal = FALSE)
  names(patch_all) <- "patch_id"

  cell_area <- terra::cellSize(mapped_trim, unit = "km")
  patch_area <- as.data.frame(terra::zonal(cell_area, patch_all, "sum", na.rm = TRUE))
  names(patch_area) <- c("patch_id_raw", "area_km2")
  patch_area <- patch_area[!is.na(patch_area$patch_id_raw), , drop = FALSE]
  n_raw_patches <- nrow(patch_area)

  keep_raw_patches <- patch_area$patch_id_raw[
    area_exceeds_threshold(patch_area$area_km2, target_row$min_patch_km2)
  ]
  n_kept_patches <- length(keep_raw_patches)
  assert(n_kept_patches > 0L, "Target species has no patches above the minimum patch area threshold.")

  patch_kept <- terra::classify(
    patch_all,
    rcl = cbind(as.integer(keep_raw_patches), seq_along(keep_raw_patches)),
    others = NA_integer_
  )
  names(patch_kept) <- "patch_id"

  patch_filter_status <- terra::ifel(!is.na(patch_kept), 2L, terra::ifel(!is.na(patch_all), 1L, NA_integer_))
  names(patch_filter_status) <- "patch_status"

  patch_table <- data.frame(
    patch_id = seq_along(keep_raw_patches),
    area_km2 = patch_area$area_km2[match(keep_raw_patches, patch_area$patch_id_raw)]
  )

  patch_polygons <- sf::st_as_sf(terra::as.polygons(patch_kept, values = TRUE, dissolve = TRUE, na.rm = TRUE))
  value_col <- setdiff(names(patch_polygons), attr(patch_polygons, "sf_column"))[1]
  patch_polygons$patch_id <- as.integer(patch_polygons[[value_col]])
  patch_polygons <- patch_polygons[!is.na(patch_polygons$patch_id), ]

  neighbors <- sf::st_is_within_distance(patch_polygons, patch_polygons, dist = target_row$disp_km * 1000)
  patch_graph <- igraph::graph_from_adj_list(neighbors, mode = "out") |>
    igraph::as_undirected("collapse") |>
    igraph::simplify(remove.multiple = TRUE, remove.loops = TRUE)

  patch_polygons$pu_prefilter <- as.integer(factor(igraph::components(patch_graph)$membership))
  n_candidate_pus <- dplyr::n_distinct(patch_polygons$pu_prefilter)

  patch_area_with_pu <- dplyr::left_join(
    patch_table,
    dplyr::select(sf::st_drop_geometry(patch_polygons), patch_id, pu_prefilter),
    by = "patch_id"
  )

  pu_area <- patch_area_with_pu |>
    dplyr::group_by(pu_prefilter) |>
    dplyr::summarise(pu_area_km2 = sum(area_km2, na.rm = TRUE), .groups = "drop")

  keep_pu <- sort(pu_area$pu_prefilter[
    area_exceeds_threshold(pu_area$pu_area_km2, target_row$min_pop_km2)
  ])
  n_final_pus <- length(keep_pu)
  assert(n_final_pus > 0L, "Target species has no population units above the minimum population-unit area threshold.")

  patch_polygons$pu_final <- match(patch_polygons$pu_prefilter, keep_pu)

  patch_to_candidate_pu <- dplyr::select(sf::st_drop_geometry(patch_polygons), patch_id, pu_prefilter)
  pu_candidate <- terra::classify(
    patch_kept,
    rcl = cbind(as.integer(patch_to_candidate_pu$patch_id), as.integer(patch_to_candidate_pu$pu_prefilter)),
    others = NA_integer_
  )
  names(pu_candidate) <- "pu_prefilter"

  patch_to_final_pu <- dplyr::select(sf::st_drop_geometry(patch_polygons), patch_id, pu_final)
  patch_to_final_pu <- patch_to_final_pu[!is.na(patch_to_final_pu$pu_final), , drop = FALSE]
  pu_final <- terra::classify(
    patch_kept,
    rcl = cbind(as.integer(patch_to_final_pu$patch_id), as.integer(patch_to_final_pu$pu_final)),
    others = NA_integer_
  )
  names(pu_final) <- "pu_id"

  final_units_sf <- patch_polygons |>
    dplyr::filter(!is.na(pu_final)) |>
    dplyr::group_by(pu_final) |>
    dplyr::summarise(geometry = sf::st_union(geometry), .groups = "drop")

  mg_sf <- madagascar_outline(mapped$habitat01)
  plot_ext <- tight_extent_from_sf(mg_sf, pad_x = 0.035, pad_y = 0.025)

  list(
    target_row = target_row,
    labels_mixed = mapped$labels,
    habitat01 = mapped$habitat01,
    presence01 = mapped$presence01,
    mapped_habitat01 = mapped$mapped_habitat,
    patch_filter_status = patch_filter_status,
    pu_candidate = pu_candidate,
    pu_final = pu_final,
    final_units_sf = final_units_sf,
    mg_sf = mg_sf,
    plot_ext = plot_ext,
    n_raw_patches = n_raw_patches,
    n_kept_patches = n_kept_patches,
    n_candidate_pus = n_candidate_pus,
    n_final_pus = n_final_pus,
    keep_pu = keep_pu,
    habitat_c = terra::crop(mapped$habitat01, plot_ext, snap = "out"),
    presence_c = terra::crop(mapped$presence01, plot_ext, snap = "out"),
    mapped_habitat_c = terra::crop(mapped$mapped_habitat, plot_ext, snap = "out"),
    patch_filter_c = terra::crop(patch_filter_status, plot_ext, snap = "out"),
    pu_candidate_c = terra::crop(pu_candidate, plot_ext, snap = "out"),
    pu_final_c = terra::crop(pu_final, plot_ext, snap = "out")
  )
}

make_single_species_process_figure <- function(layers) {
  col_hab <- "#88AE94"
  col_dist <- "#CC86B3"
  col_aoh <- "#335C67"
  col_keep <- "#4C9B7A"
  col_drop <- "grey74"
  land_fill <- "white"

  final_cols <- stats::setNames(final_pu_pal(layers$n_final_pus), paste0("PU_", seq_len(layers$n_final_pus)))
  candidate_cols <- stats::setNames(candidate_pu_pal(layers$n_candidate_pus), paste0("PU_", seq_len(layers$n_candidate_pus)))

  candidate_keep_labels <- paste0("PU_", layers$keep_pu)
  final_keep_labels <- paste0("PU_", seq_len(layers$n_final_pus))
  candidate_cols[candidate_keep_labels] <- unname(final_cols[final_keep_labels])

  panels <- list(
    bin_plot(layers$habitat_c, "Habitat", fmt_habitat_subtitle(layers$labels_mixed), col_hab, layers$mg_sf, layers$plot_ext, land_fill = land_fill),
    bin_plot(layers$presence_c, "Distribution", fmt_species_subtitle(layers$target_row$scientificName), col_dist, layers$mg_sf, layers$plot_ext, land_fill = land_fill),
    bin_plot(layers$mapped_habitat_c, "Mapped habitat", fmt_area_km2(raster_area_km2(layers$mapped_habitat_c)), col_aoh, layers$mg_sf, layers$plot_ext, land_fill = land_fill),
    patch_filter_plot(layers$patch_filter_c, "Filtered patches", paste0("retained patches = ", fmt_int(layers$n_kept_patches)), layers$mg_sf, layers$plot_ext, col_drop, col_keep, land_fill = land_fill),
    unit_plot(layers$pu_candidate_c, "Candidate population units", paste0("before PU-area filter: n = ", fmt_int(layers$n_candidate_pus)), layers$mg_sf, layers$plot_ext, candidate_cols, land_fill = land_fill),
    unit_plot(layers$pu_final_c, "Final population units", paste0("after PU-area filter: n = ", fmt_int(layers$n_final_pus)), layers$mg_sf, layers$plot_ext, final_cols, outline_sf = layers$final_units_sf, land_fill = land_fill)
  )

  cowplot::plot_grid(
    plotlist = panels,
    ncol = 3,
    labels = c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"),
    label_fontface = "bold",
    label_size = 12,
    label_x = 0.028,
    label_y = 0.978,
    hjust = 0,
    vjust = 1,
    align = "hv",
    axis = "tblr"
  )
}
