landcover_class_map <- function() {
  class_map <- list(
    forest                      = c(50L, 60L, 61L, 62L, 70L, 71L, 72L, 90L, 160L),
    savanna                     = c(120L, 121L, 122L),
    shrubland                   = c(120L, 121L, 122L, 150L, 151L, 152L, 153L, 200L, 201L, 202L),
    grassland                   = c(130L, 140L, 150L, 151L, 152L, 153L),
    wetlands_inland             = c(20L, 80L, 81L, 82L, 160L, 170L, 180L, 210L),
    rocky_areas                 = c(70L, 71L, 72L, 130L, 150L, 151L, 152L, 153L, 200L, 201L, 202L),
    desert                      = c(150L, 151L, 152L, 153L, 200L, 201L, 202L),
    arable_pastureland          = c(11L, 20L, 190L),
    plantations_degraded_forest = integer(0),
    urban_rural_gardens         = c(190L),
    artificial_aquatic          = c(20L, 160L, 170L, 180L, 190L, 210L)
  )
  class_map$artificial_terrestrial <- unique(c(
    class_map$arable_pastureland,
    class_map$urban_rural_gardens,
    class_map$plantations_degraded_forest
  ))
  class_map
}

habitat_label_map <- function() {
  c(
    "Forest"                                       = "forest",
    "Savanna"                                      = "savanna",
    "Shrubland"                                    = "shrubland",
    "Grassland"                                    = "grassland",
    "Wetlands (inland)"                            = "wetlands_inland",
    "Rocky Areas"                                  = "rocky_areas",
    "Desert"                                       = "desert",
    "Arable & Pastureland"                         = "arable_pastureland",
    "Plantations & Heavily Degraded Former Forest" = "plantations_degraded_forest",
    "Urban & Rural Gardens"                        = "urban_rural_gardens",
    "Artificial - Aquatic"                         = "artificial_aquatic",
    "Artificial - Terrestrial"                     = "artificial_terrestrial"
  )
}

parse_habitats_mixed <- function(x) {
  out <- unlist(strsplit(as.character(x), ",", fixed = TRUE))
  out <- stringr::str_squish(out)
  out <- out[nzchar(out)]
  unique(out)
}

build_habitat_masks_from_landcover <- function(landcover, class_map = landcover_class_map()) {
  base <- landcover[[1]]
  layers <- lapply(names(class_map), function(mask_name) {
    codes <- class_map[[mask_name]]
    r <- if (length(codes) == 0L || all(is.na(codes))) {
      base * NA_integer_
    } else {
      terra::ifel(base %in% codes, 1L, NA_integer_)
    }
    names(r) <- mask_name
    r
  })
  terra::rast(layers)
}

load_habitat_masks <- function(landcover_tif, roi) {
  landcover <- terra::rast(landcover_tif)
  on.exit(rm(landcover), add = TRUE)
  landcover_roi <- terra::crop(landcover, roi, snap = "out")
  build_habitat_masks_from_landcover(landcover_roi, landcover_class_map())
}

load_presence_aligned01 <- function(path, template) {
  r <- terra::rast(path)
  if (!terra::compareGeom(r, template, stopOnError = FALSE)) {
    if (!terra::same.crs(r, template)) {
      r <- terra::project(r, template, method = "near")
    }
    r <- terra::resample(r, template, method = "near")
  }
  terra::ifel(r == 1L, 1L, NA_integer_)
}

habitat_union <- function(masks, mask_names) {
  x <- masks[[mask_names]]
  if (terra::nlyr(x) == 1L) return(x)

  terra::app(x, fun = function(v) {
    if (is.matrix(v)) {
      out <- rep.int(NA_integer_, nrow(v))
      out[rowSums(v == 1L, na.rm = TRUE) > 0] <- 1L
      out
    } else {
      if (any(v == 1L, na.rm = TRUE)) 1L else NA_integer_
    }
  })
}

mapped_habitat_for_species <- function(row, habitat_masks, template, label_map = habitat_label_map()) {
  habitat_labels <- parse_habitats_mixed(row$habitats_mixed)
  mask_names <- unname(label_map[habitat_labels])
  mask_names <- unique(mask_names[!is.na(mask_names)])
  mask_names <- mask_names[mask_names %in% names(habitat_masks)]

  if (!length(mask_names)) {
    return(list(status = "skipped_no_mapped_habitat_labels", labels = habitat_labels, mask_names = mask_names))
  }

  habitat01 <- habitat_union(habitat_masks, mask_names)
  presence01 <- load_presence_aligned01(row$raster_path, template)
  mapped_habitat <- terra::ifel(habitat01 == 1L & presence01 == 1L, 1L, NA_integer_)

  list(
    status = "ok",
    labels = habitat_labels,
    mask_names = mask_names,
    habitat01 = habitat01,
    presence01 = presence01,
    mapped_habitat = mapped_habitat
  )
}

has_non_na <- function(r) {
  x <- terra::global(!is.na(r), "sum", na.rm = TRUE)[1, 1]
  is.finite(x) && x > 0
}

raster_area_km2 <- function(r01) {
  cell_area <- terra::cellSize(r01, unit = "km")
  area <- terra::global(terra::ifel(is.na(r01), NA_real_, cell_area), "sum", na.rm = TRUE)[1, 1]
  as.numeric(area)
}

bin01 <- function(r) terra::ifel(r == 1L, 1L, NA_integer_)
