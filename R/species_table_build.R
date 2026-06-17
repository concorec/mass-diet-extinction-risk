# Final assembly and validation helpers for species_table.csv.

species_table_keep_cols <- function(curves) {
  c(
    "scientificName", "className", "orderName", "familyName", "genusName", "speciesName",
    "redlistCategory",
    "raster_stem", "name_raster", "match_raster", "sdm_method", "raster_dir", "raster_file", "raster_path",
    "name_trait", "match_trait",
    "BodyMass.Value", "Diet", "diet_combo",
    "density", "home_range_size", "dispersal_dist", "min_patch_size", "min_pop_size",
    as.vector(rbind(paste0("alpha_", curves), paste0("beta_", curves))),
    "habitat_codes_suitable", "habitats_level1", "habitats_mixed"
  )
}

validate_species_table <- function(species_table, curves, use_gompertz_mammals, use_gompertz_birds) {
  if (!nrow(species_table)) {
    species_abort("No species rows remained after resolving rasters and trait records.")
  }

  required_positive <- c("BodyMass.Value", "density", "home_range_size", "dispersal_dist", "min_patch_size", "min_pop_size")
  for (col in required_positive) {
    bad <- !is.finite(species_table[[col]]) | species_table[[col]] <= 0
    if (any(bad, na.rm = TRUE)) species_abort(col, " contains missing, non-finite, or non-positive values.")
  }

  if (any(is.na(species_table$habitats_mixed) | !nzchar(species_table$habitats_mixed))) {
    species_abort("habitats_mixed contains missing or empty values; downstream habitat masking requires at least one suitable habitat label per retained species.")
  }

  if (isTRUE(use_gompertz_mammals) || isTRUE(use_gompertz_birds)) {
    for (curve in curves) {
      for (prefix in c("alpha_", "beta_")) {
        col <- paste0(prefix, curve)
        bad <- !is.finite(species_table[[col]]) | species_table[[col]] <= 0
        if (any(bad, na.rm = TRUE)) species_abort(col, " contains missing, non-finite, or non-positive values.")
      }
    }
  }

  invisible(TRUE)
}

build_species_table <- function(paths, curves, use_gompertz_mammals, use_gompertz_birds,
                                iucn_pause_seconds) {
  core <- read_species_core_inputs(paths, use_gompertz_mammals, use_gompertz_birds)
  resolved <- resolve_species_inputs(core)

  species_inputs <- add_trait_covariates(
    resolved$species_inputs,
    random_effects = core$random_effects,
    models = core$models,
    curves = curves,
    use_gompertz_mammals = use_gompertz_mammals,
    use_gompertz_birds = use_gompertz_birds
  )

  habitats <- query_species_habitats(species_inputs, iucn_pause_seconds)

  species_inputs <- species_inputs |>
    dplyr::left_join(habitats, by = c("genusName", "speciesName"))

  species_table <- species_inputs |>
    dplyr::select(dplyr::any_of(species_table_keep_cols(curves)))

  validate_species_table(species_table, curves, use_gompertz_mammals, use_gompertz_birds)

  species_table
}
