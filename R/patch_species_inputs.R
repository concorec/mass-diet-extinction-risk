required_patch_species_cols <- function(extra = character()) {
  unique(c(
    "scientificName",
    "className",
    "habitats_mixed",
    "raster_path",
    "sdm_method",
    "raster_dir",
    "raster_file",
    "density",
    "min_patch_size",
    "min_pop_size",
    "dispersal_dist",
    extra
  ))
}

read_patch_species_base <- function(species_csv) {
  required_cols <- required_patch_species_cols()
  species_header <- readr::read_csv(species_csv, show_col_types = FALSE, n_max = 0)
  need_cols(species_header, required_cols, "species_table.csv")

  readr::read_csv(
    species_csv,
    show_col_types = FALSE,
    progress = FALSE,
    col_select = dplyr::all_of(required_cols)
  ) |>
    dplyr::mutate(
      scientificName = stringr::str_squish(as.character(scientificName)),
      class_lc       = stringr::str_to_lower(stringr::str_squish(as.character(className))),
      habitats_mixed = stringr::str_squish(as.character(habitats_mixed)),
      raster_path    = stringr::str_squish(as.character(raster_path)),
      raster_dir     = stringr::str_squish(as.character(raster_dir)),
      raster_file    = stringr::str_squish(as.character(raster_file)),
      sdm_method     = stringr::str_squish(as.character(sdm_method)),
      disp_km        = suppressWarnings(readr::parse_number(as.character(dispersal_dist)))
    ) |>
    add_canonical_area_thresholds(
      validate = TRUE,
      label = "stage-5 species-table area thresholds"
    ) |>
    dplyr::mutate(
      min_patch_km2 = min_patch_area_km2,
      min_pop_km2 = min_population_area_km2
    )
}

read_patch_species_table <- function(species_csv, do_mammals = TRUE, do_birds = TRUE) {
  species_table <- read_patch_species_base(species_csv)

  selected_species <- species_table |>
    dplyr::filter(
      (class_lc == "mammalia" & isTRUE(do_mammals)) |
        (class_lc == "aves" & isTRUE(do_birds))
    )

  validate_selected_patch_species(selected_species)
  selected_species
}

read_single_patch_species <- function(species_csv, target_species) {
  species_table <- read_patch_species_base(species_csv)

  target_row <- species_table |>
    dplyr::filter(scientificName == target_species) |>
    dplyr::slice(1)

  assert(nrow(target_row) == 1L, paste0("Target species not found in species_table.csv: ", target_species))
  validate_selected_patch_species(target_row)
  target_row
}

validate_selected_patch_species <- function(selected_species) {
  assert(nrow(selected_species) > 0L, "No species selected after applying filters.")
  assert(all(nzchar(selected_species$scientificName)), "Selected species contain blank scientificName values.")
  assert(all(nzchar(selected_species$habitats_mixed)), "Selected species contain blank habitats_mixed values.")
  assert(all(nzchar(selected_species$raster_path)), "Selected species contain blank raster_path values.")
  assert(
    all(file.exists(selected_species$raster_path)),
    paste0(
      "At least one selected raster_path does not exist. First missing paths:\n",
      paste(head(selected_species$raster_path[!file.exists(selected_species$raster_path)], 20), collapse = "\n")
    )
  )
  assert(
    all(is.finite(selected_species$min_patch_km2) & selected_species$min_patch_km2 > 0),
    "Selected species contain non-positive or non-numeric min_patch_size values."
  )
  assert(
    all(is.finite(selected_species$min_pop_km2) & selected_species$min_pop_km2 > 0),
    "Selected species contain non-positive or non-numeric min_pop_size values."
  )
  assert(
    all(is.finite(selected_species$disp_km) & selected_species$disp_km >= 0),
    "Selected species contain negative or non-numeric dispersal_dist values."
  )

  patch_files <- patch_filename_from_scientific(selected_species$scientificName)
  duplicate_patch_files <- unique(patch_files[duplicated(patch_files)])
  assert(
    length(duplicate_patch_files) == 0L,
    paste0("Duplicate patch-raster filenames would be created: ", paste(duplicate_patch_files, collapse = ", "))
  )

  invisible(TRUE)
}
