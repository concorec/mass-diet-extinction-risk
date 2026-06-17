# Input readers and name-resolution helpers for species table construction.

clean_name <- function(x) stringr::str_squish(stringr::str_to_lower(as.character(x)))

as_number <- function(x) suppressWarnings(readr::parse_number(as.character(x)))

species_to_stem <- function(x) {
  paste0(
    stringr::str_replace_all(stringr::str_squish(stringr::str_to_lower(x)), "\\s+", "_"),
    "_bin"
  )
}

read_synonyms <- function(path) {
  tryCatch(
    readr::read_csv(path, show_col_types = FALSE),
    error = function(e) readr::read_delim(path, delim = "\t", show_col_types = FALSE)
  ) |>
    dplyr::transmute(
      scientificName = stringr::str_squish(scientificName),
      synonym = stringr::str_squish(paste(genusName, speciesName))
    ) |>
    dplyr::filter(nzchar(scientificName), nzchar(synonym)) |>
    dplyr::distinct()
}

read_traits <- function(path, taxon_class) {
  readr::read_tsv(
    path,
    col_types = readr::cols(.default = readr::col_character()),
    na = c("", "NA", "NaN", "NULL"),
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    dplyr::mutate(
      taxon_class = taxon_class,
      trait_key = clean_name(Scientific),
      BodyMass.Value = as_number(`BodyMass-Value`)
    )
}

list_raster_manifest <- function(raster_sources) {
  out <- lapply(seq_len(nrow(raster_sources)), function(i) {
    src <- raster_sources[i, , drop = FALSE]
    files <- list.files(src$dir_path, full.names = TRUE)
    if (!length(files)) return(NULL)

    tibble::tibble(
      taxon_class = src$taxon_class,
      sdm_method = src$sdm_method,
      raster_dir = src$dir_path,
      raster_path = files,
      raster_file = basename(files),
      raster_stem = stringr::str_to_lower(tools::file_path_sans_ext(basename(files))),
      sdm_priority = ifelse(src$sdm_method == "PPM", 0L, 1L)
    ) |>
      dplyr::filter(stringr::str_ends(raster_stem, "_bin"))
  })

  dplyr::bind_rows(out) |>
    dplyr::arrange(taxon_class, sdm_priority, raster_file) |>
    dplyr::distinct(taxon_class, sdm_method, raster_stem, .keep_all = TRUE)
}

read_random_effects <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) |>
    dplyr::mutate(
      Class = clean_name(Class),
      Order = clean_name(Order),
      Family = clean_name(Family),
      Species = clean_name(Species)
    )
}

read_species_core_inputs <- function(paths, use_gompertz_mammals, use_gompertz_birds) {
  summary_df <- readr::read_csv(paths$raw$summary, show_col_types = FALSE, progress = FALSE) |>
    dplyr::mutate(
      row_id = dplyr::row_number(),
      scientificName = stringr::str_squish(scientificName),
      taxon_class = dplyr::case_when(
        clean_name(className) == "mammalia" ~ "Mammalia",
        clean_name(className) == "aves" ~ "Aves",
        TRUE ~ NA_character_
      )
    )

  list(
    summary_df = summary_df,
    synonyms = read_synonyms(paths$raw$synonyms),
    raster_manifest = list_raster_manifest(paths$rasters),
    traits = dplyr::bind_rows(
      read_traits(paths$raw$mammal_traits, "Mammalia"),
      read_traits(paths$raw$bird_traits, "Aves")
    ),
    random_effects = read_random_effects(paths$raw$random_effects),
    models = if (isTRUE(use_gompertz_mammals) || isTRUE(use_gompertz_birds)) {
      readRDS(paths$clean$gompertz_models)
    } else NULL
  )
}

build_name_candidates <- function(summary_df, synonyms, traits) {
  trait_index <- unique(paste(traits$taxon_class, traits$trait_key, sep = "||"))

  dplyr::bind_rows(
    summary_df |>
      dplyr::transmute(row_id, taxon_class, candidate_name = scientificName, match_type = "original"),
    summary_df |>
      dplyr::select(row_id, taxon_class, scientificName) |>
      dplyr::left_join(synonyms, by = "scientificName") |>
      dplyr::transmute(row_id, taxon_class, candidate_name = synonym, match_type = "synonym")
  ) |>
    dplyr::filter(!is.na(taxon_class), !is.na(candidate_name), nzchar(candidate_name)) |>
    dplyr::mutate(
      candidate_stem = species_to_stem(candidate_name),
      candidate_key = clean_name(candidate_name),
      name_priority = ifelse(match_type == "original", 0L, 1L),
      trait_available = paste(taxon_class, candidate_key, sep = "||") %in% trait_index
    )
}

pick_rasters <- function(name_candidates, raster_manifest) {
  name_candidates |>
    dplyr::left_join(
      dplyr::select(
        raster_manifest,
        taxon_class, raster_stem, sdm_method, raster_dir, raster_path, raster_file, sdm_priority
      ),
      by = c("taxon_class", "candidate_stem" = "raster_stem")
    ) |>
    dplyr::filter(!is.na(raster_path)) |>
    dplyr::group_by(row_id) |>
    dplyr::arrange(sdm_priority, name_priority, candidate_name, .by_group = TRUE) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      row_id,
      name_raster = candidate_name,
      match_raster = match_type,
      raster_stem = candidate_stem,
      sdm_method,
      raster_dir,
      raster_file,
      raster_path
    )
}

pick_traits <- function(name_candidates) {
  name_candidates |>
    dplyr::filter(trait_available) |>
    dplyr::group_by(row_id) |>
    dplyr::arrange(name_priority, candidate_name, .by_group = TRUE) |>
    dplyr::slice(1) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      row_id,
      name_trait = candidate_name,
      match_trait = match_type,
      trait_key = candidate_key
    )
}

resolve_species_inputs <- function(core) {
  name_candidates <- build_name_candidates(core$summary_df, core$synonyms, core$traits)
  raster_pick <- pick_rasters(name_candidates, core$raster_manifest)
  trait_pick <- pick_traits(name_candidates)

  species_inputs <- core$summary_df |>
    dplyr::left_join(raster_pick, by = "row_id") |>
    dplyr::left_join(trait_pick, by = "row_id") |>
    dplyr::filter(!is.na(raster_path), !is.na(trait_key)) |>
    dplyr::left_join(core$traits, by = c("taxon_class", "trait_key")) |>
    dplyr::select(-row_id)

  list(
    species_inputs = species_inputs,
    name_candidates = name_candidates,
    raster_pick = raster_pick,
    trait_pick = trait_pick
  )
}
