# IUCN suitable-habitat query helpers.

iucn_level1 <- c(
  `1` = "Forest", `2` = "Savanna", `3` = "Shrubland", `4` = "Grassland",
  `5` = "Wetlands (inland)", `6` = "Rocky Areas", `7` = "Caves & Subterranean (non-aquatic)",
  `8` = "Desert", `9` = "Marine Neritic", `10` = "Marine Oceanic", `11` = "Marine Deep Ocean Floor",
  `12` = "Marine Intertidal", `13` = "Marine Coastal/Supratidal", `14` = "Artificial - Terrestrial",
  `15` = "Artificial - Aquatic", `16` = "Introduced Vegetation", `17` = "Other", `18` = "Unknown"
)

normalize_habitat_code <- function(x) gsub("_", ".", trimws(as.character(x)), fixed = TRUE)

summarize_artificial_terrestrial <- function(codes) {
  codes_14 <- codes[stringr::str_detect(codes, "^14(\\.|$)")]
  if (!length(codes_14)) return(character(0))
  level2 <- suppressWarnings(as.integer(stringr::str_match(codes_14, "^14\\.(\\d+)")[, 2]))
  level2 <- unique(level2[!is.na(level2)])

  out <- character(0)
  if (any(level2 %in% c(1, 2))) out <- c(out, "Arable & Pastureland")
  if (any(level2 %in% c(3, 6))) out <- c(out, "Plantations & Heavily Degraded Former Forest")
  if (any(level2 %in% c(4, 5))) out <- c(out, "Urban & Rural Gardens")
  unique(out)
}

empty_habitat_row <- function() {
  tibble::tibble(
    habitat_codes_suitable = NA_character_,
    habitats_level1 = NA_character_,
    habitats_mixed = NA_character_
  )
}

iucn_log <- function(...) {
  message(
    "[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] IUCN | ",
    paste0(..., collapse = "")
  )
}

query_iucn_habitats <- function(genus, species) {
  res <- tryCatch(
    rredlist::rl_species_latest(genus = genus, species = species, scope = "1", parse = TRUE),
    error = function(e) species_abort("IUCN API call failed for ", genus, " ", species, ": ", conditionMessage(e))
  )

  habitats <- res$habitats
  if (is.null(habitats) || !nrow(habitats)) return(empty_habitat_row())

  suitable <- habitats |>
    dplyr::mutate(suitability = stringr::str_to_lower(stringr::str_squish(as.character(suitability)))) |>
    dplyr::filter(suitability == "suitable")

  if (!nrow(suitable)) return(empty_habitat_row())

  codes <- sort(unique(normalize_habitat_code(suitable$code)))
  codes <- codes[!is.na(codes) & nzchar(codes)]
  level1_id <- suppressWarnings(as.integer(stringr::str_extract(codes, "^\\d+")))
  level1_labels <- unique(unname(iucn_level1[as.character(level1_id)]))
  level1_labels <- level1_labels[!is.na(level1_labels)]
  mixed_labels <- unique(c(level1_labels, summarize_artificial_terrestrial(codes)))

  tibble::tibble(
    habitat_codes_suitable = if (length(codes)) paste(codes, collapse = ", ") else NA_character_,
    habitats_level1 = if (length(level1_labels)) paste(level1_labels, collapse = ", ") else NA_character_,
    habitats_mixed = if (length(mixed_labels)) paste(mixed_labels, collapse = ", ") else NA_character_
  )
}

query_species_habitats <- function(species_inputs, pause_seconds) {
  query_species <- species_inputs |>
    dplyr::distinct(genusName, speciesName) |>
    dplyr::filter(!is.na(genusName), nzchar(genusName), !is.na(speciesName), nzchar(speciesName))

  dplyr::bind_rows(Map(
    f = function(genus, species, idx) {
      genus <- stringr::str_squish(genus)
      species <- stringr::str_squish(species)
      species_label <- paste(genus, species)

      t0 <- Sys.time()

      out <- query_iucn_habitats(stringr::str_squish(genus), stringr::str_squish(species))
      elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 2)

      iucn_log(
        idx, "/", nrow(query_species), " | ", species_label,
        " | elapsed=", elapsed, "s",
        " | suitable_codes=", dplyr::coalesce(out$habitat_codes_suitable[[1]], "none")
      )

      if (pause_seconds > 0) {
        Sys.sleep(pause_seconds)
      }

      dplyr::mutate(out, genusName = genus, speciesName = species)
    },
    query_species$genusName,
    query_species$speciesName,
    seq_len(nrow(query_species))
  ))
}
