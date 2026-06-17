# Configuration and precondition helpers for 4_build_species_table.Rmd.

source(file.path("R", "project_utils.R"))
source(file.path("R", "analysis_contract.R"))

species_abort <- function(...) stop(paste0(...), call. = FALSE)

species_table_curves <- function() {
  persistence_curves()
}

required_species_table_packages <- function() {
  c("readr", "dplyr", "stringr", "tibble", "tools", "rredlist")
}

load_species_table_packages <- function() {
  load_packages(required_species_table_packages())
}

trait_header <- function(path, label) {
  hdr <- trimws(strsplit(readLines(path, n = 1, warn = FALSE), "\t", fixed = TRUE)[[1]])
  missing <- setdiff(c("Scientific", "BodyMass-Value"), hdr)
  if (length(missing)) species_abort(label, " missing column(s): ", paste(missing, collapse = ", "))
  hdr
}

validate_species_table_config <- function(params, paths, use_gompertz_mammals,
                                          use_gompertz_birds, write_output,
                                          iucn_pause_seconds) {
  if (!is.finite(iucn_pause_seconds) || iucn_pause_seconds < 0) {
    species_abort("params$iucn_pause_seconds must be a non-negative number.")
  }
  if (!nzchar(Sys.getenv("IUCN_REDLIST_KEY"))) {
    species_abort("IUCN_REDLIST_KEY is not set in the environment; it is required for habitat queries.")
  }
  if (isTRUE(write_output)) {
    dir.create(dirname(paths$clean$species_table), recursive = TRUE, showWarnings = FALSE)
  }

  lapply(paths$raw, need_file)
  for (d in paths$rasters$dir_path) need_dir(d)
  if (isTRUE(use_gompertz_mammals) || isTRUE(use_gompertz_birds)) {
    need_file(paths$clean$gompertz_models)
  }

  summary_hdr <- readr::read_csv(paths$raw$summary, n_max = 0, show_col_types = FALSE)
  need_cols(
    summary_hdr,
    c("scientificName", "className", "orderName", "familyName", "genusName", "speciesName", "redlistCategory"),
    "simple_summary.csv"
  )

  syn_hdr <- tryCatch(
    readr::read_csv(paths$raw$synonyms, n_max = 0, show_col_types = FALSE),
    error = function(e) readr::read_delim(paths$raw$synonyms, delim = "\t", n_max = 0, show_col_types = FALSE)
  )
  need_cols(syn_hdr, c("scientificName", "genusName", "speciesName"), "synonyms.csv")

  headers <- list(
    mammal_traits = trait_header(paths$raw$mammal_traits, "mammal_data.txt"),
    bird_traits   = trait_header(paths$raw$bird_traits, "bird_data.txt")
  )
  if (!any(startsWith(headers$bird_traits, "Diet-"))) {
    species_abort("bird_data.txt must contain Diet-* columns.")
  }

  random_effects_hdr <- readr::read_csv(paths$raw$random_effects, n_max = 0, show_col_types = FALSE)
  need_cols(
    random_effects_hdr,
    c("Class", "Order", "Family", "Species", "Effect_Order", "Effect_Family", "Effect_Species"),
    "random_effects.csv"
  )

  invisible(TRUE)
}
