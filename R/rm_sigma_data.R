# Data preparation helpers for demographic calibration models.

MAMMAL_MARINE_RMAX_SPECIES <- c(
  "Phoca fasciata", "Phoca groenlandica", "Balaena mysticetus", "Eubalaena glacialis",
  "Balaenoptera acutorostrata", "Balaenoptera borealis", "Balaenoptera musculus",
  "Balaenoptera physalus", "Megaptera novaeangliae", "Globicephala melas",
  "Lagenorhynchus acutus", "Orcinus orca", "Stenella attenuata", "Stenella coeruleoalba",
  "Stenella longirostris", "Tursiops truncatus", "Eschrichtius robustus",
  "Delphinapterus leucas", "Monodon monoceros", "Phocoena phocoena",
  "Phocoenoides dalli", "Physeter catodon", "Pontoporia blainvillei", "Berardius bairdii",
  "Trichechus manatus"
)

MAMMAL_MARINE_SIGMA_SPECIES <- c(
  "Phoca groenlandica", "Stenella attenuata", "Eschrichtius robustus", "Trichechus manatus"
)

bird_lambda <- tibble::tibble(
  Species = c(
    "Fulmarus glacialis", "Fratercula arctica", "Gyps fulvus", "Rissa tridactyla",
    "Larus argentatus", "Chen caerulescens", "Branta leucopsis", "Phalacrocorax carbo",
    "Larus ridibundus", "Ciconia ciconia", "Sterna caspia", "Parus major", "Petronia petronia"
  ),
  lambda = c(1.06, 1.09, 1.09, 1.12, 1.13, 1.17, 1.18, 1.19, 1.19, 1.21, 1.29, 1.64, 2.15)
)

as_num <- function(x) {
  if (is.numeric(x)) return(x)
  suppressWarnings(readr::parse_number(x))
}

read_calibration_inputs <- function(paths) {
  raw_mammal_rmax <- readr::read_tsv(paths$raw$mammal_rmax, show_col_types = FALSE)
  raw_bird_traits <- readr::read_tsv(
    paths$raw$bird_traits,
    show_col_types = FALSE,
    col_types = readr::cols(.default = readr::col_character())
  )
  raw_sigma <- readr::read_csv(
    paths$raw$sigma,
    show_col_types = FALSE,
    col_select = dplyr::any_of(c("Class", "Genus", "Species", "Vr", "Mass"))
  )
  raw_synonyms <- readr::read_csv(paths$raw$bird_synonyms, show_col_types = FALSE)

  validate_raw_inputs(raw_mammal_rmax, raw_bird_traits, raw_sigma, raw_synonyms)

  list(
    mammal_rmax = raw_mammal_rmax,
    bird_traits = raw_bird_traits,
    sigma = raw_sigma,
    synonyms = raw_synonyms
  )
}

validate_raw_inputs <- function(raw_mammal_rmax, raw_bird_traits, raw_sigma, raw_synonyms) {
  assert_has_cols(raw_mammal_rmax, c("log10(M)", "log10(rm)", "Order", "Species"), "mammal_rmax.txt")
  assert_has_cols(
    raw_bird_traits,
    c(
      "Scientific", "BodyMass-Value", "Diet-5Cat",
      "Diet-Vend", "Diet-Vect", "Diet-Vfish", "Diet-Vunk", "Diet-Scav",
      "Diet-Fruit", "Diet-Nect", "Diet-PlantO", "Diet-Seed", "Diet-Inv"
    ),
    "bird_data.txt"
  )
  assert_has_cols(raw_sigma, c("Genus", "Species", "Class", "Mass", "Vr"), "sigma.csv")
  assert_has_cols(raw_synonyms, c("scientificName", "EltonTraits_scientificName"), "bird_synonyms.csv")

  if (anyDuplicated(raw_synonyms$scientificName) > 0) {
    stop("bird_synonyms.csv must map each scientificName at most once.", call. = FALSE)
  }

  invisible(TRUE)
}

prepare_bird_traits <- function(raw_bird_traits) {
  diet_cols <- c(
    "Diet-Inv", "Diet-Vend", "Diet-Vect", "Diet-Vfish", "Diet-Vunk", "Diet-Scav",
    "Diet-Fruit", "Diet-Nect", "Diet-PlantO", "Diet-Seed"
  )

  bird_traits_num <- raw_bird_traits |>
    dplyr::transmute(
      trait_scientific = Scientific,
      Mass_g = as_num(`BodyMass-Value`),
      Diet_5Cat = `Diet-5Cat`,
      dplyr::across(dplyr::all_of(diet_cols), as_num)
    )

  mass_parse_rate <- mean(is.finite(bird_traits_num$Mass_g) & bird_traits_num$Mass_g > 0, na.rm = TRUE)
  if (mass_parse_rate < 0.90) {
    stop("Fewer than 90% of bird BodyMass-Value records parse as positive numeric values.", call. = FALSE)
  }

  bad_diet_cols <- diet_cols[vapply(diet_cols, function(nm) {
    x <- bird_traits_num[[nm]]
    any(is.finite(x) & (x < 0 | x > 100))
  }, logical(1))]
  if (length(bad_diet_cols) > 0) {
    stop("Bird diet percentages must be between 0 and 100. Problem column(s): ",
         paste(bad_diet_cols, collapse = ", "), call. = FALSE)
  }

  bird_traits_num
}

build_bird_sigma_join <- function(bird_sigma, raw_synonyms, bird_traits_num) {
  synonyms <- raw_synonyms |>
    dplyr::transmute(
      sigma_scientific = scientificName,
      trait_scientific = EltonTraits_scientificName
    )

  traits <- bird_traits_num |>
    dplyr::transmute(
      trait_scientific, Mass_g, Diet_5Cat,
      Diet_Inv = `Diet-Inv`,
      Diet_Vend = `Diet-Vend`,
      Diet_Vect = `Diet-Vect`,
      Diet_Vfish = `Diet-Vfish`,
      Diet_Vunk = `Diet-Vunk`,
      Diet_Scav = `Diet-Scav`,
      Diet_Fruit = `Diet-Fruit`,
      Diet_Nect = `Diet-Nect`,
      Diet_PlantO = `Diet-PlantO`,
      Diet_Seed = `Diet-Seed`
    )

  bird_sigma |>
    dplyr::left_join(synonyms, by = "sigma_scientific") |>
    dplyr::mutate(trait_scientific = dplyr::coalesce(trait_scientific, sigma_scientific)) |>
    dplyr::inner_join(traits, by = "trait_scientific")
}

prepare_calibration_data <- function(inputs) {
  raw_sigma <- inputs$sigma |>
    dplyr::mutate(Mass = as_num(Mass), Vr = as_num(Vr)) |>
    dplyr::filter(!is.na(Mass), !is.na(Vr))

  if (!all(is.finite(raw_sigma$Mass) & raw_sigma$Mass > 0 &
           is.finite(raw_sigma$Vr) & raw_sigma$Vr > 0)) {
    stop("sigma.csv Mass and Vr must parse as finite positive values.", call. = FALSE)
  }

  bird_traits_num <- prepare_bird_traits(inputs$bird_traits)

  mammal_rmax <- inputs$mammal_rmax |>
    dplyr::rename(log10_M = `log10(M)`, log10_rm = `log10(rm)`) |>
    dplyr::filter(
      Order != "Chiroptera",
      !(Species %in% MAMMAL_MARINE_RMAX_SPECIES),
      is.finite(log10_M), is.finite(log10_rm)
    ) |>
    dplyr::transmute(Mass_g = 10^log10_M, rm = 10^log10_rm)

  sigma_std <- raw_sigma |>
    dplyr::transmute(
      Class = Class,
      Species = paste(Genus, Species),
      Mass_g = Mass,
      sigma = sqrt(Vr)
    )

  mammal_sigma <- sigma_std |>
    dplyr::filter(Class == "Mammalia", !(Species %in% MAMMAL_MARINE_SIGMA_SPECIES)) |>
    dplyr::transmute(Mass_g, sigma)

  bird_sigma <- sigma_std |>
    dplyr::filter(Class == "Aves") |>
    dplyr::transmute(sigma_scientific = Species, sigma)

  bird_mass <- bird_traits_num |>
    dplyr::transmute(Species = trait_scientific, Mass_g = Mass_g) |>
    dplyr::filter(is.finite(Mass_g), Mass_g > 0)

  bird_rmax <- bird_lambda |>
    dplyr::mutate(rm = log(lambda)) |>
    dplyr::inner_join(bird_mass, by = "Species") |>
    dplyr::transmute(Mass_g, rm)

  bird_sigma_join <- build_bird_sigma_join(bird_sigma, inputs$synonyms, bird_traits_num) |>
    dplyr::filter(is.finite(Mass_g), Mass_g > 0, is.finite(sigma), sigma > 0)

  bird_sigma_data <- bird_sigma_join |>
    dplyr::transmute(
      Mass_g, sigma, Diet_5Cat,
      Diet_Inv = Diet_Inv,
      Diet_VertFishScav = Diet_Vend + Diet_Vect + Diet_Vfish + Diet_Vunk + Diet_Scav,
      Diet_AllPlants = Diet_Fruit + Diet_Nect + Diet_PlantO + Diet_Seed
    )

  bird_sigma_model_compare_data <- bird_sigma_join |>
    dplyr::transmute(
      Mass_g, sigma, Diet_5Cat,
      Diet_Inv = Diet_Inv,
      Diet_VertFishScav = Diet_Vend + Diet_Vect + Diet_Vfish + Diet_Vunk + Diet_Scav,
      Diet_FruiNect = Diet_Fruit + Diet_Nect,
      Diet_PlantSeed = Diet_PlantO + Diet_Seed
    )

  validate_calibration_data(mammal_rmax, mammal_sigma, bird_rmax, bird_sigma_data)

  list(
    mammal_rmax = mammal_rmax,
    mammal_sigma = mammal_sigma,
    bird_rmax = bird_rmax,
    bird_sigma_data = bird_sigma_data,
    bird_sigma_model_compare_data = bird_sigma_model_compare_data
  )
}

validate_calibration_data <- function(mammal_rmax, mammal_sigma, bird_rmax, bird_sigma_data) {
  expected_n <- c(mammal_rmax = 265L, mammal_sigma = 148L, bird_rmax = 13L, bird_sigma = 225L)
  observed_n <- c(
    mammal_rmax = nrow(mammal_rmax),
    mammal_sigma = nrow(mammal_sigma),
    bird_rmax = nrow(bird_rmax),
    bird_sigma = nrow(bird_sigma_data)
  )

  bad_n <- names(observed_n)[observed_n != expected_n]
  if (length(bad_n) > 0) {
    stop(
      "Calibration sample size mismatch: ",
      paste(sprintf("%s expected %d observed %d", bad_n, expected_n[bad_n], observed_n[bad_n]),
            collapse = "; "),
      call. = FALSE
    )
  }

  diet_total <- rowSums(bird_sigma_data[, c("Diet_Inv", "Diet_AllPlants", "Diet_VertFishScav")], na.rm = FALSE)
  if (any(!is.finite(diet_total)) || any(abs(diet_total - 100) > 1e-8)) {
    stop("Grouped bird diet percentages must sum to 100 for every bird sigma calibration record.", call. = FALSE)
  }

  invisible(TRUE)
}
