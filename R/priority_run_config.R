# Canonical configuration helpers for scripts 6 and 7.

if (!exists("assert", mode = "function")) {
  source(file.path("R", "project_utils.R"))
}
source(file.path("R", "analysis_contract.R"))

priority_flags_from_tags <- function(taxa = "mammals_birds", sdm = "ppm_rangebag") {
  taxa <- as.character(taxa)
  sdm <- as.character(sdm)
  assert(taxa %in% c("mammals", "birds", "mammals_birds"),
         "taxa must be one of: mammals, birds, mammals_birds.")
  assert(sdm %in% c("ppm", "rangebag", "ppm_rangebag"),
         "sdm must be one of: ppm, rangebag, ppm_rangebag.")
  list(
    do_mammals = taxa %in% c("mammals", "mammals_birds"),
    do_birds = taxa %in% c("birds", "mammals_birds"),
    do_ppm = sdm %in% c("ppm", "ppm_rangebag"),
    do_rangebag = sdm %in% c("rangebag", "ppm_rangebag")
  )
}

validate_priority_count <- function(x, label) {
  assert(
    is.numeric(x) && length(x) == 1L && is.finite(x) && x > 0,
    paste0(label, " must be a single positive integer-like value.")
  )
  as.integer(x)
}

validate_max_stages <- function(x, label = "max_stages") {
  assert(
    is.numeric(x) && length(x) == 1L && !is.na(x) && x > 0,
    paste0(label, " must be a single positive number or Inf.")
  )
  if (is.finite(x)) as.integer(x) else Inf
}

priority_selection_from_flags <- function(
  curve_label = "q025",
  do_mammals = TRUE,
  do_birds = TRUE,
  do_ppm = TRUE,
  do_rangebag = TRUE
) {
  curve_label <- validate_persistence_curve(curve_label, "CURVE_LABEL")

  flags <- list(
    do_mammals = isTRUE(do_mammals),
    do_birds = isTRUE(do_birds),
    do_ppm = isTRUE(do_ppm),
    do_rangebag = isTRUE(do_rangebag)
  )

  assert(flags$do_mammals || flags$do_birds,
         "At least one taxonomic gate must be TRUE: DO_MAMMALS and/or DO_BIRDS.")
  assert(flags$do_ppm || flags$do_rangebag,
         "At least one SDM-source gate must be TRUE: DO_PPM and/or DO_RANGEBAG.")

  taxa <- make_selection_tag(c(
    if (flags$do_mammals) "mammals" else "",
    if (flags$do_birds) "birds" else ""
  ))

  sdm <- make_selection_tag(c(
    if (flags$do_ppm) "ppm" else "",
    if (flags$do_rangebag) "rangebag" else ""
  ))

  c(
    list(curve_label = curve_label, taxa = taxa, sdm = sdm),
    flags
  )
}

priority_run_id <- function(
  curve = "q025",
  taxa = "mammals_birds",
  sdm = "ppm_rangebag",
  cells_to_remove_per_iteration = 1000L,
  pruning_iterations_per_stage = 100L
) {
  curve <- validate_persistence_curve(curve)
  sprintf(
    "curve_%s_taxa_%s_sdm_%s_remove_%d_k_%d",
    as.character(curve),
    as.character(taxa),
    as.character(sdm),
    as.integer(cells_to_remove_per_iteration),
    as.integer(pruning_iterations_per_stage)
  )
}

priority_bundle_name <- function(
  curve = "q025",
  taxa = "mammals_birds",
  sdm = "ppm_rangebag"
) {
  curve <- validate_persistence_curve(curve)
  sprintf(
    "pruning_inputs_curve_%s_taxa_%s_sdm_%s.rds",
    as.character(curve),
    as.character(taxa),
    as.character(sdm)
  )
}

benchmark_artifact_paths <- function(
  run_paths,
  rank_method,
  rank_map_subdir = "benchmark_rank_maps"
) {
  rank_method <- validate_benchmark_rank_method(rank_method)
  lut_root <- file.path(run_paths$ana, paste0("rank_lut_", rank_method))
  list(
    rank_method = rank_method,
    rank_label = benchmark_rank_label(rank_method),
    rank_path = file.path(
      run_paths$out_curve,
      rank_map_subdir,
      paste0("rankmap_", rank_method, ".tif")
    ),
    lut_root = lut_root,
    target_summary_csv = file.path(
      lut_root,
      paste0("rank_stage_targets_", rank_method, ".csv")
    )
  )
}

persistence_comparison_paths <- function(run_paths, rank_method) {
  rank_method <- validate_benchmark_rank_method(rank_method)
  out_root <- file.path(run_paths$ana, paste0("persist_cmp_", rank_method))
  list(
    out_root = out_root,
    pu = file.path(out_root, paste0("pu_long_", rank_method, ".csv")),
    species = file.path(out_root, paste0("sp_long_", rank_method, ".csv")),
    stage = file.path(out_root, paste0("stage_sum_", rank_method, ".csv")),
    compare = file.path(out_root, paste0("stage_compare_", rank_method, ".csv"))
  )
}

rank_lut_columns <- function() {
  c(
    "stage", "method", "scientificName", "species", "patch_id", "pu_id",
    "patch_area_km2", "patch_n_cells"
  )
}

empty_rank_lut <- function() {
  data.table::data.table(
    stage = integer(),
    method = character(),
    scientificName = character(),
    species = character(),
    patch_id = integer(),
    pu_id = integer(),
    patch_area_km2 = numeric(),
    patch_n_cells = integer()
  )
}

normalized_rank_lut_columns <- function() {
  c("method", "stage", "scientificName", "species", "patch_id", "pu_id", "patch_area_km2")
}

empty_normalized_rank_lut <- function() {
  empty_rank_lut()[, normalized_rank_lut_columns(), with = FALSE]
}

standardize_rank_lut <- function(x, method, stage, species_params) {
  dt <- data.table::as.data.table(x)

  if (!nrow(dt)) {
    return(empty_normalized_rank_lut())
  }

  if (!"scientificName" %in% names(dt)) {
    assert("species" %in% names(dt),
           "LUT must contain either scientificName or species.")

    species_raw <- as.character(dt$species)

    if (all(species_raw %in% species_params$scientificName)) {
      dt[, scientificName := species_raw]
    } else if (all(species_raw %in% species_params$species)) {
      map <- species_params[, .(species, scientificName)]
      dt[, species := species_raw]
      dt <- merge(dt, map, by = "species", all.x = TRUE, sort = FALSE)
      assert(all(!is.na(dt$scientificName)),
             "Could not map all species IDs in LUT to scientificName.")
    } else {
      stop(
        "LUT species column is neither scientific names nor sid(scientificName) IDs.",
        call. = FALSE
      )
    }
  }

  assert("pu_id" %in% names(dt), "LUT is missing pu_id.")
  assert("patch_area_km2" %in% names(dt), "LUT is missing patch_area_km2.")

  if (!"patch_id" %in% names(dt)) {
    dt[, patch_id := seq_len(.N)]
  }

  dt[, scientificName := trimws(as.character(scientificName))]
  dt[, species := species_id(scientificName)]
  dt <- dt[species %in% species_params$species]

  if (!nrow(dt)) {
    return(empty_normalized_rank_lut())
  }

  dt[, patch_id := as.integer(patch_id)]
  dt[, pu_id := as.integer(pu_id)]
  dt[, patch_area_km2 := as.numeric(patch_area_km2)]

  bad <- dt[
    !is.finite(pu_id) |
      !is.finite(patch_area_km2) |
      patch_area_km2 <= 0
  ]

  assert(nrow(bad) == 0L,
         paste0("LUT for method=", method, ", stage=", stage,
                " contains invalid pu_id or patch_area_km2."))

  dt[, .(
    method = as.character(method),
    stage = as.integer(stage),
    scientificName,
    species,
    patch_id,
    pu_id,
    patch_area_km2
  )]
}

priority_paths <- function(
  curve = "q025",
  taxa = "mammals_birds",
  sdm = "ppm_rangebag",
  cells_to_remove_per_iteration = 1000L,
  pruning_iterations_per_stage = 100L,
  pruning_input_bundle_dir = file.path("Data", "Clean", "PriorityInputs"),
  pipeline_output_root = file.path("Data", "Results", "PriorityRuns"),
  ana_subdir = "ana"
) {
  run_id <- priority_run_id(
    curve = curve,
    taxa = taxa,
    sdm = sdm,
    cells_to_remove_per_iteration = cells_to_remove_per_iteration,
    pruning_iterations_per_stage = pruning_iterations_per_stage
  )
  out_curve <- file.path(pipeline_output_root, run_id)
  ana <- file.path(out_curve, ana_subdir)
  list(
    run_id = run_id,
    bundle_path = file.path(pruning_input_bundle_dir, priority_bundle_name(curve, taxa, sdm)),
    out_curve = out_curve,
    ana = ana
  )
}
