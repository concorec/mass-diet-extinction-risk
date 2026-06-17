source(file.path("R", "project_utils.R"))
source(file.path("R", "priority_run_config.R"))
source(file.path("R", "priority_outputs.R"))

run_priority_pipeline_from_saved_bundle <- function(
  curve_label = "q025",
  do_mammals = TRUE,
  do_birds = TRUE,
  do_ppm = TRUE,
  do_rangebag = TRUE,
  cells_to_remove_per_iteration = 1000L,
  pruning_iterations_per_stage = 100L,
  max_stages = 100L,
  pruning_input_bundle_dir = file.path("Data", "Clean", "PriorityInputs"),
  patch_raster_dir = file.path("Data", "Clean", "Patches"),
  pipeline_output_root = file.path("Data", "Results", "PriorityRuns")
) {
# ---------------------------------------------------------------------
# 1) User settings
# ---------------------------------------------------------------------
selection <- priority_selection_from_flags(
  curve_label = curve_label,
  do_mammals = do_mammals,
  do_birds = do_birds,
  do_ppm = do_ppm,
  do_rangebag = do_rangebag
)

CURVE_LABEL <- selection$curve_label
DO_MAMMALS <- selection$do_mammals
DO_BIRDS <- selection$do_birds
DO_PPM <- selection$do_ppm
DO_RANGEBAG <- selection$do_rangebag
taxon_selection_tag <- selection$taxa
sdm_selection_tag <- selection$sdm


# ---------------------------------------------------------------------
# 2) Validate user settings
# ---------------------------------------------------------------------
cells_to_remove_per_iteration <- validate_priority_count(
  cells_to_remove_per_iteration,
  "cells_to_remove_per_iteration"
)
pruning_iterations_per_stage <- validate_priority_count(
  pruning_iterations_per_stage,
  "pruning_iterations_per_stage"
)
max_stages <- validate_max_stages(max_stages)

# ---------------------------------------------------------------------
# 4) Paths
# ---------------------------------------------------------------------
run_paths <- priority_paths(
  curve = CURVE_LABEL,
  taxa = taxon_selection_tag,
  sdm = sdm_selection_tag,
  cells_to_remove_per_iteration = cells_to_remove_per_iteration,
  pruning_iterations_per_stage = pruning_iterations_per_stage,
  pruning_input_bundle_dir = pruning_input_bundle_dir,
  pipeline_output_root = pipeline_output_root
)

pruning_input_bundle_path <- run_paths$bundle_path
pipeline_output_dir <- run_paths$out_curve

# ---------------------------------------------------------------------
# 5) Package checks
# ---------------------------------------------------------------------
load_packages(c("data.table", "terra", "sf", "Rfast", "fastmatch"))


# ---------------------------------------------------------------------
# 6) Load the saved pruning-input bundle
# ---------------------------------------------------------------------
if (!file.exists(pruning_input_bundle_path)) {
  stop(
    "Saved pruning-input bundle not found: ", pruning_input_bundle_path, "\n",
    "Rebuild the input bundle with matching values of ",
    "CURVE_LABEL, DO_MAMMALS, DO_BIRDS, DO_PPM, and DO_RANGEBAG."
  )
}

pruning_input_bundle <- readRDS(pruning_input_bundle_path)

if (is.null(pruning_input_bundle$metadata)) {
  stop("Saved pruning-input bundle does not contain a metadata element.")
}


# ---------------------------------------------------------------------
# 7) Validate that bundle metadata matches the requested gates
# ---------------------------------------------------------------------
bundle_metadata <- pruning_input_bundle$metadata

required_metadata_fields <- c(
  "schema_version",
  "curve_label",
  "do_mammals",
  "do_birds",
  "do_ppm",
  "do_rangebag",
  "retained_species"
)

validate_priority_bundle_schema(pruning_input_bundle)

missing_metadata_fields <- setdiff(required_metadata_fields, names(bundle_metadata))

if (length(missing_metadata_fields)) {
  stop(
    "Saved pruning-input bundle is missing metadata field(s): ",
    paste(missing_metadata_fields, collapse = ", "), ".\n",
    "Regenerate the bundle with script 6."
  )
}

if (!identical(as.character(bundle_metadata$curve_label), as.character(CURVE_LABEL))) {
  stop(
    "Bundle CURVE_LABEL mismatch. Requested ", CURVE_LABEL,
    " but bundle contains ", bundle_metadata$curve_label, "."
  )
}

if (!identical(as.logical(bundle_metadata$do_mammals), DO_MAMMALS)) {
  stop("Bundle DO_MAMMALS mismatch. Regenerate or load the matching bundle.")
}

if (!identical(as.logical(bundle_metadata$do_birds), DO_BIRDS)) {
  stop("Bundle DO_BIRDS mismatch. Regenerate or load the matching bundle.")
}

if (!identical(as.logical(bundle_metadata$do_ppm), DO_PPM)) {
  stop("Bundle DO_PPM mismatch. Regenerate or load the matching bundle.")
}

if (!identical(as.logical(bundle_metadata$do_rangebag), DO_RANGEBAG)) {
  stop("Bundle DO_RANGEBAG mismatch. Regenerate or load the matching bundle.")
}

# ---------------------------------------------------------------------
# 8) Restore plain objects directly
# ---------------------------------------------------------------------
patch_table <- as.data.table(pruning_input_bundle$patch_table)
pu_graphs_by_key <- pruning_input_bundle$pu_graphs_by_key
alive_species_count_by_cell <- pruning_input_bundle$alive_species_count_by_cell
cell_area_by_cell <- pruning_input_bundle$cell_area_by_cell
species_params <- as.data.table(pruning_input_bundle$species_params)
rook_neighbor_pairs <- pruning_input_bundle$rook_neighbor_pairs


# ---------------------------------------------------------------------
# 9) Validate required restored objects
# ---------------------------------------------------------------------
if (!nrow(patch_table)) {
  stop("Loaded patch_table is empty.")
}

if (!is.list(pu_graphs_by_key)) {
  stop("Loaded pu_graphs_by_key is not a list.")
}

if (!length(alive_species_count_by_cell)) {
  stop("Loaded alive_species_count_by_cell is empty.")
}

if (!length(cell_area_by_cell)) {
  stop("Loaded cell_area_by_cell is empty.")
}

if (!nrow(species_params)) {
  stop("Loaded species_params is empty.")
}

if (!nrow(rook_neighbor_pairs)) {
  stop("Loaded rook_neighbor_pairs is empty.")
}

if (!("dispersal_distance_km" %in% names(species_params))) {
  stop(
    "The saved pruning-input bundle does not contain species_params$dispersal_distance_km. ",
    "Regenerate the bundle after patching the script-6 initialization code."
  )
}

if (any(!is.finite(species_params$dispersal_distance_km) |
        species_params$dispersal_distance_km <= 0)) {
  stop("species_params$dispersal_distance_km contains missing, non-finite, or non-positive values.")
}

required_species_param_columns <- c(
  "species",
  "density",
  "dispersal_distance_km",
  "min_patch_area_km2",
  "min_population_area_km2",
  "a_pred",
  "b_pred",
  "taxon_class",
  "sdm_method",
  "redlist_category"
)

missing_species_param_columns <- setdiff(required_species_param_columns, names(species_params))
if (length(missing_species_param_columns)) {
  stop(
    "Loaded species_params is missing required column(s): ",
    paste(missing_species_param_columns, collapse = ", ")
  )
}

validate_density_area_thresholds(
  density = species_params$density,
  min_patch_area_km2 = species_params$min_patch_area_km2,
  min_population_area_km2 = species_params$min_population_area_km2,
  label = "loaded species_params area thresholds"
)

required_patch_columns <- c("species", "patch_id", "pu_id", "patch_area_km2")
missing_patch_columns <- setdiff(required_patch_columns, names(patch_table))
if (length(missing_patch_columns)) {
  stop(
    "Loaded patch_table is missing required column(s): ",
    paste(missing_patch_columns, collapse = ", ")
  )
}


# ---------------------------------------------------------------------
# 10) Rebuild patch_id_by_species_env
# ---------------------------------------------------------------------
if (is.null(pruning_input_bundle$patch_id_by_species_list) ||
    !is.list(pruning_input_bundle$patch_id_by_species_list)) {
  stop("Bundle does not contain a valid patch_id_by_species_list.")
}

patch_id_by_species_env <- new.env(parent = emptyenv())

for (species_name in names(pruning_input_bundle$patch_id_by_species_list)) {
  assign(
    species_name,
    pruning_input_bundle$patch_id_by_species_list[[species_name]],
    envir = patch_id_by_species_env
  )
}


# ---------------------------------------------------------------------
# 11) Rebuild patch_cell_index_by_species_env
# ---------------------------------------------------------------------
if (is.null(pruning_input_bundle$patch_cell_index_by_species_list) ||
    !is.list(pruning_input_bundle$patch_cell_index_by_species_list)) {
  stop("Bundle does not contain a valid patch_cell_index_by_species_list.")
}

patch_cell_index_by_species_env <- new.env(parent = emptyenv())

for (species_name in names(pruning_input_bundle$patch_cell_index_by_species_list)) {
  assign(
    species_name,
    pruning_input_bundle$patch_cell_index_by_species_list[[species_name]],
    envir = patch_cell_index_by_species_env
  )
}


# ---------------------------------------------------------------------
# 12) Reconstruct raster templates
# ---------------------------------------------------------------------
retained_species <- bundle_metadata$retained_species

if (is.null(retained_species) || !length(retained_species)) {
  stop("The saved bundle metadata does not contain retained_species.")
}

retained_species <- as.character(retained_species)

template_species <- retained_species[[1]]

mask_template_raster_path <- file.path(
  patch_raster_dir,
  patch_filename_from_scientific(template_species)
)

if (!file.exists(mask_template_raster_path)) {
  stop(
    "Could not reconstruct mask template raster. Missing file: ",
    mask_template_raster_path
  )
}

mask_template_raster <- terra::rast(mask_template_raster_path)

species_template_layers <- lapply(retained_species, function(species_name) {
  raster_path <- file.path(
    patch_raster_dir,
    patch_filename_from_scientific(species_name)
  )

  if (!file.exists(raster_path)) {
    stop("Missing patch raster for retained species: ", raster_path)
  }

  r <- terra::rast(raster_path)
  names(r) <- species_name
  r
})

species_raster_template_stack <- terra::rast(species_template_layers)


# ---------------------------------------------------------------------
# 13) Optional post-load consistency checks
# ---------------------------------------------------------------------
species_in_env <- ls(envir = patch_id_by_species_env, all.names = TRUE)

if (!setequal(retained_species, species_in_env)) {
  stop(
    "Retained species in metadata do not match species stored in patch_id_by_species_list."
  )
}

if (!setequal(retained_species, names(species_raster_template_stack))) {
  stop(
    "Retained species in metadata do not match reconstructed species raster template stack."
  )
}

if (!all(unique(patch_table$species) %in% species_params$species)) {
  stop("Some species in patch_table are missing from species_params.")
}

if (!all(unique(patch_table$species) %in% retained_species)) {
  stop("Some species in patch_table are missing from retained_species metadata.")
}

if (length(alive_species_count_by_cell) != terra::ncell(mask_template_raster)) {
  stop(
    "alive_species_count_by_cell length does not match the number of cells in the mask template raster."
  )
}

if (length(cell_area_by_cell) != terra::ncell(mask_template_raster)) {
  stop(
    "cell_area_by_cell length does not match the number of cells in the mask template raster."
  )
}


# ---------------------------------------------------------------------
# 14) Create the root output directory
# ---------------------------------------------------------------------
dir.create(pipeline_output_dir, recursive = TRUE, showWarnings = FALSE)


# ---------------------------------------------------------------------
# 15) Optional convenience bundle in memory
# ---------------------------------------------------------------------
priority_pipeline_inputs <- list(
  curve_label                      = CURVE_LABEL,
  do_mammals                      = DO_MAMMALS,
  do_birds                        = DO_BIRDS,
  do_ppm                          = DO_PPM,
  do_rangebag                     = DO_RANGEBAG,
  cells_to_remove_per_iteration   = cells_to_remove_per_iteration,
  pruning_iterations_per_stage    = pruning_iterations_per_stage,
  max_stages                      = max_stages,
  patch_table                     = patch_table,
  pu_graphs_by_key                = pu_graphs_by_key,
  alive_species_count_by_cell     = alive_species_count_by_cell,
  patch_id_by_species_env         = patch_id_by_species_env,
  patch_cell_index_by_species_env = patch_cell_index_by_species_env,
  cell_area_by_cell               = cell_area_by_cell,
  rook_neighbor_pairs             = rook_neighbor_pairs,
  species_params                  = species_params,
  species_raster_template_stack   = species_raster_template_stack,
  mask_template_raster            = mask_template_raster,
  output_dir                      = pipeline_output_dir
)


# ---------------------------------------------------------------------
# 16) Run the full priority pipeline
# ---------------------------------------------------------------------
priority_pipeline_result <- run_priority_pipeline(
  cells_to_remove_per_iteration   = cells_to_remove_per_iteration,
  pruning_iterations_per_stage    = pruning_iterations_per_stage,
  patch_table                     = patch_table,
  pu_graphs_by_key                = pu_graphs_by_key,
  alive_species_count_by_cell     = alive_species_count_by_cell,
  patch_id_by_species_env         = patch_id_by_species_env,
  patch_cell_index_by_species_env = patch_cell_index_by_species_env,
  cell_area_by_cell               = cell_area_by_cell,
  rook_neighbor_pairs             = rook_neighbor_pairs,
  species_params                  = species_params,
  species_raster_template_stack   = species_raster_template_stack,
  mask_template_raster            = mask_template_raster,
  output_dir                      = pipeline_output_dir,
  max_stages                      = max_stages
)


# ---------------------------------------------------------------------
# 17) Small summary prints
# ---------------------------------------------------------------------
message("Priority pipeline completed.")
message("  Loaded input bundle:      ", pruning_input_bundle_path)
message("  Pipeline output dir:      ", pipeline_output_dir)
message(
  "  Removal order raster:     ",
  normalizePath(priority_pipeline_result$removal_order_path, mustWork = FALSE)
)
message(
  "  Zonation-style rank map:  ",
  normalizePath(priority_pipeline_result$rankmap_path, mustWork = FALSE)
)
message(
  "  Removal events table:     ",
  normalizePath(priority_pipeline_result$removal_events_path, mustWork = FALSE)
)
message(
  "  Patch lookup directory:   ",
  normalizePath(priority_pipeline_result$patch_lookup_output_dir, mustWork = FALSE)
)
message("  Curve label:              ", CURVE_LABEL)
message("  Taxa selected:            ", taxon_selection_tag)
message("  SDM sources selected:     ", sdm_selection_tag)
message("  Completed stages:         ", priority_pipeline_result$completed_stages)
message("  Frontier exhausted:       ", priority_pipeline_result$frontier_exhausted)
message("  Removal steps recorded:  ", priority_pipeline_result$removal_step_count)
message("  Initially alive cells:    ", priority_pipeline_result$initial_alive_cell_count)
message("  Removed initial cells:    ", priority_pipeline_result$removed_initial_alive_cell_count)
message("  Unremoved initial cells:  ", priority_pipeline_result$unremoved_initial_alive_cell_count)
message("  Final remaining patches:  ", nrow(priority_pipeline_result$patch_table))
message(
  "  Final remaining PUs:      ",
  data.table::uniqueN(
    priority_pipeline_result$patch_table,
    by = c("species", "pu_id")
  )
)
message(
  "  Final alive cells:        ",
  sum(priority_pipeline_result$alive_species_count_by_cell > 0L)
)

invisible(list(
  pruning_input_bundle_path = pruning_input_bundle_path,
  pipeline_output_dir = pipeline_output_dir,
  priority_pipeline_inputs = priority_pipeline_inputs,
  priority_pipeline_result = priority_pipeline_result
))
}
