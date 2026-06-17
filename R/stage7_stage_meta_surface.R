
suppressPackageStartupMessages({
  library(terra)
  library(data.table)
})

source(file.path("R", "project_utils.R"))
source(file.path("R", "priority_run_config.R"))

# ============================================================
# Parameters
# ============================================================

curve <- as.character(params$curve)
taxa <- as.character(params$taxa)
sdm <- as.character(params$sdm)

curve <- validate_persistence_curve(curve, "params$curve")

priority_flags <- priority_flags_from_tags(taxa, sdm)

cells_to_remove_per_iteration <- as.integer(params$cells_to_remove_per_iteration)
pruning_iterations_per_stage <- as.integer(params$pruning_iterations_per_stage)

assert(is.finite(cells_to_remove_per_iteration) &&
         cells_to_remove_per_iteration > 0L,
       "params$cells_to_remove_per_iteration must be a positive integer.")

assert(is.finite(pruning_iterations_per_stage) &&
         pruning_iterations_per_stage > 0L,
       "params$pruning_iterations_per_stage must be a positive integer.")

run_paths <- priority_paths(
  curve = curve,
  taxa = taxa,
  sdm = sdm,
  cells_to_remove_per_iteration = cells_to_remove_per_iteration,
  pruning_iterations_per_stage = pruning_iterations_per_stage,
  ana_subdir = as.character(params$ana_subdir)
)

run_id <- run_paths$run_id
bundle_path <- run_paths$bundle_path
out_curve <- run_paths$out_curve
removal_order_path <- file.path(out_curve, "removal_order.tif")
rankmap_path <- file.path(out_curve, "rankmap.tif")
removal_events_path <- file.path(out_curve, "removal_events.csv")
patch_lookup_dir <- file.path(out_curve, "patch_lookup_tables")

overwrite  <- isTRUE(params$overwrite)

need_file(bundle_path, "pruning input bundle")
need_dir(out_curve, "optimized pipeline output directory")
need_file(removal_order_path, "removal_order.tif")
need_file(rankmap_path, "rankmap.tif")
need_file(removal_events_path, "removal_events.csv")
need_dir(patch_lookup_dir, "stage patch lookup directory")

ana <- run_paths$ana
dir.create(ana, recursive = TRUE, showWarnings = FALSE)

stage_meta_csv <- file.path(ana, as.character(params$stage_meta_name))
priority_surface_tif <- rankmap_path

cat("\n=== 7.1_stage_meta_and_priority_surface ===\n")
cat("run_id      :", run_id, "\n")
cat("bundle_path :", bundle_path, "\n")
cat("out_curve   :", out_curve, "\n\n")

# ============================================================
# Load pruning-input bundle
# ============================================================

bundle <- readRDS(bundle_path)
validate_priority_bundle_schema(bundle)

assert(!is.null(bundle$metadata),
       "Bundle is missing metadata.")

assert(!is.null(bundle$alive_species_count_by_cell),
       "Bundle is missing alive_species_count_by_cell.")

assert(!is.null(bundle$cell_area_by_cell),
       "Bundle is missing cell_area_by_cell.")

assert(!is.null(bundle$metadata$retained_species),
       "Bundle is missing metadata$retained_species.")

# Optional metadata validation. This catches mismatches such as using a
# mammals+birds bundle with a mammals-only priority run.
if (isTRUE(params$validate_bundle_metadata)) {
  required_metadata_fields <- c(
    "curve_label",
    "do_mammals",
    "do_birds",
    "do_ppm",
    "do_rangebag",
    "retained_species"
  )

  missing_metadata_fields <- setdiff(required_metadata_fields, names(bundle$metadata))
  assert(!length(missing_metadata_fields),
         paste("Bundle metadata is missing field(s):",
               paste(missing_metadata_fields, collapse = ", ")))

  assert(identical(as.character(bundle$metadata$curve_label), curve),
         paste0("Bundle curve_label mismatch. Expected ", curve,
                " but bundle contains ", bundle$metadata$curve_label, "."))

  assert(identical(as.logical(bundle$metadata$do_mammals), priority_flags$do_mammals),
         "Bundle do_mammals does not match params$taxa.")

  assert(identical(as.logical(bundle$metadata$do_birds), priority_flags$do_birds),
         "Bundle do_birds does not match params$taxa.")

  assert(identical(as.logical(bundle$metadata$do_ppm), priority_flags$do_ppm),
         "Bundle do_ppm does not match params$sdm.")

  assert(identical(as.logical(bundle$metadata$do_rangebag), priority_flags$do_rangebag),
         "Bundle do_rangebag does not match params$sdm.")

}

# ============================================================
# Initial alive-cell domain from bundle
# ============================================================

alive0 <- as.integer(bundle$alive_species_count_by_cell) > 0L
cell_area_by_cell <- as.numeric(bundle$cell_area_by_cell)

template <- terra::rast(removal_order_path)

assert(length(alive0) == terra::ncell(template),
       "alive_species_count_by_cell length does not match removal_order.tif ncell.")

assert(length(cell_area_by_cell) == terra::ncell(template),
       "cell_area_by_cell length does not match removal_order.tif ncell.")

alive0_n <- sum(alive0)
alive0_area_km2 <- sum(cell_area_by_cell[alive0], na.rm = TRUE)

assert(alive0_n > 0,
       "No initially alive cells found in bundle.")

assert(is.finite(alive0_area_km2) && alive0_area_km2 > 0,
       "Initial alive area is zero or non-finite.")

assert(all(is.finite(cell_area_by_cell[alive0]) & cell_area_by_cell[alive0] > 0),
       "Some initially alive cells have missing, non-finite, or non-positive area.")

# ============================================================
# Load removal_events.csv. The final_retained row is a rank-assignment event,
# not a true ecological removal stage.
# ============================================================

rank_events <- fread(removal_events_path)

required_event_cols <- c(
  "removal_step",
  "stage",
  "event_type",
  "pruning_iteration",
  "cells_removed",
  "area_removed_km2",
  "cum_cells_removed",
  "cum_area_removed_km2",
  "cum_prop_cells_removed",
  "cum_prop_area_removed",
  "cells_retained",
  "area_retained_km2"
)

missing_event_cols <- setdiff(required_event_cols, names(rank_events))
assert(!length(missing_event_cols),
       paste("removal_events.csv is missing column(s):",
             paste(missing_event_cols, collapse = ", ")))

assert(sum(rank_events$event_type == "final_retained", na.rm = TRUE) == 1L,
       "removal_events.csv must contain exactly one final_retained event. Regenerate script 6 outputs.")

if (nrow(rank_events)) {
  setorder(rank_events, removal_step)

  assert(all(is.finite(rank_events$removal_step)),
         "removal_events.csv contains non-finite removal_step values.")

  assert(all(is.finite(rank_events$stage)),
         "removal_events.csv contains non-finite stage values.")

  assert(anyDuplicated(rank_events$removal_step) == 0L,
         "removal_events.csv contains duplicated removal_step values.")

  assert(all(rank_events$cells_removed >= 0),
         "removal_events.csv contains negative cells_removed values.")

  assert(all(rank_events$area_removed_km2 >= 0),
         "removal_events.csv contains negative area_removed_km2 values.")

  assert(all(rank_events$cum_prop_cells_removed >= 0 &
               rank_events$cum_prop_cells_removed <= 1),
         "cum_prop_cells_removed must be in [0, 1].")

  assert(all(rank_events$cum_prop_area_removed >= 0 &
               rank_events$cum_prop_area_removed <= 1),
         "cum_prop_area_removed must be in [0, 1].")
}

removal_events <- rank_events[event_type != "final_retained"]

# ============================================================
# Parse stage patch lookup tables
# ============================================================

stage_files <- list.files(
  patch_lookup_dir,
  pattern = "^stage_patch_lookup_stage_[0-9]{4}\\.csv$",
  full.names = TRUE
)

assert(length(stage_files) > 0,
       paste0("No stage_patch_lookup_stage_####.csv files found in: ", patch_lookup_dir))

st <- data.table(path = stage_files, base = basename(stage_files))
st[, stage := as.integer(sub(
  "^stage_patch_lookup_stage_([0-9]{4})\\.csv$",
  "\\1",
  base
))]
st <- st[is.finite(stage)]
setorder(st, stage)

assert(nrow(st) > 0,
       "Could not parse any stage numbers from patch lookup tables.")

assert(anyDuplicated(st$stage) == 0L,
       "Duplicate stage patch lookup files detected.")

assert(identical(st$stage, seq_len(max(st$stage))),
       "Stage patch lookup tables must be contiguous from stage 1.")

if (nrow(removal_events)) {
  assert(max(removal_events$stage, na.rm = TRUE) <= max(st$stage),
         "removal_events.csv contains events for stages without patch lookup tables.")
}

# ============================================================
# Run summary
# ============================================================

cat("\n=== 7.1 stage_meta_and_priority_surface ===\n")
cat("curve                    :", curve, "\n")
cat("bundle_path              :", bundle_path, "\n")
cat("pipeline_output_dir      :", out_curve, "\n")
cat("ana                      :", ana, "\n")
cat("alive0_n                 :", alive0_n, "\n")
cat("alive0_area_km2          :", alive0_area_km2, "\n")
cat("n_stage_lookup_tables    :", nrow(st), "\n")
cat("priority_surface_units   : cells (Zonation-style rank order)\n")
cat("priority_surface_path    :", priority_surface_tif, "\n")
cat("overwrite                :", overwrite, "\n\n")

# ============================================================
# 01_stage_meta - stage_meta.csv from removal_events.csv
#                 and stage patch lookup files
#
# Includes explicit stage 0 baseline row.
# ============================================================

if (file.exists(stage_meta_csv) && !overwrite) {
  cat("[stage_meta] Exists -> loading: ", stage_meta_csv, "\n", sep = "")
  stage_meta <- fread(stage_meta_csv)

  assert("stage" %in% names(stage_meta),
         "Existing stage_meta.csv is missing stage column.")

  assert(stage_meta[stage == 0L, .N] == 1L,
         "Existing stage_meta.csv must contain exactly one stage 0 row.")

} else {
  cat("[stage_meta] Building from removal_events.csv + stage patch lookup files ...\n")

  stage0 <- data.table(
    stage = 0L,
    patch_lookup_path = NA_character_,

    last_removal_step = 0L,

    alive_end = as.numeric(alive0_n),
    removed_end = 0,

    area_retained_km2 = as.numeric(alive0_area_km2),
    area_removed_km2 = 0,

    cells_removed_stage = 0,
    area_removed_stage_km2 = 0,
    events_in_stage = 0,

    prop_cells_removed_end = 0,
    prop_area_removed_end = 0,

    pct_cells_removed_end = 0,
    pct_area_removed_end = 0
  )

  stage_rows <- rbindlist(lapply(st$stage, function(s) {
    ev_to_stage <- removal_events[stage <= s]
    ev_this_stage <- removal_events[stage == s]

    this_lookup <- st[stage == s, path][1]

    if (nrow(ev_to_stage)) {
      last_ev <- ev_to_stage[.N]

      last_removal_step <- as.integer(last_ev$removal_step)
      alive_end <- as.numeric(last_ev$cells_retained)
      removed_end <- as.numeric(last_ev$cum_cells_removed)

      area_retained_km2 <- as.numeric(last_ev$area_retained_km2)
      area_removed_km2 <- as.numeric(last_ev$cum_area_removed_km2)

      prop_cells_removed_end <- as.numeric(last_ev$cum_prop_cells_removed)
      prop_area_removed_end <- as.numeric(last_ev$cum_prop_area_removed)
    } else {
      last_removal_step <- 0L
      alive_end <- as.numeric(alive0_n)
      removed_end <- 0

      area_retained_km2 <- as.numeric(alive0_area_km2)
      area_removed_km2 <- 0

      prop_cells_removed_end <- 0
      prop_area_removed_end <- 0
    }

    data.table(
      stage = as.integer(s),
      patch_lookup_path = this_lookup,

      last_removal_step = last_removal_step,

      alive_end = alive_end,
      removed_end = removed_end,

      area_retained_km2 = area_retained_km2,
      area_removed_km2 = area_removed_km2,

      cells_removed_stage = sum(ev_this_stage$cells_removed, na.rm = TRUE),
      area_removed_stage_km2 = sum(ev_this_stage$area_removed_km2, na.rm = TRUE),
      events_in_stage = nrow(ev_this_stage),

      prop_cells_removed_end = prop_cells_removed_end,
      prop_area_removed_end = prop_area_removed_end,

      pct_cells_removed_end = 100 * prop_cells_removed_end,
      pct_area_removed_end = 100 * prop_area_removed_end
    )
  }), use.names = TRUE, fill = TRUE)

  stage_meta <- rbindlist(
    list(stage0, stage_rows),
    use.names = TRUE,
    fill = TRUE
  )

  setorder(stage_meta, stage)

  assert(stage_meta[stage == 0L, .N] == 1L,
         "stage_meta must contain exactly one stage 0 row.")

  assert(identical(stage_meta$stage, seq.int(0L, max(stage_meta$stage))),
         "stage_meta stages must be contiguous from stage 0.")

  assert(all(stage_meta$alive_end >= 0),
         "stage_meta contains negative alive_end values.")

  assert(all(stage_meta$removed_end >= 0),
         "stage_meta contains negative removed_end values.")

  assert(all(stage_meta$area_retained_km2 >= 0),
         "stage_meta contains negative retained areas.")

  assert(all(stage_meta$area_removed_km2 >= 0),
         "stage_meta contains negative removed areas.")

  assert(all(stage_meta$prop_cells_removed_end >= 0 &
               stage_meta$prop_cells_removed_end <= 1),
         "stage_meta prop_cells_removed_end must be in [0, 1].")

  assert(all(stage_meta$prop_area_removed_end >= 0 &
               stage_meta$prop_area_removed_end <= 1),
         "stage_meta prop_area_removed_end must be in [0, 1].")

  fwrite(stage_meta, stage_meta_csv)
  cat("[stage_meta] Wrote: ", stage_meta_csv, " rows=", nrow(stage_meta), "\n", sep = "")
}

# ============================================================
# 02_priority_surface - validate canonical rankmap.tif
#
# rankmap.tif:
#   Per-cell cumulative proportion of analysis-domain cells removed when the
#   cell was lost. This follows the Zonation rankmap convention most closely:
#   values are ordinal ranks spread over the ranked raster-cell sequence.
#   Cells retained to the terminal layer get 1. Cells outside the initial
#   alive-cell domain remain NA.
# ============================================================

cat("[priority_surface] Validating canonical rankmap.tif ...\n")

removal_order <- terra::rast(removal_order_path)
priority_surface <- terra::rast(priority_surface_tif)

assert(
  terra::compareGeom(removal_order, priority_surface, stopOnError = FALSE),
  "rankmap.tif geometry does not match removal_order.tif."
)

n_removed_raster <- as.numeric(terra::global(
  terra::ifel(!is.na(removal_order), 1, 0),
  "sum",
  na.rm = TRUE
)[1, 1])

n_ranked_events <- if (nrow(rank_events)) {
  sum(rank_events$cells_removed, na.rm = TRUE)
} else {
  0
}

assert(
  as.integer(n_removed_raster) == as.integer(n_ranked_events),
  paste0(
    "Mismatch between removal_order.tif and removal_events.csv: raster has ",
    n_removed_raster,
    " ranked cells, events table has ",
    n_ranked_events,
    "."
  )
)

alive0_r <- terra::rast(template)
terra::values(alive0_r) <- as.integer(alive0)

missing_ranked_alive_count <- as.numeric(terra::global(
  terra::ifel(alive0_r == 1L & is.na(priority_surface), 1, 0),
  "sum",
  na.rm = TRUE
)[1, 1])

assert(missing_ranked_alive_count == 0,
       "Some initially alive cells have no rank value. Regenerate script 6 outputs.")

outside_alive_count <- as.numeric(terra::global(
  terra::ifel(alive0_r != 1L & !is.na(priority_surface), 1, 0),
  "sum",
  na.rm = TRUE
)[1, 1])

assert(outside_alive_count == 0,
       "rankmap.tif has non-NA values outside the initial alive-cell domain.")

rng <- terra::global(priority_surface, c("min", "max"), na.rm = TRUE)
min_val <- as.numeric(rng[1, "min"])
max_val <- as.numeric(rng[1, "max"])

assert(is.na(min_val) || min_val >= 0,
       "rankmap.tif contains values below 0.")

assert(is.na(max_val) || max_val <= 1,
       "rankmap.tif contains values above 1.")

cat("[priority_surface] Using existing rankmap.tif:\n")
cat("  - ", priority_surface_tif, "\n", sep = "")

rm(removal_order, alive0_r, priority_surface)
gc(FALSE)
