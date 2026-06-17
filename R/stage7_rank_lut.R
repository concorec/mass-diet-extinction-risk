
suppressPackageStartupMessages({
  library(terra)
  library(sf)
  library(igraph)
  library(data.table)
  library(units)
})

source(file.path("R", "project_utils.R"))
source(file.path("R", "priority_run_config.R"))

sid <- species_id

# ============================================================
# Params and paths
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

rank_method <- validate_benchmark_rank_method(
  params$rank_method,
  "params$rank_method"
)

overwrite <- isTRUE(params$overwrite)

run_paths <- priority_paths(
  curve = curve,
  taxa = taxa,
  sdm = sdm,
  cells_to_remove_per_iteration = cells_to_remove_per_iteration,
  pruning_iterations_per_stage = pruning_iterations_per_stage,
  ana_subdir = as.character(params$ana_subdir)
)

run_id <- run_paths$run_id
out_curve <- run_paths$out_curve
bundle_path <- run_paths$bundle_path
ana <- run_paths$ana
stage_meta_csv <- file.path(ana, as.character(params$stage_meta_name))

benchmark_paths <- benchmark_artifact_paths(
  run_paths,
  rank_method,
  rank_map_subdir = as.character(params$rank_map_subdir)
)
rank_path <- benchmark_paths$rank_path

species_csv <- as.character(params$species_csv)
patch_dir <- as.character(params$patch_dir)
grass_dir <- resolve_optional_grass_dir(params$grass_dir)

lut_root <- benchmark_paths$lut_root
dir.create(lut_root, recursive = TRUE, showWarnings = FALSE)

target_summary_csv <- benchmark_paths$target_summary_csv

lut_prefix <- as.character(params$lut_prefix)

need_dir(out_curve, "pipeline_output_dir")
need_file(bundle_path, "pruning input bundle")
need_file(stage_meta_csv, "stage_meta.csv from 7.1")
need_file(rank_path, "run-specific benchmark rank map")
need_file(species_csv, "species_table.csv")
need_dir(patch_dir, "patch raster directory")
if (!is.na(optional_path(grass_dir))) {
  need_dir(grass_dir, "GRASS GIS installation directory")
}
need_dir(lut_root, "rank_lut output directory")
init_faster_raster(grass_dir)

cat("\n=== 7.2_build_rank_lut ===\n")
cat("run_id                :", run_id, "\n")
cat("rank_method           :", rank_method, "\n")
cat("pipeline_output_dir   :", out_curve, "\n")
cat("ana                   :", ana, "\n")
cat("rank_path             :", rank_path, "\n")
cat("lut_root              :", lut_root, "\n")
cat("grass_dir             :", grass_dir, "\n")
cat("match_basis           : cells (Zonation-style rank order)\n")
cat("rank_values           : higher values retained later\n\n")

# ============================================================
# 01_load_inputs - stage_meta, bundle, species table, rank map
# ============================================================

bundle <- readRDS(bundle_path)
validate_priority_bundle_schema(bundle)

required_bundle_fields <- c(
  "metadata",
  "alive_species_count_by_cell",
  "cell_area_by_cell",
  "patch_table"
)

missing_bundle_fields <- setdiff(required_bundle_fields, names(bundle))
assert(!length(missing_bundle_fields),
       paste("Bundle is missing field(s):", paste(missing_bundle_fields, collapse = ", ")))

assert(!is.null(bundle$metadata$retained_species),
       "Bundle metadata is missing retained_species.")

retained_species <- as.character(bundle$metadata$retained_species)

alive0 <- as.integer(bundle$alive_species_count_by_cell) > 0L
cell_area_by_cell <- as.numeric(bundle$cell_area_by_cell)

stage_meta <- fread(stage_meta_csv)

required_stage_cols <- c(
  "stage",
  "alive_end",
  "area_retained_km2"
)

missing_stage_cols <- setdiff(required_stage_cols, names(stage_meta))
assert(!length(missing_stage_cols),
       paste("stage_meta.csv is missing column(s):", paste(missing_stage_cols, collapse = ", ")))

assert(stage_meta[stage == 0L, .N] == 1L,
       "stage_meta.csv must contain exactly one stage 0 row.")

setorder(stage_meta, stage)

species_table <- fread(species_csv)

required_species_cols <- c(
  "scientificName",
  "density",
  "min_patch_size",
  "min_pop_size",
  "dispersal_dist"
)

missing_species_cols <- setdiff(required_species_cols, names(species_table))
assert(!length(missing_species_cols),
       paste("species_table.csv is missing column(s):",
             paste(missing_species_cols, collapse = ", ")))

sp <- species_table[scientificName %in% retained_species]
assert(nrow(sp) == length(retained_species),
       "species_table does not contain exactly the retained species in the bundle.")

sp[, species := sid(scientificName)]
sp <- as.data.table(add_canonical_area_thresholds(
  sp,
  validate = FALSE,
  label = "stage-7 species-table area thresholds"
))
validate_density_area_thresholds(
  density = suppressWarnings(as.numeric(sp$density)),
  min_patch_area_km2 = sp$min_patch_area_km2,
  min_population_area_km2 = sp$min_population_area_km2,
  label = "stage-7 species-table area thresholds"
)
sp[, patch_path := patch_file_from_name(scientificName, patch_dir)]

missing_patch <- sp[!file.exists(patch_path)]
assert(nrow(missing_patch) == 0L,
       paste("Missing patch raster(s):",
             paste(missing_patch$patch_path, collapse = "\n")))

# Use removal_order.tif as the template because it is guaranteed to match
# the optimized priority run.
template_path <- file.path(out_curve, "removal_order.tif")
need_file(template_path, "removal_order.tif template")
template <- terra::rast(template_path)

assert(length(alive0) == terra::ncell(template),
       "Bundle alive_species_count_by_cell length does not match template ncell.")

assert(length(cell_area_by_cell) == terra::ncell(template),
       "Bundle cell_area_by_cell length does not match template ncell.")

assert(all(is.finite(cell_area_by_cell[alive0]) & cell_area_by_cell[alive0] > 0),
       "Some initially alive cells have missing/non-positive cell areas.")

cell_area_raster <- terra::rast(template)
terra::values(cell_area_raster) <- cell_area_by_cell
names(cell_area_raster) <- "cell_area_km2"

# Load and align rank map.
rank <- terra::rast(rank_path)

same_geom <- terra::compareGeom(rank, template, stopOnError = FALSE)
assert(same_geom,
       paste0(
         "Benchmark rank map geometry must match removal_order.tif exactly. ",
         "Reproject/resample the rank map before running 7.2."
       ))

assert(terra::ncell(rank) == terra::ncell(template),
       "Rank map ncell does not match template.")

cat("retained species:", nrow(sp), "\n")
cat("stage_meta rows :", nrow(stage_meta), "\n")
cat("alive0 cells    :", sum(alive0), "\n")
cat("alive0 area km2 :", sum(cell_area_by_cell[alive0]), "\n\n")

# ============================================================
# 02_rank_stage_targets - match benchmark rank map to stage_meta
# ============================================================

rank_values <- terra::values(rank, mat = FALSE)

assert(length(rank_values) == terra::ncell(template),
       "rank value vector length does not match template ncell.")

rank_dt <- data.table(
  cell = seq_along(rank_values),
  rank_value = as.numeric(rank_values),
  cell_area_km2 = cell_area_by_cell,
  alive0 = alive0
)

rank_dt <- rank_dt[alive0 == TRUE]

# Strong default: fail if any initially alive cell has a missing benchmark rank.
# That prevents silent mismatch between the benchmark and the persistence run.
missing_rank_n <- rank_dt[!is.finite(rank_value), .N]
assert(missing_rank_n == 0L,
       paste0("Benchmark rank map has missing/non-finite values for ",
              missing_rank_n, " initially alive cells."))

setorder(rank_dt, -rank_value, cell)

rank_dt[, cum_cells := seq_len(.N)]
rank_dt[, cum_area_km2 := cumsum(cell_area_km2)]

total_ranked_cells <- nrow(rank_dt)
total_ranked_area <- rank_dt[.N, cum_area_km2]

assert(total_ranked_cells == sum(alive0),
       "Ranked alive-cell count does not equal bundle alive0 count.")

all_stages <- sort(unique(stage_meta$stage))
assert(length(all_stages) > 0L,
       "stage_meta.csv contains no stages for rank-map LUT reconstruction.")

stage_targets <- copy(stage_meta)

find_keep_n_by_cells <- function(target_cells, n_available) {
  if (!is.finite(target_cells) || target_cells <= 0) return(0L)
  as.integer(max(0L, min(n_available, round(target_cells))))
}

stage_targets[, keep_n := vapply(
  alive_end,
  find_keep_n_by_cells,
  integer(1),
  n_available = total_ranked_cells
)]

stage_targets[, actual_cells_retained := keep_n]
stage_targets[, actual_area_retained_km2 := vapply(
  keep_n,
  function(k) {
    if (is.finite(k) && k > 0L) rank_dt$cum_area_km2[k] else 0
  },
  numeric(1)
)]

stage_targets[, target_cells_retained := alive_end]
stage_targets[, target_area_retained_km2 := area_retained_km2]

stage_targets[, area_match_error_km2 :=
                actual_area_retained_km2 - target_area_retained_km2]

stage_targets[, abs_area_match_error_km2 := abs(area_match_error_km2)]

stage_targets[, rank_cutoff := vapply(
  keep_n,
  function(k) {
    if (is.finite(k) && k > 0L) rank_dt$rank_value[k] else NA_real_
  },
  numeric(1)
)]

s0 <- stage_targets[stage == 0L]
assert(nrow(s0) == 1L && s0$keep_n == total_ranked_cells,
       "Stage 0 target must retain all initially alive cells. Check stage_meta or rank map.")

fwrite(stage_targets, target_summary_csv)
cat("[targets] Wrote: ", target_summary_csv, "\n", sep = "")

# ============================================================
# 03_reconstruction_helpers - clumping, area, connectivity, LUTs
# ============================================================

make_binary_raster_from_cells <- function(template, cells) {
  r <- terra::rast(template)
  terra::values(r) <- NA_integer_

  if (length(cells)) {
    r[cells] <- 1L
  }

  r
}

clump_rook <- function(binary_raster, trim = TRUE) {
  r01 <- terra::ifel(!is.na(binary_raster) & binary_raster != 0, 1L, NA_integer_)

  if (isTRUE(trim)) {
    r01 <- tryCatch(
      terra::trim(r01),
      error = function(e) NULL
    )
  }

  if (is.null(r01) || terra::ncell(r01) == 0L) {
    return(NULL)
  }

  has_cells <- tryCatch({
    as.numeric(terra::global(
      terra::ifel(!is.na(r01), 1, 0),
      "sum",
      na.rm = TRUE
    )[1, 1])
  }, error = function(e) 0)

  if (!is.finite(has_cells) || has_cells <= 0) {
    return(NULL)
  }

  cl <- clump_binary_raster(r01, diagonal = FALSE)

  # Defensive normalization: patch IDs must be positive integers, outside = NA.
  cl <- terra::ifel(!is.na(cl) & cl > 0, cl, NA_integer_)
  names(cl) <- "patch_id"

  cl
}


patch_areas_from_clump <- function(clump_raster, cell_area_raster) {
  if (is.null(clump_raster) || terra::ncell(clump_raster) == 0L) {
    return(data.table(
      patch_id = integer(),
      patch_area_km2 = numeric(),
      patch_n_cells = integer()
    ))
  }

  # Crop/resample cell-area raster to the clump raster's extent/geometry.
  ca <- terra::crop(cell_area_raster, clump_raster, snap = "out")

  if (!terra::compareGeom(ca, clump_raster, stopOnError = FALSE)) {
    ca <- terra::resample(ca, clump_raster, method = "near")
  }

  area_df <- as.data.frame(
    terra::zonal(ca, clump_raster, "sum", na.rm = TRUE)
  )

  if (!nrow(area_df) || ncol(area_df) < 2L) {
    return(data.table(
      patch_id = integer(),
      patch_area_km2 = numeric(),
      patch_n_cells = integer()
    ))
  }

  area_dt <- as.data.table(area_df)
  setnames(area_dt, names(area_dt)[1:2], c("patch_id", "patch_area_km2"))

  area_dt <- area_dt[
    is.finite(patch_id) &
      is.finite(patch_area_km2) &
      patch_area_km2 > 0
  ]

  if (!nrow(area_dt)) {
    return(data.table(
      patch_id = integer(),
      patch_area_km2 = numeric(),
      patch_n_cells = integer()
    ))
  }
  
  freq_dt <- as.data.table(terra::freq(clump_raster, bylayer = FALSE))
  
  if (!nrow(freq_dt)) {
    area_dt[, patch_n_cells := NA_integer_]
    area_dt[, patch_id := as.integer(patch_id)]
    return(area_dt[])
  }
  
  assert(all(c("value", "count") %in% names(freq_dt)),
         "Unexpected terra::freq() output schema.")

  n_dt <- freq_dt[
    is.finite(value),
    .(
      patch_id = as.integer(value),
      patch_n_cells = as.integer(count)
    )
  ]

  area_dt[, patch_id := as.integer(patch_id)]

  out <- merge(area_dt, n_dt, by = "patch_id", all.x = TRUE)
  out[is.na(patch_n_cells), patch_n_cells := NA_integer_]

  out[]
}

polygonize_kept_patches <- function(clump_raster, keep_patch_ids) {
  if (!length(keep_patch_ids)) return(NULL)

  r <- terra::ifel(clump_raster %in% keep_patch_ids, clump_raster, NA)

  pol <- terra::as.polygons(r, dissolve = TRUE, na.rm = TRUE)
  if (is.null(pol) || terra::nrow(pol) == 0) return(NULL)

  names(pol)[1] <- "patch_id"

  sf::st_as_sf(pol)
}

assign_population_units <- function(patch_sf, patch_area_dt, disp_km, min_pop_km2) {
  min_pop_km2 <- validate_area_threshold(min_pop_km2, "min_pop_km2")

  if (is.null(patch_sf) || !nrow(patch_sf)) {
    return(data.table(patch_id = integer(), pu_id = integer()))
  }

  patch_sf$patch_id <- as.integer(patch_sf$patch_id)

  # sf distance is in meters for lon/lat with s2; in CRS units otherwise.
  # For implementation, verify CRS behavior once with your rasters.
  patch_sf <- sf::st_make_valid(patch_sf)
  
  assert(is.finite(disp_km) && disp_km >= 0,
       "disp_km must be finite and non-negative.")
  
  if (disp_km <= 0) {
    pu_dt <- patch_area_dt[
      patch_id %in% patch_sf$patch_id,
      .(patch_id, patch_area_km2)
    ]
  
    pu_dt <- pu_dt[area_exceeds_threshold(patch_area_km2, min_pop_km2)]
  
    if (!nrow(pu_dt)) {
      return(data.table(patch_id = integer(), pu_id = integer()))
    }
  
    pu_dt[, pu_id := seq_len(.N)]
    return(pu_dt[, .(patch_id, pu_id)])
  }
  
  is_longlat <- sf::st_is_longlat(patch_sf)

  assert(!is.na(is_longlat),
       "Patch polygons have missing CRS; distance-based PU linking cannot be evaluated safely.")
  
  if (is_longlat) {
    old_s2 <- sf::sf_use_s2(TRUE)
    on.exit(sf::sf_use_s2(old_s2), add = TRUE)
    dist_arg <- units::set_units(disp_km, "km")
  } else {
    dist_arg <- units::set_units(disp_km * 1000, "m")
  }
  
  neighbors <- sf::st_is_within_distance(
    patch_sf,
    patch_sf,
    dist = dist_arg
  )

  edge_list <- lapply(seq_along(neighbors), function(i) {
    j <- as.integer(neighbors[[i]])
    j <- j[j > i]
    if (!length(j)) return(NULL)
    data.table(from = i, to = j)
  })
  
  edge_list <- Filter(Negate(is.null), edge_list)
  
  edge_dt <- if (length(edge_list)) {
    rbindlist(edge_list, use.names = TRUE)
  } else {
    data.table(from = integer(), to = integer())
  }

  g <- igraph::make_empty_graph(n = nrow(patch_sf), directed = FALSE)
  
  if (nrow(edge_dt)) {
    g <- igraph::add_edges(g, as.vector(t(as.matrix(edge_dt[, .(from, to)]))))
  }

  comp <- igraph::components(g)$membership

  pu_dt <- data.table(
    patch_id = patch_sf$patch_id,
    pu_tmp = as.integer(comp)
  )

  pu_dt <- merge(
    pu_dt,
    patch_area_dt[, .(patch_id, patch_area_km2)],
    by = "patch_id",
    all.x = TRUE
  )

  pu_area <- pu_dt[, .(
    pu_area_km2 = sum(patch_area_km2, na.rm = TRUE)
  ), by = pu_tmp]

  keep_tmp <- pu_area[
    area_exceeds_threshold(pu_area_km2, min_pop_km2),
    pu_tmp
  ]

  pu_dt <- pu_dt[pu_tmp %in% keep_tmp]

  if (!nrow(pu_dt)) {
    return(data.table(patch_id = integer(), pu_id = integer()))
  }

  pu_map <- data.table(
    pu_tmp = sort(unique(pu_dt$pu_tmp)),
    pu_id = seq_along(sort(unique(pu_dt$pu_tmp)))
  )

  merge(pu_dt[, .(patch_id, pu_tmp)], pu_map, by = "pu_tmp")[, .(patch_id, pu_id)]
}

build_lut_one_species <- function(
  sci,
  sp_id,
  retained_species_cells,
  min_patch_km2,
  min_pop_km2,
  disp_km,
  template,
  cell_area_raster
) {
  retained_species_cells <- as.integer(retained_species_cells)

  if (!length(retained_species_cells)) {
    return(list(
      n_retained_cells = 0L,
      lut = empty_rank_lut()
    ))
  }

  binary <- make_binary_raster_from_cells(template, retained_species_cells)

  clumped <- clump_rook(binary)

  if (is.null(clumped)) {
    return(list(
      n_retained_cells = length(retained_species_cells),
      lut = empty_rank_lut()
    ))
  }

  patch_area_dt <- patch_areas_from_clump(clumped, cell_area_raster)

  patch_area_dt <- patch_area_dt[
    area_exceeds_threshold(patch_area_km2, min_patch_km2)
  ]

  if (!nrow(patch_area_dt)) {
    return(list(
      n_retained_cells = length(retained_species_cells),
      lut = empty_rank_lut()
    ))
  }

  patch_sf <- polygonize_kept_patches(clumped, patch_area_dt$patch_id)

  pu_dt <- assign_population_units(
    patch_sf = patch_sf,
    patch_area_dt = patch_area_dt,
    disp_km = disp_km,
    min_pop_km2 = min_pop_km2
  )

  if (!nrow(pu_dt)) {
    return(list(
      n_retained_cells = length(retained_species_cells),
      lut = empty_rank_lut()
    ))
  }

  out <- merge(
    patch_area_dt,
    pu_dt,
    by = "patch_id",
    all = FALSE
  )

  setorder(out, pu_id, patch_id)

  out[, patch_id := seq_len(.N)]

  out[, `:=`(
    stage = NA_integer_,
    method = "rank",
    scientificName = sci,
    species = sp_id
  )]

  setcolorder(out, c(
    "stage",
    "method",
    "scientificName",
    "species",
    "patch_id",
    "pu_id",
    "patch_area_km2",
    "patch_n_cells"
  ))

  list(
    n_retained_cells = length(retained_species_cells),
    lut = out[]
  )
}

# ============================================================
# 04_species_cell_cache - initial habitat cells by species
# ============================================================

get_species_cells_from_patch_raster <- function(path) {
  r <- terra::rast(path)
  vals <- terra::values(r, mat = FALSE)
  which(is.finite(vals) & vals > 0)
}

species_cells <- vector("list", nrow(sp))
names(species_cells) <- sp$scientificName

if (!is.null(bundle$patch_cell_index_by_species_list)) {
  cat("[species_cells] Using bundle$patch_cell_index_by_species_list where possible.\n")

  for (i in seq_len(nrow(sp))) {
    sci <- sp$scientificName[i]

    x <- bundle$patch_cell_index_by_species_list[[sci]]

    if (is.null(x)) {
      # Fallback in case names are species IDs rather than scientific names.
      x <- bundle$patch_cell_index_by_species_list[[sid(sci)]]
    }

    if (is.null(x)) {
      species_cells[[sci]] <- get_species_cells_from_patch_raster(sp$patch_path[i])
    } else {
      if (is.list(x) && !is.null(x$cell)) {
        cells <- x$cell
      } else if (is.data.frame(x) && "cell" %in% names(x)) {
        cells <- x$cell
      } else {
        cells <- x
      }
    
      species_cells[[sci]] <- sort(unique(as.integer(cells)))
    }
  }

} else {
  cat("[species_cells] Bundle has no patch_cell_index_by_species_list; reading patch rasters.\n")

  for (i in seq_len(nrow(sp))) {
    species_cells[[sp$scientificName[i]]] <- get_species_cells_from_patch_raster(sp$patch_path[i])
  }
}

cell_counts <- vapply(species_cells, length, integer(1))
assert(all(cell_counts > 0),
       "At least one retained species has zero initial patch cells.")

cat("[species_cells] Cached species cell vectors for ", length(species_cells), " species.\n", sep = "")
summary(cell_counts)

n_template_cells <- terra::ncell(template)

bad_species_cells <- names(which(vapply(
  species_cells,
  function(x) {
    any(!is.finite(x) | x < 1L | x > n_template_cells)
  },
  logical(1)
)))

assert(!length(bad_species_cells),
       paste("Some species cell vectors contain invalid cell indices:",
             paste(bad_species_cells, collapse = ", ")))

bad_outside_alive0 <- names(which(vapply(
  species_cells,
  function(x) {
    any(!alive0[x])
  },
  logical(1)
)))

assert(!length(bad_outside_alive0),
       paste("Some species cell vectors include cells outside the initial alive0 domain:",
             paste(bad_outside_alive0, collapse = ", ")))

# ============================================================
# 05_stage0_lut - shared baseline
# ============================================================

baseline_lut_path <- file.path(lut_root, sprintf("%s%04d.rds", lut_prefix, 0L))

baseline_lut <- NULL

if (overwrite || !file.exists(baseline_lut_path)) {
  baseline_lut <- as.data.table(bundle$patch_table)
  
  if (!"scientificName" %in% names(baseline_lut)) {
    if ("species" %in% names(baseline_lut)) {
      baseline_lut[, scientificName := as.character(species)]
    } else {
      stop("bundle$patch_table must contain either scientificName or species.", call. = FALSE)
    }
  }
  
  required_lut_cols <- c("scientificName", "patch_id", "pu_id", "patch_area_km2")
  missing_lut_cols <- setdiff(required_lut_cols, names(baseline_lut))
  assert(!length(missing_lut_cols),
         paste("bundle$patch_table is missing column(s):",
               paste(missing_lut_cols, collapse = ", ")))

  baseline_lut <- baseline_lut[scientificName %in% retained_species]

  baseline_lut[, species := sid(scientificName)]
  baseline_lut[, stage := 0L]
  baseline_lut[, method := "rank"]

  if (!"patch_n_cells" %in% names(baseline_lut)) {
    baseline_lut[, patch_n_cells := NA_integer_]
  }

  baseline_lut <- baseline_lut[, .(
    stage,
    method,
    scientificName,
    species,
    patch_id = as.integer(patch_id),
    pu_id = as.integer(pu_id),
    patch_area_km2 = as.numeric(patch_area_km2),
    patch_n_cells = as.integer(patch_n_cells)
  )]

  saveRDS(baseline_lut, baseline_lut_path)

  cat("[stage0] Wrote: ", baseline_lut_path, "\n", sep = "")

} else {
  cat("[stage0] Exists -> skipping: ", baseline_lut_path, "\n", sep = "")
}

# ============================================================
# 06_build_stage_luts - direct rank-map retained landscape -> LUTs
# ============================================================

# Exclude stage 0 here; it was handled above.
build_targets <- stage_targets[stage != 0L]
if (nrow(build_targets) > 1L) {
  assert(all(diff(build_targets$keep_n) <= 0),
         "Rank-map keep_n should be non-increasing across increasing stages.")
}
setorder(build_targets, stage)

# Species-level monotone cache.
prev_n_cells <- setNames(rep(NA_integer_, nrow(sp)), sp$scientificName)
prev_lut <- setNames(vector("list", nrow(sp)), sp$scientificName)
for (nm in names(prev_lut)) prev_lut[[nm]] <- empty_rank_lut()

rank_cells_sorted <- as.integer(rank_dt$cell)
n_template_cells <- terra::ncell(template)

retained_mask <- rep(FALSE, n_template_cells)
current_keep_n <- 0L

update_retained_mask <- function(next_keep_n) {
  next_keep_n <- as.integer(next_keep_n)

  assert(is.finite(next_keep_n) &&
           next_keep_n >= 0L &&
           next_keep_n <= length(rank_cells_sorted),
         "Invalid keep_n while updating retained mask.")

  if (next_keep_n > current_keep_n) {
    idx <- seq.int(current_keep_n + 1L, next_keep_n)
    retained_mask[rank_cells_sorted[idx]] <<- TRUE

  } else if (next_keep_n < current_keep_n) {
    idx <- seq.int(next_keep_n + 1L, current_keep_n)
    retained_mask[rank_cells_sorted[idx]] <<- FALSE
  }

  current_keep_n <<- next_keep_n

  invisible(NULL)
}

for (rr in seq_len(nrow(build_targets))) {
  s <- as.integer(build_targets$stage[rr])
  keep_n <- as.integer(build_targets$keep_n[rr])

  out_path <- file.path(lut_root, sprintf("%s%04d.rds", lut_prefix, s))

  if (file.exists(out_path) && !overwrite) {
    cat("[stage ", s, "] Exists -> skipping: ", out_path, "\n", sep = "")
    next
  }

  cat("\n[stage ", s, "] Building rank-map LUT; keep_n = ", keep_n, "\n", sep = "")

  update_retained_mask(keep_n)

  stage_luts <- vector("list", nrow(sp))

  for (i in seq_len(nrow(sp))) {
    sci <- sp$scientificName[i]
    sp_id <- sp$species[i]
    
    cat("[stage ", s, "] Processing species ", i, "/", nrow(sp), ": ", sci, "\n", sep = "")
    
    cells_i <- species_cells[[sci]]
    retained_species_cells <- cells_i[retained_mask[cells_i]]
    species_retained_n <- length(retained_species_cells)

    if (
      is.finite(prev_n_cells[[sci]]) &&
      identical(as.integer(prev_n_cells[[sci]]), as.integer(species_retained_n))
    ) {
      stage_luts[[i]] <- copy(prev_lut[[sci]])
      if (nrow(stage_luts[[i]])) stage_luts[[i]][, stage := s]
      next
    }

    res <- build_lut_one_species(
      sci = sci,
      sp_id = sp_id,
      retained_species_cells = retained_species_cells,
      min_patch_km2 = as.numeric(sp$min_patch_area_km2[i]),
      min_pop_km2 = as.numeric(sp$min_population_area_km2[i]),
      disp_km = as.numeric(sp$dispersal_dist[i]),
      template = template,
      cell_area_raster = cell_area_raster
    )

    lut_i <- res$lut
    if (nrow(lut_i)) lut_i[, stage := s]

    stage_luts[[i]] <- lut_i

    prev_n_cells[[sci]] <- as.integer(res$n_retained_cells)
    prev_lut[[sci]] <- copy(lut_i)
  }

  stage_lut <- rbindlist(stage_luts, use.names = TRUE, fill = TRUE)

  if (!nrow(stage_lut)) {
    stage_lut <- empty_rank_lut()
  }

  saveRDS(stage_lut, out_path)

  cat("[stage ", s, "] Wrote: ", out_path, " rows=", nrow(stage_lut), "\n", sep = "")

  rm(stage_luts, stage_lut)
  gc(FALSE)
}
