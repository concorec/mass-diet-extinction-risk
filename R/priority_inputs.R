source(file.path("R", "project_utils.R"))
source(file.path("R", "priority_run_config.R"))


initialize_priority_inputs <- function(
  curve_label = "q025",
  do_mammals = TRUE,
  do_birds = TRUE,
  do_ppm = TRUE,
  do_rangebag = TRUE,
  species_table_csv_path = file.path("Data", "Clean", "species_table.csv"),
  patch_raster_dir = file.path("Data", "Clean", "Patches"),
  patch_lookup_rds_path = file.path("Data", "Clean", "all_patch_lookup.rds"),
  connectivity_rds_path = file.path("Data", "Clean", "all_connectivity.rds"),
  pruning_input_bundle_dir = file.path("Data", "Clean", "PriorityInputs")
) {
# ---------------------------------------------------------------------
# 1) User-facing settings
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
# 3) Small helpers used only in this chunk
# ---------------------------------------------------------------------

normalize_taxon_class <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0 <- gsub("[^a-z]+", "", x0)

  out <- rep(NA_character_, length(x0))
  out[x0 %in% c("mammalia", "mammal", "mammals")] <- "mammalia"
  out[x0 %in% c("aves", "bird", "birds")]        <- "aves"

  out
}

normalize_sdm_method <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0 <- gsub("[^a-z0-9]+", "", x0)

  out <- rep(NA_character_, length(x0))
  out[x0 %in% c("ppm")] <- "ppm"
  out[x0 %in% c("rangebag", "rangebags", "range")] <- "rangebag"

  out
}

# Standardize Red List category labels from species_table.csv.
normalize_redlist_category <- function(x) {
  x0 <- tolower(trimws(as.character(x)))
  x0 <- gsub("[^a-z0-9]+", "", x0)

  out <- rep(NA_character_, length(x0))

  out[x0 %in% c("lc", "leastconcern")] <- "least_concern"
  out[x0 %in% c("nt", "nearthreatened")] <- "near_threatened"
  out[x0 %in% c("vu", "vulnerable")] <- "vulnerable"
  out[x0 %in% c("en", "endangered")] <- "endangered"
  out[x0 %in% c("cr", "criticallyendangered")] <- "critically_endangered"
  out[x0 %in% c("dd", "datadeficient")] <- "data_deficient"

  out
}

# Enforce undirected symmetry on a CSR graph.
#
# Input assumptions:
#   - row_ptr is 0-based cumulative offsets
#   - col_idx is 1-based local node indices in 1..N
#
# Output:
#   - symmetric adjacency
#   - no self-loops
#   - unique, sorted neighbors in each row
csr_symmetrize_undirected <- function(row_ptr, col_idx) {
  row_ptr <- as.integer(row_ptr)
  col_idx <- as.integer(col_idx)

  n_nodes <- length(row_ptr) - 1L

  if (n_nodes <= 0L) {
    stop("CSR graph has no rows.")
  }

  row_degree <- as.integer(diff(row_ptr))

  if (sum(row_degree) != length(col_idx)) {
    stop("Malformed CSR graph: sum(diff(row_ptr)) != length(col_idx).")
  }

  from <- rep.int(seq_len(n_nodes), row_degree)
  to   <- col_idx

  valid <- is.finite(to) & to >= 1L & to <= n_nodes & (from != to)
  from  <- from[valid]
  to    <- to[valid]

  # Add reverse edges.
  from2 <- c(from, to)
  to2   <- c(to, from)

  # Deduplicate directed edges.
  edge_key <- as.numeric(from2) + (as.numeric(to2) - 1) * as.numeric(n_nodes)
  ordering <- order(edge_key)

  edge_key <- edge_key[ordering]
  from2    <- from2[ordering]
  to2      <- to2[ordering]

  keep <- !duplicated(edge_key)
  from2 <- from2[keep]
  to2   <- to2[keep]

  adjacency_by_row <- split(to2, from2)

  full_adjacency <- vector("list", n_nodes)
  row_lengths    <- integer(n_nodes)

  for (u in seq_len(n_nodes)) {
    neighbors <- adjacency_by_row[[as.character(u)]]

    if (is.null(neighbors)) {
      full_adjacency[[u]] <- integer(0L)
      row_lengths[u] <- 0L
    } else {
      neighbors <- sort(unique(as.integer(neighbors)))
      full_adjacency[[u]] <- neighbors
      row_lengths[u] <- length(neighbors)
    }
  }

  row_ptr_new <- as.integer(c(0L, cumsum(row_lengths)))
  col_idx_new <- as.integer(unlist(full_adjacency, use.names = FALSE))

  list(
    row_ptr = row_ptr_new,
    col_idx = col_idx_new
  )
}

# Check whether a CSR graph is symmetric.
#
# Input assumptions:
#   - row_ptr is 0-based cumulative offsets
#   - col_idx is 1-based local node indices
check_csr_symmetric <- function(row_ptr, col_idx) {
  row_ptr <- as.integer(row_ptr)
  col_idx <- as.integer(col_idx)

  n_nodes <- length(row_ptr) - 1L
  if (n_nodes <= 0L) return(TRUE)

  neighbors_of <- function(u) {
    first_edge_index <- row_ptr[u] + 1L
    last_edge_index  <- row_ptr[u + 1L]

    if (last_edge_index < first_edge_index) {
      integer(0L)
    } else {
      as.integer(col_idx[first_edge_index:last_edge_index])
    }
  }

  for (u in seq_len(n_nodes)) {
    neighbors_u <- neighbors_of(u)

    if (!length(neighbors_u)) next

    for (v in neighbors_u) {
      neighbors_v <- neighbors_of(v)

      if (!(u %in% neighbors_v)) {
        return(FALSE)
      }
    }
  }

  TRUE
}


# ---------------------------------------------------------------------
# 4) Package checks
# ---------------------------------------------------------------------
load_packages(c("data.table", "terra", "stringr"))


# ---------------------------------------------------------------------
# 5) Existence checks for upstream outputs
# ---------------------------------------------------------------------
if (!file.exists(species_table_csv_path)) {
  stop("Missing species table output from script 4: ", species_table_csv_path)
}

if (!dir.exists(patch_raster_dir)) {
  stop("Missing patch raster directory from script 5: ", patch_raster_dir)
}

if (!file.exists(patch_lookup_rds_path)) {
  stop("Missing patch lookup output from script 5: ", patch_lookup_rds_path)
}

if (!file.exists(connectivity_rds_path)) {
  stop("Missing connectivity output from script 5: ", connectivity_rds_path)
}


# ---------------------------------------------------------------------
# 6) Read species table from script 4 and build species_params
# ---------------------------------------------------------------------
species_table_raw <- data.table::fread(species_table_csv_path)

required_species_columns <- c(
  "scientificName",
  "className",
  "redlistCategory",
  "sdm_method",
  "density",
  "dispersal_dist",
  "min_patch_size",
  "min_pop_size",
  paste0("alpha_", CURVE_LABEL),
  paste0("beta_", CURVE_LABEL)
)

missing_species_columns <- setdiff(required_species_columns, names(species_table_raw))
if (length(missing_species_columns)) {
  stop(
    "species_table.csv is missing required column(s): ",
    paste(missing_species_columns, collapse = ", ")
  )
}

# Standardize taxonomic class and SDM-method labels before filtering.
species_table_raw[
  ,
  taxon_class_std := normalize_taxon_class(className)
]

species_table_raw[
  ,
  sdm_method_std := normalize_sdm_method(sdm_method)
]

species_table_raw[
  ,
  redlist_category_std := normalize_redlist_category(redlistCategory)
]

unknown_taxon_values <- sort(unique(
  as.character(species_table_raw$className[
    is.na(species_table_raw$taxon_class_std) &
      !is.na(species_table_raw$className) &
      nzchar(trimws(as.character(species_table_raw$className)))
  ])
))

if (length(unknown_taxon_values)) {
  stop(
    "Unrecognized className value(s) in species_table.csv: ",
    paste(unknown_taxon_values, collapse = ", ")
  )
}

unknown_sdm_values <- sort(unique(
  as.character(species_table_raw$sdm_method[
    is.na(species_table_raw$sdm_method_std) &
      !is.na(species_table_raw$sdm_method) &
      nzchar(trimws(as.character(species_table_raw$sdm_method)))
  ])
))

unknown_redlist_values <- sort(unique(
  as.character(species_table_raw$redlistCategory[
    is.na(species_table_raw$redlist_category_std) &
      !is.na(species_table_raw$redlistCategory) &
      nzchar(trimws(as.character(species_table_raw$redlistCategory)))
  ])
))

if (length(unknown_redlist_values)) {
  stop(
    "Unrecognized redlistCategory value(s) in species_table.csv: ",
    paste(unknown_redlist_values, collapse = ", ")
  )
}

if (length(unknown_sdm_values)) {
  stop(
    "Unrecognized sdm_method value(s) in species_table.csv: ",
    paste(unknown_sdm_values, collapse = ", "),
    ". Expected values equivalent to PPM or RangeBag."
  )
}

selected_taxon_rows <- (
  (DO_MAMMALS & species_table_raw$taxon_class_std == "mammalia") |
  (DO_BIRDS   & species_table_raw$taxon_class_std == "aves")
)

selected_sdm_rows <- (
  (DO_PPM      & species_table_raw$sdm_method_std == "ppm") |
  (DO_RANGEBAG & species_table_raw$sdm_method_std == "rangebag")
)

selected_species_rows <- selected_taxon_rows & selected_sdm_rows

species_table_selected <- species_table_raw[selected_species_rows]

if (!nrow(species_table_selected)) {
  stop(
    "After applying taxon and SDM-source gates, no species remain. ",
    "Current gates: DO_MAMMALS=", DO_MAMMALS,
    ", DO_BIRDS=", DO_BIRDS,
    ", DO_PPM=", DO_PPM,
    ", DO_RANGEBAG=", DO_RANGEBAG, "."
  )
}

selection_summary <- species_table_selected[
  ,
  .N,
  by = .(taxon_class_std, sdm_method_std, redlist_category_std)
][order(taxon_class_std, sdm_method_std, redlist_category_std)]

message(
  "Species selected before spatial-output filtering: ",
  nrow(species_table_selected),
  " species across ",
  nrow(selection_summary),
  " taxon/SDM/Red List metadata groups."
)

species_table_selected <- add_canonical_area_thresholds(
  species_table_selected,
  validate = TRUE,
  label = "species_params area thresholds"
)

species_params <- species_table_selected[
  ,
  .(
    species               = as.character(scientificName),
    taxon_class           = as.character(taxon_class_std),
    redlist_category      = as.character(redlist_category_std),
    sdm_method            = as.character(sdm_method_std),
    density               = as.numeric(density),
    dispersal_distance_km = as.numeric(dispersal_dist),
    min_patch_area_km2      = as.numeric(min_patch_area_km2),
    min_population_area_km2 = as.numeric(min_population_area_km2),
    a_pred                = as.numeric(get(paste0("alpha_", CURVE_LABEL))),
    b_pred                = as.numeric(get(paste0("beta_", CURVE_LABEL)))
  )
]

if (any(!is.finite(species_params$dispersal_distance_km) |
        species_params$dispersal_distance_km <= 0)) {
  stop("species_params contains non-positive or non-finite dispersal_distance_km values.")
}

if (any(!is.finite(species_params$density) | species_params$density <= 0)) {
  stop("species_params contains non-positive or non-finite density values.")
}

if (any(!is.finite(species_params$min_patch_area_km2) |
        species_params$min_patch_area_km2 <= 0)) {
  stop("species_params contains non-positive or non-finite min_patch_area_km2 values.")
}

if (any(!is.finite(species_params$min_population_area_km2) |
        species_params$min_population_area_km2 <= 0)) {
  stop("species_params contains non-positive or non-finite min_population_area_km2 values.")
}

if (any(!is.finite(species_params$a_pred) | species_params$a_pred <= 0)) {
  stop("species_params contains non-positive or non-finite a_pred values.")
}

if (any(!is.finite(species_params$b_pred) | species_params$b_pred <= 0)) {
  stop("species_params contains non-positive or non-finite b_pred values.")
}

setkey(species_params, species)


# ---------------------------------------------------------------------
# 7) Read patch lookup from script 5 and standardize it to patch_table
# ---------------------------------------------------------------------
patch_lookup_raw <- readRDS(patch_lookup_rds_path)
patch_lookup_dt  <- as.data.table(patch_lookup_raw)

required_patch_lookup_columns <- c("scientificName", "patch_id", "pu_id", "patch_area_km2")
missing_patch_lookup_columns <- setdiff(required_patch_lookup_columns, names(patch_lookup_dt))

if (length(missing_patch_lookup_columns)) {
  stop(
    "all_patch_lookup.rds is missing required column(s): ",
    paste(missing_patch_lookup_columns, collapse = ", ")
  )
}

patch_table <- patch_lookup_dt[
  ,
  .(
    species        = as.character(scientificName),
    patch_id       = as.integer(patch_id),
    pu_id          = as.integer(pu_id),
    patch_area_km2 = as.numeric(patch_area_km2)
  )
]

patch_table <- patch_table[species %in% species_params$species]

if (!nrow(patch_table)) {
  stop("patch_table is empty after filtering to the selected species.")
}

setkey(patch_table, species, patch_id)


# ---------------------------------------------------------------------
# 8) Restrict to species that actually have final spatial outputs
# ---------------------------------------------------------------------
candidate_species <- species_params$species

patch_raster_file_by_species <- setNames(
  file.path(
    patch_raster_dir,
    vapply(candidate_species, patch_filename_from_scientific, character(1L))
  ),
  candidate_species
)

species_with_patch_raster <- names(patch_raster_file_by_species)[file.exists(patch_raster_file_by_species)]
species_with_patch_lookup <- unique(patch_table$species)

retained_test_species <- intersect(species_with_patch_raster, species_with_patch_lookup)

if (!length(retained_test_species)) {
  stop("No species remain after restricting to those with both final patch rasters and patch lookup rows.")
}

species_params <- species_params[species %in% retained_test_species]
patch_table    <- patch_table[species %in% retained_test_species]

patch_raster_file_by_species <- patch_raster_file_by_species[species_params$species]

final_selection_summary <- species_params[
  ,
  .N,
  by = .(taxon_class, sdm_method, redlist_category)
][order(taxon_class, sdm_method, redlist_category)]

message(
  "Species retained after requiring patch rasters and patch lookup rows: ",
  nrow(species_params),
  " species across ",
  nrow(final_selection_summary),
  " taxon/SDM/Red List metadata groups."
)


# ---------------------------------------------------------------------
# 9) Read final patch rasters from script 5 and stack them
# ---------------------------------------------------------------------
patch_layers <- lapply(seq_along(patch_raster_file_by_species), function(i) {
  species_name <- names(patch_raster_file_by_species)[i]
  raster_path  <- patch_raster_file_by_species[[i]]

  r <- terra::rast(raster_path)
  names(r) <- species_name
  r
})

patch_stack <- terra::rast(patch_layers)
template_raster <- patch_stack[[1]]

if (terra::nlyr(patch_stack) == 0L) {
  stop("patch_stack has zero layers.")
}


# ---------------------------------------------------------------------
# 10) Build patch_id_by_species_env
# ---------------------------------------------------------------------
patch_id_by_species_env <- new.env(parent = emptyenv())

for (species_name in names(patch_stack)) {
  patch_ids <- as.integer(terra::values(patch_stack[[species_name]], mat = FALSE))
  assign(species_name, patch_ids, envir = patch_id_by_species_env)
}


# ---------------------------------------------------------------------
# 11) Build patch_cell_index_by_species_env
# ---------------------------------------------------------------------
patch_cell_index_by_species_env <- new.env(parent = emptyenv())

for (species_name in names(patch_stack)) {
  patch_ids <- get(species_name, envir = patch_id_by_species_env, inherits = FALSE)

  occupied_cells <- which(!is.na(patch_ids))

  if (!length(occupied_cells)) {
    assign(species_name, NULL, envir = patch_cell_index_by_species_env)
    next
  }

  occupied_patch_ids <- patch_ids[occupied_cells]
  ordering <- order(occupied_patch_ids)

  assign(
    species_name,
    list(
      pid  = as.integer(occupied_patch_ids[ordering]),
      cell = as.integer(occupied_cells[ordering])
    ),
    envir = patch_cell_index_by_species_env
  )
}


# ---------------------------------------------------------------------
# 12) Derive global static raster objects
# ---------------------------------------------------------------------
cell_area_by_cell <- as.numeric(
  terra::values(
    terra::cellSize(template_raster, unit = "km"),
    mat = FALSE
  )
)

rook_neighbor_pairs <- terra::adjacent(
  template_raster,
  cells      = seq_len(terra::ncell(template_raster)),
  directions = 4,
  pairs      = TRUE
)


# ---------------------------------------------------------------------
# 13) Derive the initial alive-species count vector
# ---------------------------------------------------------------------
patch_value_matrix <- terra::values(patch_stack, mat = TRUE)

alive_species_count_by_cell <- rowSums(!is.na(patch_value_matrix))


# ---------------------------------------------------------------------
# 14) Read connectivity from script 5 and convert it to pu_graphs_by_key
# ---------------------------------------------------------------------
all_connectivity_raw <- readRDS(connectivity_rds_path)

if (!is.list(all_connectivity_raw)) {
  stop("all_connectivity.rds is not a list.")
}

csr_version   <- attr(all_connectivity_raw, "csr_version", exact = TRUE)
col_idx_base  <- attr(all_connectivity_raw, "col_idx_base", exact = TRUE)
col_idx_space <- attr(all_connectivity_raw, "col_idx_space", exact = TRUE)

if (!is.null(csr_version) && csr_version != 2L) {
  stop("Unexpected csr_version attribute in all_connectivity.rds: ", csr_version)
}

if (!is.null(col_idx_base) && col_idx_base != "0-based") {
  stop("Unexpected col_idx_base attribute in all_connectivity.rds: ", col_idx_base)
}

if (!is.null(col_idx_space) && col_idx_space != "pu_local_index") {
  stop("Unexpected col_idx_space attribute in all_connectivity.rds: ", col_idx_space)
}

all_connectivity_filtered <- Filter(
  f = function(x) {
    is.list(x) &&
      !is.null(x$species) &&
      x$species %in% species_params$species
  },
  x = all_connectivity_raw
)

pu_graphs_by_key <- list()

for (e in all_connectivity_filtered) {
  id2patch <- as.integer(e$patch_ids)

  # Stored col_idx is 0-based PU-local indexing; convert to 1-based for R.
  col_idx <- as.integer(e$col_idx) + 1L

  # Keep row_ptr as 0-based cumulative offsets.
  row_ptr <- as.integer(e$row_ptr)

  sym <- csr_symmetrize_undirected(
    row_ptr = row_ptr,
    col_idx = col_idx
  )

  graph_key <- paste0(as.character(e$species), "|", as.integer(e$pu_id))

  pu_graphs_by_key[[graph_key]] <- list(
    species  = as.character(e$species),
    pu_id    = as.integer(e$pu_id),
    id2patch = id2patch,
    row_ptr  = sym$row_ptr,
    col_idx  = sym$col_idx
  )
}


# ---------------------------------------------------------------------
# 15) Internal validation checks
# ---------------------------------------------------------------------
if (anyDuplicated(patch_table[, .(species, patch_id)]) > 0L) {
  stop("patch_table contains duplicate (species, patch_id) rows.")
}

if (!all(unique(patch_table$species) %in% species_params$species)) {
  stop("Some species in patch_table are missing from species_params.")
}

if (!all(unique(patch_table$species) %in% names(patch_stack))) {
  stop("Some species in patch_table are missing from patch_stack.")
}

if (!all(alive_species_count_by_cell == rowSums(!is.na(patch_value_matrix)))) {
  stop("alive_species_count_by_cell does not match the non-NA count in patch_stack.")
}

for (species_name in names(patch_stack)) {
  raster_patch_ids <- sort(unique(as.integer(patch_value_matrix[, species_name][!is.na(patch_value_matrix[, species_name])])))
  table_patch_ids  <- sort(unique(as.integer(patch_table$patch_id[patch_table$species == species_name])))

  if (!identical(raster_patch_ids, table_patch_ids)) {
    stop(
      "Patch ID mismatch for species '", species_name,
      "': raster IDs do not match patch_table IDs."
    )
  }
}

for (graph_key in names(pu_graphs_by_key)) {
  g <- pu_graphs_by_key[[graph_key]]

  valid_patch_ids <- patch_table$patch_id[patch_table$species == g$species]

  if (!all(g$id2patch %in% valid_patch_ids)) {
    stop("Graph ", graph_key, " references patch IDs not found in patch_table.")
  }

  if (!check_csr_symmetric(g$row_ptr, g$col_idx)) {
    stop("Graph ", graph_key, " is still asymmetric after symmetrization.")
  }
}


# ---------------------------------------------------------------------
# 16) Save a reusable pruning-input bundle BEFORE any pruning call
# ---------------------------------------------------------------------
#
# This bundle captures a clean starting state for repeated tests in future
# R sessions. The two environments are saved as named lists, then rebuilt
# on load. That is simpler and more robust than saving the environments
# directly.
#
dir.create(pruning_input_bundle_dir, recursive = TRUE, showWarnings = FALSE)

pruning_input_bundle_path <- file.path(
  pruning_input_bundle_dir,
  priority_bundle_name(
    curve = CURVE_LABEL,
    taxa = taxon_selection_tag,
    sdm = sdm_selection_tag
  )
)

patch_id_by_species_list <- setNames(
  lapply(ls(envir = patch_id_by_species_env, all.names = TRUE), function(species_name) {
    get(species_name, envir = patch_id_by_species_env, inherits = FALSE)
  }),
  ls(envir = patch_id_by_species_env, all.names = TRUE)
)

patch_cell_index_by_species_list <- setNames(
  lapply(ls(envir = patch_cell_index_by_species_env, all.names = TRUE), function(species_name) {
    get(species_name, envir = patch_cell_index_by_species_env, inherits = FALSE)
  }),
  ls(envir = patch_cell_index_by_species_env, all.names = TRUE)
)

pruning_input_bundle <- list(
  metadata = list(
    schema_version            = priority_bundle_schema_version(),
    created_at                = as.character(Sys.time()),
    curve_label               = CURVE_LABEL,
    do_mammals                = DO_MAMMALS,
    do_birds                  = DO_BIRDS,
    do_ppm                    = DO_PPM,
    do_rangebag               = DO_RANGEBAG,
    retained_species          = species_params$species,
    retained_species_table    = species_params[
      ,
      .(
        species          = species,
        taxon_class      = taxon_class,
        redlist_category = redlist_category,
        sdm_method       = sdm_method
      )
    ]
  ),
  patch_table                      = patch_table,
  pu_graphs_by_key                 = pu_graphs_by_key,
  alive_species_count_by_cell      = alive_species_count_by_cell,
  cell_area_by_cell                = cell_area_by_cell,
  species_params                   = species_params,
  rook_neighbor_pairs              = rook_neighbor_pairs,
  patch_id_by_species_list         = patch_id_by_species_list,
  patch_cell_index_by_species_list = patch_cell_index_by_species_list
)

saveRDS(pruning_input_bundle, pruning_input_bundle_path)

message("Saved pruning-input bundle to: ", pruning_input_bundle_path)

invisible(list(
  pruning_input_bundle_path = pruning_input_bundle_path,
  pruning_input_bundle = pruning_input_bundle,
  patch_table = patch_table,
  pu_graphs_by_key = pu_graphs_by_key,
  alive_species_count_by_cell = alive_species_count_by_cell,
  cell_area_by_cell = cell_area_by_cell,
  species_params = species_params,
  rook_neighbor_pairs = rook_neighbor_pairs,
  patch_id_by_species_env = patch_id_by_species_env,
  patch_cell_index_by_species_env = patch_cell_index_by_species_env
))
}
