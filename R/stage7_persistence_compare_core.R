# Core setup and persistence-trajectory computation for 7.3_persist_cmp.Rmd.
# Produces pu_long, sp_long, stage_sum, and stage_compare outputs.


suppressPackageStartupMessages({
  library(data.table)
})

source(file.path("R", "project_utils.R"))
source(file.path("R", "priority_run_config.R"))

sid <- species_id

# ============================================================
# Parameters and paths
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

rank_map_subdir <- as.character(params$rank_map_subdir)

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
species_csv <- as.character(params$species_csv)

ana <- run_paths$ana
stage_meta_csv <- file.path(ana, as.character(params$stage_meta_name))

pipe_lookup_dir <- file.path(out_curve, as.character(params$pipe_patch_lookup_dir))
pipe_prefix <- as.character(params$pipe_patch_lookup_prefix)

benchmark_paths <- benchmark_artifact_paths(
  run_paths,
  rank_method,
  rank_map_subdir = rank_map_subdir
)
rank_label <- benchmark_paths$rank_label
rank_lut_root <- benchmark_paths$lut_root
rank_prefix <- as.character(params$rank_lut_prefix)

rank_path <- benchmark_paths$rank_path
target_summary_csv <- benchmark_paths$target_summary_csv

comparison_paths <- persistence_comparison_paths(run_paths, rank_method)
out_root <- comparison_paths$out_root
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

out_pu <- comparison_paths$pu
out_sp <- comparison_paths$species
out_stage <- comparison_paths$stage
out_compare <- comparison_paths$compare

need_dir(out_curve, "pipeline_output_dir")
need_file(bundle_path, "pruning input bundle")
need_file(species_csv, "species_table.csv")
need_file(stage_meta_csv, "stage_meta.csv from 7.1")
need_dir(pipe_lookup_dir, "optimized script-6 patch_lookup_tables")
need_dir(rank_lut_root, "rank-map LUT output from 7.2")
need_dir(out_root, "persist_cmp output directory")

cat("\n=== 7.3_persist_cmp ===\n")
cat("run_id      :", run_id, "\n")
cat("rank_method :", rank_method, "\n")
cat("rank_label  :", rank_label, "\n")
cat("out_curve   :", out_curve, "\n")
cat("rank_lut    :", rank_lut_root, "\n")
cat("out_root    :", out_root, "\n\n")

# ============================================================
# Load bundle and stage metadata
# ============================================================

bundle <- readRDS(bundle_path)
validate_priority_bundle_schema(bundle)

required_bundle_fields <- c("metadata", "patch_table")
missing_bundle_fields <- setdiff(required_bundle_fields, names(bundle))
assert(!length(missing_bundle_fields),
       paste("Bundle is missing field(s):",
             paste(missing_bundle_fields, collapse = ", ")))

assert(!is.null(bundle$metadata$retained_species),
       "Bundle metadata is missing retained_species.")

retained_species <- as.character(bundle$metadata$retained_species)

stage_meta <- fread(stage_meta_csv)

required_stage_cols <- c(
  "stage",
  "alive_end",
  "area_retained_km2",
  "pct_cells_removed_end"
)

missing_stage_cols <- setdiff(required_stage_cols, names(stage_meta))
assert(!length(missing_stage_cols),
       paste("stage_meta.csv is missing column(s):",
             paste(missing_stage_cols, collapse = ", ")))

assert(stage_meta[stage == 0L, .N] == 1L,
       "stage_meta.csv must contain exactly one stage 0 row.")

stage_meta[, stage := as.integer(stage)]
setorder(stage_meta, stage)

stage_meta[, stage_order := seq_len(.N) - 1L]
stage_meta[, x := as.numeric(pct_cells_removed_end)]

assert(all(is.finite(stage_meta$x)),
       "Selected stage metadata has non-finite x values.")

# ============================================================
# Load species parameters
# ============================================================

a_col <- paste0("alpha_", curve)
b_col <- paste0("beta_", curve)

species_header <- names(fread(species_csv, nrows = 0))

required_species_cols <- c(
  "scientificName",
  "className",
  "density",
  "min_patch_size",
  "min_pop_size",
  a_col,
  b_col
)

missing_species_cols <- setdiff(required_species_cols, species_header)
assert(!length(missing_species_cols),
       paste("species_table.csv is missing column(s):",
             paste(missing_species_cols, collapse = ", ")))

sp0 <- fread(species_csv, select = required_species_cols)

sp0[, scientificName := trimws(as.character(scientificName))]
sp0 <- sp0[scientificName %in% retained_species]

assert(nrow(sp0) == length(retained_species),
       paste0(
         "species_table.csv does not contain exactly the retained species from the bundle. ",
         "Expected ", length(retained_species), ", got ", nrow(sp0), "."
       ))

sp0[, species := sid(scientificName)]

sp0 <- as.data.table(add_canonical_area_thresholds(
  sp0,
  validate = FALSE,
  label = "stage-7 species-table area thresholds"
))
sp0[, density := suppressWarnings(as.numeric(density))]
sp0[, c_th := min_population_area_km2]
sp0[, a_pred := suppressWarnings(as.numeric(get(a_col)))]
sp0[, b_pred := suppressWarnings(as.numeric(get(b_col)))]
sp0[, className := as.character(className)]

bad_species <- sp0[
  !is.finite(density) | density <= 0 |
    !is.finite(min_patch_area_km2) | min_patch_area_km2 <= 0 |
    !is.finite(c_th) | c_th <= 0 |
    !is.finite(a_pred) | a_pred <= 0 |
    !is.finite(b_pred) | b_pred <= 0
]

assert(
  nrow(bad_species) == 0L,
  paste0(
    "Found ", nrow(bad_species),
    " species with invalid density/min_pop_size/Gompertz parameters."
  )
)

validate_density_area_thresholds(
  density = sp0$density,
  min_patch_area_km2 = sp0$min_patch_area_km2,
  min_population_area_km2 = sp0$c_th,
  label = "stage-7 species-table area thresholds"
)

assert(nrow(sp0) > 0L,
       "No valid retained species remain after filtering invalid parameters.")

par <- sp0[, .(
  scientificName,
  species,
  className,
  density,
  a_pred,
  b_pred,
  c_th
)]

assert(anyDuplicated(par$species) == 0L,
       "Duplicate species IDs after sid(scientificName).")

setkey(par, species)

cat("\n=== 7.3_persist_cmp ===\n")
cat("curve              :", curve, "\n")
cat("pipeline_output_dir:", out_curve, "\n")
cat("ana                :", ana, "\n")
cat("stages             :", nrow(stage_meta), "\n")
cat("retained species   :", nrow(par), "\n")
cat("x_axis             : pct_removed\n")
cat("out_root           :", out_root, "\n\n")

# ============================================================
# 01_discover_method_inputs
#
# Persistence-based method:
#   stage 0: bundle$patch_table
#   stage >0: optimized script-6 patch_lookup_tables/stage_patch_lookup_stage_####.csv
#
# Rank-map benchmark:
#   stage >=0: 7.2 ana/rank_lut/lut_stage_####.rds
# ============================================================

pipe_stage_path <- function(stage) {
  stage <- as.integer(stage)

  if (stage == 0L) {
    return(NA_character_)
  }

  # Prefer the exact path recorded by 7.1 stage_meta if present.
  if ("patch_lookup_path" %in% names(stage_meta)) {
    p <- stage_meta[stage == !!stage, patch_lookup_path][1]
    if (!is.na(p) && nzchar(p) && file.exists(p)) {
      return(p)
    }
  }

  file.path(
    pipe_lookup_dir,
    sprintf("%s%04d.csv", pipe_prefix, stage)
  )
}

rank_stage_path <- function(stage) {
  file.path(
    rank_lut_root,
    sprintf("%s%04d.rds", rank_prefix, as.integer(stage))
  )
}

pipe_st <- copy(stage_meta)[
  ,
  .(
    method = "pipe",
    stage,
    stage_order,
    x,
    path = vapply(stage, pipe_stage_path, character(1))
  )
]

rank_st <- copy(stage_meta)[
  ,
  .(
    method = "rank",
    stage,
    stage_order,
    x,
    path = vapply(stage, rank_stage_path, character(1))
  )
]

missing_pipe <- pipe_st[stage != 0L & !file.exists(path)]
assert(nrow(missing_pipe) == 0L,
       paste("Missing persistence-based stage patch lookup file(s):",
             paste(missing_pipe$path, collapse = "\n")))

missing_rank <- rank_st[!file.exists(path)]
assert(nrow(missing_rank) == 0L,
       paste("Missing rank-map LUT file(s) from 7.2:",
             paste(missing_rank$path, collapse = "\n")))

cat("[inputs] Persistence-based stages:", nrow(pipe_st), "\n")
cat("[inputs] Rank-map stages         :", nrow(rank_st), "\n")

# ============================================================
# 02_persistence_helpers
# ============================================================

load_pipe_lut <- function(stage, path) {
  stage <- as.integer(stage)

  if (stage == 0L) {
    return(standardize_rank_lut(
      x = bundle$patch_table,
      method = "pipe",
      stage = 0L,
      species_params = par
    ))
  }

  standardize_rank_lut(
    x = fread(path),
    method = "pipe",
    stage = stage,
    species_params = par
  )
}

load_rank_lut <- function(stage, path) {
  standardize_rank_lut(
    x = readRDS(path),
    method = "rank",
    stage = as.integer(stage),
    species_params = par
  )
}

# PU-level persistence probability.
# The curve is P(K) = exp[-alpha * (K - K0)^(-beta)].
# Since K = density * area and c_th = K0 / density,
# this is equivalent to:
# P(A) = exp[-alpha * density^(-beta) * (A - c_th)^(-beta)].
pu_persist <- function(area_km2, density, a_pred, b_pred, c_th) {
  A <- as.numeric(area_km2)
  dens <- as.numeric(density)
  a <- as.numeric(a_pred)
  b <- as.numeric(b_pred)
  c <- as.numeric(c_th)

  delta <- A - c
  out <- numeric(length(A))

  ok <- is.finite(delta) & delta > 0 &
    is.finite(dens) & dens > 0 &
    is.finite(a) & a > 0 &
    is.finite(b) & b > 0

  out[!ok] <- 0

  if (any(ok)) {
    log_inner <- log(a[ok]) -
      b[ok] * log(dens[ok]) -
      b[ok] * log(delta[ok])

    inner <- exp(log_inner)
    out[ok] <- exp(-inner)
  }

  out[!is.finite(out)] <- 0
  out[out < 0] <- 0
  out[out > 1] <- 1

  out
}

sp_persist <- function(P_pu) {
  P <- as.numeric(P_pu)
  P <- P[is.finite(P)]

  if (!length(P)) return(0)

  P <- pmin(pmax(P, 0), 1)

  if (any(P >= 1)) return(1)

  log_ext <- sum(log1p(-P))
  S <- 1 - exp(log_ext)

  if (!is.finite(S)) S <- 0

  min(max(S, 0), 1)
}

# ============================================================
# 03_compute_method
# ============================================================

compute_method <- function(stage_table, method, load_fn, par) {
  pu_list <- vector("list", nrow(stage_table))
  sp_list <- vector("list", nrow(stage_table))

  par_small <- par[, .(
    scientificName,
    species,
    className,
    density,
    a_pred,
    b_pred,
    c_th
  )]

  for (i in seq_len(nrow(stage_table))) {
    s <- as.integer(stage_table$stage[i])
    so <- as.integer(stage_table$stage_order[i])
    x <- as.numeric(stage_table$x[i])
    path <- stage_table$path[i]

    cat(sprintf("[%s] stage %04d (%d/%d)\n", method, s, i, nrow(stage_table)))

    lut <- load_fn(s, path)

    if (nrow(lut)) {
      pu <- lut[
        ,
        .(
          pu_area_km2 = sum(patch_area_km2, na.rm = TRUE),
          n_patches = .N
        ),
        by = .(scientificName, species, pu_id)
      ]

      pu <- merge(
        pu,
        par_small,
        by = c("scientificName", "species"),
        all.x = TRUE,
        sort = FALSE
      )

      assert(all(!is.na(pu$density)),
             paste0("Missing species parameters after joining LUT for method=",
                    method, ", stage=", s, "."))

      pu[, P_pu := pu_persist(
        area_km2 = pu_area_km2,
        density = density,
        a_pred = a_pred,
        b_pred = b_pred,
        c_th = c_th
      )]

      pu_keep <- pu[, .(
        method = method,
        stage = s,
        stage_order = so,
        x = x,
        scientificName,
        species,
        pu_id,
        pu_area_km2,
        n_patches,
        P_pu
      )]

      sp_stage <- pu_keep[
        ,
        .(
          sp_persist = sp_persist(P_pu),
          n_pu = .N,
          total_pu_area_km2 = sum(pu_area_km2, na.rm = TRUE)
        ),
        by = .(scientificName, species)
      ]

    } else {
      pu_keep <- data.table(
        method = character(),
        stage = integer(),
        stage_order = integer(),
        x = numeric(),
        scientificName = character(),
        species = character(),
        pu_id = integer(),
        pu_area_km2 = numeric(),
        n_patches = integer(),
        P_pu = numeric()
      )

      sp_stage <- data.table(
        scientificName = character(),
        species = character(),
        sp_persist = numeric(),
        n_pu = integer(),
        total_pu_area_km2 = numeric()
      )
    }

    sp_full <- merge(
      par[, .(scientificName, species, className)],
      sp_stage,
      by = c("scientificName", "species"),
      all.x = TRUE,
      sort = FALSE
    )

    sp_full[is.na(sp_persist), sp_persist := 0]
    sp_full[is.na(n_pu), n_pu := 0L]
    sp_full[is.na(total_pu_area_km2), total_pu_area_km2 := 0]

    sp_full[, `:=`(
      method = method,
      stage = s,
      stage_order = so,
      x = x
    )]

    setcolorder(sp_full, c(
      "method",
      "stage",
      "stage_order",
      "x",
      "scientificName",
      "species",
      "className",
      "sp_persist",
      "n_pu",
      "total_pu_area_km2"
    ))

    pu_list[[i]] <- pu_keep
    sp_list[[i]] <- sp_full

    rm(lut, pu_keep, sp_stage, sp_full)
    if (exists("pu")) rm(pu)
    gc(FALSE)
  }

  list(
    pu_long = rbindlist(pu_list, use.names = TRUE, fill = TRUE),
    sp_long = rbindlist(sp_list, use.names = TRUE, fill = TRUE)
  )
}

# ============================================================
# 04_run_write
# ============================================================

all_outputs_exist <- all(file.exists(c(out_pu, out_sp, out_stage, out_compare)))

if (!overwrite && all_outputs_exist) {
  cat("[persist_cmp] Outputs exist and overwrite=FALSE -> loading existing outputs.\n")
  pu_long <- fread(out_pu)
  sp_long <- fread(out_sp)
  stage_sum <- fread(out_stage)
  stage_compare <- fread(out_compare)

} else {
  res_pipe <- compute_method(
    stage_table = pipe_st,
    method = "pipe",
    load_fn = load_pipe_lut,
    par = par
  )

  res_rank <- compute_method(
    stage_table = rank_st,
    method = "rank",
    load_fn = load_rank_lut,
    par = par
  )

  pu_long <- rbindlist(
    list(res_pipe$pu_long, res_rank$pu_long),
    use.names = TRUE,
    fill = TRUE
  )

  sp_long <- rbindlist(
    list(res_pipe$sp_long, res_rank$sp_long),
    use.names = TRUE,
    fill = TRUE
  )

  # Fixed denominator: all retained species appear in sp_long at every stage.
  stage_sum <- sp_long[
    ,
    .(
      n_species = .N,
      n_species_with_pu = sum(n_pu > 0),
      mean_persist = mean(sp_persist),
      median_persist = median(sp_persist),
      q10_persist = as.numeric(quantile(sp_persist, 0.10, names = FALSE)),
      q90_persist = as.numeric(quantile(sp_persist, 0.90, names = FALSE)),
      min_persist = min(sp_persist),
      max_persist = max(sp_persist),
      frac_gt_05 = mean(sp_persist > 0.5),
      total_pu_area_km2 = sum(total_pu_area_km2, na.rm = TRUE)
    ),
    by = .(method, stage, stage_order, x)
  ]

  setorder(stage_sum, method, stage_order)

  stage_compare <- dcast(
    stage_sum[
      ,
      .(
        stage,
        stage_order,
        x,
        method,
        mean_persist,
        median_persist,
        frac_gt_05,
        total_pu_area_km2
      )
    ],
    stage + stage_order + x ~ method,
    value.var = c(
      "mean_persist",
      "median_persist",
      "frac_gt_05",
      "total_pu_area_km2"
    )
  )

  if (all(c("mean_persist_pipe", "mean_persist_rank") %in% names(stage_compare))) {
    stage_compare[
      ,
      mean_persist_diff_pipe_minus_rank :=
        mean_persist_pipe - mean_persist_rank
    ]
  }

  if (all(c("median_persist_pipe", "median_persist_rank") %in% names(stage_compare))) {
    stage_compare[
      ,
      median_persist_diff_pipe_minus_rank :=
        median_persist_pipe - median_persist_rank
    ]
  }

  if (all(c("total_pu_area_km2_pipe", "total_pu_area_km2_rank") %in% names(stage_compare))) {
    stage_compare[
      ,
      total_pu_area_diff_pipe_minus_rank :=
        total_pu_area_km2_pipe - total_pu_area_km2_rank
    ]
  }

  setorder(stage_compare, stage_order)

  fwrite(pu_long, out_pu)
  fwrite(sp_long, out_sp)
  fwrite(stage_sum, out_stage)
  fwrite(stage_compare, out_compare)

  cat("\n[persist_cmp] Wrote:\n")
  cat("  - ", out_pu, "\n", sep = "")
  cat("  - ", out_sp, "\n", sep = "")
  cat("  - ", out_stage, "\n", sep = "")
  cat("  - ", out_compare, "\n", sep = "")

  rm(res_pipe, res_rank)
  gc(FALSE)
}

