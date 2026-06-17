# Fixed scientific and artifact contracts shared across pipeline stages.

if (!exists("assert", mode = "function")) {
  source(file.path("R", "project_utils.R"))
}

analysis_constants <- function() {
  list(
    quasi_extinction_abundance = 500L,
    minimum_patch_abundance = 10L,
    persistence_horizon_years = 100L,
    priority_bundle_schema_version = 2L
  )
}

quasi_extinction_abundance <- function() {
  analysis_constants()$quasi_extinction_abundance
}

minimum_patch_abundance <- function() {
  analysis_constants()$minimum_patch_abundance
}

persistence_horizon_years <- function() {
  analysis_constants()$persistence_horizon_years
}

priority_bundle_schema_version <- function() {
  analysis_constants()$priority_bundle_schema_version
}

persistence_quantiles <- function() {
  c(q50 = 0.50, q16 = 0.16, q025 = 0.025)
}

persistence_curves <- function() {
  names(persistence_quantiles())
}

validate_persistence_curve <- function(curve, label = "curve") {
  curve <- as.character(curve)
  assert(
    length(curve) == 1L && !is.na(curve) && curve %in% persistence_curves(),
    paste0(
      label, " must be one of: ",
      paste(persistence_curves(), collapse = ", "), "."
    )
  )
  curve
}

benchmark_rank_labels <- function() {
  c(
    abf = "Zonation ABF",
    caz1 = "Zonation CAZ1",
    caz2 = "Zonation CAZ2",
    cazmax = "Zonation CAZMAX"
  )
}

benchmark_rank_methods <- function() {
  names(benchmark_rank_labels())
}

normalize_benchmark_rank_method <- function(method) {
  tolower(trimws(as.character(method)))
}

validate_benchmark_rank_method <- function(method, label = "rank_method") {
  method <- normalize_benchmark_rank_method(method)
  assert(
    length(method) == 1L && !is.na(method) && method %in% benchmark_rank_methods(),
    paste0(
      label, " must be one of: ",
      paste(benchmark_rank_methods(), collapse = ", "), "."
    )
  )
  method
}

benchmark_rank_label <- function(method) {
  method <- validate_benchmark_rank_method(method)
  unname(benchmark_rank_labels()[method])
}

area_exceeds_threshold <- function(area, threshold) {
  is.finite(area) & is.finite(threshold) & area > threshold
}

validate_area_threshold <- function(threshold, label = "area threshold") {
  assert(
    length(threshold) == 1L && is.finite(threshold) && threshold > 0,
    paste0(label, " must be one finite positive value.")
  )
  as.numeric(threshold)
}

validate_density_area_thresholds <- function(
  density,
  min_patch_area_km2,
  min_population_area_km2,
  label = "species thresholds",
  tolerance = sqrt(.Machine$double.eps)
) {
  expected_patch <- minimum_patch_abundance() / density
  expected_population <- quasi_extinction_abundance() / density

  ok <- all(is.finite(density) & density > 0) &&
    all(is.finite(min_patch_area_km2) & min_patch_area_km2 > 0) &&
    all(is.finite(min_population_area_km2) & min_population_area_km2 > 0) &&
    isTRUE(all.equal(
      as.numeric(min_patch_area_km2),
      as.numeric(expected_patch),
      tolerance = tolerance,
      check.attributes = FALSE
    )) &
    isTRUE(all.equal(
      as.numeric(min_population_area_km2),
      as.numeric(expected_population),
      tolerance = tolerance,
      check.attributes = FALSE
    ))

  assert(
    ok,
    paste0(
      label, " are inconsistent with density and the canonical abundance ",
      "thresholds (", minimum_patch_abundance(), " and ",
      quasi_extinction_abundance(), " individuals)."
    )
  )
  invisible(TRUE)
}

add_canonical_area_thresholds <- function(
  x,
  density_col = "density",
  min_patch_col = "min_patch_size",
  min_population_col = "min_pop_size",
  validate = TRUE,
  label = "species-table area thresholds"
) {
  need_cols(x, c(density_col, min_patch_col, min_population_col), label)

  x$min_patch_area_km2 <- suppressWarnings(as.numeric(x[[min_patch_col]]))
  x$min_population_area_km2 <- suppressWarnings(as.numeric(x[[min_population_col]]))

  if (isTRUE(validate)) {
    validate_density_area_thresholds(
      density = suppressWarnings(as.numeric(x[[density_col]])),
      min_patch_area_km2 = x$min_patch_area_km2,
      min_population_area_km2 = x$min_population_area_km2,
      label = label
    )
  }

  x
}

validate_priority_bundle_schema <- function(bundle, label = "script-6 pruning-input bundle") {
  found <- suppressWarnings(as.integer(bundle$metadata$schema_version))
  expected <- priority_bundle_schema_version()

  assert(
    length(found) == 1L && !is.na(found) && found == expected,
    paste0(
      "Obsolete ", label, " schema. Expected version ", expected, "; found ",
      if (length(found) && !is.na(found)) found else "an unversioned legacy bundle",
      ". Regenerate the bundle with 6_spatial_prioritization_pipeline.Rmd."
    )
  )
  invisible(TRUE)
}
