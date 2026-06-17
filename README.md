# Persistence-Based Spatial Conservation Prioritization

This repository contains the complete analysis workflow for a persistence-based spatial conservation prioritization study in Madagascar. The workflow links demographic allometry, stochastic population simulation, shifted-Gompertz persistence curves, species traits, species distribution rasters, habitat suitability, spatial patches, dispersal-connected population units, a reverse-removal prioritization algorithm, and comparisons with four benchmark rank maps.

The project is organized as numbered R Markdown drivers and focused helper files in `R/`. The R Markdown files define the scientific order of operations and expose the parameters that are intended to be changed between runs. The helper files contain the reusable calculations, validation rules, spatial data structures, simulation routines, output writers, reporting code, and figure builders.

This README is both a user guide and a technical description of how the implementation works. It follows the pipeline chronologically, describes the helper functions called by each driver, and records the contracts that connect one stage to the next.

## Contents

1. [Execution Model](#execution-model)
2. [Repository Layout](#repository-layout)
3. [Software Requirements](#software-requirements)
4. [Pipeline Overview](#pipeline-overview)
5. [Shared Conventions](#shared-conventions)
6. [Stage 1: Demographic Calibration](#stage-1-demographic-calibration)
7. [Stage 2: Persistence-Point Simulation](#stage-2-persistence-point-simulation)
8. [Stage 3: Gompertz and LOESS Models](#stage-3-gompertz-and-loess-models)
9. [Stage 4: Canonical Species Table](#stage-4-canonical-species-table)
10. [Stage 5: Habitat Patches and Connectivity](#stage-5-habitat-patches-and-connectivity)
11. [Stage 5.1: Single-Species Process Figure](#stage-51-single-species-process-figure)
12. [Stage 5.2: Zonation Benchmark Inputs](#stage-52-zonation-benchmark-inputs)
13. [Stage 6: Spatial Prioritization](#stage-6-spatial-prioritization)
14. [Stage 7.1: Stage Metadata and Priority Surface](#stage-71-stage-metadata-and-priority-surface)
15. [Stage 7.2: Benchmark Rank Lookup Tables](#stage-72-benchmark-rank-lookup-tables)
16. [Stage 7.3: Persistence Comparison](#stage-73-persistence-comparison)
17. [Manuscript Figures](#manuscript-figures)
18. [Output Contracts](#output-contracts)

## Execution Model

### Working directory

Run every R Markdown driver from the repository root. Open the `.Rmd` file in
RStudio and run it from the project root. The YAML
`params:` block at the top of each Rmd is the intended place to edit run
settings before knitting or running chunks.

Running a driver from another working directory will generally cause file-existence checks to fail.

### Drivers and helpers

The numbered `.Rmd` files are the public entry points. Each driver:

1. Reads YAML parameters.
2. Sources configuration and utility helpers.
3. Validates packages, files, schemas, and parameter values.
4. Calls the stage-specific computational helpers.
5. Writes canonical outputs.
6. Produces a concise knitted summary and, where applicable, manuscript figures.

Most scientific logic is in `R/`. This separation is deliberate: the Rmd files show the chronological workflow, while helper files make the calculations testable and reusable.

Helper files are grouped by stage and by shared contract:

- `project_utils.R` provides general assertions, package loading, path checks,
  species-safe IDs, and small numerical helpers.
- `analysis_contract.R` holds fixed scientific constants, curve labels,
  benchmark labels, threshold semantics, and bundle schema validation.
- `priority_run_config.R` centralizes script-6/7 run IDs, benchmark artifact
  paths, comparison paths, and rank-LUT schemas.
- Stage-specific prefixes (`rm_sigma_`, `persist_`, `gompertz_`,
  `species_table_`, `patch_`, `priority_`, and `stage7_`) keep most logic near
  the driver that uses it.

When reviewing code, start from the numbered Rmd, then follow the sourced files
in order. That mirrors the runtime order and avoids mistaking helper definition
order for scientific execution order.

### Parameter handling

The Rmd YAML blocks are the parameter source of truth. Helper files may convert
those visible values into structured objects, but helper files no longer carry
alternate defaults for stage-specific run parameters.

Interactive chunk execution deserves care. If a chunk is run outside the normal
RStudio Rmd context, define `params` from the YAML block first.

### Expensive and destructive operations

Several stages overwrite canonical outputs when their overwrite or run flags permit it:

- Stage 1 can refit JAGS models and replace posterior draw CSVs.
- Stage 2 can rebuild the C++ cache and replace or append persistence-point output.
- Stage 5 can overwrite species patch rasters and the combined lookup/connectivity RDS files.
- Stage 5.2 can overwrite binary patch rasters, Zonation feature lists, and staged benchmark rank maps.
- Stage 6 writes a complete run directory for the selected run ID.
- Stages 7.1-7.3 can replace analysis products when `overwrite = true`.

Review YAML parameters before rendering. In particular, `6_spatial_prioritization_pipeline.Rmd` currently defaults to:

```yaml
initialize_inputs: false
run_pipeline: false
```

Therefore, rendering script 6 with no parameter override reports the selected
configuration without rebuilding the saved input bundle or running the full
priority pipeline.

### Generated data and version consistency

Downstream outputs are not interchangeable across arbitrary upstream versions. A
saved bundle, rank map, patch lookup, or persistence-point table can have the
correct filename while representing an older domain or parameterization. The
safest workflow is to regenerate downstream products after changing their inputs
or schema contracts rather than mixing artifacts from different runs. Existing
files are not deleted automatically.

## Repository Layout

```text
.
|-- 1_rm_sigma_models.Rmd
|-- 2_persist_points_crn.Rmd
|-- 3_gompertz_loess_models.Rmd
|-- 4_build_species_table.Rmd
|-- 5_build_patches_and_connectivity.Rmd
|-- 5.1_single_species_aoh_patches_pu_process.Rmd
|-- 5.2_build_binary_patch_rasters_zonation_feature_list.Rmd
|-- 6_spatial_prioritization_pipeline.Rmd
|-- 7.1_stage_meta_and_priority_surface.Rmd
|-- 7.2_build_rank_lut.Rmd
|-- 7.3_persist_cmp.Rmd
|-- Main_Text.tex
|-- Supplementary_Methods.tex
|-- Title_Page.tex
|-- references.bib
|-- R/
|   |-- project_utils.R
|   |-- analysis_contract.R
|   |-- figure_utils.R
|   |-- rm_sigma_*.R
|   |-- persist_*.R
|   |-- gompertz_*.R
|   |-- species_table_*.R
|   |-- patch_*.R
|   |-- priority_*.R
|   `-- stage7_*.R
|-- src/
|   `-- simulate_persist_probs_cpp.cpp
|-- Data/
|   |-- Raw/
|   |-- Clean/
|   `-- Results/
|-- Figures/
|-- bird_ppm_bin/
|-- bird_rangebag_bin/
|-- mammal_ppm_bin/
`-- mammal_rangebag_bin/
```

### Directory roles

| Path | Role |
|---|---|
| `Data/Raw/` | Source tabular and raster inputs that are not generated by this pipeline. |
| `Data/Clean/` | Canonical cleaned tables, fitted model objects, patch rasters, connectivity objects, saved priority-input bundles, and benchmark staging products. |
| `Data/Results/` | Simulation outputs and run-specific prioritization/comparison results. |
| `R/` | Stage-specific and shared R helpers. |
| `src/` | Rcpp simulation implementation. |
| `Figures/` | Code-generated manuscript figures. |
| `*_ppm_bin/`, `*_rangebag_bin/` | Binary species distribution inputs separated by taxon and SDM source. |
| `Main_Text.tex`, `Supplementary_Methods.tex`, `Title_Page.tex`, `references.bib` | Manuscript and bibliography source files. |

Large raster directories and generated result trees may not be practical to store in ordinary Git history. Restore all required local inputs before beginning a full run.

### Expected Data Objects

This section lists the main data objects expected by the scripts or produced by
them. It is intended as a practical editing map: when changing a filename,
directory, or object schema, search for the object here first and then follow
the stage-specific description later in the README.

#### External inputs

These objects are inputs to the repository rather than products of the numbered
pipeline:

```text
Data/Raw/
|-- mammal_rmax.txt
|-- mammal_data.txt
|-- bird_data.txt
|-- sigma.csv
|-- bird_synonyms.csv
|-- simple_summary.csv
|-- synonyms.csv
|-- random_effects.csv
`-- esacci_2022_pfts.tif

bird_ppm_bin/
bird_rangebag_bin/
mammal_ppm_bin/
mammal_rangebag_bin/
```

The four species-distribution directories contain binary species rasters used
by script 5. Their names encode taxon and SDM source, and script 5 selects among
them using its visible YAML parameters.

#### Generated clean objects

The workflow uses these reusable cleaned or staged objects in `Data/Clean/`:

```text
Data/Clean/
|-- post_mammal_rm_coefs.csv
|-- post_mammal_sigma_coefs.csv
|-- post_bird_r_coefs.csv
|-- post_bird_sigma_ilr_coefs.csv
|-- bird_ilr_combo_coords.csv
|-- bird_mass_grid.csv
|-- mammal_mass_grid.csv
|-- gompertz_loess_models.rds
|-- species_table.csv
|-- Patches/
|   |-- <species>.tif
|   `-- ...
|-- Patches_binary/
|   |-- <species>.tif
|   `-- ...
|-- all_patch_lookup.rds
|-- all_connectivity.rds
|-- PriorityInputs/
|   `-- pruning_inputs_curve_<curve>_taxa_<taxa>_sdm_<sdm>.rds
`-- ZonationOutputs/
    |-- zonation_binary_patch_feature_list.txt
    |-- settings.z5.txt
    |-- ABF/
    |-- CAZ1/
    |-- CAZ2/
    |-- CAZMAX/
    `-- benchmark_rank_maps/
        `-- rankmap_<rank_method>.tif
```

`species_table.csv` intentionally keeps the source-facing threshold column
names `min_patch_size` and `min_pop_size`. When script 6 constructs
`species_params`, those columns are renamed internally to
`min_patch_area_km2` and `min_population_area_km2` so later code is explicit
that the values are areas.

#### Generated result objects

Script 2 writes persistence simulations:

```text
Data/Results/
|-- persist_points_mammals.csv
`-- persist_points_birds.csv
```

Script 6 writes one run directory per selected curve, taxon set, SDM source,
and pruning schedule:

```text
Data/Results/PriorityRuns/<run_id>/
|-- removal_order.tif
|-- rankmap.tif
|-- removal_events.csv
`-- patch_lookup_tables/
    |-- stage_patch_lookup_stage_0001.csv
    |-- stage_patch_lookup_stage_0002.csv
    `-- ...
```

Stage 7.1 adds stage metadata to that run directory and validates the existing
stage-6 `rankmap.tif` as the canonical priority surface:

```text
Data/Results/PriorityRuns/<run_id>/ana/
`-- stage_meta.csv
```

Stage 7.2 expects run-local benchmark rank maps for each benchmark method that
will be compared:

```text
Data/Results/PriorityRuns/<run_id>/benchmark_rank_maps/
|-- rankmap_abf.tif
|-- rankmap_caz1.tif
|-- rankmap_caz2.tif
`-- rankmap_cazmax.tif
```

For each selected benchmark method, stage 7.2 writes:

```text
Data/Results/PriorityRuns/<run_id>/ana/rank_lut_<rank_method>/
|-- rank_stage_targets_<rank_method>.csv
|-- lut_stage_0000.rds
|-- lut_stage_0001.rds
`-- ...
```

Stage 7.3 writes paired persistence-comparison tables for each benchmark
method:

```text
Data/Results/PriorityRuns/<run_id>/ana/persist_cmp_<rank_method>/
|-- pu_long_<rank_method>.csv
|-- sp_long_<rank_method>.csv
|-- stage_sum_<rank_method>.csv
`-- stage_compare_<rank_method>.csv
```

Manuscript figure scripts write PNGs to `Figures/`. The manuscript source files
reference these generated assets but are not automatically synchronized with
them.

## Software Requirements

The workflow is written in R and R Markdown. Package requirements vary by stage.

### Rendering

- `knitr`
- `rmarkdown`

### General data handling

- `data.table`
- `readr`
- `dplyr`
- `tibble`
- `stringr`
- `tools`
- `scales`

### Statistical modeling and compositional data

- `rjags`
- `coda`
- JAGS installed on the operating system
- `compositions`
- `zCompositions`

### Simulation

- `Rcpp`
- `fastmatch`
- `Rfast`
- A working C++ toolchain compatible with the installed R version

### Spatial processing

- `terra`
- `sf`
- `igraph`
- `units`
- `fasterRaster`
- GRASS GIS for faster `fasterRaster` clumping

### Plotting

- `ggplot2`
- `cowplot`
- `ggtern`
- `viridisLite`
- `tidyterra`
- `rnaturalearth`

### IUCN queries

- `rredlist`
- An IUCN Red List API key in `IUCN_REDLIST_KEY`

Set the key before running script 4:

```r
Sys.setenv(IUCN_REDLIST_KEY = "your_key_here")
```

Scripts 5, 5.1, and 7.2 attempt to use `fasterRaster` for connected-patch
clumping. If `fasterRaster` or GRASS cannot be initialized, shared utilities
fall back to `terra::patches()` with a warning. For production runs, provide a
GRASS GIS installation directory through the Rmd `grass_dir` parameter or leave
it blank to auto-detect GRASS from common install locations, from an activated
conda environment, from a `grass` executable on `PATH`, or from one of these
environment variables:

- `FASTER_RASTER_GRASS_DIR`
- `GRASS_DIR`
- `GISBASE`

On a server, install GRASS in the conda environment that will render the Rmds:

```bash
conda install -c conda-forge grass
```

Then activate the environment before launching R or rendering scripts:

```bash
conda activate <env-name>
```

If auto-detection fails, set `FASTER_RASTER_GRASS_DIR` explicitly. In many
conda installations this directory is under `$CONDA_PREFIX/lib`, for example:

```bash
export FASTER_RASTER_GRASS_DIR="$CONDA_PREFIX/lib/grass84"
```

You can also pass the directory directly when rendering, for example:

```r
rmarkdown::render(
  "5_build_patches_and_connectivity.Rmd",
  params = list(grass_dir = Sys.getenv("FASTER_RASTER_GRASS_DIR"))
)
```

## Pipeline Overview

Run the canonical workflow in this order:

1. `1_rm_sigma_models.Rmd`
2. `2_persist_points_crn.Rmd`
3. `3_gompertz_loess_models.Rmd`
4. `4_build_species_table.Rmd`
5. `5_build_patches_and_connectivity.Rmd`
6. `5.2_build_binary_patch_rasters_zonation_feature_list.Rmd`
7. `6_spatial_prioritization_pipeline.Rmd`
8. `7.1_stage_meta_and_priority_surface.Rmd`
9. `7.2_build_rank_lut.Rmd`
10. `7.3_persist_cmp.Rmd`

`5.1_single_species_aoh_patches_pu_process.Rmd` is a figure-only branch. It uses the stage-4 species table and stage-5 helper logic but does not produce canonical stage-5 inputs.

Stage 5.2 is required only for benchmark analyses that use Zonation rank maps.
Run it after stage 5 to create binary patch rasters and the Zonation feature
list. After Zonation has been run, place the completed benchmark rank maps in
the stage-6 run's `benchmark_rank_maps/` directory before stages 7.2 and 7.3.

### Data flow

```text
Raw demographic data
  -> Stage 1 posterior demographic coefficients and mass/diet grids
  -> Stage 2 stochastic persistence points
  -> Stage 3 shifted-Gompertz and LOESS models

Species metadata + traits + SDM rasters + IUCN habitats
  -> Stage 4 canonical species table
  -> Stage 5 patches and connected population units
  -> Stage 5.2 binary patch rasters and Zonation benchmark inputs

Stage 3 models + Stage 4 species parameters + Stage 5 spatial objects
  -> Stage 6 reverse-removal prioritization
  -> Stage 7.1 stage metadata and priority-surface validation

Run-local benchmark rank maps + Stage 7.1 targets + Stage 5 baseline patches
  -> Stage 7.2 benchmark stage lookup tables

Stage 6 stage lookups + Stage 7.2 benchmark lookups
  -> Stage 7.3 persistence comparison, statistics, and figures
```

### What each stage contributes

The pipeline is intentionally sequential, but each stage has a distinct
contract. A useful way to decide what to rerun is to ask which contract has changed:

| Stage | Main question answered | Canonical handoff |
|---:|---|---|
| 1 | What demographic allometries and uncertainty distributions describe mammals and birds? | Posterior coefficient CSVs and mass/diet grids in `Data/Clean/`. |
| 2 | For each body mass or bird diet/mass combination, what abundance is needed to meet fixed persistence probabilities? | `persist_points_mammals.csv` and `persist_points_birds.csv`. |
| 3 | How do shifted-Gompertz parameters vary with body mass and diet? | `gompertz_loess_models.rds`. |
| 4 | Which species are analyzed, and what are their density, dispersal, habitat, and persistence parameters? | `species_table.csv`. |
| 5 | Where are each species' habitat patches and dispersal-connected population units? | Patch rasters, `all_patch_lookup.rds`, and `all_connectivity.rds`. |
| 5.2 | What binary patch rasters and feature lists are needed for Zonation benchmarks? | `Patches_binary/`, `ZonationOutputs/`, and copied `rankmap_<method>.tif` files when available. |
| 6 | In what order should cells be removed under the persistence-based objective? | A run directory containing rank rasters, events, and stage patch lookups. |
| 7.1 | How should the completed run be summarized stage by stage? | `stage_meta.csv`; `rankmap.tif` is validated as the canonical priority surface. |
| 7.2 | What patch/PU landscapes would the benchmark rank maps retain at the same stage targets? | `rank_lut_<method>/lut_stage_####.rds`. |
| 7.3 | How do persistence trajectories compare between the pipeline and each benchmark? | Long-form comparison CSVs, statistics, report text, and the comparison figure. |

### Regeneration guide

Do not delete artifacts automatically. Instead, regenerate only the downstream
products whose inputs changed:

| If this changes | Regenerate |
|---|---|
| Raw demographic data, JAGS settings, or stage-1 posterior files | Stages 1-4, then all spatial/downstream stages that depend on the species table. |
| Stage-2 simulation settings or persistence-point CSVs | Stages 2-4, then spatial/downstream stages if `species_table.csv` changes. |
| Gompertz fitting settings or `gompertz_loess_models.rds` | Stages 3-4, then spatial/downstream stages if `species_table.csv` changes. |
| Species metadata, IUCN habitat interpretation, density, dispersal, or threshold columns | Stages 4-7. |
| SDM rasters, habitat rasters, GRASS/fasterRaster behavior, or patch construction code | Stages 5-7, including stage 5.2 if benchmark rank maps will be rebuilt. |
| Script-6 run settings such as curve, taxa, SDM source, cells removed per iteration, or iterations per stage | Rebuild the script-6 bundle if selection changed; rerun stage 6 and all stage-7 analyses. |
| Benchmark rank maps | Rerun stages 7.2 and 7.3 for the affected benchmark methods. |

When an upstream contract or spatial domain changes, treat existing downstream
files as stale until the affected stages have been regenerated under the same
settings.

## Shared Conventions

### Quasi-extinction threshold

`R/analysis_contract.R` is the single source of truth for fixed scientific
constants and cross-stage vocabularies. `quasi_extinction_abundance()` returns
500 individuals, `minimum_patch_abundance()` returns 10 individuals, and
`persistence_horizon_years()` returns 100 years.

Stage 2 uses the quasi-extinction abundance for both the simulator threshold
and shifted-Gompertz grid construction. Stage 3 obtains `K0` from the same
helper rather than exposing a separate parameter. Stage 4 converts the two
abundance thresholds to area:

```text
min_patch_size = minimum_patch_abundance() / density
min_pop_size   = quasi_extinction_abundance() / density
```

Stages 6 and 7 validate these derived areas against density with floating-point
tolerance. The source CSV retains `min_patch_size` and `min_pop_size`; the
script-6 bundle uses `min_patch_area_km2` and `min_population_area_km2`.

### Strict area-threshold rule

Every spatial stage uses `area_exceeds_threshold()`. An area survives only when
it is finite and strictly greater than its threshold. Areas equal to, below, or
non-finite relative to the threshold fail. This applies to initial patch and
PU construction, process-figure reconstruction, pruning, fragmentation,
distance repair, and benchmark LUT reconstruction.

### Curve labels

Three demographic uncertainty curves are propagated:

- `q50`: median persistence behavior.
- `q16`: lower-tail behavior at the 16th percentile.
- `q025`: lower-tail behavior at the 2.5th percentile.

The selected stage-6 `curve_label` determines which fitted `alpha_<curve>` and `beta_<curve>` columns become `a_pred` and `b_pred` in the priority-input bundle.

`persistence_quantiles()` defines the named probabilities once,
`persistence_curves()` returns their ordered names, and
`validate_persistence_curve()` is used at configuration boundaries.

### Taxon and SDM tags

Canonical selection tags are:

- Taxa: `mammals`, `birds`, or `mammals_birds`.
- SDM sources: `ppm`, `rangebag`, or `ppm_rangebag`.

`priority_flags_from_tags()` converts tags back to logical gates for stages 7.1-7.3. `make_selection_tag()` creates stable underscore-separated tags during stage 6.

### Species identifiers

Scientific names remain the authoritative biological identifiers. `species_id()` converts a scientific name to a filename/table-safe ID by replacing spaces with underscores. Functions that read stage lookup tables accept scientific names where possible and reconstruct this ID consistently.

### Raster-cell domain

The stage-6 initial analysis domain is:

```r
alive_species_count_by_cell > 0L
```

Every initially alive cell must receive a removal order and rank value. Cells outside this domain remain `NA`. Benchmark rank maps used in stage 7.2 must have a finite value for every cell in this same domain and no geometric mismatch with `removal_order.tif`.

### Patch and population-unit definitions

- A patch is a rook-contiguous group of mapped-habitat cells.
- Patches at or below `min_patch_size` are removed.
- Retained patches within `dispersal_dist` are connected.
- Connected components form candidate population units.
- Population units at or below `min_pop_size` are removed.

The stage-6 algorithm updates these structures after cell loss. Stage 7.2 reconstructs the same concepts independently for each retained benchmark landscape. In both cases, equality does not survive.

## Stage 1: Demographic Calibration

### Driver and helpers

Driver:

- `1_rm_sigma_models.Rmd`

Helpers:

- `R/project_utils.R`
- `R/figure_utils.R`
- `R/rm_sigma_config.R`
- `R/compositions.R`
- `R/rm_sigma_data.R`
- `R/rm_sigma_models.R`
- `R/rm_sigma_figures.R`

### Purpose

Stage 1 calibrates body-mass and diet relationships for the demographic parameters used by the stochastic simulator:

- Maximum intrinsic growth rate for mammals.
- Environmental variability for mammals.
- Maximum intrinsic growth rate for birds.
- Environmental variability for birds, including diet composition.

It writes posterior coefficient draws, deterministic mass grids, and the canonical 66-combination bird diet grid.

### Parameters

| Parameter | Rmd default | Meaning |
|---|---:|---|
| `run_fits` | `false` | Refit JAGS models when `true`; validate and reuse posterior CSVs when `false`. |
| `seed` | `123` | Base seed for model initialization and reproducibility. |
| `n_chains` | `3` | Number of JAGS chains. |
| `n_adapt` | `5000` | Adaptation iterations per chain. |
| `n_iter` | `5000` | Retained sampling iterations before thinning. |
| `thin` | `10` | Thinning interval. |
| `tau_beta` | `0.001` | Regression-coefficient prior precision passed to JAGS. |
| `sigma_min` | `0` | Lower bound for residual sigma priors. |
| `sigma_max` | `1` | Upper bound for residual sigma priors. |

With the defaults, the expected posterior row count is:

```text
3 chains * 5000 iterations / 10 = 1500 draws
```

### Inputs

- `Data/Raw/mammal_rmax.txt`
- `Data/Raw/bird_data.txt`
- `Data/Raw/sigma.csv`
- `Data/Raw/bird_synonyms.csv`

### Chronological execution

#### 1. Configuration and preflight checks

The setup chunk sources `project_utils.R`, `figure_utils.R`, and `rm_sigma_config.R`.

The Rmd YAML block defines the public parameters and every stage-1 input/output path. `bayes_config()` converts the YAML values into the JAGS configuration, including prior precision, sigma bounds, chain settings, and seed.

`required_rm_sigma_packages()` adds `rjags` and `coda` only when `run_fits = true`. `load_packages()` fails before analysis if any required package is unavailable.

`assert_raw_inputs_exist()` verifies the four raw files. `assert_jags_available()` checks the system JAGS installation when fitting is requested.

When `run_fits = false`, `assert_posterior_inputs_exist()` requires all four posterior files and `assert_posterior_schema()` checks:

- Required coefficient columns.
- Exactly the expected number of posterior rows.

This means the default render is a validation-and-reuse run, not a silent partial fit.

#### 2. Reading and validating calibration data

`read_calibration_inputs()` reads the four source files and calls `validate_raw_inputs()`. The validation checks required columns and enforces a unique synonym mapping.

`prepare_calibration_data()` then constructs the model datasets:

1. `prepare_bird_traits()` parses bird mass and diet values, requires a high successful mass-parse rate, and rejects percentages outside 0-100.
2. Mammal growth data exclude bats and a fixed list of marine mammals before back-transforming log mass and log growth.
3. Environmental variability is computed as `sqrt(Vr)`.
4. A separate marine-mammal exclusion list is applied to mammal sigma data.
5. Bird growth values are derived as `log(lambda)` and joined to bird body masses.
6. `build_bird_sigma_join()` resolves sigma names through the bird synonym table and joins the trait data.
7. Bird diet is grouped into invertebrates, all plants, and vertebrate/fish/scavenging components.
8. `validate_calibration_data()` enforces the expected calibration sample sizes and requires every grouped bird diet composition to sum to 100.

The expected sample sizes encoded in the helper are:

| Dataset | Expected rows |
|---|---:|
| Mammal growth | 265 |
| Mammal sigma | 148 |
| Bird growth | 13 |
| Bird sigma | 225 |

#### 3. Constructing bird diet and mass grids

`ternary_diet_grid_10pct()` enumerates all non-negative three-part diets in 10% increments that sum to 100. This produces exactly 66 combinations.

`bird_ilr_combo_design()` calls `to_ilr3()` to:

1. Normalize the components.
2. Replace structural zeros with `zero_replace()` using the count-zero multiplicative method.
3. Apply the isometric log-ratio transform.
4. Validate that all `ilr1` and `ilr2` coordinates are finite.

The Rmd writes:

- `Data/Clean/bird_ilr_combo_coords.csv`
- A 31-point bird mass grid from 2 g to 11,236 g.
- A 31-point mammal mass grid from 2 g to 4,750,000 g.

Both mass grids are generated by `log_space()`.

#### 4. Fitting the demographic models

When `run_fits = true`, `fit_all_demographic_models()` fits four regressions.

`fit_loglog_allometry()` wraps `fit_jags_lm()` for:

- Mammal `rm`.
- Mammal `sigma`.
- Bird `r`.

Each response is modeled on the `log10` scale against `log10(body mass)`.

For bird sigma, `fit_all_demographic_models()` computes three-part ILR coordinates and calls `fit_jags_lm()` with predictors:

```text
log10(body mass), ilr1, ilr2
```

`fit_jags_lm()`:

1. Builds JAGS data and chain-specific initial values.
2. Fits the Gaussian regression in `jags_lm_model`.
3. Samples `alpha`, each beta coefficient, and residual `sigma`.
4. Converts the combined MCMC object to a table.
5. Renames beta columns to stable semantic names.
6. Writes the posterior CSV.

After fitting, `assert_posterior_schema()` is run again.

#### 5. Figure construction

The remaining chunks use the in-memory calibration objects and posterior files to build:

- The bird body-mass and diet sigma figure.
- The 2-by-2 allometric calibration figure.
- The bird sigma candidate-model comparison.

`save_manuscript_figure()` maps a stable figure ID to its canonical filename and writes a white-background PNG.

### Outputs

| Output | Description |
|---|---|
| `Data/Clean/post_mammal_rm_coefs.csv` | Mammal growth posterior draws. |
| `Data/Clean/post_mammal_sigma_coefs.csv` | Mammal environmental-variability posterior draws. |
| `Data/Clean/post_bird_r_coefs.csv` | Bird growth posterior draws. |
| `Data/Clean/post_bird_sigma_ilr_coefs.csv` | Bird sigma posterior draws with body mass and ILR diet effects. |
| `Data/Clean/bird_ilr_combo_coords.csv` | Canonical 66-combination diet design. |
| `Data/Clean/bird_mass_grid.csv` | 31-point bird body-mass grid. |
| `Data/Clean/mammal_mass_grid.csv` | 31-point mammal body-mass grid. |

Stable posterior schemas:

```text
post_mammal_rm_coefs.csv:
alpha, beta_logM, sigma

post_mammal_sigma_coefs.csv:
alpha, beta_logM, sigma

post_bird_r_coefs.csv:
alpha, beta_logM, sigma

post_bird_sigma_ilr_coefs.csv:
alpha, beta_logM, beta_ilr1, beta_ilr2, sigma
```

## Stage 2: Persistence-Point Simulation

### Driver and helpers

Driver:

- `2_persist_points_crn.Rmd`

Helpers:

- `R/persist_config.R`
- `R/persist_inputs.R`
- `R/persist_simulation.R`
- `src/simulate_persist_probs_cpp.cpp`

### Purpose

Stage 2 samples demographic parameter uncertainty from stage-1 posterior draws and evaluates stochastic population persistence across body masses, bird diet combinations, and equilibrium abundance values. It uses common random numbers (CRN) so comparisons across `K`, mass, and uncertainty curves are not dominated by unrelated Monte Carlo noise.

### Parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `run_mammals` | `true` | Generate mammal persistence points. |
| `run_birds` | `false` | Generate bird persistence points. |
| `force_rerun_bird_combos` | `false` | Recompute complete bird combinations. |
| `rebuild_cpp` | `true` | Force recompilation in a fresh cache directory. |
| `verbose` | `false` | Print detailed simulation progress. |
| `run_crn_self_test` | `true` | Verify exact CRN reproducibility and threshold behavior before simulation. |
| `cap_factor` | `1.1` | Carrying-capacity multiplier used when constructing simulation abundances. |
| `r_buffer` | `0.8` | Multiplicative buffer applied to growth-rate draws before simulation. |
| `n_draws` | `760` | Posterior draw indices sampled for each scenario. |
| `reps` | `2500` | Monte Carlo replicates per posterior draw. |
| `chunk_size` | `25` | Number of posterior draws processed per C++ chunk. |
| `base_seed` | `123` | Seed used to create the CRN context and posterior draw samples. |
| `p_anchor` | `[0.025, 0.16, 0.50, 0.84, 0.975]` | Anchor probabilities used while fitting the inversion curve. |
| `p_grid_lower` | `0.025` | First probability evaluated on the regular inversion grid. |
| `p_grid_step` | `0.025` | Step size for the regular inversion grid. |
| `p_grid_upper` | `0.975` | Last probability evaluated on the regular inversion grid. |
| `p_grid_extra` | `0.99` | Extra high-persistence target appended to the regular grid. |
| `k_search_start_k` | `512` | Starting abundance for inversion searches. |
| `k_search_rel_tol` | `0.01` | Relative tolerance for abundance inversion. |
| `k_search_round_k` | `true` | Round inverted abundances to integer `K` values. |

### Chronological execution

#### 1. Configuration

The Rmd YAML block defines the public run settings, simulation dimensions, CRN settings, target probability grids, stage-1 input files, and stage-2 output files. The setup chunk converts those visible values into `SIM`, `GRID`, and `PATHS` lists.

`validate_persist_config()` requires:

- At least one taxon.
- Positive simulation dimensions.
- A starting `K` above the threshold.
- Strictly increasing probability targets in `(0, 1)`.

The quasi-extinction threshold and time horizon are obtained from
`analysis_contract.R`; the Rmd contains no independent shift or horizon
setting.

#### 2. Compile and verify the C++ simulator

`load_persist_inputs()` begins with `check_persist_files()` and `compile_persist_cpp()`.

`compile_persist_cpp()` calls `Rcpp::sourceCpp()` and verifies four exported symbols:

- `crn_context_create`
- `simulate_persist_qvec_cpp`
- `simulate_persist_probs_cpp`
- `persist_cpp_contract`

The contract string must equal `initial-threshold-v2`. This check catches an old compiled object even if the cache exists.

#### 3. Load demographic inputs

`load_mammal_demographic_inputs()` reads the mammal mass grid and the `alpha`/`beta_logM` columns from the mammal growth and sigma posteriors.

`load_bird_demographic_inputs()` additionally reads `combo`, `ilr1`, and `ilr2`. `validate_bird_ilr_grid()` requires:

- Exactly 66 rows.
- Unique combo labels.
- Three integer components separated by underscores.
- Non-negative 10% increments.
- A row sum of 100.

`sample_posterior_indices()` samples growth and sigma posterior rows independently with deterministic seed offsets. The sampled indices are reused across scenarios.

#### 4. Create the CRN context

`create_crn_context()` allocates the C++ random-number context for all posterior draws, replicates, years, and chunks.

When enabled, `run_crn_self_test()`:

1. Creates a tiny context.
2. Calls the simulator twice with identical inputs and requires bit-for-bit identical results.
3. Checks that `K <= 500` always produces persistence zero.

#### 5. Construct demographic samplers

`make_mammal_sampler()` predicts sampled growth and sigma values from body mass.

`make_bird_sampler()` does the same but adds the selected diet combination's `ilr1` and `ilr2` effects to bird sigma.

Growth is multiplied by the fixed `r_buffer` before simulation.

#### 6. Generate persistence points

`run_one_scenario()` loops over body masses and the three uncertainty curves.

For each mass:

1. The taxon-specific sampler produces vectors of `r` and `sigma`.
2. `eval_qvec_crn()` evaluates persistence quantiles through the C++ simulator.
3. Results are cached by `K` so repeated searches do not repeat simulations.
4. `find_K_for_target()` uses bracket expansion and bisection to locate `K` values for five anchor probabilities.
5. `fit_gompertz_ab()` fits a temporary shifted-Gompertz approximation to those anchors.
6. `invert_gompertz()` converts the full probability grid to proposed `K` values.
7. The simulator evaluates those proposed values directly.
8. Rows are appended to the output CSV.

The temporary Gompertz fit is a grid-construction device. Stage 3 subsequently refits the final curves to the generated persistence points.

#### 7. Mammal output behavior

The mammal chunk calls `run_one_scenario()` once with `combo = NA`. The output is replaced because `append = false`.

`validate_points_schema()` then requires the standard columns and verifies that every `K` is finite and strictly greater than 500.

#### 8. Bird resume behavior

Bird simulations are managed by `run_bird_combo_set()`:

1. `expected_rows_per_combo()` calculates the required row count.
2. `remove_partial_combos()` removes incomplete combinations from an existing output.
3. `completed_combos()` identifies combinations safe to skip.
4. Each requested combination is written first to a temporary CSV.
5. The temporary output is appended only after its row count is validated.
6. `finalize_bird_output()` filters to the canonical 66 combinations, sorts the file, and rejects missing or incomplete combinations.

This design allows an interrupted bird run to resume without treating partial output as complete.

### Outputs

- `Data/Results/persist_points_mammals.csv`
- `Data/Results/persist_points_birds.csv`
- `Data/.rcpp_cache/`

Schema:

```text
combo, mass_idx, mass_g, curve, K, p_eval
```

Under the current grids:

| Output | Expected size |
|---|---:|
| Mammals | 3,720 rows |
| Each bird combo | 3,720 rows |
| All 66 bird combos | 245,520 rows |

## Stage 3: Gompertz and LOESS Models

### Driver and helpers

Driver:

- `3_gompertz_loess_models.Rmd`

Helpers:

- `R/project_utils.R`
- `R/figure_utils.R`
- `R/gompertz_config.R`
- `R/gompertz_core.R`
- `R/gompertz_figures.R`

### Purpose

Stage 3 converts the simulated persistence points into compact models that can be evaluated for any retained species:

1. A shifted two-parameter Gompertz curve is fitted for every mass, curve, and bird diet combination.
2. LOESS models interpolate positive Gompertz parameters over body mass.
3. The fitted LOESS objects are saved for stage 4.

The persistence function is:

```text
P(K) = exp[-alpha * (K - K0)^(-beta)]
```

with `K0 = 500`.

### Parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `run_mammals` | `true` | Fit mammal curves. |
| `run_birds` | `true` | Fit bird curves. |
| `write_models` | `true` | Save the model RDS. |
| `loess_span` | `0.75` | LOESS smoothing span. |
| `loess_z` | `1.96` | Multiplier for displayed prediction bands. |

### Chronological execution

#### 1. Setup and validation

`gompertz_curves()` delegates to `persistence_curves()` in
`analysis_contract.R`. `K0` is obtained from
`quasi_extinction_abundance()` and passed explicitly to the mathematical
Gompertz helpers.

`validate_gompertz_config()` requires:

- At least one taxon.
- A LOESS span in `(0, 1]`.
- A positive `loess_z`.
- The canonical three curve labels.

`check_gompertz_inputs()` validates the persistence-point schemas and, for birds, the diet-coordinate and sigma-posterior files used in figures.

#### 2. Read persistence points

`load_all_persistence_points()` calls `read_persistence_points()` for each requested taxon. It standardizes types and discards rows outside the selected curves, invalid masses, and `K <= 500`.

`validate_bird_combo_coverage()` compares the bird combos in the persistence table against `bird_ilr_combo_coords.csv`. Missing or extra combinations are fatal.

#### 3. Fit shifted-Gompertz curves

`fit_all_gompertz_parameters()` groups by:

```text
group, combo, mass_idx, mass_g, curve
```

and calls `fit_gompertz_ab()` for each group.

`fit_gompertz_ab()`:

1. Clamps probabilities away from exact zero and one.
2. Uses a linearized relationship to obtain initial `alpha` and `beta`.
3. Tries multiple positive beta starting values.
4. Fits the nonlinear model on log-parameter scales.
5. Back-transforms the fitted parameters.
6. Calculates a linearized RMSE.
7. Returns `fit_ok = false` and `NA` parameters instead of aborting the entire run when one fit fails.

#### 4. Fit body-mass interpolation models

`build_gompertz_loess_models()` separates successful mammal and bird fits.

`fit_curve_loess()` calls `fit_loess_logparam()` separately for `alpha` and `beta`. Both axes are transformed:

```text
x = log10(body mass)
y = log(parameter)
```

Mammals receive one LOESS set per uncertainty curve. Birds receive a separate LOESS set for every diet combo and uncertainty curve.

`predict_loess_exp()` and `predict_loess_set()` produce a dense prediction grid for plotting, including exponentiated mid, lower, and upper values.

The saved model object includes metadata describing:

- Creation time.
- `K0`.
- Curve labels.
- LOESS span.
- Predictor and response transforms.
- Source persistence-point files.

#### 5. Figures

The Rmd builds:

- A shifted-Gompertz versus Wolff comparison for representative mammal masses.
- Mammal LOESS `alpha` and `beta` curves.
- Bird LOESS parameter curves colored by the diet-dependent sigma effect.

The bird figure reads the stage-1 ILR design and sigma posterior to reproduce the same color mapping used in the bird sigma figure.

### Output

- `Data/Clean/gompertz_loess_models.rds`

Top-level structure:

```text
meta
mammals
  q50
  q16
  q025
birds
  <diet_combo>
    q50
    q16
    q025
```

Each curve contains separate positive-scale LOESS fits for `alpha` and `beta`.

## Stage 4: Canonical Species Table

### Driver and helpers

Driver:

- `4_build_species_table.Rmd`

Helpers:

- `R/project_utils.R`
- `R/species_table_config.R`
- `R/species_table_inputs.R`
- `R/species_table_traits.R`
- `R/species_table_iucn.R`
- `R/species_table_build.R`

### Purpose

Stage 4 creates the canonical row-per-species input for all spatial stages. It resolves species names across metadata, traits, and raster filenames; computes density and movement-derived thresholds; predicts Gompertz parameters; and queries suitable IUCN habitats.

### Parameters

| Parameter | Default | Meaning |
|---|---:|---|
| `write_output` | `true` | Write `species_table.csv`. |
| `gompertz_mammals` | `true` | Predict mammal Gompertz parameters. |
| `gompertz_birds` | `true` | Predict bird Gompertz parameters. |
| `iucn_pause_seconds` | `2` | Delay between IUCN species queries. |

### Inputs

- `Data/Raw/simple_summary.csv`
- `Data/Raw/synonyms.csv`
- `Data/Raw/mammal_data.txt`
- `Data/Raw/bird_data.txt`
- `Data/Raw/random_effects.csv`
- `Data/Clean/gompertz_loess_models.rds`
- `mammal_ppm_bin/`
- `mammal_rangebag_bin/`
- `bird_ppm_bin/`
- `bird_rangebag_bin/`
- `IUCN_REDLIST_KEY`

### Chronological execution

#### 1. Configuration validation

`validate_species_table_config()` checks:

- A non-negative API pause.
- A non-empty IUCN key.
- Every raw file and raster directory.
- Required columns in the summary, synonym, trait, and random-effect files.
- Presence of bird `Diet-*` columns.
- The Gompertz RDS when either taxon uses fitted models.

#### 2. Read core inputs

`read_species_core_inputs()`:

1. Reads the species summary and standardizes taxon class.
2. Calls `read_synonyms()`.
3. Calls `list_raster_manifest()` across all four raster directories.
4. Calls `read_traits()` for mammal and bird tables.
5. Calls `read_random_effects()`.
6. Loads the stage-3 model object when requested.

The raster manifest retains only files whose lowercase stem ends in `_bin`. PPM receives a lower numeric priority than RangeBag when both are available for the same candidate name.

#### 3. Resolve names, rasters, and traits

`build_name_candidates()` creates original-name and synonym candidates. Each candidate records:

- Candidate scientific name.
- Match type.
- Expected raster stem.
- Normalized trait key.
- Whether a matching trait record exists.

`pick_rasters()` joins candidates to the raster manifest and chooses one raster per source species using:

1. SDM priority.
2. Original name before synonym.
3. Candidate name as a stable tie breaker.

`pick_traits()` independently chooses the best trait match.

`resolve_species_inputs()` retains only rows with both a raster and trait match, joins the selected trait record, and preserves columns describing how each match was made.

#### 4. Derive diet, density, movement, and thresholds

`add_trait_covariates()` first uses `diet_aggregates()` to combine raw diet columns.

For mammals, `classify_mammal_diet()` assigns carnivore, herbivore, or omnivore. For birds, `get_bird_diet5()` uses an existing five-category field when available and otherwise derives a dominant category.

`make_bird_combo()` converts the three grouped bird diet percentages to the stage-1 combo label.

`attach_random_effects()` joins order, family, and species effects; missing effects are set to zero.

`compute_density_and_dispersal()` then calculates:

- Taxon-specific density from body mass, diet, and random effects.
- Taxon-specific home-range size.
- Dispersal distance as a function of home range.
- Minimum patch area as `minimum_patch_abundance() / density`.
- Minimum population-unit area as `quasi_extinction_abundance() / density`.

The helper return values are 10 and 500 individuals, respectively.

#### 5. Predict species-specific Gompertz parameters

`add_gompertz_parameters()` loops over the ordered values returned by
`persistence_curves()`.

For mammals, it evaluates the taxon-wide LOESS fits. For birds, it chooses the LOESS set matching `diet_combo`.

`predict_positive_loess()` clamps `log10(body mass)` to the training range before prediction and exponentiates the result. This avoids uncontrolled extrapolation and preserves positive `alpha` and `beta`.

#### 6. Query IUCN habitats

`query_species_habitats()` queries each unique genus/species pair sequentially and logs elapsed time and suitable codes.

`query_iucn_habitats()`:

1. Calls `rredlist::rl_species_latest()`.
2. Keeps habitat rows marked suitable.
3. Normalizes IUCN habitat codes.
4. Maps top-level codes to descriptive labels.
5. Uses `summarize_artificial_terrestrial()` to refine code 14 into arable/pastureland, plantations/degraded forest, and urban/rural gardens.

The output columns are:

- `habitat_codes_suitable`
- `habitats_level1`
- `habitats_mixed`

#### 7. Assemble and validate the final table

`build_species_table()` joins the habitat results, selects `species_table_keep_cols()`, and calls `validate_species_table()`.

Validation requires:

- At least one retained row.
- Positive body mass, density, home range, dispersal, patch threshold, and population threshold.
- A non-empty `habitats_mixed` field for every retained row.
- Positive finite Gompertz parameters when those models are enabled.

### Output

- `Data/Clean/species_table.csv`

Important downstream columns:

```text
scientificName
className
redlistCategory
sdm_method
raster_path
raster_dir
raster_file
BodyMass.Value
Diet
diet_combo
density
home_range_size
dispersal_dist
min_patch_size
min_pop_size
alpha_q50, beta_q50
alpha_q16, beta_q16
alpha_q025, beta_q025
habitat_codes_suitable
habitats_level1
habitats_mixed
```

## Stage 5: Habitat Patches and Connectivity

### Driver and helpers

Driver:

- `5_build_patches_and_connectivity.Rmd`

Helpers:

- `R/project_utils.R`
- `R/patch_config.R`
- `R/patch_species_inputs.R`
- `R/patch_habitat.R`
- `R/patch_processing.R`

### Purpose

Stage 5 converts each selected species' distribution and habitat associations into:

- A final patch-ID raster.
- A table linking patches to population units.
- A sparse graph describing patch connectivity within each population unit.

### Parameters

The Rmd exposes taxon gates, all input/output paths, GRASS location, ROI bounds, raster overwrite behavior, and whether a species error stops the full run.

Current ROI defaults cover Madagascar:

```text
x: 43.18 to 50.56
y: -25.64 to -11.89
```

### Chronological execution

#### 1. Configure packages, paths, and ROI

`load_patch_packages()` loads spatial dependencies and enables `sf` s2 geometry.

`patch_paths()`, `patch_run_options()`, and `patch_roi()` convert parameters into structured objects.

`validate_patch_config()` verifies the species table, land-cover raster, GRASS directory, output directories, and taxon gates. It initializes `fasterRaster` with the selected GRASS installation.

#### 2. Read selected species

`read_patch_species_table()` reads only the columns needed by the patch pipeline, parses thresholds, normalizes text, filters to the requested taxa, and calls `validate_selected_patch_species()`.

The validation checks:

- Non-empty names, habitat labels, and raster paths.
- Existence of every species raster.
- Positive patch/population thresholds.
- Non-negative dispersal.
- Unique output patch filenames.

#### 3. Build reusable land-cover masks

`load_habitat_masks()` crops the ESA CCI raster to the ROI and calls `build_habitat_masks_from_landcover()`.

`landcover_class_map()` defines the ESA class-code crosswalk for:

- Forest
- Savanna
- Shrubland
- Grassland
- Inland wetlands
- Rocky areas
- Desert
- Arable/pastureland
- Plantations/degraded forest
- Urban/rural gardens
- Artificial aquatic habitat
- Artificial terrestrial habitat

Each mask is represented as `1` for included cells and `NA` elsewhere.

#### 4. Process each species

`run_patch_pipeline()` loops over selected species and calls `process_patch_species()`. Errors are either recorded as status `error` or rethrown depending on `stop_on_species_error`.

For each species, `process_patch_species()` calls `mapped_habitat_for_species()`:

1. `parse_habitats_mixed()` separates the species' suitable labels.
2. `habitat_label_map()` maps labels to mask-layer names.
3. `habitat_union()` unions all suitable masks.
4. `load_presence_aligned01()` reads the distribution raster, projects/resamples it with nearest-neighbor methods when needed, and converts value 1 to presence.
5. The habitat and presence rasters are intersected to form mapped habitat.

If mapped habitat is non-empty, `patch_components_for_species()` performs the spatial construction:

1. Trim mapped habitat to its occupied extent.
2. Use `fasterRaster::clump(..., diagonal = FALSE)` to label rook-contiguous raw patches.
3. Calculate patch area with `terra::cellSize()` and `terra::zonal()`.
4. Remove patches at or below `min_patch_km2`.
5. Renumber retained patches.
6. Polygonize retained patches.
7. Link polygons within `disp_km` using `sf::st_is_within_distance()`.
8. Build an undirected patch graph with `igraph`.
9. Assign candidate PUs from graph components.
10. Sum patch area by candidate PU.
11. Remove PUs at or below `min_pop_km2`.
12. Renumber final patches and PUs.
13. Build a patch lookup table.
14. Call `build_pu_connectivity()` to encode each final PU as a sparse CSR graph.

`process_patch_species()` writes the final patch raster as an unsigned integer GeoTIFF and returns its tables and connectivity objects.

#### 5. Combine species outputs

After the loop, `run_patch_pipeline()`:

- Keeps results with status `retained`.
- Row-binds all patch lookup tables.
- Flattens all PU connectivity entries.
- Attaches CSR metadata attributes.
- Builds a status summary with `patch_status_summary()`.

The Rmd writes the two combined RDS objects.

### Outputs

- `Data/Clean/Patches/<species_id>.tif`
- `Data/Clean/all_patch_lookup.rds`
- `Data/Clean/all_connectivity.rds`

Patch lookup schema:

```text
scientificName, patch_id, pu_id, patch_area_km2
```

Connectivity entry:

```text
species
pu_id
patch_ids
row_ptr
col_idx
```

Connectivity attributes:

```text
csr_version = 2
row_ptr_base = "0-based"
col_idx_space = "pu_local_index"
col_idx_base = "0-based"
```

## Stage 5.1: Single-Species Process Figure

### Driver and helpers

Driver:

- `5.1_single_species_aoh_patches_pu_process.Rmd`

Additional helper:

- `R/patch_process_figure.R`

Stage 5.1 also sources the stage-5 configuration, species-input, habitat, and processing helpers.

### Purpose

This script reproduces the stage-5 logic for one species and exposes intermediate rasters for a six-panel process figure. It does not write canonical patch rasters, patch lookups, or connectivity objects.

### Chronological execution

1. `validate_patch_config(..., production = FALSE)` checks inputs and initializes GRASS without creating production outputs.
2. `read_single_patch_species()` finds and validates the requested species row.
3. `build_single_species_layers()` rebuilds habitat masks, distribution presence, and mapped habitat.
4. It clumps all raw patches and records their count.
5. It applies the minimum-patch threshold and builds a retained/discarded status raster.
6. It polygonizes retained patches and links them by dispersal distance.
7. It labels candidate PUs and applies the minimum-PU-area filter.
8. It creates candidate and final PU rasters, dissolves final PU outlines, and crops all layers to a shared plotting extent.
9. `make_single_species_process_figure()` builds panels for habitat, distribution, mapped habitat, patch filtering, candidate PUs, and final PUs.
10. `save_manuscript_figure()` writes the canonical figure.

### Output

- `Figures/spatial_process.png`

## Stage 5.2: Zonation Benchmark Inputs

### Driver

- `5.2_build_binary_patch_rasters_zonation_feature_list.Rmd`

### Purpose

Stage 5.2 prepares the stage-5 patch rasters for Zonation benchmark analyses.
It converts final patch-ID rasters into binary retained-habitat rasters, writes a
Zonation feature list with equal feature weights, creates expected Zonation
output directories, and copies completed Zonation rank maps into a shared
benchmark staging directory when they are available.

### Chronological execution

1. The driver reads final patch rasters from `Data/Clean/Patches/`.
2. Each non-`NA` patch cell is converted to value 1 and written to
   `Data/Clean/Patches_binary/` with the original raster filename.
3. `Data/Clean/ZonationOutputs/zonation_binary_patch_feature_list.txt` is
   written with one equal-weight row per binary raster.
4. Output directories are created for `ABF`, `CAZ1`, `CAZ2`, and `CAZMAX`.
5. If a completed Zonation `rankmap.tif` exists in one of those directories, it
   is copied to `Data/Clean/ZonationOutputs/benchmark_rank_maps/` as
   `rankmap_<rank_method>.tif`.

The files in `Data/Clean/ZonationOutputs/benchmark_rank_maps/` are a staging
copy. Stage 7.2 reads run-local benchmark maps from
`Data/Results/PriorityRuns/<run_id>/benchmark_rank_maps/`, so rank maps must be
copied there for each stage-6 run that will be compared.

### Outputs

```text
Data/Clean/
|-- Patches_binary/
|   |-- <species>.tif
|   `-- ...
`-- ZonationOutputs/
    |-- zonation_binary_patch_feature_list.txt
    |-- ABF/
    |-- CAZ1/
    |-- CAZ2/
    |-- CAZMAX/
    `-- benchmark_rank_maps/
        `-- rankmap_<rank_method>.tif
```

## Stage 6: Spatial Prioritization

### Driver and helper architecture

Driver:

- `6_spatial_prioritization_pipeline.Rmd`

Configuration and orchestration:

- `R/priority_run_config.R`
- `R/priority_inputs.R`
- `R/priority_outputs.R`
- `R/priority_pipeline.R`
- `R/priority_saved_run.R`

Pruning:

- `R/priority_pruning_frontier.R`
- `R/priority_pruning_scoring.R`
- `R/priority_pruning_graph.R`
- `R/priority_pruning_iteration.R`
- `R/priority_pruning_stage.R`

Fragmentation repair:

- `R/priority_fragmentation_patches.R`
- `R/priority_fragmentation_graph.R`
- `R/priority_fragmentation_stage.R`

Distance-connectivity repair:

- `R/priority_distance_geometry.R`
- `R/priority_distance_graph.R`
- `R/priority_distance_stage.R`

### Purpose

Stage 6 performs reverse-removal spatial prioritization. Low-scoring frontier cells are removed in batches. After pruning, the code repairs patch fragmentation and dispersal connectivity, removes patches or PUs that no longer satisfy biological thresholds, and records a complete removal order for the initial occupied-cell domain.

The driver supports two operations:

1. Build a reusable clean priority-input bundle.
2. Run the prioritization from a saved bundle.

These are controlled independently by `initialize_inputs` and `run_pipeline`.

### Run naming

`priority_run_id()` constructs:

```text
curve_<curve>_taxa_<taxa>_sdm_<sdm>_remove_<cells>_k_<iterations>
```

Example:

```text
curve_q025_taxa_mammals_birds_sdm_ppm_rangebag_remove_1000_k_100
```

`priority_bundle_name()` omits pruning-rate settings because the same initial bundle can support multiple removal schedules.

### Part A: Build the priority-input bundle

When `initialize_inputs = true`, the driver calls `initialize_priority_inputs()`.

#### A1. Validate selection gates

The function validates curve, taxon gates, and SDM gates. At least one taxon and one SDM source are required.

`normalize_taxon_class()`, `normalize_sdm_method()`, and `normalize_redlist_category()` standardize metadata labels. Red List category is retained as metadata but does not define a run branch or filename tag.

#### A2. Construct `species_params`

The stage-4 table is filtered to the selected taxon and SDM sources. The selected curve columns become:

```text
alpha_<curve> -> a_pred
beta_<curve>  -> b_pred
```

The table also carries density, dispersal distance, patch threshold, population threshold, taxon, SDM method, and Red List category. At this boundary the stage-4 columns are renamed:

```text
min_patch_size -> min_patch_area_km2
min_pop_size   -> min_population_area_km2
```

Positive finite values are required, and
`validate_density_area_thresholds()` verifies both areas against density and
the canonical abundance constants.

#### A3. Standardize the patch table

`all_patch_lookup.rds` is reduced to:

```text
species, patch_id, pu_id, patch_area_km2
```

Only selected species are retained.

#### A4. Require final spatial outputs

Species must have both:

- A patch raster in `Data/Clean/Patches/`.
- At least one patch lookup row.

The species parameter and patch tables are reduced to this final retained set.

#### A5. Build raster-index structures

All retained patch rasters are stacked on a common grid.

For each species:

- `patch_id_by_species_env` stores the full cell-to-patch vector.
- `patch_cell_index_by_species_env` stores a compact patch-sorted `pid`/`cell` index.

The code also derives:

- `cell_area_by_cell`
- Rook-neighbor cell pairs
- `alive_species_count_by_cell`, the number of species occupying each raster cell

#### A6. Convert connectivity graphs

The stage-5 connectivity object is validated against its CSR attributes.

Stored `col_idx` values are converted from zero-based to one-based R indexing. `csr_symmetrize_undirected()` removes loops, adds reverse edges, deduplicates edges, and rebuilds sorted CSR arrays. `check_csr_symmetric()` verifies the result.

Graphs are stored in `pu_graphs_by_key` under keys:

```text
<species>|<pu_id>
```

#### A7. Validate and save

Internal checks compare:

- Species across parameters, patch table, and raster stack.
- Patch IDs in rasters versus patch lookup rows.
- Graph patch IDs versus patch table.
- Alive-species counts versus non-`NA` raster layers.
- Graph symmetry.

Environments are converted to named lists before serialization. The saved bundle represents the state before any pruning call. Its metadata contains `schema_version = 2`, obtained from `priority_bundle_schema_version()`.

Schema version 2 is a deliberate clean break. Unversioned bundles and bundles
containing the former `min_patch_size500` field are rejected by
`validate_priority_bundle_schema()` with an instruction to regenerate them.

Bundle output:

```text
Data/Clean/PriorityInputs/
  pruning_inputs_curve_<curve>_taxa_<taxa>_sdm_<sdm>.rds
```

### Part B: Load and run a saved bundle

When `run_pipeline = true`, the driver calls `run_priority_pipeline_from_saved_bundle()`.

#### B1. Validate run settings and bundle identity

The wrapper reconstructs the taxon and SDM tags, locates the bundle, validates
its schema before its fields, and verifies that its metadata exactly matches
the requested curve and logical gates.

It restores the patch table, graph list, alive counts, areas, species parameters, and rook-neighbor pairs. The saved patch-index lists are reconstructed as live environments.

#### B2. Reconstruct raster templates

The first retained species patch raster becomes the mask template. All retained patch rasters are loaded into `species_raster_template_stack` for later distance-stage geometry.

The wrapper checks species names and raster dimensions before calling `run_priority_pipeline()`.

### Core stage loop in `run_priority_pipeline()`

#### 1. Initialize output and global state

The function creates the run directory and `patch_lookup_tables/`.

It defines:

- `initial_alive_by_cell`
- Initial alive-cell count and area
- `removal_order_by_cell`, initially `NA`
- A removal-step counter
- A list of removal-event rows
- Completed-stage count
- Frontier exhaustion state

The local `record_removed_cells()` helper assigns one removal step to newly removed initial-domain cells and records event type, stage, pruning iteration, cell count, and area.

`record_final_retained_cells()` later assigns all still-unordered initial cells to one synthetic terminal event.

#### 2. Initialize the frontier

`build_rook_neighbor_index()` converts rook-neighbor pairs into a compact symmetric neighbor index.

`initialize_frontier_state()` identifies currently alive cells adjacent to at least one empty cell. `get_frontier_cells()` retrieves the active frontier, and `update_frontier_state_after_empty_cells()` updates only affected neighborhoods after later removals.

#### 3. Pruning stage

Each global stage calls `run_pruning_stage()`, which repeatedly calls `run_pruning_iteration()` up to `pruning_iterations_per_stage` times or until the frontier is exhausted.

Within a pruning iteration:

1. `update_patch_score_cache()` computes or refreshes scores for dirty species/PUs.
2. `compute_pu_log_scores()` evaluates the persistence-based contribution of each PU.
3. Patch-level scores are cached by species and patch ID.
4. `score_frontier_cells()` combines the patch scores of species occupying each frontier cell.
5. `choose_frontier_cells_to_remove()` selects up to `cells_to_remove_per_iteration` cells.
6. Removed area is summarized with `summarize_removed_patch_area()`.
7. Species patch vectors and patch areas are updated.
8. `rebuild_pu_after_patch_loss()` uses the shared CSR helpers in `R/priority_csr_graph.R` to repair affected PU graphs after complete patch loss.
9. The patch table, graph list, alive-species counts, frontier state, and dirty sets are returned.

The outer pipeline records each pruning iteration as its own chronological removal event.

#### 4. Fragmentation stage

`run_fragmentation_stage()` examines patches changed by pruning.

Its supporting functions:

- `label_patch_cell_components()` identifies rook-connected fragments within a patch.
- `summarize_species_patch_areas()` computes fragment areas.
- `rebuild_species_patch_index()` rebuilds compact patch-to-cell indexing.
- `build_provisional_pu_graph()` integrates retained/reassigned patch nodes with an existing PU graph.
- `make_lazy_added_node_overlay_pu_graph()` avoids eagerly copying large graph structures when nodes are added.
- `rebuild_pu_after_patch_loss()` re-evaluates PU components and minimum-area viability.

Fragments and resulting PUs that fail thresholds are removed. The stage returns removed cells and patches requiring a distance recheck. The main loop updates the frontier and records one `fragmentation` event.

#### 5. Distance-connectivity stage

`run_distance_connectivity_stage()` verifies that graph edges still satisfy species dispersal distances after patch geometry changes.

The geometry helpers:

- `materialize_species_patch_raster()`
- `extract_candidate_patch_ids_for_distance()`
- `build_local_patch_polygons()`
- `build_distance_predicate_lookup()`
- `get_species_dispersal_threshold_km()`

The graph helpers, including shared CSR utilities from `R/priority_csr_graph.R`:

- `materialize_lazy_added_node_overlay_graph()`
- `filter_distance_invalid_edges_for_pu()`
- `find_csr_components()`
- `build_csr_subgraph()`
- `rebuild_pu_after_edge_filter()`

Invalid edges are removed, disconnected components are reassigned, and components at or below the population threshold are removed. The main loop updates the frontier and records one `distance` event.

#### 6. Complete the global stage

After pruning, fragmentation, and distance repair:

- `write_stage_patch_lookup_table()` writes the current patch table.
- `log_priority_pipeline_stage()` reports frontier status, removal steps, remaining patches, remaining PUs, and remaining alive cells.
- The completed-stage count advances.

The loop stops when the frontier is exhausted or `max_stages` is reached.

#### 7. Assign the terminal layer

Every initially alive cell must have a rank. `record_final_retained_cells()` assigns all cells still retained at the stopping condition to one synthetic final event at `completed_stages + 1`.

This synthetic event is not a stage patch lookup. It exists to complete the rank ordering.

#### 8. Write outputs

`write_removal_order_raster()` writes integer removal steps.

The event rows are combined and augmented with:

- Cumulative cells and area removed.
- Cumulative cell and area proportions.
- Cells and area retained.

`write_rankmap_raster()` maps each removal step to cumulative cell proportion removed. The terminal retained event maps to 1.

The script does not reread those just-written rasters as a separate
verification pass. Later stages check the files they need when they load them,
which keeps script 6 focused on the expensive prioritization work and its
canonical outputs.

### Outputs

```text
Data/Results/PriorityRuns/<run_id>/
|-- removal_order.tif
|-- rankmap.tif
|-- removal_events.csv
`-- patch_lookup_tables/
    |-- stage_patch_lookup_stage_0001.csv
    |-- stage_patch_lookup_stage_0002.csv
    `-- ...
```

## Stage 7.1: Stage Metadata and Priority Surface

### Driver and helper

- `7.1_stage_meta_and_priority_surface.Rmd`
- `R/stage7_stage_meta_surface.R`

### Purpose

Stage 7.1 converts event-level stage-6 output into one row per completed
ecological stage and validates the stage-6 `rankmap.tif` as the canonical
priority surface.

### Chronological execution

1. `priority_paths()` reconstructs the exact bundle and run directory from the curve, taxon, SDM, and pruning settings.
2. The helper loads the bundle and optionally validates its metadata against `priority_flags_from_tags()`.
3. The initial domain and cell areas are read from the bundle.
4. `removal_events.csv` is loaded and validated. Exactly one `final_retained` event is required.
5. The synthetic terminal event is excluded from ecological stage summaries.
6. Stage patch lookup filenames are parsed, required to be unique, and required to form a contiguous sequence beginning at stage 1.
7. An explicit stage-0 row is created from the initial bundle.
8. For each stage, the helper finds the last event at or before that stage and records cumulative retained/removed cells and area.
9. It also sums cells, area, and event count removed within the stage.
10. Cell-based and area-based removal proportions are stored separately.
11. The existing stage-6 `rankmap.tif` is loaded as the priority surface.
12. The rank map is checked against `removal_order.tif`, the initial alive-cell
    domain, and the expected `[0, 1]` value range.

### Outputs

- `Data/Results/PriorityRuns/<run_id>/ana/stage_meta.csv`

Important stage metadata:

```text
stage
patch_lookup_path
last_removal_step
alive_end
removed_end
area_retained_km2
area_removed_km2
cells_removed_stage
area_removed_stage_km2
events_in_stage
prop_cells_removed_end
prop_area_removed_end
pct_cells_removed_end
pct_area_removed_end
```

Canonical priority surface interpretation:

- Higher values indicate later removal and therefore higher priority.
- Terminal retained cells have value 1.
- Cells outside the initial domain are `NA`.
- Values are cell-order proportions, not area-weighted proportions.
- The priority surface path is `Data/Results/PriorityRuns/<run_id>/rankmap.tif`.

## Stage 7.2: Benchmark Rank Lookup Tables

### Driver and helper

- `7.2_build_rank_lut.Rmd`
- `R/stage7_rank_lut.R`

### Purpose

Stage 7.2 reconstructs the landscape retained by a benchmark rank map at every stage-6 cell-count target. It then rebuilds species patches and population units so benchmark and persistence-based landscapes can be evaluated with the same demographic persistence model.

Supported methods:

- `abf`
- `caz1`
- `caz2`
- `cazmax`

`benchmark_rank_methods()` and `benchmark_rank_labels()` define this registry
once. `validate_benchmark_rank_method()` normalizes and validates user input.
`benchmark_artifact_paths()` constructs the rank-map, LUT-directory,
and target-summary paths shared by stages 7.2 and 7.3.

Expected benchmark filename:

```text
Data/Results/PriorityRuns/<run_id>/
  benchmark_rank_maps/rankmap_<rank_method>.tif
```

Stage 5.2 can stage completed Zonation rank maps under
`Data/Clean/ZonationOutputs/benchmark_rank_maps/`, but stage 7.2 reads only the
run-local copy for the selected stage-6 run.

### Chronological execution

#### 1. Load and validate inputs

The helper reconstructs run paths, reads the stage-6 bundle, stage metadata, species table, and selected benchmark rank map.

The bundle supplies:

- Initial alive-cell domain.
- Cell areas.
- Retained species.
- Baseline patch table.
- Optional compact species cell indices.

The benchmark raster must match `removal_order.tif` exactly in geometry and cell count. No resampling is performed in stage 7.2.

Every initially alive cell must have a finite benchmark rank. This is a deliberate hard failure because silently dropping cells would compare different landscapes.

#### 2. Convert stage targets to benchmark cutoffs

The alive-domain benchmark cells are sorted by:

```text
descending rank value, then ascending cell index
```

The cell-index tie breaker makes equal-rank behavior deterministic.

`find_keep_n_by_cells()` converts each stage's `alive_end` to a bounded integer count. The helper records:

- Target and actual retained cells.
- Target and actual retained area.
- Area mismatch.
- Rank cutoff.

Matching is by cell count, not area. Area mismatch is retained as a diagnostic rather than optimized away.

#### 3. Cache initial species cells

`get_species_cells_from_patch_raster()` reads occupied cells from patch rasters when needed. When available, the bundle's `patch_cell_index_by_species_list` is preferred.

Each species cell vector is checked for:

- Non-zero length.
- Valid raster indices.
- Membership in the initial alive domain.

#### 4. Write the shared stage-0 LUT

The baseline LUT is derived from `bundle$patch_table`, restricted to retained species, standardized to the benchmark LUT schema, labeled `method = "rank"`, and saved as `lut_stage_0000.rds`.

Both methods therefore begin from the same stage-0 patch/PU representation.

#### 5. Reconstruct each benchmark stage

The retained benchmark mask is updated incrementally by `update_retained_mask()`. Because stage targets are non-increasing, the mask usually removes lower-ranked cells as stages advance.

For each species:

1. Intersect its initial cells with the current retained mask.
2. Reuse the previous species LUT when its retained-cell count has not changed.
3. Otherwise call `build_lut_one_species()`.

`build_lut_one_species()` calls:

- `make_binary_raster_from_cells()`
- `clump_rook()`
- `patch_areas_from_clump()`
- `polygonize_kept_patches()`
- `assign_population_units()`

The function:

1. Builds a binary retained-habitat raster.
2. Labels rook-connected patches with `fasterRaster::clump(..., diagonal = FALSE)`.
3. Calculates patch area and cell count.
4. Removes patches at or below `min_patch_size`.
5. Polygonizes retained patches.
6. Links patches within `dispersal_dist`.
7. Assigns graph components as candidate PUs.
8. Removes PUs at or below `min_pop_size`.
9. Returns a standardized patch-to-PU LUT.

The stage LUT is the row-bind of all species LUTs.

The fasterRaster/GRASS step is limited to rook clumping. Stage 7.2 still uses
`terra::zonal()` and `terra::freq()` for area and cell-count summaries,
`terra::as.polygons()` for polygonization, `sf::st_is_within_distance()` for
dispersal-neighbor detection, and `igraph::components()` for PU labels. This
means fasterRaster can reduce the connected-patch labeling cost without
changing the downstream representation or the persistence-comparison contract.

Two caches reduce repeated work:

- `species_cells` stores the initial occupied cell indices for each retained
  species, preferably from the saved script-6 bundle.
- `prev_n_cells` and `prev_lut` reuse a species' previous LUT when the current
  benchmark stage retains the same number of cells for that species.

#### 6. Finish benchmark LUT construction

The target summary and stage LUT files are the complete stage-7.2 contract.
Stage 7.3 checks the target summary and LUT schema when it loads them, so stage
7.2 does not maintain a separate inventory file or reread every output at the
end of the run.

### Outputs

```text
Data/Results/PriorityRuns/<run_id>/ana/rank_lut_<rank_method>/
|-- rank_stage_targets_<rank_method>.csv
|-- lut_stage_0000.rds
|-- lut_stage_0001.rds
`-- ...
```

LUT schema:

```text
stage
method
scientificName
species
patch_id
pu_id
patch_area_km2
patch_n_cells
```

## Stage 7.3: Persistence Comparison

### Driver and helpers

Driver:

- `7.3_persist_cmp.Rmd`

Standalone single-species wrapper:

- `R/stage7_persistence_compare.R`

Implementation:

- `R/stage7_persistence_compare_core.R`
- `R/stage7_persistence_compare_stats.R`
- `R/stage7_persistence_compare_figures.R`
- `R/stage7_persistence_compare_report.R`

### Purpose

Stage 7.3 applies the same PU- and species-level persistence model to:

- The persistence-based stage-6 landscapes (`method = "pipe"`).
- One benchmark method's stage-7.2 landscapes (`method = "rank"`).

It writes long-form trajectories, stage summaries, paired differences, detailed
manuscript statistics, a console report, and four-panel comparison figures for
the configured illustrative species.

### Orchestration order

The Rmd sources the implementation helpers directly. It first sources:

1. Core calculations, quietly.
2. Statistics, quietly.

It then loops over `params$illustrative_species`, sets each species as the
current `focus_species`, and sources the figure helper for each one. The console
report is printed only for the first illustrative species. If a four-panel
figure exists, the Rmd sources `figure_utils.R` and prints it in the knitted
document; only the first illustrative-species figure is saved to
`Figures/results_comparison_4panel.png`.

`R/stage7_persistence_compare.R` remains available as a single-species wrapper
that sources the same core, statistics, figure, and report helpers in that
order.

### Core chronological execution

#### 1. Resolve run paths and species parameters

The core helper validates `rank_method` with
`validate_benchmark_rank_method()` against:

```text
abf, caz1, caz2, cazmax
```

It loads the bundle, stage metadata, and species table. The selected curve's alpha/beta columns are standardized as `a_pred` and `b_pred`. The area threshold is:

```text
c_th = quasi_extinction_abundance() / density
```

The species-table areas are validated against density and the canonical
abundance constants before persistence is calculated.
`benchmark_artifact_paths()` and `persistence_comparison_paths()` provide the
shared stage-7 path contract.

#### 2. Discover method inputs

`pipe_stage_path()` resolves:

- Stage 0 from `bundle$patch_table`.
- Later pipeline stages from the exact path in `stage_meta.csv` when available, otherwise from the canonical filename.

`rank_stage_path()` resolves each stage-7.2 RDS.

All expected stage files must exist before calculations begin.

#### 3. Standardize LUTs

`standardize_rank_lut()` in `R/priority_run_config.R` accepts scientific names or safe species IDs, requires `pu_id` and `patch_area_km2`, creates a patch ID if absent, filters to retained species, validates positive areas, and returns one common schema.

`load_pipe_lut()` and `load_rank_lut()` wrap CSV and RDS reading respectively.

#### 4. Calculate PU persistence

`pu_persist()` evaluates:

```text
P(K) = exp[-alpha * (K - 500)^(-beta)]
```

with `K = density * area`. In area form:

```text
P(A) = exp[
  -alpha * density^(-beta) * (A - 500/density)^(-beta)
]
```

PUs at or below the threshold area receive persistence zero. Invalid values are set to zero and final results are bounded to `[0, 1]`.

#### 5. Calculate species persistence

`sp_persist()` combines PU probabilities under the complement-of-joint-extinction expression:

```text
P(species persists) = 1 - product(1 - P(PU persists))
```

The implementation uses `log1p()` for numerical stability. A species with no PUs receives persistence zero.

#### 6. Compute both methods

`compute_method()` loops over every stage:

1. Loads and standardizes the stage LUT.
2. Sums patch area within species-PU combinations.
3. Joins demographic parameters.
4. Computes `P_pu`.
5. Aggregates PU persistence to species persistence.
6. Adds zero-valued rows for species absent from the stage landscape.

This fixed-denominator design guarantees one species row per retained species, method, and stage.

The function is called once for `pipe` and once for `rank`.

#### 7. Write core comparison tables

`stage_sum` calculates:

- Species count.
- Species with at least one PU.
- Mean, median, 10th percentile, and 90th percentile persistence.
- Minimum and maximum persistence.
- Fraction above 0.5.
- Total PU area.

`stage_compare` reshapes selected summaries side by side and adds pipeline-minus-rank differences.

When `overwrite = false` and all four core files exist, they are loaded instead of recomputed.

#### 8. Validate output completeness

The core helper requires:

- Both methods.
- Exactly one species row for every method-stage-species combination.
- The expected stage count for each method.
- PU and species persistence values in `[0, 1]`.

### Statistics

`stage7_persistence_compare_stats.R` computes six groups of reported metrics:

1. Late-stage mean persistence and persistence-loss reduction.
2. Species-level signed area between method trajectories.
3. Area-normalized persistence advantage.
4. Interpolated threshold-crossing asymmetry.
5. Checkpoint tail-risk counts.
6. Timing of PU redundancy loss.

Important helpers include:

- `trapz_xy()` for signed trajectory integration.
- `area_norm_one_species()` for comparison on a common retained-area grid.
- `interp_first_crossing()` for continuous threshold crossing.
- `first_condition_time()` for redundancy-loss timing.
- `resolve_species_id()` for stable focal-species lookup.

### Figures

`stage7_persistence_compare_figures.R` builds:

- Panel A: mean persistence across the removal sequence.
- Panel B: focal-species persistence.
- Panel C: focal-species PU area composition.
- Panel D: focal-species persistence versus total retained area.

It also calculates signed species AUC values used for ranking examples. `panel_header()`, `add_header_strip()`, and `tag_panel()` assemble the final publication layout.

### Console report

`stage7_persistence_compare_report.R` formats the major statistics and prints selected best/worst examples. Helpers such as `print_auc_examples()`, `print_area_examples()`, `print_cross_examples()`, `print_tail_examples()`, and `print_redundancy_examples()` keep the report compact and consistent.

### Outputs

```text
Data/Results/PriorityRuns/<run_id>/ana/persist_cmp_<rank_method>/
|-- pu_long_<rank_method>.csv
|-- sp_long_<rank_method>.csv
|-- stage_sum_<rank_method>.csv
`-- stage_compare_<rank_method>.csv
```

Figure:

- `Figures/results_comparison_4panel.png`

## Manuscript Figures

`figure_utils.R` contains `manuscript_figure_manifest()`, which maps stable IDs used in code to output filenames. `save_manuscript_figure()` creates `Figures/`, chooses the configured dimensions and DPI, and writes a white-background PNG.

| File | Stage | Content |
|---|---:|---|
| `Figures/bird_sigma.png` | 1 | Bird environmental variability versus body mass and diet. |
| `Figures/S2-allometric-calibration-relationships.png` | 1 | Mammal and bird demographic allometries. |
| `Figures/S3-bird-sigma-model-comparison.png` | 1 | Candidate bird sigma model comparison. |
| `Figures/area_curve.png` | 4 | Focal-species abundance-persistence curve in area units. |
| `Figures/gompertz_wolff.png` | 3 | Shifted-Gompertz and Wolff comparison. |
| `Figures/S3-mammal-loess-parameters-all-curves.png` | 3 | Mammal alpha/beta LOESS relationships. |
| `Figures/S3-bird-loess-parameters-by-diet.png` | 3 | Bird alpha/beta LOESS relationships by diet. |
| `Figures/spatial_process.png` | 5.1 | Habitat-to-population-unit workflow. |
| `Figures/results_comparison_4panel.png` | 7.3 | Persistence-based versus benchmark comparison. |

The manuscript source files reference generated figure assets, but figure
generation does not automatically update manuscript text.

## Output Contracts

### Stage 5 patch lookup

```text
scientificName, patch_id, pu_id, patch_area_km2
```

Each row describes one patch. Multiple patches can share a `pu_id` within a species.

### Stage 6 removal events

```text
removal_step
stage
event_type
pruning_iteration
cells_removed
area_removed_km2
cum_cells_removed
cum_area_removed_km2
cum_prop_cells_removed
cum_prop_area_removed
cells_retained
area_retained_km2
```

`event_type` can be:

- `prune`
- `fragmentation`
- `distance`
- `final_retained`

### Removal order and rank map

`removal_order.tif` stores the chronological removal step. Larger values were retained longer.

`rankmap.tif` stores cumulative cell-removal proportion. It is in `[0, 1]`, with the terminal retained layer equal to 1.

Both rasters:

- Have identical geometry.
- Are finite over exactly the initial alive-cell domain.
- Are `NA` outside that domain.

### Stage lookup tables

Pipeline stage CSVs and benchmark stage RDS files both represent patch-to-PU structure at a completed stage. Stage 7.3 standardizes them before persistence calculation.

### Persistence comparison tables

`pu_long_<method>.csv` contains one row per represented PU and stage, including area, patch count, and `P_pu`.

`sp_long_<method>.csv` contains one row per retained species, stage, and method, including species persistence, PU count, and total PU area. Species absent from a landscape are explicitly represented with zeros.

`stage_sum_<method>.csv` contains method-level summaries at each stage.

`stage_compare_<method>.csv` places pipeline and benchmark summaries side by side and includes pipeline-minus-benchmark differences.

