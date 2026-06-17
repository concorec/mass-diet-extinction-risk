# Shared helpers used by the numbered analysis drivers.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

assert <- function(ok, msg) {
  if (!isTRUE(ok)) stop(msg, call. = FALSE)
}

need_file <- function(path, label = path) {
  assert(file.exists(path), paste0("Missing file: ", label, "\nPath: ", path))
}

need_dir <- function(path, label = path) {
  assert(dir.exists(path), paste0("Missing directory: ", label, "\nPath: ", path))
}

need_cols <- function(x, cols, label) {
  missing_cols <- setdiff(cols, names(x))
  assert(
    length(missing_cols) == 0L,
    paste0(label, " is missing required column(s): ", paste(missing_cols, collapse = ", "))
  )
}

assert_has_cols <- need_cols

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  need_dir(path)
  invisible(path)
}

ensure_writable_dir <- function(path, label = path) {
  ensure_dir(path)
  probe <- tempfile("write_test_", tmpdir = path, fileext = ".tmp")
  writable <- tryCatch(file.create(probe), error = function(e) FALSE)
  if (isTRUE(writable)) unlink(probe)
  assert(writable, paste0("Directory is not writable: ", label, "\nPath: ", path))
  invisible(path)
}

optional_path <- function(path) {
  if (is.null(path) || !length(path)) return(NA_character_)
  path <- trimws(as.character(path[1L]))
  if (is.na(path) || !nzchar(path) || toupper(path) %in% c("NA", "NULL")) return(NA_character_)
  path.expand(path)
}

grass_config_dir <- function() {
  grass_bin <- Sys.which(c("grass", "grass84", "grass83", "grass82"))
  grass_bin <- grass_bin[nzchar(grass_bin)]
  if (!length(grass_bin)) return(NA_character_)

  for (bin in grass_bin) {
    out <- tryCatch(
      system2(bin, c("--config", "path"), stdout = TRUE, stderr = FALSE),
      error = function(e) character(0)
    )
    out <- optional_path(out[1L])
    if (!is.na(out) && dir.exists(out)) return(out)
  }

  NA_character_
}

resolve_grass_dir <- function(grass_dir = NULL) {
  explicit <- optional_path(grass_dir)
  env_candidates <- vapply(
    c("FASTER_RASTER_GRASS_DIR", "GRASS_DIR", "GISBASE"),
    function(x) optional_path(Sys.getenv(x, unset = NA_character_)),
    character(1)
  )
  path_candidate <- grass_config_dir()

  conda_prefix <- optional_path(Sys.getenv("CONDA_PREFIX", unset = NA_character_))
  conda_candidates <- character(0)
  if (!is.na(conda_prefix)) {
    conda_candidates <- c(
      Sys.glob(file.path(conda_prefix, "grass*")),
      Sys.glob(file.path(conda_prefix, "lib", "grass*"))
    )
  }

  os_candidates <- c(
    Sys.glob("C:/Program Files/GRASS GIS *"),
    Sys.glob("/usr/local/grass*"),
    Sys.glob("/usr/lib/grass*"),
    Sys.glob("/opt/grass*")
  )

  candidates <- unique(c(explicit, env_candidates, path_candidate, conda_candidates, os_candidates))
  candidates <- candidates[!is.na(candidates) & nzchar(candidates)]
  existing <- candidates[dir.exists(candidates)]

  if (!length(existing)) {
    stop(
      paste(
        "Could not find a GRASS GIS installation directory for fasterRaster.",
        "Set params$grass_dir, FASTER_RASTER_GRASS_DIR, GRASS_DIR, or GISBASE.",
        "For conda, activate the environment before rendering; if auto-detection fails,",
        "set FASTER_RASTER_GRASS_DIR to the GRASS directory, often $CONDA_PREFIX/lib/grass84.",
        sep = "\n"
      ),
      call. = FALSE
    )
  }

  normalizePath(existing[1L], winslash = "/", mustWork = TRUE)
}

resolve_optional_grass_dir <- function(grass_dir = NULL) {
  if (!requireNamespace("fasterRaster", quietly = TRUE)) {
    return(optional_path(grass_dir))
  }

  tryCatch(
    resolve_grass_dir(grass_dir),
    error = function(e) {
      warning(
        paste0(
          "Could not find a GRASS GIS installation for fasterRaster; ",
          "raster clumping will use terra::patches() fallback if GRASS clumping fails. ",
          "Original error: ", conditionMessage(e)
        ),
        call. = FALSE
      )
      NA_character_
    }
  )
}

init_faster_raster <- function(grass_dir = NULL, ...) {
  if (!requireNamespace("fasterRaster", quietly = TRUE)) {
    warning(
      "fasterRaster is not installed; raster clumping will use terra::patches() fallback.",
      call. = FALSE
    )
    return(invisible(NA_character_))
  }

  resolved_grass_dir <- tryCatch(
    resolve_grass_dir(grass_dir),
    error = function(e) {
      warning(
        paste0(
          "Could not initialize fasterRaster because GRASS GIS was not found; ",
          "raster clumping will use terra::patches() fallback. ",
          "Original error: ", conditionMessage(e)
        ),
        call. = FALSE
      )
      NA_character_
    }
  )

  if (is.na(resolved_grass_dir)) {
    return(invisible(NA_character_))
  }

  configure_grass_path(resolved_grass_dir)
  fasterRaster::faster(grassDir = resolved_grass_dir, ...)
  invisible(resolved_grass_dir)
}

configure_grass_path <- function(grass_dir) {
  grass_dir <- normalizePath(grass_dir, winslash = "/", mustWork = TRUE)
  grass_path_dirs <- file.path(grass_dir, c("bin", "extrabin", "scripts", "lib"))
  grass_path_dirs <- normalizePath(
    grass_path_dirs[dir.exists(grass_path_dirs)],
    winslash = "/",
    mustWork = TRUE
  )

  current_path <- strsplit(Sys.getenv("PATH"), .Platform$path.sep, fixed = TRUE)[[1]]
  current_path_norm <- normalizePath(
    current_path[nzchar(current_path) & dir.exists(current_path)],
    winslash = "/",
    mustWork = FALSE
  )
  add_path_dirs <- grass_path_dirs[!grass_path_dirs %in% current_path_norm]

  if (length(add_path_dirs)) {
    Sys.setenv(PATH = paste(c(add_path_dirs, Sys.getenv("PATH")), collapse = .Platform$path.sep))
  }
  Sys.setenv(GISBASE = grass_dir)

  invisible(grass_dir)
}

clump_binary_raster <- function(x, diagonal = FALSE) {
  directions <- if (isTRUE(diagonal)) 8L else 4L

  if (!requireNamespace("fasterRaster", quietly = TRUE)) {
    warning(
      "fasterRaster is not installed; using terra::patches() fallback.",
      call. = FALSE
    )
    return(terra::patches(x, directions = directions))
  }

  out <- tryCatch(
    {
      terra::rast(fasterRaster::clump(fasterRaster::fast(x), diagonal = isTRUE(diagonal)))
    },
    error = function(e) {
      warning(
        paste0(
          "fasterRaster::clump() failed; using terra::patches() fallback. ",
          "Original error: ", conditionMessage(e)
        ),
        call. = FALSE
      )
      terra::patches(x, directions = directions)
    }
  )

  out
}

check_packages <- function(packages) {
  packages <- unique(as.character(packages))
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  assert(
    length(missing) == 0L,
    paste0(
      "Missing required R package(s): ", paste(missing, collapse = ", "),
      "\nInstall them before running this script."
    )
  )
  invisible(packages)
}

load_packages <- function(packages) {
  packages <- check_packages(packages)
  suppressPackageStartupMessages(
    invisible(lapply(packages, library, character.only = TRUE))
  )
  invisible(packages)
}

species_id <- function(scientific_name) {
  x <- trimws(as.character(scientific_name))
  x <- gsub("\\s+", "_", x)
  x <- gsub("[/\\\\:<>\"|?*]+", "_", x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  gsub("^_|_$", "", x)
}

patch_filename_from_scientific <- function(scientific_name) {
  paste0(species_id(scientific_name), ".tif")
}

patch_file_from_name <- function(scientific_name, patch_dir) {
  file.path(patch_dir, patch_filename_from_scientific(scientific_name))
}

make_selection_tag <- function(parts) {
  parts <- as.character(parts)
  parts <- parts[nzchar(parts)]
  if (length(parts)) paste(parts, collapse = "_") else "none"
}

log_msg <- function(..., flush = FALSE) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), paste(..., collapse = " ")))
  if (isTRUE(flush)) flush.console()
  invisible(NULL)
}

log_space <- function(x_min, x_max, n) {
  assert(
    all(is.finite(c(x_min, x_max, n))) && x_min > 0 && x_max > x_min && n >= 2,
    "log_space() requires 0 < x_min < x_max and n >= 2."
  )
  10^seq(log10(x_min), log10(x_max), length.out = as.integer(n))
}

theme_pub <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.margin = ggplot2::margin(8, 8, 8, 8),
      axis.title = ggplot2::element_text(face = "bold", colour = "grey10"),
      axis.text = ggplot2::element_text(colour = "grey20"),
      axis.line = ggplot2::element_line(linewidth = 0.45, colour = "grey20"),
      axis.ticks = ggplot2::element_line(linewidth = 0.45, colour = "grey20"),
      legend.title = ggplot2::element_text(face = "bold", colour = "grey10"),
      legend.text = ggplot2::element_text(colour = "grey20"),
      legend.key = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, colour = "grey10")
    )
}
