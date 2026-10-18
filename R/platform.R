#' Supported detection platforms and their configurations
#'
#' Each platform defines: id (lowercase key), display_name, signal_label,
#' signal_unit, typical_range, recommended_methods, description.

.platform_env <- new.env(parent = emptyenv())

.platform_env$registry <- list(
  msd = list(
    id                = "msd",
    display_name      = "MSD",
    full_name         = "Meso Scale Discovery",
    signal_label      = "Signal",
    signal_unit       = "electrochemiluminescence units",
    typical_range     = c(10, 200000),
    recommended_methods = c("4pl", "5pl", "log"),
    description       = "Electrochemiluminescence-based multiplex immunoassay"
  ),
  luminex = list(
    id                = "luminex",
    display_name      = "Luminex",
    full_name         = "Luminex xMAP",
    signal_label      = "MFI",
    signal_unit       = "Median Fluorescence Intensity",
    typical_range     = c(10, 25000),
    recommended_methods = c("5pl", "4pl", "cubic"),
    description       = "Bead-based multiplex immunoassay (xMAP technology)"
  ),
  elisa = list(
    id                = "elisa",
    display_name      = "ELISA",
    full_name         = "Enzyme-Linked Immunosorbent Assay",
    signal_label      = "OD",
    signal_unit       = "Optical Density",
    typical_range     = c(0.01, 4.0),
    recommended_methods = c("4pl", "log", "quadratic"),
    description       = "Plate-based colorimetric/fluorometric immunoassay"
  ),
  simoa = list(
    id                = "simoa",
    display_name      = "Simoa",
    full_name         = "Single Molecule Array (Quanterix)",
    signal_label      = "AEB",
    signal_unit       = "Average Enzymes per Bead",
    typical_range     = c(0.001, 100),
    recommended_methods = c("4pl", "linear", "log"),
    description       = "Digital ELISA with ultra-high sensitivity"
  ),
  cba = list(
    id                = "cba",
    display_name      = "CBA",
    full_name         = "Cytometric Bead Array (BD)",
    signal_label      = "MFI",
    signal_unit       = "Median Fluorescence Intensity",
    typical_range     = c(5, 50000),
    recommended_methods = c("4pl", "5pl", "log"),
    description       = "Flow cytometry-based bead array immunoassay"
  ),
  ella = list(
    id                = "ella",
    display_name      = "ELLA",
    full_name         = "Simple Plex / ELLA (ProteinSimple)",
    signal_label      = "Signal",
    signal_unit       = "chemiluminescence units",
    typical_range     = c(100, 500000),
    recommended_methods = c("4pl", "5pl", "log"),
    description       = "Automated microfluidic immunoassay"
  )
)

#' Minimum data points required per fitting method
.method_min_points <- c(
  linear = 3, quadratic = 3, cubic = 4, log = 3,
  power = 3, exponential = 3, "4pl" = 5, "5pl" = 6
)

#' List all supported detection platforms
#'
#' @return character vector of platform ids (invisibly); prints summary to console
list_platforms <- function() {
  cat("Supported detection platforms:\n")
  for (p in .platform_env$registry) {
    cat(sprintf("  %-10s %-35s signal: %-10s range: [%s, %s]\n",
                p$id, p$full_name, p$signal_label,
                format(p$typical_range[1], scientific = FALSE),
                format(p$typical_range[2], scientific = FALSE)))
  }
  cat("\nRecommended methods per platform:\n")
  for (p in .platform_env$registry) {
    cat(sprintf("  %-10s %s\n", p$id, paste(p$recommended_methods, collapse = ", ")))
  }
  invisible(names(.platform_env$registry))
}

#' Get platform configuration
#'
#' @param platform character, platform name (case-insensitive)
#' @return list with platform configuration
get_platform <- function(platform) {
  platform <- tolower(trimws(platform))
  if (!platform %in% names(.platform_env$registry)) {
    stop(sprintf("Unknown platform '%s'. Available: %s",
                 platform, paste(names(.platform_env$registry), collapse = ", ")))
  }
  .platform_env$registry[[platform]]
}

#' Register a custom platform
#'
#' @param name character, short identifier (will be lowercased)
#' @param full_name character, full display name
#' @param signal_label character, column/axis label for signal values
#' @param signal_unit character, unit description
#' @param typical_range numeric vector of length 2 (min, max)
#' @param recommended_methods character vector of method names
#' @param description character, platform description
#' @return the platform id, invisibly
register_platform <- function(name, full_name, signal_label, signal_unit = "",
                              typical_range = c(1, 1e6),
                              recommended_methods = c("4pl", "linear"),
                              description = "") {
  name <- tolower(trimws(name))
  if (nchar(name) == 0) stop("Platform name must not be empty")
  if (!is.numeric(typical_range) || length(typical_range) != 2 ||
      typical_range[1] >= typical_range[2])
    stop("typical_range must be a numeric vector of length 2 with min < max")
  if (!is.character(recommended_methods) || length(recommended_methods) == 0)
    stop("recommended_methods must be a non-empty character vector")
  if (!is.character(signal_label) || nchar(signal_label) == 0)
    stop("signal_label must be a non-empty string")
  if (!is.character(full_name) || nchar(full_name) == 0)
    stop("full_name must be a non-empty string")

  if (name %in% names(.platform_env$registry)) {
    warning(sprintf("Platform '%s' already exists; overwriting", name))
  }
  .platform_env$registry[[name]] <- list(
    id                = name,
    display_name      = full_name,
    full_name         = full_name,
    signal_label      = signal_label,
    signal_unit       = signal_unit,
    typical_range     = typical_range,
    recommended_methods = recommended_methods,
    description       = description
  )
  invisible(name)
}