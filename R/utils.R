# -----------------------------------------------------------------------------
# utils.R — small helpers
# -----------------------------------------------------------------------------

#' Parse command-line arguments of the form --key=value
#'
#' @param defaults named list of default values. Every key in `defaults` is a
#'   valid option; unknown options raise an error. A default of NULL marks the
#'   option as required.
#' @return named list with the parsed values (all character).
parse_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  opts <- defaults
  for (a in args) {
    if (a %in% c("-h", "--help")) {
      cat("Options (default):\n")
      for (k in names(defaults)) {
        d <- defaults[[k]]
        cat(sprintf("  --%s=%s\n", k, if (is.null(d)) "<required>" else d))
      }
      quit(save = "no", status = 0)
    }
    if (!grepl("^--[^=]+=", a)) stop("Arguments must look like --key=value, got: ", a)
    key <- sub("^--([^=]+)=.*$", "\\1", a)
    val <- sub("^--[^=]+=", "", a)
    if (!key %in% names(defaults)) stop("Unknown option --", key)
    opts[[key]] <- val
  }
  missing <- names(defaults)[vapply(names(defaults), function(k) is.null(opts[[k]]), logical(1))]
  if (length(missing) > 0) {
    stop("Missing required option(s): ", paste0("--", missing, collapse = ", "),
         "\nRun with --help to see all options.")
  }
  opts
}

#' Split a comma-separated string into a character vector
split_csv <- function(x) trimws(strsplit(x, ",")[[1]])

#' Interpret "1"/"TRUE"/"yes" as TRUE
as_flag <- function(x) tolower(as.character(x)) %in% c("1", "true", "yes", "y")

#' Create a directory (recursively) if it does not exist and return its path
ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    message("Directory created: ", path)
  }
  path
}

#' Stop with an informative message if a file is missing
require_file <- function(path, hint = NULL) {
  if (!file.exists(path)) {
    stop("File not found: ", path, if (!is.null(hint)) paste0("\n", hint))
  }
  path
}

#' Number of parallel workers: "auto" -> all cores minus one
resolve_workers <- function(x) {
  if (identical(x, "auto")) max(1L, parallel::detectCores() - 1L) else as.integer(x)
}
