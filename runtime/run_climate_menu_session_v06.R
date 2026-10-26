suppressPackageStartupMessages(library(targets))

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1L) {
  stop("Expected exactly one argument: session_spec.rds")
}

session_spec_file <- normalizePath(
  args[[1]],
  winslash = "/",
  mustWork = TRUE
)

spec <- readRDS(session_spec_file)

required <- c(
  "project_root",
  "request_file",
  "store",
  "result_file",
  "output_dir",
  "session_id"
)

missing <- setdiff(required, names(spec))

if (length(missing)) {
  stop(
    "Session specification is missing: ",
    paste(missing, collapse = ", ")
  )
}

setwd(spec$project_root)

Sys.setenv(
  CLIMATE_MENU_REQUEST_FILE = spec$request_file,
  CLIMATE_MENU_SESSION_ID = spec$session_id,
  CLIMATE_MENU_OUTPUT_DIR = spec$output_dir
)

dir.create(spec$store, recursive = TRUE, showWarnings = FALSE)
dir.create(spec$output_dir, recursive = TRUE, showWarnings = FALSE)

cat("CLIMATE_MENU_SESSION_RUNNER_STARTED\n")
cat("SESSION_ID=", spec$session_id, "\n", sep = "")
cat("REQUEST_FILE=", spec$request_file, "\n", sep = "")
cat("TARGET_STORE=", spec$store, "\n", sep = "")
cat("OUTPUT_DIR=", spec$output_dir, "\n", sep = "")

targets_to_make <- c(
  "raster_hist_mean_file",
  "raster_proj_mean_file",
  "raster_delta_mean_file",
  "plot_mean_comparison"
)

targets::tar_make(
  names = tidyselect::any_of(targets_to_make),
  script = "Climate_Menu_targets_20261026_v06.R",
  store = spec$store,
  callr_function = NULL
)

read_target <- function(name) {
  tryCatch(
    targets::tar_read_raw(name, store = spec$store),
    error = function(e) NULL
  )
}

plots <- read_target("plot_mean_comparison")

if (is.null(plots)) {
  plots <- list()
}

plot_paths <- unlist(plots, use.names = TRUE)

if (length(plot_paths)) {
  plot_paths <- plot_paths[file.exists(plot_paths)]
}

find_plot <- function(pattern) {

  if (!length(plot_paths)) {
    return(NULL)
  }

  hit <- grepl(
    pattern,
    basename(plot_paths),
    ignore.case = TRUE
  )

  if (!any(hit) && length(names(plot_paths))) {
    hit <- grepl(
      pattern,
      names(plot_paths),
      ignore.case = TRUE
    )
  }

  if (!any(hit)) {
    return(NULL)
  }

  unname(plot_paths[which(hit)[1]])
}

result <- list(
  success = TRUE,
  request = readRDS(spec$request_file),
  rasters = list(
    historical = read_target("raster_hist_mean_file"),
    future = read_target("raster_proj_mean_file"),
    change = read_target("raster_delta_mean_file")
  ),
  plots = list(
    historical = find_plot("hist"),
    future = find_plot("proj|future"),
    change = find_plot("delta|change")
  ),
  completed_at = Sys.time()
)

saveRDS(result, spec$result_file)

cat("CLIMATE_MENU_SESSION_RUN=PASS\n")
