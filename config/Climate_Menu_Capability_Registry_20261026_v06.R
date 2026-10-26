###############################################################################
# CLIMATE MENU CAPABILITY REGISTRY
# Version: 20261026_v06
#
# PURPOSE
# -------
# User-facing scientific capability registry for Climate Menu.
#
# This separates:
#   1. the scientific question shown to the user;
#   2. the analysis Climate Menu performs; and
#   3. the underlying NARCliM source variable/product.
#
# The Shiny interface should read this registry rather than hard-code
# NARCliM variable names or scientific definitions.
###############################################################################

climate_menu_capability_registry <- data.frame(

  capability_id = c(
    "tasmax_mean",
    "tasmin_mean",
    "tasmax_percentile",
    "hot_days_threshold",
    "hot_event_duration",

    "precip_mean",
    "precip_percentile",
    "heavy_rain_days",

    "ffdi_mean",
    "ffdi_threshold_days",
    "ffdi_event_duration"
  ),

  exposure_family = c(
    rep("Temperature and heat", 5),
    rep("Precipitation", 3),
    rep("Fire weather", 3)
  ),

  display_name = c(
    "Mean maximum temperature",
    "Mean minimum temperature",
    "Extreme temperature percentile",
    "Hot days above threshold",
    "Consecutive hot-day events",

    "Mean precipitation",
    "Extreme precipitation percentile",
    "Heavy-rain days above threshold",

    "Mean fire weather danger",
    "Very high fire-weather days",
    "Consecutive high fire-weather events"
  ),

  description = c(
    "Average daily maximum near-surface air temperature over the selected period.",
    "Average daily minimum near-surface air temperature over the selected period.",
    "Selected percentile of daily maximum temperature over the selected period.",
    "Number of days on which daily maximum temperature reaches or exceeds a selected threshold.",
    "Frequency and duration of consecutive days on which daily maximum temperature reaches or exceeds a selected threshold.",

    "Average daily precipitation over the selected period.",
    "Selected percentile of daily precipitation over the selected period.",
    "Number of days on which daily precipitation reaches or exceeds a selected threshold.",

    "Average Forest Fire Danger Index over the selected period.",
    "Number of days on which FFDI reaches or exceeds a selected threshold.",
    "Frequency and duration of consecutive days on which FFDI reaches or exceeds a selected threshold."
  ),

  analysis_type = c(
    "mean",
    "mean",
    "percentile",
    "threshold",
    "events",

    "mean",
    "percentile",
    "threshold",

    "mean",
    "threshold",
    "events"
  ),

  preferred_product = c(
    rep("bias_adjusted", 8),
    rep("raw", 3)
  ),

  preferred_variable = c(
    "tasmaxAdjust",
    "tasminAdjust",
    "tasmaxAdjust",
    "tasmaxAdjust",
    "tasmaxAdjust",

    "prAdjust",
    "prAdjust",
    "prAdjust",

    "FFDI",
    "FFDI",
    "FFDI"
  ),

  raw_fallback_variable = c(
    "tasmax",
    "tasmin",
    "tasmax",
    "tasmax",
    "tasmax",

    "pr",
    "pr",
    "pr",

    "FFDI",
    "FFDI",
    "FFDI"
  ),

  threshold_required = c(
    FALSE,
    FALSE,
    FALSE,
    TRUE,
    TRUE,

    FALSE,
    FALSE,
    TRUE,

    FALSE,
    TRUE,
    TRUE
  ),

  percentile_required = c(
    FALSE,
    FALSE,
    TRUE,
    FALSE,
    FALSE,

    FALSE,
    TRUE,
    FALSE,

    FALSE,
    FALSE,
    FALSE
  ),

  event_duration_required = c(
    FALSE,
    FALSE,
    FALSE,
    FALSE,
    TRUE,

    FALSE,
    FALSE,
    FALSE,

    FALSE,
    FALSE,
    TRUE
  ),

  backend_status = c(
    "ready",
    "ready",
    "ready",
    "ready",
    "ready",

    "ready",
    "ready",
    "ready",

    "ready_local",
    "ready_local",
    "ready_local"
  ),

  stringsAsFactors = FALSE
)


###############################################################################
# SPATIAL CAPABILITY
###############################################################################

climate_menu_spatial_registry <- data.frame(

  spatial_mode = c(
    "named_geography",
    "multiple_geographies",
    "bbox"
  ),

  display_name = c(
    "Single named area",
    "Multiple named areas",
    "Custom bounding box"
  ),

  description = c(
    "Analyse one selected administrative or study geography.",
    "Analyse multiple selected areas separately or combine them into one study region.",
    "Analyse a user-defined longitude/latitude bounding box."
  ),

  backend_status = c(
    "ready",
    "ready",
    "ready"
  ),

  stringsAsFactors = FALSE
)


###############################################################################
# COMPARISON CAPABILITY
###############################################################################

climate_menu_comparison_registry <- data.frame(

  comparison_mode = c(
    "historical",
    "future",
    "historical_future",
    "historical_future_change"
  ),

  display_name = c(
    "Historical climate",
    "Future climate",
    "Historical and future",
    "Future change from historical"
  ),

  backend_status = "ready",

  stringsAsFactors = FALSE
)


###############################################################################
# MODEL CAPABILITY
###############################################################################

climate_menu_model_registry <- data.frame(

  model_mode = c(
    "individual",
    "ensemble_R3",
    "ensemble_R5",
    "ensemble_all"
  ),

  display_name = c(
    "Individual climate model",
    "R3 ensemble",
    "R5 ensemble",
    "All available models"
  ),

  backend_status = "ready",

  stringsAsFactors = FALSE
)


###############################################################################
# DATA SOURCE CAPABILITY
#
# Only combinations already validated for the Climate Menu should be enabled
# here. SEAus-04 bias-adjusted is deliberately retained but not yet enabled
# until its correct catalogue path is wired into discovery.
###############################################################################

climate_menu_source_registry <- data.frame(

  domain = c(
    "AUS-18",
    "AUS-18",
    "SEAus-04",
    "SEAus-04"
  ),

  product = c(
    "raw",
    "bias_adjusted",
    "raw",
    "bias_adjusted"
  ),

  approximate_resolution = c(
    "~20 km",
    "~20 km",
    "~4 km",
    "~4 km"
  ),

  enabled = c(
    TRUE,
    TRUE,
    TRUE,
    FALSE
  ),

  status = c(
    "validated",
    "validated",
    "validated",
    "known_to_exist_path_pending",
    "known_to_exist_path_pending"
  )[1:4],

  stringsAsFactors = FALSE
)


###############################################################################
# USER-FACING HELPERS
###############################################################################

climate_menu_exposure_families <- function() {

  unique(
    climate_menu_capability_registry$exposure_family[
      climate_menu_capability_registry$backend_status %in%
        c("ready", "ready_local")
    ]
  )
}


climate_menu_metrics <- function(exposure_family) {

  x <- climate_menu_capability_registry[
    climate_menu_capability_registry$exposure_family ==
      exposure_family &
    climate_menu_capability_registry$backend_status %in%
      c("ready", "ready_local"),
    ,
    drop = FALSE
  ]

  stats::setNames(
    x$capability_id,
    x$display_name
  )
}


climate_menu_capability <- function(capability_id) {

  x <- climate_menu_capability_registry[
    climate_menu_capability_registry$capability_id ==
      capability_id,
    ,
    drop = FALSE
  ]

  if (nrow(x) != 1L) {

    stop(
      "Unknown or duplicated Climate Menu capability: ",
      capability_id
    )
  }

  x
}


###############################################################################
# CONTRACT CHECKS
###############################################################################

stopifnot(
  !anyDuplicated(
    climate_menu_capability_registry$capability_id
  )
)

stopifnot(
  all(
    climate_menu_capability_registry$analysis_type %in%
      c(
        "mean",
        "percentile",
        "threshold",
        "events"
      )
  )
)

stopifnot(
  !anyDuplicated(
    climate_menu_source_registry[
      c("domain", "product")
    ]
  )
)

###############################################################################
# END
###############################################################################