
#----
# Climate menu Smoke Test
#----
# ----------------------
# Climate Menu recovery smoke test
# ----------------------

data_dir <- "data"

gpkg_path <-
        "data/boundaries/GCCSA_2016/GCCSA_2016_AUST.shp"

variables <- c("FFDI")

scenarios <- c(
        "historical",
        "ssp370"
)

models <- c(
        "EC-Earth3-Veg"
)

# ----------------------
# NARCliM model selection
# ----------------------

gcms <- c(
        "EC-Earth3-Veg"
)

rcms <- c(
        "WRF412R3"
)

model_mode <- "ensemble_R3" #individual"

# Allowed:
# "individual"
# "ensemble_R3"
# "ensemble_R5"
# "ensemble_all"

hist_start_year <- 1990
hist_end_year   <- 1990

proj_start_year <- 2020
proj_end_year   <- 2020

Percentile <- 99

threshold_value <- 50
threshold_var   <- "TXge"
min_duration    <- 3

abs_filter_column <- "GCC_NAME16"
abs_filter_values <- c("Greater Perth")

# FIRST RECOVERY TEST:
# deliberately exclude annual/monthly until temporal aggregation is repaired.
selected_targets <- c(
        "total",
        "mean",
        "annual",
        "monthly",
        "percentile",
        "delta_mean",
        "threshold",
        "events",
        "plots",
        "shapes"
)

precip_in_mm_day <- FALSE

test_mode <- FALSE

# Enough layers for a smoke test, but results are TEST results only.
sample_layers <- NULL

# =============================================================================
# GEOGRAPHY / BOUNDARY SETTINGS
# =============================================================================

boundary_path <-
        "data/boundaries/GCCSA_2016/GCCSA_2016_AUST.shp"

# For shapefiles leave NULL.
# For GeoPackages specify the layer, e.g. "LGA_2021_AUST_GDA2020"
boundary_layer <- NULL

# Columns identifying each geography
boundary_id_column <- "GCC_CODE16"
boundary_name_column <- "GCC_NAME16"

# Column and values used to select the user's area(s)
boundary_filter_column <- "GCC_NAME16"
boundary_filter_values <- c(
        "Greater Perth"
)

output_dir <- "outputs/20261026_v06"

if (!dir.exists(output_dir)) {
        dir.create(
                output_dir,
                recursive = TRUE
        )
}

# =============================================================================
# V04 DATA SOURCE AND STUDY-AREA CONTRACT
# =============================================================================
#
# This layer separates:
#
#   1. WHERE climate data are acquired from
#   2. WHAT geographical area is analysed
#   3. HOW downstream climate metrics are calculated
#
# Existing Climate Menu scientific targets remain responsible for (3).
#
# IMPORTANT:
# Human-readable geography names are metadata only.
# They must not be used as functional target or raster filenames.
# =============================================================================


# -----------------------------------------------------------------------------
# DATA SOURCE
# -----------------------------------------------------------------------------

# Allowed:
#   "local"
#   "narclim_remote"
#
# Keep local as the default until the remote pathway has passed the
# full historical/future/change regression tests.

data_source <- "local"


# -----------------------------------------------------------------------------
# STUDY AREA / ROI
# -----------------------------------------------------------------------------

# Stable internal identifier.
#
# This protects downstream filenames and targets from long or changing
# administrative-area names.

roi_id <- "ROI_001"


# Allowed:
#   "named_geography"
#   "bbox"
#
# A future Shiny release can also construct named_geography from an
# uploaded user boundary without changing the downstream ROI contract.

roi_mode <- "named_geography"


# Number of NARCliM native-grid cells added around the initial retrieval
# window.
#
# Padding is an ACQUISITION property, not an expansion of the final
# analysis geography.

roi_padding_cells <- 1L


# When multiple polygons are selected:
#
# TRUE  = dissolve them into one study region for analysis
# FALSE = retain polygon identities for separate summaries
#
# A single climate retrieval window may still cover all selected polygons.

combine_geographies <- TRUE


# -----------------------------------------------------------------------------
# NAMED GEOGRAPHY
# -----------------------------------------------------------------------------

# Initial validated geography remains GCCSA.
#
# V04 will add a geography catalogue so the same interface can support:
#
#   GCCSA
#   LGA
#   SA1
#   SA2
#   SA3
#   SA4
#   State/Territory
#   custom boundaries
#
# Exact field names will be obtained from the supplied boundary datasets
# rather than guessed.

geography_type <- "GCCSA"

geography_ids <- NULL

geography_names <- c(
        "Greater Perth"
)


# -----------------------------------------------------------------------------
# OPTIONAL EXPLICIT COORDINATE BOUNDING BOX
# -----------------------------------------------------------------------------
#
# Used only when:
#
#   roi_mode == "bbox"
#
# Coordinates are WGS84 longitude/latitude.
#
# Leave NULL for named-geography mode.

roi_bbox <- list(
        lon_min = NULL,
        lon_max = NULL,
        lat_min = NULL,
        lat_max = NULL
)


# -----------------------------------------------------------------------------
# HISTORICAL / FUTURE COMPARISON
# -----------------------------------------------------------------------------
#
# These switches explicitly preserve the central Climate Menu workflow:
#
#   historical climate
#   future climate
#   future-minus-historical change
#
# Existing hist_start_year/hist_end_year and
# proj_start_year/proj_end_year remain authoritative.

comparison_mode <- "historical_future_change"

allowed_comparison_modes <- c(
        "historical",
        "future",
        "historical_future",
        "historical_future_change"
)


# -----------------------------------------------------------------------------
# V04 DEVELOPMENT SAFETY
# -----------------------------------------------------------------------------

climate_menu_version <- "20261026_v06"

enable_remote_targets <- FALSE


###############################################################################
# CLIMATE MENU v04 - NARCliM product / domain contract
#
# Scientific choices are represented separately:
#
#   domain:
#     AUS-18   = Australia-wide NARCliM2 domain (~20 km)
#     SEAus-04 = finer southeast/eastern NARCliM2 domain (~4 km)
#
#   product:
#     raw            = original NARCliM2 model output
#     bias_adjusted  = bias-adjusted NARCliM2 output
#
# Currently exposed Climate Menu combinations:
#
#   AUS-18   + raw
#   AUS-18   + bias_adjusted
#   SEAus-04 + raw
#
# SEAus-04 + bias_adjusted is deliberately NOT enabled unless/until
# availability is established from the NCI catalogue.
#
# These settings describe the requested climate product. They do not alter
# the established NARCliM streaming implementation.
###############################################################################

narclim_domain <- "AUS-18"
narclim_product <- "bias_adjusted"

allowed_narclim_domains <- c(
  "AUS-18",
  "SEAus-04"
)

allowed_narclim_products <- c(
  "raw",
  "bias_adjusted"
)

allowed_narclim_product_domain <- data.frame(
  narclim_domain = c(
    "AUS-18",
    "AUS-18",
    "SEAus-04"
  ),
  narclim_product = c(
    "raw",
    "bias_adjusted",
    "raw"
  ),
  approximate_resolution = c(
    "~20 km",
    "~20 km",
    "~4 km"
  ),
  stringsAsFactors = FALSE
)

validate_narclim_product_domain <- function(
  domain = narclim_domain,
  product = narclim_product
) {

  if (!domain %in% allowed_narclim_domains) {
    stop(
      "Unsupported NARCliM domain: ",
      domain
    )
  }

  if (!product %in% allowed_narclim_products) {
    stop(
      "Unsupported NARCliM product: ",
      product
    )
  }

  valid <- allowed_narclim_product_domain[
    allowed_narclim_product_domain$narclim_domain == domain &
      allowed_narclim_product_domain$narclim_product == product,
    ,
    drop = FALSE
  ]

  if (nrow(valid) != 1L) {
    stop(
      "NARCliM product/domain combination is not currently enabled: ",
      domain,
      " + ",
      product
    )
  }

  invisible(valid)
}

validate_narclim_product_domain()
