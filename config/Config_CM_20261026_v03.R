
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

output_dir <- "outputs/20261024_v02"

if (!dir.exists(output_dir)) {
        dir.create(
                output_dir,
                recursive = TRUE
        )
}
