# ----------------------
# User Settings
# ----------------------
#
#
data_dir <- "D:/NARCLIM2.0_Data/Fire_Index/WA/all_data"
#data_dir <- "D:/NARCLIM2.0_Data/Fire_Index/WA/testing"
#gpkg_path<- "D:/NARCLIM2.0_Data/Shape_Files/ASGS_2021_MAIN_STRUCTURE_GPKG_GDA2020/ASGS_2021_Main_Structure_GDA2020.gpkg"
gpkg_path<-"C:\\Users\\300888D\\Nextcloud\\Environment_General\\ABS_data\\ABS_Census_2016\\abs_gcc_2016_data_provided\\GCCSA_2016_AUST.shp"

variables <- c("FFDI")
scenarios <- c("historical","ssp370")
models <- c("EC-Earth3-Veg")

hist_start_year <- 1990
hist_end_year <- 2005
proj_start_year <- 2020
proj_end_year <- 2039

Percentile <- 99

# heatwave settings
threshold_value <- 35  #308.15
threshold_var <- "TXge"
min_duration <- 3

# Geopackage containing multiple shape layers eg GCCSA, SA2 etc
#abs_level <- "GCCSA" 
#abs_filter_column <- "GCCSA_NAME_2021"
#abs_filter_value <- c("Greater Perth","Greater Sydney") # example multiple locations

# file with only GCC
abs_filter_column <- "GCC_NAME16"
abs_filter_values <- c("Greater Perth")#,"Greater Sydney") # example multiple locations


#selected_targets <- c("mean","annual","monthly","delta_mean","plots","percentile","total","days35","precip","HWD","HWF","HWN"#"shapes")
selected_targets <- c("mean","annual","monthly","plots","total", "percentile","shapes")


precip_in_mm_day <- FALSE

test_mode <- TRUE
#crop_extent <- terra::ext(112,120,-35,-25)
sample_layers <- 12 
extent <- c(112,129,-35.2,-13.8)
output_dir <- "outputs_test"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive=TRUE)
