library(targets)
library(terra)
library(sf)
library(stringr)
library(dplyr)
library(ggplot2)
library(lubridate)

source("config/Config_CM.R")
source("functions/_functions_unified.R")

tar_option_set(packages = c("terra","sf","dplyr","ggplot2","tidyr","stringr"), format = "rds")

list(
        # Step 1: Load rasters
        tar_target(nc_files, filter_nc_files(data_dir, variables, scenarios, models)),
        tar_target(nc_hist_files, nc_files[str_detect(tolower(nc_files), tolower(scenarios[1]))]),
        tar_target(nc_proj_files, nc_files[str_detect(tolower(nc_files), tolower(scenarios[2]))]),
        tar_target(city_shape, {
                if (!"shapes" %in% selected_targets) return(NULL)
                get_abs_geography(gpkg_path, abs_filter_column, abs_filter_values[1])
        })
        # ,
        # 
        # tar_target(t_raster_hist_file, {
        #         if (!"total" %in% selected_targets) return(NULL)
        #         r <- load_and_crop_rasters(nc_hist_files, variables[1], hist_start_year, hist_end_year,
        #                                    test_mode, if ("shapes" %in% selected_targets) city_shape else NULL,
        #                                    if (isTRUE(test_mode)) sample_layers else NULL)
        #         if (is.null(r)) return(NULL)
        #         path <- file.path(output_dir, paste0("hist_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
        #         safe_write_raster(r, path)
        # }, format = "file"),
        # 
        # tar_target(t_raster_proj_file, {
        #         if (!"total" %in% selected_targets) return(NULL)
        #         r <- load_and_crop_rasters(nc_proj_files, variables[1], proj_start_year, proj_end_year,
        #                                    test_mode, if ("shapes" %in% selected_targets) city_shape else NULL,
        #                                    if (isTRUE(test_mode)) sample_layers else NULL)
        #         if (is.null(r)) return(NULL)
        #         path <- file.path(output_dir, paste0("proj_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
        #         safe_write_raster(r, path)
        # }, format = "file"),
        # 
        # # Step 2: Aggregations
        # tar_target(raster_hist_mean_file, {
        #         if (!"mean" %in% selected_targets) return(NULL)
        #         if (is.null(t_raster_hist_file)) return(NULL)
        #         r <- terra::rast(t_raster_hist_file)
        #         r_mean <- deep_mean(r)
        #         path <- file.path(output_dir, paste0("hist_mean_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
        #         safe_write_raster(r_mean, path)
        # }, format = "file"),
        # 
        # tar_target(raster_proj_mean_file, {
        #         if (!"mean" %in% selected_targets) return(NULL)
        #         if (is.null(t_raster_proj_file)) return(NULL)
        #         r <- terra::rast(t_raster_proj_file)
        #         r_mean <- deep_mean(r)
        #         path <- file.path(output_dir, paste0("proj_mean_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
        #         safe_write_raster(r_mean, path)
        # }, format = "file"),
        # 
        # tar_target(raster_delta_mean_file, {
        #         if (!"delta_mean" %in% selected_targets) return(NULL)
        #         if (is.null(raster_hist_mean_file) || is.null(raster_proj_mean_file)) return(NULL)
        #         delta <- compute_delta_raster(raster_proj_mean_file, raster_hist_mean_file)
        #         path <- file.path(output_dir, paste0("delta_mean_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
        #         safe_write_raster(delta, path)
        # }, format = "file"),
        # 
        # tar_target(raster_monthly_hist_file, {
        #         if (!"monthly" %in% selected_targets) return(NULL)
        #         safe_aggregate_raster(t_raster_hist_file, hist_start_year, hist_end_year, "monthly", output_dir, "monthly_hist")
        # }, format = "file"),
        # 
        # tar_target(raster_monthly_proj_file, {
        #         if (!"monthly" %in% selected_targets) return(NULL)
        #         safe_aggregate_raster(t_raster_proj_file, proj_start_year, proj_end_year, "monthly", output_dir, "monthly_proj")
        # }, format = "file"),
        # 
        # tar_target(raster_annual_hist_file, {
        #         if (!"annual" %in% selected_targets) return(NULL)
        #         safe_aggregate_raster(t_raster_hist_file, hist_start_year, hist_end_year, "annual", output_dir, "annual_hist")
        # }, format = "file"),
        # 
        # tar_target(raster_annual_proj_file, {
        #         if (!"annual" %in% selected_targets) return(NULL)
        #         safe_aggregate_raster(t_raster_proj_file, proj_start_year, proj_end_year, "annual", output_dir, "annual_proj")
        # }, format = "file"),
        # 
        # # Percentile rasters
        # tar_target(raster_percentile_hist_file, {
        #         if (!"percentile" %in% selected_targets) return(NULL)
        #         r <- terra::rast(t_raster_hist_file)
        #         perc <- percentile_map(r, Percentile / 100)
        #         if (is.null(perc)) return(NULL)
        #         path <- file.path(output_dir, paste0("hist_percentile_", Percentile, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
        #         safe_write_raster(perc, path)
        # }, format = "file"),
        # 
        # tar_target(raster_percentile_proj_file, {
        #         if (!"percentile" %in% selected_targets) return(NULL)
        #         r <- terra::rast(t_raster_proj_file)
        #         perc <- percentile_map(r, Percentile / 100)
        #         if (is.null(perc)) return(NULL)
        #         path <- file.path(output_dir, paste0("proj_percentile_", Percentile, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
        #         safe_write_raster(perc, path)
        # }, format = "file"),
        # 
        # # Percentile plots
        # tar_target(plot_percentile_hist, {
        #         if (!"plots" %in% selected_targets|| !"percentile" %in% selected_targets) 
        #                 return(NULL)
        #         if (is.null(raster_percentile_hist_file) || 
        #             !file.exists(raster_percentile_hist_file)) return(NULL)
        #         if(file.info(raster_percentile_hist_file)$size==0){
        #                 message("Skipping percentile plot: raster file is empty.")
        #                 return(NULL)}
        #         perc_df <- to_df(terra::rast(raster_percentile_hist_file))
        #         plot_extent <- as.vector(terra::ext(terra::rast(raster_percentile_hist_file)))
        #         plot_percentile_map(perc_df, plot_extent, variables[1], Percentile,
        #                             paste(hist_start_year, hist_end_year), output_dir)
        # }),
        # 
        # tar_target(plot_percentile_proj, {
        #         if (!"plots" %in% selected_targets || !"percentile" %in% selected_targets)
        #                 return(NULL)
        #         if (is.null(raster_percentile_proj_file) || 
        #             !file.exists(raster_percentile_proj_file)) return(NULL)
        #         if(file.info(raster_percentile_proj_file)$size==0){
        #                 message("Skipping percentile plot: raster file is empty.")
        #                 return(NULL)}
        #         perc_df <- to_df(terra::rast(raster_percentile_proj_file))
        #         plot_extent <- as.vector(terra::ext(terra::rast(raster_percentile_proj_file)))
        #         plot_percentile_map(perc_df, plot_extent, variables[1], Percentile,
        #                             paste(proj_start_year, proj_end_year), output_dir)
        # }),
        # 
        # # Polygon summaries
        # tar_target(polygon_summary_hist, {
        #         if (!"shapes" %in% selected_targets) return(NULL)
        #         if (is.null(city_shape) || is.null(raster_hist_mean_file)) return(NULL)
        #         r <- terra::rast(raster_hist_mean_file)
        #         aggregate_by_polygon(r, city_shape, paste0(variables[1], "_hist_mean"), output_dir)
        # }),
        # 
        # tar_target(polygon_summary_proj, {
        #         if (!"shapes" %in% selected_targets) return(NULL)
        #         if (is.null(city_shape) || is.null(raster_proj_mean_file)) return(NULL)
        #         r <- terra::rast(raster_proj_mean_file)
        #         aggregate_by_polygon(r, city_shape, paste0(variables[1], "_proj_mean"), output_dir)
        # }),
        # 
        # 
        # # Step 3: Plots
        # tar_target(plot_mean_comparison, {
        #         if (!"plots" %in% selected_targets) return(NULL)
        #         
        #         # Check files exist
        #         if (is.null(raster_hist_mean_file) || is.null(raster_proj_mean_file)) return(NULL)
        #         if (!file.exists(raster_hist_mean_file) || !file.exists(raster_proj_mean_file)) return(NULL)
        #         
        #         # Load rasters
        #         hist_r <- terra::rast(raster_hist_mean_file)
        #         proj_r <- terra::rast(raster_proj_mean_file)
        #         
        #         # Compute delta directly for plotting
        #         delta_r <- tryCatch({
        #                 if (!is.null(raster_delta_mean_file) && file.exists(raster_delta_mean_file)) {
        #                         terra::rast(raster_delta_mean_file)
        #                 } else {
        #                         proj_r - hist_r
        #                 }
        #         }, error = function(e) {
        #                 message("Delta calculation failed: ", e$message)
        #                 return(NULL)
        #         })
        #         
        #         # Convert rasters to data frames
        #         hist_df <- to_df(hist_r)
        #         proj_df <- to_df(proj_r)
        #         delta_df <- if (!is.null(delta_r)) to_df(delta_r) else NULL
        #         
        #         # Compute plot extent (fallback to raster extent)
        #         plot_extent <- tryCatch({
        #                 if ("shapes" %in% selected_targets && !is.null(city_shape)) {
        #                         as.vector(sf::st_bbox(sf::st_transform(city_shape, crs = terra::crs(hist_r))))
        #                 } else {
        #                         as.vector(terra::ext(hist_r))
        #                 }
        #         }, error = function(e) {
        #                 message("Extent calculation failed, using default.")
        #                 c(-180,180,-90,90)
        #         })
        #         
        #         # Generate plot
        #         plot_raster_comparison(hist_df, proj_df, delta_df, plot_extent, variables[1],
        #                                paste(hist_start_year, hist_end_year),
        #                                paste(proj_start_year, proj_end_year), output_dir)
        # }),
        # 
        # tar_target(plot_monthly_hist, {
        #         if (!"plots" %in% selected_targets) return(NULL)
        #         if (is.null(raster_monthly_hist_file) || !file.exists(raster_monthly_hist_file)) return(NULL)
        #         
        #         r <- terra::rast(raster_monthly_hist_file)
        #         names(r) <- paste0("monthly_", seq_len(nlyr(r)))  # Force correct names
        #         monthly_df <- to_df(r)
        #         monthly_df$month <- factor(monthly_df$layer, levels = paste0("monthly_", 1:12), labels = month.abb)
        #         plot_extent <- as.vector(terra::ext(r))
        #         plot_monthly_climatology(monthly_df, plot_extent, variables[1],
        #                                  paste(hist_start_year, hist_end_year), output_dir)
        # }),
        # 
        # tar_target(plot_monthly_proj, {
        #         if (!"plots" %in% selected_targets) return(NULL)
        #         if (is.null(raster_monthly_proj_file) || !file.exists(raster_monthly_proj_file)) return(NULL)
        #         
        #         r <- terra::rast(raster_monthly_proj_file)
        #         names(r) <- paste0("monthly_", seq_len(nlyr(r)))
        #         monthly_df <- to_df(r)
        #         monthly_df$month <- factor(monthly_df$layer, levels = paste0("monthly_", 1:12), labels = month.abb)
        #         plot_extent <- as.vector(terra::ext(r))
        #         plot_monthly_climatology(monthly_df, plot_extent, variables[1],
        #                                  paste(proj_start_year, proj_end_year), output_dir)
        # }),
        # 
        # tar_target(plot_annual_hist, {
        #         if (!"plots" %in% selected_targets) return(NULL)
        #         if (is.null(raster_annual_hist_file) || !file.exists(raster_annual_hist_file)) return(NULL)
        #         
        #         r <- terra::rast(raster_annual_hist_file)
        #         names(r) <- paste0("annual_", seq_len(nlyr(r)))
        #         if (is.null(names(r)) || !all(grepl("annual", names(r)))) {
        #                 names(r) <- paste0("annual_", seq_len(nlyr(r)))
        #         }
        #         annual_df <- to_df(r)
        #         plot_extent <- as.vector(terra::ext(r))
        #         plot_annual_or_monthly(annual_df, plot_extent, variables[1], "annual", output_dir)
        # }),
        # 
        # tar_target(plot_annual_proj, {
        #         if (!"plots" %in% selected_targets) return(NULL)
        #         if (is.null(raster_annual_proj_file) || !file.exists(raster_annual_proj_file)) return(NULL)
        #         
        #         r <- terra::rast(raster_annual_proj_file)
        #         names(r) <- paste0("annual_", seq_len(nlyr(r)))
        #         if (is.null(names(r)) || !all(grepl("annual", names(r)))) {
        #                 names(r) <- paste0("annual_", seq_len(nlyr(r)))
        #         }
        #         annual_df <- to_df(r)
        #         plot_extent <- as.vector(terra::ext(r))
        #         plot_annual_or_monthly(annual_df, plot_extent, variables[1], "annual", output_dir)
        # })
)