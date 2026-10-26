library(targets)
library(terra)
library(sf)
library(stringr)
library(dplyr)
library(ggplot2)
library(lubridate)

#source("config/Config_CM.R")
source("config/Config_CM_20261026_v03.R") # Test config menu
source("functions/_functions_20261026_v03.R")

tar_option_set(packages = c("terra","sf","dplyr","ggplot2","tidyr","stringr"), format = "rds")

list(
        # Step 1: Load rasters
        tar_target(nc_files, filter_nc_files(data_dir=data_dir, variables=variables, scenarios=scenarios, models=NULL)),
        #====
        # NARCliM member discovery and user selection
        #====
        tar_target(narclim_member_table, parse_narclim_member(nc_files)),
        
        tar_target(selected_member_table, select_narclim_members(member_table = narclim_member_table,
                        gcms = if (model_mode == "individual") {gcms} else {NULL},
                        rcms = if (model_mode == "individual") {rcms} else {NULL},
                        model_mode = model_mode)),
        tar_target(selected_projection_members,
                selected_member_table[
                        tolower(selected_member_table$scenario) ==
                                tolower(scenarios[2]),
                        ,
                        drop = FALSE]),
        
        # ==========================================================
        # Projection member-level processing
        #
        # Each selected GCM x RCM projection member is processed
        # independently before any ensemble statistic is calculated.
        # ==========================================================
        
        tar_target(
                projection_member_index,
                seq_len(
                        nrow(selected_projection_members)
                )
        ),
        
        tar_target(
                projection_member_mean_file,
                {
                        
                        member_row <- selected_projection_members[
                                projection_member_index,
                                ,
                                drop = FALSE
                        ]
                        
                        process_member_mean(
                                member_row = member_row,
                                start_year = proj_start_year,
                                end_year = proj_end_year,
                                
                                city_shape =
                                        if ("shapes" %in% selected_targets) {
                                                city_shape
                                        } else {
                                                NULL
                                        },
                                
                                test_mode = test_mode,
                                
                                sample_layers =
                                        if (isTRUE(test_mode)) {
                                                sample_layers
                                        } else {
                                                NULL
                                        },
                                
                                output_dir = output_dir
                        )
                },
                pattern = map(projection_member_index),
                format = "file"
        ),
        
        tar_target(
                projection_ensemble_mean_file,
                {
                        
                        if (model_mode == "individual") {
                                return(NULL)
                        }
                        
                        if (
                                is.null(projection_member_mean_file) ||
                                length(projection_member_mean_file) == 0L
                        ) {
                                return(NULL)
                        }
                        
                        ensemble <- compute_ensemble_summary(
                                raster_files =
                                        projection_member_mean_file
                        )
                        
                        if (is.null(ensemble)) {
                                return(NULL)
                        }
                        
                        path <- file.path(
                                output_dir,
                                paste0(
                                        "ensemble_mean_",
                                        safe_name(scenarios[2]),
                                        "_",
                                        safe_name(model_mode),
                                        ".tif"
                                )
                        )
                        
                        safe_write_raster(
                                ensemble,
                                path
                        )
                },
                format = "file"
        ),
                
        tar_target(
                nc_hist_files,
                selected_member_table$file[
                        tolower(selected_member_table$scenario) ==
                                tolower(scenarios[1])]),
        
        tar_target(
                nc_proj_files,
                selected_member_table$file[
                        tolower(selected_member_table$scenario) ==
                                tolower(scenarios[2])
                ]
        ),
        #====
        # Resolved metric metadata from actual raster products
        #====
        tar_target(hist_metric_metadata,{if (!isTRUE(hist_metric_resolution$available)) {return(NULL)}
                        inspect_metric_raster(raster_file = t_raster_hist_file,
                                dataset = "historical",
                                metric = variables[1],
                                source_variable = if (is.null(hist_metric_resolution$source_variable)) {
                                                "derived"} else {
                                                hist_metric_resolution$source_variable},
                                source_mode = hist_metric_resolution$mode,
                                definition_id = paste0(
                                        variables[1],
                                        "_NARCliM2_definition"),
                                units_override = NULL)}),
        
        tar_target(
                proj_metric_metadata,
                {
                        
                        if (!isTRUE(proj_metric_resolution$available)) {
                                return(NULL)
                        }
                        
                        inspect_metric_raster(
                                raster_file = t_raster_proj_file,
                                dataset = "NARCliM2",
                                metric = variables[1],
                                
                                source_variable =
                                        if (is.null(proj_metric_resolution$source_variable)) {
                                                "derived"
                                        } else {
                                                proj_metric_resolution$source_variable
                                        },
                                
                                source_mode =
                                        proj_metric_resolution$mode,
                                
                                definition_id = paste0(
                                        variables[1],
                                        "_NARCliM2_definition"
                                ),
                                
                                units_override = NULL
                        )
                }
        ),
        
        #====
        # Metric planning / provenance
        #====
        
        tar_target(
                metric_derivation_requirements,
                list(
                        FFDI = c(
                                "tasmax",
                                "relative_humidity",
                                "wind_speed",
                                "drought_factor"
                        )
                )
        ),
        
        tar_target(
                hist_metric_resolution,
                resolve_metric_source(
                        requested_metric = variables[1],
                        available_variables =
                                extract_available_variable_names(
                                        nc_hist_files
                                ),
                        direct_metric_variables = variables[1],
                        derivation_requirements =
                                metric_derivation_requirements
                )
        ),
        
        tar_target(
                proj_metric_resolution,
                resolve_metric_source(
                        requested_metric = variables[1],
                        available_variables =
                                extract_available_variable_names(
                                        nc_proj_files
                                ),
                        direct_metric_variables = variables[1],
                        derivation_requirements =
                                metric_derivation_requirements
                )
        ),
        
        tar_target(
                metric_comparison_plan,
                {
                        
                        decide_metric_comparison(
                                hist_resolution =
                                        hist_metric_resolution,
                                
                                proj_resolution =
                                        proj_metric_resolution,
                                
                                hist_meta =
                                        hist_metric_metadata,
                                
                                proj_meta =
                                        proj_metric_metadata
                        )
                }
        ),
        
        tar_target(city_shape, {
                if (!"shapes" %in% selected_targets) return(NULL)
                get_abs_geography(gpkg_path, abs_filter_column, abs_filter_values[1])
        }),
        
        tar_target(t_raster_hist_file, {
                if (!"total" %in% selected_targets) return(NULL)
                r <- load_and_crop_rasters(nc_hist_files, variables[1], hist_start_year, hist_end_year,
                                           test_mode, if ("shapes" %in% selected_targets) city_shape else NULL,
                                           if (isTRUE(test_mode)) sample_layers else NULL)
                if (is.null(r)) return(NULL)
                path <- file.path(output_dir, paste0("hist_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
                safe_write_raster(r, path)
        }, format = "file"),
        
        tar_target(t_raster_proj_file, {
                if (!"total" %in% selected_targets) return(NULL)
                r <- load_and_crop_rasters(nc_proj_files, variables[1], proj_start_year, proj_end_year,
                                           test_mode, if ("shapes" %in% selected_targets) city_shape else NULL,
                                           if (isTRUE(test_mode)) sample_layers else NULL)
                if (is.null(r)) return(NULL)
                path <- file.path(output_dir, paste0("proj_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
                safe_write_raster(r, path)
        }, format = "file"),
        
        # Step 2: Aggregations with defensive checks
        tar_target(raster_hist_mean_file, {
                if (!"mean" %in% selected_targets) return(NULL)
                if (is.null(t_raster_hist_file) || length(t_raster_hist_file) == 0 || !file.exists(t_raster_hist_file)) return(NULL)
                r <- terra::rast(t_raster_hist_file)
                if (nlyr(r) == 0) return(NULL)
                r_mean <- deep_mean(r)
                if (is.null(r_mean)) return(NULL)
                path <- file.path(output_dir, paste0("hist_mean_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
                safe_write_raster(r_mean, path)
        }, format = "file"),
        
        tar_target(raster_proj_mean_file, {
                if (!"mean" %in% selected_targets) return(NULL)
                if (is.null(t_raster_proj_file) || length(t_raster_proj_file) == 0 || !file.exists(t_raster_proj_file)) return(NULL)
                r <- terra::rast(t_raster_proj_file)
                if (nlyr(r) == 0) return(NULL)
                r_mean <- deep_mean(r)
                if (is.null(r_mean)) return(NULL)
                path <- file.path(output_dir, paste0("proj_mean_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
                safe_write_raster(r_mean, path)
        }, format = "file"),
        
        tar_target(
                raster_delta_mean_file,
                {
                        
                        # User did not request a historical-projection mean delta.
                        if (!"delta_mean" %in% selected_targets) {
                                return(NULL)
                        }
                        
                        # Scientific compatibility gate.
                        if (
                                is.null(metric_comparison_plan) ||
                                !isTRUE(metric_comparison_plan$delta_allowed)
                        ) {
                                
                                message(
                                        "Skipping mean delta: historical/projection comparison ",
                                        "is not scientifically permitted. ",
                                        if (!is.null(metric_comparison_plan$reason)) {
                                                paste0(
                                                        "Reason: ",
                                                        metric_comparison_plan$reason
                                                )
                                        } else {
                                                ""
                                        }
                                )
                                
                                return(NULL)
                        }
                        
                        # Required products must exist.
                        if (
                                is.null(raster_hist_mean_file) ||
                                length(raster_hist_mean_file) == 0 ||
                                !file.exists(raster_hist_mean_file)
                        ) {
                                return(NULL)
                        }
                        
                        if (
                                is.null(raster_proj_mean_file) ||
                                length(raster_proj_mean_file) == 0 ||
                                !file.exists(raster_proj_mean_file)
                        ) {
                                return(NULL)
                        }
                        
                        delta <- compute_delta_raster(
                                raster_proj_mean_file,
                                raster_hist_mean_file
                        )
                        
                        if (is.null(delta)) {
                                return(NULL)
                        }
                        
                        path <- file.path(
                                output_dir,
                                paste0(
                                        "delta_mean_",
                                        format(
                                                Sys.time(),
                                                "%Y%m%d_%H%M%S"
                                        ),
                                        ".tif"
                                )
                        )
                        
                        safe_write_raster(
                                delta,
                                path
                        )
                },
                format = "file"
        ),
        
        tar_target(raster_monthly_hist_file, {
                if (!"monthly" %in% selected_targets) return(NULL)
                if (is.null(t_raster_hist_file) || length(t_raster_hist_file) == 0 || !file.exists(t_raster_hist_file)) return(NULL)
                safe_aggregate_raster(t_raster_hist_file, hist_start_year, hist_end_year, "monthly", output_dir, "monthly_hist")
        }, format = "file"),
        
        tar_target(raster_monthly_proj_file, {
                if (!"monthly" %in% selected_targets) return(NULL)
                if (is.null(t_raster_proj_file) || length(t_raster_proj_file) == 0 || !file.exists(t_raster_proj_file)) return(NULL)
                safe_aggregate_raster(t_raster_proj_file, proj_start_year, proj_end_year, "monthly", output_dir, "monthly_proj")
        }, format = "file"),
        
        tar_target(raster_annual_hist_file, {
                if (!"annual" %in% selected_targets) return(NULL)
                if (is.null(t_raster_hist_file) || length(t_raster_hist_file) == 0 || !file.exists(t_raster_hist_file)) return(NULL)
                safe_aggregate_raster(t_raster_hist_file, hist_start_year, hist_end_year, "annual", output_dir, "annual_hist")
        }, format = "file"),
        
        tar_target(raster_annual_proj_file, {
                if (!"annual" %in% selected_targets) return(NULL)
                if (is.null(t_raster_proj_file) || length(t_raster_proj_file) == 0 || !file.exists(t_raster_proj_file)) return(NULL)
                safe_aggregate_raster(t_raster_proj_file, proj_start_year, proj_end_year, "annual", output_dir, "annual_proj")
        }, format = "file"),
        
        # Percentile rasters
        tar_target(raster_percentile_hist_file, {
                if (!"percentile" %in% selected_targets) return(NULL)
                if (is.null(t_raster_hist_file) || length(t_raster_hist_file) == 0 || !file.exists(t_raster_hist_file)) return(NULL)
                r <- terra::rast(t_raster_hist_file)
                if (nlyr(r) == 0) return(NULL)
                perc <- percentile_map(r, Percentile / 100)
                if (is.null(perc)) return(NULL)
                path <- file.path(output_dir, paste0("hist_percentile_", Percentile, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
                safe_write_raster(perc, path)
        }, format = "file"),
        
        tar_target(raster_percentile_proj_file, {
                if (!"percentile" %in% selected_targets) return(NULL)
                if (is.null(t_raster_proj_file) || length(t_raster_proj_file) == 0 || !file.exists(t_raster_proj_file)) return(NULL)
                r <- terra::rast(t_raster_proj_file)
                if (nlyr(r) == 0) return(NULL)
                perc <- percentile_map(r, Percentile / 100)
                if (is.null(perc)) return(NULL)
                path <- file.path(output_dir, paste0("proj_percentile_", Percentile, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
                safe_write_raster(perc, path)
        }, format = "file"),
        
        #====
        # Threshold Exceedance Rasters
        #====
        tar_target(
                raster_threshold_hist_file,
                {
                        
                        if (!"threshold" %in% selected_targets) {
                                return(NULL)
                        }
                        
                        if (
                                is.null(t_raster_hist_file) ||
                                length(t_raster_hist_file) == 0 ||
                                !file.exists(t_raster_hist_file)
                        ) {
                                return(NULL)
                        }
                        
                        r <- terra::rast(
                                t_raster_hist_file
                        )
                        
                        if (
                                terra::nlyr(r) == 0 ||
                                !terra::hasValues(r)
                        ) {
                                return(NULL)
                        }
                        
                        threshold_r <- compute_threshold_days(
                                r = r,
                                threshold = threshold_value
                        )
                        
                        if (is.null(threshold_r)) {
                                return(NULL)
                        }
                        
                        path <- file.path(
                                output_dir,
                                paste0(
                                        "hist_threshold_",
                                        safe_name(threshold_value),
                                        "_",
                                        format(
                                                Sys.time(),
                                                "%Y%m%d_%H%M%S"
                                        ),
                                        ".tif"
                                )
                        )
                        
                        safe_write_raster(
                                threshold_r,
                                path
                        )
                },
                format = "file"
        ),
        
        tar_target(
                raster_threshold_proj_file,
                {
                        
                        if (!"threshold" %in% selected_targets) {
                                return(NULL)
                        }
                        
                        if (
                                is.null(t_raster_proj_file) ||
                                length(t_raster_proj_file) == 0 ||
                                !file.exists(t_raster_proj_file)
                        ) {
                                return(NULL)
                        }
                        
                        r <- terra::rast(
                                t_raster_proj_file
                        )
                        
                        if (
                                terra::nlyr(r) == 0 ||
                                !terra::hasValues(r)
                        ) {
                                return(NULL)
                        }
                        
                        threshold_r <- compute_threshold_days(
                                r = r,
                                threshold = threshold_value
                        )
                        
                        if (is.null(threshold_r)) {
                                return(NULL)
                        }
                        
                        path <- file.path(
                                output_dir,
                                paste0(
                                        "proj_threshold_",
                                        safe_name(threshold_value),
                                        "_",
                                        format(
                                                Sys.time(),
                                                "%Y%m%d_%H%M%S"
                                        ),
                                        ".tif"
                                )
                        )
                        
                        safe_write_raster(
                                threshold_r,
                                path
                        )
                },
                format = "file"
        ),
        
        ###############################################################################
        # Consecutive threshold-event statistics
        ###############################################################################
        
        tar_target(
                raster_event_metrics_hist_file,
                {
                        
                        if (!"events" %in% selected_targets) {
                                return(NULL)
                        }
                        
                        if (
                                is.null(t_raster_hist_file) ||
                                length(t_raster_hist_file) == 0 ||
                                !file.exists(t_raster_hist_file)
                        ) {
                                return(NULL)
                        }
                        
                        r <- terra::rast(
                                t_raster_hist_file
                        )
                        
                        if (
                                terra::nlyr(r) == 0 ||
                                !terra::hasValues(r)
                        ) {
                                return(NULL)
                        }
                        
                        metrics <- compute_consecutive_event_metrics(
                                r = r,
                                threshold = threshold_value,
                                min_duration = min_duration
                        )
                        
                        if (is.null(metrics)) {
                                return(NULL)
                        }
                        
                        path <- file.path(
                                output_dir,
                                paste0(
                                        "hist_event_metrics_",
                                        safe_name(threshold_value),
                                        "_dur",
                                        min_duration,
                                        "_",
                                        format(
                                                Sys.time(),
                                                "%Y%m%d_%H%M%S"
                                        ),
                                        ".tif"
                                )
                        )
                        
                        safe_write_raster(
                                metrics,
                                path
                        )
                },
                format = "file"
        ),
        
        tar_target(
                raster_event_metrics_proj_file,
                {
                        
                        if (!"events" %in% selected_targets) {
                                return(NULL)
                        }
                        
                        if (
                                is.null(t_raster_proj_file) ||
                                length(t_raster_proj_file) == 0 ||
                                !file.exists(t_raster_proj_file)
                        ) {
                                return(NULL)
                        }
                        
                        r <- terra::rast(
                                t_raster_proj_file
                        )
                        
                        if (
                                terra::nlyr(r) == 0 ||
                                !terra::hasValues(r)
                        ) {
                                return(NULL)
                        }
                        
                        metrics <- compute_consecutive_event_metrics(
                                r = r,
                                threshold = threshold_value,
                                min_duration = min_duration
                        )
                        
                        if (is.null(metrics)) {
                                return(NULL)
                        }
                        
                        path <- file.path(
                                output_dir,
                                paste0(
                                        "proj_event_metrics_",
                                        safe_name(threshold_value),
                                        "_dur",
                                        min_duration,
                                        "_",
                                        format(
                                                Sys.time(),
                                                "%Y%m%d_%H%M%S"
                                        ),
                                        ".tif"
                                )
                        )
                        
                        safe_write_raster(
                                metrics,
                                path
                        )
                },
                format = "file"
        ),
        #==========================================================
        ## Plotting Targets
        #Percentile plots
        tar_target(plot_percentile_hist, {
                if (!"plots" %in% selected_targets|| !"percentile" %in% selected_targets) 
                        return(NULL)
                if (is.null(raster_percentile_hist_file) || 
                    !file.exists(raster_percentile_hist_file)) return(NULL)
                if(file.info(raster_percentile_hist_file)$size==0){
                        message("Skipping percentile plot: raster file is empty.")
                        return(NULL)}
                perc_df <- to_df(terra::rast(raster_percentile_hist_file))
                plot_extent <- as.vector(terra::ext(terra::rast(raster_percentile_hist_file)))
                plot_percentile_map(perc_df, plot_extent, variables[1], Percentile,
                                    paste(hist_start_year, hist_end_year), scenario = scenarios[1], output_dir)
        }),
        
        tar_target(plot_percentile_proj, {
                if (!"plots" %in% selected_targets || !"percentile" %in% selected_targets)
                        return(NULL)
                if (is.null(raster_percentile_proj_file) || 
                    !file.exists(raster_percentile_proj_file)) return(NULL)
                if(file.info(raster_percentile_proj_file)$size==0){
                        message("Skipping percentile plot: raster file is empty.")
                        return(NULL)}
                perc_df <- to_df(terra::rast(raster_percentile_proj_file))
                plot_extent <- as.vector(terra::ext(terra::rast(raster_percentile_proj_file)))
                plot_percentile_map(perc_df, plot_extent, variables[1], Percentile,
                                    paste(proj_start_year, proj_end_year), scenario = scenarios[2], output_dir)
        }),
        
        # Polygon summaries
        tar_target(polygon_summary_hist, {
                if (!"shapes" %in% selected_targets) return(NULL)
                if (is.null(city_shape) || is.null(raster_hist_mean_file)) return(NULL)
                r <- terra::rast(raster_hist_mean_file)
                aggregate_by_polygon(r, city_shape, paste0(variables[1], "_hist_mean"), output_dir)
        }),
        
        tar_target(polygon_summary_proj, {
                if (!"shapes" %in% selected_targets) return(NULL)
                if (is.null(city_shape) || is.null(raster_proj_mean_file)) return(NULL)
                r <- terra::rast(raster_proj_mean_file)
                aggregate_by_polygon(r, city_shape, paste0(variables[1], "_proj_mean"), output_dir)
        }),
        
        
        # Step 3: Plots
        tar_target(
                plot_mean_comparison,
                {
                        
                        if (!"plots" %in% selected_targets) {
                                return(NULL)
                        }
                        if (
                                is.null(metric_comparison_plan) ||
                                !isTRUE(metric_comparison_plan$comparison_allowed)
                        ) {
                                
                                message(
                                        "Skipping historical-projection comparison plot: ",
                                        "comparison is not scientifically permitted. ",
                                        if (!is.null(metric_comparison_plan$reason)) {
                                                paste0(
                                                        "Reason: ",
                                                        metric_comparison_plan$reason
                                                )
                                        } else {
                                                ""
                                        }
                                )
                                
                                return(NULL)
                        }
                        
                        
                        if (
                                is.null(raster_hist_mean_file) ||
                                is.null(raster_proj_mean_file) ||
                                is.null(raster_delta_mean_file)
                        ) {
                                return(NULL)
                        }
                        
                        if (
                                !file.exists(raster_hist_mean_file) ||
                                !file.exists(raster_proj_mean_file) ||
                                !file.exists(raster_delta_mean_file)
                        ) {
                                return(NULL)
                        }
                        
                        # Load the already-computed Targets products.
                        hist_r <- terra::rast(
                                raster_hist_mean_file
                        )
                        
                        proj_r <- terra::rast(
                                raster_proj_mean_file
                        )
                        
                        delta_r <- terra::rast(
                                raster_delta_mean_file
                        )
                        
                        # Defensive checks.
                        if (
                                !terra::hasValues(hist_r) ||
                                !terra::hasValues(proj_r) ||
                                !terra::hasValues(delta_r)
                        ) {
                                stop(
                                        "One or more mean-comparison rasters have no values."
                                )
                        }
                        
                        # These should normally already match because the delta target
                        # has completed successfully.
                        if (
                                !terra::compareGeom(
                                        hist_r,
                                        proj_r,
                                        stopOnError = FALSE
                                )
                        ) {
                                hist_r <- align_and_resample(
                                        ref = proj_r,
                                        target = hist_r
                                )
                        }
                        
                        # Convert existing validated raster products to plotting data.
                        hist_df <- to_df(hist_r)
                        proj_df <- to_df(proj_r)
                        delta_df <- to_df(delta_r)
                        
                        if (
                                nrow(hist_df) == 0 ||
                                nrow(proj_df) == 0 ||
                                nrow(delta_df) == 0
                        ) {
                                stop(
                                        "Raster-to-data-frame conversion returned an empty result."
                                )
                        }
                        
                        plot_extent <- c(
                                terra::xmin(hist_r),
                                terra::xmax(hist_r),
                                terra::ymin(hist_r),
                                terra::ymax(hist_r)
                        )
                        
                        plot_raster_comparison(
                                hist_df,
                                proj_df,
                                delta_df,
                                plot_extent,
                                variables[1],
                                paste(
                                        hist_start_year,
                                        hist_end_year
                                ),
                                paste(
                                        proj_start_year,
                                        proj_end_year
                                ),
                                output_dir
                        )
                }
        ),
        
        tar_target(plot_monthly_hist, {
                if (!"plots" %in% selected_targets) return(NULL)
                if (is.null(raster_monthly_hist_file) || !file.exists(raster_monthly_hist_file)) return(NULL)
                
                r <- terra::rast(raster_monthly_hist_file)
                names(r) <- paste0("monthly_", seq_len(nlyr(r)))  # Force correct names
                monthly_df <- to_df(r)
                monthly_df$month <- factor(monthly_df$layer, levels = paste0("monthly_", 1:12), labels = month.abb)
                plot_extent <- as.vector(terra::ext(r))
                plot_monthly_climatology(monthly_df, plot_extent, variables[1],
                                         paste(hist_start_year, hist_end_year), scenario = scenarios[1], output_dir)
        }),
        
        tar_target(plot_monthly_proj, {
                if (!"plots" %in% selected_targets) return(NULL)
                if (is.null(raster_monthly_proj_file) || !file.exists(raster_monthly_proj_file)) return(NULL)
                
                r <- terra::rast(raster_monthly_proj_file)
                names(r) <- paste0("monthly_", seq_len(nlyr(r)))
                monthly_df <- to_df(r)
                monthly_df$month <- factor(monthly_df$layer, levels = paste0("monthly_", 1:12), labels = month.abb)
                plot_extent <- as.vector(terra::ext(r))
                plot_monthly_climatology(monthly_df, plot_extent, variables[1],
                                         paste(proj_start_year, proj_end_year), scenario = scenarios[2], output_dir)
        }),
        
        tar_target(
                plot_annual_hist,
                {
                        
                        if (
                                !"plots" %in% selected_targets ||
                                !"annual" %in% selected_targets
                        ) {
                                return(NULL)
                        }
                        
                        if (
                                is.null(raster_annual_hist_file) ||
                                !file.exists(raster_annual_hist_file)
                        ) {
                                return(NULL)
                        }
                        
                        r <- terra::rast(
                                raster_annual_hist_file
                        )
                        
                        if (
                                terra::nlyr(r) == 0 ||
                                !terra::hasValues(r)
                        ) {
                                return(NULL)
                        }
                        
                        # Preserve calendar-aware names such as annual_1990.
                        annual_df <- to_df(r)
                        
                        if (nrow(annual_df) == 0) {
                                return(NULL)
                        }
                        
                        plot_extent <- c(
                                terra::xmin(r),
                                terra::xmax(r),
                                terra::ymin(r),
                                terra::ymax(r)
                        )
                        
                        plot_annual_or_monthly(
                                r_df = annual_df,
                                plot_extent = plot_extent,
                                variable = variables[1],
                                period = "annual",
                                range_label = paste(
                                        hist_start_year,
                                        hist_end_year
                                ),
                                scenario = scenarios[1],
                                output_dir = output_dir
                        )
                }
        ),
        
        tar_target(
                plot_annual_proj,
                {
                        
                        if (
                                !"plots" %in% selected_targets ||
                                !"annual" %in% selected_targets
                        ) {
                                return(NULL)
                        }
                        
                        if (
                                is.null(raster_annual_proj_file) ||
                                !file.exists(raster_annual_proj_file)
                        ) {
                                return(NULL)
                        }
                        
                        r <- terra::rast(
                                raster_annual_proj_file
                        )
                        
                        if (
                                terra::nlyr(r) == 0 ||
                                !terra::hasValues(r)
                        ) {
                                return(NULL)
                        }
                        
                        # Preserve calendar-aware names such as annual_2020.
                        annual_df <- to_df(r)
                        
                        if (nrow(annual_df) == 0) {
                                return(NULL)
                        }
                        
                        plot_extent <- c(
                                terra::xmin(r),
                                terra::xmax(r),
                                terra::ymin(r),
                                terra::ymax(r)
                        )
                        
                        plot_annual_or_monthly(
                                r_df = annual_df,
                                plot_extent = plot_extent,
                                variable = variables[1],
                                period = "annual",
                                range_label = paste(
                                        proj_start_year,
                                        proj_end_year
                                ),
                                scenario = scenarios[2],
                                output_dir = output_dir
                        )
                }
        )
)


