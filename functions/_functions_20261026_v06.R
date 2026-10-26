# ----------------------
# Raster Loading & Cropping
# ----------------------
filter_nc_files <- function(
                data_dir,
                variables,
                scenarios,
                models = NULL
) {
        
        all_files <- list.files(
                data_dir,
                pattern = "\\.nc$",
                full.names = TRUE,
                recursive = TRUE
        )
        
        if (length(all_files) == 0) {
                return(character(0))
        }
        
        x <- tolower(all_files)
        
        # Variable filter
        if (!is.null(variables) && length(variables) > 0) {
                
                pat_var <- paste(
                        tolower(variables),
                        collapse = "|"
                )
                
                keep_var <- grepl(
                        pat_var,
                        x
                )
                
                all_files <- all_files[keep_var]
                x <- x[keep_var]
        }
        
        # Scenario filter
        if (!is.null(scenarios) && length(scenarios) > 0) {
                
                pat_scen <- paste(
                        tolower(scenarios),
                        collapse = "|"
                )
                
                keep_scen <- grepl(
                        pat_scen,
                        x
                )
                
                all_files <- all_files[keep_scen]
                x <- x[keep_scen]
        }
        
        # GCM/model filter
        if (!is.null(models) && length(models) > 0) {
                
                pat_model <- paste(
                        tolower(models),
                        collapse = "|"
                )
                
                keep_model <- grepl(
                        pat_model,
                        x
                )
                
                all_files <- all_files[keep_model]
        }
        
        sort(all_files)
}

# ----------------------
# Polygon aggregation
# ----------------------

aggregate_by_polygon <- function(r, polygons, metric_name, output_dir) {
        if (is.null(r) || is.null(polygons)) return(NULL)
        
        # Force CRS match
        polygons <- sf::st_transform(polygons, terra::crs(r))
        
        extracted <- terra::extract(r, terra::vect(polygons), fun = mean, na.rm = TRUE)
        df_ex <- as.data.frame(extracted)
        if ("ID" %in% names(df_ex)) df_ex$ID <- NULL
        
        # Keep polygons as sf object
        tab <- polygons
        tab <- cbind(tab, df_ex)
        
        csv_path <- file.path(output_dir, paste0("polygon_summary_", safe_name(metric_name), ".csv"))
        tryCatch(sf::st_write(tab, csv_path, layer_options = "GEOMETRY=AS_WKT", delete_dsn = TRUE, quiet = TRUE),
                 error = function(e) write.csv(as.data.frame(tab), csv_path, row.names = FALSE))
        
        # Plot
        map_path <- NULL
        numcols <- names(tab)[sapply(tab, is.numeric)]
        if (length(numcols) > 0) {
                fillc <- numcols[1]
                p <- ggplot2::ggplot(tab) +
                        ggplot2::geom_sf(aes(fill = !!sym(fillc)), color = "black") +
                        safe_scale_fill_viridis(na.value = "grey80") +
                        ggplot2::labs(title = paste("Mean", metric_name)) +
                        ggplot2::theme_minimal()
                map_path <- file.path(output_dir, paste0("polygon_map_", safe_name(metric_name), ".png"))
                ggplot2::ggsave(map_path, p, width = 10, height = 8)
        }
        list(table = tab, csv = csv_path, map = map_path)
}

#---------------------------
# Safe name 
#---------------------------
safe_name <- function(x) {
        x <- as.character(x)
        x <- gsub("[[:space:]/\\\\]+", "_", x)
        x <- gsub("[^A-Za-z0-9_.:-]", "", x)
        x
}
# ----------------------
# Raster Loading & Cropping
# Load, time-subset and spatially crop climate rasters
#
# IMPORTANT:
# User-selected dates are applied BEFORE downstream statistics.
#
# Processing order:
#   identify files
#   -> open
#   -> retain requested dates
#   -> combine
#   -> optional test subset
#   -> crop/mask
#----------------------

load_and_crop_rasters <- function(
                files,
                var_token,
                start_year,
                end_year,
                test_mode = FALSE,
                city_shape = NULL,
                sample_layers = NULL,
                threshold_var = NULL
) {
        
        # -------------------------------------------------------------------------
        # Basic input checks
        # -------------------------------------------------------------------------
        
        if (
                is.null(files) ||
                length(files) == 0
        ) {
                return(NULL)
        }
        
        # -------------------------------------------------------------------------
        # Select requested variable
        # -------------------------------------------------------------------------
        
        keep_var <- grepl(
                tolower(var_token),
                tolower(basename(files)),
                fixed = TRUE
        )
        
        files_var <- files[keep_var]
        
        if (!is.null(threshold_var)) {
                
                keep_threshold <- grepl(
                        tolower(threshold_var),
                        tolower(basename(files_var)),
                        fixed = TRUE
                )
                
                files_var <- files_var[
                        keep_threshold
                ]
        }
        
        if (length(files_var) == 0) {
                return(NULL)
        }
        
        # -------------------------------------------------------------------------
        # Process each file independently.
        #
        # This is still the LOCAL-file implementation.
        # We will later replace the data-access layer with the proven NARCliM
        # streaming implementation without changing downstream functionality.
        # -------------------------------------------------------------------------
        
        raster_list <- vector(
                "list",
                length(files_var)
        )
        
        for (i in seq_along(files_var)) {
                
                f <- files_var[i]
                
                message(
                        "Loading climate file [",
                        i,
                        "/",
                        length(files_var),
                        "]: ",
                        basename(f)
                )
                
                r <- tryCatch(
                        terra::rast(f),
                        error = function(e) {
                                
                                warning(
                                        "Could not open ",
                                        basename(f),
                                        ": ",
                                        e$message
                                )
                                
                                NULL
                        }
                )
                
                if (is.null(r)) {
                        next
                }
                
                # ---------------------------------------------------------------------
                # Use actual NetCDF time coordinate.
                # Do NOT infer years from layer positions.
                # ---------------------------------------------------------------------
                
                time_vals <- tryCatch(
                        terra::time(r),
                        error = function(e) NULL
                )
                
                if (
                        is.null(time_vals) ||
                        length(time_vals) != terra::nlyr(r)
                ) {
                        
                        stop(
                                "No valid time coordinate for file:\n",
                                f,
                                "\nCannot safely apply user-selected date range."
                        )
                }
                
                years <- lubridate::year(
                        time_vals
                )
                
                keep_idx <- which(
                        years >= start_year &
                                years <= end_year
                )
                
                if (length(keep_idx) == 0) {
                        
                        message(
                                "No data in requested period ",
                                start_year,
                                "-",
                                end_year,
                                " for ",
                                basename(f)
                        )
                        
                        next
                }
                
                r <- r[[keep_idx]]
                
                message(
                        "  retained ",
                        terra::nlyr(r),
                        " layers for ",
                        start_year,
                        "-",
                        end_year
                )
                
                raster_list[[i]] <- r
        }
        
        # -------------------------------------------------------------------------
        # Remove files that contributed no requested dates
        # -------------------------------------------------------------------------
        
        raster_list <- raster_list[
                !vapply(
                        raster_list,
                        is.null,
                        logical(1)
                )
        ]
        
        if (length(raster_list) == 0) {
                return(NULL)
        }
        
        # -------------------------------------------------------------------------
        # Combine only the temporally relevant data
        # -------------------------------------------------------------------------
        
        if (length(raster_list) == 1) {
                
                r_all <- raster_list[[1]]
                
        } else {
                
                r_all <- do.call(
                        c,
                        raster_list
                )
        }
        
        # -------------------------------------------------------------------------
        # Ensure chronological order
        # -------------------------------------------------------------------------
        
        time_all <- tryCatch(
                terra::time(r_all),
                error = function(e) NULL
        )
        
        if (
                !is.null(time_all) &&
                length(time_all) == terra::nlyr(r_all)
        ) {
                
                ord <- order(time_all)
                
                if (!identical(
                        ord,
                        seq_along(ord)
                )) {
                        r_all <- r_all[[ord]]
                }
        }
        
        # -------------------------------------------------------------------------
        # Test mode
        #
        # Use a CONTIGUOUS subset rather than random days.
        # Random sampling destroys temporal continuity and would invalidate
        # later heatwave-duration calculations.
        # -------------------------------------------------------------------------
        
        if (
                isTRUE(test_mode) &&
                !is.null(sample_layers) &&
                sample_layers < terra::nlyr(r_all)
        ) {
                
                idx <- seq_len(
                        sample_layers
                )
                
                r_all <- r_all[[idx]]
                
                message(
                        "Test mode: retained first ",
                        sample_layers,
                        " consecutive layers."
                )
        }
        
        # -------------------------------------------------------------------------
        # Spatial crop and mask
        # -------------------------------------------------------------------------
        
        if (!is.null(city_shape)) {
                
                city_shape_use <- city_shape
                
                city_shape_use <- tryCatch(
                        sf::st_transform(
                                city_shape_use,
                                crs = terra::crs(r_all)
                        ),
                        error = function(e) {
                                
                                stop(
                                        "Could not transform selected geography to climate CRS: ",
                                        e$message
                                )
                        }
                )
                
                city_vect <- terra::vect(
                        city_shape_use
                )
                
                r_all <- terra::crop(
                        r_all,
                        city_vect
                )
                
                r_all <- terra::mask(
                        r_all,
                        city_vect
                )
        }
        
        # -------------------------------------------------------------------------
        # Final date-range QA
        # -------------------------------------------------------------------------
        
        final_time <- tryCatch(
                terra::time(r_all),
                error = function(e) NULL
        )
        
        if (
                !is.null(final_time) &&
                length(final_time) > 0
        ) {
                
                final_years <- lubridate::year(
                        final_time
                )
                
                if (
                        min(final_years) < start_year ||
                        max(final_years) > end_year
                ) {
                        
                        stop(
                                "Internal date-range QA failed. ",
                                "Output raster contains dates outside ",
                                start_year,
                                "-",
                                end_year,
                                "."
                        )
                }
        }
        
        r_all
}
# ----------------------
# Safe Write Raster
# ----------------------
safe_write_raster <- function(r, path, overwrite = TRUE) {
        if (is.null(r) || nlyr(r) == 0 || all(is.na(r[]))) {
                message("Raster is empty, not writing: ", path)
                return(NULL)
        }
        terra::writeRaster(r, path, overwrite = overwrite)
        return(path)
}

# ----------------------
# Compute Mean
# ----------------------
deep_mean <- function(r) {
        if (is.null(r) || nlyr(r) == 0) return(NULL)
        terra::app(r, mean, na.rm = TRUE)
}

# ----------------------
# Aggregate Time Period
# ----------------------
aggregate_time_period <- function(
                r,
                start_year,
                end_year,
                period = c("annual", "monthly"),
                fun = mean
) {
        
        if (is.null(r) || terra::nlyr(r) == 0) {
                return(NULL)
        }
        
        period <- match.arg(period)
        
        # ------------------------------------------------------------
        # Read actual NetCDF time coordinate
        # ------------------------------------------------------------
        
        time_vals <- tryCatch(
                terra::time(r),
                error = function(e) NULL
        )
        
        if (
                is.null(time_vals) ||
                length(time_vals) != terra::nlyr(r)
        ) {
                stop(
                        "Raster has no valid time coordinate. ",
                        "Annual/monthly aggregation will not guess dates from layer position."
                )
        }
        
        dates <- as.Date(time_vals)
        
        if (anyNA(dates)) {
                stop(
                        "One or more raster time values could not be converted to Date."
                )
        }
        
        years  <- lubridate::year(dates)
        months <- lubridate::month(dates)
        
        # ------------------------------------------------------------
        # Restrict to user-selected years
        # ------------------------------------------------------------
        
        keep <- which(
                years >= start_year &
                        years <= end_year
        )
        
        if (length(keep) == 0) {
                return(NULL)
        }
        
        r <- r[[keep]]
        dates  <- dates[keep]
        years  <- years[keep]
        months <- months[keep]
        
        # ------------------------------------------------------------
        # Calendar grouping
        # ------------------------------------------------------------
        
        if (period == "annual") {
                
                groups <- as.character(years)
                
        } else {
                
                # Existing Climate Menu meaning:
                # monthly = 12-layer monthly climatology
                groups <- sprintf(
                        "%02d",
                        months
                )
        }
        
        grouped <- split(
                seq_len(terra::nlyr(r)),
                groups
        )
        
        # ------------------------------------------------------------
        # Aggregate each calendar group
        # ------------------------------------------------------------
        
        out_list <- lapply(
                grouped,
                function(idx) {
                        
                        terra::app(
                                r[[idx]],
                                fun = fun,
                                na.rm = TRUE
                        )
                }
        )
        
        out <- do.call(
                c,
                out_list
        )
        
        if (!inherits(out, "SpatRaster")) {
                out <- terra::rast(out_list)
        }
        
        # ------------------------------------------------------------
        # Stable layer names
        # ------------------------------------------------------------
        
        if (period == "annual") {
                
                names(out) <- paste0(
                        "annual_",
                        names(grouped)
                )
                
        } else {
                
                names(out) <- paste0(
                        "monthly_",
                        as.integer(names(grouped))
                )
        }
        
        out
}
# ----------------------
# Safe Aggregate Wrapper
# ----------------------
safe_aggregate_raster <- function(raster_file, start_year, end_year, period, output_dir, prefix) {
        r <- terra::rast(raster_file)
        if (nlyr(r) == 0) {
                message("No layers in raster, skipping aggregation.")
                return(NULL)
        }
        
        agg_list <- aggregate_time_period(r, start_year, end_year, period)
        if (is.null(agg_list) || nlyr(agg_list) == 0) {
                message("No layers after aggregation, skipping raster.")
                return(NULL)
        }
        
        outfile <- file.path(output_dir, paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif"))
        safe_write_raster(agg_list, outfile, overwrite = TRUE)
        return(outfile)
}

###############################################################################
# Raster geometry alignment
###############################################################################

align_and_resample <- function(
                ref,
                target,
                method = "bilinear",
                tol = 1e-6
) {
        
        if (is.null(ref) || is.null(target)) {
                return(NULL)
        }
        
        if (terra::crs(ref) != terra::crs(target)) {
                
                target <- terra::project(
                        target,
                        terra::crs(ref)
                )
        }
        
        ref_ext <- c(
                terra::xmin(ref),
                terra::xmax(ref),
                terra::ymin(ref),
                terra::ymax(ref)
        )
        
        target_ext <- c(
                terra::xmin(target),
                terra::xmax(target),
                terra::ymin(target),
                terra::ymax(target)
        )
        
        ext_diff <- max(
                abs(ref_ext - target_ext)
        )
        
        res_diff <- max(
                abs(
                        terra::res(ref) -
                                terra::res(target)
                )
        )
        
        if (
                ext_diff > tol ||
                res_diff > tol
        ) {
                
                target <- terra::resample(
                        target,
                        ref,
                        method = method
                )
        }
        
        target
}


# ----------------------
# Compute Delta Raster
# ----------------------
compute_delta_raster <- function(rast_proj, rast_hist) {
        if (is.null(rast_proj) || is.null(rast_hist)) return(NULL)
        rast_proj <- terra::rast(rast_proj)
        rast_hist <- terra::rast(rast_hist)
        if (!terra::compareGeom(rast_proj, rast_hist, stopOnError = FALSE)) {
                rast_hist <- terra::resample(rast_hist, rast_proj)
        }
        delta <- rast_proj - rast_hist
        return(delta)
}

# ----------------------
# Convert Raster to Data Frame
# ----------------------
to_df <- function(r) {
        
        if (is.null(r)) {
                return(data.frame())
        }
        
        if (
                !inherits(r, "SpatRaster") ||
                terra::nlyr(r) == 0 ||
                !terra::hasValues(r)
        ) {
                return(data.frame())
        }
        
        df <- terra::as.data.frame(
                r,
                xy = TRUE,
                na.rm = FALSE
        )
        
        if (
                is.null(df) ||
                nrow(df) == 0 ||
                !all(c("x", "y") %in% names(df))
        ) {
                return(data.frame())
        }
        
        value_cols <- setdiff(
                names(df),
                c("x", "y")
        )
        
        if (length(value_cols) == 0) {
                return(data.frame())
        }
        
        # Keep only numeric raster-value columns.
        value_cols <- value_cols[
                vapply(
                        df[value_cols],
                        is.numeric,
                        logical(1)
                )
        ]
        
        if (length(value_cols) == 0) {
                return(data.frame())
        }
        
        if (length(value_cols) == 1) {
                
                out <- data.frame(
                        x = df$x,
                        y = df$y,
                        layer = value_cols[1],
                        value = as.numeric(
                                df[[value_cols[1]]]
                        )
                )
                
        } else {
                
                out <- tidyr::pivot_longer(
                        df[
                                ,
                                c(
                                        "x",
                                        "y",
                                        value_cols
                                ),
                                drop = FALSE
                        ],
                        cols = tidyselect::all_of(
                                value_cols
                        ),
                        names_to = "layer",
                        values_to = "value"
                )
                
                out$value <- as.numeric(
                        out$value
                )
        }
        
        out
}

#-------------------------
# get shape files
#-------------------------
get_abs_geography <- function(path, filter_column = NULL, filter_value = NULL) {
        if (is.null(path) || !file.exists(path)) {
                warning("Geometry file missing: ", path)
                return(NULL)
        }
        sf_obj <- tryCatch(sf::st_read(path, quiet = TRUE), error = function(e) { warning(e$message); NULL })
        if (is.null(sf_obj)) return(NULL)
        if (!is.null(filter_column) && !is.null(filter_value)) {
                if (!(filter_column %in% names(sf_obj))) {
                        warning("Filter column not found: ", filter_column)
                        return(NULL)
                }
                sf_obj <- sf_obj[sf_obj[[filter_column]] %in% filter_value, , drop = FALSE]
                if (nrow(sf_obj) == 0) return(NULL)
        }
        sf_obj
}

safe_scale_fill_viridis <- function(...) {
        ggplot2::scale_fill_viridis_c(...)
}

# ----------------------
# Plot Mean Comparison
# ----------------------
# plot_raster_comparison <- function(hist_df, proj_df = NULL, delta_df = NULL, plot_extent, variable, hist_range, proj_range, output_dir) {
#         if (is.null(plot_extent) || length(plot_extent) < 4) plot_extent <- c(-180,180,-90,90)
#         vals <- c(if (nrow(hist_df)>0) hist_df$value else NA,
#                   if (!is.null(proj_df) && nrow(proj_df)>0) proj_df$value else NA,
#                   if (!is.null(delta_df) && nrow(delta_df)>0) delta_df$value else NA)
#         vals <- vals[!is.na(vals)]
#         if (length(vals) == 0) return(NULL)
#         minv <- min(vals); maxv <- max(vals); mid <- (minv+maxv)/2
#         timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
#         
#         p_hist <- ggplot2::ggplot(hist_df, ggplot2::aes(x=x,y=y,fill=value)) +
#                 ggplot2::geom_raster() +
#                 ggplot2::coord_quickmap(xlim = plot_extent[1:2], ylim = plot_extent[3:4]) +
#                 ggplot2::scale_fill_gradient2(limits=c(minv,maxv), low="blue", mid="white", high="red", midpoint=mid, na.value="palegoldenrod") +
#                 ggplot2::labs(title=paste(variable,"Historical",hist_range)) +
#                 ggplot2::theme_minimal(base_size=14)+
#                 ggplot2::theme(panel.background = element_rect(fill="white"),
#                                plot.background = element_rect(fill="white"))
#         f <- file.path(output_dir, paste0(variable,"_comparison_",timestamp,".png"))
#         ggplot2::ggsave(f, p_hist, width=8, height=6)
#         return(f)
# }

# 
# 
# 
# 
# plot_monthly_climatology <- function(monthly_df, plot_extent, variable, range_label, output_dir) {
#         if (nrow(monthly_df) == 0) return(NULL)
#         monthly_df$month <- factor(monthly_df$layer, levels = paste0("monthly_",1:12), labels = month.abb)
#         p <- ggplot2::ggplot(monthly_df, ggplot2::aes(x=x,y=y,fill=value)) +
#                 ggplot2::geom_raster() +
#                 ggplot2::coord_quickmap(xlim=plot_extent[1:2], ylim=plot_extent[3:4]) +
#                 ggplot2::facet_wrap(~month, ncol=4) +
#                 ggplot2::labs(title=paste(variable,"Monthly Climatology",range_label)) +
#                 ggplot2::theme_minimal(base_size=14)+
#                 ggplot2::theme(panel.background = element_rect(fill="white"),
#                                plot.background = element_rect(fill="white"))
#         f <- file.path(output_dir, paste0(safe_name(variable),"_monthly_climatology_",format(Sys.time(),"%Y%m%d_%H%M%S"),".png"))
#         ggplot2::ggsave(f, p, width=12, height=8)
#         f
# }
# 
# plot_annual_or_monthly <- function(r_df, plot_extent, variable, period = c("annual", "monthly"), output_dir) {
#         period <- match.arg(period)
#         if (nrow(r_df) == 0) return(NULL)
#         r_df$period <- factor(r_df$layer, levels = unique(r_df$layer))
#         p <- ggplot2::ggplot(r_df, ggplot2::aes(x = x, y = y, fill = value)) +
#                 ggplot2::geom_raster() +
#                 ggplot2::coord_quickmap(xlim = plot_extent[1:2], ylim = plot_extent[3:4]) +
#                 ggplot2::facet_wrap(~period) +
#                 ggplot2::labs(title = paste(toupper(substr(period, 1, 1)), substr(period, 2, nchar(period)), "mean for", variable)) +
#                 ggplot2::theme_minimal(base_size=14)+
#                 ggplot2::theme(panel.background = element_rect(fill="white"),
#                                plot.background = element_rect(fill="white"))
#         f <- file.path(output_dir, paste0(safe_name(variable), "_", period, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"))
#         ggplot2::ggsave(f, p, width = 12, height = 8)
#         f
# }



# plot_raster_comparison <- function(hist_df, proj_df = NULL, delta_df = NULL, plot_extent, variable, hist_range, proj_range, output_dir) {
#         if (nrow(hist_df) == 0) return(NULL)
#         vals <- c(hist_df$value, if (!is.null(proj_df)) proj_df$value, if (!is.null(delta_df)) delta_df$value)
#         vals <- vals[!is.na(vals)]
#         if (length(vals) == 0) return(NULL)
#         minv <- min(vals); maxv <- max(vals); mid <- (minv + maxv) / 2
#         timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
#         
#         title_hist <- paste(variable, "Historical", hist_range, "-", scenarios[1], "-", models[1])
#         title_proj <- paste(variable, "Projected", proj_range, "-", scenarios[2], "-", models[1])
#         title_delta <- paste(variable, "Delta (Proj - Hist)", "-", models[1])
#         
#         # Historical
#         p_hist <- ggplot2::ggplot(hist_df, ggplot2::aes(x=x,y=y,fill=value)) +
#                 ggplot2::geom_raster() +
#                 ggplot2::coord_quickmap(xlim=plot_extent[1:2], ylim=plot_extent[3:4]) +
#                 ggplot2::scale_fill_gradient2(limits=c(minv,maxv), low="blue", mid="white", high="red", midpoint=mid) +
#                 ggplot2::labs(title=title_hist, fill="Value") +
#                 ggplot2::theme_minimal(base_size=14)+
#                 ggplot2::theme(panel.background = element_rect(fill="white"),
#                                  plot.background = element_rect(fill="white"))
#         f_hist <- file.path(output_dir, paste0(safe_name(variable),"_hist_comparison_",timestamp,".png"))
#         ggplot2::ggsave(f_hist, p_hist, width=8, height=6)
#         
#         # Projected
#         if (!is.null(proj_df)) {
#                 p_proj <- ggplot2::ggplot(proj_df, ggplot2::aes(x=x,y=y,fill=value)) +
#                         ggplot2::geom_raster() +
#                         ggplot2::coord_quickmap(xlim=plot_extent[1:2], ylim=plot_extent[3:4]) +
#                         ggplot2::scale_fill_gradient2(limits=c(minv,maxv), low="blue", mid="white", high="red", midpoint=mid) +
#                         ggplot2::labs(title=title_proj, fill="Value") +
#                         ggplot2::theme_minimal(base_size=14) +
#                         ggplot2::theme(panel.background=element_rect(fill="white"))
#                 f_proj <- file.path(output_dir, paste0(safe_name(variable),"_proj_comparison_",timestamp,".png"))
#                 ggplot2::ggsave(f_proj, p_proj, width=8, height=6)
#         }
#         
#         # Delta
#         if (!is.null(delta_df)) {
#                 dl <- max(abs(delta_df$value), na.rm=TRUE)
#                 p_delta <- ggplot2::ggplot(delta_df, ggplot2::aes(x=x,y=y,fill=value)) +
#                         ggplot2::geom_raster() +
#                         ggplot2::coord_quickmap(xlim=plot_extent[1:2], ylim=plot_extent[3:4]) +
#                         ggplot2::scale_fill_gradient2(limits=c(-dl,dl), low="blue", mid="white", high="red", midpoint=0) +
#                         ggplot2::labs(title=title_delta, fill="Delta") +
#                         ggplot2::theme_minimal(base_size=14)+
#                         ggplot2::theme(panel.background = element_rect(fill="white"),
#                                        plot.background = element_rect(fill="white"))
#                 f_delta <- file.path(output_dir, paste0(safe_name(variable),"_delta_",timestamp,".png"))
#                 ggplot2::ggsave(f_delta, p_delta, width=8, height=6)
#         }
#         
#         return(list(hist=f_hist, proj=f_proj, delta=f_delta))
# }

plot_raster_comparison <- function(hist_df, proj_df = NULL, delta_df = NULL, plot_extent, variable, hist_range, proj_range, output_dir) {
        if (nrow(hist_df) == 0) return(NULL)
        vals <- c(hist_df$value, if (!is.null(proj_df)) proj_df$value, if (!is.null(delta_df)) delta_df$value)
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) return(NULL)
        minv <- min(vals); maxv <- max(vals); mid <- (minv + maxv) / 2
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        
        title_hist <- paste(variable, "Historical", hist_range, "-", scenarios[1], "-", models[1])
        title_proj <- paste(variable, "Projected", proj_range, "-", scenarios[2], "-", models[1])
        title_delta <- paste(variable, "Delta (Proj - Hist)", "-", models[1])
        
        # Historical
        p_hist <- ggplot2::ggplot(hist_df, ggplot2::aes(x=x,y=y,fill=value)) +
                ggplot2::geom_raster() +
                ggplot2::coord_quickmap(xlim=plot_extent[1:2], ylim=plot_extent[3:4]) +
                ggplot2::scale_fill_gradient2(limits=c(minv,maxv), low="blue", mid="white", high="red", midpoint=mid) +
                ggplot2::labs(title=title_hist, fill="Value") +
                ggplot2::theme_minimal(base_size=14) +
                ggplot2::labs(x = "Longitude", y = "Latitude")+
                ggplot2::theme(panel.background=element_rect(fill="white"))
        f_hist <- file.path(output_dir, paste0(safe_name(variable),"_hist_comparison_",timestamp,".png"))
        ggplot2::ggsave(f_hist, p_hist, width=8, height=6, bg = "white")
        
        # Projected
        f_proj <- NULL
        if (!is.null(proj_df)) {
                p_proj <- ggplot2::ggplot(proj_df, ggplot2::aes(x=x,y=y,fill=value)) +
                        ggplot2::geom_raster() +
                        ggplot2::coord_quickmap(xlim=plot_extent[1:2], ylim=plot_extent[3:4]) +
                        ggplot2::scale_fill_gradient2(limits=c(minv,maxv), low="blue", mid="white", high="red", midpoint=mid) +
                        ggplot2::labs(title=title_proj, fill="Value") +
                        ggplot2::theme_minimal(base_size=14) +
                        ggplot2::labs(x = "Longitude", y = "Latitude")+
                        ggplot2::theme(panel.background=element_rect(fill="white"))
                f_proj <- file.path(output_dir, paste0(safe_name(variable),"_proj_comparison_",timestamp,".png"))
                ggplot2::ggsave(f_proj, p_proj, width=8, height=6, bg = "white")
        }
        
        # Delta
        f_delta <- NULL
        if (!is.null(delta_df) && nrow(delta_df) > 0) {
                dl <- max(abs(delta_df$value), na.rm=TRUE)
                p_delta <- ggplot2::ggplot(delta_df, ggplot2::aes(x=x,y=y,fill=value)) +
                        ggplot2::geom_raster() +
                        ggplot2::coord_quickmap(xlim=plot_extent[1:2], ylim=plot_extent[3:4]) +
                        ggplot2::scale_fill_gradient2(limits=c(-dl,dl), low="blue", mid="white", high="red", midpoint=0) +
                        ggplot2::labs(title=title_delta, fill="Delta") +
                        ggplot2::theme_minimal(base_size=14) +
                        ggplot2::labs(x = "Longitude", y = "Latitude")+
                        ggplot2::theme(panel.background=element_rect(fill="white"))
                f_delta <- file.path(output_dir, paste0(safe_name(variable),"_delta_",timestamp,".png"))
                ggplot2::ggsave(f_delta, p_delta, width=8, height=6, bg = "white")
        }
        
        return(list(hist=f_hist, proj=f_proj, delta=f_delta))
}

plot_monthly_climatology <- function(
                monthly_df,
                plot_extent,
                variable,
                range_label,
                scenario,
                output_dir
) {
        
        if (nrow(monthly_df) == 0) {
                return(NULL)
        }
        
        monthly_df$month <- factor(
                monthly_df$layer,
                levels = paste0(
                        "monthly_",
                        1:12
                ),
                labels = month.abb
        )
        
        title_text <- paste(
                variable,
                "Monthly Climatology",
                range_label,
                "-",
                scenario,
                "-",
                models[1]
        )
        
        p <- ggplot2::ggplot(
                monthly_df,
                ggplot2::aes(
                        x = x,
                        y = y,
                        fill = value
                )
        ) +
                ggplot2::geom_raster() +
                ggplot2::coord_quickmap(
                        xlim = plot_extent[1:2],
                        ylim = plot_extent[3:4]
                ) +
                ggplot2::facet_wrap(
                        ~month,
                        ncol = 4
                ) +
                ggplot2::scale_fill_viridis_c(
                        option = "C"
                ) +
                ggplot2::labs(
                        title = title_text,
                        fill = "Value",
                        x = "Longitude",
                        y = "Latitude"
                ) +
                ggplot2::theme_minimal(
                        base_size = 14
                ) +
                ggplot2::theme(
                        panel.background =
                                element_rect(fill = "white"),
                        plot.background =
                                element_rect(fill = "white")
                )
        
        f <- file.path(
                output_dir,
                paste0(
                        safe_name(variable),
                        "_",
                        safe_name(scenario),
                        "_monthly_",
                        format(
                                Sys.time(),
                                "%Y%m%d_%H%M%S"
                        ),
                        ".png"
                )
        )
        
        ggplot2::ggsave(
                f,
                p,
                width = 12,
                height = 8,
                bg = "white"
        )
        
        f
}

plot_annual_or_monthly <- function(
                r_df,
                plot_extent,
                variable,
                period,
                range_label,
                scenario,
                output_dir
) {
        
        if (nrow(r_df) == 0) {
                return(NULL)
        }
        
        r_df$period <- factor(
                r_df$layer,
                levels = unique(r_df$layer)
        )
        
        title_text <- paste(
                variable,
                toupper(period),
                "Climatology",
                range_label,
                "-",
                scenario,
                "-",
                models[1]
        )
        
        p <- ggplot2::ggplot(
                r_df,
                ggplot2::aes(
                        x = x,
                        y = y,
                        fill = value
                )
        ) +
                ggplot2::geom_raster() +
                ggplot2::coord_quickmap(
                        xlim = plot_extent[1:2],
                        ylim = plot_extent[3:4]
                ) +
                ggplot2::facet_wrap(~period) +
                ggplot2::scale_fill_viridis_c(
                        option = "C"
                ) +
                ggplot2::labs(
                        title = title_text,
                        fill = "Value",
                        x = "Longitude",
                        y = "Latitude"
                ) +
                ggplot2::theme_minimal(
                        base_size = 14
                ) +
                ggplot2::theme(
                        panel.background =
                                element_rect(fill = "white"),
                        plot.background =
                                element_rect(fill = "white")
                )
        
        f <- file.path(
                output_dir,
                paste0(
                        safe_name(variable),
                        "_",
                        safe_name(scenario),
                        "_",
                        period,
                        "_",
                        format(
                                Sys.time(),
                                "%Y%m%d_%H%M%S"
                        ),
                        ".png"
                )
        )
        
        ggplot2::ggsave(
                f,
                p,
                width = 12,
                height = 8,
                bg = "white"
        )
        
        f
}

plot_percentile_map <- function(
                percentile_df,
                plot_extent,
                variable,
                percentile,
                range_label,
                scenario,
                output_dir,
                common_limits = NULL
) {
        
        if (nrow(percentile_df) == 0) {
                return(NULL)
        }
        
        title_text <- paste(
                variable,
                percentile,
                "Percentile",
                range_label,
                "-",
                scenario,
                "-",
                models[1]
        )
        
        p <- ggplot2::ggplot(
                percentile_df,
                ggplot2::aes(
                        x = x,
                        y = y,
                        fill = value
                )
        ) +
                ggplot2::geom_raster() +
                ggplot2::coord_quickmap(
                        xlim = plot_extent[1:2],
                        ylim = plot_extent[3:4]
                ) +
                ggplot2::scale_fill_viridis_c(
                        option = "C",
                        limits = common_limits
                ) +
                ggplot2::labs(
                        title = title_text,
                        fill = "Value",
                        x = "Longitude",
                        y = "Latitude"
                ) +
                ggplot2::theme_minimal(
                        base_size = 14
                ) +
                ggplot2::theme(
                        panel.background =
                                element_rect(fill = "white"),
                        plot.background =
                                element_rect(fill = "white")
                )
        
        f <- file.path(
                output_dir,
                paste0(
                        safe_name(variable),
                        "_",
                        safe_name(scenario),
                        "_percentile_",
                        percentile,
                        "_",
                        format(
                                Sys.time(),
                                "%Y%m%d_%H%M%S"
                        ),
                        ".png"
                )
        )
        
        ggplot2::ggsave(
                f,
                p,
                width = 8,
                height = 6,
                bg = "white"
        )
        
        f
}

percentile_map <- function(r, prob) {
        if (is.null(r)) return(NULL)
        tryCatch(terra::app(r, fun = function(x) stats::quantile(x, probs = prob, na.rm = TRUE)), error = function(e) NULL)
}

###############################################################################
# Threshold exceedance count
#
# Returns one raster containing the number of time steps at or above
# the selected threshold.
###############################################################################

compute_threshold_days <- function(
                r,
                threshold,
                inclusive = TRUE
) {
        
        if (
                is.null(r) ||
                terra::nlyr(r) == 0
        ) {
                return(NULL)
        }
        
        if (isTRUE(inclusive)) {
                
                out <- terra::app(
                        r,
                        fun = function(x) {
                                sum(
                                        x >= threshold,
                                        na.rm = TRUE
                                )
                        }
                )
                
        } else {
                
                out <- terra::app(
                        r,
                        fun = function(x) {
                                sum(
                                        x > threshold,
                                        na.rm = TRUE
                                )
                        }
                )
        }
        
        names(out) <- paste0(
                "days_ge_",
                safe_name(threshold)
        )
        
        out
}

###############################################################################
# Consecutive threshold-event metrics
#
# Calculates period-wide event statistics from a DAILY raster:
#
#   event_days   = total number of days belonging to qualifying events
#   max_duration = longest qualifying event
#   event_count  = number of qualifying events
#
# A qualifying event is >= min_duration consecutive days meeting the
# threshold criterion.
#
# IMPORTANT:
# threshold must be expressed in the SAME UNITS as the input raster.
###############################################################################

compute_consecutive_event_metrics <- function(
                r,
                threshold,
                min_duration = 3,
                inclusive = TRUE
) {
        
        if (
                is.null(r) ||
                terra::nlyr(r) == 0
        ) {
                return(NULL)
        }
        
        if (
                !is.numeric(threshold) ||
                length(threshold) != 1L ||
                is.na(threshold)
        ) {
                stop("threshold must be one finite numeric value.")
        }
        
        if (
                !is.numeric(min_duration) ||
                length(min_duration) != 1L ||
                is.na(min_duration) ||
                min_duration < 1
        ) {
                stop("min_duration must be a positive integer.")
        }
        
        min_duration <- as.integer(
                min_duration
        )
        
        eligibility <- check_metric_temporal_eligibility(
                r = r,
                metric = "consecutive_daily_event"
        )
        
        if (!isTRUE(eligibility$eligible)) {
                
                stop(
                        "Consecutive-event metric cannot be calculated: ",
                        eligibility$reason
                )
        }
        # -------------------------------------------------------------------------
        # Validate time coordinate
        # -------------------------------------------------------------------------
        
        tt <- tryCatch(
                terra::time(r),
                error = function(e) NULL
        )
        
        if (
                is.null(tt) ||
                length(tt) != terra::nlyr(r)
        ) {
                stop(
                        "Consecutive-event calculation requires a valid daily ",
                        "time coordinate."
                )
        }
        
        dates <- as.Date(tt)
        
        if (anyNA(dates)) {
                stop(
                        "Consecutive-event calculation found invalid dates."
                )
        }
        
        if (anyDuplicated(dates)) {
                stop(
                        "Consecutive-event calculation found duplicated dates."
                )
        }
        
        if (length(dates) > 1L) {
                
                step_days <- as.integer(
                        diff(dates)
                )
                
                if (any(step_days != 1L)) {
                        
                        stop(
                                "Consecutive-event calculation requires consecutive ",
                                "daily observations. Missing/non-daily dates detected."
                        )
                }
        }
        
        # -------------------------------------------------------------------------
        # Cell-level event summary
        # -------------------------------------------------------------------------
        
        event_fun <- function(x) {
                
                if (all(is.na(x))) {
                        
                        return(
                                c(
                                        event_days = NA_real_,
                                        max_duration = NA_real_,
                                        event_count = NA_real_
                                )
                        )
                }
                
                if (isTRUE(inclusive)) {
                        exceed <- x >= threshold
                } else {
                        exceed <- x > threshold
                }
                
                # Missing observations break an event.
                exceed[is.na(exceed)] <- FALSE
                
                rr <- rle(
                        as.logical(exceed)
                )
                
                qualifying_lengths <- rr$lengths[
                        rr$values &
                                rr$lengths >= min_duration
                ]
                
                if (length(qualifying_lengths) == 0L) {
                        
                        return(
                                c(
                                        event_days = 0,
                                        max_duration = 0,
                                        event_count = 0
                                )
                        )
                }
                
                c(
                        event_days =
                                sum(qualifying_lengths),
                        
                        max_duration =
                                max(qualifying_lengths),
                        
                        event_count =
                                length(qualifying_lengths)
                )
        }
        
        out <- terra::app(
                r,
                fun = event_fun
        )
        
        names(out) <- c(
                "event_days",
                "max_duration",
                "event_count"
        )
        
        out
}

###############################################################################
# Detect temporal resolution
#
# Uses the actual raster time coordinate.
#
# Returns a small metadata list rather than changing the raster.
###############################################################################

detect_temporal_resolution <- function(r) {
        
        if (
                is.null(r) ||
                !inherits(r, "SpatRaster") ||
                terra::nlyr(r) == 0
        ) {
                
                return(
                        list(
                                valid = FALSE,
                                resolution = "unknown",
                                step_seconds = NA_real_,
                                regular = FALSE,
                                n_times = 0L,
                                reason = "Raster is NULL or has no layers."
                        )
                )
        }
        
        tt <- tryCatch(
                terra::time(r),
                error = function(e) NULL
        )
        
        if (
                is.null(tt) ||
                length(tt) != terra::nlyr(r)
        ) {
                
                return(
                        list(
                                valid = FALSE,
                                resolution = "unknown",
                                step_seconds = NA_real_,
                                regular = FALSE,
                                n_times = terra::nlyr(r),
                                reason = "No valid time coordinate."
                        )
                )
        }
        
        if (length(tt) == 1L) {
                
                return(
                        list(
                                valid = TRUE,
                                resolution = "single_time",
                                step_seconds = NA_real_,
                                regular = NA,
                                n_times = 1L,
                                reason = NULL
                        )
                )
        }
        
        # Convert all time differences to seconds.
        step_seconds <- as.numeric(
                diff(tt),
                units = "secs"
        )
        
        # Use median so a DST/calendar representation issue does not dominate.
        median_step <- stats::median(
                step_seconds,
                na.rm = TRUE
        )
        
        tolerance <- max(
                1,
                abs(median_step) * 1e-6
        )
        
        regular <- all(
                abs(step_seconds - median_step) <= tolerance
        )
        
        # Classification boundaries.
        resolution <- dplyr::case_when(
                
                abs(median_step - 3600) <= 60 ~
                        "hourly",
                
                abs(median_step - 3 * 3600) <= 60 ~
                        "3_hourly",
                
                abs(median_step - 6 * 3600) <= 60 ~
                        "6_hourly",
                
                abs(median_step - 12 * 3600) <= 60 ~
                        "12_hourly",
                
                abs(median_step - 24 * 3600) <= 60 ~
                        "daily",
                
                abs(median_step - 7 * 24 * 3600) <= 3600 ~
                        "weekly",
                
                median_step >= 27 * 24 * 3600 &
                        median_step <= 32 * 24 * 3600 ~
                        "monthly",
                
                median_step >= 360 * 24 * 3600 &
                        median_step <= 370 * 24 * 3600 ~
                        "annual",
                
                TRUE ~
                        "other"
        )
        
        list(
                valid = TRUE,
                resolution = resolution,
                step_seconds = median_step,
                regular = regular,
                n_times = length(tt),
                first_time = min(tt),
                last_time = max(tt),
                reason = NULL
        )
}

#====
# Check whether a requested climate statistic is valid for the temporal
# resolution supplied.
#====

check_metric_temporal_eligibility <- function(
                r,
                metric,
                required_resolution = NULL
) {
        
        info <- detect_temporal_resolution(r)
        
        if (!isTRUE(info$valid)) {
                
                return(
                        list(
                                eligible = FALSE,
                                temporal = info,
                                reason = info$reason
                        )
                )
        }
        
        metric <- tolower(metric)
        
        # -------------------------------------------------------------------------
        # Default metric requirements
        # -------------------------------------------------------------------------
        
        default_requirements <- list(
                
                consecutive_daily_event =
                        c("daily"),
                
                daily_threshold_count =
                        c("daily"),
                
                heatwave =
                        c("daily"),
                
                annual_mean =
                        c(
                                "hourly",
                                "3_hourly",
                                "6_hourly",
                                "12_hourly",
                                "daily",
                                "weekly",
                                "monthly"
                        ),
                
                monthly_climatology =
                        c(
                                "hourly",
                                "3_hourly",
                                "6_hourly",
                                "12_hourly",
                                "daily",
                                "weekly",
                                "monthly"
                        ),
                
                percentile =
                        c(
                                "hourly",
                                "3_hourly",
                                "6_hourly",
                                "12_hourly",
                                "daily",
                                "weekly",
                                "monthly"
                        )
        )
        
        if (!is.null(required_resolution)) {
                
                allowed <- required_resolution
                
        } else if (metric %in% names(default_requirements)) {
                
                allowed <- default_requirements[[metric]]
                
        } else {
                
                return(
                        list(
                                eligible = FALSE,
                                temporal = info,
                                reason = paste0(
                                        "No temporal-resolution rule has been defined for metric '",
                                        metric,
                                        "'."
                                )
                        )
                )
        }
        
        if (!isTRUE(info$regular)) {
                
                return(
                        list(
                                eligible = FALSE,
                                temporal = info,
                                reason =
                                        "Time coordinate is irregular; temporal metric not calculated."
                        )
                )
        }
        
        if (!info$resolution %in% allowed) {
                
                return(
                        list(
                                eligible = FALSE,
                                temporal = info,
                                reason = paste0(
                                        "Metric '",
                                        metric,
                                        "' requires ",
                                        paste(
                                                allowed,
                                                collapse = " / "
                                        ),
                                        " data, but input is ",
                                        info$resolution,
                                        "."
                                )
                        )
                )
        }
        
        list(
                eligible = TRUE,
                temporal = info,
                reason = NULL
        )
}

#====
# Check whether historical and projection metric products may be compared.
#====

check_metric_compatibility <- function(
                hist_meta,
                proj_meta
) {
        
        required_fields <- c(
                "metric",
                "units",
                "temporal_resolution",
                "definition_id"
        )
        
        missing_hist <- setdiff(
                required_fields,
                names(hist_meta)
        )
        
        missing_proj <- setdiff(
                required_fields,
                names(proj_meta)
        )
        
        if (
                length(missing_hist) > 0 ||
                length(missing_proj) > 0
        ) {
                
                return(
                        list(
                                compatible = FALSE,
                                reason =
                                        "Historical/projection metric metadata are incomplete."
                        )
                )
        }
        
        checks <- c(
                
                metric =
                        identical(
                                hist_meta$metric,
                                proj_meta$metric
                        ),
                
                units =
                        identical(
                                hist_meta$units,
                                proj_meta$units
                        ),
                
                temporal_resolution =
                        identical(
                                hist_meta$temporal_resolution,
                                proj_meta$temporal_resolution
                        ),
                
                definition =
                        identical(
                                hist_meta$definition_id,
                                proj_meta$definition_id
                        )
        )
        
        if (!all(checks)) {
                
                failed <- names(checks)[
                        !checks
                ]
                
                return(
                        list(
                                compatible = FALSE,
                                reason = paste0(
                                        "Historical and projection metrics cannot be compared: ",
                                        paste(
                                                failed,
                                                collapse = ", "
                                        ),
                                        " differ."
                                )
                        )
                )
        }
        
        list(
                compatible = TRUE,
                reason = NULL
        )
}

#====
# Metric source / provenance helpers
#====

make_metric_metadata <- function(
                dataset,
                metric,
                source_variable,
                source_mode = c("precalculated", "derived"),
                units,
                temporal_resolution,
                definition_id,
                aggregation = NULL,
                threshold = NULL,
                min_duration = NULL
) {
        
        source_mode <- match.arg(source_mode)
        
        list(
                dataset = dataset,
                metric = metric,
                source_variable = source_variable,
                source_mode = source_mode,
                units = units,
                temporal_resolution = temporal_resolution,
                definition_id = definition_id,
                aggregation = aggregation,
                threshold = threshold,
                min_duration = min_duration
        )
}


check_metric_compatibility <- function(
                hist_meta,
                proj_meta
) {
        
        required_fields <- c(
                "metric",
                "units",
                "temporal_resolution",
                "definition_id"
        )
        
        if (
                !all(required_fields %in% names(hist_meta)) ||
                !all(required_fields %in% names(proj_meta))
        ) {
                
                return(
                        list(
                                compatible = FALSE,
                                reason =
                                        "Historical/projection metric metadata are incomplete."
                        )
                )
        }
        
        checks <- c(
                metric =
                        identical(
                                hist_meta$metric,
                                proj_meta$metric
                        ),
                
                units =
                        identical(
                                hist_meta$units,
                                proj_meta$units
                        ),
                
                temporal_resolution =
                        identical(
                                hist_meta$temporal_resolution,
                                proj_meta$temporal_resolution
                        ),
                
                definition =
                        identical(
                                hist_meta$definition_id,
                                proj_meta$definition_id
                        )
        )
        
        if (!all(checks)) {
                
                failed <- names(checks)[!checks]
                
                return(
                        list(
                                compatible = FALSE,
                                reason = paste0(
                                        "Historical and projection metrics are not directly comparable: ",
                                        paste(
                                                failed,
                                                collapse = ", "
                                        ),
                                        " differ."
                                )
                        )
                )
        }
        
        list(
                compatible = TRUE,
                reason = NULL
        )
}

#====
# Resolve whether a metric is direct, derived, or unavailable
#
# Resolve whether a requested metric is:
#   1. available as an approved pre-calculated product,
#   2. derivable from required source variables,
#   3. unavailable.
#
# IMPORTANT:
# Presence of a variable with the same name is NOT automatically treated
# as an approved pre-calculated metric.
#====

resolve_metric_source <- function(
                requested_metric,
                available_variables,
                direct_metric_variables = character(0),
                derivation_requirements = list()
) {
        
        requested_metric_lc <- tolower(
                requested_metric
        )
        
        available_lc <- tolower(
                available_variables
        )
        
        direct_lc <- tolower(
                direct_metric_variables
        )
        
        # -------------------------------------------------------------------------
        # 1. Approved pre-calculated product
        # -------------------------------------------------------------------------
        
        if (
                requested_metric_lc %in% direct_lc &&
                requested_metric_lc %in% available_lc
        ) {
                
                matched_name <- available_variables[
                        match(
                                requested_metric_lc,
                                available_lc
                        )
                ]
                
                return(
                        list(
                                available = TRUE,
                                mode = "precalculated",
                                metric = requested_metric,
                                source_variable = matched_name,
                                required_variables = matched_name,
                                reason = NULL
                        )
                )
        }
        

        # -------------------------------------------------------------------------
        # 2. Derivation rule
        # Match derivation-rule names case-insensitively.
        # -------------------------------------------------------------------------
        
        if (length(derivation_requirements) > 0L) {
                
                derivation_names <- names(derivation_requirements)
                
                if (
                        !is.null(derivation_names) &&
                        length(derivation_names) > 0L
                ) {
                        
                        derivation_names_lc <- tolower(derivation_names)
                        
                        derivation_match <- match(
                                requested_metric_lc,
                                derivation_names_lc
                        )
                        
                        if (!is.na(derivation_match)) {
                                
                                req <- derivation_requirements[[derivation_match]]
                                
                                req_lc <- tolower(req)
                                
                                present <- req_lc %in% available_lc
                                
                                if (all(present)) {
                                        
                                        matched_variables <- available_variables[
                                                match(
                                                        req_lc,
                                                        available_lc
                                                )
                                        ]
                                        
                                        return(
                                                list(
                                                        available = TRUE,
                                                        mode = "derived",
                                                        metric = requested_metric,
                                                        source_variable = NULL,
                                                        required_variables = matched_variables,
                                                        missing_variables = character(0),
                                                        reason = NULL
                                                )
                                        )
                                }
                                
                                missing_variables <- req[!present]
                                
                                return(
                                        list(
                                                available = FALSE,
                                                mode = "unavailable",
                                                metric = requested_metric,
                                                source_variable = NULL,
                                                required_variables = req,
                                                missing_variables = missing_variables,
                                                reason = paste0(
                                                        "Metric '",
                                                        requested_metric,
                                                        "' cannot be derived because required source variables are missing: ",
                                                        paste(
                                                                missing_variables,
                                                                collapse = ", "
                                                        )
                                                )
                                        )
                                )
                        }
                }
        }
        # -------------------------------------------------------------------------
        # 3. No approved direct product and no derivation rule
        # -------------------------------------------------------------------------
        
        list(
                available = FALSE,
                mode = "unavailable",
                metric = requested_metric,
                source_variable = NULL,
                required_variables = character(0),
                reason = paste0(
                        "Metric '",
                        requested_metric,
                        "' is not available as an approved pre-calculated product ",
                        "and no derivation rule has been defined."
                )
        )
}


#====
# Decide whether historical/projection metric comparison is permitted
#
# Rules:
#   - both unavailable          -> nothing to process
#   - projection only          -> projection-only output
#   - historical only          -> historical-only output
#   - both available + same definition -> comparison/delta allowed
#   - both available + incompatible definitions -> separate outputs only
#====

decide_metric_comparison <- function(
                hist_resolution,
                proj_resolution,
                hist_meta = NULL,
                proj_meta = NULL
) {
        
        hist_available <- isTRUE(hist_resolution$available)
        proj_available <- isTRUE(proj_resolution$available)
        
        # Neither available
        if (!hist_available && !proj_available) {
                
                return(
                        list(
                                mode = "unavailable",
                                historical_allowed = FALSE,
                                projection_allowed = FALSE,
                                comparison_allowed = FALSE,
                                delta_allowed = FALSE,
                                reason = paste(
                                        "Metric unavailable for both historical and projection data.",
                                        "Historical:",
                                        hist_resolution$reason,
                                        "Projection:",
                                        proj_resolution$reason
                                )
                        )
                )
        }
        
        # Projection only
        if (!hist_available && proj_available) {
                
                return(
                        list(
                                mode = "projection_only",
                                historical_allowed = FALSE,
                                projection_allowed = TRUE,
                                comparison_allowed = FALSE,
                                delta_allowed = FALSE,
                                reason = paste0(
                                        "Projection metric is available, but no scientifically ",
                                        "equivalent historical metric can currently be produced. ",
                                        "Projection-only output is permitted."
                                )
                        )
                )
        }
        
        # Historical only
        if (hist_available && !proj_available) {
                
                return(
                        list(
                                mode = "historical_only",
                                historical_allowed = TRUE,
                                projection_allowed = FALSE,
                                comparison_allowed = FALSE,
                                delta_allowed = FALSE,
                                reason = paste0(
                                        "Historical metric is available, but projection metric ",
                                        "is unavailable. Historical-only output is permitted."
                                )
                        )
                )
        }
        
        # Both products exist, but metadata are required before comparison.
        if (is.null(hist_meta) || is.null(proj_meta)) {
                
                return(
                        list(
                                mode = "separate_only",
                                historical_allowed = TRUE,
                                projection_allowed = TRUE,
                                comparison_allowed = FALSE,
                                delta_allowed = FALSE,
                                reason = paste0(
                                        "Historical and projection metrics are both available, ",
                                        "but comparison metadata are incomplete."
                                )
                        )
                )
        }
        
        compatibility <- check_metric_compatibility(
                hist_meta,
                proj_meta
        )
        
        if (!isTRUE(compatibility$compatible)) {
                
                return(
                        list(
                                mode = "separate_only",
                                historical_allowed = TRUE,
                                projection_allowed = TRUE,
                                comparison_allowed = FALSE,
                                delta_allowed = FALSE,
                                reason = compatibility$reason
                        )
                )
        }
        
        # Both available and scientifically compatible
        list(
                mode = "comparison",
                historical_allowed = TRUE,
                projection_allowed = TRUE,
                comparison_allowed = TRUE,
                delta_allowed = TRUE,
                reason = NULL
        )
}

extract_available_variable_names <- function(files) {
        
        if (
                is.null(files) ||
                length(files) == 0
        ) {
                return(character(0))
        }
        
        basename_files <- basename(files)
        
        unique(
                sub(
                        "_.*$",
                        "",
                        basename_files
                )
        )
}

###############################################################################
# Build metadata from an actual resolved metric raster
#
# Extracts properties that can safely be determined from the data itself:
#   temporal resolution
#   regularity
#   number of observations
#   time range
#   raster units where available
#
# Scientific definition_id is supplied explicitly because it should NOT
# be guessed from a filename.
###############################################################################

inspect_metric_raster <- function(
                raster_file,
                dataset,
                metric,
                source_variable,
                source_mode,
                definition_id,
                units_override = NULL
) {
        
        if (
                is.null(raster_file) ||
                length(raster_file) == 0L ||
                !file.exists(raster_file)
        ) {
                return(NULL)
        }
        
        r <- terra::rast(raster_file)
        
        if (
                terra::nlyr(r) == 0L ||
                !terra::hasValues(r)
        ) {
                return(NULL)
        }
        
        temporal <- detect_temporal_resolution(r)
        
        # Try to read units from the raster metadata.
        raster_units <- tryCatch(
                terra::units(r),
                error = function(e) NULL
        )
        
        if (
                !is.null(units_override) &&
                length(units_override) == 1L &&
                !is.na(units_override)
        ) {
                
                units_final <- as.character(
                        units_override
                )
                
        } else if (
                !is.null(raster_units) &&
                length(raster_units) > 0L &&
                !all(is.na(raster_units))
        ) {
                
                units_final <- unique(
                        as.character(
                                raster_units[
                                        !is.na(raster_units)
                                ]
                        )
                )
                
                units_final <- paste(
                        units_final,
                        collapse = ";"
                )
                
        } else {
                
                units_final <- NA_character_
        }
        
        list(
                dataset = dataset,
                metric = metric,
                source_variable = source_variable,
                source_mode = source_mode,
                
                units = units_final,
                
                temporal_resolution =
                        temporal$resolution,
                
                temporal_regular =
                        temporal$regular,
                
                n_times =
                        temporal$n_times,
                
                first_time =
                        if (!is.null(temporal$first_time)) {
                                temporal$first_time
                        } else {
                                NA
                        },
                
                last_time =
                        if (!is.null(temporal$last_time)) {
                                temporal$last_time
                        } else {
                                NA
                        },
                
                definition_id = definition_id
        )
}

#====
# Parse NARCliM member metadata from filenames
#
# Parse NARCliM member metadata from filenames
#
# Expected current filename structure:
#
# variable_scenario_GCM_NARCliM...RCM_DOMAIN.nc
#
# Example:
# FFDI_ssp370_EC-Earth3-Veg_NARCliM2-0-WRF412R3_AUS-20i.nc
#====

parse_narclim_member <- function(files) {
        
        if (
                is.null(files) ||
                length(files) == 0L
        ) {
                
                return(
                        data.frame(
                                file = character(0),
                                variable = character(0),
                                scenario = character(0),
                                gcm = character(0),
                                rcm = character(0),
                                domain = character(0),
                                stringsAsFactors = FALSE
                        )
                )
        }
        
        b <- basename(files)
        
        # Remove extension only.
        stem <- sub(
                "\\.nc$",
                "",
                b,
                ignore.case = TRUE
        )
        
        # -------------------------------------------------------------------------
        # Variable and scenario are the first two underscore-separated fields.
        # -------------------------------------------------------------------------
        
        parts <- strsplit(
                stem,
                "_",
                fixed = TRUE
        )
        
        variable <- vapply(
                parts,
                function(x) {
                        if (length(x) >= 1L) x[1] else NA_character_
                },
                character(1)
        )
        
        scenario <- vapply(
                parts,
                function(x) {
                        if (length(x) >= 2L) x[2] else NA_character_
                },
                character(1)
        )
        
        # -------------------------------------------------------------------------
        # GCM
        #
        # Remove each file's OWN variable + scenario prefix.
        # This must be done row-by-row rather than with a vector regex.
        # -------------------------------------------------------------------------
        
        gcm <- vapply(
                seq_along(stem),
                function(i) {
                        
                        prefix <- paste0(
                                variable[i],
                                "_",
                                scenario[i],
                                "_"
                        )
                        
                        if (!startsWith(
                                tolower(stem[i]),
                                tolower(prefix)
                        )) {
                                return(NA_character_)
                        }
                        
                        remainder <- substring(
                                stem[i],
                                nchar(prefix) + 1L
                        )
                        
                        if (!grepl(
                                "_NARCliM",
                                remainder,
                                ignore.case = TRUE
                        )) {
                                return(NA_character_)
                        }
                        
                        sub(
                                "_NARCliM.*$",
                                "",
                                remainder,
                                ignore.case = TRUE
                        )
                },
                character(1)
        )
        
        # -------------------------------------------------------------------------
        # RCM
        # -------------------------------------------------------------------------
        
        rcm <- ifelse(
                grepl(
                        "WRF412R3",
                        b,
                        ignore.case = TRUE
                ),
                "WRF412R3",
                ifelse(
                        grepl(
                                "WRF412R5",
                                b,
                                ignore.case = TRUE
                        ),
                        "WRF412R5",
                        NA_character_
                )
        )
        
        # -------------------------------------------------------------------------
        # Domain
        # -------------------------------------------------------------------------
        
        domain <- ifelse(
                grepl(
                        "AUS-20",
                        b,
                        ignore.case = TRUE
                ),
                "AUS20",
                ifelse(
                        grepl(
                                "SEAus-04",
                                b,
                                ignore.case = TRUE
                        ),
                        "SEAus04",
                        NA_character_
                )
        )
        
        out <- data.frame(
                file = files,
                variable = variable,
                scenario = scenario,
                gcm = gcm,
                rcm = rcm,
                domain = domain,
                stringsAsFactors = FALSE
        )
        
        rownames(out) <- NULL
        
        out
}

#====
# Select NARCliM members according to user model mode
#====

select_narclim_members <- function(
                member_table,
                gcms = NULL,
                rcms = NULL,
                model_mode = c(
                        "individual",
                        "ensemble_R3",
                        "ensemble_R5",
                        "ensemble_all"
                )
) {
        
        model_mode <- match.arg(model_mode)
        
        if (
                is.null(member_table) ||
                nrow(member_table) == 0L
        ) {
                return(member_table)
        }
        
        x <- member_table
        
        if (
                !is.null(gcms) &&
                length(gcms) > 0L
        ) {
                x <- x[
                        tolower(x$gcm) %in% tolower(gcms),
                        ,
                        drop = FALSE
                ]
        }
        
        if (model_mode == "individual") {
                
                if (
                        !is.null(rcms) &&
                        length(rcms) > 0L
                ) {
                        x <- x[
                                tolower(x$rcm) %in% tolower(rcms),
                                ,
                                drop = FALSE
                        ]
                }
                
        } else if (model_mode == "ensemble_R3") {
                
                x <- x[
                        tolower(x$rcm) == "wrf412r3",
                        ,
                        drop = FALSE
                ]
                
        } else if (model_mode == "ensemble_R5") {
                
                x <- x[
                        tolower(x$rcm) == "wrf412r5",
                        ,
                        drop = FALSE
                ]
                
        } else if (model_mode == "ensemble_all") {
                
                x <- x[
                        tolower(x$rcm) %in%
                                c(
                                        "wrf412r3",
                                        "wrf412r5"
                                ),
                        ,
                        drop = FALSE
                ]
        }
        
        rownames(x) <- NULL
        
        x
}

#====
# Compute ensemble summary from already-comparable member rasters
#
# Returns:
#   median
#   p10
#   p90
#   min
#   max
#====

###############################################################################
# Ensemble summary across member-level metric rasters
#
# Each input raster represents the metric already calculated independently
# for one climate-model member.
#
# Ensemble statistics are calculated cell-by-cell across members.
###############################################################################

###############################################################################
# Ensemble summary across member-level metric rasters
#
# Each input raster must already represent the requested metric calculated
# independently for one climate-model member.
#
# Ensemble statistics are then calculated CELL-BY-CELL across members.
#
# Outputs:
#   ensemble_median
#   ensemble_p10
#   ensemble_p90
#   ensemble_min
#   ensemble_max
###############################################################################

compute_ensemble_summary <- function(
                raster_files,
                align_method = "bilinear"
) {
        
        # ---------------------------------------------------------------------
        # Normalise file input
        # ---------------------------------------------------------------------
        
        if (
                is.null(raster_files) ||
                length(raster_files) == 0L
        ) {
                return(NULL)
        }
        
        raster_files <- unlist(
                raster_files,
                use.names = FALSE
        )
        
        raster_files <- raster_files[
                !is.na(raster_files) &
                        nzchar(raster_files) &
                        file.exists(raster_files)
        ]
        
        if (length(raster_files) == 0L) {
                return(NULL)
        }
        
        # ---------------------------------------------------------------------
        # First member defines reference geometry
        # ---------------------------------------------------------------------
        
        ref <- terra::rast(
                raster_files[[1]]
        )
        
        if (
                terra::nlyr(ref) != 1L
        ) {
                stop(
                        "Each member metric raster must contain exactly one layer."
                )
        }
        
        # ---------------------------------------------------------------------
        # ONE-MEMBER ENSEMBLE
        #
        # This is an important explicit edge case.
        #
        # With only one member:
        #
        # median = p10 = p90 = min = max = member value
        #
        # Avoid terra::app() because no across-member reduction is required.
        # ---------------------------------------------------------------------
        
        if (length(raster_files) == 1L) {
                
                out <- c(
                        ref,
                        ref,
                        ref,
                        ref,
                        ref
                )
                
                names(out) <- c(
                        "ensemble_median",
                        "ensemble_p10",
                        "ensemble_p90",
                        "ensemble_min",
                        "ensemble_max"
                )
                
                return(out)
        }
        
        # ---------------------------------------------------------------------
        # MULTI-MEMBER ENSEMBLE
        #
        # Align every member to the first member before stacking.
        # ---------------------------------------------------------------------
        
        ensemble_stack <- ref
        
        for (i in 2:length(raster_files)) {
                
                r <- terra::rast(
                        raster_files[[i]]
                )
                
                if (
                        terra::nlyr(r) != 1L
                ) {
                        stop(
                                paste0(
                                        "Member ",
                                        i,
                                        " contains ",
                                        terra::nlyr(r),
                                        " layers. Expected exactly one ",
                                        "metric layer per ensemble member."
                                )
                        )
                }
                
                same_geometry <- terra::compareGeom(
                        ref,
                        r,
                        stopOnError = FALSE
                )
                
                if (!isTRUE(same_geometry)) {
                        
                        r <- align_and_resample(
                                ref = ref,
                                target = r,
                                method = align_method
                        )
                }
                
                ensemble_stack <- c(
                        ensemble_stack,
                        r
                )
        }
        
        if (
                terra::nlyr(ensemble_stack) !=
                length(raster_files)
        ) {
                stop(
                        paste0(
                                "Ensemble stack contains ",
                                terra::nlyr(ensemble_stack),
                                " layers but ",
                                length(raster_files),
                                " member rasters were supplied."
                        )
                )
        }
        
        # ---------------------------------------------------------------------
        # Cell-wise statistics across climate-model members
        #
        # terra::app() supplies a matrix for multi-layer input.
        # Each row represents one grid cell; columns are ensemble members.
        # ---------------------------------------------------------------------
        
        ensemble_median <- terra::app(
                ensemble_stack,
                fun = function(x) {
                        
                        if (is.matrix(x)) {
                                
                                apply(
                                        x,
                                        1,
                                        stats::median,
                                        na.rm = TRUE
                                )
                                
                        } else {
                                
                                stats::median(
                                        x,
                                        na.rm = TRUE
                                )
                        }
                }
        )
        
        ensemble_p10 <- terra::app(
                ensemble_stack,
                fun = function(x) {
                        
                        if (is.matrix(x)) {
                                
                                apply(
                                        x,
                                        1,
                                        function(v) {
                                                stats::quantile(
                                                        v,
                                                        probs = 0.10,
                                                        na.rm = TRUE,
                                                        names = FALSE
                                                )
                                        }
                                )
                                
                        } else {
                                
                                stats::quantile(
                                        x,
                                        probs = 0.10,
                                        na.rm = TRUE,
                                        names = FALSE
                                )
                        }
                }
        )
        
        ensemble_p90 <- terra::app(
                ensemble_stack,
                fun = function(x) {
                        
                        if (is.matrix(x)) {
                                
                                apply(
                                        x,
                                        1,
                                        function(v) {
                                                stats::quantile(
                                                        v,
                                                        probs = 0.90,
                                                        na.rm = TRUE,
                                                        names = FALSE
                                                )
                                        }
                                )
                                
                        } else {
                                
                                stats::quantile(
                                        x,
                                        probs = 0.90,
                                        na.rm = TRUE,
                                        names = FALSE
                                )
                        }
                }
        )
        
        ensemble_min <- terra::app(
                ensemble_stack,
                fun = function(x) {
                        
                        if (is.matrix(x)) {
                                
                                apply(
                                        x,
                                        1,
                                        min,
                                        na.rm = TRUE
                                )
                                
                        } else {
                                
                                min(
                                        x,
                                        na.rm = TRUE
                                )
                        }
                }
        )
        
        ensemble_max <- terra::app(
                ensemble_stack,
                fun = function(x) {
                        
                        if (is.matrix(x)) {
                                
                                apply(
                                        x,
                                        1,
                                        max,
                                        na.rm = TRUE
                                )
                                
                        } else {
                                
                                max(
                                        x,
                                        na.rm = TRUE
                                )
                        }
                }
        )
        
        # ---------------------------------------------------------------------
        # Combine ensemble products
        # ---------------------------------------------------------------------
        
        out <- c(
                ensemble_median,
                ensemble_p10,
                ensemble_p90,
                ensemble_min,
                ensemble_max
        )
        
        names(out) <- c(
                "ensemble_median",
                "ensemble_p10",
                "ensemble_p90",
                "ensemble_min",
                "ensemble_max"
        )
        
        out
}
#====
# Create safe member identifier
#====

narclim_member_id <- function(
                gcm,
                rcm,
                scenario,
                variable
) {
        
        safe_name(
                paste(
                        variable,
                        scenario,
                        gcm,
                        rcm,
                        sep = "_"
                )
        )
}

#====
# Build one mean raster for one NARCliM member row
#====

process_member_mean <- function(
                member_row,
                start_year,
                end_year,
                city_shape = NULL,
                test_mode = FALSE,
                sample_layers = NULL,
                output_dir
) {
        
        if (
                is.null(member_row) ||
                nrow(member_row) != 1L
        ) {
                stop("process_member_mean() requires exactly one member row.")
        }
        
        r <- load_and_crop_rasters(
                files = member_row$file,
                var_token = member_row$variable,
                start_year = start_year,
                end_year = end_year,
                test_mode = test_mode,
                city_shape = city_shape,
                sample_layers = sample_layers
        )
        
        if (is.null(r)) {
                return(NULL)
        }
        
        r_mean <- deep_mean(r)
        
        if (is.null(r_mean)) {
                return(NULL)
        }
        
        member_id <- narclim_member_id(
                gcm = member_row$gcm,
                rcm = member_row$rcm,
                scenario = member_row$scenario,
                variable = member_row$variable
        )
        
        path <- file.path(
                output_dir,
                paste0(
                        "member_mean_",
                        member_id,
                        ".tif"
                )
        )
        
        safe_write_raster(
                r_mean,
                path
        )
}

split_member_table <- function(member_table) {
        
        if (
                is.null(member_table) ||
                nrow(member_table) == 0L
        ) {
                return(list())
        }
        
        split(
                member_table,
                seq_len(nrow(member_table))
        )
}
###############################################################################
# NARCliM2 REMOTE SOURCE CONFIGURATION
###############################################################################

narclim_source_config <- function(
                product = c(
                        "raw",
                        "bias_adjusted"
                )
) {
        
        product <- match.arg(product)
        
        if (product == "raw") {
                
                return(
                        list(
                                product = "raw",
                                
                                catalogue_base = paste0(
                                        "https://thredds.nci.org.au/thredds/catalog/",
                                        "zz63/NARCliM2-0"
                                ),
                                
                                opendap_base = paste0(
                                        "https://thredds.nci.org.au/thredds/dodsC/"
                                ),
                                
                                mip_era = "output-CMIP6",
                                activity_id = "DD",
                                institution_id = "NSW-Government",
                                version_realisation = "v1-r1",
                                version = "latest"
                        )
                )
        }
        
        if (product == "bias_adjusted") {
                
                return(
                        list(
                                product = "bias_adjusted",
                                
                                catalogue_base = paste0(
                                        "https://thredds.nci.org.au/thredds/catalog/",
                                        "zz63/NARCliM2-0-derived/output-CMIP6/",
                                        "bias-adjusted-output"
                                ),
                                
                                opendap_base = paste0(
                                        "https://thredds.nci.org.au/thredds/dodsC/",
                                        "zz63/NARCliM2-0-derived/output-CMIP6/",
                                        "bias-adjusted-output"
                                ),
                                
                                institution_id = "NSW-Government",
                                version_realisation =
                                        "v1-r1-NSWGovernment-CDF-AGCDv1-1990-2009",
                                version = "latest"
                        )
                )
        }
        
        stop(
                "Unsupported NARCliM product: ",
                product
        )
}

###############################################################################
# NARCliM2 GCM / variant configuration
###############################################################################

narclim_model_table <- function() {
        
        data.frame(
                gcm = c(
                        "ACCESS-ESM1-5",
                        "EC-Earth3-Veg",
                        "MPI-ESM1-2-HR",
                        "NorESM2-MM",
                        "UKESM1-0-LL"
                ),
                
                variant = c(
                        "r6i1p1f1",
                        "r1i1p1f1",
                        "r1i1p1f1",
                        "r1i1p1f1",
                        "r1i1p1f2"
                ),
                
                stringsAsFactors = FALSE
        )
}
###############################################################################
# Discover NARCliM2 remote resources
#
# Returns one row per available annual NetCDF resource matching the requested:
# product, variable, domain, scenario, GCM, RCM, frequency and year range.
###############################################################################

discover_narclim_resources <- function(
                product = c("raw", "bias_adjusted"),
                variable,
                domain,
                scenario,
                gcm,
                rcm,
                frequency = "day",
                start_year = NULL,
                end_year = NULL
) {
        
        product <- match.arg(product)
        
        cfg <- narclim_source_config(
                product = product
        )
        
        model_tbl <- narclim_model_table()
        
        model_match <- model_tbl[
                model_tbl$gcm == gcm,
                ,
                drop = FALSE
        ]
        
        if (nrow(model_match) != 1L) {
                stop(
                        "Could not resolve unique NARCliM variant for GCM: ",
                        gcm
                )
        }
        
        variant <- model_match$variant[[1]]
        
        if (!rcm %in% c(
                "NARCliM2-0-WRF412R3",
                "NARCliM2-0-WRF412R5"
        )) {
                stop(
                        "Unsupported NARCliM RCM: ",
                        rcm
                )
        }
        
        # ----------------------------------------------------------
        # Build catalogue path
        # ----------------------------------------------------------
        
        if (product == "raw") {
                
                catalogue_url <- paste(
                        cfg$catalogue_base,
                        cfg$mip_era,
                        cfg$activity_id,
                        if (identical(domain, "SEAus-04"))
                                "NARCliM2-0-SEAus-04"
                        else
                                domain,
                        cfg$institution_id,
                        gcm,
                        scenario,
                        variant,
                        rcm,
                        cfg$version_realisation,
                        frequency,
                        variable,
                        cfg$version,
                        "catalog.xml",
                        sep = "/"
                )
                
        } else {
                
                catalogue_url <- paste(
                        cfg$catalogue_base,
                        domain,
                        cfg$institution_id,
                        gcm,
                        scenario,
                        variant,
                        rcm,
                        cfg$version_realisation,
                        frequency,
                        variable,
                        cfg$version,
                        "catalog.xml",
                        sep = "/"
                )
        }
        
        cat(
                "\nTesting NARCliM catalogue:\n",
                catalogue_url,
                "\n",
                sep = ""
        )
        
        catalogue <- tryCatch(
                thredds::get_catalog(
                        catalogue_url
                ),
                error = function(e) {
                        
                        message(
                                "Catalogue unavailable: ",
                                conditionMessage(e)
                        )
                        
                        return(NULL)
                }
        )
        
        if (is.null(catalogue)) {
                
                return(
                        data.frame(
                                product = product,
                                domain = domain,
                                variable = variable,
                                scenario = scenario,
                                gcm = gcm,
                                variant = variant,
                                rcm = rcm,
                                frequency = frequency,
                                year = integer(0),
                                filename = character(0),
                                dataset_path = character(0),
                                opendap_url = character(0),
                                stringsAsFactors = FALSE
                        )
                )
        }
        
        dataset_names <- catalogue$get_dataset_names()
        datasets <- catalogue$get_datasets()
        
        if (
                length(dataset_names) == 0L ||
                length(datasets) == 0L
        ) {
                
                return(
                        data.frame(
                                product = product,
                                domain = domain,
                                variable = variable,
                                scenario = scenario,
                                gcm = gcm,
                                variant = variant,
                                rcm = rcm,
                                frequency = frequency,
                                year = integer(0),
                                filename = character(0),
                                dataset_path = character(0),
                                opendap_url = character(0),
                                stringsAsFactors = FALSE
                        )
                )
        }
        
        rows <- lapply(
                seq_along(datasets),
                function(i) {
                        
                        ds <- datasets[[i]]
                        
                        dataset_path <- NULL
                        
                        if (!is.null(ds$url)) {
                                dataset_path <- ds$url
                        }
                        
                        if (
                                is.null(dataset_path) &&
                                !is.null(ds$urlPath)
                        ) {
                                dataset_path <- ds$urlPath
                        }
                        
                        if (
                                is.null(dataset_path) ||
                                length(dataset_path) == 0L ||
                                is.na(dataset_path[[1]])
                        ) {
                                return(NULL)
                        }
                        
                        dataset_path <- as.character(
                                dataset_path[[1]]
                        )
                        
                        filename <- basename(
                                dataset_path
                        )
                        
                        date_match <- stringr::str_match(
                                filename,
                                "_([0-9]{8})-([0-9]{8})\\.nc$"
                        )
                        
                        year <- suppressWarnings(
                                as.integer(
                                        substr(
                                                date_match[, 2],
                                                1,
                                                4
                                        )
                                )
                        )
                        
                        # ------------------------------------------
                        # Convert catalogue/file path to OPeNDAP
                        # ------------------------------------------
                        
                        clean_path <- sub(
                                "^https?://[^/]+/thredds/(dodsC|fileServer)/",
                                "",
                                dataset_path
                        )
                        
                        clean_path <- sub(
                                "^/+",
                                "",
                                clean_path
                        )
                        
                        if (product == "raw") {
                                
                                opendap_url <- paste0(
                                        sub("/$", "", cfg$opendap_base),
                                        "/",
                                        clean_path
                                )
                                
                        } else {
                                
                                expected_prefix <- paste0(
                                        "zz63/NARCliM2-0-derived/output-CMIP6/",
                                        "bias-adjusted-output/"
                                )
                                
                                if (
                                        startsWith(
                                                clean_path,
                                                expected_prefix
                                        )
                                ) {
                                        
                                        opendap_url <- paste0(
                                                "https://thredds.nci.org.au/thredds/dodsC/",
                                                clean_path
                                        )
                                        
                                } else {
                                        
                                        opendap_url <- paste0(
                                                sub("/$", "", cfg$opendap_base),
                                                "/",
                                                clean_path
                                        )
                                }
                        }
                        
                        data.frame(
                                product = product,
                                domain = domain,
                                variable = variable,
                                scenario = scenario,
                                gcm = gcm,
                                variant = variant,
                                rcm = rcm,
                                frequency = frequency,
                                year = year,
                                filename = filename,
                                dataset_path = dataset_path,
                                opendap_url = opendap_url,
                                stringsAsFactors = FALSE
                        )
                }
        )
        
        rows <- rows[
                !vapply(
                        rows,
                        is.null,
                        logical(1)
                )
        ]
        
        if (length(rows) == 0L) {
                
                return(
                        data.frame(
                                product = product,
                                domain = domain,
                                variable = variable,
                                scenario = scenario,
                                gcm = gcm,
                                variant = variant,
                                rcm = rcm,
                                frequency = frequency,
                                year = integer(0),
                                filename = character(0),
                                dataset_path = character(0),
                                opendap_url = character(0),
                                stringsAsFactors = FALSE
                        )
                )
        }
        
        out <- do.call(
                rbind,
                rows
        )
        
        # ----------------------------------------------------------
        # Restrict to user-requested years
        # ----------------------------------------------------------
        
        if (!is.null(start_year)) {
                
                out <- out[
                        !is.na(out$year) &
                                out$year >= start_year,
                        ,
                        drop = FALSE
                ]
        }
        
        if (!is.null(end_year)) {
                
                out <- out[
                        !is.na(out$year) &
                                out$year <= end_year,
                        ,
                        drop = FALSE
                ]
        }
        
        out <- out[
                order(out$year),
                ,
                drop = FALSE
        ]
        
        rownames(out) <- NULL
        
        out
}

###############################################################################
# Inspect one NARCliM OPeNDAP resource
#
# Reads metadata and coordinates only.
# Does NOT read the full climate variable.
###############################################################################

inspect_narclim_resource <- function(
                opendap_url,
                variable
) {
        
        nc <- ncdf4::nc_open(
                opendap_url,
                readunlim = FALSE
        )
        
        on.exit(
                ncdf4::nc_close(nc),
                add = TRUE
        )
        
        if (!variable %in% names(nc$var)) {
                stop(
                        "Variable '",
                        variable,
                        "' not found in remote NARCliM resource."
                )
        }
        
        required_coord_vars <- c(
                "lon",
                "lat"
        )
        
        missing_coords <- setdiff(
                required_coord_vars,
                names(nc$var)
        )
        
        if (length(missing_coords) > 0L) {
                stop(
                        "Remote NARCliM resource is missing coordinate variables: ",
                        paste(
                                missing_coords,
                                collapse = ", "
                        )
                )
        }
        
        if (!all(c("rlon", "rlat", "time") %in% names(nc$dim))) {
                stop(
                        "Expected rlon, rlat and time dimensions were not all found."
                )
        }
        
        dimension_order <- vapply(
                nc$var[[variable]]$dim,
                function(x) x$name,
                character(1)
        )
        
        lon2d <- ncdf4::ncvar_get(
                nc,
                "lon"
        )
        
        lat2d <- ncdf4::ncvar_get(
                nc,
                "lat"
        )
        
        rlon <- nc$dim$rlon$vals
        rlat <- nc$dim$rlat$vals
        time_values <- nc$dim$time$vals
        
        time_units <- tryCatch(
                ncdf4::ncatt_get(
                        nc,
                        "time",
                        "units"
                )$value,
                error = function(e) NA_character_
        )
        
        calendar <- tryCatch(
                ncdf4::ncatt_get(
                        nc,
                        "time",
                        "calendar"
                )$value,
                error = function(e) "standard"
        )
        
        variable_units <- tryCatch(
                ncdf4::ncatt_get(
                        nc,
                        variable,
                        "units"
                )$value,
                error = function(e) NA_character_
        )
        
        list(
                variable = variable,
                dimension_order = dimension_order,
                
                n_rlon = length(rlon),
                n_rlat = length(rlat),
                n_time = length(time_values),
                
                rlon = rlon,
                rlat = rlat,
                lon2d = lon2d,
                lat2d = lat2d,
                
                time_values = time_values,
                time_units = time_units,
                calendar = calendar,
                
                variable_units = variable_units
        )
}

###############################################################################
# Convert a WGS84 bounding box to NARCliM native-grid index bounds
###############################################################################

narclim_bbox_to_indices <- function(
                lon2d,
                lat2d,
                lon_min,
                lon_max,
                lat_min,
                lat_max
) {
        
        if (!identical(
                dim(lon2d),
                dim(lat2d)
        )) {
                stop(
                        "lon2d and lat2d have different dimensions."
                )
        }
        
        inside <- which(
                lon2d >= lon_min &
                        lon2d <= lon_max &
                        lat2d >= lat_min &
                        lat2d <= lat_max,
                arr.ind = TRUE
        )
        
        if (nrow(inside) == 0L) {
                stop(
                        "No NARCliM grid cells intersect the requested bounding box."
                )
        }
        
        x_min <- min(inside[, 1])
        x_max <- max(inside[, 1])
        
        y_min <- min(inside[, 2])
        y_max <- max(inside[, 2])
        
        list(
                x_start = x_min,
                x_count = x_max - x_min + 1L,
                
                y_start = y_min,
                y_count = y_max - y_min + 1L,
                
                x_min = x_min,
                x_max = x_max,
                y_min = y_min,
                y_max = y_max,
                
                n_matching_cells = nrow(inside)
        )
}

###############################################################################
# Stream a small NARCliM native-grid subset from one OPeNDAP resource
#
# This reads only:
#   selected x indices
#   selected y indices
#   selected time indices
#
# No full-file download.
###############################################################################

stream_narclim_subset <- function(
                opendap_url,
                variable,
                indices,
                t_start = 1L,
                t_count = 3L
) {
        
        nc <- ncdf4::nc_open(
                opendap_url,
                readunlim = FALSE
        )
        
        on.exit(
                ncdf4::nc_close(nc),
                add = TRUE
        )
        
        if (!variable %in% names(nc$var)) {
                stop(
                        "Variable '",
                        variable,
                        "' not found in remote resource."
                )
        }
        
        dims <- vapply(
                nc$var[[variable]]$dim,
                function(x) x$name,
                character(1)
        )
        
        if (!identical(
                dims,
                c("rlon", "rlat", "time")
        )) {
                stop(
                        "Unexpected variable dimension order: ",
                        paste(
                                dims,
                                collapse = ", "
                        )
                )
        }
        
        n_time <- length(
                nc$dim$time$vals
        )
        
        if (
                t_start < 1L ||
                t_start > n_time
        ) {
                stop(
                        "t_start is outside the resource time dimension."
                )
        }
        
        if (t_count == -1L) {
                t_count <- n_time - t_start + 1L
        }
        
        t_count <- min(
                as.integer(t_count),
                n_time - t_start + 1L
        )
        
        vals <- ncdf4::ncvar_get(
                nc,
                variable,
                start = c(
                        indices$x_start,
                        indices$y_start,
                        t_start
                ),
                count = c(
                        indices$x_count,
                        indices$y_count,
                        t_count
                ),
                collapse_degen = FALSE
        )
        
        rlon <- ncdf4::ncvar_get(
                nc,
                "rlon",
                start = indices$x_start,
                count = indices$x_count
        )
        
        rlat <- ncdf4::ncvar_get(
                nc,
                "rlat",
                start = indices$y_start,
                count = indices$y_count
        )
        
        lon2d <- ncdf4::ncvar_get(
                nc,
                "lon",
                start = c(
                        indices$x_start,
                        indices$y_start
                ),
                count = c(
                        indices$x_count,
                        indices$y_count
                )
        )
        
        lat2d <- ncdf4::ncvar_get(
                nc,
                "lat",
                start = c(
                        indices$x_start,
                        indices$y_start
                ),
                count = c(
                        indices$x_count,
                        indices$y_count
                )
        )
        
        time_values <- ncdf4::ncvar_get(
                nc,
                "time",
                start = t_start,
                count = t_count
        )
        
        time_units <- ncdf4::ncatt_get(
                nc,
                "time",
                "units"
        )$value
        
        calendar <- ncdf4::ncatt_get(
                nc,
                "time",
                "calendar"
        )$value
        
        time_table <- parse_cf_time(
                values = time_values,
                units = time_units,
                calendar = calendar
        )
        
        variable_units <- ncdf4::ncatt_get(
                nc,
                variable,
                "units"
        )$value
        
        list(
                values = vals,
                rlon = rlon,
                rlat = rlat,
                lon2d = lon2d,
                lat2d = lat2d,
                time_values = time_values,
                time_units = time_units,
                calendar = calendar,
                time_table = time_table,
                variable_units = variable_units,
                dimension_order = dims
        )
}

###############################################################################
# Calendar-aware day of year
###############################################################################

calendar_day_of_year <- function(
                year,
                month,
                day,
                calendar
) {
        
        calendar <- tolower(
                trimws(
                        as.character(calendar)[1L]
                )
        )
        
        if (
                is.na(calendar) ||
                !nzchar(calendar)
        ) {
                stop(
                        "Missing calendar in calendar_day_of_year()."
                )
        }
        
        if (identical(
                calendar,
                "360_day"
        )) {
                
                return(
                        as.integer(
                                (month - 1L) * 30L + day
                        )
                )
        }
        
        if (
                calendar %in%
                c(
                        "noleap",
                        "365_day"
                )
        ) {
                
                month_lengths <- c(
                        31L, 28L, 31L, 30L,
                        31L, 30L, 31L, 31L,
                        30L, 31L, 30L, 31L
                )
                
                return(
                        as.integer(
                                vapply(
                                        seq_along(month),
                                        function(i) {
                                                
                                                if (month[i] == 1L) {
                                                        
                                                        day[i]
                                                        
                                                } else {
                                                        
                                                        sum(
                                                                month_lengths[
                                                                        seq_len(
                                                                                month[i] - 1L
                                                                        )
                                                                ]
                                                        ) + day[i]
                                                }
                                        },
                                        integer(1)
                                )
                        )
                )
        }
        
        date_value <- as.Date(
                sprintf(
                        "%04d-%02d-%02d",
                        year,
                        month,
                        day
                )
        )
        
        as.integer(
                format(
                        date_value,
                        "%j"
                )
        )
}


###############################################################################
# Parse CF-compliant NetCDF time values
###############################################################################

parse_cf_time <- function(
                values,
                units,
                calendar = "standard"
) {
        
        calendar <- tolower(
                trimws(
                        as.character(calendar)[1L]
                )
        )
        
        units <- trimws(
                as.character(units)[1L]
        )
        
        if (
                is.na(calendar) ||
                !nzchar(calendar)
        ) {
                calendar <- "standard"
        }
        
        if (identical(
                calendar,
                "365_day"
        )) {
                calendar <- "noleap"
        }
        
        supported <- c(
                "standard",
                "gregorian",
                "proleptic_gregorian",
                "noleap",
                "360_day"
        )
        
        if (!calendar %in% supported) {
                stop(
                        "Unsupported calendar: ",
                        calendar
                )
        }
        
        if (
                is.na(units) ||
                !nzchar(units)
        ) {
                stop(
                        "Missing CF time units."
                )
        }
        
        values <- as.numeric(
                values
        )
        
        if (
                length(values) == 0L ||
                anyNA(values) ||
                any(!is.finite(values))
        ) {
                stop(
                        "CF time coordinate contains missing or non-finite values."
                )
        }
        
        cf <- CFtime::CFtime(
                units,
                calendar,
                values
        )
        
        cf_timestamp <- CFtime::as_timestamp(
                cf,
                format = "date",
                asPOSIX = FALSE
        )
        
        cf_date <- substr(
                as.character(cf_timestamp),
                1L,
                10L
        )
        
        year <- as.integer(
                substr(
                        cf_date,
                        1L,
                        4L
                )
        )
        
        month <- as.integer(
                substr(
                        cf_date,
                        6L,
                        7L
                )
        )
        
        day <- as.integer(
                substr(
                        cf_date,
                        9L,
                        10L
                )
        )
        
        data.frame(
                resource_time_index =
                        seq_along(values),
                
                time_value =
                        values,
                
                calendar =
                        calendar,
                
                cf_date =
                        cf_date,
                
                calendar_year =
                        year,
                
                calendar_month =
                        month,
                
                calendar_day =
                        day,
                
                calendar_yday =
                        calendar_day_of_year(
                                year,
                                month,
                                day,
                                calendar
                        ),
                
                stringsAsFactors = FALSE
        )
}

###############################################################################
# Stream NARCliM resources over a selected period and calculate
# a cell-by-cell temporal mean without retaining the full daily cube.
#
# resources:
#   output from discover_narclim_resources()
#
# indices:
#   output from narclim_bbox_to_indices()
#
# Returns:
#   mean_values     matrix [x, y]
#   count_values    matrix [x, y]
#   metadata
###############################################################################

stream_narclim_period_mean <- function(
                resources,
                variable,
                indices,
                start_date,
                end_date,
                chunk_days = 31L
) {
        
        if (
                is.null(resources) ||
                nrow(resources) == 0L
        ) {
                stop(
                        "No NARCliM resources supplied."
                )
        }
        
        start_date <- as.Date(start_date)
        end_date   <- as.Date(end_date)
        
        if (
                is.na(start_date) ||
                is.na(end_date) ||
                start_date > end_date
        ) {
                stop(
                        "Invalid start_date/end_date."
                )
        }
        
        chunk_days <- as.integer(chunk_days)
        
        if (
                is.na(chunk_days) ||
                chunk_days < 1L
        ) {
                stop(
                        "chunk_days must be >= 1."
                )
        }
        
        # ---------------------------------------------------------------------
        # Running accumulators
        # ---------------------------------------------------------------------
        
        running_sum <- matrix(
                0,
                nrow = indices$x_count,
                ncol = indices$y_count
        )
        
        running_count <- matrix(
                0L,
                nrow = indices$x_count,
                ncol = indices$y_count
        )
        
        total_selected_days <- 0L
        
        observed_units <- character(0)
        observed_calendars <- character(0)
        
        first_selected_date <- NULL
        last_selected_date <- NULL
        
        # ---------------------------------------------------------------------
        # Process annual resources
        # ---------------------------------------------------------------------
        
        resources <- resources[
                order(resources$year),
                ,
                drop = FALSE
        ]
        
        for (resource_i in seq_len(nrow(resources))) {
                
                resource_url <- resources$opendap_url[[resource_i]]
                
                cat(
                        "\nOpening resource ",
                        resource_i,
                        " of ",
                        nrow(resources),
                        ":\n",
                        basename(resources$filename[[resource_i]]),
                        "\n",
                        sep = ""
                )
                
                info <- inspect_narclim_resource(
                        opendap_url = resource_url,
                        variable = variable
                )
                
                full_time_table <- parse_cf_time(
                        values = info$time_values,
                        units = info$time_units,
                        calendar = info$calendar
                )
                
                # -------------------------------------------------------------
                # Date filtering
                #
                # For standard/Gregorian dates this can use as.Date safely.
                # We will generalise the range comparison for 360_day/noleap
                # after validating the standard-calendar path.
                # -------------------------------------------------------------
                
                if (
                        !info$calendar %in%
                        c(
                                "standard",
                                "gregorian",
                                "proleptic_gregorian"
                        )
                ) {
                        stop(
                                paste0(
                                        "Initial period-mean implementation ",
                                        "currently accepts standard/Gregorian ",
                                        "calendar resources only. Observed: ",
                                        info$calendar
                                )
                        )
                }
                
                resource_dates <- as.Date(
                        full_time_table$cf_date
                )
                
                selected_index <- which(
                        resource_dates >= start_date &
                                resource_dates <= end_date
                )
                
                if (length(selected_index) == 0L) {
                        next
                }
                
                # The selected part of an annual resource should be contiguous
                # for one requested date range.
                first_index <- min(selected_index)
                last_index  <- max(selected_index)
                
                if (
                        !identical(
                                selected_index,
                                seq.int(
                                        first_index,
                                        last_index
                                )
                        )
                ) {
                        stop(
                                "Selected time indices are not contiguous."
                        )
                }
                
                selected_dates <- resource_dates[
                        selected_index
                ]
                
                if (is.null(first_selected_date)) {
                        first_selected_date <- min(
                                selected_dates
                        )
                }
                
                last_selected_date <- max(
                        selected_dates
                )
                
                observed_units <- unique(
                        c(
                                observed_units,
                                info$variable_units
                        )
                )
                
                observed_calendars <- unique(
                        c(
                                observed_calendars,
                                info$calendar
                        )
                )
                
                # -------------------------------------------------------------
                # Chunk only the selected time indices
                # -------------------------------------------------------------
                
                chunk_starts <- seq.int(
                        first_index,
                        last_index,
                        by = chunk_days
                )
                
                for (
                        chunk_i in seq_along(
                                chunk_starts
                        )
                ) {
                        
                        t_start <- chunk_starts[[chunk_i]]
                        
                        t_count <- min(
                                chunk_days,
                                last_index - t_start + 1L
                        )
                        
                        cat(
                                "  streaming time indices ",
                                t_start,
                                "-",
                                t_start + t_count - 1L,
                                "\n",
                                sep = ""
                        )
                        
                        chunk <- stream_narclim_subset(
                                opendap_url = resource_url,
                                variable = variable,
                                indices = indices,
                                t_start = t_start,
                                t_count = t_count
                        )
                        
                        vals <- chunk$values
                        
                        if (
                                !identical(
                                        dim(vals),
                                        c(
                                                indices$x_count,
                                                indices$y_count,
                                                t_count
                                        )
                                )
                        ) {
                                stop(
                                        "Unexpected streamed chunk dimensions."
                                )
                        }
                        
                        # -----------------------------------------------------
                        # Cell-by-cell running sum/count
                        # -----------------------------------------------------
                        
                        valid <- is.finite(
                                vals
                        )
                        
                        chunk_sum <- apply(
                                vals,
                                c(1, 2),
                                function(x) {
                                        sum(
                                                x,
                                                na.rm = TRUE
                                        )
                                }
                        )
                        
                        chunk_count <- apply(
                                valid,
                                c(1, 2),
                                sum
                        )
                        
                        running_sum <- running_sum +
                                chunk_sum
                        
                        running_count <- running_count +
                                chunk_count
                        
                        total_selected_days <-
                                total_selected_days +
                                t_count
                }
        }
        
        if (total_selected_days == 0L) {
                stop(
                        "No time steps matched the requested period."
                )
        }
        
        mean_values <- running_sum /
                running_count
        
        mean_values[
                running_count == 0L
        ] <- NA_real_
        
        list(
                mean_values = mean_values,
                count_values = running_count,
                
                variable = variable,
                units = observed_units,
                calendars = observed_calendars,
                
                start_date =
                        as.character(first_selected_date),
                
                end_date =
                        as.character(last_selected_date),
                
                selected_time_steps =
                        total_selected_days,
                
                x_count =
                        indices$x_count,
                
                y_count =
                        indices$y_count,
                
                chunk_days =
                        chunk_days
        )
}


###############################################################################
# Read and select arbitrary polygon geography
#
# Works with:
#   LGA
#   SA2
#   GCCSA
#   State
#   custom polygon layers
#
# Returned polygons always contain the standard fields:
#   GEO_ID
#   GEO_NAME
###############################################################################

read_geography_boundaries <- function(
                boundary_path,
                boundary_layer = NULL,
                id_column,
                name_column,
                filter_column = NULL,
                filter_values = NULL,
                target_crs = 4326
) {
        
        if (
                is.null(boundary_path) ||
                !file.exists(boundary_path)
        ) {
                stop(
                        "Boundary file does not exist:\n",
                        boundary_path
                )
        }
        
        # ---------------------------------------------------------------------
        # Read spatial data
        # ---------------------------------------------------------------------
        
        if (
                is.null(boundary_layer) ||
                !nzchar(boundary_layer)
        ) {
                
                boundaries <- sf::st_read(
                        boundary_path,
                        quiet = TRUE
                )
                
        } else {
                
                boundaries <- sf::st_read(
                        boundary_path,
                        layer = boundary_layer,
                        quiet = TRUE
                )
        }
        
        if (nrow(boundaries) == 0L) {
                stop(
                        "Boundary file contains no features."
                )
        }
        
        # ---------------------------------------------------------------------
        # Validate required columns
        # ---------------------------------------------------------------------
        
        required_columns <- c(
                id_column,
                name_column
        )
        
        if (!is.null(filter_column)) {
                required_columns <- unique(
                        c(
                                required_columns,
                                filter_column
                        )
                )
        }
        
        missing_columns <- setdiff(
                required_columns,
                names(boundaries)
        )
        
        if (length(missing_columns) > 0L) {
                
                stop(
                        paste0(
                                "Boundary file is missing required column(s): ",
                                paste(
                                        missing_columns,
                                        collapse = ", "
                                ),
                                "\nAvailable columns:\n",
                                paste(
                                        names(boundaries),
                                        collapse = ", "
                                )
                        )
                )
        }
        
        # ---------------------------------------------------------------------
        # Optional user selection
        # ---------------------------------------------------------------------
        
        if (
                !is.null(filter_column) &&
                !is.null(filter_values)
        ) {
                
                filter_values <- as.character(
                        filter_values
                )
                
                boundaries <- boundaries[
                        as.character(
                                boundaries[[filter_column]]
                        ) %in% filter_values,
                        ,
                        drop = FALSE
                ]
        }
        
        if (nrow(boundaries) == 0L) {
                
                stop(
                        paste0(
                                "No polygons matched the requested boundary selection."
                        )
                )
        }
        
        # ---------------------------------------------------------------------
        # Standardise identifiers
        # ---------------------------------------------------------------------
        
        boundaries$GEO_ID <- as.character(
                boundaries[[id_column]]
        )
        
        boundaries$GEO_NAME <- as.character(
                boundaries[[name_column]]
        )
        
        if (
                anyNA(boundaries$GEO_ID) ||
                any(!nzchar(boundaries$GEO_ID))
        ) {
                stop(
                        "Selected geography contains missing GEO_ID values."
                )
        }
        
        # ---------------------------------------------------------------------
        # Geometry QA
        # ---------------------------------------------------------------------
        
        boundaries <- sf::st_make_valid(
                boundaries
        )
        
        empty_geometry <- sf::st_is_empty(
                boundaries
        )
        
        if (any(empty_geometry)) {
                
                stop(
                        "One or more selected geography polygons have empty geometry."
                )
        }
        
        # ---------------------------------------------------------------------
        # NARCliM lon/lat arrays are WGS84
        # ---------------------------------------------------------------------
        
        boundaries_wgs84 <- sf::st_transform(
                boundaries,
                target_crs
        )
        
        boundaries_wgs84
}

###############################################################################
# Convert arbitrary geography polygons to a NARCliM streaming window
#
# IMPORTANT:
# The bounding rectangle is used ONLY to minimise the OPeNDAP read.
# It does not define the analytical geography.
#
# Exact polygon/native-cell intersection is performed later.
###############################################################################

geography_to_narclim_indices <- function(
                geography_sf,
                lon2d,
                lat2d,
                padding_cells = 1L
) {
        
        if (
                is.null(geography_sf) ||
                nrow(geography_sf) == 0L
        ) {
                stop(
                        "No geography polygons supplied."
                )
        }
        
        geography_wgs84 <- sf::st_transform(
                geography_sf,
                4326
        )
        
        bb <- sf::st_bbox(
                geography_wgs84
        )
        
        indices <- narclim_bbox_to_indices(
                lon2d = lon2d,
                lat2d = lat2d,
                
                lon_min = as.numeric(bb["xmin"]),
                lon_max = as.numeric(bb["xmax"]),
                lat_min = as.numeric(bb["ymin"]),
                lat_max = as.numeric(bb["ymax"])
        )
        
        # ---------------------------------------------------------------------
        # One-cell padding around the bounding window.
        #
        # Useful because later exact polygon/native-cell intersection can
        # legitimately involve cells whose centre lies just outside the polygon
        # bounding extent.
        # ---------------------------------------------------------------------
        
        padding_cells <- as.integer(
                padding_cells
        )
        
        if (
                is.na(padding_cells) ||
                padding_cells < 0L
        ) {
                stop(
                        "padding_cells must be >= 0."
                )
        }
        
        x_min <- max(
                1L,
                indices$x_min - padding_cells
        )
        
        x_max <- min(
                nrow(lon2d),
                indices$x_max + padding_cells
        )
        
        y_min <- max(
                1L,
                indices$y_min - padding_cells
        )
        
        y_max <- min(
                ncol(lon2d),
                indices$y_max + padding_cells
        )
        
        list(
                x_start = x_min,
                x_count = x_max - x_min + 1L,
                
                y_start = y_min,
                y_count = y_max - y_min + 1L,
                
                x_min = x_min,
                x_max = x_max,
                
                y_min = y_min,
                y_max = y_max,
                
                bounding_box = c(
                        xmin = as.numeric(bb["xmin"]),
                        xmax = as.numeric(bb["xmax"]),
                        ymin = as.numeric(bb["ymin"]),
                        ymax = as.numeric(bb["ymax"])
                ),
                
                geography_count =
                        nrow(geography_sf),
                
                geography_ids =
                        geography_sf$GEO_ID,
                
                geography_names =
                        geography_sf$GEO_NAME,
                
                padding_cells =
                        padding_cells
        )
}

###############################################################################
# Build native NARCliM cell polygons for a streamed subset
#
# Uses the authoritative 2-D WGS84 lon/lat cell centres.
# Cell boundaries are approximated from neighbouring cell centres.
#
# Output:
#   one sf polygon per native NARCliM cell
#   with native x/y indices preserved
###############################################################################

build_narclim_cell_polygons <- function(
                lon2d,
                lat2d,
                indices
) {
        
        if (!identical(dim(lon2d), dim(lat2d))) {
                stop("lon2d and lat2d dimensions do not match.")
        }
        
        x_seq <- indices$x_min:indices$x_max
        y_seq <- indices$y_min:indices$y_max
        
        centres <- expand.grid(
                x_index = x_seq,
                y_index = y_seq
        )
        
        centres$lon <- mapply(
                function(x, y) lon2d[x, y],
                centres$x_index,
                centres$y_index
        )
        
        centres$lat <- mapply(
                function(x, y) lat2d[x, y],
                centres$x_index,
                centres$y_index
        )
        
        if (
                anyNA(centres$lon) ||
                anyNA(centres$lat)
        ) {
                stop("Missing NARCliM cell-centre coordinates.")
        }
        
        # ----------------------------------------------------------
        # Estimate cell edges from midpoint spacing
        # ----------------------------------------------------------
        
        x_unique <- sort(unique(centres$x_index))
        y_unique <- sort(unique(centres$y_index))
        
        polygons <- vector(
                "list",
                nrow(centres)
        )
        
        for (i in seq_len(nrow(centres))) {
                
                x <- centres$x_index[i]
                y <- centres$y_index[i]
                
                lon_c <- lon2d[x, y]
                lat_c <- lat2d[x, y]
                
                # neighbouring centres in x direction
                lon_left <- if (x > 1L) {
                        lon2d[x - 1L, y]
                } else {
                        NA_real_
                }
                
                lon_right <- if (x < nrow(lon2d)) {
                        lon2d[x + 1L, y]
                } else {
                        NA_real_
                }
                
                # neighbouring centres in y direction
                lat_down <- if (y > 1L) {
                        lat2d[x, y - 1L]
                } else {
                        NA_real_
                }
                
                lat_up <- if (y < ncol(lat2d)) {
                        lat2d[x, y + 1L]
                } else {
                        NA_real_
                }
                
                lon_min <- if (is.finite(lon_left)) {
                        (lon_left + lon_c) / 2
                } else if (is.finite(lon_right)) {
                        lon_c - (lon_right - lon_c) / 2
                } else {
                        stop("Cannot determine longitude cell bounds.")
                }
                
                lon_max <- if (is.finite(lon_right)) {
                        (lon_c + lon_right) / 2
                } else if (is.finite(lon_left)) {
                        lon_c + (lon_c - lon_left) / 2
                } else {
                        stop("Cannot determine longitude cell bounds.")
                }
                
                lat_min <- if (is.finite(lat_down)) {
                        (lat_down + lat_c) / 2
                } else if (is.finite(lat_up)) {
                        lat_c - (lat_up - lat_c) / 2
                } else {
                        stop("Cannot determine latitude cell bounds.")
                }
                
                lat_max <- if (is.finite(lat_up)) {
                        (lat_c + lat_up) / 2
                } else if (is.finite(lat_down)) {
                        lat_c + (lat_c - lat_down) / 2
                } else {
                        stop("Cannot determine latitude cell bounds.")
                }
                
                xy <- matrix(
                        c(
                                lon_min, lat_min,
                                lon_max, lat_min,
                                lon_max, lat_max,
                                lon_min, lat_max,
                                lon_min, lat_min
                        ),
                        ncol = 2,
                        byrow = TRUE
                )
                
                polygons[[i]] <- sf::st_polygon(
                        list(xy)
                )
        }
        
        out <- sf::st_sf(
                x_index = centres$x_index,
                y_index = centres$y_index,
                cell_id = paste(
                        centres$x_index,
                        centres$y_index,
                        sep = "_"
                ),
                lon = centres$lon,
                lat = centres$lat,
                geometry = sf::st_sfc(
                        polygons,
                        crs = 4326
                )
        )
        
        sf::st_make_valid(out)
}

###############################################################################
# Intersect native NARCliM cells with arbitrary selected geography
#
# Works for:
#   LGA
#   SA2
#   GCCSA
#   State
#   custom polygons
###############################################################################

intersect_narclim_cells_geography <- function(
                narclim_cells,
                geography_sf
) {
        
        if (
                is.null(narclim_cells) ||
                nrow(narclim_cells) == 0L
        ) {
                stop("No NARCliM cells supplied.")
        }
        
        if (
                is.null(geography_sf) ||
                nrow(geography_sf) == 0L
        ) {
                stop("No geography polygons supplied.")
        }
        
        geography_sf <- sf::st_transform(
                geography_sf,
                sf::st_crs(narclim_cells)
        )
        
        geography_keep <- geography_sf[
                ,
                c(
                        "GEO_ID",
                        "GEO_NAME"
                )
        ]
        
        intersections <- suppressWarnings(
                sf::st_intersection(
                        narclim_cells,
                        geography_keep
                )
        )
        
        if (nrow(intersections) == 0L) {
                stop(
                        "No NARCliM cells intersect the selected geography."
                )
        }
        
        # Use an equal-area CRS for areas.
        intersections_area <- sf::st_transform(
                intersections,
                3577
        )
        
        cell_area_source <- sf::st_transform(
                narclim_cells[
                        match(
                                intersections$cell_id,
                                narclim_cells$cell_id
                        ),
                ],
                3577
        )
        
        intersections$intersection_area_m2 <-
                as.numeric(
                        sf::st_area(
                                intersections_area
                        )
                )
        
        intersections$cell_area_m2 <-
                as.numeric(
                        sf::st_area(
                                cell_area_source
                        )
                )
        
        intersections$cell_fraction <-
                intersections$intersection_area_m2 /
                intersections$cell_area_m2
        
        if (
                any(
                        !is.finite(intersections$cell_fraction)
                ) ||
                any(
                        intersections$cell_fraction <= 0
                )
        ) {
                stop(
                        "Invalid native-cell intersection fractions."
                )
        }
        
        intersections
}

###############################################################################
# Build authoritative native NARCliM cell table from remote grid metadata
#
# Uses:
#   rlon / rlat logical coordinates
#   2-D lon / lat WGS84 cell centres
#
# Climate data are NOT read here.
###############################################################################

build_narclim_native_cell_table <- function(
                remote_info,
                earth_radius_m = 6370000
) {
        
        rlon <- as.numeric(remote_info$rlon)
        rlat <- as.numeric(remote_info$rlat)
        
        lon <- remote_info$lon2d
        lat <- remote_info$lat2d
        
        n_rlon <- length(rlon)
        n_rlat <- length(rlat)
        
        if (!identical(dim(lon), dim(lat))) {
                stop(
                        "NARCliM lon and lat arrays have different dimensions."
                )
        }
        
        if (identical(dim(lon), c(n_rlon, n_rlat))) {
                
                lon_matrix <- lon
                lat_matrix <- lat
                
        } else if (identical(dim(lon), c(n_rlat, n_rlon))) {
                
                lon_matrix <- t(lon)
                lat_matrix <- t(lat)
                
        } else {
                
                stop(
                        paste0(
                                "Unexpected NARCliM lon/lat dimensions: ",
                                paste(dim(lon), collapse = " x ")
                        )
                )
        }
        
        native_cells <- tidyr::expand_grid(
                rlat_index = seq_len(n_rlat),
                rlon_index = seq_len(n_rlon)
        ) |>
                dplyr::arrange(
                        rlat_index,
                        rlon_index
                ) |>
                dplyr::mutate(
                        native_cell_id = dplyr::row_number(),
                        
                        rlon = rep(
                                rlon,
                                times = n_rlat
                        ),
                        
                        rlat = rep(
                                rlat,
                                each = n_rlon
                        ),
                        
                        lon = as.vector(lon_matrix),
                        lat = as.vector(lat_matrix)
                ) |>
                dplyr::select(
                        native_cell_id,
                        rlon_index,
                        rlat_index,
                        rlon,
                        rlat,
                        lon,
                        lat
                )
        
        if (
                any(!is.finite(native_cells$lon)) ||
                any(!is.finite(native_cells$lat))
        ) {
                stop(
                        "Non-finite NARCliM lon/lat values detected."
                )
        }
        
        # ----------------------------------------------------------
        # Nominal native cell area
        #
        # Retained from the validated Heat-Mortality implementation.
        # ----------------------------------------------------------
        
        rlon_steps <- abs(
                diff(rlon)
        )
        
        ordinary_rlon_steps <- rlon_steps[
                rlon_steps <
                        stats::quantile(
                                rlon_steps,
                                0.95,
                                na.rm = TRUE,
                                names = FALSE
                        )
        ]
        
        if (length(ordinary_rlon_steps) == 0L) {
                ordinary_rlon_steps <- rlon_steps
        }
        
        delta_rlon_deg <- stats::median(
                ordinary_rlon_steps,
                na.rm = TRUE
        )
        
        delta_rlat_deg <- stats::median(
                abs(diff(rlat)),
                na.rm = TRUE
        )
        
        lower_lat_rad <- pmax(
                -90,
                native_cells$rlat -
                        delta_rlat_deg / 2
        ) * pi / 180
        
        upper_lat_rad <- pmin(
                90,
                native_cells$rlat +
                        delta_rlat_deg / 2
        ) * pi / 180
        
        delta_lon_rad <-
                delta_rlon_deg *
                pi / 180
        
        native_cells <- native_cells |>
                dplyr::mutate(
                        cell_area_km2 =
                                (
                                        earth_radius_m^2 *
                                                delta_lon_rad *
                                                (
                                                        sin(upper_lat_rad) -
                                                                sin(lower_lat_rad)
                                                )
                                ) / 1e6
                )
        
        native_cells
}

###############################################################################
# Pad coordinate matrix for native-cell corner construction
###############################################################################

pad_coordinate_matrix <- function(m) {
        
        nr <- nrow(m)
        nc <- ncol(m)
        
        out <- matrix(
                NA_real_,
                nrow = nr + 2L,
                ncol = nc + 2L
        )
        
        out[
                2:(nr + 1L),
                2:(nc + 1L)
        ] <- m
        
        out[
                1,
                2:(nc + 1L)
        ] <-
                2 * m[1, ] -
                m[2, ]
        
        out[
                nr + 2L,
                2:(nc + 1L)
        ] <-
                2 * m[nr, ] -
                m[nr - 1L, ]
        
        out[
                2:(nr + 1L),
                1
        ] <-
                2 * m[, 1] -
                m[, 2]
        
        out[
                2:(nr + 1L),
                nc + 2L
        ] <-
                2 * m[, nc] -
                m[, nc - 1L]
        
        out[1, 1] <-
                2 * out[2, 2] -
                out[3, 3]
        
        out[1, nc + 2L] <-
                2 * out[2, nc + 1L] -
                out[3, nc]
        
        out[nr + 2L, 1] <-
                2 * out[nr + 1L, 2] -
                out[nr, 3]
        
        out[nr + 2L, nc + 2L] <-
                2 * out[nr + 1L, nc + 1L] -
                out[nr, nc]
        
        out
}

###############################################################################
# Derive native NARCliM cell-corner arrays
###############################################################################

build_native_cell_corner_arrays <- function(
                native_cells
) {
        
        n_rlon <- max(
                native_cells$rlon_index
        )
        
        n_rlat <- max(
                native_cells$rlat_index
        )
        
        lon_matrix <- matrix(
                native_cells$lon,
                nrow = n_rlon,
                ncol = n_rlat
        )
        
        lat_matrix <- matrix(
                native_cells$lat,
                nrow = n_rlon,
                ncol = n_rlat
        )
        
        lon_pad <- pad_coordinate_matrix(
                lon_matrix
        )
        
        lat_pad <- pad_coordinate_matrix(
                lat_matrix
        )
        
        lon_corner <- matrix(
                NA_real_,
                nrow = n_rlon + 1L,
                ncol = n_rlat + 1L
        )
        
        lat_corner <- matrix(
                NA_real_,
                nrow = n_rlon + 1L,
                ncol = n_rlat + 1L
        )
        
        for (i in seq_len(n_rlon + 1L)) {
                
                for (j in seq_len(n_rlat + 1L)) {
                        
                        lon_corner[i, j] <- mean(
                                lon_pad[
                                        i:(i + 1L),
                                        j:(j + 1L)
                                ],
                                na.rm = TRUE
                        )
                        
                        lat_corner[i, j] <- mean(
                                lat_pad[
                                        i:(i + 1L),
                                        j:(j + 1L)
                                ],
                                na.rm = TRUE
                        )
                }
        }
        
        list(
                lon_corner = lon_corner,
                lat_corner = lat_corner
        )
}

###############################################################################
# Build exact native NARCliM polygons potentially overlapping one geography
###############################################################################

build_candidate_native_polygons_generic <- function(
                native_cells,
                corner_arrays,
                geography_one,
                domain = "AUS-18",
                buffer_degrees = NULL
) {
        
        geography_wgs84 <- sf::st_transform(
                geography_one,
                4326
        )
        
        bbox <- sf::st_bbox(
                geography_wgs84
        )
        
        # Conservative candidate-cell buffer only.
        # Exact st_intersection() determines the final overlap.
        
        if (is.null(buffer_degrees)) {
                
                buffer_degrees <- if (
                        grepl(
                                "AUS.?18|20km",
                                domain,
                                ignore.case = TRUE
                        )
                ) {
                        1.0
                } else {
                        0.25
                }
        }
        
        candidate_cells <- native_cells |>
                dplyr::filter(
                        lon >= bbox[["xmin"]] -
                                buffer_degrees,
                        lon <= bbox[["xmax"]] +
                                buffer_degrees,
                        lat >= bbox[["ymin"]] -
                                buffer_degrees,
                        lat <= bbox[["ymax"]] +
                                buffer_degrees
                )
        
        if (nrow(candidate_cells) == 0L) {
                return(NULL)
        }
        
        lon_corner <-
                corner_arrays$lon_corner
        
        lat_corner <-
                corner_arrays$lat_corner
        
        polygon_list <- vector(
                "list",
                nrow(candidate_cells)
        )
        
        for (
                k in seq_len(
                        nrow(candidate_cells)
                )
        ) {
                
                i <-
                        candidate_cells$rlon_index[k]
                
                j <-
                        candidate_cells$rlat_index[k]
                
                coords <- matrix(
                        c(
                                lon_corner[i, j],
                                lat_corner[i, j],
                                
                                lon_corner[i + 1L, j],
                                lat_corner[i + 1L, j],
                                
                                lon_corner[i + 1L, j + 1L],
                                lat_corner[i + 1L, j + 1L],
                                
                                lon_corner[i, j + 1L],
                                lat_corner[i, j + 1L],
                                
                                lon_corner[i, j],
                                lat_corner[i, j]
                        ),
                        ncol = 2,
                        byrow = TRUE
                )
                
                polygon_list[[k]] <-
                        sf::st_polygon(
                                list(coords)
                        )
        }
        
        polygons <- sf::st_sf(
                candidate_cells,
                geometry = sf::st_sfc(
                        polygon_list,
                        crs = 4326
                )
        )
        
        invalid <- !sf::st_is_valid(
                polygons
        )
        
        if (any(invalid)) {
                
                polygons[invalid, ] <-
                        sf::st_make_valid(
                                polygons[invalid, ]
                        )
        }
        
        polygons
}

###############################################################################
# Exact native NARCliM cell x arbitrary geography intersections
#
# Required geography fields:
#   GEO_ID
#   GEO_NAME
###############################################################################

create_cell_geography_intersections <- function(
                native_cells,
                geography_sf,
                domain = "AUS-18"
) {
        
        required_fields <- c(
                "GEO_ID",
                "GEO_NAME"
        )
        
        missing_fields <- setdiff(
                required_fields,
                names(geography_sf)
        )
        
        if (length(missing_fields) > 0L) {
                
                stop(
                        "Geography is missing required field(s): ",
                        paste(
                                missing_fields,
                                collapse = ", "
                        )
                )
        }
        
        corner_arrays <-
                build_native_cell_corner_arrays(
                        native_cells
                )
        
        geography_wgs84 <-
                sf::st_transform(
                        geography_sf,
                        4326
                )
        
        intersection_list <- vector(
                "list",
                nrow(geography_wgs84)
        )
        
        for (
                g in seq_len(
                        nrow(geography_wgs84)
                )
        ) {
                
                geography_one <-
                        geography_wgs84[g, ]
                
                candidate_polygons <-
                        build_candidate_native_polygons_generic(
                                native_cells =
                                        native_cells,
                                
                                corner_arrays =
                                        corner_arrays,
                                
                                geography_one =
                                        geography_one,
                                
                                domain =
                                        domain
                        )
                
                if (is.null(candidate_polygons)) {
                        next
                }
                
                piece <- suppressWarnings(
                        sf::st_intersection(
                                candidate_polygons,
                                
                                geography_one[
                                        ,
                                        c(
                                                "GEO_ID",
                                                "GEO_NAME"
                                        )
                                ]
                        )
                )
                
                if (nrow(piece) == 0L) {
                        next
                }
                
                piece <-
                        piece[
                                !sf::st_is_empty(piece),
                        ]
                
                if (nrow(piece) == 0L) {
                        next
                }
                
                intersection_list[[g]] <-
                        piece[
                                ,
                                c(
                                        "native_cell_id",
                                        "rlon_index",
                                        "rlat_index",
                                        "rlon",
                                        "rlat",
                                        "lon",
                                        "lat",
                                        "cell_area_km2",
                                        "GEO_ID",
                                        "GEO_NAME",
                                        "geometry"
                                )
                        ]
        }
        
        intersections <- dplyr::bind_rows(
                intersection_list
        )
        
        if (nrow(intersections) == 0L) {
                
                stop(
                        "No native NARCliM cell ? geography intersections were produced."
                )
        }
        
        intersections$intersection_area_km2 <-
                as.numeric(
                        sf::st_area(
                                sf::st_transform(
                                        intersections,
                                        3577
                                )
                        )
                ) / 1e6
        
        intersections
}

###############################################################################
# Stream NARCliM daily data into a Climate Menu SpatRaster
#
# The raster remains in the native NARCliM rotated-pole grid.
# CRS parameters are taken from the NARCliM CF grid mapping.
#
# Scientific rules:
#   - preserve daily temporal resolution;
#   - preserve actual CF-derived dates;
#   - preserve missing values;
#   - do not aggregate during acquisition;
#   - do not label rotated coordinates as ordinary WGS84 longitude/latitude.
###############################################################################

stream_narclim_daily_raster <- function(
        opendap_url,
        variable,
        indices,
        t_start = 1L,
        t_count = -1L
) {

    # Use the established NARCliM streaming interface unchanged.
    chunk <- stream_narclim_subset(
        opendap_url = opendap_url,
        variable = variable,
        indices = indices,
        t_start = t_start,
        t_count = t_count
    )

    vals <- chunk$values

    if (length(dim(vals)) != 3L) {
        stop(
            "Remote NARCliM values must have dimensions ",
            "rlon x rlat x time."
        )
    }

    expected_dim <- c(
        as.integer(indices$x_count),
        as.integer(indices$y_count),
        nrow(chunk$time_table)
    )

    if (!identical(
        as.integer(dim(vals)),
        expected_dim
    )) {
        stop(
            "Unexpected remote NARCliM array dimensions. Observed: ",
            paste(dim(vals), collapse = " x "),
            "; expected: ",
            paste(expected_dim, collapse = " x "),
            "."
        )
    }

    if (
        length(chunk$rlon) < 1L ||
        length(chunk$rlat) < 1L
    ) {
        stop("Remote NARCliM spatial coordinates are empty.")
    }

    # rlon/rlat are native NARCliM rotated-grid coordinates.
    # They are NOT ordinary geographic longitude/latitude.
    #
    # The raster therefore remains in native-grid coordinate space
    # with CRS deliberately unset. Geographic operations must use
    # the authoritative NARCliM lon2d/lat2d native-cell pathway.

    xres <- if (length(chunk$rlon) > 1L) {
        stats::median(abs(diff(chunk$rlon)))
    } else {
        1
    }

    yres <- if (length(chunk$rlat) > 1L) {
        stats::median(abs(diff(chunk$rlat)))
    } else {
        1
    }

    if (
        !is.finite(xres) ||
        !is.finite(yres) ||
        xres <= 0 ||
        yres <= 0
    ) {
        stop(
            "Could not determine valid NARCliM native-grid resolution."
        )
    }

    nx <- dim(vals)[1L]
    ny <- dim(vals)[2L]
    nt <- dim(vals)[3L]

    r <- terra::rast(
        nrows = ny,
        ncols = nx,
        nlyrs = nt,
        xmin = min(chunk$rlon) - xres / 2,
        xmax = max(chunk$rlon) + xres / 2,
        ymin = min(chunk$rlat) - yres / 2,
        ymax = max(chunk$rlat) + yres / 2,
        crs = ""
    )

    # ncdf4 returns [x, y, time].
    # terra raster matrices are [row(y), column(x)] and row 1 is ymax.
    y_decreasing <- length(chunk$rlat) > 1L &&
        chunk$rlat[2L] < chunk$rlat[1L]

    for (k in seq_len(nt)) {

        layer_matrix <- t(
            vals[, , k, drop = TRUE]
        )

        if (!y_decreasing) {
            layer_matrix <- layer_matrix[
                nrow(layer_matrix):1L,
                ,
                drop = FALSE
            ]
        }

        terra::values(r[[k]]) <-
            as.vector(t(layer_matrix))
    }

    # Preserve the actual parsed CF dates from the established streamer.
    calendar <- tolower(
        trimws(
            as.character(chunk$calendar)[1L]
        )
    )

    if (!calendar %in% c(
        "standard",
        "gregorian",
        "proleptic_gregorian"
    )) {
        stop(
            "Remote daily-raster adapter currently requires a ",
            "standard/Gregorian calendar. Observed calendar: ",
            calendar,
            "."
        )
    }

    dates <- as.Date(chunk$time_table$cf_date)

    if (
        length(dates) != terra::nlyr(r) ||
        anyNA(dates)
    ) {
        stop(
            "Could not assign a complete CF-derived date vector ",
            "to the remote NARCliM raster."
        )
    }

    terra::time(r) <- dates

    names(r) <- paste0(
        variable,
        "_",
        format(dates, "%Y%m%d")
    )

    attr(r, "narclim_native_grid") <- TRUE
    attr(r, "narclim_spatial_reference") <-
        "rlon_rlat_native_grid_crs_unset"
    attr(r, "narclim_geographic_coordinates") <-
        "use_lon2d_lat2d_for_geographic_operations"
    attr(r, "narclim_calendar") <- chunk$calendar
    attr(r, "narclim_variable") <- variable
    attr(r, "narclim_variable_units") <- chunk$variable_units
    attr(r, "narclim_opendap_url") <- opendap_url

    r
}

# =============================================================================
# CLIMATE MENU V04 - STANDARD STUDY-AREA / ROI LAYER
# =============================================================================
#
# Purpose:
#
# Convert user geography choices into one stable internal object BEFORE
# climate-data acquisition.
#
# The ROI object deliberately separates:
#
#   analysis geography
#   retrieval extent
#   human-readable labels
#   stable internal identifiers
#
# Retrieval padding is NOT part of the analysis geography.
#
# =============================================================================


validate_climate_menu_bbox <- function(bbox) {

  required <- c(
    "lon_min",
    "lon_max",
    "lat_min",
    "lat_max"
  )

  if (!is.list(bbox)) {
    stop("roi_bbox must be a list.")
  }

  if (!all(required %in% names(bbox))) {
    stop(
      "roi_bbox must contain: ",
      paste(required, collapse = ", ")
    )
  }

  values <- unlist(
    bbox[required],
    use.names = TRUE
  )

  if (
    length(values) != 4L ||
    anyNA(values) ||
    !all(is.finite(values))
  ) {
    stop(
      "All four roi_bbox coordinates must be finite numeric values."
    )
  }

  values <- as.numeric(values)
  names(values) <- required

  if (values["lon_min"] >= values["lon_max"]) {
    stop("roi_bbox lon_min must be less than lon_max.")
  }

  if (values["lat_min"] >= values["lat_max"]) {
    stop("roi_bbox lat_min must be less than lat_max.")
  }

  if (
    values["lon_min"] < -180 ||
    values["lon_max"] > 180 ||
    values["lat_min"] < -90 ||
    values["lat_max"] > 90
  ) {
    stop("roi_bbox coordinates fall outside valid WGS84 bounds.")
  }

  values
}


build_climate_menu_roi <- function(
    roi_id = "ROI_001",
    roi_mode = c(
      "named_geography",
      "bbox"
    ),
    geography_type = NULL,
    geography_names = NULL,
    geography_ids = NULL,
    combine_geographies = TRUE,
    boundary_path = NULL,
    boundary_layer = NULL,
    boundary_id_column = NULL,
    boundary_name_column = NULL,
    boundary_filter_column = NULL,
    roi_bbox = NULL,
    padding_cells = 1L
) {

  roi_mode <- match.arg(roi_mode)

  if (
    length(roi_id) != 1L ||
    is.na(roi_id) ||
    !nzchar(roi_id)
  ) {
    stop("roi_id must be one non-empty value.")
  }

  if (
    length(combine_geographies) != 1L ||
    is.na(combine_geographies) ||
    !is.logical(combine_geographies)
  ) {
    stop("combine_geographies must be TRUE or FALSE.")
  }

  if (
    length(padding_cells) != 1L ||
    is.na(padding_cells) ||
    !is.numeric(padding_cells) ||
    padding_cells < 0 ||
    padding_cells != as.integer(padding_cells)
  ) {
    stop("padding_cells must be a non-negative integer.")
  }

  padding_cells <- as.integer(padding_cells)


  # ---------------------------------------------------------------------------
  # EXPLICIT WGS84 BOUNDING BOX
  # ---------------------------------------------------------------------------

  if (identical(roi_mode, "bbox")) {

    bbox_values <- validate_climate_menu_bbox(
      roi_bbox
    )

    return(
      structure(
        list(
          roi_id = roi_id,
          roi_mode = "bbox",
          geography_type = "custom_bbox",

          geography_ids = NULL,
          geography_names = NULL,

          combine_geographies = TRUE,

          geometry = NULL,

          bbox_wgs84 = bbox_values,

          retrieval_padding_cells =
            padding_cells,

          analysis_definition =
            "explicit_wgs84_bbox"
        ),
        class = c(
          "climate_menu_roi",
          "list"
        )
      )
    )
  }


  # ---------------------------------------------------------------------------
  # NAMED POLYGON GEOGRAPHY
  # ---------------------------------------------------------------------------

  if (
    is.null(boundary_path) ||
    length(boundary_path) != 1L ||
    !nzchar(boundary_path)
  ) {
    stop(
      "boundary_path is required for named geography."
    )
  }

  if (!file.exists(boundary_path)) {
    stop(
      "Boundary file does not exist: ",
      boundary_path
    )
  }

  if (
    is.null(boundary_name_column) ||
    length(boundary_name_column) != 1L ||
    !nzchar(boundary_name_column)
  ) {
    stop(
      "boundary_name_column is required for named geography."
    )
  }

  geography <-
    if (
      is.null(boundary_layer) ||
      length(boundary_layer) == 0L ||
      is.na(boundary_layer) ||
      !nzchar(boundary_layer)
    ) {

      sf::st_read(
        dsn = boundary_path,
        quiet = TRUE
      )

    } else {

      sf::st_read(
        dsn = boundary_path,
        layer = boundary_layer,
        quiet = TRUE
      )
    }

  if (nrow(geography) == 0L) {
    stop("Boundary dataset contains no features.")
  }

  if (
    is.na(sf::st_crs(geography))
  ) {
    stop(
      "Boundary dataset has no defined CRS."
    )
  }

  if (
    !boundary_name_column %in%
      names(geography)
  ) {
    stop(
      "Boundary name column not found: ",
      boundary_name_column
    )
  }

  if (
    !is.null(boundary_id_column) &&
    !boundary_id_column %in%
      names(geography)
  ) {
    stop(
      "Boundary ID column not found: ",
      boundary_id_column
    )
  }

  # Prefer explicit IDs when supplied.
  # Otherwise select using human-readable names.

  if (
    !is.null(geography_ids) &&
    length(geography_ids) > 0L
  ) {

    if (is.null(boundary_id_column)) {
      stop(
        "geography_ids were supplied but boundary_id_column is NULL."
      )
    }

    keep <- as.character(
      geography[[boundary_id_column]]
    ) %in% as.character(
      geography_ids
    )

  } else {

    if (
      is.null(geography_names) ||
      length(geography_names) == 0L
    ) {
      stop(
        "Supply geography_ids or geography_names."
      )
    }

    filter_column <- boundary_filter_column

    if (
      is.null(filter_column) ||
      !nzchar(filter_column)
    ) {
      filter_column <- boundary_name_column
    }

    if (!filter_column %in% names(geography)) {
      stop(
        "Boundary filter column not found: ",
        filter_column
      )
    }

    keep <- as.character(
      geography[[filter_column]]
    ) %in% as.character(
      geography_names
    )
  }

  selected <- geography[
    keep,
    ,
    drop = FALSE
  ]

  if (nrow(selected) == 0L) {
    stop(
      "No boundary features matched the requested geography."
    )
  }


  # Check that every requested name/ID matched.

  if (
    !is.null(geography_ids) &&
    length(geography_ids) > 0L
  ) {

    found_ids <- unique(
      as.character(
        selected[[boundary_id_column]]
      )
    )

    missing_ids <- setdiff(
      as.character(geography_ids),
      found_ids
    )

    if (length(missing_ids) > 0L) {
      stop(
        "Requested geography IDs not found: ",
        paste(
          missing_ids,
          collapse = ", "
        )
      )
    }

  } else {

    filter_column <- boundary_filter_column

    if (
      is.null(filter_column) ||
      !nzchar(filter_column)
    ) {
      filter_column <- boundary_name_column
    }

    found_names <- unique(
      as.character(
        selected[[filter_column]]
      )
    )

    missing_names <- setdiff(
      as.character(geography_names),
      found_names
    )

    if (length(missing_names) > 0L) {
      stop(
        "Requested geography names not found: ",
        paste(
          missing_names,
          collapse = ", "
        )
      )
    }
  }


  # Preserve the selected polygon identities for later summaries.

  selected_ids <-
    if (!is.null(boundary_id_column)) {
      as.character(
        selected[[boundary_id_column]]
      )
    } else {
      as.character(seq_len(nrow(selected)))
    }

  selected_names <- as.character(
    selected[[boundary_name_column]]
  )


  # Geometry used for the overall analysis/retrieval envelope.
  #
  # If combine_geographies == FALSE, individual polygons remain in
  # `selected_features`; the union is used only to obtain one efficient
  # retrieval envelope.

  selected_wgs84 <- sf::st_transform(
    selected,
    4326
  )

  analysis_geometry <-
    if (combine_geographies) {

      sf::st_as_sf(
        data.frame(
          roi_id = roi_id
        ),
        geometry = sf::st_sfc(
          sf::st_union(
            sf::st_geometry(
              selected_wgs84
            )
          ),
          crs = 4326
        )
      )

    } else {

      selected_wgs84
    }


  envelope_geometry <- sf::st_union(
    sf::st_geometry(
      selected_wgs84
    )
  )

  bbox <- sf::st_bbox(
    envelope_geometry
  )

  bbox_values <- c(
    lon_min = unname(bbox["xmin"]),
    lon_max = unname(bbox["xmax"]),
    lat_min = unname(bbox["ymin"]),
    lat_max = unname(bbox["ymax"])
  )


  structure(
    list(
      roi_id = roi_id,

      roi_mode =
        "named_geography",

      geography_type =
        geography_type,

      geography_ids =
        selected_ids,

      geography_names =
        selected_names,

      combine_geographies =
        combine_geographies,

      geometry =
        analysis_geometry,

      selected_features =
        selected_wgs84,

      bbox_wgs84 =
        bbox_values,

      retrieval_padding_cells =
        padding_cells,

      boundary_source =
        boundary_path,

      boundary_id_column =
        boundary_id_column,

      boundary_name_column =
        boundary_name_column,

      analysis_definition =
        if (combine_geographies) {
          "combined_named_geography"
        } else {
          "separate_named_geographies"
        }
    ),
    class = c(
      "climate_menu_roi",
      "list"
    )
  )
}


print.climate_menu_roi <- function(
    x,
    ...
) {

  cat(
    "Climate Menu ROI\n",
    "  ID: ",
    x$roi_id,
    "\n",
    "  mode: ",
    x$roi_mode,
    "\n",
    "  geography type: ",
    x$geography_type,
    "\n",
    sep = ""
  )

  if (
    length(x$geography_names) > 0L
  ) {
    cat(
      "  areas: ",
      paste(
        x$geography_names,
        collapse = " | "
      ),
      "\n",
      sep = ""
    )
  }

  cat(
    "  analysis: ",
    x$analysis_definition,
    "\n",
    "  retrieval padding: ",
    x$retrieval_padding_cells,
    " cell(s)\n",
    sep = ""
  )

  if (
    !is.null(x$bbox_wgs84)
  ) {
    cat(
      "  WGS84 bbox: ",
      paste(
        round(
          x$bbox_wgs84,
          5
        ),
        collapse = ", "
      ),
      "\n",
      sep = ""
    )
  }

  invisible(x)
}

###############################################################################
# Climate Menu remote daily-period acquisition
# Added 2026-10-26 v06 production integration
###############################################################################

stream_narclim_period_daily_file <- function(
    opendap_url,
    variable,
    indices,
    remote_info,
    start_year,
    end_year,
    output_file,
    chunk_days = 31L
) {

    required <- c(
        'time_values',
        'time_units',
        'calendar'
    )

    missing_fields <- setdiff(
        required,
        names(remote_info)
    )

    if (length(missing_fields) > 0L) {
        stop(
            'Remote resource is missing required CF time metadata: ',
            paste(missing_fields, collapse = ', ')
        )
    }

    calendar <- tolower(
        as.character(remote_info$calendar[[1]])
    )

    if (!calendar %in% c(
        'standard',
        'gregorian',
        'proleptic_gregorian'
    )) {
        stop(
            'Remote calendar is not currently supported by this ',
            'production helper: ',
            remote_info$calendar[[1]],
            '. Refusing to infer dates.'
        )
    }

    units <- as.character(
        remote_info$time_units[[1]]
    )

    origin_match <- regexec(
        '^days since[[:space:]]+([0-9]{4}-[0-9]{2}-[0-9]{2})',
        units,
        ignore.case = TRUE
    )

    origin_parts <- regmatches(
        units,
        origin_match
    )[[1]]

    if (length(origin_parts) < 2L) {
        stop(
            'Unsupported CF time units: ',
            units,
            '. Expected days since YYYY-MM-DD.'
        )
    }

    origin <- as.Date(origin_parts[[2]])

    time_values <- as.numeric(
        remote_info$time_values
    )

    # NARCliM daily timestamps are centred on the day (e.g. .5).
    # floor() maps those timestamps to the corresponding calendar day.
    dates <- origin + floor(time_values)

    if (length(dates) != as.integer(remote_info$n_time)) {
        stop(
            'CF date count does not match remote n_time.'
        )
    }

    years <- as.integer(
        format(dates, '%Y')
    )

    wanted <- which(
        years >= as.integer(start_year) &
        years <= as.integer(end_year)
    )

    if (length(wanted) == 0L) {
        stop(
            'No remote daily timesteps found for requested period ',
            start_year, '-', end_year, '.'
        )
    }

    if (
        length(wanted) > 1L &&
        any(diff(wanted) != 1L)
    ) {
        stop(
            'Requested remote time indices are not contiguous.'
        )
    }

    chunk_days <- as.integer(chunk_days)

    if (
        !is.finite(chunk_days) ||
        chunk_days < 1L
    ) {
        stop(
            'chunk_days must be a positive integer.'
        )
    }

    starts <- seq.int(
        from = wanted[[1]],
        to = wanted[[length(wanted)]],
        by = chunk_days
    )

    chunks <- vector(
        'list',
        length(starts)
    )

    for (i in seq_along(starts)) {

        first <- starts[[i]]

        remaining <-
            wanted[[length(wanted)]] - first + 1L

        count <- min(
            chunk_days,
            remaining
        )

        message(
            'Streaming ',
            variable,
            ': ',
            as.character(dates[[first]]),
            ' to ',
            as.character(dates[[first + count - 1L]]),
            ' (', count, ' daily layers)'
        )

        chunks[[i]] <- stream_narclim_daily_raster(
            opendap_url = opendap_url,
            variable = variable,
            indices = indices,
            t_start = first,
            t_count = count
        )
    }

    r <- do.call(c, chunks)

    if (
        is.null(r) ||
        terra::nlyr(r) != length(wanted)
    ) {
        stop(
            'Remote period stream returned an unexpected layer count.'
        )
    }

    selected_dates <- dates[wanted]

    names(r) <- paste0(
        'day_',
        format(selected_dates, '%Y_%m_%d')
    )

    # Preserve actual calendar dates as SpatRaster time metadata as well
    # as layer names, so downstream calendar-aware processing has both.
    terra::time(r) <- selected_dates

    # Convert temperature only when the source metadata says Kelvin.
    # This avoids relying solely on a numerical magnitude heuristic.
    source_units <-
        if ('variable_units' %in% names(remote_info)) {
            as.character(remote_info$variable_units[[1]])
        } else {
            NA_character_
        }

    if (
        grepl('^tas', variable, ignore.case = TRUE) &&
        !is.na(source_units) &&
        tolower(trimws(source_units)) %in% c(
            'k',
            'kelvin'
        )
    ) {
        r <- r - 273.15
    }

    dir.create(
        dirname(output_file),
        recursive = TRUE,
        showWarnings = FALSE
    )

    safe_write_raster(
        r,
        output_file
    )

    if (!file.exists(output_file)) {
        stop(
            'Remote raster output was not created.'
        )
    }

    output_file
}

