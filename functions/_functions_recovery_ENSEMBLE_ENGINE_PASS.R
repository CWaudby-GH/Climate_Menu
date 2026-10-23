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