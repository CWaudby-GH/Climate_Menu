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