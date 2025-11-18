# _targets_functions.R
# -------------------
# Load required libraries


# -------------------
# Function: aggregate_time_period
aggregate_time_period <- function(data, time_var = "time", period = "month", fun = mean) {
        # Aggregate raster or vector data over a time period
        if (inherits(data, "RasterStack") | inherits(data, "SpatRaster")) {
                # For raster/SpatRaster
                time_index <- lubridate::floor_date(as.Date(names(data)), unit = period)
                aggregated <- terra::aggregate(data, by = time_index, fun = fun)
        } else if (inherits(data, "data.frame") | inherits(data, "data.table")) {
                # For data frames
                data[[time_var]] <- lubridate::floor_date(as.Date(data[[time_var]]), unit = period)
                aggregated <- data %>%
                        group_by(across(-matches(time_var))) %>%
                        summarise(value = fun(value, na.rm = TRUE), .groups = "drop")
        } else {
                stop("Unsupported data type for aggregation.")
        }
        return(aggregated)
}

# -------------------
# Function: convert_precip_units
convert_precip_units <- function(precip_data, from = "kg/m2/s", to = "mm/day") {
        # Convert precipitation units
        if (from == "kg/m2/s" & to == "mm/day") {
                converted <- precip_data * 86400
        } else if (from == "mm/day" & to == "kg/m2/s") {
                converted <- precip_data / 86400
        } else {
                stop("Unsupported unit conversion.")
        }
        return(converted)
}

# -------------------
# Function: polygon_aggregate
polygon_aggregate <- function(raster_data, polygons_sf, fun = mean) {
        # Aggregate raster over polygons
        if (!inherits(polygons_sf, "sf")) stop("polygons_sf must be an sf object")
        result <- exactextractr::exact_extract(raster_data, polygons_sf, fun = fun)
        return(result)
}

# -------------------
# Function: percentile_days
percentile_days <- function(raster_data, threshold = 35) {
        # Count days exceeding threshold
        if (!inherits(raster_data, "RasterStack") & !inherits(raster_data, "SpatRaster")) {
                stop("raster_data must be RasterStack or SpatRaster")
        }
        exceed <- raster_data > threshold
        count <- terra::app(exceed, sum)
        return(count)
}

# -------------------
# Function: heatwave_days
heatwave_days <- function(raster_data, threshold = 35, consecutive = 3) {
        # Identify consecutive days exceeding threshold
        if (!inherits(raster_data, "RasterStack") & !inherits(raster_data, "SpatRaster")) {
                stop("raster_data must be RasterStack or SpatRaster")
        }
        exceed <- raster_data > threshold
        # Apply run-length encoding to count consecutive days
        heatwave_count <- terra::app(exceed, function(x) {
                rle_x <- rle(x)
                max_streak <- max(rle_x$lengths[rle_x$values], na.rm = TRUE)
                return(ifelse(max_streak >= consecutive, max_streak, 0))
        })
        return(heatwave_count)
}

# -------------------
# Function: raster_stack_mean
raster_stack_mean <- function(raster_data) {
        if (!inherits(raster_data, "RasterStack") & !inherits(raster_data, "SpatRaster")) {
                stop("raster_data must be RasterStack or SpatRaster")
        }
        mean_raster <- terra::app(raster_data, mean, na.rm = TRUE)
        return(mean_raster)
}

# -------------------
# Function: time_series_extract
time_series_extract <- function(raster_data, points_sf) {
        if (!inherits(points_sf, "sf")) stop("points_sf must be sf object")
        ts_data <- terra::extract(raster_data, vect(points_sf))
        return(ts_data)
}

# -------------------
# Function: combine_targets
combine_targets <- function(..., fun = sum) {
        rasters <- list(...)
        combined <- do.call(terra::app, c(rasters, list(fun = fun)))
        return(combined)
}

# -------------------
# Function: save_raster_png
save_raster_png <- function(raster_data, file_path, zlim = NULL) {
        if (!inherits(raster_data, "Raster") & !inherits(raster_data, "SpatRaster")) {
                stop("raster_data must be Raster or SpatRaster")
        }
        png(filename = file_path, width = 1200, height = 900)
        plot(raster_data, zlim = zlim)
        dev.off()
}

# -------------------
# Any other functions from your previous pipeline can be added here
# Keep names and functionality exactly as in your working _targets_14g_11.R
# File: functions/_targets_functions.R
# ----------------------
# Helper functions for Targets pipeline
# ----------------------

# ----------------------
# Utilities
# ----------------------
# Helper: convert expand.grid to list of lists for tar_map
map_city_var <- function(cities, vars) {
        expand.grid(city = cities, var = vars, stringsAsFactors = FALSE) %>%
                split(seq(nrow(.))) %>%
                purrr::map(~list(city = .$city, var = .$var))
}


# Helper function: safe aggregation and writing
# safe_aggregate_raster <- function(raster_file, start_year, end_year, period = c("monthly", "annual"), 
#                                   output_dir, prefix = "FFD", overwrite = TRUE) {
#         period <- match.arg(period)
#         
#         # Load raster
#         r <- terra::rast(raster_file)
#         
#         # Skip if too few layers (monthly requires >=12)
#         if (period == "monthly" && nlyr(r) < 12) {
#                 message("Skipping monthly aggregation: fewer than 12 layers.")
#                 return(NULL)
#         }
#         
#         # Aggregate
#         r_agg <- aggregate_time_period(r, start_year, end_year, period)
#         
#         # Check if aggregation produced valid values
#         if (is.null(r_agg) || nlyr(r_agg) == 0 || all(is.na(terra::values(r_agg)))) {
#                 message("Skipping ", period, " aggregation: no valid raster cells.")
#                 return(NULL)
#         }
#         
#         # Build output path
#         path <- file.path(
#                 output_dir,
#                 paste0(prefix, "_", period, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif")
#         )
#         
#         # Write raster safely
#         safe_write_raster(r_agg, path, overwrite = overwrite)
#         return(path)
# }

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
        
        r_agg <- agg_list
        
        outfile <- file.path(
                output_dir,
                paste0(prefix, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif")
        )
        
        safe_write_raster(r_agg, outfile, overwrite = TRUE)
        outfile
}  # <- this closes safe_aggregate_raster


safe_write_raster <- function(r, path, overwrite = FALSE) {
        if (is.null(r) || nlyr(r) == 0 || all(is.na(r[]))) {
                message("Raster is empty, not writing: ", path)
                return(NULL)
        }
        terra::writeRaster(r, filename = path, overwrite = overwrite)
        path
}

aggregated <- aggregate_time_period(r, start_year, end_year, period)
if (is.null(aggregated)) {
        message("Aggregation returned NULL: skipping.")
        return(NULL)
}

all_na <- all(sapply(1:nlyr(aggregated), function(i) all(is.na(aggregated[[i]][]))))
if (all_na) {
        message("All aggregated layers are NA: skipping write.")
        return(NULL)
}

r_agg <- terra::rast(aggregated)  # convert list of rasters to SpatRaster
path <- file.path(
        output_dir, 
        paste0(prefix, "_", period, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".tif")
)
safe_write_raster(r_agg, path, overwrite = TRUE)
return(path) 

safe_name <- function(x) {
        x <- as.character(x)
        x <- gsub("[[:space:]/\\\\]+", "_", x)
        x <- gsub("[^A-Za-z0-9_.:-]", "", x)
        x
}



# safe_scale_fill_viridis <- function(...) {
#         
#         # If viridis is available → use the correct function
#         if ("viridis" %in% loadedNamespaces()) {
#                 return(viridis::scale_fill_viridis_c(...))
#         }
#         
#         # If viridisLite is available → also use *the correct function*
#         if ("viridisLite" %in% loadedNamespaces()) {
#                 return(ggplot2::scale_fill_viridis_c(...))   # wrapper provided by ggplot2
#         }
#         
#         # Fallback (if neither is available)
#         ggplot2::scale_fill_gradient(low = "white", high = "black", ...)
# }
# safe_scale_fill_viridis <- function(...){
#         # prefer viridis package (exports scale_fill_viridis_c)
#         if (requireNamespace("viridis", quietly = TRUE)) {
#                 return(viridis::scale_fill_viridis_c(...))
#         }
#         # otherwise fall back to a simple gradient from ggplot2
#         ggplot2::scale_fill_gradient(low = "white", high = "black", ...)
# }

safe_scale_fill_viridis <- function(...) {
        ggplot2::scale_fill_viridis_c(...)
}


safe_ggsave <- function(df, csv_file, png_file, x_col = "Month", y_col = "Mean") {
        write.csv(df, csv_file, row.names = FALSE)
        p <- ggplot2::ggplot(df, ggplot2::aes_string(x = x_col, y = y_col, fill = y_col)) +
                ggplot2::geom_col() +
                safe_scale_fill_viridis() +
                ggplot2::theme_minimal()
        ggplot2::ggsave(filename = png_file, plot = p)
        return(csv_file)
}

# ----------------------
# File & geometry handling
# ----------------------
filter_nc_files <- function(data_dir, variables, scenarios, models = NULL) {
        all_files <- list.files(data_dir, pattern = "\\.nc$", full.names = TRUE, recursive = TRUE)
        if (length(all_files) == 0) return(character(0))
        pat_var <- paste(tolower(variables), collapse = "|")
        pat_scen <- paste(tolower(scenarios), collapse = "|")
        keep <- all_files %>% purrr::keep(~ stringr::str_detect(tolower(.x), pat_var)) %>%
                purrr::keep(~ stringr::str_detect(tolower(.x), pat_scen))
        if (!is.null(models) && length(models) > 0) {
                pat_model <- paste(tolower(models), collapse = "|")
                keep <- keep %>% purrr::keep(~ stringr::str_detect(tolower(.x), pat_model))
        }
        sort(keep)
}

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

# ----------------------
# Raster loading & manipulation
# ----------------------
convert_precip_to_mmday <- function(r) {
        if (is.null(r)) return(NULL)
        r * 86400
}

load_and_crop_rasters <- function(files, var_token, start_year, end_year,
                                  test_mode = FALSE, city_shape = NULL,
                                  sample_layers = NULL, threshold_var = NULL) {
        if (is.null(files) || length(files) == 0) return(NULL)
        
        # Filter files by variable token
        files_var <- files[stringr::str_detect(tolower(basename(files)), tolower(var_token))]
        if (!is.null(threshold_var)) {
                files_var <- files_var[stringr::str_detect(tolower(basename(files_var)), tolower(threshold_var))]
        }
        if (length(files_var) == 0) return(NULL)
        
        # Load rasters
        rasters <- lapply(files_var, function(f) tryCatch(terra::rast(f), error = function(e) NULL))
        rasters <- purrr::compact(rasters)
        if (length(rasters) == 0) return(NULL)
        
        # Combine into one SpatRaster
        r_all <- do.call(c, rasters)
        
        # Time subset
        time_vals <- tryCatch(terra::time(r_all), error = function(e) NULL)
        if (!is.null(time_vals)) {
                yrs <- lubridate::year(time_vals)
                keep_idx <- which(yrs >= start_year & yrs <= end_year)
                if (length(keep_idx) == 0) return(NULL)
                r_all <- r_all[[keep_idx]]
        }
        
        # Sample layers if test_mode is TRUE
        if (!is.null(sample_layers) && sample_layers < nlyr(r_all)) {
                set.seed(42)
                idx <- sample(seq_len(nlyr(r_all)), sample_layers)
                r_all <- r_all[[idx]]
        }
        
        # Crop/mask only if city_shape is provided
        if (!is.null(city_shape)) {
                try({
                        if (!is.na(sf::st_crs(city_shape)$epsg) && !is.na(terra::crs(r_all))) {
                                if (sf::st_crs(city_shape) != sf::st_crs(r_all)) {
                                        city_shape <- sf::st_transform(city_shape, crs = terra::crs(r_all))
                                }
                        }
                }, silent = TRUE)
                try({
                        r_all <- terra::crop(r_all, terra::vect(city_shape))
                        r_all <- terra::mask(r_all, terra::vect(city_shape))
                }, silent = TRUE)
        }
        
        r_all
}

align_and_resample <- function(ref, target, method = "bilinear", tol = 1e-6) {
        if (is.null(ref) || is.null(target)) return(NULL)
        if (terra::crs(ref) != terra::crs(target)) {
                target <- terra::project(target, terra::crs(ref))
        }
        ext_diff <- max(abs(as.numeric(terra::ext(ref)) - as.numeric(terra::ext(target))))
        res_diff <- max(abs(as.numeric(terra::res(ref)) - as.numeric(terra::res(target))))
        if (ext_diff > tol || res_diff > tol) {
                target <- terra::resample(target, ref, method = method)
        }
        target
}

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

# safe_write_raster <- function(r, filename) {
#         if (is.null(r)) return(NULL)
#         filename <- as.character(filename)
#         
#         # Ensure directory exists
#         dir.create(dirname(filename), showWarnings = FALSE, recursive = TRUE)
#         
#         # Ensure SpatRaster
#         if (!inherits(r, "SpatRaster")) {
#                 r <- terra::rast(r)
#         }
#         
#         terra::writeRaster(
#                 r,
#                 filename,
#                 overwrite = TRUE,
#                 wopt = list(datatype = "FLT4S", gdal = c("COMPRESS=LZW"))
#         )
#         
#         filename
# }

# safe_write_raster <- function(r, path, overwrite = TRUE) {
#         # r: SpatRaster
#         # path: output path
#         # overwrite: whether to overwrite existing file
#         terra::writeRaster(r, path, overwrite = overwrite)
#         return(path)
# }
# ----------------------
# Aggregation & statistics
# ----------------------
percentile_map <- function(r, prob) {
        if (is.null(r)) return(NULL)
        tryCatch(terra::app(r, fun = function(x) stats::quantile(x, probs = prob, na.rm = TRUE)), error = function(e) NULL)
}

deep_mean <- function(r) {
        if (is.null(r)) return(NULL)
        if (nlyr(r) == 0) return(NULL)
        terra::app(r, fun = mean, na.rm = TRUE)
}

# aggregate_time_period <- function(r, start_year, end_year, period = c("annual", "monthly"), var_is_precip = FALSE) {
#         if (is.null(r) || nlyr(r) == 0) return(NULL)
#         period <- match.arg(period)
#         
#         # Try to extract time metadata
#         time_vals <- tryCatch(terra::time(r), error = function(e) NULL)
#         
#         if (!is.null(time_vals)) {
#                 yrs <- lubridate::year(time_vals)
#                 mons <- lubridate::month(time_vals)
#         } else {
#                 # Fallback: assume sequential layers represent time
#                 yrs <- rep(seq(start_year, end_year, length.out = nlyr(r)), each = 1)
#                 mons <- rep(1:12, length.out = nlyr(r))
#         }
#         
#         # Filter layers by year range
#         keep_idx <- which(yrs >= start_year & yrs <= end_year)
#         if (length(keep_idx) == 0) return(NULL)
#         r <- r[[keep_idx]]
#         
#         # Group layers by year or month
#         groups <- if (period == "annual") yrs[keep_idx] else mons[keep_idx]
#         grouped <- split(seq_along(keep_idx), groups)
#         
#         # Aggregate each group
#         out_list <- lapply(grouped, function(idx) {
#                 if (var_is_precip && period == "annual") {
#                         terra::app(r[[idx]], fun = sum, na.rm = TRUE)
#                 } else {
#                         terra::app(r[[idx]], fun = mean, na.rm = TRUE)
#                 }
#         })
#         
#         # Combine into one SpatRaster
#         out <- do.call(c, out_list)
#         if (!inherits(out, "SpatRaster")) out <- terra::rast(out_list)
#         
#         names(out) <- paste0(period, "_", names(grouped))
#         terra::ext(out) <- terra::ext(r)
#         terra::crs(out) <- terra::crs(r)
#         
#         return(out)
# }

aggregate_time_period <- function(r, start_year, end_year, period = c("annual", "monthly"), var_is_precip = FALSE) {
        if (is.null(r) || nlyr(r) == 0) return(NULL)
        period <- match.arg(period)
        
        n_layers <- nlyr(r)
        
        # Grouping logic
        if (period == "monthly") {
                # If fewer than 12 layers, each layer becomes its own month
                groups <- if (n_layers < 12) seq_len(n_layers) else rep(1:12, length.out = n_layers)
        } else if (period == "annual") {
                # Group layers into blocks of 12, or treat each layer as a year if fewer than 12
                n_years <- if (n_layers < 12) n_layers else ceiling(n_layers / 12)
                groups <- if (n_layers < 12) seq_len(n_layers) else rep(1:n_years, each = 12)[1:n_layers]
        }
        
        grouped <- split(seq_len(n_layers), groups)
        
        # Aggregate each group
        out_list <- lapply(grouped, function(idx) {
                terra::app(r[[idx]], fun = mean, na.rm = TRUE)
        })
        
        # Combine into one SpatRaster
        out <- do.call(c, out_list)
        if (!inherits(out, "SpatRaster")) out <- terra::rast(out_list)
        
        names(out) <- paste0(period, "_", seq_along(grouped))
        terra::ext(out) <- terra::ext(r)
        terra::crs(out) <- terra::crs(r)
        
        return(out)
}

safe_aggregate_time_period <- function(r, start_year, end_year, period) {
        if (is.null(r) || nlyr(r) == 0) {
                message("Input raster is NULL or has no layers.")
                return(terra::rast())  # return empty SpatRaster
        }
        
        # Call your original aggregation function
        aggregated <- tryCatch({
                aggregate_time_period(r, start_year, end_year, period)
        }, error = function(e) {
                message("Aggregation error: ", e$message)
                return(NULL)
        })
        
        # If aggregation returned NULL, replace with empty raster
        if (is.null(aggregated)) {
                message("Aggregation returned NULL, creating empty raster with same dimensions.")
                aggregated <- terra::rast(nrows = nrow(r), ncols = ncol(r), nlyr = nlyr(r), ext = ext(r), crs = crs(r))
                aggregated[] <- NA
        }
        
        # Remove layers that are completely NA
        keep_layers <- sapply(1:nlyr(aggregated), function(i) !all(is.na(aggregated[[i]][])))
        if (!any(keep_layers)) {
                message("All aggregated layers are NA, returning empty raster with one NA layer.")
                aggregated <- aggregated[[1]]
                aggregated[] <- NA
        } else if (!all(keep_layers)) {
                aggregated <- aggregated[[keep_layers]]
        }
        
        terra::rast(aggregated)
}

# ----------------------
# Heatwave & threshold
# ----------------------
compute_hw_metrics <- function(temp_raster, threshold = 308.15, min_duration = 3) {
        if (is.null(temp_raster)) return(NULL)
        hot_days <- temp_raster > (threshold - 273.15)
        HWF <- terra::app(hot_days, fun = function(x) {
                if (all(is.na(x))) return(NA_real_)
                r <- rle(as.logical(x))
                sum(r$lengths[r$values == TRUE & r$lengths >= min_duration])
        })
        HWD <- terra::app(hot_days, fun = function(x) {
                if (all(is.na(x))) return(NA_real_)
                r <- rle(as.logical(x))
                max(c(0, r$lengths[r$values == TRUE & r$lengths >= min_duration]))
        })
        HWN <- terra::app(hot_days, fun = function(x) {
                if (all(is.na(x))) return(NA_real_)
                r <- rle(as.logical(x))
                length(r$lengths[r$values == TRUE & r$lengths >= min_duration])
        })
        list(HWF = HWF, HWD = HWD, HWN = HWN)
}

compute_threshold_days <- function(hist_raster = NULL, proj_raster = NULL, threshold) {
        out <- list()
        if (!is.null(hist_raster)) out$hist <- terra::app(hist_raster, fun = function(x) sum(x >= threshold, na.rm = TRUE))
        if (!is.null(proj_raster)) out$proj <- terra::app(proj_raster, fun = function(x) sum(x >= threshold, na.rm = TRUE))
        out
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

# ----------------------
# Convert SpatRaster -> long df
# ----------------------
to_df <- function(r) {
        if (is.null(r)) return(data.frame())
        df <- as.data.frame(r, xy = TRUE)
        if (ncol(df) <= 2) return(df)
        tidyr::pivot_longer(df, cols = -c(x, y), names_to = "layer", values_to = "value")
}

# ----------------------
# Plotting functions
# ----------------------
# ----------------------
# Plotting functions
# ----------------------

plot_raster_comparison <- function(hist_df, proj_df = NULL, delta_df = NULL, plot_extent, variable, hist_range, proj_range, output_dir) {
        if (is.null(plot_extent) || length(plot_extent) < 4) plot_extent <- c(-180,180,-90,90)
        vals <- c(if (nrow(hist_df)>0) hist_df$value else NA,
                  if (!is.null(proj_df) && nrow(proj_df)>0) proj_df$value else NA,
                  if (!is.null(delta_df) && nrow(delta_df)>0) delta_df$value else NA)
        vals <- vals[!is.na(vals)]
        if (length(vals) == 0) return(NULL)
        minv <- min(vals); maxv <- max(vals); mid <- (minv+maxv)/2
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        
        # Historical
        if (nrow(hist_df) > 0) {
                p_hist <- ggplot2::ggplot(hist_df, ggplot2::aes(x=x,y=y,fill=value)) +
                        ggplot2::geom_raster() +
                        ggplot2::coord_quickmap(xlim = plot_extent[1:2], ylim = plot_extent[3:4]) +
                        ggplot2::scale_fill_gradient2(limits=c(minv,maxv), low="blue", mid="white", high="red", midpoint=mid, na.value="palegoldenrod") +
                        ggplot2::labs(title=paste(variable,"Historical",hist_range)) +
                        ggplot2::theme_minimal()
                f <- file.path(output_dir, paste0(safe_name(variable),"_historical_",timestamp,".png"))
                ggplot2::ggsave(f, p_hist, width=8, height=6)
        }
        
        # Projected
        if (!is.null(proj_df) && nrow(proj_df) > 0) {
                p_proj <- ggplot2::ggplot(proj_df, ggplot2::aes(x=x,y=y,fill=value)) +
                        ggplot2::geom_raster() +
                        ggplot2::coord_quickmap(xlim = plot_extent[1:2], ylim = plot_extent[3:4]) +
                        ggplot2::scale_fill_gradient2(limits=c(minv,maxv), low="blue", mid="white", high="red", midpoint=mid, na.value="palegoldenrod") +
                        ggplot2::labs(title=paste(variable,"Projected",proj_range)) +
                        ggplot2::theme_minimal()
                f2 <- file.path(output_dir, paste0(safe_name(variable),"_projected_",timestamp,".png"))
                ggplot2::ggsave(f2, p_proj, width=8, height=6)
        }
        
        # Delta
        if (!is.null(delta_df) && nrow(delta_df) > 0) {
                dl <- max(abs(delta_df$value), na.rm = TRUE)
                p_d <- ggplot2::ggplot(delta_df, ggplot2::aes(x=x,y=y,fill=value)) +
                        ggplot2::geom_raster() +
                        ggplot2::coord_quickmap(xlim=plot_extent[1:2], ylim=plot_extent[3:4]) +
                        ggplot2::scale_fill_gradient2(limits=c(-dl,dl), low="blue", mid="white", high="red", midpoint=0, na.value="palegoldenrod") +
                        ggplot2::labs(title=paste(variable,"Delta (proj-hist)")) +
                        ggplot2::theme_minimal()
                f3 <- file.path(output_dir, paste0(safe_name(variable),"_delta_",timestamp,".png"))
                ggplot2::ggsave(f3, p_d, width=8, height=6)
        }
        
        TRUE
}

plot_monthly_climatology <- function(monthly_df, plot_extent, variable, range_label, output_dir) {
        if (nrow(monthly_df) == 0) return(NULL)
        monthly_df$month <- factor(monthly_df$layer, levels = paste0("monthly_",1:12), labels = month.abb)
        p <- ggplot2::ggplot(monthly_df, ggplot2::aes(x=x,y=y,fill=value)) +
                ggplot2::geom_raster() +
                ggplot2::coord_quickmap(xlim=plot_extent[1:2], ylim=plot_extent[3:4]) +
                ggplot2::facet_wrap(~month, ncol=4) +
                ggplot2::labs(title=paste(variable,"Monthly Climatology",range_label)) +
                ggplot2::theme_minimal()
        f <- file.path(output_dir, paste0(safe_name(variable),"_monthly_climatology_",format(Sys.time(),"%Y%m%d_%H%M%S"),".png"))
        ggplot2::ggsave(f, p, width=12, height=8)
        f
}

plot_percentile_map <- function(percentile_df, plot_extent, variable, percentile, range_label = "Percentile", output_dir) {
        if (nrow(percentile_df) == 0) return(NULL)
        p <- ggplot2::ggplot(percentile_df, ggplot2::aes(x = x, y = y, fill = value)) +
                ggplot2::geom_raster() +
                ggplot2::coord_quickmap(xlim = plot_extent[1:2], ylim = plot_extent[3:4]) +
                safe_scale_fill_viridis() +
                ggplot2::labs(title = paste(variable, percentile, "Percentile", range_label)) +
                ggplot2::theme_minimal()
        f <- file.path(output_dir, paste0(safe_name(variable), "_percentile_", percentile, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"))
        ggplot2::ggsave(f, p, width = 8, height = 6)
        f
}

plot_heatwave_metrics <- function(hw_hist_df, hw_proj_df = NULL, plot_extent, variable, output_dir) {
        df <- hw_hist_df
        if (!is.null(hw_proj_df)) df <- dplyr::bind_rows(df, hw_proj_df)
        if (nrow(df) == 0) return(NULL)
        p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
                ggplot2::geom_raster() +
                ggplot2::coord_quickmap(xlim = plot_extent[1:2], ylim = plot_extent[3:4]) +
                ggplot2::facet_wrap(~layer) +
                ggplot2::labs(title = paste("Heatwave metrics for", variable)) +
                ggplot2::theme_minimal()
        f <- file.path(output_dir, paste0(safe_name(variable), "_heatwave_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"))
        ggplot2::ggsave(f, p, width = 10, height = 8)
        f
}

plot_threshold_days <- function(th_hist_df, th_proj_df = NULL, plot_extent, variable, threshold, output_dir) {
        df <- th_hist_df
        if (!is.null(th_proj_df)) df <- dplyr::bind_rows(df, th_proj_df)
        if (nrow(df) == 0) return(NULL)
        threshold <- as.numeric(threshold)
        p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, fill = value)) +
                ggplot2::geom_raster() +
                ggplot2::coord_quickmap(xlim = plot_extent[1:2], ylim = plot_extent[3:4]) +
                ggplot2::facet_wrap(~layer) +
                safe_scale_fill_viridis() +
                ggplot2::labs(title = paste("Threshold days >=", round(threshold - 273.15), "°C for", variable)) +
                ggplot2::theme_minimal()
        f <- file.path(output_dir, paste0(safe_name(variable), "_threshold_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"))
        ggplot2::ggsave(f, p, width = 10, height = 8)
        f
}

plot_annual_or_monthly <- function(r_df, plot_extent, variable, period = c("annual", "monthly"), output_dir) {
        period <- match.arg(period)
        if (nrow(r_df) == 0) return(NULL)
        r_df$period <- factor(r_df$layer, levels = unique(r_df$layer))
        p <- ggplot2::ggplot(r_df, ggplot2::aes(x = x, y = y, fill = value)) +
                ggplot2::geom_raster() +
                ggplot2::coord_quickmap(xlim = plot_extent[1:2], ylim = plot_extent[3:4]) +
                ggplot2::facet_wrap(~period) +
                ggplot2::labs(title = paste(toupper(substr(period,1,1)), substr(period,2,nchar(period)), "mean for", variable)) +
                ggplot2::theme_minimal()
        f <- file.path(output_dir, paste0(safe_name(variable), "_", period, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".png"))
        ggplot2::ggsave(f, p, width = 12, height = 8)
        f
}
