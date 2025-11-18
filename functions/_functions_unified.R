# ----------------------
# Raster Loading & Cropping
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
# ----------------------
load_and_crop_rasters <- function(files, var_token, start_year, end_year,
                                  test_mode = FALSE, city_shape = NULL,
                                  sample_layers = NULL) {
        if (is.null(files) || length(files) == 0) return(NULL)
        files_var <- files[stringr::str_detect(tolower(basename(files)), tolower(var_token))]
        if (length(files_var) == 0) return(NULL)
        
        rasters <- lapply(files_var, function(f) tryCatch(terra::rast(f), error = function(e) NULL))
        rasters <- purrr::compact(rasters)
        if (length(rasters) == 0) return(NULL)
        
        r_all <- do.call(c, rasters)
        
        # Sample layers if test_mode
        if (!is.null(sample_layers) && sample_layers < nlyr(r_all)) {
                set.seed(42)
                idx <- sample(seq_len(nlyr(r_all)), sample_layers)
                r_all <- r_all[[idx]]
        }
        
        # Crop/mask only if city_shape provided
        if (!is.null(city_shape)) {
                city_shape <- sf::st_transform(city_shape, crs = terra::crs(r_all))
                r_all <- terra::crop(r_all, terra::vect(city_shape))
                r_all <- terra::mask(r_all, terra::vect(city_shape))
        }
        
        return(r_all)
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
aggregate_time_period <- function(r, start_year, end_year, period = c("annual", "monthly")) {
        if (is.null(r) || nlyr(r) == 0) return(NULL)
        period <- match.arg(period)
        n_layers <- nlyr(r)
        
        if (n_layers < 12) {
                message("Skipping aggregation: fewer than 12 layers.")
                return(NULL)
        }
        
        if (period == "monthly") {
                groups <- rep(1:12, length.out = n_layers)
        } else {
                n_years <- ceiling(n_layers / 12)
                groups <- rep(1:n_years, each = 12)[1:n_layers]
        }
        
        grouped <- split(seq_len(n_layers), groups)
        out_list <- lapply(grouped, function(idx) terra::app(r[[idx]], mean, na.rm = TRUE))
        
        # Combine layers into one SpatRaster
        out <- do.call(c, out_list)
        if (!inherits(out, "SpatRaster")) out <- terra::rast(out_list)
        
        names(out) <- paste0(period, "_", seq_along(grouped))
        terra::ext(out) <- terra::ext(r)
        terra::crs(out) <- terra::crs(r)
        return(out)
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
        if (is.null(r)) return(data.frame())
        df <- as.data.frame(r, xy = TRUE)
        if (ncol(df) <= 2) return(df)
        tidyr::pivot_longer(df, cols = -c(x, y), names_to = "layer", values_to = "value")
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

plot_monthly_climatology <- function(monthly_df, plot_extent, variable, range_label, output_dir) {
        if (nrow(monthly_df) == 0) return(NULL)
        monthly_df$month <- factor(monthly_df$layer, levels=paste0("monthly_",1:12), labels=month.abb)
        title_text <- paste(variable, "Monthly Climatology", range_label, "-", scenarios[1], "-", models[1])
        p <- ggplot2::ggplot(monthly_df, ggplot2::aes(x=x,y=y,fill=value)) +
                ggplot2::geom_raster() +
                ggplot2::coord_quickmap(xlim=plot_extent[1:2], ylim=plot_extent[3:4]) +
                ggplot2::facet_wrap(~month, ncol=4) +
                ggplot2::scale_fill_viridis_c(option="C") +
                ggplot2::labs(title=title_text, fill="Value") +
                ggplot2::theme_minimal(base_size=14)+
                ggplot2::labs(x = "Longitude", y = "Latitude")+
                ggplot2::theme(panel.background = element_rect(fill="white"),
                               plot.background = element_rect(fill="white"))
        f <- file.path(output_dir, paste0(safe_name(variable),"_monthly_climatology_",format(Sys.time(),"%Y%m%d_%H%M%S"),".png"))
        ggplot2::ggsave(f, p, width=12, height=8, bg = "white")
        f
}

plot_annual_or_monthly <- function(r_df, plot_extent, variable, period=c("annual","monthly"), output_dir) {
        period <- match.arg(period)
        if (nrow(r_df) == 0) return(NULL)
        r_df$period <- factor(r_df$layer, levels=unique(r_df$layer))
        title_text <- paste(variable, toupper(period), "Climatology", "-", scenarios[2], "-", models[1])
        p <- ggplot2::ggplot(r_df, ggplot2::aes(x=x,y=y,fill=value)) +
                ggplot2::geom_raster() +
                ggplot2::coord_quickmap(xlim=plot_extent[1:2], ylim=plot_extent[3:4]) +
                ggplot2::facet_wrap(~period) +
                ggplot2::scale_fill_viridis_c(option="C") +
                ggplot2::labs(title=title_text, fill="Value") +
                ggplot2::theme_minimal(base_size=14)+
                ggplot2::labs(x = "Longitude", y = "Latitude")+
                ggplot2::theme(panel.background = element_rect(fill="white"),
                               plot.background = element_rect(fill="white"))
        f <- file.path(output_dir, paste0(safe_name(variable),"_",period,"_",format(Sys.time(),"%Y%m%d_%H%M%S"),".png"))
        ggplot2::ggsave(f, p, width=12, height=8, bg = "white")
        f
}

plot_percentile_map <- function(percentile_df, plot_extent, variable, percentile, range_label, output_dir) {
        if (nrow(percentile_df) == 0) return(NULL)
        title_text <- paste(variable, percentile, "Percentile", range_label, "-", models[1])
        p <- ggplot2::ggplot(percentile_df, ggplot2::aes(x=x,y=y,fill=value)) +
                ggplot2::geom_raster() +
                ggplot2::coord_quickmap(xlim=plot_extent[1:2], ylim=plot_extent[3:4]) +
                ggplot2::scale_fill_viridis_c(option="C") +
                ggplot2::labs(title=title_text, fill="Value") +
                ggplot2::theme_minimal(base_size=14)+
                ggplot2::labs(x = "Longitude", y = "Latitude")+
                ggplot2::theme(panel.background = element_rect(fill="white"),
                               plot.background = element_rect(fill="white"))
        f <- file.path(output_dir, paste0(safe_name(variable),"_percentile_",percentile,"_",format(Sys.time(),"%Y%m%d_%H%M%S"),".png"))
        ggplot2::ggsave(f, p, width=8, height=6, bg = "white")
        f
}

percentile_map <- function(r, prob) {
        if (is.null(r)) return(NULL)
        tryCatch(terra::app(r, fun = function(x) stats::quantile(x, probs = prob, na.rm = TRUE)), error = function(e) NULL)
}