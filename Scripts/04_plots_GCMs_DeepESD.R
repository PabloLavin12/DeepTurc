# ==============================================================================
# MULTI-MODEL HISTORICAL EVALUATION & GWL PROJECTIONS (DeepESD)
# 
# Description:
# This script evaluates the historical performance (1981-2005) of Global Climate 
# Models (GCMs) downscaled via DeepESD against ERA5 reanalysis and observations. 
# It computes spatial Climatologies, Mean Bias, and Interannual Correlations.
# 
# Additionally, it extracts and visualizes future anomalies (Deltas) for specific 
# Global Warming Levels (GWL: +2.0ºC, +3.0ºC, and +4.0ºC) under the RCP8.5 scenario.
# Finally, it computes and plots the Multi-Model Ensemble Mean for both historical 
# baselines and future projections, handling sub-daily data irregularities and 
# missing GWL thresholds (e.g., models that never reach +4.0ºC).
# ==============================================================================

# ==============================================================================
# SECTION 0: INITIAL SETUP & LIBRARIES
# ==============================================================================
library(visualizeR)
library(magrittr)
library(RColorBrewer)
library(gridExtra)
library(grid)

# List of models to process
gcms <- c("ACCESS1-0", "CanESM2", "MPI-ESM-MR", "NorESM1-M")

# Central Years Table for Global Warming Levels (GWL) under RCP8.5
gwl_data <- data.frame(
  GCM = gcms,
  GWL_2 = c(2041, 2026, 2039, 2048),
  GWL_3 = c(2062, 2049, 2060, 2072),
  GWL_4 = c(2081, 2068, 2082, NA) # Note: NorESM1-M does not reach +4.0ºC
)

# Standardized limits for visualizations (to allow visual comparison across models)
limite_bias  <- 2
limite_delta <- 5

# ==============================================================================
# SECTION 1: HELPER FUNCTIONS (AGGREGATION, CLEANING & METRICS)
# ==============================================================================

calcular_climatologia <- function(grid, mask_oceanos, nombre_dataset = "Dataset") {
  
  message(sprintf("\n--- Processing: %s ---", nombre_dataset))
  
  raw_dates <- as.character(unlist(grid$Dates))
  d_str <- grep("^[0-9]{4}-[0-9]{2}-[0-9]{2}", raw_dates, value = TRUE)
  
  n_time <- dim(grid$Data)[1]
  if (length(d_str) > n_time) {
    d_str <- d_str[1:n_time]
  }
  
  # Filter duplicates of sub-daily data of GCMs
  unique_days <- unique(d_str)
  if (length(unique_days) < length(d_str)) {
    message(sprintf("Detected %d records for %d unique days.", length(d_str), length(unique_days)))
    message("Collapsing duplicates by calculating the daily MEAN...")
    
    daily_data <- array(NA, dim = c(length(unique_days), dim(grid$Data)[2], dim(grid$Data)[3]))
    
    for (d in seq_len(length(unique_days))) {
      idx_day <- which(d_str == unique_days[d])
      if (length(idx_day) > 1) {
        # If there are multiple records for the same day, average them
        daily_data[d, , ] <- colMeans(grid$Data[idx_day, , , drop = FALSE], na.rm = TRUE, dims = 1)
      } else {
        daily_data[d, , ] <- grid$Data[idx_day, , ]
      }
    }
    
    # Replace raw data with the clean daily aggregated data
    grid$Data <- daily_data
    d_str <- unique_days
  }
  # =================================================================
  
  # -- Monthly Aggregation --
  message("Summing daily values to monthly resolution...")
  ym <- substr(d_str, 1, 7) 
  meses_unicos <- unique(ym)
  n_meses <- length(meses_unicos)
  
  mensual_data <- array(NA, dim = c(n_meses, dim(grid$Data)[2], dim(grid$Data)[3]))
  
  for (m in seq_len(n_meses)) {
    idx <- which(ym == meses_unicos[m])
    mensual_data[m, , ] <- colSums(grid$Data[idx, , , drop = FALSE], na.rm = TRUE, dims = 1)
  }
  
  grid_monthly <- grid
  grid_monthly$Data <- mensual_data
  grid_monthly$Dates$start <- unname(as.character(tapply(d_str, ym, function(x) x[1])))
  grid_monthly$Dates$end <- unname(as.character(tapply(d_str, ym, function(x) x[length(x)])))
  
  # -- Annual Aggregation --
  message("Summing monthly values to annual resolution...")
  years <- substr(grid_monthly$Dates$start, 1, 4) 
  years_unicos <- unique(years)
  n_years <- length(years_unicos)
  
  anual_data <- array(NA, dim = c(n_years, dim(grid_monthly$Data)[2], dim(grid_monthly$Data)[3]))
  
  for (y in seq_len(n_years)) {
    idx <- which(years == years_unicos[y])
    anual_data[y, , ] <- colSums(grid_monthly$Data[idx, , , drop = FALSE], na.rm = TRUE, dims = 1)
  }
  
  grid_annual <- grid_monthly
  grid_annual$Data <- anual_data
  grid_annual$Dates$start <- unname(as.character(tapply(grid_monthly$Dates$start, years, function(x) x[1])))
  grid_annual$Dates$end <- unname(as.character(tapply(grid_monthly$Dates$end, years, function(x) x[length(x)])))
  attr(grid_annual$Data, "dimensions") <- c("time", "lat", "lon")
  
  # -- Climatology --
  message("Calculating the climatological mean of the annual accumulations...")
  clim_data <- colMeans(grid_annual$Data, na.rm = TRUE, dims = 1)
  
  grid_clim <- grid_annual
  grid_clim$Data <- array(clim_data, dim = c(1, dim(grid_annual$Data)[2], dim(grid_annual$Data)[3]))
  grid_clim$Dates$start <- grid_annual$Dates$start[1]
  grid_clim$Dates$end <- grid_annual$Dates$end[length(grid_annual$Dates$end)]
  
  # -- Restore Ocean Mask --
  capa_clim <- grid_clim$Data[1, , ]
  capa_clim[mask_oceanos] <- NA
  grid_clim$Data[1, , ] <- capa_clim
  attr(grid_clim$Data, "dimensions") <- c("time", "lat", "lon")
  
  if(is.null(grid_clim$Variable$varName)) grid_clim$Variable$varName <- "Turc"
  
  return(list(clim = grid_clim, annual = grid_annual))
}

# ===================================================================
# FUNCTION TO CLEAN SUB-DAILY DATA AT THE SOURCE
# ===================================================================
limpiar_subdiario <- function(grid) {
  raw_dates <- as.character(unlist(grid$Dates))
  d_str <- grep("^[0-9]{4}-[0-9]{2}-[0-9]{2}", raw_dates, value = TRUE)
  n_time <- dim(grid$Data)[1]
  if (length(d_str) > n_time) d_str <- d_str[1:n_time]
  
  unique_days <- unique(d_str)
  
  if (length(unique_days) < length(d_str)) {
    message("Cleaning duplicate/sub-daily records to enforce daily resolution...")
    daily_data <- array(NA, dim = c(length(unique_days), dim(grid$Data)[2], dim(grid$Data)[3]))
    
    for (d in seq_len(length(unique_days))) {
      idx_day <- which(d_str == unique_days[d])
      if (length(idx_day) > 1) {
        daily_data[d, , ] <- colMeans(grid$Data[idx_day, , , drop = FALSE], na.rm = TRUE, dims = 1)
      } else {
        daily_data[d, , ] <- grid$Data[idx_day, , ]
      }
    }
    
    # Overwrite data and rebuild dates so intersectGrid operates smoothly
    grid$Data <- daily_data
    attr(grid$Data, "dimensions") <- c("time", "lat", "lon")
    grid$Dates$start <- paste0(unique_days, " 00:00:00 GMT")
    grid$Dates$end   <- paste0(unique_days, " 24:00:00 GMT")
  }
  return(grid)
}

calc_cor_interanual <- function(obs_annual, pred_annual, mask_oceanos, threshold = 0.05) {
  message("--- Calculating Interannual Correlation ---")
  years_obs <- substr(obs_annual$Dates$start, 1, 4)
  years_pred <- substr(pred_annual$Dates$start, 1, 4)
  
  common_years <- intersect(years_obs, years_pred)
  idx_obs <- match(common_years, years_obs)
  idx_pred <- match(common_years, years_pred)
  
  data_obs_matched <- obs_annual$Data[idx_obs, , , drop = FALSE]
  data_pred_matched <- pred_annual$Data[idx_pred, , , drop = FALSE]
  lat_n <- dim(obs_annual$Data)[2]; lon_n <- dim(obs_annual$Data)[3]
  
  cor_array <- matrix(NA, nrow = lat_n, ncol = lon_n)
  pval_array <- matrix(NA, nrow = lat_n, ncol = lon_n)
  
  for (i in 1:lat_n) {
    for (j in 1:lon_n) {
      if (mask_oceanos[i, j]) next
      obs_series <- data_obs_matched[, i, j]
      pred_series <- data_pred_matched[, i, j]
      valid_idx <- complete.cases(obs_series, pred_series)
      
      if (sum(valid_idx) >= 10) {
        test <- cor.test(pred_series[valid_idx], obs_series[valid_idx], method = "pearson")
        cor_array[i, j] <- test$estimate
        pval_array[i, j] <- test$p.value
      }
    }
  }
  cor_grid <- list(Data = cor_array, xyCoords = obs_annual$xyCoords, Variable = list(varName="Correlation"), Dates = obs_annual$Dates)
  attr(cor_grid$Data, "dimensions") <- c("lat", "lon")
  class(cor_grid) <- "grid"
  
  pval_grid <- list(Data = pval_array, xyCoords = obs_annual$xyCoords, Variable = list(varName="p-values"), Dates = obs_annual$Dates)
  attr(pval_grid$Data, "dimensions") <- c("lat", "lon")
  class(pval_grid) <- "grid"
  
  pts <- map.stippling(climatology(pval_grid), threshold = threshold, condition = "LT", pch = 19, col = "black", cex = 0.03) %>% suppressMessages() %>% suppressWarnings()
  return(list(cor_grid = climatology(cor_grid), pts = pts))
}


# =================================================================
# SIZE VARIABLES FOR TITLES AND LEGENDS (ADJUSTED FOR LATEX)
# =================================================================
cex_main <- 2.2
cex_na   <- 20
cex_top  <- 26
cex_legend   <- 1.6
width_legend <- 1.8

my_colorkey <- list(labels = list(cex = cex_legend), width = width_legend)

# ===================================================================
# SECTION 2: HISTORICAL COMPARISON (ERA5 + ALL GCMS)
# ===================================================================
message("\n=========================================")
message(" INITIATING MULTI-MODEL HISTORICAL BLOCK ")
message("=========================================\n")

# Load generic historical baselines
data_obs_raw <- readRDS('../../data/indexTurc_daily.rds')
attr(data_obs_raw$Data, "dimensions") <- c("time", "lat", "lon")
data_era5_raw <- readRDS('../../data/pred_DeepESD_CV.rds')

# Filter to standard historical period
years_target <- 1981:2005
data_obs_base <- subsetGrid(data_obs_raw, years = years_target)
data_era5_base <- subsetGrid(data_era5_raw, years = years_target)

# Intersect Obs with ERA5
data_obs_sync <- intersectGrid(data_era5_base, data_obs_base, type = "temporal", which.return = 2)
data_era5_sync <- intersectGrid(data_era5_base, data_obs_base, type = "temporal", which.return = 1)

# Generate Ocean Mask
suma_total <- colSums(data_era5_sync$Data, na.rm = TRUE, dims = 1)
mask_oceanos <- is.na(suma_total) | suma_total == 0

# Process ERA5
res_obs_era5 <- calcular_climatologia(data_obs_sync, mask_oceanos, "Obs (Sync with ERA5)")
res_era5 <- calcular_climatologia(data_era5_sync, mask_oceanos, "ERA5")

# ERA5 Bias and Correlation
clim_bias_era5 <- res_era5$clim
clim_bias_era5$Data <- res_era5$clim$Data - res_obs_era5$clim$Data
attr(clim_bias_era5$Data, "dimensions") <- c("time", "lat", "lon")
cor_era5 <- calc_cor_interanual(res_obs_era5$annual, res_era5$annual, mask_oceanos)


# Generate ERA5 plots
plot_list_hist <- list()
plot_list_hist[[1]] <- spatialPlot(res_era5$clim, backdrop.theme="countries", at=seq(0,40,length.out=21), main=list(label="ERA5", cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
plot_list_hist[[2]] <- spatialPlot(clim_bias_era5, backdrop.theme="countries", at=seq(-limite_bias, limite_bias, length.out=21), set.min=-limite_bias, set.max=limite_bias, main=list(label="ERA5 Bias", cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
plot_list_hist[[3]] <- spatialPlot(cor_era5$cor_grid, sp.layout=list(cor_era5$pts), backdrop.theme="countries", at=seq(-1, 1, length.out=21), main=list(label="ERA5 Cor", cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)

# Dictionary to store historical GCM climatologies
hist_clim_dict <- list()
hist_bias_list <- list()
hist_cor_list  <- list()

for (gcm in gcms) {
  ruta_hist <- paste0('../../data/proyecciones_ESD/Turc_DeepESD_', gcm, '_historical.rds')
  data_hist_raw <- readRDS(ruta_hist)
  data_hist_raw <- limpiar_subdiario(data_hist_raw)
  available_years <- unique(getYearsAsINDEX(data_hist_raw))
  safe_years <- intersect(years_target, available_years)
  if(length(safe_years) == 0) {
    message(paste("WARNING: No common years for", gcm, "in the historical period. Skipping..."))
    next
  }
  
  message(sprintf("Safe years applied for %s: %d to %d", gcm, min(safe_years), max(safe_years)))
  data_hist_sub <- subsetGrid(data_hist_raw, years = safe_years)
  
  # Temporal synchronization with GCM
  data_obs_gcm <- intersectGrid(data_hist_sub, data_obs_base, type="temporal", which.return=2)
  data_hist_gcm <- intersectGrid(data_hist_sub, data_obs_base, type="temporal", which.return=1)
  
  res_obs_gcm <- calcular_climatologia(data_obs_gcm, mask_oceanos, paste("Obs_sync", gcm))
  res_hist_gcm <- calcular_climatologia(data_hist_gcm, mask_oceanos, gcm)
  
  # Bias
  clim_bias_gcm <- res_hist_gcm$clim
  clim_bias_gcm$Data <- res_hist_gcm$clim$Data - res_obs_gcm$clim$Data
  attr(clim_bias_gcm$Data, "dimensions") <- c("time", "lat", "lon")
  
  # Correlation
  cor_gcm <- calc_cor_interanual(res_obs_gcm$annual, res_hist_gcm$annual, mask_oceanos)
  
  # Save pure matrices for the Ensemble Mean
  hist_clim_dict[[gcm]] <- res_hist_gcm$clim
  hist_bias_list[[gcm]] <- clim_bias_gcm$Data
  hist_cor_list[[gcm]]  <- cor_gcm$cor_grid$Data
  
  # GCM Plots
  p_pred <- spatialPlot(res_hist_gcm$clim, backdrop.theme="countries", at=seq(0,40,length.out=21), main=list(label=paste(gcm), cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
  p_bias <- spatialPlot(clim_bias_gcm, backdrop.theme="countries", at=seq(-limite_bias, limite_bias, length.out=21), set.min=-limite_bias, set.max=limite_bias, main=list(label=paste(gcm, "Bias"), cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
  p_cor <- spatialPlot(cor_gcm$cor_grid, sp.layout=list(cor_gcm$pts), backdrop.theme="countries", at=seq(-1, 1, length.out=21), main=list(label=paste(gcm, "Cor"), cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
  
  plot_list_hist <- c(plot_list_hist, list(p_pred, p_bias, p_cor))
}

# =================================================================
# HISTORICAL ENSEMBLE MEAN CALCULATION
# =================================================================
message("Calculating Historical Ensemble Mean...")
# Calculate the mean of the matrices
ens_clim_data <- Reduce("+", lapply(hist_clim_dict, function(x) x$Data)) / length(hist_clim_dict)
ens_bias_data <- Reduce("+", hist_bias_list) / length(hist_bias_list)
ens_cor_data  <- Reduce("+", hist_cor_list) / length(hist_cor_list)

# Use the first GCM as a "template" to inject the ensemble data
ens_template <- hist_clim_dict[[1]]

ens_clim_grid <- ens_template; ens_clim_grid$Data <- ens_clim_data
ens_bias_grid <- ens_template; ens_bias_grid$Data <- ens_bias_data
ens_cor_grid  <- ens_template; ens_cor_grid$Data  <- ens_cor_data

p_ens_pred <- spatialPlot(ens_clim_grid, backdrop.theme="countries", at=seq(0,40,length.out=21), main=list(label="Ensemble Mean", cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
p_ens_bias <- spatialPlot(ens_bias_grid, backdrop.theme="countries", at=seq(-limite_bias, limite_bias, length.out=21), set.min=-limite_bias, set.max=limite_bias, main=list(label="Ensemble Mean Bias", cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
p_ens_cor  <- spatialPlot(ens_cor_grid, backdrop.theme="countries", at=seq(-1, 1, length.out=21), main=list(label="Ensemble Mean Cor", cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)

plot_list_hist <- c(plot_list_hist, list(p_ens_pred, p_ens_bias, p_ens_cor))

# Export Historical Comparison PDF
message(">>> Exporting Historical Comparison PDF...")
pdf(file = "Comparativa_Turc_Hist_Deep.pdf", width = 16, height = 25)
grid.arrange(grobs = plot_list_hist, ncol = 3, nrow = length(gcms) + 2, 
             top = textGrob("Historical Comparison (1981-2005) of the Turc Index: ERA5 and CMIP5 GCMs Models", gp = gpar(fontsize = cex_top, font = 2)))
dev.off()

png(file = "Comparativa_Turc_Hist_Deep.png", 
    width = 16, height = 22, units = "in", res = 300)

grid.arrange(grobs = plot_list_hist, ncol = 3, nrow = length(gcms) + 2, 
             top = textGrob("Historical Comparison (1981-2005) of the Turc Index: ERA5 and CMIP5 GCMs Models", gp = gpar(fontsize = cex_top, font = 2)))
dev.off()

# ==============================================================================
# SECTION 3: FUTURE GWL PROJECTIONS & ENSEMBLES
# ==============================================================================
message("\n=========================================")
message(" INITIATING GWL PROJECTIONS BLOCK        ")
message("=========================================\n")

plot_list_proj <- list()
niveles_gwl <- c("GWL_2", "GWL_3", "GWL_4")
grados <- c("2ºC", "3ºC", "4ºC")

# List to store anomaly matrices per GWL level
deltas_gwl_list <- list(GWL_2 = list(), GWL_3 = list(), GWL_4 = list())

for (i in seq_along(gcms)) {
  gcm <- gcms[i]
  ruta_rcp <- paste0('../../data/proyecciones_ESD/Turc_DeepESD_', gcm, '_rcp85.rds')
  
  if (file.exists(ruta_rcp)) {
    data_rcp_raw <- readRDS(ruta_rcp)
    clim_hist_gcm <- hist_clim_dict[[gcm]] 
    
    for (j in seq_along(niveles_gwl)) {
      lvl <- niveles_gwl[j]
      año_central <- gwl_data[i, lvl]
      
      if (is.na(año_central)) {
        p_na <- textGrob(paste(gcm, "\nDoes not reach", grados[j]), gp = gpar(fontsize=cex_na, col="gray50", fontface="italic"))
        plot_list_proj <- c(plot_list_proj, list(p_na))
      } else {
        gwl_years <- (año_central - 9):(año_central + 10)
        data_gwl <- subsetGrid(data_rcp_raw, years = gwl_years)
        
        res_gwl <- calcular_climatologia(data_gwl, mask_oceanos, paste(gcm, "-", grados[j]))
        
        delta_gwl <- res_gwl$clim
        delta_gwl$Data <- res_gwl$clim$Data - clim_hist_gcm$Data
        attr(delta_gwl$Data, "dimensions") <- c("time", "lat", "lon")
        
        # Save matrix for the ensemble of this GWL
        deltas_gwl_list[[lvl]][[gcm]] <- delta_gwl$Data
        
        titulo_mapa <- paste(gcm, "- Delta", grados[j], "\n(", gwl_years[1], "-", gwl_years[20], ")")
        
        # Added COLORKEY
        p_delta <- spatialPlot(delta_gwl, backdrop.theme = "countries",
                               at = seq(-limite_delta, limite_delta, length.out = 21), 
                               set.min = -limite_delta, set.max = limite_delta,
                               main = list(label=titulo_mapa, cex=cex_main),
                               color.theme = "RdYlBu", rev.colors = TRUE,
                               colorkey=my_colorkey)
        
        plot_list_proj <- c(plot_list_proj, list(p_delta))
      }
    }
    rm(data_rcp_raw)
    gc()
  } else {
    message(paste("File not found for:", gcm))
  }
}

# =================================================================
# FUTURE ENSEMBLE MEAN CALCULATION (GWL)
# =================================================================
message("Calculating Projections Ensemble Mean...")

for (j in seq_along(niveles_gwl)) {
  lvl <- niveles_gwl[j]
  datos_nivel <- deltas_gwl_list[[lvl]]
  
  if (length(datos_nivel) > 0) {
    # Average models that actually reached this GWL
    ens_delta_data <- Reduce("+", datos_nivel) / length(datos_nivel)
    
    # Use a template from the first part
    ens_delta_grid <- hist_clim_dict[[1]]
    ens_delta_grid$Data <- ens_delta_data
    
    titulo_ens <- paste("Ensemble Mean - Delta", grados[j], "\n(", length(datos_nivel), "models)")
    
    # Added COLORKEY
    p_ens_delta <- spatialPlot(ens_delta_grid, backdrop.theme = "countries",
                               at = seq(-limite_delta, limite_delta, length.out = 21), 
                               set.min = -limite_delta, set.max = limite_delta,
                               main = list(label=titulo_ens, cex=cex_main),
                               color.theme = "RdYlBu", rev.colors = TRUE,
                               colorkey=my_colorkey)
    plot_list_proj <- c(plot_list_proj, list(p_ens_delta))
  } else {
    p_na <- textGrob(paste("Ensemble\nNo data for", grados[j]), gp = gpar(fontsize=cex_na, col="gray50", fontface="italic"))
    plot_list_proj <- c(plot_list_proj, list(p_na))
  }
}

# Export Projections PDF
message(">>> Exporting GWL Projections PDF...")
pdf(file = "Proyecciones_Deltas_Deep.pdf", width = 16, height = 22)
grid.arrange(grobs = plot_list_proj, ncol = 3, nrow = length(gcms) + 1, 
             top = textGrob("Turc Index Anomalies by Global Warming Levels (GWL) in the RCP8.5 vs Historical Scenario", gp = gpar(fontsize = cex_top, font = 2)))
dev.off()

png(file = "Proyecciones_Deltas_Deep.png", 
    width = 16, height = 22, units = "in", res = 300)

grid.arrange(grobs = plot_list_proj, ncol = 3, nrow = length(gcms) + 1, 
             top = textGrob("Turc Index Anomalies by Global Warming Levels (GWL) in the RCP8.5 vs Historical Scenario", gp = gpar(fontsize = cex_top, font = 2)))
dev.off()

message("\n>>> PROCESS COMPLETED!")