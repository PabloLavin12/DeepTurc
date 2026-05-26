# ==============================================================================
# COMPREHENSIVE TURC INDEX ANALYSIS PIPELINE
# 
# Description:
# This unified script performs a complete end-to-end analysis of the Turc 
# Agricultural Potential Index using observational data (PTI), reanalysis (ERA5), 
# and downscaled GCM predictions (DeepESD & KNN).
# 
# The workflow is divided into three main blocks:
#   PART 1: Index Calculation & Spatial Methodological Comparison 
#           (Adapted Daily vs Classic Monthly formulations).
#   PART 2: Regional Time Series Validation
#           (Extracting and evaluating DeepESD and KNN predictions across 6 
#           specific Spanish regions with embedded performance metrics).
#   PART 3: Global Climate Drivers Analysis
#           (Evaluating the temporal evolution of the Turc Index and its 
#           seasonal correlation with Mean Temperature and Precipitation).
# ==============================================================================

# ==============================================================================
# SECTION 0: INITIAL SETUP & LIBRARIES
# ==============================================================================
options(java.parameters = "-Xmx16g")

library(dplyr)
library(abind)
library(loadeR)
library(transformeR)
library(visualizeR)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)
library(grid)

# Create output directories if they do not exist
dir.create("series", showWarnings = FALSE)

# ==============================================================================
# SECTION 1: CORE FUNCTIONS
# ==============================================================================

#' Calculate Daily Turc Index (Adapted for daily climate data)
#' @param data_list List containing grids: tmax, tmin, hr, pr, ssrd
#' @param wfc Soil field capacity (default 100 mm)
#' @param source_type String to define radiation conversion: "PTI", "ERA5", or "GCM"
calc_turc_daily <- function(data_list, wfc = 100, source_type = "PTI") {
  ntime <- dim(data_list$tmax$Data)[1]
  nlat  <- dim(data_list$tmax$Data)[2]
  nlon  <- dim(data_list$tmax$Data)[3]
  lats  <- data_list$tmax$xyCoords$y
  
  CA_3d <- array(NA, dim = c(ntime, nlat, nlon))
  
  lat_2d <- matrix(rep(lats, nlon), nrow = nlat, ncol = nlon)
  lat_rad <- lat_2d * pi / 180
  
  fechas <- as.Date(data_list$tmax$Dates$start)
  jdays <- as.numeric(format(fechas, "%j"))
  
  R <- matrix(wfc, nrow = nlat, ncol = nlon) 
  R[is.na(data_list$pr$Data[1,,])] <- NA     
  ET <- matrix(0, nrow = nlat, ncol = nlon)  
  
  for (i in 1:ntime) {
    tmax_i <- data_list$tmax$Data[i,,]
    tmin_i <- data_list$tmin$Data[i,,]
    hr_i   <- data_list$hr$Data[i,,]
    pr_i   <- data_list$pr$Data[i,,]
    ssrd_i <- data_list$ssrd$Data[i,,]
    
    # Radiation adjustment based on data source
    if (source_type == "GCM") {
      rs_i <- ssrd_i * 2.065 # GCM: W/m2 to cal/cm2/day (Langleys)
    } else if (source_type == "ERA5") {
      rs_i <- ((ssrd_i * 24) / 1000) / 41.84 # ERA5: Hourly J to cal/cm2/day
    } else {
      rs_i <- ssrd_i / 41.84 # PTI: standard conversion
    }
    
    # Thermal Factor (Ft)
    tmean_i <- (tmax_i + tmin_i) / 2
    ft1 <- (tmean_i * (60 - tmean_i) / 1000) * ((tmin_i - 1) / 4)
    ft2 <- (tmean_i * (60 - tmean_i) / 1000)
    ft <- ifelse(tmin_i <= 1, 0, ifelse(tmin_i < 5, ft1, ft2))
    ft <- pmax(ft, 0, na.rm = FALSE)
    
    # Solar Factor (Fh)
    delta <- 0.409 * sin(0.0172 * jdays[i] - 1.39)
    arg_acos <- pmin(1, pmax(-1, -tan(lat_rad) * tan(delta)))
    omegas <- acos(arg_acos)
    N <- 7.64 * omegas
    fh1 <- N - 5 - ((lat_2d / 40)^2)
    fh2 <- 0.03 * (rs_i - 100)
    fh <- pmax(0, pmin(fh1, fh2), na.rm = FALSE)
    
    # Potential Evapotranspiration (ETo)
    eto1 <- 0.0133 * (tmean_i / (tmean_i + 15)) * (rs_i + 50)
    eto2 <- eto1 * (1 + (50 - hr_i) / 70)
    eto <- pmax(0, ifelse(hr_i > 50, eto1, eto2), na.rm = FALSE)
    
    # Water balance updates
    if (i > 1) {
      pr_prev <- data_list$pr$Data[i-1,,]
      R <- pmax(0, pmin(100, R + pr_prev - ET))
    }
    DIF <- R + pr_i - eto
    ET <- ifelse(DIF > 0, eto, R + pr_i)
    
    # Dryness Factor (Fs)
    X <- pmin(eto, eto * 0.3 + 1.64)
    Deficit <- eto - ET
    fs_i <- ifelse(X == 0 & Deficit == 0, 0, pmax((X - Deficit) / X, 0))
    
    CA_3d[i,,] <- (ft * fh * fs_i) / 30.4
  }
  
  turc_grid <- data_list$tmax
  turc_grid$Data <- CA_3d
  turc_grid$Variable$varName <- "turc_daily"
  turc_grid$Variable$description <- "Daily Turc Agricultural Potential Index"
  turc_grid$Variable$units <- "index"
  attr(turc_grid$Data, "dimensions") <- c("time", "lat", "lon")
  
  return(turc_grid)
}

#' Calculate Classic Monthly Turc Index
calc_turc_monthly_classic <- function(data_list, wfc = 100) {
  ntime <- dim(data_list$tmax$Data)[1] 
  nlat  <- dim(data_list$tmax$Data)[2]
  nlon  <- dim(data_list$tmax$Data)[3]
  lats  <- data_list$tmax$xyCoords$y
  
  CA_3d <- array(NA, dim = c(ntime, nlat, nlon))
  lat_2d <- matrix(lats, nrow = nlat, ncol = nlon, byrow = FALSE)
  lat_rad <- lat_2d * pi / 180
  
  fechas <- as.Date(data_list$tmax$Dates$start)
  meses <- as.numeric(format(fechas, "%m"))
  midpoints_ref <- c(15, 46, 74, 104, 135, 165, 196, 227, 257, 288, 318, 349)
  jdays_mensuales <- midpoints_ref[meses]
  
  R <- matrix(wfc, nrow = nlat, ncol = nlon) 
  R[is.na(data_list$pr$Data[1,,])] <- NA     
  ET <- matrix(0, nrow = nlat, ncol = nlon)  
  
  for (i in 1:ntime) {
    tmax_i <- data_list$tmax$Data[i,,]
    tmin_i <- data_list$tmin$Data[i,,]
    hr_i   <- data_list$hr$Data[i,,]
    pr_i   <- data_list$pr$Data[i,,]
    ssrd_i <- data_list$ssrd$Data[i,,]
    
    rs_i <- ssrd_i / 41.84 
    
    tmean_i <- (tmax_i + tmin_i) / 2
    ft1 <- (tmean_i * (60 - tmean_i) / 1000) * ((tmin_i - 1) / 4)
    ft2 <- (tmean_i * (60 - tmean_i) / 1000)
    ft <- ifelse(tmin_i <= 1, 0, ifelse(tmin_i < 5, ft1, ft2))
    ft <- pmax(ft, 0, na.rm = FALSE)
    
    delta <- 0.409 * sin(0.0172 * jdays_mensuales[i] - 1.39)
    arg_acos <- pmin(1, pmax(-1, -tan(lat_rad) * tan(delta)))
    omegas <- acos(arg_acos)
    N <- 7.64 * omegas
    fh1 <- N - 5 - ((lat_2d / 40)^2)
    fh2 <- 0.03 * (rs_i - 100)
    fh <- pmax(0, pmin(fh1, fh2), na.rm = FALSE)
    
    eto1 <- 0.4 * (tmean_i / (tmean_i + 15)) * (rs_i + 50)
    eto2 <- eto1 * (1 + (50 - hr_i) / 70)
    eto <- pmax(0, ifelse(hr_i > 50, eto1, eto2), na.rm = FALSE)
    
    if (i > 1) {
      pr_prev <- data_list$pr$Data[i-1,,]
      R <- pmax(0, pmin(100, R + pr_prev - ET))
    }
    
    DIF <- R + pr_i - eto
    ET <- ifelse(DIF > 0, eto, R + pr_i)
    
    X <- pmin(eto, eto * 0.3 + 50)
    Deficit <- eto - ET
    fs_i <- ifelse(X == 0 & Deficit == 0, 0, pmax((X - Deficit) / X, 0))
    
    CA_3d[i,,] <- (ft * fh * fs_i)
  }
  
  turc_grid <- data_list$tmax
  turc_grid$Data <- CA_3d
  turc_grid$Variable$varName <- "Turc_Original"
  turc_grid$Variable$description <- "Monthly Turc Agricultural Potential Index (Original)"
  attr(turc_grid$Data, "dimensions") <- c("time", "lat", "lon") 
  
  return(turc_grid)
}

# GCM unit conversion
convert_gcm_units <- function(grid, value, operation, new_unit) {
  new_grid <- gridArithmetics(grid, value, operator = operation)
  new_grid$Variable <- grid$Variable
  new_grid$Variable$units <- new_unit
  return(new_grid)
}

# Universal Sea Mask
apply_sea_mask <- function(grid, mask) {
  dims <- length(dim(grid$Data))
  if (dims == 2) {
    grid$Data[mask] <- NA
  } else if (dims == 3) {
    for (i in 1:dim(grid$Data)[1]) {
      temp_layer <- grid$Data[i, , ]  
      temp_layer[mask] <- NA          
      grid$Data[i, , ] <- temp_layer  
    }
  } else if (dims == 4) {
    for (m in 1:dim(grid$Data)[1]) {
      for (t in 1:dim(grid$Data)[2]) {
        temp_layer <- grid$Data[m, t, , ]
        temp_layer[mask] <- NA
        grid$Data[m, t, , ] <- temp_layer
      }
    }
  }
  return(grid)
}

# Fast Annual Aggregation
agregacion_anual_rapida <- function(grid_diario) {
  fechas <- as.Date(grid_diario$Dates$start)
  years <- format(fechas, "%Y")
  years_unicos <- unique(years)
  n_years <- length(years_unicos)
  n_lat <- dim(grid_diario$Data)[2]; n_lon <- dim(grid_diario$Data)[3]
  
  data_anual <- array(NA, dim = c(n_years, n_lat, n_lon))
  for (y in 1:n_years) {
    idx <- which(years == years_unicos[y])
    data_anual[y, , ] <- colSums(grid_diario$Data[idx, , , drop = FALSE], na.rm = TRUE, dims = 1)
  }
  
  grid_anual <- grid_diario
  grid_anual$Data <- data_anual
  attr(grid_anual$Data, "dimensions") <- c("time", "lat", "lon")
  grid_anual$Dates$start <- as.character(tapply(grid_diario$Dates$start, years, function(x) x[1]))
  grid_anual$Dates$end <- as.character(tapply(grid_diario$Dates$end, years, function(x) x[length(x)]))
  return(grid_anual)
}

# Coordinate Time Series Extraction
extraer_serie_punto <- function(grid_obs, grid_pred, lon_target, lat_target, nombre_region) {
  lons <- grid_obs$xyCoords$x
  lats <- grid_obs$xyCoords$y
  
  idx_lon <- which.min(abs(lons - lon_target))
  idx_lat <- which.min(abs(lats - lat_target))
  
  lon_real <- lons[idx_lon]
  lat_real <- lats[idx_lat]
  
  fechas <- as.Date(grid_obs$Dates$start)
  
  ts_obs <- grid_obs$Data[, idx_lat, idx_lon]
  ts_pred <- grid_pred$Data[, idx_lat, idx_lon]
  
  df <- data.frame(
    Fecha = fechas,
    Region = nombre_region,
    Lon_Target = lon_target, Lat_Target = lat_target,
    Lon_Real = lon_real, Lat_Real = lat_real,
    Observado = ts_obs,
    Modelo = ts_pred,
    stringsAsFactors = FALSE
  )
  return(df)
}

# Fast Spatial Mean Calculator
get_spatial_mean <- function(grid_data, mask) {
  if (length(dim(grid_data)) == 3) {
    apply(grid_data, 1, function(layer) mean(layer[!mask], na.rm = TRUE))
  } else if (length(dim(grid_data)) == 4) {
    apply(grid_data[1, , , ], 1, function(layer) mean(layer[!mask], na.rm = TRUE))
  }
}


# ==============================================================================
# PART 1: INDEX CALCULATION & METHODOLOGICAL COMPARISON
# ==============================================================================
message("\n--- PART 1: INDEX CALCULATION & SPATIAL CLIMATOLOGIES ---")

message("Processing PTI observational data...")
data_list_pti <- readRDS('../../data/PTI_grid_5km_1981_2020.rds')
turc_daily_pti <- calc_turc_daily(data_list = data_list_pti, source_type = "PTI")
turc_monthly_pti <- aggregateGrid(turc_daily_pti, aggr.m = list(FUN = "sum", na.rm = TRUE))
turc_annual_pti  <- aggregateGrid(turc_daily_pti, aggr.y = list(FUN = "sum", na.rm = TRUE))
pti_grid_ref <- getGrid(turc_daily_pti) 

message("Processing ERA5 reanalysis variables...")
data_list_era5 <- list(
  tmax = readRDS("../../data/ERA5-Turc/t2mx_ERA5_daily_1981_2020_5km.rds"),
  tmin = readRDS("../../data/ERA5-Turc/t2mn_ERA5_daily_1981_2020_5km.rds"),
  hr   = readRDS("../../data/ERA5-Turc/hurs_ERA5_daily_1981_2020_5km.rds"),
  pr   = readRDS("../../data/ERA5-Turc/tp_ERA5_daily_1981_2020_5km.rds"),
  ssrd = readRDS("../../data/ERA5-Turc/ssrd_ERA5_daily_1981_2020_5km.rds")
)
turc_daily_era5 <- calc_turc_daily(data_list = data_list_era5, source_type = "ERA5")
turc_monthly_era5 <- aggregateGrid(turc_daily_era5, aggr.m = list(FUN = "sum", na.rm = TRUE))
turc_annual_era5  <- aggregateGrid(turc_daily_era5, aggr.y = list(FUN = "sum", na.rm = TRUE))

message("Processing GCM variables (CanESM2)...")
gcm_raw <- readRDS("../../data/GCMs/GCM_data_Surface_CanESM2.rds")
hist_data <- gcm_raw$historical
data_list_gcm_5km <- list(
  tmax = interpGrid(convert_gcm_units(hist_data$tasmax, 273.15, "-", "degC"), new.coordinates = pti_grid_ref, method = "bilinear"),
  tmin = interpGrid(convert_gcm_units(hist_data$tasmin, 273.15, "-", "degC"), new.coordinates = pti_grid_ref, method = "bilinear"),
  hr   = interpGrid(hist_data$hur, new.coordinates = pti_grid_ref, method = "bilinear"),
  pr   = interpGrid(convert_gcm_units(hist_data$pr, 86400, "*", "mm/day"), new.coordinates = pti_grid_ref, method = "bilinear"),
  ssrd = interpGrid(hist_data$rsds, new.coordinates = pti_grid_ref, method = "bilinear")
)
turc_daily_gcm <- calc_turc_daily(data_list = data_list_gcm_5km, source_type = "GCM")
turc_annual_gcm <- aggregateGrid(turc_daily_gcm, aggr.y = list(FUN = "sum", na.rm = TRUE))

message("Calculating Classic Monthly baseline using aggregated PTI data...")
data_list_monthly_pti <- list(
  tmax = aggregateGrid(data_list_pti$tmax, aggr.m = list(FUN = "mean", na.rm = TRUE)),
  tmin = aggregateGrid(data_list_pti$tmin, aggr.m = list(FUN = "mean", na.rm = TRUE)),
  hr   = aggregateGrid(data_list_pti$hr,   aggr.m = list(FUN = "mean", na.rm = TRUE)),
  ssrd = aggregateGrid(data_list_pti$ssrd, aggr.m = list(FUN = "mean", na.rm = TRUE)),
  pr   = aggregateGrid(data_list_pti$pr,   aggr.m = list(FUN = "sum",  na.rm = TRUE)) 
)
turc_original_monthly <- calc_turc_monthly_classic(data_list = data_list_monthly_pti)

message("Calculating Climatologies and Spatial Bias...")
clim_pti  <- climatology(turc_annual_pti)
clim_era5 <- climatology(turc_annual_era5)
clim_gcm  <- climatology(turc_annual_gcm)
mask_na <- is.na(clim_pti$Data[1, , ]) | clim_pti$Data[1, , ] == 0 | is.na(clim_gcm$Data[1, , ])

clim_pti  <- apply_sea_mask(clim_pti, mask_na)
clim_era5 <- apply_sea_mask(clim_era5, mask_na)
clim_gcm  <- apply_sea_mask(clim_gcm, mask_na)

bias_annual_era5 <- gridArithmetics(clim_era5, clim_pti, operator = "-")
bias_annual_gcm  <- gridArithmetics(clim_gcm, clim_pti, operator = "-")

# Methodological Comparisons
annual_orig <- aggregateGrid(turc_original_monthly, aggr.y = list(FUN = "sum", na.rm = TRUE))
annual_adap <- aggregateGrid(turc_monthly_pti, aggr.y = list(FUN = "sum", na.rm = TRUE))
clim_annual_orig <- apply_sea_mask(climatology(annual_orig), mask_na)
clim_annual_adap <- apply_sea_mask(climatology(annual_adap), mask_na)
clim_annual_adap$Dates <- clim_annual_orig$Dates 
bias_annual_method <- gridArithmetics(clim_annual_adap, clim_annual_orig, operator = "-")

seasons_list <- list(DJF = c(12, 1, 2), MAM = c(3, 4, 5), JJA = c(6, 7, 8), SON = c(9, 10, 11))
list_clim_orig <- list(); list_clim_adap <- list(); list_bias_est <- list()

for (est in names(seasons_list)) {
  months <- seasons_list[[est]]
  c_orig <- apply_sea_mask(climatology(subsetGrid(turc_original_monthly, season = months)), mask_na)
  c_adap <- apply_sea_mask(climatology(subsetGrid(turc_monthly_pti, season = months)), mask_na)
  list_clim_orig[[paste("Original", est)]] <- c_orig
  list_clim_adap[[paste("Adapted", est)]]  <- c_adap
  list_bias_est[[paste("Bias", est)]] <- gridArithmetics(c_adap, c_orig, operator = "-")
}

list_seasonal_abs <- c(list_clim_orig, list_clim_adap)
dates_ref_est  <- list_seasonal_abs[[1]]$Dates
dates_ref_bias <- list_bias_est[[1]]$Dates
list_seasonal_abs <- lapply(list_seasonal_abs, function(g) { g$Dates <- dates_ref_est; return(g) })
list_bias_est     <- lapply(list_bias_est, function(g) { g$Dates <- dates_ref_bias; return(g) })

multigrid_annual   <- bindGrid(clim_annual_orig, clim_annual_adap, dimension = "member")
names_annual       <- c("Classic Monthly", "Adapted Daily")
multigrid_seasonal <- do.call(bindGrid, c(list_seasonal_abs, list(dimension = "member")))
names_seasonal     <- names(list_seasonal_abs)
multigrid_bias_est <- do.call(bindGrid, c(list_bias_est, list(dimension = "member")))
names_bias_est     <- names(list_bias_est)

# Mapping exports
message("Exporting Maps to PDF/PNG...")
limit_bias_models <- 5
cex_main <- 2.2 
cex_strip <- 1.6 
cex_legend <- 1.6   
width_legend <- 1.8 
my_colorkey <- list(labels = list(cex = cex_legend), width = width_legend)

plot_era5 <- function() {
  p1 <- spatialPlot(grid = clim_pti, backdrop.theme = "countries", main = list(label = "Annual Climatology Turc (PTI)", cex = cex_main), color.theme = "RdYlBu", rev.colors = TRUE, at = seq(0, 40, length.out = 21), colorkey = my_colorkey) 
  p2 <- spatialPlot(grid = clim_era5, backdrop.theme = "countries", main = list(label = "Annual Climatology Turc (ERA5)", cex = cex_main), color.theme = "RdYlBu", rev.colors = TRUE, at = seq(0, 40, length.out = 21), colorkey = my_colorkey) 
  p3 <- spatialPlot(grid = bias_annual_era5, backdrop.theme = "countries", main = list(label = "Bias (ERA5 - PTI)", cex = cex_main), color.theme = "RdYlBu", rev.colors = TRUE, at = seq(-limit_bias_models, limit_bias_models, length.out = 21), set.min = -limit_bias_models, set.max = limit_bias_models, colorkey = my_colorkey) 
  print(p1, split = c(1, 1, 3, 1), more = TRUE); print(p2, split = c(2, 1, 3, 1), more = TRUE); print(p3, split = c(3, 1, 3, 1), more = FALSE)
}
pdf(file = "Comparativa_Turc_ERA5_vs_PTI.pdf", width = 20, height = 7); plot_era5(); dev.off()
png(filename = "Comparativa_Turc_ERA5_vs_PTI.png", width = 20, height = 7, units = "in", res = 300); plot_era5(); dev.off()

plot_combined_method <- function() {
  max_val_ann <- max(abs(bias_annual_method$Data), na.rm = TRUE)
  lim_bias_ann <- seq(-max_val_ann, max_val_ann, length.out = 11)
  max_val_est <- max(abs(multigrid_bias_est$Data), na.rm = TRUE)
  lim_bias_est <- seq(-max_val_est, max_val_est, length.out = 11)
  
  p_ann_abs <- spatialPlot(multigrid_annual, names.attr = names_annual, backdrop.theme = "countries", main = list(label = "Annual Climatology", cex = cex_main), par.strip.text = list(cex = cex_strip), color.theme = "RdYlBu", rev.colors = TRUE, layout = c(2, 1), colorkey = my_colorkey) 
  p_ann_bias <- spatialPlot(bias_annual_method, backdrop.theme = "countries", main = list(label = "Annual Bias (Adapted - Classic)", cex = cex_main), color.theme = "RdYlBu", rev.colors = TRUE, at = lim_bias_ann, colorkey = my_colorkey) 
  p_est_abs <- spatialPlot(multigrid_seasonal, names.attr = names_seasonal, backdrop.theme = "countries", main = list(label = "Seasonal Climatology", cex = cex_main), par.strip.text = list(cex = cex_strip), color.theme = "RdYlBu", rev.colors = TRUE, layout = c(4, 2), colorkey = my_colorkey) 
  p_est_bias <- spatialPlot(multigrid_bias_est, names.attr = names_bias_est, backdrop.theme = "countries", main = list(label = "Seasonal Bias (Adapted - Classic)", cex = cex_main), par.strip.text = list(cex = cex_strip), color.theme = "RdYlBu", rev.colors = TRUE, layout = c(4, 1), at = lim_bias_est, colorkey = my_colorkey) 
  
  top_plot <- arrangeGrob(p_est_abs, p_est_bias, nrow = 2, heights = c(2, 1))
  bottom_plot <- arrangeGrob(p_ann_abs, p_ann_bias, ncol = 2, widths = c(2, 1))
  grid.arrange(top_plot, bottom_plot, nrow = 2, heights = c(2.5, 1))
}
pdf("Turc_Method_Comparison_Combined.pdf", width = 18, height = 18); plot_combined_method(); dev.off()
png("Turc_Method_Comparison_Combined.png", width = 18, height = 18, units = "in", res = 300); plot_combined_method(); dev.off()


# ==============================================================================
# PART 2: REGIONAL TIME SERIES VALIDATION (DEEP ESD vs KNN)
# ==============================================================================
message("\n--- PART 2: REGIONAL TIME SERIES VALIDATION ---")

path_obs <- "../../data/indexTurc_daily.rds"
path_pred_deep <- "../../data/pred_DeepESD_CV.rds"
path_pred_knn <- "../../data/pred_KNN_PCA_CV_FINAL.rds"

data_list_orig <- readRDS(file = path_obs)
data_list_pred_deep <- readRDS(file = path_pred_deep)
data_list_pred_knn <- readRDS(file = path_pred_knn)

obs_anual <- agregacion_anual_rapida(data_list_orig)
pred_deep_anual <- agregacion_anual_rapida(data_list_pred_deep)
pred_knn_anual <- agregacion_anual_rapida(data_list_pred_knn)

puntos_espana <- data.frame(
  Region = c("Galicia (Northwest)", "Catalonia (Northeast)", "Mediterranean (Levante)", 
             "Center (Plateau)", "South (Andalusia)", "North (Cantabria)"),
  Lon = c(-8.0, 2.0, -0.5, -3.7, -5.9, -4.0),
  Lat = c(42.5, 41.8, 39.5, 40.4, 37.4, 43.3),
  stringsAsFactors = FALSE
)

formatear_df_con_metricas <- function(obs_grid, pred_grid) {
  df_todas <- data.frame()
  df_metricas <- data.frame()
  
  for(i in 1:nrow(puntos_espana)) {
    df_punto <- extraer_serie_punto(obs_grid, pred_grid, puntos_espana$Lon[i], puntos_espana$Lat[i], puntos_espana$Region[i])
    df_todas <- rbind(df_todas, df_punto)
    
    bias <- mean(df_punto$Modelo - df_punto$Observado, na.rm = TRUE)
    corr <- cor(df_punto$Modelo, df_punto$Observado, use = "complete.obs")
    var_ratio <- var(df_punto$Modelo, na.rm = TRUE) / var(df_punto$Observado, na.rm = TRUE)
    label_text <- sprintf("Bias: %.2f\nCorr: %.2f\nVar Ratio: %.2f", bias, corr, var_ratio)
    
    df_metricas <- rbind(df_metricas, data.frame(Region = puntos_espana$Region[i], label = label_text, stringsAsFactors = FALSE))
  }
  
  df_obs <- data.frame(Fecha = df_todas$Fecha, Region = df_todas$Region, Tipo = "Observed", Indice_Turc = df_todas$Observado)
  df_mod <- data.frame(Fecha = df_todas$Fecha, Region = df_todas$Region, Tipo = "Model", Indice_Turc = df_todas$Modelo)
  df_long <- rbind(df_obs, df_mod)
  df_long <- df_long[!is.na(df_long$Indice_Turc), ]
  
  return(list(long = df_long, raw = df_todas, metricas = df_metricas))
}

datos_deep <- formatear_df_con_metricas(obs_anual, pred_deep_anual)
datos_knn <- formatear_df_con_metricas(obs_anual, pred_knn_anual)

tema_series <- theme_minimal(base_size = 20) +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    strip.text = element_text(face = "bold", size = 18, margin = margin(t = 10, b = 10)),
    strip.background = element_rect(fill = "gray90", color = NA),
    plot.title = element_text(face = "bold", size = 26), 
    axis.title = element_text(face = "bold", size = 22),
    axis.text = element_text(size = 18) 
  )

crear_grafico_anual <- function(datos, titulo, archivo_salida) {
  fecha_min <- min(datos$long$Fecha)
  p <- ggplot(datos$long, aes(x = Fecha, y = Indice_Turc, color = Tipo)) +
    geom_line(alpha = 0.7, linewidth = 1.2) +
    geom_point(size = 3) + 
    facet_wrap(~ Region, scales = "free_y", ncol = 2) +
    scale_color_manual(values = c("Observed" = "black", "Model" = "firebrick")) +
    geom_label(data = datos$metricas, aes(x = fecha_min, y = Inf, label = label), inherit.aes = FALSE, hjust = 0, vjust = 1.1, fill = "white", alpha = 0.85, label.size = 0.5, fontface = "bold", size = 6) +
    labs(title = titulo, x = "Year", y = "Accumulated Annual Turc Index", color = "Dataset") +
    tema_series
  ggsave(archivo_salida, plot = p, width = 18, height = 14) 
  return(p)
}

message("Generating Time Series PDF outputs...")
crear_grafico_anual(datos_deep, "Annual Turc Index Time Series - DeepESD", "series/Series_Turc_DeepESD.pdf")
crear_grafico_anual(datos_knn, "Annual Turc Index Time Series - KNN", "series/Series_Turc_KNN.pdf")


# ==============================================================================
# PART 3: GLOBAL CLIMATE DRIVERS ANALYSIS
# ==============================================================================
message("\n--- PART 3: GLOBAL CLIMATE DRIVERS (TURC vs CLIMATE VARIABLES) ---")

message("Calculating Daily Mean Temperature (Tmean)...")
tmean_daily <- gridArithmetics(data_list_pti$tmax, data_list_pti$tmin, operator = "+")
tmean_daily <- gridArithmetics(tmean_daily, 2, operator = "/")

message("Extracting global spatial daily means...")
fechas  <- as.Date(turc_daily_pti$Dates$start)
years   <- as.numeric(format(fechas, "%Y"))

ts_turc_daily  <- get_spatial_mean(turc_daily_pti$Data, mask_na)
ts_tmean_daily <- get_spatial_mean(tmean_daily$Data, mask_na)
ts_pr_daily    <- get_spatial_mean(data_list_pti$pr$Data, mask_na)
ts_ssrd_daily  <- get_spatial_mean(data_list_pti$ssrd$Data, mask_na)
ts_hr_daily    <- get_spatial_mean(data_list_pti$hr$Data, mask_na)

df_daily <- data.frame(
  Date = fechas, Year = years, Turc = ts_turc_daily,
  Tmean = ts_tmean_daily, Pr = ts_pr_daily, Rad = ts_ssrd_daily, Hr = ts_hr_daily
)

df_annual <- df_daily %>%
  group_by(Year) %>%
  summarise(
    Turc_Sum = sum(Turc, na.rm = TRUE),
    Tmean_Mean = mean(Tmean, na.rm = TRUE),
    Pr_Sum = sum(Pr, na.rm = TRUE),
    Rad_Mean = mean(Rad, na.rm = TRUE),
    Hr_Mean = mean(Hr, na.rm = TRUE)
  )

plot_turc_trend <- ggplot(df_annual, aes(x = Year, y = Turc_Sum)) +
  geom_line(color = "black", linewidth = 1.2) +
  geom_point(color = "black", size = 3) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 20), axis.text.x = element_text(angle = 45, hjust = 1), axis.title = element_text(face = "bold")) +
  scale_x_continuous(breaks = seq(min(df_annual$Year), max(df_annual$Year), by = 2)) +
  labs(title = "Annual Evolution of the Turc Index (1981 - 2020)", x = "Year", y = "Turc Index")

ggsave("Annual_Evolution_Turc.png", plot = plot_turc_trend, width = 12, height = 7, dpi = 300)

message("Calculating seasonal aggregations...")
df_daily_season <- df_daily %>%
  mutate(Month = as.numeric(format(Date, "%m")), Season = case_when(Month %in% c(12, 1, 2) ~ "Winter", Month %in% c(3, 4, 5) ~ "Spring", Month %in% c(6, 7, 8) ~ "Summer", Month %in% c(9, 10, 11) ~ "Autumn")) %>%
  mutate(Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn")))

df_season <- df_daily_season %>%
  group_by(Year, Season) %>%
  summarise(Turc_Sum = sum(Turc, na.rm = TRUE), Tmean_Mean = mean(Tmean, na.rm = TRUE), Pr_Sum = sum(Pr, na.rm = TRUE), .groups = "drop")

cor_tmean_season <- df_season %>% group_by(Season) %>% summarise(cor_val = cor(Tmean_Mean, Turc_Sum, use = "complete.obs")) %>% mutate(cor_label = sprintf("r = %.2f", cor_val))
plot_season_tmean <- ggplot(df_season, aes(x = Tmean_Mean, y = Turc_Sum)) +
  geom_point(color = "black", size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "red", fill = "gray70", se = TRUE) +
  geom_text(data = cor_tmean_season, aes(x = -Inf, y = Inf, label = cor_label), hjust = -0.2, vjust = 1.5, size = 6, fontface = "bold", color = "black") +
  facet_wrap(~ Season, scales = "free_x", ncol = 4) + 
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), strip.background = element_rect(fill = "#fdf6e3", color = "black"), strip.text = element_text(face = "bold", color = "black")) +
  labs(title = "Effect of Temperature on the Turc Index by Season", y = "Seasonal Turc Index", x = "Mean Seasonal Temperature (°C)")

cor_pr_season <- df_season %>% group_by(Season) %>% summarise(cor_val = cor(Pr_Sum, Turc_Sum, use = "complete.obs")) %>% mutate(cor_label = sprintf("r = %.2f", cor_val))
plot_season_pr <- ggplot(df_season, aes(x = Pr_Sum, y = Turc_Sum)) +
  geom_point(color = "black", size = 2.5, alpha = 0.8) +
  geom_smooth(method = "lm", color = "blue", fill = "gray70", se = TRUE) + 
  geom_text(data = cor_pr_season, aes(x = -Inf, y = Inf, label = cor_label), hjust = -0.2, vjust = 1.5, size = 6, fontface = "bold", color = "black") +
  facet_wrap(~ Season, scales = "free_x", ncol = 4) +
  theme_bw(base_size = 18) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5), strip.background = element_rect(fill = "#e0f3db", color = "black"), strip.text = element_text(face = "bold", color = "black")) +
  labs(title = "Effect of Precipitation on the Turc Index by Season", y = "Seasonal Turc Index", x = "Total Seasonal Precipitation (mm)")

message("Exporting Combined Seasonal Drivers Panel...")
png("Combined_Seasonal_Panel.png", width = 16, height = 12, units = "in", res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(2, 1)))
print(plot_season_tmean, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(plot_season_pr, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

message("\n>>> FULL WORKFLOW COMPLETED SUCCESSFULLY!")