# ==============================================================================
# GLOBAL CLIMATE MODEL (GCM) DOWNSCALING & GWL PROJECTIONS (kNN)
# 
# Description:
# This script applies a globally trained K-Nearest Neighbors (kNN) analogue 
# model with PCA to downscale a suite of Global Climate Models (GCMs). 
# It evaluates the historical period (1981-2005) against ERA5 and observations, 
# and computes future anomalies (deltas) at specific Global Warming Levels (GWL):
# +2.0ºC, +3.0ºC, and +4.0ºC based on the RCP8.5 scenario.
#
# The workflow includes:
#   1. Training a global kNN model on the full historical ERA5 dataset.
#   2. GCM downscaling utilizing Scaling Delta Mapping for bias correction.
#   3. Multi-model Historical Evaluation (Bias & Correlation mapping).
#   4. Future GWL extraction and multi-model Ensemble Mean calculations.
#   5. High-resolution multi-panel plotting (Historical baselines & Future Deltas).
# ==============================================================================

# ==============================================================================
# SECTION 0: INITIAL SETUP & LIBRARIES
# ==============================================================================
options(java.parameters = "-Xmx32g")

library(magrittr)
library(yaml)
library(downscaleR)
library(transformeR)
library(visualizeR)
library(RColorBrewer)
library(gridExtra)
library(grid)

# --- Model List & GWL Target Years Table ---
gcms <- c("ACCESS1-0", "CanESM2", "MPI-ESM-MR", "NorESM1-M")

gwl_data <- data.frame(
  GCM = gcms,
  GWL_2 = c(2041, 2026, 2039, 2048),
  GWL_3 = c(2062, 2049, 2060, 2072),
  GWL_4 = c(2081, 2068, 2082, NA) # Note: NorESM1-M does not reach +4.0ºC
)

# --- General Parameters & Hyperparameters ---
n_vecinos <- 1        # Number of analogues (k)
varianza_pca <- 0.90  # Retained variance for PCA (90%)
years_target <- 1981:2005

# --- Visualization Parameters ---
limite_bias  <- 2
limite_delta <- 5

cex_main   <- 2.2
cex_na     <- 20   
cex_top    <- 26   
cex_legend <- 1.6 
width_legend <- 1.8 

my_colorkey <- list(labels = list(cex = cex_legend), width = width_legend)

# ==============================================================================
# SECTION 1: HELPER FUNCTIONS (HARMONIZATION & VISUALIZATION)
# ==============================================================================

# Scaling Delta Mapping for bias correction
scalingDeltaMapping <- function(grid, base, ref) {
  grid_detrended <- scaleGrid(grid, base = grid, ref = base, type = "center", spatial.frame = "gridbox", time.frame = "monthly", skip.season.check = TRUE)  
  grid_detrended_corrected <- scaleGrid(grid_detrended, base = base, ref = ref, type = "standardize", spatial.frame = "gridbox", time.frame = "monthly", skip.season.check = TRUE)    
  grid_corrected <- scaleGrid(grid_detrended_corrected, base = base, ref = grid, type = "center", spatial.frame = "gridbox", time.frame = "monthly", skip.season.check = TRUE)    
  return(grid_corrected)
}

# Standard scaling
standardize <- function(grid, base){
  scaleGrid(grid, base = base, ref = NULL, type = "standardize", skip.season.check = TRUE) 
} 

# Fill spatial NAs via nearest neighbor (crucial for GCM topologies)
fill_all_nas_multigrid <- function(multigrid) {
  multigrid <- transformeR::redim(multigrid, drop = TRUE)
  data_array <- multigrid$Data
  n_vars <- dim(data_array)[1]
  n_lat <- dim(data_array)[3]
  n_lon <- dim(data_array)[4]
  lat_matrix <- matrix(rep(1:n_lat, n_lon), nrow = n_lat, ncol = n_lon)
  lon_matrix <- matrix(rep(1:n_lon, each = n_lat), nrow = n_lat, ncol = n_lon)
  
  for (v in 1:n_vars) {
    capa_ref <- data_array[v, 1, , ]
    nas <- is.na(capa_ref)
    if (any(nas) && !all(nas)) {
      valid_idx <- which(!nas)
      na_idx <- which(nas)
      nearest_valid_idx <- integer(length(na_idx))
      for (k in seq_along(na_idx)) {
        i <- na_idx[k]
        distancias <- (lat_matrix[i] - lat_matrix[valid_idx])^2 + (lon_matrix[i] - lon_matrix[valid_idx])^2
        nearest_valid_idx[k] <- valid_idx[which.min(distancias)]
      }
      for (t in 1:dim(data_array)[2]) {
        capa_dia <- data_array[v, t, , ]
        capa_dia[na_idx] <- capa_dia[nearest_valid_idx]
        data_array[v, t, , ] <- capa_dia
      }
    }
  }
  multigrid$Data <- data_array
  return(multigrid)
}

# Climatology (Handles messy sub-daily aggregation in GCMs)
calcular_climatologia <- function(grid, mask_oceanos, nombre_dataset = "Dataset") {
  message(sprintf("\n--- Processing climatology: %s ---", nombre_dataset))
  raw_dates <- as.character(unlist(grid$Dates))
  d_str <- grep("^[0-9]{4}-[0-9]{2}-[0-9]{2}", raw_dates, value = TRUE)
  n_time <- dim(grid$Data)[1]
  if (length(d_str) > n_time) d_str <- d_str[1:n_time]
  
  unique_days <- unique(d_str)
  if (length(unique_days) < length(d_str)) {
    message(sprintf("  [!] Detected %d records for %d unique days. Collapsing (MEAN)...", length(d_str), length(unique_days)))
    daily_data <- array(NA, dim = c(length(unique_days), dim(grid$Data)[2], dim(grid$Data)[3]))
    for (d in seq_len(length(unique_days))) {
      idx_day <- which(d_str == unique_days[d])
      if (length(idx_day) > 1) {
        daily_data[d, , ] <- colMeans(grid$Data[idx_day, , , drop = FALSE], na.rm = TRUE, dims = 1)
      } else {
        daily_data[d, , ] <- grid$Data[idx_day, , ]
      }
    }
    grid$Data <- daily_data
    d_str <- unique_days
  }
  
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
  
  clim_data <- colMeans(grid_annual$Data, na.rm = TRUE, dims = 1)
  grid_clim <- grid_annual
  grid_clim$Data <- array(clim_data, dim = c(1, dim(grid_annual$Data)[2], dim(grid_annual$Data)[3]))
  grid_clim$Dates$start <- grid_annual$Dates$start[1]
  grid_clim$Dates$end <- grid_annual$Dates$end[length(grid_annual$Dates$end)]
  
  capa_clim <- grid_clim$Data[1, , ]
  capa_clim[mask_oceanos] <- NA
  grid_clim$Data[1, , ] <- capa_clim
  attr(grid_clim$Data, "dimensions") <- c("time", "lat", "lon")
  if(is.null(grid_clim$Variable$varName)) grid_clim$Variable$varName <- "Turc"
  
  return(list(clim = grid_clim, annual = grid_annual))
}

# Forces sub-daily dates into strict YYYY-MM-DD formats for consistency
limpiar_subdiario <- function(grid) {
  raw_dates <- as.character(unlist(grid$Dates))
  d_str <- grep("^[0-9]{4}-[0-9]{2}-[0-9]{2}", raw_dates, value = TRUE)
  n_time <- dim(grid$Data)[1]
  if (length(d_str) > n_time) d_str <- d_str[1:n_time]
  unique_days <- unique(d_str)
  
  if (length(unique_days) < length(d_str)) {
    daily_data <- array(NA, dim = c(length(unique_days), dim(grid$Data)[2], dim(grid$Data)[3]))
    for (d in seq_len(length(unique_days))) {
      idx_day <- which(d_str == unique_days[d])
      if (length(idx_day) > 1) {
        daily_data[d, , ] <- colMeans(grid$Data[idx_day, , , drop = FALSE], na.rm = TRUE, dims = 1)
      } else {
        daily_data[d, , ] <- grid$Data[idx_day, , ]
      }
    }
    grid$Data <- daily_data
    attr(grid$Data, "dimensions") <- c("time", "lat", "lon")
    grid$Dates$start <- paste0(unique_days, " 00:00:00 GMT")
    grid$Dates$end   <- paste0(unique_days, " 24:00:00 GMT")
  }
  return(grid)
}

# Interannual Correlation Calculation with Stippling
calc_cor_interanual <- function(obs_annual, pred_annual, mask_oceanos, threshold = 0.05) {
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

# ==============================================================================
# SECTION 2: DATA LOADING (ERA5 PREDICTORS & PREDICTAND)
# ==============================================================================
message("\n>>> Loading observations and ERA5 predictors...")
path_data <- "../../data/ERA5/"
variables_base <- c("u", "v", "q", "zg", "t")
levels <- c("850", "500")
variables_list <- as.vector(outer(variables_base, levels, paste0))

grid_list <- lapply(variables_list, function(var) {
  file_path <- paste0(path_data, var, "_ERA5_daily_1981_2020.rds")
  if (file.exists(file_path)) return(readRDS(file_path)) else return(NULL)
})

x_ERA5 <- do.call(makeMultiGrid, grid_list) %>% redim(drop = TRUE)

y_obs = readRDS(file = '../../data/indexTurc_daily.rds')
attr(y_obs$Data, "dimensions") <- c("time", "lat", "lon")
y_obs <- filterNA(y_obs)
y_obs[['Data']][is.na(y_obs[['Data']])] <- -999

xT = intersectGrid(x_ERA5, y_obs, which.return = 1)
yT = intersectGrid(x_ERA5, y_obs, which.return = 2)

rm(x_ERA5, grid_list)
gc()

# ==============================================================================
# SECTION 3: CROSS-VALIDATION (ERA5 BASELINE)
# ==============================================================================
dir.create("predicciones_tmp", showWarnings = FALSE)
years_all = unique(getYearsAsINDEX(xT)) 
k = 4
folds = split(years_all, cut(seq_along(years_all), breaks = k, labels = FALSE))

vars_to_keep = c(ls(), "vars_to_keep", "i", "folds", "k")

for (i in 1:k) {
    years_test = folds[[i]]
    years_train = setdiff(years_all, years_test)
    
    f_out = file.path("predicciones_tmp", paste0("prediccion_KNN_PCA_fold_", i, ".rds"))
    if (file.exists(f_out)) next
    
    message("--- Starting Fold ", i, " (KNN + PCA 90%) ---")
    xT_train = subsetGrid(xT, years = years_train)
    yT_train = subsetGrid(yT, years = years_train)
    xT_test = subsetGrid(xT, years = years_test)

    xT_train_scaled = scaleGrid(xT_train, type = "standardize", spatial.frame = "gridbox") %>% redim(drop = TRUE) %>% suppressMessages()
    xT_train_scaled[['Data']][is.na(xT_train_scaled[['Data']])] = 0

    xT_test_scaled = scaleGrid(xT_test, base = xT_train, type = "standardize", spatial.frame = "gridbox") %>% redim(drop = TRUE) %>% suppressMessages()
    xT_test_scaled[['Data']][is.na(xT_test_scaled[['Data']])] = 0
    
    yT_train_clean <- yT_train %>% redim(drop = TRUE)
    
    data_train <- prepareData(x = xT_train_scaled, y = yT_train_clean, spatial.predictors = list(v.exp = varianza_pca))
    model <- downscaleTrain(data_train, method = "analogs", n.analogs = n_vecinos)
    
    newdata <- prepareNewData(xT_test_scaled, data_train)
    pred <- downscalePredict(newdata, model)

    yT_test_template <- subsetGrid(yT, years = years_test)
    pred[['Data']][yT_test_template[['Data']] < 0] <- NA
    
    saveRDS(pred, file = f_out)
    rm(list = setdiff(ls(), vars_to_keep))
    gc()
}

if (!file.exists('../../data/pred_KNN_PCA_CV_FINAL.rds')) {
  message("Merging Cross-Validation predictions...")
  files = list.files("predicciones_tmp", pattern = "KNN_PCA", full.names = TRUE)
  folds_in_files = as.numeric(gsub(".*_fold_(\\d+)\\.rds", "\\1", basename(files)))
  prediction_final = bindGrid(lapply(files[order(folds_in_files)], readRDS), dimension = "time")
  saveRDS(prediction_final, '../../data/pred_KNN_PCA_CV_FINAL.rds')
}

# ==============================================================================
# SECTION 4: GLOBAL MODEL TRAINING
# ==============================================================================
message("\n--- Training Global Model (1981-2020) for GCM projections ---")
xT_scaled_full = scaleGrid(xT, type = "standardize", spatial.frame = "gridbox") %>% redim(drop = TRUE) %>% suppressMessages()
xT_scaled_full[['Data']][is.na(xT_scaled_full[['Data']])] = 0

data_train_global <- prepareData(x = xT_scaled_full, y = yT %>% redim(drop = TRUE), spatial.predictors = list(v.exp = varianza_pca))
model_global <- downscaleTrain(data_train_global, method = "analogs", n.analogs = n_vecinos)
rm(xT_scaled_full); gc()

# ==============================================================================
# SECTION 5: GCM DOWNSCALING LOOP
# ==============================================================================
dir.create("../../data/proyecciones_kNN", showWarnings = FALSE, recursive = TRUE)

for (model_name in gcms) {
  message(sprintf("\n>>> Processing Downscaling for GCM: %s ...", model_name))
  file_path <- paste0("../../data/GCMs/GCM_data_", model_name, ".rds")
  
  if(!file.exists(file_path)) {
    message("   [!] File not found for ", model_name, ". Skipping prediction...")
    next
  }
  
  gcm_data <- readRDS(file_path)
  base_vars_gcm <- c("ua", "va", "hus", "zg", "ta") 
  levels_gcm <- c("85000", "50000")
  var_list_gcm_target <- as.vector(outer(base_vars_gcm, levels_gcm, paste, sep="_"))
  var_list_era5_target <- as.vector(outer(c("u", "v", "q", "zg", "t"), c("850", "500"), paste0))
  
  harmonize.base <- NULL 
  
  for (period_name in c("historical", "rcp85")) {
    output.file <- sprintf("../../data/proyecciones_kNN/Turc_KNN_PCA_%s_%s.rds", model_name, period_name)
    
    if (!file.exists(output.file)) {
      grid_list_gcm <- list()
      for (i in seq_along(var_list_gcm_target)) {
        k <- var_list_gcm_target[i]; new_name <- var_list_era5_target[i]
        grid <- gcm_data[[period_name]][[k]]
        grid$Data[grid$Data > 1e10] <- NA; grid$Data[grid$Data < -1e10] <- NA
        grid$Variable$varName <- new_name
        grid <- interpGrid(grid, new.coordinates = getGrid(xT), method = "bilinear")
        grid_list_gcm[[new_name]] <- grid
      }
      
      # Strict temporal subsetting
      common_dates <- Reduce(intersect, lapply(grid_list_gcm, function(g) substr(as.character(g$Dates$start), 1, 10)))
      grid_list_gcm_harmonized <- lapply(names(grid_list_gcm), function(name) {
        g <- grid_list_gcm[[name]]
        idx <- match(common_dates, substr(as.character(g$Dates$start), 1, 10))
        transformeR::subsetDimension(g, dimension = "time", indices = idx[!is.na(idx)])
      })
      names(grid_list_gcm_harmonized) <- names(grid_list_gcm)
      
      args_multigrid <- unname(grid_list_gcm_harmonized)
      args_multigrid$skip.temporal.check <- TRUE
      x_raw <- do.call(makeMultiGrid, args_multigrid) %>% redim(drop = TRUE)
      x_raw$Variable <- xT$Variable
      x_raw <- fill_all_nas_multigrid(x_raw)
      
      # Extract bias correction reference period from the historical scenario
      if (period_name == "historical") {
        ref.period <- intersect(unique(getYearsAsINDEX(xT)), unique(getYearsAsINDEX(x_raw))) 
        harmonize.base <- subsetGrid(x_raw, years = ref.period)
      }
      
      x_harm <- scalingDeltaMapping(grid = x_raw, base = harmonize.base, ref = subsetGrid(xT, years = ref.period) %>% suppressMessages())
      xn <- standardize(grid = x_harm, base = xT)
      xn[['Data']][is.na(xn[['Data']])] <- 0
      
      newdata_gcm <- prepareNewData(newdata = xn, data.structure = data_train_global)
      pred_grid <- downscalePredict(newdata = newdata_gcm, model = model_global)
      
      saveRDS(pred_grid, file = output.file, compress = "xz")
      rm(grid_list_gcm, x_raw, x_harm, xn, newdata_gcm, pred_grid); gc()
    } else {
      message("   -> Prediction ", period_name, " already exists for ", model_name)
    }
  }
  rm(gcm_data); gc()
}

# ==============================================================================
# SECTION 6: HISTORICAL MULTI-MODEL ANALYSIS & ENSEMBLES
# ==============================================================================
message("\n=========================================")
message(" INITIATING MULTI-MODEL HISTORICAL BLOCK ")
message("=========================================\n")

data_obs_raw <- readRDS('../../data/indexTurc_daily.rds')
attr(data_obs_raw$Data, "dimensions") <- c("time", "lat", "lon")
data_era5_raw <- readRDS('../../data/pred_KNN_PCA_CV_FINAL.rds')

data_obs_base <- subsetGrid(data_obs_raw, years = years_target)
data_era5_base <- subsetGrid(data_era5_raw, years = years_target)

data_obs_sync <- intersectGrid(data_era5_base, data_obs_base, type = "temporal", which.return = 2)
data_era5_sync <- intersectGrid(data_era5_base, data_obs_base, type = "temporal", which.return = 1)

suma_total <- colSums(data_era5_sync$Data, na.rm = TRUE, dims = 1)
mask_oceanos <- is.na(suma_total) | suma_total == 0

res_obs_era5 <- calcular_climatologia(data_obs_sync, mask_oceanos, "Obs (Sync with ERA5)")
res_era5 <- calcular_climatologia(data_era5_sync, mask_oceanos, "ERA5")

clim_bias_era5 <- res_era5$clim
clim_bias_era5$Data <- res_era5$clim$Data - res_obs_era5$clim$Data
attr(clim_bias_era5$Data, "dimensions") <- c("time", "lat", "lon")
cor_era5 <- calc_cor_interanual(res_obs_era5$annual, res_era5$annual, mask_oceanos)

plot_list_hist <- list()
plot_list_hist[[1]] <- spatialPlot(res_era5$clim, backdrop.theme="countries", at=seq(0,40,length.out=21), main=list(label="ERA5", cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
plot_list_hist[[2]] <- spatialPlot(clim_bias_era5, backdrop.theme="countries", at=seq(-limite_bias, limite_bias, length.out=21), set.min=-limite_bias, set.max=limite_bias, main=list(label="ERA5 Bias", cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
plot_list_hist[[3]] <- spatialPlot(cor_era5$cor_grid, sp.layout=list(cor_era5$pts), backdrop.theme="countries", at=seq(-1, 1, length.out=21), main=list(label="ERA5 Cor", cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)

hist_clim_dict <- list()
hist_bias_list <- list()
hist_cor_list  <- list()

for (gcm in gcms) {
  ruta_hist <- paste0('../../data/proyecciones_kNN/Turc_KNN_PCA_', gcm, '_historical.rds')
  if(!file.exists(ruta_hist)) next
  
  data_hist_raw <- readRDS(ruta_hist) %>% limpiar_subdiario()
  safe_years <- intersect(years_target, unique(as.numeric(substr(data_hist_raw$Dates$start, 1, 4))))
  
  if(length(safe_years) == 0) {
    message(paste("  [!] No common years for", gcm))
    next
  }
  
  data_hist_sub <- subsetGrid(data_hist_raw, years = safe_years)
  data_obs_gcm <- intersectGrid(data_hist_sub, data_obs_base, type="temporal", which.return=2)
  data_hist_gcm <- intersectGrid(data_hist_sub, data_obs_base, type="temporal", which.return=1)
  
  res_obs_gcm <- calcular_climatologia(data_obs_gcm, mask_oceanos, paste("Obs_sync", gcm))
  res_hist_gcm <- calcular_climatologia(data_hist_gcm, mask_oceanos, gcm)
  
  clim_bias_gcm <- res_hist_gcm$clim
  clim_bias_gcm$Data <- res_hist_gcm$clim$Data - res_obs_gcm$clim$Data
  attr(clim_bias_gcm$Data, "dimensions") <- c("time", "lat", "lon")
  
  cor_gcm <- calc_cor_interanual(res_obs_gcm$annual, res_hist_gcm$annual, mask_oceanos)
  
  hist_clim_dict[[gcm]] <- res_hist_gcm$clim
  hist_bias_list[[gcm]] <- clim_bias_gcm$Data
  hist_cor_list[[gcm]]  <- cor_gcm$cor_grid$Data
  
  p_pred <- spatialPlot(res_hist_gcm$clim, backdrop.theme="countries", at=seq(0,40,length.out=21), main=list(label=paste(gcm), cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
  p_bias <- spatialPlot(clim_bias_gcm, backdrop.theme="countries", at=seq(-limite_bias, limite_bias, length.out=21), set.min=-limite_bias, set.max=limite_bias, main=list(label=paste(gcm, "Bias"), cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
  p_cor <- spatialPlot(cor_gcm$cor_grid, sp.layout=list(cor_gcm$pts), backdrop.theme="countries", at=seq(-1, 1, length.out=21), main=list(label=paste(gcm, "Cor"), cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
  
  plot_list_hist <- c(plot_list_hist, list(p_pred, p_bias, p_cor))
}

# --- Historical Ensemble ---
message("  -> Calculating Historical Ensemble Mean...")
ens_clim_data <- Reduce("+", lapply(hist_clim_dict, function(x) x$Data)) / length(hist_clim_dict)
ens_bias_data <- Reduce("+", hist_bias_list) / length(hist_bias_list)
ens_cor_data  <- Reduce("+", hist_cor_list) / length(hist_cor_list)

ens_template <- hist_clim_dict[[1]]
ens_clim_grid <- ens_template; ens_clim_grid$Data <- ens_clim_data
ens_bias_grid <- ens_template; ens_bias_grid$Data <- ens_bias_data
ens_cor_grid  <- ens_template; ens_cor_grid$Data  <- ens_cor_data

p_ens_pred <- spatialPlot(ens_clim_grid, backdrop.theme="countries", at=seq(0,40,length.out=21), main=list(label="Ensemble Mean", cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
p_ens_bias <- spatialPlot(ens_bias_grid, backdrop.theme="countries", at=seq(-limite_bias, limite_bias, length.out=21), set.min=-limite_bias, set.max=limite_bias, main=list(label="Ensemble Mean Bias", cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)
p_ens_cor  <- spatialPlot(ens_cor_grid, backdrop.theme="countries", at=seq(-1, 1, length.out=21), main=list(label="Ensemble Mean Cor", cex=cex_main), color.theme="RdYlBu", rev.colors=TRUE, colorkey=my_colorkey)

plot_list_hist <- c(plot_list_hist, list(p_ens_pred, p_ens_bias, p_ens_cor))

message(">>> Exporting Historical Comparison PDF...")
pdf(file = "Comparativa_Turc_KNN_Hist.pdf", width = 16, height = 25)
grid.arrange(grobs = plot_list_hist, ncol = 3, nrow = length(gcms) + 2, 
             top = textGrob("Historical Comparison (1981-2005) of Turc Index: ERA5 & KNN Downscaled GCMs", gp = gpar(fontsize = cex_top, font = 2)))
dev.off()

message(">>> Exporting Historical Comparison PNG...")
png(file = "Comparativa_Turc_KNN_Hist.png", width = 16, height = 22, units = "in", res = 300)
grid.arrange(grobs = plot_list_hist, ncol = 3, nrow = length(gcms) + 2, 
             top = textGrob("Historical Comparison (1981-2005) of Turc Index: ERA5 & KNN Downscaled GCMs", gp = gpar(fontsize = cex_top, font = 2)))
dev.off()

# ==============================================================================
# SECTION 7: FUTURE GWL PROJECTIONS & ENSEMBLES
# ==============================================================================
message("\n=========================================")
message(" INITIATING GWL PROJECTIONS BLOCK        ")
message("=========================================\n")

plot_list_proj <- list()
niveles_gwl <- c("GWL_2", "GWL_3", "GWL_4")
grados <- c("2ºC", "3ºC", "4ºC")
deltas_gwl_list <- list(GWL_2 = list(), GWL_3 = list(), GWL_4 = list())

for (i in seq_along(gcms)) {
  gcm <- gcms[i]
  ruta_rcp <- paste0('../../data/proyecciones_kNN/Turc_KNN_PCA_', gcm, '_rcp85.rds')
  
  if (file.exists(ruta_rcp) && !is.null(hist_clim_dict[[gcm]])) {
    data_rcp_raw <- readRDS(ruta_rcp) %>% limpiar_subdiario()
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
        
        deltas_gwl_list[[lvl]][[gcm]] <- delta_gwl$Data
        
        titulo_mapa <- paste(gcm, "- Delta", grados[j], "\n(", gwl_years[1], "-", gwl_years[20], ")")
        
        p_delta <- spatialPlot(delta_gwl, backdrop.theme = "countries",
                               at = seq(-limite_delta, limite_delta, length.out = 21), 
                               set.min = -limite_delta, set.max = limite_delta,
                               main = list(label=titulo_mapa, cex=cex_main),
                               color.theme = "RdYlBu", rev.colors = TRUE,
                               colorkey=my_colorkey)
        
        plot_list_proj <- c(plot_list_proj, list(p_delta))
      }
    }
    rm(data_rcp_raw); gc()
  } else {
    message(paste("File not found or historical baseline failed for:", gcm))
  }
}

# --- Projections Ensemble ---
message("  -> Calculating Projections Ensemble Mean...")
for (j in seq_along(niveles_gwl)) {
  lvl <- niveles_gwl[j]
  datos_nivel <- deltas_gwl_list[[lvl]]
  
  if (length(datos_nivel) > 0) {
    ens_delta_data <- Reduce("+", datos_nivel) / length(datos_nivel)
    ens_delta_grid <- hist_clim_dict[[1]]
    ens_delta_grid$Data <- ens_delta_data
    
    titulo_ens <- paste("Ensemble Mean - Delta", grados[j], "\n(", length(datos_nivel), "models)")
    
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

message(">>> Exporting GWL Projections PDF...")
pdf(file = "Proyecciones_Deltas_KNN.pdf", width = 16, height = 22)
grid.arrange(grobs = plot_list_proj, ncol = 3, nrow = length(gcms) + 1, 
             top = textGrob("Turc Index Anomalies (KNN Downscaled) by GWL in RCP8.5 vs Historical", gp = gpar(fontsize = cex_top, font = 2)))
dev.off()

message(">>> Exporting GWL Projections PNG in high resolution...")
png(file = "Proyecciones_Deltas_KNN.png", width = 16, height = 22, units = "in", res = 300)
grid.arrange(grobs = plot_list_proj, ncol = 3, nrow = length(gcms) + 1, 
             top = textGrob("Turc Index Anomalies (KNN Downscaled) by GWL in RCP8.5 vs Historical", gp = gpar(fontsize = cex_top, font = 2)))
dev.off()

message("\n>>> PROCESS COMPLETED!")