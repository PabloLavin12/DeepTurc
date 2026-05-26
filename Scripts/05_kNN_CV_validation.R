# ==============================================================================
# K-NEAREST NEIGHBORS (KNN) BASELINE MODEL & VALIDATION
# 
# Description:
# This script establishes a baseline Empirical Statistical Downscaling (ESD) 
# model using the k-Nearest Neighbors (kNN) analogue method, coupled with 
# Principal Component Analysis (PCA). It evaluates the KNN model using the 
# same rigorous cross-validation (4-Fold) and metrics applied to the DeepESD 
# model, ensuring a fair 1:1 comparison.
#
# The workflow includes:
#   1. Helper functions for metrics (Bias, Correlation, Intra-annual distribution).
#   2. Data loading and temporal intersection (climate4R framework).
#   3. 4-Fold cross-validation using Analogs (kNN) with 90% PCA variance.
#   4. Generation of high-resolution spatial maps and boxplots for KNN.
#   5. Visual Comparison (2x2 facet grids) of Spatial Correlation: DeepESD vs KNN.
# ==============================================================================

# ==============================================================================
# SECTION 0: INITIAL SETUP & LIBRARY LOADING
# ==============================================================================
options(java.parameters = "-Xmx32g")

library(magrittr)
library(yaml)
library(downscaleR)
library(transformeR)
library(visualizeR)
library(RColorBrewer)
library(ggplot2)
library(gridExtra)
library(grid)

# Create output directories if they do not exist
dir.create("mapas_metricas_knn", showWarnings = FALSE)
dir.create("boxplot_knn", showWarnings = FALSE)
dir.create("series_knn", showWarnings = FALSE)

# --- General Parameters & Hyperparameters ---
n_vecinos <- 1        # Number of analogues (k)
varianza_pca <- 0.90  # Retained variance for PCA (90%)

# ==============================================================================
# SECTION 1: HELPER FUNCTIONS FOR VALIDATION
# ==============================================================================

# Efficiently aggregates daily grid data into monthly totals.
agregacion_mensual_rapida <- function(grid_diario) {
  fechas <- as.Date(grid_diario$Dates$start)
  ym <- format(fechas, "%Y-%m")
  meses_unicos <- unique(ym)
  n_meses <- length(meses_unicos)
  n_lat <- dim(grid_diario$Data)[2]; n_lon <- dim(grid_diario$Data)[3]
  
  data_mensual <- array(NA, dim = c(n_meses, n_lat, n_lon))
  for (m in 1:n_meses) {
    idx <- which(ym == meses_unicos[m])
    data_mensual[m, , ] <- colSums(grid_diario$Data[idx, , , drop = FALSE], na.rm = TRUE, dims = 1)
  }
  
  grid_mensual <- grid_diario
  grid_mensual$Data <- data_mensual
  attr(grid_mensual$Data, "dimensions") <- c("time", "lat", "lon")
  grid_mensual$Dates$start <- as.character(tapply(grid_diario$Dates$start, ym, function(x) x[1]))
  grid_mensual$Dates$end <- as.character(tapply(grid_diario$Dates$end, ym, function(x) x[length(x)]))
  return(grid_mensual)
}

# Applies a logical mask to replace ocean pixels with NA across 2D or 3D grids.
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
  }
  return(grid)
}

# Calculates Pearson correlation and p-values pixel-by-pixel, applying stippling.
calc_cor_interanual <- function(obs_annual, pred_annual, threshold = 0.05) {
  lat_n <- dim(obs_annual$Data)[2]; lon_n <- dim(obs_annual$Data)[3]
  cor_array <- matrix(NA, nrow = lat_n, ncol = lon_n)
  pval_array <- matrix(NA, nrow = lat_n, ncol = lon_n)
  
  for (i in 1:lat_n) {
    for (j in 1:lon_n) {
      obs_series <- obs_annual$Data[, i, j]
      pred_series <- pred_annual$Data[, i, j]
      if (all(is.na(obs_series))) next
      
      valid_idx <- complete.cases(obs_series, pred_series)
      if (sum(valid_idx) >= 10) {
        suppressWarnings({
          test <- try(cor.test(pred_series[valid_idx], obs_series[valid_idx], method = "pearson"), silent = TRUE)
          if (!inherits(test, "try-error")) {
            cor_array[i, j] <- test$estimate
            pval_array[i, j] <- test$p.value
          }
        })
      }
    }
  }
  cor_grid <- list(Data = cor_array, xyCoords = obs_annual$xyCoords, Variable = list(varName="Correlation"), Dates = obs_annual$Dates)
  attr(cor_grid$Data, "dimensions") <- c("lat", "lon"); class(cor_grid) <- "grid"
  
  pval_grid <- list(Data = pval_array, xyCoords = obs_annual$xyCoords, Variable = list(varName="p-values"), Dates = obs_annual$Dates)
  attr(pval_grid$Data, "dimensions") <- c("lat", "lon"); class(pval_grid) <- "grid"
  
  pts <- map.stippling(climatology(pval_grid), threshold = threshold, condition = "LT", pch = 19, col = adjustcolor("black", alpha.f = 0.2), cex = 0.01) %>% suppressMessages() %>% suppressWarnings()
  return(list(cor_grid = climatology(cor_grid), pts = pts))
}

# Calculates global spatial bias and spatial correlation for each year.
calc_metricas_espaciales <- function(obs_annual, pred_annual) {
  fechas <- as.integer(format(as.Date(obs_annual$Dates$start), "%Y"))
  n_years <- length(fechas)
  
  df_spat <- data.frame(Year = fechas, Spatial_Bias = numeric(n_years), Spatial_Corr = numeric(n_years))
  df_bias_pixels <- data.frame() 
  
  for (i in 1:n_years) {
    obs_map <- as.vector(obs_annual$Data[i,,])
    pred_map <- as.vector(pred_annual$Data[i,,])
    
    valid_idx <- !is.na(obs_map) & !is.na(pred_map)
    if (sum(valid_idx) > 10) {
      df_spat$Spatial_Bias[i] <- mean(pred_map[valid_idx] - obs_map[valid_idx])
      df_spat$Spatial_Corr[i] <- cor(obs_map[valid_idx], pred_map[valid_idx])
      
      tmp_df <- data.frame(Year = fechas[i], Bias = pred_map[valid_idx] - obs_map[valid_idx])
      df_bias_pixels <- rbind(df_bias_pixels, tmp_df)
    } else {
      df_spat[i, 2:3] <- NA
    }
  }
  return(list(spat = df_spat, bias_pixels = df_bias_pixels))
}

# INTRA-ANNUAL Daily Metrics (Correlation Only)
calc_metricas_diarias_intra <- function(obs_daily, pred_daily, mask_na, s_months, s_name) {
  fechas_todas <- as.Date(obs_daily$Dates$start)
  meses_todos <- as.numeric(format(fechas_todas, "%m"))
  
  idx_estacion <- which(meses_todos %in% s_months)
  fechas_filtradas <- fechas_todas[idx_estacion]
  obs_data <- obs_daily$Data[idx_estacion, , , drop = FALSE]
  pred_data <- pred_daily$Data[idx_estacion, , , drop = FALSE]
  
  years_filtrados <- as.numeric(format(fechas_filtradas, "%Y"))
  if (s_name == "DJF") {
    # Shift December into the following year for proper winter grouping
    meses_filtrados <- as.numeric(format(fechas_filtradas, "%m"))
    idx_dic <- which(meses_filtrados == 12)
    years_filtrados[idx_dic] <- years_filtrados[idx_dic] + 1
  }
  
  unique_years <- unique(years_filtrados)
  df_intra <- data.frame()
  
  for (y in unique_years) {
    idx_y <- which(years_filtrados == y)
    if (length(idx_y) < 10) next 
    
    obs_y <- obs_data[idx_y, , , drop = FALSE]
    pred_y <- pred_data[idx_y, , , drop = FALSE]
    
    dim(obs_y) <- c(dim(obs_y)[1], dim(obs_y)[2] * dim(obs_y)[3])
    dim(pred_y) <- c(dim(pred_y)[1], dim(pred_y)[2] * dim(pred_y)[3])
    
    valid_pixels <- which(!as.vector(mask_na))
    
    if (length(valid_pixels) > 0) {
      obs_land <- obs_y[, valid_pixels, drop = FALSE]
      pred_land <- pred_y[, valid_pixels, drop = FALSE]
      
      cor_vals <- sapply(1:ncol(obs_land), function(k) {
        o_vec <- obs_land[, k]
        p_vec <- pred_land[, k]
        valid_idx <- complete.cases(o_vec, p_vec)
        if (sum(valid_idx) >= 10) {
          res <- try(cor(o_vec[valid_idx], p_vec[valid_idx], method = "pearson"), silent = TRUE)
          if (!inherits(res, "try-error")) return(res)
        }
        return(NA) 
      })
      tmp_df <- data.frame(Year = as.integer(y), Intra_Corr = cor_vals)
      df_intra <- rbind(df_intra, na.omit(tmp_df))
    }
  }
  return(df_intra)
}

# OVERLEAF GIANT FONT THEME
tema_overleaf <- theme_minimal(base_size = 24) +
  theme(plot.title = element_text(face = "bold", size = 28, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 22, margin = margin(b = 15)),
        axis.title = element_text(face = "bold", size = 24),
        axis.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", size = 24),
        legend.position = "bottom")

lattice_margins <- list(layout.heights = list(top.padding = 0, bottom.padding = 0, main.key.padding = 0),
                        layout.widths = list(left.padding = 0, right.padding = 0))

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

# Convert NAs to -999 to allow the modeling framework to process the grid
y_obs[['Data']][is.na(y_obs[['Data']])] <- -999

xT = intersectGrid(x_ERA5, y_obs, which.return = 1)
yT = intersectGrid(x_ERA5, y_obs, which.return = 2)

rm(x_ERA5, grid_list)
gc()

# ==============================================================================
# SECTION 3: K-FOLD CROSS-VALIDATION (ERA5)
# ==============================================================================
dir.create("predicciones_tmp", showWarnings = FALSE)
years_all = unique(getYearsAsINDEX(xT)) 
k = 4
folds = split(years_all, cut(seq_along(years_all), breaks = k, labels = FALSE))

# Variables to preserve during garbage collection to prevent memory leaks
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
    
    # Train KNN Analogs model with PCA reduction
    data_train <- prepareData(x = xT_train_scaled, y = yT_train_clean, spatial.predictors = list(v.exp = varianza_pca))
    model <- downscaleTrain(data_train, method = "analogs", n.analogs = n_vecinos)
    
    # Downscale predictions
    newdata <- prepareNewData(xT_test_scaled, data_train)
    pred <- downscalePredict(newdata, model)

    # Post-process: Mask out oceans based on the real test template
    yT_test_template <- subsetGrid(yT, years = years_test)
    pred[['Data']][yT_test_template[['Data']] < 0] <- NA
    
    saveRDS(pred, file = f_out)
    rm(list = setdiff(ls(), vars_to_keep))
    gc()
}

# Bind cross-validated folds into a single seamless historical prediction grid
if (!file.exists('../../data/pred_KNN_PCA_CV_FINAL.rds')) {
  message("Merging Cross-Validation predictions...")
  files = list.files("predicciones_tmp", pattern = "KNN_PCA", full.names = TRUE)
  folds_in_files = as.numeric(gsub(".*_fold_(\\d+)\\.rds", "\\1", basename(files)))
  prediction_final = bindGrid(lapply(files[order(folds_in_files)], readRDS), dimension = "time")
  saveRDS(prediction_final, '../../data/pred_KNN_PCA_CV_FINAL.rds')
} else {
  prediction_final <- readRDS('../../data/pred_KNN_PCA_CV_FINAL.rds')
}

# ==============================================================================
# SECTION 4: KNN VALIDATION (SYMMETRIC TO DEEP ESD METRICS)
# ==============================================================================
message("\n=========================================")
message(" INITIATING KNN MODEL VALIDATION         ")
message("=========================================\n")

# Extract variables for metric calculation
data_list_orig <- y_obs
data_list_pred <- prediction_final
data_list_orig[['Data']][data_list_orig[['Data']] < 0] <- NA

message("Pre-aggregating monthly data...")
obs_mensual <- agregacion_mensual_rapida(data_list_orig)
pred_mensual <- agregacion_mensual_rapida(data_list_pred)

message("Generating master Ocean Mask...")
suma_total_obs <- colSums(obs_mensual$Data, na.rm = TRUE, dims = 1)
mask_maestra <- is.na(suma_total_obs) | suma_total_obs == 0 | is.na(obs_mensual$Data[1,,]) 

data_list_orig <- apply_sea_mask(data_list_orig, mask_maestra)
data_list_pred <- apply_sea_mask(data_list_pred, mask_maestra)
obs_mensual <- apply_sea_mask(obs_mensual, mask_maestra)
pred_mensual <- apply_sea_mask(pred_mensual, mask_maestra)

estaciones <- list("Annual" = 1:12, "DJF" = c(12, 1, 2), "MAM" = c(3, 4, 5), "JJA" = c(6, 7, 8), "SON" = c(9, 10, 11))

mapas_todos <- list()
df_corr_espacial_knn <- data.frame()
df_bias_pixels_knn <- data.frame()
df_intra_corr_knn <- data.frame()

for (s_name in names(estaciones)) {
  message("   Processing season: ", s_name)
  s_months <- estaciones[[s_name]]
  
  if (s_name == "Annual") {
    obs_ann <- aggregateGrid(obs_mensual, aggr.y = list(FUN = "sum", na.rm = TRUE))
    pred_ann <- aggregateGrid(pred_mensual, aggr.y = list(FUN = "sum", na.rm = TRUE))
    map_max <- 40 # Limit for Annual plots
  } else {
    obs_ann <- aggregateGrid(subsetGrid(obs_mensual, season = s_months), aggr.y = list(FUN = "sum", na.rm = TRUE))
    pred_ann <- aggregateGrid(subsetGrid(pred_mensual, season = s_months), aggr.y = list(FUN = "sum", na.rm = TRUE))
    
    # Dynamic map limits based on the specific season
    if (s_name == "DJF") map_max <- 7
    if (s_name == "MAM") map_max <- 12
    if (s_name == "JJA") map_max <- 17
    if (s_name == "SON") map_max <- 9
  }
  
  obs_clim <- climatology(obs_ann, clim.fun = list(FUN = "mean", na.rm = TRUE))
  pred_clim <- climatology(pred_ann, clim.fun = list(FUN = "mean", na.rm = TRUE))
  
  bias_map <- gridArithmetics(pred_clim, obs_clim, operator = "-")
  res_cor <- calc_cor_interanual(obs_ann, pred_ann)
  cor_map <- res_cor$cor_grid; puntos_sig <- res_cor$pts
  
  # Apply the master ocean mask
  obs_clim <- apply_sea_mask(obs_clim, mask_maestra)
  pred_clim <- apply_sea_mask(pred_clim, mask_maestra)
  bias_map <- apply_sea_mask(bias_map, mask_maestra)
  cor_map <- apply_sea_mask(cor_map, mask_maestra)
  
  obs_clim$Data[obs_clim$Data <= 0] <- NA
  pred_clim$Data[pred_clim$Data <= 0] <- NA
  
  # Calculate spatial mean metrics for map annotations
  mean_bias <- mean(bias_map$Data, na.rm = TRUE)
  mean_cor <- mean(cor_map$Data, na.rm = TRUE)
  
  x_coords <- bias_map$xyCoords$x
  y_coords <- bias_map$xyCoords$y
  lon_pos <- max(x_coords) - 0.02 * diff(range(x_coords))
  lat_pos <- min(y_coords) + 0.02 * diff(range(y_coords))
  
  txt_bias <- list("sp.text", c(lon_pos, lat_pos), sprintf("Mean: %.2f", mean_bias), cex = 1.6, adj = c(1, 0), font = 2)
  txt_cor <- list("sp.text", c(lon_pos, lat_pos), sprintf("Mean: %.2f", mean_cor), cex = 1.6, adj = c(1, 0), font = 2)
  
  sp_layout_bias <- list(txt_bias)
  sp_layout_cor <- list(puntos_sig, txt_cor)
  
  # Map generation
  p_obs <- spatialPlot(obs_clim, backdrop.theme = "countries", main = list(label=paste(s_name, "- Obs"), cex=2.5), 
                       color.theme = "RdYlBu", rev.colors = TRUE, set.min = 0, set.max = map_max, at = seq(0, map_max, length.out=21), par.settings = lattice_margins)
                       
  p_pred <- spatialPlot(pred_clim, backdrop.theme = "countries", main = list(label=paste(s_name, "- KNN"), cex=2.5), 
                        color.theme = "RdYlBu", rev.colors = TRUE, set.min = 0, set.max = map_max, at = seq(0, map_max, length.out=21), par.settings = lattice_margins)
                        
  p_bias <- spatialPlot(bias_map, backdrop.theme = "countries", main = list(label=paste(s_name, "- Bias"), cex=2.5), 
                        color.theme = "RdYlBu", set.min = -2, set.max = 2, at = seq(-2, 2, length.out=21), rev.colors = TRUE, 
                        sp.layout = sp_layout_bias, par.settings = lattice_margins)
                        
  p_cor <- spatialPlot(cor_map, backdrop.theme = "countries", main = list(label=paste(s_name, "- Cor"), cex=2.5), 
                       color.theme = "RdYlBu", set.min = -1, set.max = 1, at = seq(-1, 1, length.out=21), rev.colors = TRUE, 
                       sp.layout = sp_layout_cor, par.settings = lattice_margins)
  
  mapas_todos[[paste0(s_name, "_obs")]] <- p_obs
  mapas_todos[[paste0(s_name, "_pred")]] <- p_pred
  mapas_todos[[paste0(s_name, "_bias")]] <- p_bias
  mapas_todos[[paste0(s_name, "_cor")]] <- p_cor
  
  # Seasonal Metrics Calculation
  if (s_name != "Annual") {
    res_metricas <- calc_metricas_espaciales(obs_ann, pred_ann)
    
    df_spat <- res_metricas$spat
    df_spat$Season <- s_name
    df_corr_espacial_knn <- rbind(df_corr_espacial_knn, df_spat[, c("Year", "Spatial_Corr", "Season")])
    
    df_bias_pix <- res_metricas$bias_pixels
    df_bias_pix$Season <- s_name
    df_bias_pixels_knn <- rbind(df_bias_pixels_knn, df_bias_pix)
    
    df_diario_pixels <- calc_metricas_diarias_intra(data_list_orig, data_list_pred, mask_maestra, s_months, s_name)
    df_diario_pixels$Season <- s_name
    df_intra_corr_knn <- rbind(df_intra_corr_knn, df_diario_pixels)
  }
}

# Export KNN Maps
message("Saving combined KNN maps (PDF)...")
pdf(file = "mapas_metricas_knn/Mapas_Turc_KNN_Combinados.pdf", width = 24, height = 28)
grid.arrange(grobs = mapas_todos, ncol = 4, nrow = 5)
dev.off()

message("Saving combined high-resolution KNN maps (PNG)...")
png(file = "mapas_metricas_knn/Mapas_Turc_KNN_Combinados.png", width = 24, height = 28, units = "in", res = 300)
grid.arrange(grobs = mapas_todos, ncol = 4, nrow = 5)
dev.off()

# Export KNN Boxplots
df_bias_pixels_knn$Season <- factor(df_bias_pixels_knn$Season, levels = c("DJF", "MAM", "JJA", "SON"))
df_intra_corr_knn$Season <- factor(df_intra_corr_knn$Season, levels = c("DJF", "MAM", "JJA", "SON"))

bx_bias_global_knn <- ggplot(df_bias_pixels_knn, aes(x = as.factor(Year), y = Bias)) + 
  geom_boxplot(fill = "firebrick", alpha = 0.5, outlier.size = 0.5, outlier.alpha = 0.3) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1.2) + 
  stat_summary(fun = mean, geom = "point", color = "black", size = 2) + 
  facet_grid(Season ~ .) + 
  labs(title = "KNN: Pixel-by-Pixel Annual Bias Distribution", x = "Year", y = "Bias (Pred - Obs)") + 
  tema_overleaf

pdf(file = "boxplot_knn/Boxplot_Global_Bias_Estaciones_KNN.pdf", width = 18, height = 16)
print(bx_bias_global_knn)
dev.off()

ggsave(filename = "boxplot_knn/Boxplot_Global_Bias_Estaciones_KNN.png", 
       plot = bx_bias_global_knn, 
       width = 18, height = 16, units = "in", dpi = 300)

bx_corr_global_knn <- ggplot(df_intra_corr_knn, aes(x = as.factor(Year), y = Intra_Corr)) + 
  geom_boxplot(fill = "dodgerblue", alpha = 0.5, outlier.size = 0.5, outlier.alpha = 0.3) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey", linewidth = 1.2) + 
  stat_summary(fun = mean, geom = "point", color = "black", size = 2) + 
  facet_grid(Season ~ .) + 
  labs(title = "KNN: Daily Intra-Annual Correlation Distribution", x = "Year", y = "Correlation (r)") + 
  tema_overleaf

pdf(file = "boxplot_knn/Boxplot_Global_IntraCorr_Estaciones_KNN.pdf", width = 18, height = 16)
print(bx_corr_global_knn)
dev.off()

ggsave(filename = "boxplot_knn/Boxplot_Global_IntraCorr_Estaciones_KNN.png", 
       plot = bx_corr_global_knn, 
       width = 18, height = 16, units = "in", dpi = 300)

# ==============================================================================
# SECTION 5: COMPARISON: KNN vs DEEP ESD (2x2 SPATIAL CORRELATION PLOT)
# ==============================================================================
message("\n=========================================")
message(" GENERATING COMPARISON: KNN vs DEEP ESD  ")
message("=========================================\n")

# Load DeepESD data generated from the previous script
archivo_deep_esd <- "../DeepESD/series/Datos_Correlacion_Espacial_Estaciones.rds"

if (file.exists(archivo_deep_esd)) {
  df_corr_deep <- readRDS(archivo_deep_esd)
  df_corr_deep$Model <- "DeepESD"
  
  # Prepare KNN data
  df_corr_knn <- df_corr_espacial_knn
  df_corr_knn$Model <- "KNN"
  
  # Combine datasets for ggplot grouping
  df_compare <- rbind(df_corr_deep, df_corr_knn)
  df_compare$Season <- factor(df_compare$Season, levels = c("DJF", "MAM", "JJA", "SON"))
  
  # Create 2x2 facet grid comparison plot
  p_comparativa <- ggplot(df_compare, aes(x = Year, y = Spatial_Corr, color = Model)) +
    geom_line(linewidth = 1.5, alpha = 0.8) +
    geom_point(size = 3) +
    facet_wrap(~Season, nrow = 2, ncol = 2) +
    scale_color_manual(values = c("DeepESD" = "dodgerblue", "KNN" = "firebrick")) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "black", linewidth = 1) +
    labs(title = "Spatial Correlation Comparison: DeepESD vs KNN (By Season)",
         x = "Year", 
         y = "Pearson Correlation (r)") +
    tema_overleaf
  
  pdf(file = "series_knn/Comparativa_Correlacion_Espacial_Modelos.pdf", width = 16, height = 12)
  print(p_comparativa)
  dev.off()

  ggsave(filename = "series_knn/Comparativa_Correlacion_Espacial_Modelos.png", 
         plot = p_comparativa, 
         width = 16, height = 12, units = "in", dpi = 300)
  
} else {
  message("DeepESD file not found at: ", archivo_deep_esd)
}

message("\n>>> PROCESS COMPLETED!")