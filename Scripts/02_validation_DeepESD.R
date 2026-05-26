# ==============================================================================
# CLIMATE METRICS AND VISUALIZATION (DeepESD vs OBSERVATIONS)
# 
# Description:
# This script evaluates the performance of the DeepESD downscaling model by 
# comparing its predictions against observational data (e.g., Turc Index). 
# It computes and visualizes spatial and temporal performance metrics across 
# different seasons (Annual, DJF, MAM, JJA, SON).
# 
# The workflow includes:
#   1. Fast monthly aggregations and Sea/Ocean masking.
#   2. Calculation of Climatologies, Mean Bias, and Interannual Correlation.
#   3. Calculation of Spatial Aggregation Metrics and Intra-annual Correlation.
#   4. Generation of high-resolution spatial maps (5x4 grids) for visualization.
#   5. Generation of global boxplots comparing pixel-by-pixel bias and correlation.
# ==============================================================================

# ==============================================================================
# SECTION 1: LIBRARIES & DIRECTORIES
# ==============================================================================
library(visualizeR)
library(transformeR)   
library(magrittr)
library(RColorBrewer)
library(ggplot2)       
library(gridExtra)     

# Create output directories if they do not exist
dir.create("mapas_metricas", showWarnings = FALSE)
dir.create("boxplot", showWarnings = FALSE)
dir.create("series", showWarnings = FALSE)

# ==============================================================================
# SECTION 2: CORE FUNCTIONS
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

# Calculates Pearson correlation and p-values pixel-by-pixel. Includes stippling for significance.
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
  
  pts <- map.stippling(climatology(pval_grid), threshold = threshold, condition = "LT", pch = 19, col = adjustcolor("black", alpha.f = 0.3), cex = 0.01) %>% suppressMessages() %>% suppressWarnings()
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

# Calculates temporal correlation based on daily data within specific seasons.
calc_metricas_diarias_intra <- function(obs_daily, pred_daily, mask_na, s_months, s_name) {
  message("   Extracting intra-annual daily correlation (pixel by pixel)...")
  
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
          o_valid <- o_vec[valid_idx]
          p_valid <- p_vec[valid_idx]
          sd_o <- sd(o_valid)
          sd_p <- sd(p_valid)
          if (!is.na(sd_o) && !is.na(sd_p) && sd_o > 0 && sd_p > 0) {
            res <- try(cor(o_valid, p_valid, method = "pearson"), silent = TRUE)
            if (!inherits(res, "try-error")) return(res)
          }
        }
        return(NA) 
      })
      
      tmp_df <- data.frame(Year = as.integer(y), Intra_Corr = cor_vals)
      df_intra <- rbind(df_intra, na.omit(tmp_df))
    }
  }
  return(df_intra)
}

# ==============================================================================
# SECTION 3: DATA LOADING AND PRE-PROCESSING
# ==============================================================================

message("Loading original daily data...")
data_list_orig <- readRDS(file = '../../data/indexTurc_daily.rds')
data_list_pred <- readRDS(file = '../../data/pred_DeepESD_CV.rds')

message("Pre-aggregating monthly data...")
obs_mensual <- agregacion_mensual_rapida(data_list_orig)
pred_mensual <- agregacion_mensual_rapida(data_list_pred)

message("Generating master Ocean Mask...")
suma_total_obs <- colSums(obs_mensual$Data, na.rm = TRUE, dims = 1)
mask_maestra <- is.na(suma_total_obs) | suma_total_obs == 0 | is.na(obs_mensual$Data[1,,]) 

# Apply the master mask to all grid instances
data_list_orig <- apply_sea_mask(data_list_orig, mask_maestra)
data_list_pred <- apply_sea_mask(data_list_pred, mask_maestra)
obs_mensual <- apply_sea_mask(obs_mensual, mask_maestra)
pred_mensual <- apply_sea_mask(pred_mensual, mask_maestra)

# Define target seasons
estaciones <- list("Annual" = 1:12, "DJF" = c(12, 1, 2), "MAM" = c(3, 4, 5), "JJA" = c(6, 7, 8), "SON" = c(9, 10, 11))

# Custom publication-ready ggplot theme (tailored for Overleaf/LaTeX integration)
tema_overleaf <- theme_minimal(base_size = 24) +
  theme(plot.title = element_text(face = "bold", size = 28, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 22, margin = margin(b = 15)),
        axis.title = element_text(face = "bold", size = 24),
        axis.text = element_text(size = 20, color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(face = "bold", size = 24), 
        legend.position = "bottom")

lattice_margins <- list(
  layout.heights = list(top.padding = 0, bottom.padding = 0, main.key.padding = 0),
  layout.widths = list(left.padding = 0, right.padding = 0)
)

# ==============================================================================
# SECTION 4: GLOBAL DATA STRUCTURES
# ==============================================================================
# Initialize structures to collect plot outputs and dataframes for global plotting
mapas_todos <- list()
df_corr_espacial_todos <- data.frame()
df_bias_pixels_todos <- data.frame()
df_intra_corr_todos <- data.frame()

# ==============================================================================
# SECTION 5: MAIN SEASONAL LOOP
# ==============================================================================

for (s_name in names(estaciones)) {
  message("\n--- PROCESSING SEASON: ", s_name, " ---")
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
  
  # Re-apply Ocean Mask to results
  obs_clim <- apply_sea_mask(obs_clim, mask_maestra)
  pred_clim <- apply_sea_mask(pred_clim, mask_maestra)
  bias_map <- apply_sea_mask(bias_map, mask_maestra)
  cor_map <- apply_sea_mask(cor_map, mask_maestra)
  
  # Calculate spatial mean metrics for annotations
  mean_bias <- mean(bias_map$Data, na.rm = TRUE)
  mean_cor <- mean(cor_map$Data, na.rm = TRUE)
  
  # Calculate bottom-right corner coordinates for text placement
  x_coords <- bias_map$xyCoords$x
  y_coords <- bias_map$xyCoords$y
  lon_pos <- max(x_coords) - 0.02 * diff(range(x_coords))
  lat_pos <- min(y_coords) + 0.02 * diff(range(y_coords))
  
  # Text objects for corner annotations
  txt_bias <- list("sp.text", c(lon_pos, lat_pos), sprintf("Mean: %.2f", mean_bias), cex = 1.6, adj = c(1, 0), font = 2)
  txt_cor <- list("sp.text", c(lon_pos, lat_pos), sprintf("Mean: %.2f", mean_cor), cex = 1.6, adj = c(1, 0), font = 2)
  
  # Package for layout visualization
  sp_layout_bias <- list(txt_bias)
  sp_layout_cor <- list(puntos_sig, txt_cor)

  # -----------------------------------------------------------------
  # 5.1 MAP GENERATION (All seasons including Annual)
  # -----------------------------------------------------------------
  message("   Generating maps for ", s_name, "...")
  
  p_obs <- spatialPlot(obs_clim, backdrop.theme = "countries", main = list(label=paste(s_name, "- Obs"), cex=2.5), 
                       color.theme = "RdYlBu", rev.colors = TRUE, set.min = 0, set.max = map_max, at = seq(0, map_max, length.out=21), par.settings = lattice_margins)
  
  p_pred <- spatialPlot(pred_clim, backdrop.theme = "countries", main = list(label=paste(s_name, "- Pred"), cex=2.5), 
                        color.theme = "RdYlBu", rev.colors = TRUE, set.min = 0, set.max = map_max, at = seq(0, map_max, length.out=21), par.settings = lattice_margins)
  
  p_bias <- spatialPlot(bias_map, backdrop.theme = "countries", main = list(label=paste(s_name, "- Bias"), cex=2.5), 
                        color.theme = "RdYlBu", set.min = -2, set.max = 2, at = seq(-2, 2, length.out=21), rev.colors = TRUE, 
                        sp.layout = sp_layout_bias, par.settings = lattice_margins)
  
  p_cor <- spatialPlot(cor_map, backdrop.theme = "countries", main = list(label=paste(s_name, "- Cor"), cex=2.5), 
                       color.theme = "RdYlBu", set.min = -1, set.max = 1, at = seq(-1, 1, length.out=21), rev.colors = TRUE, 
                       sp.layout = sp_layout_cor, par.settings = lattice_margins)
  
  # Store in the global list
  mapas_todos[[paste0(s_name, "_obs")]] <- p_obs
  mapas_todos[[paste0(s_name, "_pred")]] <- p_pred
  mapas_todos[[paste0(s_name, "_bias")]] <- p_bias
  mapas_todos[[paste0(s_name, "_cor")]] <- p_cor
  
  # -----------------------------------------------------------------
  # 5.2 METRICS CALCULATION (Seasonal only, excludes Annual)
  # -----------------------------------------------------------------
  if (s_name != "Annual") {

    # Spatial Metrics
    res_metricas <- calc_metricas_espaciales(obs_ann, pred_ann)
    df_spat <- res_metricas$spat
    df_spat$Season <- s_name
    
    # Collect spatial correlation for later export
    df_corr_espacial_todos <- rbind(df_corr_espacial_todos, df_spat[, c("Year", "Spatial_Corr", "Season")])
    
    # Collect pixel-by-pixel Bias
    df_bias_pix <- res_metricas$bias_pixels
    df_bias_pix$Season <- s_name
    df_bias_pixels_todos <- rbind(df_bias_pixels_todos, df_bias_pix)
    
    # Intra-annual Metrics
    df_diario_pixels <- calc_metricas_diarias_intra(data_list_orig, data_list_pred, mask_maestra, s_months, s_name)
    df_diario_pixels$Season <- s_name
    df_intra_corr_todos <- rbind(df_intra_corr_todos, df_diario_pixels)
  }
}

# ==============================================================================
# SECTION 6: JOINT FINAL EXPORT
# ==============================================================================
message("\n--- EXPORTING GLOBAL RESULTS ---")

# COMBINED MAPS (Grid: 5 Rows x 4 Columns)
message("Saving combined maps...")
pdf(file = "mapas_metricas/Mapas_Metricas.pdf", width = 24, height = 28)
grid.arrange(grobs = mapas_todos, ncol = 4, nrow = 5)
dev.off()

message("Saving high-resolution PNG maps...")
png(file = "mapas_metricas/Mapas_Metricas.png", width = 24, height = 28, units = "in", res = 300)
grid.arrange(grobs = mapas_todos, ncol = 4, nrow = 5)
dev.off()

# EXPORT SPATIAL CORRELATION DATA (Raw dataframe)
message("Exporting seasonal spatial correlation data...")
saveRDS(df_corr_espacial_todos, file = "series/Datos_Correlacion_Espacial_Estaciones.rds")

# COMBINED ROW-BASED BOXPLOTS (Using facet_grid)
message("Generating global seasonal boxplot figures...")

# Ensure strict chronological order for factors in the plot
df_bias_pixels_todos$Season <- factor(df_bias_pixels_todos$Season, levels = c("DJF", "MAM", "JJA", "SON"))
df_intra_corr_todos$Season <- factor(df_intra_corr_todos$Season, levels = c("DJF", "MAM", "JJA", "SON"))

# Global Bias Boxplot (4 rows stacked)
bx_bias_global <- ggplot(df_bias_pixels_todos, aes(x = as.factor(Year), y = Bias)) + 
  geom_boxplot(fill = "firebrick", alpha = 0.5, outlier.size = 0.5, outlier.alpha = 0.3) + 
  geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1.2) + 
  stat_summary(fun = mean, geom = "point", color = "black", size = 2) + 
  facet_grid(Season ~ .) + 
  labs(title = "Pixel-by-Pixel Annual Bias Distribution (By Season)", 
       x = "Year", y = "Bias (Pred - Obs)") + 
  tema_overleaf

pdf(file = "boxplot/Boxplot_Global_Bias_Seasons.pdf", width = 18, height = 16)
print(bx_bias_global)
dev.off()

ggsave(filename = "boxplot/Boxplot_Global_Bias_Seasons.png", 
       plot = bx_bias_global, 
       width = 18, height = 16, units = "in", dpi = 300)

# Global Intra-annual Correlation Boxplot
bx_corr_global <- ggplot(df_intra_corr_todos, aes(x = as.factor(Year), y = Intra_Corr)) + 
  geom_boxplot(fill = "dodgerblue", alpha = 0.5, outlier.size = 0.5, outlier.alpha = 0.3) + 
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey", linewidth = 1.2) + 
  stat_summary(fun = mean, geom = "point", color = "black", size = 2) + 
  facet_grid(Season ~ .) + 
  labs(title = "Daily Intra-Annual Correlation Distribution (By Season)", 
       x = "Year", y = "Correlation (r)") + 
  tema_overleaf

pdf(file = "boxplot/Boxplot_Global_IntraCorr_Seasons.pdf", width = 18, height = 16)
print(bx_corr_global)
dev.off()

ggsave(filename = "boxplot/Boxplot_Global_IntraCorr_Seasons.png", 
       plot = bx_corr_global, 
       width = 18, height = 16, units = "in", dpi = 300)

message("\nSUCCESS: All processes completed!")