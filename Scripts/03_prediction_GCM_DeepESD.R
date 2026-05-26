# ==============================================================================
# DEEP EMPIRICAL STATISTICAL DOWNSCALING (DeepESD) - GCM PROJECTIONS
# 
# Description:
# This script executes a complete DeepESD pipeline for climate change projections.
# It trains a Convolutional Neural Network (CNN) on historical reanalysis data 
# (ERA5) to downscale a climate target (Turc Index). It then applies the trained 
# model to Global Climate Model (GCM) outputs for historical and future scenarios.
# 
# The workflow includes:
#   1. ERA5 & Predictand data loading and temporal intersection.
#   2. CNN Architecture definition and custom loss function (masking oceans).
#   3. Model training on the complete ERA5 series and re-substitution prediction.
#   4. Advanced GCM pre-processing: Spatial NA filling via nearest-neighbor 
#      and bias correction via Scaling Delta Mapping.
#   5. GCM Downscaling: Date harmonization across variables and final projection.
# ==============================================================================

# ==============================================================================
# SECTION 0: INITIAL SETUP & LIBRARIES
# ==============================================================================
# Increase Java memory BEFORE loading climate4R packages to handle 32GB RAM
options(java.parameters = "-Xmx32g") 
library(reticulate)

# CONNECT TO PYTHON
# (Note for GitHub: Update this path to your local Python environment)
path_to_python = "your_path_to_python"
use_python(path_to_python, required = TRUE)

# LIBRARY LOADING
library(tensorflow)
library(keras)
library(magrittr)
library(yaml)
library(downscaleR)
library(transformeR)
library(downscaleR.keras)

# Read configuration parameters from YAML file 
config = yaml.load_file("config.yaml")

# ==============================================================================
# SECTION 1: LOAD ERA5 PREDICTORS AND OBSERVATIONAL TARGET
# ==============================================================================
message(">>> Loading observations and ERA5 predictors...")
path_data <- "../../data/ERA5/"
variables_base <- c("u", "v", "q", "zg", "t")
levels <- c("850", "500")
variables_list <- as.vector(outer(variables_base, levels, paste0))

# Load ERA5 variables into a list
grid_list <- lapply(variables_list, function(var) {
  file_path <- paste0(path_data, var, "_ERA5_daily_1981_2020.rds")
  if (file.exists(file_path)) {
    return(readRDS(file_path))
  } else {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
})

# Create the ERA5 MultiGrid
x_ERA5 <- do.call(makeMultiGrid, grid_list) %>% redim(drop = TRUE)
message("Process finished successfully. 'x_ERA5' object is ready.")

# Load Predictand (Turc Index)
y_obs = readRDS(file = '../../data/indexTurc_daily.rds')
attr(y_obs$Data, "dimensions") <- c("time", "lat", "lon")

# Filter out empty days and replace remaining NAs (oceans) with -999
y_obs <- filterNA(y_obs)
y_obs[['Data']][is.na(y_obs[['Data']])] <- -999

# Intersection to ensure ERA5 and observations match perfectly in time
xT = intersectGrid(x_ERA5, y_obs, which.return = 1)
yT = intersectGrid(x_ERA5, y_obs, which.return = 2)

# Free up the memory of the massive un-intersected grids
rm(x_ERA5, y_obs, grid_list)
gc()

# ==============================================================================
# SECTION 2: CNN MODEL DEFINITION
# ==============================================================================

architectures = function(input_shape, output_shape) {
    inputs = layer_input(shape = input_shape)
    x = inputs
    l1 = layer_conv_2d(x , filters = 50, kernel_size = c(3,3), activation = 'relu', padding = "valid")
    l2 = layer_conv_2d(l1, filters = 25, kernel_size = c(3,3), activation = 'relu', padding = "valid")
    l3 = layer_conv_2d(l2, filters = 10, kernel_size = c(3,3), activation = 'relu', padding = "valid")
    l4 = layer_flatten(l3)
    outputs = layer_dense(l4, units = output_shape, activation = 'relu')
    model = keras_model(inputs = inputs, outputs = outputs)
    return(model)
}

# Custom Loss function to ignore ocean pixels (values >= 0)
loss_custom <- function(y_true, y_pred) {
  mask <- k_cast(k_greater_equal(y_true, 0), dtype = k_floatx())
  squared_error <- k_square(y_true - y_pred) * mask
  sum_squared_error <- k_sum(squared_error, axis = -1L)
  sum_mask <- k_sum(mask, axis = -1L) + k_epsilon()
  mse_batch <- sum_squared_error / sum_mask
  return(k_mean(mse_batch))
}
attr(loss_custom, "name") = "loss_custom"

# Directory preparation for saving weights
dir.create("modelos_guardados", showWarnings = FALSE)

# ==============================================================================
# SECTION 3: STANDARDIZATION AND SINGLE MODEL TRAINING (ERA5 BASELINE)
# ==============================================================================
model_file = file.path("modelos_guardados", "model_DeepESD_FINAL.h5")

message("--- Preparing base data structure (ERA5) ---")
# Standardize ERA5. CRITICAL: Keep `xT` intact as it will be used as the 
# reference for the GCM Delta Mapping later.
xT_scaled = scaleGrid(grid = xT, type = "standardize", spatial.frame = "gridbox") %>% redim(drop = TRUE) %>% suppressMessages() %>% suppressWarnings()
xT_scaled[['Data']][is.na(xT_scaled[['Data']])] = 0

yT_clean <- yT %>% redim(drop = TRUE)

# Prepare Keras Matrices. This 'data' structure MUST be kept alive 
# to act as a template for the GCM later.
data = prepareData.keras(x = xT_scaled, y = yT_clean,
                         spatial.predictors = NULL,
                         first.connection = "conv",
                         last.connection = "dense",
                         channels = "last")

# Free the scaled version now safely inside the 'data' list
rm(xT_scaled)
gc()

# --- MODEL EXISTENCE CHECK ---
if (file.exists(model_file)) {
    message("\nThe model already exists at '", model_file, "'!")
    message("Loading the model and proceeding directly to GCM predictions.\n")
    model <- load_model_hdf5(model_file, custom_objects = c("loss_custom" = loss_custom))
} else {
    message("\n--- Training DeepESD model using the full ERA5 series ---")
    
    # Instantiate the model
    model = architectures(input_shape = dim(data$x.global)[-1],
                          output_shape = dim(data$y$Data)[2])

    model %>% compile(
        loss = loss_custom,
        optimizer = optimizer_adam(learning_rate = config[['params']][['learning_rate']])
    )

    # Train the model
    hist = model %>% fit(
        x = data$x.global,
        y = data$y$Data,
        batch_size = config[['params']][['batch']],
        epochs = config[['params']][['epochs']],
        validation_split = 0.2, 
        verbose = 1,
        callbacks = list(
            callback_early_stopping(patience = config[['params']][['patience']]), 
            callback_model_checkpoint(
                filepath = model_file,
                monitor = 'val_loss',
                save_best_only = TRUE,
                save_weights_only = FALSE
            )
        )
    )
    message("Model trained and saved to: ", model_file, "!")
    rm(hist)
    gc() 
}

# ==============================================================================
# SECTION 3.5: FINAL PREDICTION ON ERA5 (RE-SUBSTITUTION)
# ==============================================================================
file_pred_era5 <- "modelos_guardados/pred_DeepESD_ERA5_Final.rds"

if (!file.exists(file_pred_era5)) {
  message("\n--- Running final prediction on ERA5 predictors ---")
  
  # SCALE DATA (Crucial for Keras)
  message("  -> Standardizing ERA5 predictors...")
  xT_scaled_pred <- scaleGrid(grid = xT, type = "standardize", spatial.frame = "gridbox") %>% 
    transformeR::redim(drop = TRUE) %>% 
    suppressMessages() %>% suppressWarnings()
  
  # Fill NAs with 0 (mean value in a standardized distribution)
  xT_scaled_pred[['Data']][is.na(xT_scaled_pred[['Data']])] = 0
  
  # Prepare data into Keras tensor format using our previous template
  message("  -> Preparing tensors...")
  newdata_era5 <- prepareNewData.keras(
    newdata = xT_scaled_pred, 
    data.structure = data 
  )
  
  # Configure saved model access
  model_args_era5 <- list(
    "filepath" = model_file,
    "custom_objects" = list("loss_custom" = loss_custom)
  )
  
  # Execute prediction
  message("  -> Executing model...")
  pred_era5_final <- downscalePredict.keras(
    newdata = newdata_era5,
    model = model_args_era5,
    C4R.template = yT_clean,
    clear.session = TRUE
  )
  
  # Save output and clean RAM
  saveRDS(pred_era5_final, file = file_pred_era5, compress = "xz")
  message(">>> ERA5 final prediction successfully saved to: ", file_pred_era5)
  
  rm(xT_scaled_pred, newdata_era5, pred_era5_final)
  gc()
} else {
  message("\nThe ERA5 final prediction file already exists at '", file_pred_era5, "'! Skipping...")
}

# ==============================================================================
# SECTION 4: DELTA MAPPING AND SPATIAL INTERPOLATION FUNCTIONS
# ==============================================================================
message(">>> Loading harmonization functions...")

scalingDeltaMapping <- function(grid, base, ref) {
  # Remove the seasonal trend
  grid_detrended <- scaleGrid(grid, 
                              base = grid, 
                              ref = base, 
                              type = "center", 
                              spatial.frame = "gridbox", 
                              time.frame = "monthly",
                              skip.season.check = TRUE)  
  
  # Bias correct the mean and variance
  grid_detrended_corrected <- scaleGrid(grid_detrended, 
                                        base = base, 
                                        ref = ref, 
                                        type = "standardize", 
                                        spatial.frame = "gridbox", 
                                        time.frame = "monthly",
                                        skip.season.check = TRUE)    
  
  # Add the seasonal trend back
  grid_corrected <- scaleGrid(grid_detrended_corrected, 
                              base = base, 
                              ref = grid, 
                              type = "center", 
                              spatial.frame = "gridbox", 
                              time.frame = "monthly",
                              skip.season.check = TRUE)    
  return(grid_corrected)
}

standardize <- function(grid, base){
  scaleGrid(grid, 
            base = base,
            ref = NULL,
            type = "standardize",
            skip.season.check = TRUE
  ) 
} 

# ------------------------------------------------------------------------------
# NEAREST-NEIGHBOR SPATIAL NA FILLING (Crucial for GCM topologies)
# ------------------------------------------------------------------------------
fill_all_nas_multigrid <- function(multigrid) {
  # Ensure clean structure
  multigrid <- transformeR::redim(multigrid, drop = TRUE)
  data_array <- multigrid$Data
  
  # Dimensions: [variables, time, lat, lon]
  n_vars <- dim(data_array)[1]
  n_lat <- dim(data_array)[3]
  n_lon <- dim(data_array)[4]
  
  # Relative coordinate matrices
  lat_matrix <- matrix(rep(1:n_lat, n_lon), nrow = n_lat, ncol = n_lon)
  lon_matrix <- matrix(rep(1:n_lon, each = n_lat), nrow = n_lat, ncol = n_lon)
  
  # Clean variable by variable
  for (v in 1:n_vars) {
    # Take the first day as a reference (geographical NAs are usually static)
    capa_ref <- data_array[v, 1, , ]
    nas <- is.na(capa_ref)
    
    if (any(nas) && !all(nas)) {
      message(sprintf("Filling NAs for variable %d...", v))
      valid_idx <- which(!nas)
      na_idx <- which(nas)
      
      # Find the nearest valid neighbor for each NA pixel
      nearest_valid_idx <- integer(length(na_idx))
      for (k in seq_along(na_idx)) {
        i <- na_idx[k]
        distancias <- (lat_matrix[i] - lat_matrix[valid_idx])^2 + (lon_matrix[i] - lon_matrix[valid_idx])^2
        nearest_valid_idx[k] <- valid_idx[which.min(distancias)]
      }
      
      # Replace NA with its neighbor simultaneously across ALL time steps
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

# ==============================================================================
# SECTION 5: GCM PROCESSING & PREDICTION (HISTORICAL AND RCP8.5)
# ==============================================================================
message(">>> Processing GCM data independently (Historical and RCP8.5)...")
dir.create("proyecciones_ESD", showWarnings = FALSE, recursive = TRUE)

model_name <- "HadGEM2-ES" 
file_path <- paste0("../../data/GCMs/GCM_data_", model_name, ".rds")
gcm_data <- readRDS(file_path)

# Source nomenclature (GCM)
base_vars_gcm <- c("ua", "va", "hus", "zg", "ta") 
levels_gcm <- c("85000", "50000")
var_list_gcm_target <- as.vector(outer(base_vars_gcm, levels_gcm, paste, sep="_"))

# Target nomenclature (ERA5 - Expected by Keras)
base_vars_era5 <- c("u", "v", "q", "zg", "t")
levels_era5 <- c("850", "500")
var_list_era5_target <- as.vector(outer(base_vars_era5, levels_era5, paste0))

# Global variable to store the historical base used by Delta Mapping in both periods
harmonize.base <- NULL 

# ==============================================================================
# MAIN SCENARIO LOOP
# ==============================================================================
periods <- c("historical", "rcp85")

for (period_name in periods) {
  message(sprintf("\n--- Processing period: %s ---", period_name))
  
  output.file <- sprintf("proyecciones_ESD/Turc_DeepESD_%s_%s.rds", model_name, period_name)
  
  if (!file.exists(output.file)) {
    
    # SORTING AND RENAMING MULTIGRID CHANNELS
    message("Sorting and cleaning channels for Keras...")
    grid_list_gcm <- list()
    
    for (i in seq_along(var_list_gcm_target)) {
      k <- var_list_gcm_target[i] 
      new_name <- var_list_era5_target[i]
      grid <- gcm_data[[period_name]][[k]]
      
      if (!is.null(grid)) {
        
        # Exorcise CMIP ghost values (outliers)
        grid$Data[grid$Data > 1e10] <- NA
        grid$Data[grid$Data < -1e10] <- NA
        
        # Rename and interpolate to ERA5 grid
        grid$Variable$varName <- new_name
        grid <- interpGrid(grid, new.coordinates = getGrid(xT), method = "bilinear")
        
        grid_list_gcm[[new_name]] <- grid
      } else {
        stop(paste("Critical Error: Missing GCM data for variable:", k, "in", period_name))
      }
    }
    
    # ====================================================================
    # MANUAL DATE HARMONIZATION (WITH SANITY CHECKS)
    # ====================================================================
    message("Calculating the exact dates shared by ALL variables...")
    
    # Extract date vectors IGNORING THE HOUR (YYYY-MM-DD only)
    dates_list <- lapply(grid_list_gcm, function(g) {
      substr(as.character(g$Dates$start), 1, 10)
    })
    
    # Original length of dates per variable
    message("\n--- Original days per variable ---")
    for (name in names(grid_list_gcm)) {
      message(sprintf("Variable %s: %d original days", name, length(dates_list[[name]])))
    }
    
    # Find the pure intersection (now strictly day-based)
    common_dates <- Reduce(intersect, dates_list)
    message(sprintf("\n -> Expected final common days: %d", length(common_dates)))
    
    # Prevent crashing down the line
    if (length(common_dates) == 0) {
      stop("CRITICAL ERROR: Zero common days found. Date ranges for humidity (hus/q) and other variables do not overlap at all in this period.")
    }
    
    # Filter each variable to strictly retain these common days
    message("\n--- Subsetting variables and verifying dimensions ---")
    grid_list_gcm_harmonized <- lapply(names(grid_list_gcm), function(name) {
      g <- grid_list_gcm[[name]]
      
      # Strip hour again to ensure a correct match
      fechas_orig <- substr(as.character(g$Dates$start), 1, 10)
      
      idx <- match(common_dates, fechas_orig)
      idx <- idx[!is.na(idx)]
      
      dim_antes <- dim(g$Data)
      message(sprintf(" -> %s | Dim BEFORE: [%s] | Indices found: %d", 
                      name, paste(dim_antes, collapse=" x "), length(idx)))
      
      # Subset along the time dimension
      g_sub <- transformeR::subsetDimension(g, dimension = "time", indices = idx)
      
      dim_despues <- dim(g_sub$Data)
      message(sprintf("    %s | Dim AFTER: [%s]", 
                      name, paste(dim_despues, collapse=" x ")))
      
      return(g_sub)
    })
    
    # Restore list names
    names(grid_list_gcm_harmonized) <- names(grid_list_gcm)
    
    # Final overview of the time dimension before collapsing
    message("\n--- Final check of the 'time' dimension ---")
    time_dims <- sapply(grid_list_gcm_harmonized, function(g) {
      # Search for the time dimension based on attributes
      idx_time <- which(attr(g$Data, "dimensions") == "time")
      return(dim(g$Data)[idx_time])
    })
    print(time_dims)
    
    if (length(unique(time_dims)) > 1) {
      stop("CRITICAL ERROR: Time dimensions are still mismatched. Review the sanity checks above.")
    }
    
    # Create the MultiGrid with the harmonized list
    message("\nCreating MultiGrid...")
    args_multigrid <- unname(grid_list_gcm_harmonized)
    args_multigrid$skip.temporal.check <- TRUE
    x_raw <- do.call(makeMultiGrid, args_multigrid) %>% redim(drop = TRUE)
    
    # Fill spatial NAs
    message("Filling GCM spatial NAs (Vital for Keras to process properly)...")
    x_raw <- fill_all_nas_multigrid(x_raw)
    
    # ====================================================================
    # REFERENCE CALCULATION AND SCALING DELTA MAPPING
    # ====================================================================
    if (period_name == "historical") {
      # Calculate the REAL intersection of surviving years
      years_era5 <- unique(getYearsAsINDEX(xT))
      years_gcm_raw <- unique(getYearsAsINDEX(x_raw))
      
      # Save globally to ensure RCP8.5 uses the historical baseline
      ref.period <<- intersect(years_era5, years_gcm_raw) 
      
      message(sprintf("REAL reference period for bias correction: %d to %d", min(ref.period), max(ref.period)))
      
      harmonize.base <<- subsetGrid(x_raw, years = ref.period)
    }
    
    message("Extracting common reference period for ERA5...")
    xT_ref <- subsetGrid(xT, years = ref.period) %>% suppressMessages()
    
    message("Applying Scaling Delta Mapping...")
    x_harm <- scalingDeltaMapping(grid = x_raw, base = harmonize.base, ref = xT_ref)
    
    message("Standardizing relative to ERA5 (1981-2020)...")
    xn <- standardize(grid = x_harm, base = xT)
    xn[['Data']][is.na(xn[['Data']])] <- 0

    # KERAS MODEL PREDICTION (VECTORIZED/NO MANUAL LOOPS)
    message("Preparing tensors and performing downscaling...")
    newdata <- prepareNewData.keras(newdata = xn, data.structure = data)
    
    # Create a dummy template for Climate4R to assemble the 3D grid natively
    dummy_template <- yT_clean
    n_time_gcm <- dim(xn$Data)[2]
    dummy_template$Data <- array(NA, dim = c(n_time_gcm, dim(yT_clean$Data)[2], dim(yT_clean$Data)[3]))
    dummy_template$Dates <- xn$Dates
    attr(dummy_template$Data, "dimensions") <- c("time", "lat", "lon")
    
    model_args <- list(
      "filepath" = model_file,
      "custom_objects" = list("loss_custom" = loss_custom)
    )

    pred_grid <- downscalePredict.keras(
      newdata = newdata,
      model = model_args,
      C4R.template = dummy_template,
      clear.session = TRUE
    )

    # SAVE AND CLEANUP
    saveRDS(pred_grid, file = output.file, compress = "xz")
    message(">>> Prediction successfully saved to: ", output.file)
    
    rm(grid_list_gcm, x_raw, x_harm, xn, newdata, dummy_template, pred_grid)
    gc()
    
  } else {
    message(sprintf("The file %s already exists. Skipping...", output.file))
  }
}

message("\n>>> Downscaling workflow completed successfully!")