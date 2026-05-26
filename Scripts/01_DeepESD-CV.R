# ==============================================================================
# DEEP EMPIRICAL STATISTICAL DOWNSCALING (DeepESD) FOR CLIMATE DATA
# 
# Description:
# This script performs Deep Empirical Statistical Downscaling (DeepESD) using a 
# Convolutional Neural Network (CNN) architecture. It trains a Keras/TensorFlow 
# model to downscale ERA5 reanalysis data to a specific climate target (Turc Index).
# 
# The workflow includes:
#   1. Data loading and temporal intersection (climate4R framework).
#   2. 4-Fold cross-validation training loop with early stopping.
#   3. Custom Loss Function handling missing ocean data.
#   4. Reassembly of cross-validated predictions into a single seamless grid.
# ==============================================================================

# ==============================================================================
# SECTION 1: ENVIRONMENT SETUP & DEPENDENCIES
# ==============================================================================

# Increase Java memory BEFORE loading climate4R packages
options(java.parameters = "-Xmx8g") 
library(reticulate)

# CONNECT TO PYTHON 
path_to_python = "your_path_to_python"
use_python(path_to_python, required = TRUE)

# LIBRARY LOADING
library(tensorflow)
library(keras)
library(magrittr)
library(yaml)
library(downscaleR)
library(transformeR)
library(visualizeR)
library(downscaleR.keras)

# Read configuration parameters from YAML file
config = yaml.load_file("config.yaml")

# ==============================================================================
# SECTION 2: DATA LOADING & PREPROCESSING (PREDICTORS AND PREDICTAND)
# ==============================================================================

# Paths and variable configuration
path_data <- "../../data/ERA5/"
variables_base <- c("u", "v", "q", "zg", "t")
levels <- c("850", "500")
variables_list <- as.vector(outer(variables_base, levels, paste0))

# Create a list to store the loaded ERA5 grids
grid_list <- lapply(variables_list, function(var) {
  file_path <- paste0(path_data, var, "_ERA5_daily_1981_2020.rds")
  
  if (file.exists(file_path)) {
    return(readRDS(file_path))
  } else {
    warning(paste("File not found:", file_path))
    return(NULL)
  }
})

# Combine all variables into a single MultiGrid object
x_ERA5 <- do.call(makeMultiGrid, grid_list) %>% redim(drop = TRUE)
message("Process finished successfully. 'x_ERA5' object is ready.")

# Load Predictand (Turc Index)
y_obs = readRDS(file = '../../data/indexTurc_daily.rds')
attr(y_obs$Data, "dimensions") <- c("time", "lat", "lon")

# Filter out days entirely filled with NaNs or NAs
y_obs <- filterNA(y_obs)

# Substitute remaining NaNs (typically ocean pixels) with -999
y_obs[['Data']][is.na(y_obs[['Data']])] <- -999

# Initial intersection to ensure ERA5 predictors and observational data match perfectly in time
xT = intersectGrid(x_ERA5, y_obs, which.return = 1)
yT = intersectGrid(x_ERA5, y_obs, which.return = 2)

# ==============================================================================
# SECTION 3: MODEL ARCHITECTURE & CUSTOM LOSS FUNCTION
# ==============================================================================

# Define the DeepESD CNN Architecture
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

# Define a Custom Loss function that ignores ocean pixels (marked as -999)
loss_custom <- function(y_true, y_pred) {
  # Mask: 1 where >= 0 (land), 0 where < 0 (ocean)
  mask <- k_cast(k_greater_equal(y_true, 0), dtype = k_floatx())
  
  # Pixel-by-pixel squared error multiplied by the mask
  squared_error <- k_square(y_true - y_pred) * mask
  sum_squared_error <- k_sum(squared_error, axis = -1L)
  sum_mask <- k_sum(mask, axis = -1L) + k_epsilon()
  
  # MSE for each sample in the batch
  mse_batch <- sum_squared_error / sum_mask
  
  return(k_mean(mse_batch))
}
attr(loss_custom, "name") = "loss_custom"

# ==============================================================================
# SECTION 4: K-FOLD CROSS-VALIDATION SETUP
# ==============================================================================

# Directory preparation for outputs
dir.create("predicciones_tmp", showWarnings = FALSE)
dir.create("modelos_guardados", showWarnings = FALSE)

# Extract all available years in the dataset
years_all = unique(getYearsAsINDEX(xT)) 

# Create 4 sequential folds
k = 4
folds = split(years_all, cut(seq_along(years_all), breaks = k, labels = FALSE))

pred_list = list()
# Variables to preserve during loop cleanup to avoid memory leaks
vars_to_keep = c(ls(), "vars_to_keep", "i", "folds", "k")

# ==============================================================================
# SECTION 5: TRAINING & PREDICTION LOOP
# ==============================================================================

for (i in 1:k) {
    years_test = folds[[i]]
    years_train = setdiff(years_all, years_test)
    
    # Output file paths
    f_out = file.path("predicciones_tmp", paste0("prediccion_DeepESD_fold_", i, ".rds"))
    model_file = file.path("modelos_guardados", paste0("model_DeepESD_Fold_", i, ".h5"))
    
    # Skip if fold has already been processed
    if (file.exists(f_out) && file.exists(model_file)) {
        message("Skipping Fold ", i, ": Prediction and model already exist.")
        next
    }
    
    message("--- Starting Fold ", i, " ---")
    
    # --- DATA SPLIT (Train & Test) ---
    # Training set
    xT_train = subsetGrid(xT, years = years_train)
    yT_train = subsetGrid(yT, years = years_train)
    
    # Test set
    xT_test = subsetGrid(xT, years = years_test)

    # --- STANDARDIZATION ---
    # Standardize Train (ERA5) - Predictors only
    xT_train_scaled = scaleGrid(grid = xT_train,
                                type = "standardize",
                                spatial.frame = "gridbox") %>% redim(drop = TRUE) %>% suppressMessages() %>% suppressWarnings()

    xT_train_scaled[['Data']][is.na(xT_train_scaled[['Data']])] = 0

    # Standardize Test using Train (ERA5) mean and standard deviation
    xT_test_scaled = scaleGrid(grid = xT_test, 
                               base = xT_train,
                               type = "standardize",
                               spatial.frame = "gridbox") %>% redim(drop = TRUE) %>% suppressMessages() %>% suppressWarnings()

    xT_test_scaled[['Data']][is.na(xT_test_scaled[['Data']])] = 0
    
    # Prepare clean predictand grid for Keras
    yT_train_clean <- yT_train %>% redim(drop = TRUE)
    
    # --- PREPARE KERAS TENSORS ---
    data = prepareData.keras(x = xT_train_scaled, y = yT_train_clean,
                             spatial.predictors = NULL,
                             first.connection = "conv",
                             last.connection = "dense",
                             channels = "last")

    # SANITY CHECK: DIMENSIONS
    message("--- SHAPE VERIFICATION (FOLD ", i, ") ---")
    message("Full X tensor dimensions: ", paste(dim(data$x.global), collapse = " x "))
    message("Input shape to be passed to Keras (lat x lon x channels): ", paste(dim(data$x.global)[-1], collapse = " x "))
    
    # Instantiate the model
    model = architectures(input_shape = dim(data$x.global)[-1],
                          output_shape = dim(data$y$Data)[2])

    # --- COMPILE & TRAIN ---
    model %>% compile(
        loss = loss_custom,
        optimizer = optimizer_adam(learning_rate = config[['params']][['learning_rate']])
    )

    message("--- TRAINING SANITY CHECK ---")
    # Filter out -999 (oceans) to display valid land target metrics
    y_valida <- data$y$Data[data$y$Data > -900] 
    message("Target daily max in this Fold: ", max(y_valida, na.rm = TRUE))
    message("Target daily mean in this Fold: ", mean(y_valida, na.rm = TRUE))

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

    # --- PREDICTION (DOWNSCALING) ---
    newdata = prepareNewData.keras(xT_test_scaled, data)
    
    model_args = list(
        "filepath" = model_file,
        "custom_objects" = list("loss_custom" = loss_custom)
    )

    pred = downscalePredict.keras(
        newdata = newdata,
        model = model_args,
        C4R.template = subsetGrid(yT, years = years_test),
        clear.session = TRUE
    )

    # --- PREDICTION POST-PROCESSING ---
    # Restore NAs (oceans) based on the reference grid (Test)
    yT_test_template <- subsetGrid(yT, years = years_test)
    
    # Identify ocean masks
    mask_oceanos <- yT_test_template[['Data']] < 0
    
    # Apply mask to our continuous prediction
    pred[['Data']][mask_oceanos] <- NA

    message("--- PREDICTION SANITY CHECK ---")
    pred_valida <- pred[['Data']][pred[['Data']] > -900]
    message("Predicted daily max: ", max(pred_valida, na.rm = TRUE))
    message("Predicted daily mean: ", mean(pred_valida, na.rm = TRUE))
    
    # Save fold prediction
    saveRDS(pred, file = f_out)

    # --- CLEANUP ---
    try({ keras::k_clear_session() }, silent = TRUE)
    current_vars = ls()
    vars_to_remove = setdiff(current_vars, vars_to_keep)
    rm(list = vars_to_remove)
    gc()
    
    message("--- Fold ", i, " completed successfully ---")
}

# ==============================================================================
# SECTION 6: MERGE FINAL PREDICTIONS
# ==============================================================================

archivo_final <- 'pred_DeepESD_FINAL.rds'

if (file.exists(archivo_final)) {
  message("The file '", archivo_final, "' already exists. Skipping merging phase.")
} else {
  message("K-Fold loop finished. Merging predictions...")
  files = list.files("predicciones_tmp", pattern = "\\.rds$", full.names = TRUE)

  # Ensure numerical order of folds (1, 2, 3, 4)
  folds_in_files = as.numeric(gsub(".*_fold_(\\d+)\\.rds", "\\1", basename(files)))
  files = files[order(folds_in_files)]

  list_of_grids = lapply(files, readRDS)

  # Merge along the temporal dimension
  prediction_final = bindGrid(list_of_grids, dimension = "time")

  # Save unified result
  saveRDS(prediction_final, archivo_final)
  message("Process successfully finished! Final prediction saved.")
}

