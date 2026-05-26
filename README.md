# Turc Index & Climate Downscaling Intercomparison

This repository contains the research, code, and resources for my Master's Thesis in Data Science at the University of Cantabria.

## Project Overview

This study focuses on evaluating agricultural potential under future climate scenarios by emulating and downscaling the **Turc Index** of agricultural potential. The core of the research is a methodological intercomparison of statistical downscaling techniques, driven by 10 upper-air atmospheric variables. 

We compare two distinct approaches:
1.  **K-Nearest Neighbors (KNN):** A classic analog-based statistical method.
2.  **DeepESD:** A modern Deep Learning approach utilizing Convolutional Neural Networks (CNNs).

The models are evaluated under a "perfect prognosis" framework using reanalysis data (ERA5). Subsequently, we apply the trained models to various Global Climate Models (GCMs) to generate climate change projections under the **RCP8.5 scenario**, specifically analyzing impacts at Global Warming Levels (GWL) of +2.0ºC, +3.0ºC, and +4.0ºC.

## Repository Structure

Designed with open science and reproducibility in mind, this repository is organized as follows:

* **`Scripts/`**: Source code for data processing, model training (KNN & DeepESD), and evaluation metrics.
* **`Figures/`**: High-resolution charts, spatial maps, and visual intercomparisons presented in the final manuscript.
* **`main.ipynb/`**: Interactive Jupyter Notebook detailing the step-by-step methodology for reproducibility.

## Environment Installation

The required packages and dependencies to run the experiments are listed in `environment.yaml`. To set up the environment, please follow these steps:

1.  **Create the environment** (We recommend using Mamba for faster dependency resolution):
    ```bash
    conda create -n turc-downscaling -c conda-forge mamba
    ```

2.  **Activate the base environment:**
    ```bash
    conda activate turc-downscaling
    ```

3.  **Install the required packages using Mamba:**
    ```bash
    mamba env create -f environment.yaml
    ```

4.  **Activate the final environment:**
    ```bash
    conda activate turc-downscaling 
    ```

Once the environment is installed and activated, you will be able to open the Jupyter notebook to explore the emulation & downscaling process and climate change projections!
