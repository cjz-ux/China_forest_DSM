# QRF Soil Properties Modeling

This repository contains scripts for soil property prediction using Quantile Random Forest (QRF), including hyperparameter tuning, model training, and spatial prediction with raster data.

## Project Structure

- `1_train_tuning/train_tune_qrf.R`: Batch training and hyperparameter tuning for QRF models.
- `2_raster_prediction/raster_qrf_predict.R`: Raster-based spatial prediction using optimized QRF models.
- `requirements_R.txt`: List of required R packages.
- `example_data/`: Example datasets for demonstration (optional).

## Usage

### 1. Environment Setup

Install R and all required packages:

```r
install.packages(c("raster", "randomForest", "caret", "quantregForest", "ggplot2", "readxl", "stringr"))
remotes::install_github("mlr-org/mlr3")
remotes::install_github("mlr-org/mlr3learners")
remotes::install_github("mlr-org/mlr3tuning")
remotes::install_github("mlr-org/paradox")
```

Or install each package listed in `requirements_R.txt`.

### 2. Data Preparation

- `_combined_data.csv`: Merged tabular data; first column is the response variable, remaining columns are predictors.
- `_best_features.rds`: R serialized file containing selected feature names.
- Raster stack files and layer name `.tif` and `.rds` files for spatial predictions.

### 3. Workflow

#### Step 1: Model Training & Hyperparameter Tuning

Edit and run `1_train_tuning/train_tune_qrf.R` to train models and generate optimized parameters and performance metrics.

#### Step 2: Raster-Based Prediction

Edit and run `2_raster_prediction/raster_qrf_predict.R` to predict soil properties spatially and output mosaicked results.

### 4. Outputs

- `*_optimized_performance.csv`: Performance and tuning results for each model.
- `*_pred.tif`, `*_lower.tif`, `*_upper.tif`: Quantile prediction raster files.
- `total_optimized_model_results.csv`: Summary table for all models.

## Contributing

Issues and pull requests are welcome for improvements or suggestions.

## License

Please add a LICENSE file (MIT, GPL, etc.) as appropriate for your project.
