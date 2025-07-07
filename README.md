# NPF_ML
## MATLAB scripts

### read_AOSnc_file_extract_10282024:
This script processes .nc and .cdf files by extracting data based on file date and category, initializing data fields within a structured .mat file, and saving the consolidated data to a dated .mat file, ensuring continuity across compatible dimensions.

### read_AOSnc_file_append_10282024: 
This MATLAB script processes .nc and .cdf files by locating, categorizing, and appending their data to existing .mat files, updating and saving each .mat file with new data fields.

### all_together_10282024_plot: 
This script processes data from .mat files, performing data alignment, interpolation, flag calculation, and data visualization, including data loading and preparation, data interpolation and transformation, flag calculation and post-processing, data plotting and export, data export to CSV.

### HYSPLIT_COlor_COded_Height_interval_10282024: 
This script visualizes backward atmospheric trajectories from multiple files by plotting latitude and longitude segments colored according to height above ground level (AGL). 

## Python scripts
### plotfrequency_10282024: 
This notebook analyzes and visualizes New Particle Formation occurrences from 2018 to 2023 using environmental data. It performs data preprocessing, grouping, and visualization to understand monthly patterns and trends.

### NPF_identification_10282024: 
This notebook is designed to identify the onset time of New Particle Formation events based on a Random Forest Classifier.

### NPF_driving_factors_10282024: 
This notebook aims to identify the driving factors behind New Particle Formation events using a Random Forest Classifier and feature importance analysis.

### NPF_driving_factors_10282024_predict(latest)_Gaurav
This notebook addresses class imbalance in the NPF dataset to improve the recall of the Random Forest model. It also implements an early warning approach to predict events before their occurrence and includes a feature importance analysis.

### NPF_ML_Models
This notebook focuses on employing various widely used ML models, such as Random Forest, XGboost, LightGBM, SVM, and MLP, and compares them across key criteria such as predictive accuracy, computational cost, model complexity, and interpretability.
