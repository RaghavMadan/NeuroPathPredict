# NeuroPathPredict Pipeline (V1.1)

A modular pipeline for analyzing neuropathological data using Elastic Net and Universal Kriging approaches.

## Overview

This pipeline processes neuropathological data through several stages:
1. Data preprocessing and validation (including covariate transformation)
2. Elastic Net modeling
3. Factor analysis
4. Universal Kriging analysis

## Directory Structure

```
Model_V1.1/
├── config/
│   ├── input_variables.txt      # Input paths and parameters (paths relative to project root)
│   ├── input_variables.example.txt
│   ├── subject_list.txt         # Subject IDs (one per line)
│   └── subject_list.example.txt
├── scripts/                     # Pipeline scripts
├── input/
│   ├── raw/                     # Place your QNP and covariate CSVs here
│   └── processed/               # Processed inputs (e.g. MNI template)
└── output/                      # Generated outputs (not committed by default)
    └── {subject_id}/            # Subject-specific output
```

Run all steps from the **project root** (the directory containing `config/`, `scripts/`, `input/`, `output/`).

## Input Data Requirements

Place required input files in `input/raw/` (or set paths in `config/input_variables.txt`):
- `MFG_QNP_vox_index_combined.csv`: QNP data
- `Cov_all_consolidated_MFG_v1.csv`: Covariate data

## Pipeline Components

### 1. Preprocessing (`01_preprocessing.R`)
- Data loading and validation
- Basic cleaning and formatting
- Coordinate system alignment
- Missing value handling

### 2. Elastic Net (`02_elastic_net.R`)
- Feature selection and analysis
- Model training and validation
- Performance metrics calculation
- Feature importance visualization

### 3. Factor Analysis (`03_factor_analysis.R`)
- Predictor review and selection for downstream UK.

### 4. Universal Kriging (`04_universal_kriging.R`)
- 3D grid-based cross-validation
- Variogram modeling
- Spatial prediction
- Performance evaluation

### 5. Block Split Creation & NIfTI Visualization (`05_block_visualization.py`)
- Performs KMeans clustering on the subject's coordinate data for each specified split (e.g., 2, 4, 8, 16).
- **Mapping logic:** The input coordinates (X, Y, Z) are assumed to be in voxel space for the MNI152 2009b 0.5mm template. No affine transformation is applied; coordinates are rounded and used directly as voxel indices.
- For each split, a new column is added to the CSV (e.g., `2_block_splits`), and a NIfTI file is created with the block assignments.
- The output NIfTI files are written with the standard MNI affine in the header for compatibility.
- **Note:** If your input coordinates are in MNI mm space (not voxel indices), you must apply the affine transformation to convert to voxel indices. See the script docstring for details.

## Usage

1. From the project root, copy example configs if needed:
   - `cp config/input_variables.example.txt config/input_variables.txt`
   - `cp config/subject_list.example.txt config/subject_list.txt`
2. Edit `config/input_variables.txt`: input paths are relative to the project root (default: `input/raw/` for QNP and covariate CSVs).
3. Edit `config/subject_list.txt`: one subject ID per line.
4. Run the pipeline:
```bash
cd /path/to/Model_V1.1   # project root
Rscript scripts/run_pipeline.R
```
Output is written to `output/` and is not committed to the repo by default.

## Configuration

### Input Variables
- Subject ID
- Input/output paths
- Model parameters
- Cross-validation settings

### Subject List
- One subject ID per line
- Default: 6966

## Output

Each subject's output is organized in their respective directory:
- `preprocessing/`: Cleaned and validated data
- `covariate_analysis/`: Transformation analysis and plots
- `elastic_net/`: Model results and feature analysis
- `universal_kriging/`: Spatial analysis results
- `universal_kriging/block_visualizations/`: NIfTI files for each block split (e.g., `{subject_id}_blocks_4splits.nii.gz`).
- `universal_kriging/6966_uk_data.csv`: Updated with new block split columns for each split.

## Dependencies

**R** (install via `install.packages(...)` or run `Rscript scripts/setup.R`):
- dplyr, tidyr, readr, stringr, purrr, tibble
- data.table
- ggplot2
- caret, glmnet
- gstat, sp, raster
- logging, config
- moments, MASS
- yaml, RNifti
- furrr, tictoc

**Python** (for block visualization and some analysis scripts): numpy, pandas, scikit-learn, nibabel

## Error Handling

The pipeline includes:
- Input validation
- Process tracking
- Checkpointing
- Error logging
- Progress monitoring

## Contributing

Please read CONTRIBUTING.md for details on our code of conduct and the process for submitting pull requests.

## License

This project is licensed under the MIT License - see the LICENSE file for details. 