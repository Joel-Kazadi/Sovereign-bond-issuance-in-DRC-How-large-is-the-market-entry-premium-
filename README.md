# Eurobond DRC Project

This repository contains the full R pipeline for the empirical analysis of the inaugural Eurobond issued by the Democratic Republic of Congo (DRC). The project is organized into four sequential scripts covering preprocessing, exploratory data analysis, econometric estimation, and simulations.

The workflow is designed to be run **in order**. Each script builds on objects created by the previous one.

## Project overview

The empirical strategy is based on a cross-sectional comparison of first-time sovereign bond issuers worldwide. The main goal is to estimate the DRC’s **market-entry premium (MEP)** by comparing its observed spread with a model-implied counterfactual spread built from similar countries.

The pipeline produces:
- a cleaned analysis dataset,
- exploratory plots,
- weighted regression estimates and similarity diagnostics,
- counterfactual spread comparisons,
- bootstrap significance tests for the MEP,
- future learning-curve simulations,
- CSV and PDF outputs for reproducibility.

## Repository structure

```text
.
├── 1 Preprocessing.R
├── 2 EDA.R
├── 3 Estimation.R
├── 4 Simulations.R
├── eurobond_dataset.xlsx
├── plots/
├── simulation_outputs/
│   ├── plots/
│   └── tables/
└── README.md
```

## Data requirements

The pipeline expects an Excel workbook named `eurobond_dataset.xlsx`, with at least the following sheets:
- `cross sectional data` — main estimation dataset;
- `drc data` — future DRC projections used in the simulation stage.

The scripts assume the workbook is placed in the project root.

## Required R packages

The scripts use the following packages: `readxl`, `dplyr`, `stringr`, `lubridate`, `readr`, `tidyr`, `ggplot2`, `GGally`, `ggrepel`, `RColorBrewer`, `patchwork`, `scales`, `here`, `cluster`, `MASS`, `sandwich`, `lmtest`, `tibble`

Install any missing package before running the pipeline.

## How to run the pipeline

Run the scripts in this order:

- `1 Preprocessing.R`
- `2 EDA.R`
- `3 Estimation.R`
- `4 Simulations.R`

Each script assumes that the objects created by the previous step are already available in the R environment.

## Script-by-script guide

### 1 Preprocessing

#### Purpose

This script imports the raw Excel data and prepares the master cleaned dataset used in the rest of the analysis.

#### Outputs created

- In memory: `raw_df`, `clean_df`, `removed_df`, and `na_summary`.
- On disk: `cleaned_dataset.csv`.

#### Notes

This script is the foundation of the pipeline. The later scripts rely on `clean_df` being available in memory.

### 2 EDA

#### Purpose

This script performs the exploratory data analysis and generates the figures used in the EDA section of the paper.

#### Outputs created

PDF plots saved in plots/: `correlation_plot.pdf`, `coupon_vs_maturity.pdf`, `spread_vs_rating.pdf`, and `debt_distribution.pdf`.

#### Notes

This script is descriptive only. It does not modify the cleaned estimation dataset used later in the econometric stage.

### 3 Estimation

#### Purpose

This script estimates the locally weighted regression models used to reconstruct the DRC’s counterfactual spread and to build the similarity diagnostics.

#### Outputs created

- In memory: `estimation_store`, `estimation_summary`, and `similarity_plots`.
- On disk: `estimation_summary.csv`, and `similarity plots` in plots/, one PDF per combination.

#### Notes

This script creates the helper functions that the simulation stage later reuses. For that reason, `3 Estimation.R` should be sourced before `4 Simulations.R` if the scripts are run manually.

### 4 Simulations

#### Purpose

This script produces the simulation layer of the analysis:
- bootstrap significance tests for the MEP curve,
- comparison of the DRC’s MEP with the ten closest countries, and
- future learning-curve simulations under alternative scenarios.

#### Outputs created

- In memory: `ppe_significance_tbl`, `country_comparison_tbl`, `learning_curve_tbl`, `curve_plots`, and `learning_plots`.
- On disk: `simulation_outputs/tables/ppe_significance_summary.csv`, `simulation_outputs/tables/ppe_country_comparison_all_combos.csv`, `simulation_outputs/tables/learning_curve_all_combos.csv`, `individual comparison tables by combination in simulation_outputs/tables/`, `MEP significance plots in simulation_outputs/plots/`, and `learning-curve plots in simulation_outputs/plots/`.

#### Notes

This script depends on the estimation script because it reuses the fitted models and helper functions. If you run the scripts interactively, source `3 Estimation.R` first.

## Reproducibility

The repository is designed so that every major analytical step leaves a persistent output on disk. The main figures are saved as PDF files, and the main tables are exported as CSV files.

To reproduce the results:
1. place `eurobond_dataset.xlsx` in the repository root,
2. install the required R packages,
3. source the four scripts in order,
4. inspect the output folders `plots/` and `simulation_outputs/`.
