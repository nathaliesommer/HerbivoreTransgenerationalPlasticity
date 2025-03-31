# Population climate history and predation risk, not environmental predictability, shape transgenerational plasticity 
[Manuscript in review]

Authors: Nathalie R. Sommer, Matthew S. Baker, Annise M. Dobson, Geoffrey C. Trussell, Oswald J. Schmitz

Code by Nathalie R. Sommer

This repository contains all data, scripts, and documentation necessary to reproduce the analyses and figures presented in our manuscript. Upon acceptance, this README will be updated with the manuscript DOI.

## Table of Contents

- [Project Overview](#project-overview)
- [Data](#data)
- [Scripts](#scripts)
- [Dependencies](#dependencies)
- [Reproducibility](#reproducibility)

## Project Overview

This research investigates transgenerational plasticity response to predation and temperature predictability metrics. The study examines how  conditions experienced by parent generations of grasshoppers influence offspring behavior and physiology, with implications for ecosystem dynamics and conservation.

## Data

All data is contained within the `Data` folder in their raw, original form. Descriptions of each file are contained below, with all data processing and analyses occurring in the [Scripts](#scripts).

### Behavior and Physiology Data
- `Data/Behavior_Raw.csv`: Raw behavioral data from experiments
- `Data/Behavior_Temps_2023_F1_F2_Raw.csv`: Temperature-dependent behavioral data for F1 and F2 generations
- `Data/Respriation_Mass_2023_F1_F2.csv`: Mass measurements for respiration
- `Data/Respriation_2023_F1_F2_Raw.csv`: Raw respiration data from F1 and F2 generations
- `Data/Environmental_Changes.csv`: Environmental change data
- `Data/Temperatures/`: Directory containing temperature-related data

## Scripts

All scripts are contained in the root directory. The main analysis scripts include:

- `Trait-Assessment.R`: Main script for processing and analyzing trait data, including:
  - Behavior analysis
  - Respiration analysis
  - Figure generation
- `Temperatures.R`: Script for processing and analyzing temperature data

## Dependencies

The analyses were conducted in R with the following packages:

### Data Processing and Analysis
- `dplyr` (v.1.1.4): Data manipulation and transformation
- `tidyr` (v.1.3.0): Data tidying and reshaping
- `lubridate` (v.1.9.2): Date and time handling
- `readr` (v.2.1.4): Data import
- `zoo` (v.1.8-12): Time series analysis
- `entropy` (v.1.6.1): Information theoretic measures

### Statistical Analysis
- `lme4` (v.1.1-35.5): Linear mixed-effects models
- `car` (v.3.1-30): Companion to Applied Regression
- `DHARMa` (v.0.4.6): Residual diagnostics for hierarchical models
- `emmeans` (v.1.8.0): Estimated marginal means
- `ks` (v.1.14.0): Kernel smoothing

### Visualization
- `ggplot2` (v.3.5.1): Data visualization
- `ggforce` (v.0.4.1): Additional plotting functionality
- `gghalves` (v.0.1.4): Half-violin plots

## Reproducibility

To reproduce the analysis and figures:

1. Clone the repository from GitHub:
   ```bash
   git clone https://github.com/nathaliesommer/HerbivoreTransgenerationalPlasticity.git
   ```
2. Set your working directory to the cloned repository
3. Run the scripts in the following order:
   - Start with `Temperatures.R` to process temperature data
   - Proceed with `Trait-Assessment.R` for trait analyses and visualization
4. Results will be generated in the working directory and printed to the console