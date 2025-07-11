# Heterogeneous Treatment Effects Simulation Study

This repository contains R code for reproducing the main simulation study and analysis for the manuscript **Modelling uncertain heterogeneous treatment effects for decision analytical models: An early exploration**

## File Structure

### Main Files

- **`multi_sim_pop.R`**: Generates population data for all outcome surface scenarios. Creates trial datasets with varying control and treatment effect patterns (constant, linear, non-linear monotonic, non-linear non-monotonic).

- **`multi_sim_all.R`**: Main simulation loop. Fits five modeling approaches (unadjusted, adjusted, linear interaction, unrestricted spline, monotonic spline) to generated trial data and runs cohort simulations for cost-effectiveness analysis.

- **`MC_analysis.R`**: Monte Carlo analysis of simulation results. Computes posterior mean INMB, probability of cost-effectiveness, and log predictive density across all scenarios at willingness-to-pay threshold of £15,000.

- **`ceac_scenario.R`**: Generates cost-effectiveness acceptability curves (CEACs) across willingness-to-pay thresholds (£0-£30,000) for all scenarios and modeling approaches.

- **`appendix-plot.R`**: Creates side-by-side cost-effectiveness plane and outcome surface plots for manuscript appendix.

### Helper Functions

- **`data_generation.R`**: Functions for sampling ages, generating outcome probabilities, and creating trial datasets with specified control and treatment effect patterns.

- **`simulation_functions.R`**: Core simulation functions including Bayesian model fitting, cohort state-transition modeling, and cost-effectiveness calculations.

- **`analysis_functions.R`**: Analysis functions for computing INMB, probability of cost-effectiveness, log predictive density, and generating visualization plots.

## Workflow for Reproducing Results

### 1. Data Generation
Run `multi_sim_pop.R` to generate population data for all 12 scenarios (3 control × 4 treatment surface combinations) with 1,000 Monte Carlo replications each.

### 2. Main Simulation
Run `multi_sim_all.R` to execute the main simulation loop fitting all five modeling approaches to each scenario. Requires pre-compiled Stan models and parallel processing.

### 3. Results Analysis
Run `MC_analysis.R` to analyze simulation results and generate manuscript figures:
- Posterior mean INMB plots
- Log predictive density plots
- Probability of cost-effectiveness plots (limited data)

### 4. CEAC Analysis
Run `ceac_scenario.R` to generate cost-effectiveness acceptability curves across willingness-to-pay thresholds. Saves plots to `figures/ceac_all_limited.png` and `figures/ceac_all_extended.png`.

### 5. Appendix Plots
Run `appendix-plot.R` to create combined cost-effectiveness plane and outcome surface plots for manuscript appendix.

## Output

- **`figures/`**: Contains all plots used in the manuscript
- **`results/`**: Simulation results organized by modeling approach
- **`pop_data/`**: Generated population datasets for each scenario

## Important Notes

1. **Data Availability**: Due to file size limitations, population data and simulation results are not uploaded to this repository. The data are available upon request.

2. **Appendix Figure Compatibility**: To reproduce the appendix figure, ggplot2 version 3.4.x, ggdist 3.3.1, and cowplot 1.1.3 are required, as the relevant ggplot objects were created using older versions of these packages.
