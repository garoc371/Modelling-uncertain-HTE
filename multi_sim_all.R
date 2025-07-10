library(rstanarm)
library(tidyverse)
library(splines2)
library(cmdstanr)
library(tidybayes)
library(patchwork)
library(doParallel)
library(doSNOW)
library(tidyr)
library(dplyr)
library(here)

source(here("microsim_helper.R"))
source(here("cmdstan_fit_plot.R"))
source(here("gen_surface.R"))

source(here("sim_multi_cohort.R"))
source(here("cohortsim_helper.R"))

# Generate outcome surfaces for all scenarios
n_trial <- 500
n_extra <- 200
n_target <- 20000

age_grid <- 40:80
extended_age_grid <- 40:120
limited_probs <- uniform_probs(age_grid)
extra_probs <- half_normal_probs(age_grid, sd = 5)
# Generate ages for control and treatment groups
ages_treatment <- ages_control <- sample_age(n_trial, limited_probs, age_grid)

extra_trt <- extra_control <- sample_age(n_extra, probs = extra_probs, age_grid)
age_target <- sample(40:80, size = n_target, replace = T)


# Utility values
utility_healthy <- 1
utility_diseased <- 0.7
utility_dead <- 0

# Costs per cycle
cost_healthy <- 100
cost_diseased <- 500
cost_dead <- 0

treatment_cost_per_cycle <- 5000

# Time horizon and number of cycles
time_horizon <- 40


n_simulations <- 4000
n_patients <- n_target
n_cycles <- 40

baseline_probs = 0.4

utilities <- c(utility_healthy, utility_diseased, 0)
cost_ctrl <- c(cost_healthy, cost_diseased, 0)
cost_trt <- c(cost_healthy + treatment_cost_per_cycle, cost_diseased, 0)

# Define all combinations for control and treatment surfaces
combinations <- expand.grid(
  control_type = c(
    "non_linear_mon_increase",
    "non_linear_mon_decrease",
    "non_linear_non_mon"
  ),
  treatment_type = c(
    "constant",
    "non_linear_mon_increase",
    "non_linear_mon_decrease",
    "non_linear_non_mon"
  )
)

# Iterate through the combinations

rw_file <- file.path(cmdstan_path(), "rw_spline.stan")
mod_sep_rw <- cmdstan_model(rw_file)

sep_file <- file.path(cmdstan_path(), "sep_dir.stan")
mod_sep_dir <- cmdstan_model(sep_file)

n_rep <- 1000


# Add short versions to the combinations data frame
combinations <- combinations %>%
  mutate(
    control_short = gsub(
      "non_linear",
      "non",
      gsub("linear", "lin", as.character(control_type))
    ),
    treatment_short = gsub(
      "non_linear",
      "non",
      gsub("linear", "lin", as.character(treatment_type))
    ),
    scenario = paste0(control_short, "_ctrl_", treatment_short, "_trt")
  )


for (i in 1:nrow(combinations)) {
  # Extract the current control and treatment types
  control_type <- missing_combinations$control_type[i]
  treatment_type <- missing_combinations$treatment_type[i]

  # Generate the filename
  control_short <- gsub(
    "non_linear",
    "non",
    gsub("linear", "lin", as.character(control_type))
  )
  treatment_short <- gsub(
    "non_linear",
    "non",
    gsub("linear", "lin", as.character(treatment_type))
  )
  pop_path <- here(
    "pop_data",
    paste0(
      control_short,
      "_ctrl_",
      treatment_short,
      "_trt"
    )
  )
  plot_path <- here(
    "figures",
    paste0(
      control_short,
      "_ctrl_",
      treatment_short,
      "_trt"
    )
  )
  res_true_path <- here(
    "results",
    "true",
    paste0(
      control_short,
      "_ctrl_",
      treatment_short,
      "_trt"
    )
  )

  cohort_true <- sim_true_cohort(control_type, treatment_type)
  true_weights <- cohort_true$target_prop
  saveRDS(true_weights, paste0(res_true_path, "_weights_true.RDs"))
  saveRDS(
    cohort_true$res_true_cohort,
    paste0(res_true_path, "_cohort_true.RDs")
  )

  # call sim_cohort to simulate cohort model for all scenarios
  start_ages <- list(
    start_age1 = c(40, 50, 60, 70, 80)
  )

  true_weights_g5 <- aggregate_weights(true_weights, start_ages[[1]])

  trial_data <- readRDS(file.path(pop_path, "trial_data.RDs"))
  trial_data_extra <- readRDS(file.path(pop_path, "trial_data_extra.RDs"))

  target_probs_treatment <- generate_probabilities(
    age_target,
    baseline_probs,
    control_type,
    treatment_type
  )
  target_probs_control <- plogis(control_outcome_logit(
    baseline_probs,
    age_target,
    control_type
  ))

  unadjusted_results <- adjusted_results <-
    linear_results <- unrestricted_spline_results <- monotonic_spline_results <- vector(
      "list",
      n_rep
    )

  unadjusted_path <- here(
    "results",
    "unadjusted",
    paste0(
      control_short,
      "_ctrl_",
      treatment_short,
      "_trt"
    )
  )
  adjusted_path <- here(
    "results",
    "adjusted",
    paste0(
      control_short,
      "_ctrl_",
      treatment_short,
      "_trt"
    )
  )
  linear_path <- here(
    "results",
    "linear_interaction",
    paste0(
      control_short,
      "_ctrl_",
      treatment_short,
      "_trt"
    )
  )
  unres_spline_path <- here(
    "results",
    "unrestricted_spline",
    paste0(
      control_short,
      "_ctrl_",
      treatment_short,
      "_trt"
    )
  )
  mon_spline_path <- here(
    "results",
    "monotonic_spline",
    paste0(
      control_short,
      "_ctrl_",
      treatment_short,
      "_trt"
    )
  )

  ncores <- 32
  cl <- makeCluster(ncores, type = "SOCK")
  registerDoSNOW(cl)

  unadjusted_results <- foreach(
    j = 1:n_rep,
    .errorhandling = "pass",
    .packages = c('tidyverse', 'rstanarm', 'tidybayes')
  ) %dopar%
    {
      unadjusted_model(trial_data[[j]], trial_data_extra[[j]])
    }

  unadjusted_results <- lapply(unadjusted_results, function(sublist) {
    multi_cohort_average(sublist, true_weights_g5)
  })

  unadjusted_limited <- lapply(unadjusted_results, function(x) {
    x$unadjusted_limited
  })
  unadjusted_extended <- lapply(unadjusted_results, function(x) {
    x$unadjusted_extended
  })
  saveRDS(unadjusted_limited, paste0(unadjusted_path, "_limited.RDs"))
  saveRDS(unadjusted_extended, paste0(unadjusted_path, "_extended.RDs"))
  rm(unadjusted_limited, unadjusted_extended, unadjusted_results)
  gc()

  adjusted_results <- foreach(
    j = 1:n_rep,
    .errorhandling = "pass",
    .packages = c('tidyverse', 'rstanarm', 'tidybayes')
  ) %dopar%
    {
      adjusted_model(
        trial_data[[j]],
        trial_data_extra[[j]],
        filename = plot_path,
        plot_outcome = (j == 1)
      )
    }

  adjusted_results <- lapply(adjusted_results, function(sublist) {
    multi_cohort_average(sublist, true_weights_g5)
  })

  adjusted_limited <- lapply(adjusted_results, function(x) x$adjusted_limited)
  adjusted_extended <- lapply(adjusted_results, function(x) x$adjusted_extended)
  saveRDS(adjusted_limited, paste0(adjusted_path, "_limted.RDs"))
  saveRDS(adjusted_extended, paste0(adjusted_path, "_extended.RDs"))
  rm(adjusted_limited, adjusted_extended, adjusted_results)
  gc()

  linear_results <- foreach(
    j = 1:n_rep,
    .errorhandling = "pass",
    .packages = c('tidyverse', 'rstanarm', 'tidybayes')
  ) %dopar%
    {
      linear_interaction_model(
        trial_data[[j]],
        trial_data_extra[[j]],
        filename = plot_path,
        plot_outcome = (j == 1)
      )
    }

  linear_results <- lapply(linear_results, function(sublist) {
    multi_cohort_average(sublist, true_weights_g5)
  })

  linear_limited <- lapply(linear_results, function(x) x$linear_limited)
  linear_extended <- lapply(linear_results, function(x) x$linear_extended)
  saveRDS(linear_limited, paste0(linear_path, "_limited.RDs"))
  saveRDS(linear_extended, paste0(linear_path, "_extended.RDs"))
  rm(linear_limited, linear_extended, linear_results)
  gc()

  unrestricted_spline_results <- foreach(
    j = 1:n_rep,
    .errorhandling = "pass",
    .packages = c(
      'dplyr',
      'cmdstanr',
      'tidyr',
      'ggplot2',
      'splines2',
      'arm',
      'tidybayes'
    )
  ) %dopar%
    {
      unrestricted_spline_model(
        mod_sep_rw,
        trial_data[[j]],
        trial_data_extra[[j]],
        filename = plot_path,
        plot_outcome = (j == 1)
      )
    }

  unrestricted_spline_results <- lapply(
    unrestricted_spline_results,
    function(sublist) {
      multi_cohort_average(sublist, true_weights_g5)
    }
  )

  unrestricted_spline_limited <- lapply(
    unrestricted_spline_results,
    function(x) x$unres_spline_limited
  )
  unrestricted_spline_extended <- lapply(
    unrestricted_spline_results,
    function(x) x$unres_spline_extended
  )
  saveRDS(
    unrestricted_spline_limited,
    paste0(unres_spline_path, "_limited.RDs")
  )
  saveRDS(
    unrestricted_spline_extended,
    paste0(unres_spline_path, "_extended.RDs")
  )
  rm(
    unrestricted_spline_limited,
    unrestricted_spline_extended,
    unrestricted_spline_results
  )
  gc()

  monotonic_spline_results <- foreach(
    j = 1:n_rep,
    .errorhandling = "pass",
    .packages = c(
      'dplyr',
      'cmdstanr',
      'tidyr',
      'ggplot2',
      'splines2',
      'arm',
      'tidybayes'
    )
  ) %dopar%
    {
      monotonic_spline_model(
        mod_sep_dir,
        trial_data[[j]],
        trial_data_extra[[j]],
        filename = plot_path,
        plot_outcome = (j == 1)
      )
    }

  monotonic_spline_results <- lapply(
    monotonic_spline_results,
    function(sublist) {
      multi_cohort_average(sublist, true_weights_g5)
    }
  )

  monotonic_spline_limited <- lapply(monotonic_spline_results, function(x) {
    x$mono_spline_limited
  })
  monotonic_spline_extended <- lapply(monotonic_spline_results, function(x) {
    x$mono_spline_extended
  })
  saveRDS(monotonic_spline_limited, paste0(mon_spline_path, "_limited.RDs"))
  saveRDS(monotonic_spline_extended, paste0(mon_spline_path, "_extended.RDs"))
  rm(
    monotonic_spline_limited,
    monotonic_spline_extended,
    monotonic_spline_results
  )
  gc()

  stopCluster(cl)
}
