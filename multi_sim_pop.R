library(tidyverse)


source(here("data_generation.R"))


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

# no. of MC replications
n_rep <- 1000


for (i in 1:nrow(combinations)) {
  control_type <- combinations$control[i]
  treatment_type <- combinations$treatment[i]

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
  pop_path <- paste0(
    "pop_data/",
    control_short,
    "_ctrl_",
    treatment_short,
    "_trt",
    sep = ""
  )
  dir.create(pop_path, showWarnings = F)

  probs_control <- plogis(control_outcome_logit(
    baseline_probs,
    ages_control,
    control_type
  ))
  probs_treatment <- generate_probabilities(
    ages_treatment,
    baseline_probs,
    control_type,
    treatment_type
  )
  probs_control_extra <- plogis(control_outcome_logit(
    baseline_probs,
    extra_control,
    control_type
  ))
  probs_treatment_extra <- generate_probabilities(
    extra_trt,
    baseline_probs,
    control_type,
    treatment_type
  )

  # Generate the population data

  pop_data <- replicate(
    n = n_rep,
    expr = gen_trial_data(
      n_trial,
      n_extra,
      probs_control,
      probs_treatment,
      probs_control_extra,
      probs_treatment_extra,
      ages_control,
      ages_treatment,
      extra_control,
      extra_trt
    ),
    simplify = F
  )

  trial_data_list <- lapply(pop_data, function(x) x$trial_data)
  extra_data_list <- lapply(pop_data, function(x) x$trial_data_extra)

  # Save the population data
  saveRDS(trial_data_list, file = paste0(pop_path, "/trial_data.RDs"))
  saveRDS(extra_data_list, file = paste0(pop_path, "/trial_data_extra.RDs"))
}
