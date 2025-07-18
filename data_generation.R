sample_age <- function(n, probs, age_grid, tol = 1e-8) {
  if (length(probs) != length(age_grid)) {
    stop("The length of 'probs' and 'age_grid' must be the same.")
  }
  if (abs(sum(probs) - 1) > tol) {
    stop("The probabilities 'probs' must sum to 1.")
  }
  return(sample(age_grid, size = n, replace = TRUE, prob = probs))
}

uniform_probs <- function(age_grid) {
  probs <- ifelse(age_grid >= 40 & age_grid <= 60, 1, 0)
  return(probs / sum(probs))
}

half_normal_probs <- function(age_grid, sd = 5) {
  probs <- ifelse(
    age_grid > 60 & age_grid <= 80,
    dnorm(age_grid, mean = 70, sd = sd),
    0
  )
  return(probs / sum(probs))
}

control_outcome_logit <- function(baseline_probs, age, type = "constant") {
  baseline_logit <- log(baseline_probs / (1 - baseline_probs)) #log(baseline_probs / (1 - baseline_probs)) # Lower bound

  type <- as.character(type)

  switch(
    type,
    "constant" = {
      return(rep(baseline_logit, length(age)))
    },
    "linear_increasing" = {
      return(baseline_logit + (age - 40) * 0.025)
    },
    "linear_decreasing" = {
      return(baseline_logit - (age - 40) * 0.025)
    },
    "non_linear_mon_increase" = {
      return(baseline_logit + (age - 40) * 0.025 + log(age - 40 + 1) * 0.1)
    },
    "non_linear_mon_decrease" = {
      return(baseline_logit - (age - 40) * 0.025 - log(age - 40 + 1) * 0.1)
    },
    "non_linear_non_mon" = {
      return(-0.002 * (age - 55)^2 + baseline_logit)
    }
  )
}

treatment_effect_logit <- function(age, type = "constant") {
  type <- as.character(type)

  switch(
    type,
    "constant" = {
      return(rep(log(0.5), length(age)))
    },
    "linear_increasing" = {
      return(-(age - 39) * 0.02) # Increasing effectiveness (negative coefficient)
    },
    "linear_decreasing" = {
      return(-(80 - age) * 0.02) # Decreasing effectiveness (positive coefficient)
    },
    "non_linear_mon_increase" = {
      return(-2 * log(1 + 0.1 * sqrt((age - 39) + 1))) # Increasing effectiveness
    },
    "non_linear_mon_decrease" = {
      return(-log(1.5 - 0.2 * sqrt((age - 39)) + 1)) # Decreasing effectiveness
    },
    "non_linear_non_mon" = {
      return(-0.001 * (age - 60)**2 - 0.2)
    }
  )
}

generate_probabilities <- function(
  age_grid,
  baseline_probs,
  control_type,
  treatment_type
) {
  # Control logit
  control_logit <- control_outcome_logit(baseline_probs, age_grid, control_type)

  # Treatment logit
  treatment_logit <- treatment_effect_logit(age_grid, treatment_type)

  # Combined logit
  combined_logit <- control_logit + treatment_logit

  # Convert to probability
  prob <- exp(combined_logit) / (1 + exp(combined_logit))

  # Rescale probabilities to fit within desired range

  return(prob)
}

gen_trial_data <- function(
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
) {
  age_grid <- 40:80

  control_data <- tibble(
    probs = probs_control,
    age = ages_control,
    trt = rep(0, n_trial)
  )
  control_data$outcome <- rbinom(n_trial, 1, control_data$probs)

  treatment_data <- tibble(
    probs = probs_treatment,
    age = ages_treatment,
    trt = rep(1, n_trial)
  )
  treatment_data$outcome <- rbinom(n_trial, 1, treatment_data$probs)

  trial_data <- rbind(control_data, treatment_data)
  trial_data$age_std <- as.vector(scale(trial_data$age))

  control_data_extra <- tibble(
    probs = probs_control_extra,
    age = extra_control,
    trt = rep(0, n_extra)
  )
  control_data_extra$outcome <- rbinom(n_extra, 1, control_data_extra$probs)
  control_data_extra <- rbind(control_data, control_data_extra)

  treatment_data_extra <- tibble(
    probs = probs_treatment_extra,
    age = extra_trt,
    trt = rep(1, n_extra)
  )
  treatment_data_extra$outcome <- rbinom(n_extra, 1, treatment_data_extra$probs)
  treatment_data_extra <- rbind(treatment_data, treatment_data_extra)

  trial_data_extra <- rbind(control_data_extra, treatment_data_extra)
  trial_data_extra$age_std <- as.vector(scale(trial_data_extra$age))

  return(list(trial_data = trial_data, trial_data_extra = trial_data_extra))
}
