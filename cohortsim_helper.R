posterior_transition_probs <- function(
  model,
  trt_status,
  age_varying = TRUE,
  standardize_age = FALSE
) {
  if (age_varying) {
    pred_data <- tibble(age = 40:80, trt = rep(trt_status, length(40:80)))
    if (standardize_age) {
      pred_data$age_std = scale(pred_data$age) %>% as.vector()
    }
    pred_result <- pred_data %>%
      add_epred_draws(., model) %>%
      ungroup() %>%
      select(age, .draw, .epred) %>%
      pivot_wider(names_from = age, values_from = .epred) %>%
      select(-.draw) %>%
      as.matrix()

    # Extend to 81-120 with the last column values
    last_col <- pred_result[, ncol(pred_result)]
    extended_cols <- matrix(
      rep(last_col, times = length(81:120)),
      nrow = n_simulations,
      byrow = F
    )

    final_result <- cbind(pred_result, extended_cols)
  } else {
    pred_result <- model %>%
      add_epred_draws(., newdata = tibble(trt = trt_status)) %>%
      ungroup() %>%
      select(.epred) %>%
      as.matrix()
    final_result <- pred_result # No need to extend if it's not age-varying
  }

  return(final_result)
}

posterior_spline_matrix <- function(model, fitted_basis, newx) {
  pred_result <- cmdstan_epred_sep(model, fitted_basis, newx)
  last_col <- pred_result[, ncol(pred_result)]
  extended_cols <- matrix(
    rep(last_col, times = length(81:120)),
    nrow = n_simulations,
    byrow = F
  )
  final_result <- cbind(pred_result, extended_cols)

  return(final_result)
}

sim_true_cohort <- function(control_type, treatment_type) {
  age_grid <- 40:80
  limited_probs <- uniform_probs(age_grid)
  extra_probs <- half_normal_probs(age_grid, sd = 5)
  combined_probs <- c(0.6 * limited_probs[1:21], 0.4 * extra_probs[22:41])

  true_cohort <- target_prob_vector(control_type, treatment_type, max_age = 120)
  true_cohort_ctrl <- matrix(true_cohort$control, nrow = 1)
  true_cohort_trt <- matrix(true_cohort$treatment, nrow = 1)

  utilities <- c(utility_healthy, utility_diseased, 0)
  cost_ctrl <- c(cost_healthy, cost_diseased, 0)
  cost_trt <- c(cost_healthy + treatment_cost_per_cycle, cost_diseased, 0)

  ncores <- 32
  cl <- makeCluster(ncores, type = "SOCK")
  clusterExport(
    cl,
    c("cohort_state_transition", "n_cycles", "n_target", "compute_delta")
  )
  start_ages <- 40:80
  extended_age_grid <- 40:120
  res_true_cohort_list <- parLapply(cl, start_ages, function(start_age) {
    age_idx_range = which(
      extended_age_grid >= start_age & extended_age_grid < (start_age + 40)
    )
    cohort_state_transition(
      true_cohort_ctrl,
      true_cohort_trt,
      utilities,
      cost_ctrl,
      cost_trt,
      n_cycles,
      age_idx_range
    )
  })

  stopCluster(cl)

  cohort_true_collapsed <- bind_rows(res_true_cohort_list)
  res_true_cohort <- cohort_true_collapsed %>%
    summarise(
      delta_QALYs = weighted.mean(
        cohort_true_collapsed$delta_QALYs,
        combined_probs
      ),
      delta_costs = weighted.mean(
        cohort_true_collapsed$delta_costs,
        combined_probs
      )
    )

  return(list(res_true_cohort = res_true_cohort, target_prop = combined_probs))
}


run_multi_cohort_sequential <- function(
  start_ages,
  final_tp_ctrl,
  final_tp_trt,
  n_cycles,
  utilities,
  cost_ctrl,
  cost_trt
) {
  final_tp_ctrl <- final_tp_ctrl
  final_tp_trt <- final_tp_trt

  # Using lapply for sequential processing
  multi_cohort_stm <- lapply(start_ages, function(start_age) {
    cohort_list <- list()
    age_idx_range <- which(
      extended_age_grid >= start_age & extended_age_grid < (start_age + 40)
    )
    for (i in 1:nrow(final_tp_ctrl)) {
      ctrl_probs <- final_tp_ctrl[i, ]
      trt_probs <- final_tp_trt[i, ]
      cohort_list[[i]] <- cohort_state_transition(
        ctrl_probs,
        trt_probs,
        utilities = utilities,
        time_horizon = n_cycles,
        cost_ctrl = cost_ctrl,
        cost_trt = cost_trt,
        age_idx_range = age_idx_range
      )
    }

    return(cohort_list)
  })

  return(multi_cohort_stm)
}


target_prob_vector <- function(control_type, treatment_type, max_age) {
  age_grid_var = 40:80
  age_grid_const = 81:max_age

  # Generate probabilities using your specific functions for control and treatment
  control_probs_var = plogis(control_outcome_logit(
    baseline_probs,
    age_grid_var,
    control_type
  ))
  treatment_probs_var = generate_probabilities(
    age_grid_var,
    baseline_probs,
    control_type,
    treatment_type
  )

  # For ages beyond 80, use the last value
  control_probs_const = rep(
    control_probs_var[length(control_probs_var)],
    length(age_grid_const)
  )
  treatment_probs_const = rep(
    treatment_probs_var[length(treatment_probs_var)],
    length(age_grid_const)
  )

  # Combine
  control_probs = c(control_probs_var, control_probs_const)
  treatment_probs = c(treatment_probs_var, treatment_probs_const)

  return(list(control = control_probs, treatment = treatment_probs))
}


average_matrices <- function(mat_list) {
  avg_mat_ctrl <- Reduce("+", lapply(mat_list, `[[`, 1)) / length(mat_list)
  avg_mat_trt <- Reduce("+", lapply(mat_list, `[[`, 2)) / length(mat_list)
  list(trace_ctrl = avg_mat_ctrl, trace_trt = avg_mat_trt)
}

process_multi_cohort_stm <- function(multi_cohort_stm) {
  n_cohorts <- length(multi_cohort_stm)
  multi_cohort_trace <- vector("list", n_cohorts)

  # Extract trace matrices for each cohort and store in a list
  for (i in 1:n_cohorts) {
    cohort_list <- multi_cohort_stm[[i]]
    mat_list <- lapply(cohort_list, function(x) {
      list(x$trace_prop_ctrl, x$trace_prop_trt)
    })

    # Use average_matrices function to get averaged trace matrices
    multi_cohort_trace[[i]] <- average_matrices(mat_list)
  }

  return(multi_cohort_trace)
}


multi_cohort_average <- function(result, weights) {
  # Use lapply to loop over each 'method' in 'result'
  lapply(result, function(single_method_results) {
    # Use Reduce to sum across all cohort lists
    weighted_delta_QALYs <- Reduce(
      '+',
      lapply(seq_along(single_method_results), function(i) {
        sapply(single_method_results[[i]], function(x) x$delta_QALYs) *
          weights[i]
      })
    )
    weighted_delta_costs <- Reduce(
      '+',
      lapply(seq_along(single_method_results), function(i) {
        sapply(single_method_results[[i]], function(x) x$delta_costs) *
          weights[i]
      })
    )

    # Return list with weighted delta_QALYs and delta_costs
    list(delta_QALYs = weighted_delta_QALYs, delta_costs = weighted_delta_costs)
  })
}


aggregate_weights <- function(true_weights, starting_ages) {
  if (any(diff(starting_ages) == 10)) {
    # For increments of 10
    bins <- list(c(1:6), c(7:16), c(17:26), c(27:36), c(37:41))
  } else if (any(diff(starting_ages) == 5)) {
    # For increments of 5
    bins <- lapply(seq(1, 36, by = 5), function(x) seq(x, x + 4))
    bins <- append(bins, 41)
  }

  aggregated_weights <- sapply(bins, function(bin) {
    sum(true_weights[bin])
  })

  return(aggregated_weights)
}


compute_delta <- function(trace_ctrl, trace_trt, payoff_ctrl, payoff_trt) {
  # Initialize discount rate and vector
  discount_rate <- 0.03
  discount_vector <- (1 / (1 + discount_rate))^(1:n_cycles)

  control_vec <- cumsum((payoff_ctrl %*% t(trace_ctrl[-1, ])) * discount_vector)
  treatment_vec <- cumsum((payoff_trt %*% t(trace_trt[-1, ])) * discount_vector)

  delta_vec <- treatment_vec - control_vec
  return(delta = delta_vec)
}


extract_plot_data <- function(plot_list) {
  plot_list$plot$data
}
