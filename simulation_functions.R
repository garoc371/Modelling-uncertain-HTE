cmdstan_epred_sep <- function(fit, fitted_basis, newx) {
  new_basis <- predict(fitted_basis, newx)
  nbasis <- dim(new_basis)[2]
  beta_0 <- fit$draws("Intercept", format = "df")[, 1] %>% as.matrix()
  s <- fit$draws("s_0", format = "draws_df")[, 1:nbasis] %>% as.matrix()
  S <- s %*% t(new_basis)

  y_linpred <- sweep(S, 1, beta_0, "+")
  y_pred <- arm::invlogit(y_linpred)

  attributes(y_pred)$dimnames[[2]] <- NULL

  return(y_pred)
}


outcome_data <- function(true_outcome_matrix, treatment_label, age_vector) {
  as.data.frame(true_outcome_matrix) %>%
    mutate(age = age_vector) %>%
    pivot_longer(
      cols = -age,
      names_to = "sample",
      values_to = "probability"
    ) %>%
    mutate(treatment = treatment_label, type = "True")
}

posterior_summary <- function(
  pred,
  treatment_label,
  age_vector,
  ci_levels = c(0.5, 0.8, 0.95)
) {
  pred <- as.data.frame(t(pred)) %>%
    mutate(age = age_vector) %>%
    pivot_longer(
      cols = -age,
      names_to = "sample",
      values_to = "probability"
    ) %>%
    mutate(treatment = treatment_label) %>%
    group_by(age, treatment) %>%
    median_qi(probability, .width = ci_levels)
}

plot_outcome_surface <- function(data, true_data, trans_logit = FALSE) {
  if (trans_logit) {
    data$probability = logit(data$probability)
    true_data$probability = logit(true_data$probability)
  }

  ggplot() +
    geom_lineribbon(
      data = data,
      aes(x = as.numeric(age), y = probability, ymin = .lower, ymax = .upper),
      alpha = 0.6
    ) +
    geom_line(
      data = true_data,
      aes(x = as.numeric(age), y = probability, color = type, linewidth = type)
    ) +
    scale_fill_brewer(palette = "Blues") +
    scale_color_manual(values = c("True" = "red")) +
    scale_linewidth_manual(values = c("True" = 1.5)) +
    facet_wrap(~treatment) +
    theme_minimal() +
    labs(x = "Age", y = expression("P"))
}

plot_ce <- function(cep_data) {
  ce_plot <- ggplot(
    cep_data,
    aes(x = delta_QALYs, y = delta_costs, color = source)
  ) +
    geom_point(alpha = 0.4, size = 1) +
    scale_color_manual(values = c("blue", "red")) +
    geom_point(
      data = cep_data[cep_data$source == "true_value", ],
      color = "red",
      size = 4,
      shape = 18
    ) +
    theme_minimal() +
    theme(legend.position = "bottom") +
    labs(
      x = "Incremental QALYs",
      y = "Incremental Costs",
      color = "Source"
    )

  return(ce_plot)
}

knots_quantile <- function(data, increment) {
  knots <- quantile(data, probs = seq(.05, .95, increment))
  len <- length(knots)
  bknots <- range(data)
  if (knots[1] == bknots[1]) {
    knots[1] <- knots[1] + 0.5
  }
  if (knots[len] == bknots[2]) {
    knots[len] <- knots[len] - 0.5
  }
  return(knots)
}

cmdstan_fit_plot <- function(
  mod,
  increment,
  target_age,
  degree,
  ctrl_dat,
  ctrl_dat_ex,
  trt_dat,
  trt_dat_ex,
  target_probs_control,
  target_probs_treatment,
  ispline = TRUE,
  plot_outcome = F,
  file_name
) {
  knots1 <- knots_quantile(ctrl_dat$age, increment)
  knots2 <- knots_quantile(trt_dat$age, increment)
  knots3 <- knots_quantile(ctrl_dat_ex$age, increment)
  knots4 <- knots_quantile(trt_dat_ex$age, increment)

  mspline_degree <- degree

  pred_age_target <- target_age
  bknots <- range(pred_age_target)
  bknots[2] <- bknots[2] + 0.5

  if (ispline) {
    ctrl_spline_basis <- iSpline(
      ctrl_dat$age,
      knots = knots1,
      Boundary.knots = bknots,
      degree = mspline_degree,
      intercept = T
    )
    trt_spline_basis <- iSpline(
      trt_dat$age,
      knots = knots2,
      Boundary.knots = bknots,
      degree = mspline_degree,
      intercept = T
    )

    ex_ctrl_spline_basis <- iSpline(
      ctrl_dat_ex$age,
      knots = knots3,
      Boundary.knots = bknots,
      degree = mspline_degree,
      intercept = T
    )
    ex_trt_spline_basis <- iSpline(
      trt_dat_ex$age,
      knots = knots4,
      Boundary.knots = bknots,
      degree = mspline_degree,
      intercept = T
    )
  } else {
    ctrl_spline_basis <- bSpline(
      ctrl_dat$age,
      knots = knots1[-1],
      Boundary.knots = bknots,
      degree = mspline_degree,
      intercept = T
    )
    trt_spline_basis <- bSpline(
      trt_dat$age,
      knots = knots2[-1],
      Boundary.knots = bknots,
      degree = mspline_degree,
      intercept = T
    )

    ex_ctrl_spline_basis <- bSpline(
      ctrl_dat_ex$age,
      knots = knots3[-1],
      Boundary.knots = bknots,
      degree = mspline_degree,
      intercept = T
    )
    ex_trt_spline_basis <- bSpline(
      trt_dat_ex$age,
      knots = knots4[-1],
      Boundary.knots = bknots,
      degree = mspline_degree,
      intercept = T
    )
  }
  nbasis <- dim(ctrl_spline_basis)[2]
  # mod <- cmdstan_model("rw_spline.stan")
  ctrl_data_list <- list(
    N = nrow(ctrl_dat),
    Y = ctrl_dat$outcome,
    m = nbasis,
    S = ctrl_spline_basis
  )

  ex_ctrl_data_list <- list(
    N = nrow(ctrl_dat_ex),
    Y = ctrl_dat_ex$outcome,
    m = nbasis,
    S = ex_ctrl_spline_basis
  )

  trt_data_list <- list(
    N = nrow(trt_dat),
    Y = trt_dat$outcome,
    m = nbasis,
    S = trt_spline_basis
  )

  ex_trt_data_list <- list(
    N = nrow(trt_dat_ex),
    Y = trt_dat_ex$outcome,
    m = nbasis,
    S = ex_trt_spline_basis
  )

  fit_ctrl <- mod$sample(
    data = ctrl_data_list,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.95
  )

  fit_ctrl_ex <- mod$sample(
    data = ex_ctrl_data_list,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.95
  )

  fit_trt <- mod$sample(
    data = trt_data_list,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.95
  )

  fit_trt_ex <- mod$sample(
    data = ex_trt_data_list,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.95
  )

  ctrl_pred_basis <- predict(ctrl_spline_basis, newx = pred_age_target)
  ctrl_pred_basis_ex <- predict(ex_ctrl_spline_basis, newx = pred_age_target)
  trt_pred_basis <- predict(trt_spline_basis, newx = pred_age_target)
  trt_pred_basis_ex <- predict(ex_trt_spline_basis, newx = pred_age_target)

  if (plot_outcome) {
    pred_y0s <- cmdstan_epred_sep(fit_ctrl, ctrl_spline_basis, pred_age_target)
    pred_y1s <- cmdstan_epred_sep(fit_trt, trt_spline_basis, pred_age_target)

    pred_diff <- pred_y1s - pred_y0s

    posterior_control_sep <- posterior_summary(
      pred_y0s,
      "Control",
      pred_age_target
    )
    posterior_trt_sep <- posterior_summary(
      pred_y1s,
      "Treatment",
      pred_age_target
    )
    posterior_effect <- posterior_summary(
      pred_diff,
      "Treatment_effect",
      pred_age_target
    )

    rm(pred_y0s, pred_y1s)
    gc()

    pred_y0s_ex <- cmdstan_epred_sep(
      fit_ctrl_ex,
      ex_ctrl_spline_basis,
      pred_age_target
    )
    pred_y1s_ex <- cmdstan_epred_sep(
      fit_trt_ex,
      ex_trt_spline_basis,
      pred_age_target
    )
    pred_diff_ex <- pred_y1s_ex - pred_y0s_ex

    posterior_control_exsep <- posterior_summary(
      pred_y0s_ex,
      "Control",
      pred_age_target
    )
    posterior_trt_exsep <- posterior_summary(
      pred_y1s_ex,
      "Treatment",
      pred_age_target
    )
    posterior_effect_ex <- posterior_summary(
      pred_diff_ex,
      "Treatment_effect",
      pred_age_target
    )

    rm(pred_y0s_ex, pred_y1s_ex)

    true_outcome_control <- outcome_data(
      target_probs_control,
      "Control",
      age_target
    ) %>%
      mutate(sample = "True")
    true_outcome_trt <- outcome_data(
      target_probs_treatment,
      "Treatment",
      age_target
    ) %>%
      mutate(sample = "True")

    true_effect <- target_probs_treatment - target_probs_control
    true_effect <- outcome_data(true_effect, "Treatment_effect", age_target) %>%
      mutate(sample = "True")
    # true_effect_logit <- logit(target_probs_treatment)-logit(target_probs_control)
    # true_effect_logit <- outcome_data(true_effect_logit, "Treatment_effect(logit)", age_target) %>% mutate(sample = "True")

    sep_spline_combined_data <- bind_rows(
      posterior_control_sep,
      posterior_trt_sep
    )

    true_data <- bind_rows(true_outcome_control, true_outcome_trt)

    s1 <- plot_outcome_surface(sep_spline_combined_data, true_data = true_data)
    rm(sep_spline_combined_data)
    gc()

    ex_sep_ispline_combined_data <- bind_rows(
      posterior_control_exsep,
      posterior_trt_exsep
    )
    s2 <- plot_outcome_surface(
      ex_sep_ispline_combined_data,
      true_data = true_data
    )
    rm(ex_sep_ispline_combined_data, true_data)
    gc()

    s3 <- plot_outcome_surface(posterior_effect, true_data = true_effect)
    s4 <- plot_outcome_surface(posterior_effect_ex, true_data = true_effect)

    plots_all <- list(
      limited_outcome = s1,
      extended_outcome = s2,
      limited_effect = s3,
      extended_effect = s4
    )

    saveRDS(plots_all, file_name)
  }

  return(list(
    fit_ctrl = fit_ctrl,
    fit_trt = fit_trt,
    fit_ctrl_ex = fit_ctrl_ex,
    fit_trt_ex = fit_trt_ex,
    ctrl_spline_basis = ctrl_spline_basis,
    trt_spline_basis = trt_spline_basis,
    ex_ctrl_spline_basis = ex_ctrl_spline_basis,
    ex_trt_spline_basis = ex_trt_spline_basis
  ))
}

cohort_state_transition <- function(
  ctrl_probs,
  trt_probs,
  utilities,
  cost_ctrl,
  cost_trt,
  time_horizon,
  age_idx_range
) {
  if (length(ctrl_probs) == 1) {
    ctrl_probs <- rep(ctrl_probs, length(extended_age_grid))
  }
  if (length(trt_probs) == 1) {
    trt_probs <- rep(trt_probs, length(extended_age_grid))
  }

  ctrl_probs <- matrix(ctrl_probs, nrow = 1)
  trt_probs <- matrix(trt_probs, nrow = 1)
  ctrl_probs <- ctrl_probs[, age_idx_range]
  trt_probs <- trt_probs[, age_idx_range]

  trans_array_trt <- trans_array_ctrl <- array(0, dim = c(3, 3, time_horizon))
  for (t in 1:time_horizon) {
    p_ctrl <- ctrl_probs[t]
    p_trt <- trt_probs[t]

    trans_mat_ctrl <- matrix(
      c(1 - 0.5 * p_ctrl, 0.5 * p_ctrl, 0, 0, 1 - p_ctrl, p_ctrl, 0, 0, 1),
      nrow = 3,
      ncol = 3,
      byrow = TRUE
    )

    trans_mat_trt <- matrix(
      c(1 - 0.5 * p_trt, 0.5 * p_trt, 0, 0, 1 - p_trt, p_trt, 0, 0, 1),
      nrow = 3,
      ncol = 3,
      byrow = TRUE
    )

    trans_array_trt[,, t] <- trans_mat_trt
    trans_array_ctrl[,, t] <- trans_mat_ctrl
  }

  trace_mat_ctrl <- trace_mat_trt <- matrix(
    NA,
    nrow = 3,
    ncol = time_horizon + 1
  )
  trace_mat_ctrl[, 1] <- trace_mat_trt[, 1] <- c(n_target, 0, 0)

  delta_qaly_cohort <- delta_cost_cohort <- list()

  for (s in 1:time_horizon) {
    trace_mat_ctrl[, s + 1] <- trace_mat_ctrl[, s] %*% trans_array_ctrl[,, s]
    trace_mat_trt[, s + 1] <- trace_mat_trt[, s] %*% trans_array_trt[,, s]
  }

  discount_rate <- 0.03
  discount_vector <- (1 / (1 + discount_rate))^(1:(time_horizon))

  discounted_trace_ctrl <- sweep(
    trace_mat_ctrl[, -1],
    2,
    discount_vector,
    FUN = "*"
  )
  discounted_trace_trt <- sweep(
    trace_mat_trt[, -1],
    2,
    discount_vector,
    FUN = "*"
  )

  # Apply discounting after state transition but before summing across time
  cum_qalys_ctrl <- sum(utilities %*% discounted_trace_ctrl)
  cum_costs_ctrl <- sum(cost_ctrl %*% discounted_trace_ctrl)

  cum_qalys_trt <- sum(utilities %*% discounted_trace_trt)
  cum_costs_trt <- sum(cost_trt %*% discounted_trace_trt)

  delta_qaly_cohort <- cum_qalys_trt - cum_qalys_ctrl
  delta_cost_cohort <- cum_costs_trt - cum_costs_ctrl

  # trace_prop_ctrl <- t(trace_mat_ctrl)
  # trace_prop_trt <- t(trace_mat_trt)
  #
  # delta_trace_qaly <- compute_delta(trace_ctrl = trace_prop_ctrl, trace_trt = trace_prop_trt,
  #                                   payoff_ctrl = utilities, payoff_trt = utilities)
  # delta_trace_costs <- compute_delta(trace_ctrl = trace_prop_ctrl, trace_trt = trace_prop_trt,
  #                                    payoff_ctrl = cost_ctrl, payoff_trt = cost_trt)

  return(list(delta_QALYs = delta_qaly_cohort, delta_costs = delta_cost_cohort))
}


vis_outcome_effect <- function(
  true_outcome_control,
  true_outcome_trt,
  bayesian_model,
  bayesian_model_ex,
  age_target
) {
  # Calculate predictions and posterior data

  y0_pred <- posterior_epred(
    bayesian_model,
    newdata = data.frame(
      trt = rep(0, length(age_target)),
      age_std = scale(age_target)
    )
  )
  y1_pred <- posterior_epred(
    bayesian_model,
    newdata = data.frame(
      trt = rep(1, length(age_target)),
      age_std = scale(age_target)
    )
  )

  pred_diff <- y1_pred - y0_pred
  pred_diff_logit <- logit(y1_pred) - logit(y0_pred)

  posterior_control <- posterior_summary(y0_pred, "Control", age_target)
  posterior_trt <- posterior_summary(y1_pred, "Treatment", age_target)
  posterior_effect <- posterior_summary(
    pred_diff,
    "Treatment_effect",
    age_target
  )
  posterior_effect_logit <- posterior_summary(
    pred_diff_logit,
    "Treatment_effect(logit)",
    age_target
  )

  # Calculate predictions and posterior data for extended data
  y0_pred_ex <- posterior_epred(
    bayesian_model_ex,
    newdata = data.frame(
      trt = rep(0, length(age_target)),
      age_std = scale(age_target)
    )
  )
  y1_pred_ex <- posterior_epred(
    bayesian_model_ex,
    newdata = data.frame(
      trt = rep(1, length(age_target)),
      age_std = scale(age_target)
    )
  )
  pred_diff_ex <- y1_pred_ex - y0_pred_ex
  pred_diff_ex_logit <- logit(y1_pred_ex) - logit(y0_pred_ex)

  posterior_control_ex <- posterior_summary(y0_pred_ex, "Control", age_target)
  posterior_trt_ex <- posterior_summary(y1_pred_ex, "Treatment", age_target)
  posterior_effect_ex <- posterior_summary(
    pred_diff_ex,
    "Treatment_effect",
    age_target
  )
  posterior_effect_logit_ex <- posterior_summary(
    pred_diff_ex_logit,
    "Treatment_effect(logit)",
    age_target
  )

  rm(
    y0_pred,
    y0_pred_ex,
    pred_diff_ex,
    pred_diff,
    pred_diff_ex_logit,
    pred_diff_logit,
    y1_pred,
    y1_pred_ex
  )

  # Generate plots
  true_data <- bind_rows(true_outcome_control, true_outcome_trt)
  combined_data <- bind_rows(posterior_control, posterior_trt)

  m1 <- plot_outcome_surface(combined_data, true_data)

  true_data_ex <- bind_rows(true_outcome_control, true_outcome_trt)
  ex_combined_data <- bind_rows(posterior_control_ex, posterior_trt_ex)

  m2 <- plot_outcome_surface(ex_combined_data, true_data_ex)

  # Return plots
  list(limited_outcome = m1, extended_outcome = m2)
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

adjusted_model <- function(
  trial_data,
  trial_data_extra,
  filename,
  plot_outcome = F
) {
  adjusted_model <- stan_glm(
    outcome ~ trt + age_std,
    data = trial_data,
    family = binomial(link = "logit"),
    prior = normal(0, 2.5, autoscale = T),
    prior_intercept = normal(0, 2.5, autoscale = T),
    chains = 4,
    iter = 2000
  )

  adjusted_model_ex <- stan_glm(
    outcome ~ trt + age_std,
    data = trial_data_extra,
    family = binomial(link = "logit"),
    prior = normal(0, 2.5, autoscale = T),
    prior_intercept = normal(0, 2.5, autoscale = T),
    chains = 4,
    iter = 2000
  )

  if (plot_outcome) {
    n_samples = 50

    true_outcome_control <- outcome_data(
      target_probs_control,
      "Control",
      age_target
    ) %>%
      mutate(sample = "True")
    true_outcome_trt <- outcome_data(
      target_probs_treatment,
      "Treatment",
      age_target
    ) %>%
      mutate(sample = "True")

    adjusted_plots <- vis_outcome_effect(
      true_outcome_control,
      true_outcome_trt,
      adjusted_model,
      adjusted_model_ex,
      age_grid
    )

    saveRDS(adjusted_plots, paste0(filename, "_adjusted.RDs"))
    rm(adjusted_plots)
    gc()
  }

  adjusted_tp_ctrl <- posterior_transition_probs(
    adjusted_model,
    trt_status = 0,
    age_varying = T,
    standardize_age = T
  )
  adjusted_tp_trt <- posterior_transition_probs(
    adjusted_model,
    trt_status = 1,
    age_varying = T,
    standardize_age = T
  )

  res_adjusted_multi_g5 <- run_multi_cohort_sequential(
    start_ages = start_ages[[1]],
    adjusted_tp_ctrl,
    adjusted_tp_trt,
    n_cycles = n_cycles,
    utilities = utilities,
    cost_ctrl = cost_ctrl,
    cost_trt = cost_trt
  )

  adjusted_tp_ctrl_ex <- posterior_transition_probs(
    adjusted_model_ex,
    trt_status = 0,
    age_varying = T,
    standardize_age = T
  )
  adjusted_tp_trt_ex <- posterior_transition_probs(
    adjusted_model_ex,
    trt_status = 1,
    age_varying = T,
    standardize_age = T
  )

  res_adjusted_multi_ex_g5 <- run_multi_cohort_sequential(
    start_ages = start_ages[[1]],
    adjusted_tp_ctrl_ex,
    adjusted_tp_trt_ex,
    n_cycles = n_cycles,
    utilities = utilities,
    cost_ctrl = cost_ctrl,
    cost_trt = cost_trt
  )

  return(list(
    adjusted_limited = res_adjusted_multi_g5,
    adjusted_extended = res_adjusted_multi_ex_g5
  ))
}

linear_interaction_model <- function(
  trial_data,
  trial_data_extra,
  filename,
  plot_outcome = F
) {
  bayesian_model <- stan_glm(
    outcome ~ trt * age_std,
    data = trial_data,
    family = binomial(link = "logit"),
    prior = normal(0, 2.5, autoscale = T),
    prior_intercept = normal(0, 2.5, autoscale = T),
    chains = 4,
    iter = 2000
  )

  bayesian_model_ex <- stan_glm(
    outcome ~ trt * age_std,
    data = trial_data_extra,
    family = binomial(link = "logit"),
    prior = normal(0, 2.5, autoscale = T),
    prior_intercept = normal(0, 2.5, autoscale = T),
    chains = 4,
    iter = 2000
  )

  if (plot_outcome) {
    n_samples = 50

    true_outcome_control <- outcome_data(
      target_probs_control,
      "Control",
      age_target
    ) %>%
      mutate(sample = "True")
    true_outcome_trt <- outcome_data(
      target_probs_treatment,
      "Treatment",
      age_target
    ) %>%
      mutate(sample = "True")

    linear_plots <- vis_outcome_effect(
      true_outcome_control,
      true_outcome_trt,
      bayesian_model,
      bayesian_model_ex,
      age_grid
    )
    saveRDS(linear_plots, paste0(filename, "_linear.RDs"))
    rm(linear_plots)
    gc()
  }

  linear_tp_ctrl <- posterior_transition_probs(
    bayesian_model,
    trt_status = 0,
    age_varying = T,
    standardize_age = T
  )
  linear_tp_trt <- posterior_transition_probs(
    bayesian_model,
    trt_status = 1,
    age_varying = T,
    standardize_age = T
  )

  res_linear_multi_g5 <- run_multi_cohort_sequential(
    start_ages = start_ages[[1]],
    linear_tp_ctrl,
    linear_tp_trt,
    n_cycles = n_cycles,
    utilities = utilities,
    cost_ctrl = cost_ctrl,
    cost_trt = cost_trt
  )

  linear_tp_ctrl_ex <- posterior_transition_probs(
    bayesian_model_ex,
    trt_status = 0,
    age_varying = T,
    standardize_age = T
  )
  linear_tp_trt_ex <- posterior_transition_probs(
    bayesian_model_ex,
    trt_status = 1,
    age_varying = T,
    standardize_age = T
  )

  res_linear_multi_ex_g5 <- run_multi_cohort_sequential(
    start_ages = start_ages[[1]],
    linear_tp_ctrl_ex,
    linear_tp_trt_ex,
    n_cycles = n_cycles,
    utilities = utilities,
    cost_ctrl = cost_ctrl,
    cost_trt = cost_trt
  )

  return(list(
    linear_limited = res_linear_multi_g5,
    linear_extended = res_linear_multi_ex_g5
  ))
}

unrestricted_spline_model <- function(
  stan_model,
  trial_data,
  trial_data_extra,
  filename,
  plot_outcome = F
) {
  file_name <- paste0(filename, "_unrestricted_spline.RDs")

  control_data <- trial_data %>% filter(trt == 0)
  treatment_data <- trial_data %>% filter(trt == 1)
  control_data_extra <- trial_data_extra %>% filter(trt == 0)
  treatment_data_extra <- trial_data_extra %>% filter(trt == 1)
  bspline_rw <- cmdstan_fit_plot(
    mod = stan_model,
    increment = 0.15,
    target_age = age_grid,
    degree = 3,
    ctrl_dat = control_data,
    ctrl_dat_ex = control_data_extra,
    trt_dat = treatment_data,
    trt_dat_ex = treatment_data_extra,
    target_probs_control = target_probs_control,
    target_probs_treatment = target_probs_treatment,
    ispline = F,
    plot_outcome = plot_outcome,
    file_name = file_name
  )

  bspline_tp_ctrl <- posterior_spline_matrix(
    bspline_rw$fit_ctrl,
    bspline_rw$ctrl_spline_basis,
    newx = age_grid
  )
  bspline_tp_ctrl_ex <- posterior_spline_matrix(
    bspline_rw$fit_ctrl_ex,
    bspline_rw$ex_ctrl_spline_basis,
    newx = age_grid
  )
  bspline_tp_trt <- posterior_spline_matrix(
    bspline_rw$fit_trt,
    bspline_rw$trt_spline_basis,
    newx = age_grid
  )
  bspline_tp_trt_ex <- posterior_spline_matrix(
    bspline_rw$fit_trt_ex,
    bspline_rw$ex_trt_spline_basis,
    newx = age_grid
  )

  res_bspline_multi_g5 <- run_multi_cohort_sequential(
    start_ages = start_ages[[1]],
    bspline_tp_ctrl,
    bspline_tp_trt,
    n_cycles = n_cycles,
    utilities = utilities,
    cost_ctrl = cost_ctrl,
    cost_trt = cost_trt
  )

  res_bspline_multi_ex_g5 <- run_multi_cohort_sequential(
    start_ages = start_ages[[1]],
    bspline_tp_ctrl_ex,
    bspline_tp_trt_ex,
    n_cycles = n_cycles,
    utilities = utilities,
    cost_ctrl = cost_ctrl,
    cost_trt = cost_trt
  )

  return(list(
    unres_spline_limited = res_bspline_multi_g5,
    unres_spline_extended = res_bspline_multi_ex_g5
  ))
}

unadjusted_model <- function(trial_data, trial_data_extra) {
  unadjusted_model <- stan_glm(
    outcome ~ trt,
    data = trial_data,
    family = binomial(link = "logit"),
    prior = normal(0, 2.5, autoscale = T),
    prior_intercept = normal(0, 2.5, autoscale = T),
    chains = 4,
    iter = 2000
  )

  unadjusted_model_ex <- stan_glm(
    outcome ~ trt,
    data = trial_data_extra,
    family = binomial(link = "logit"),
    prior = normal(0, 2.5, autoscale = T),
    prior_intercept = normal(0, 2.5, autoscale = T),
    chains = 4,
    iter = 2000
  )

  # calculate the transition probability for unadjusted model
  # and run the mutli-cohort simulation loop

  unadjusted_tp_ctrl <- posterior_transition_probs(
    unadjusted_model,
    trt_status = 0,
    age_varying = F
  )
  unadjusted_tp_trt <- posterior_transition_probs(
    unadjusted_model,
    trt_status = 1,
    age_varying = F
  )

  res_unadjusted_multi_g5 <- run_multi_cohort_sequential(
    start_ages = start_ages[[1]],
    unadjusted_tp_ctrl,
    unadjusted_tp_trt,
    n_cycles = n_cycles,
    utilities = utilities,
    cost_ctrl = cost_ctrl,
    cost_trt = cost_trt
  )

  # res_unadjusted_multi_g10 <- run_multi_cohort_parallel(start_ages = start_ages[[2]], unadjusted_tp_ctrl,
  #                                                      unadjusted_tp_trt, n_cycles = n_cycles, utilities = utilities,
  #                                                      cost_ctrl = cost_ctrl, cost_trt = cost_trt)

  gc()

  # res_unadjusted_multi_all <- run_multi_cohort_parallel(start_ages = start_ages[[3]], unadjusted_tp_ctrl,
  #                                                      unadjusted_tp_trt, n_cycles = n_cycles, utilities = utilities,
  #                                                      cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  #
  # gc()

  unadjusted_tp_ctrl_ex <- posterior_transition_probs(
    unadjusted_model_ex,
    trt_status = 0,
    age_varying = F
  )
  unadjusted_tp_trt_ex <- posterior_transition_probs(
    unadjusted_model_ex,
    trt_status = 1,
    age_varying = F
  )

  res_unadjusted_multi_ex_g5 <- run_multi_cohort_sequential(
    start_ages = start_ages[[1]],
    unadjusted_tp_ctrl_ex,
    unadjusted_tp_trt_ex,
    n_cycles = n_cycles,
    utilities = utilities,
    cost_ctrl = cost_ctrl,
    cost_trt = cost_trt
  )
  gc()

  return(list(
    unadjusted_limited = res_unadjusted_multi_g5,
    unadjusted_extended = res_unadjusted_multi_ex_g5
  ))
}

monotonic_spline_model <- function(
  stan_model,
  trial_data,
  trial_data_extra,
  filename,
  plot_outcome = F
) {
  filename <- paste0(filename, "_monotonic_spline.RDs")

  control_data <- trial_data %>% filter(trt == 0)
  treatment_data <- trial_data %>% filter(trt == 1)
  control_data_extra <- trial_data_extra %>% filter(trt == 0)
  treatment_data_extra <- trial_data_extra %>% filter(trt == 1)

  smol_d3 <- cmdstan_fit_plot(
    mod = stan_model,
    increment = 0.1,
    target_age = age_grid,
    degree = 3,
    ctrl_dat = control_data,
    ctrl_dat_ex = control_data_extra,
    trt_dat = treatment_data,
    trt_dat_ex = treatment_data_extra,
    target_probs_control = target_probs_control,
    target_probs_treatment = target_probs_treatment,
    file_name = filename,
    plot_outcome = plot_outcome
  )

  spline_tp_ctrl <- posterior_spline_matrix(
    smol_d3$fit_ctrl,
    smol_d3$ctrl_spline_basis,
    newx = age_grid
  )
  spline_tp_ctrl_ex <- posterior_spline_matrix(
    smol_d3$fit_ctrl_ex,
    smol_d3$ex_ctrl_spline_basis,
    newx = age_grid
  )
  spline_tp_trt <- posterior_spline_matrix(
    smol_d3$fit_trt,
    smol_d3$trt_spline_basis,
    newx = age_grid
  )
  spline_tp_trt_ex <- posterior_spline_matrix(
    smol_d3$fit_trt_ex,
    smol_d3$ex_trt_spline_basis,
    newx = age_grid
  )

  res_spline_multi_g5 <- run_multi_cohort_sequential(
    start_ages = start_ages[[1]],
    spline_tp_ctrl,
    spline_tp_trt,
    n_cycles = n_cycles,
    utilities = utilities,
    cost_ctrl = cost_ctrl,
    cost_trt = cost_trt
  )

  res_spline_multi_ex_g5 <- run_multi_cohort_sequential(
    start_ages = start_ages[[1]],
    spline_tp_ctrl_ex,
    spline_tp_trt_ex,
    n_cycles = n_cycles,
    utilities = utilities,
    cost_ctrl = cost_ctrl,
    cost_trt = cost_trt
  )

  return(list(
    mono_spline_limited = res_spline_multi_g5,
    mono_spline_extended = res_spline_multi_ex_g5
  ))
}
