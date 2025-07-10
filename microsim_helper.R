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


# Plot the outcome surfaces using ggplot2
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
