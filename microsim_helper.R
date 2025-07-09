# this file contains functions for performing microsimulation
# based on a fitted model
# @microsim_post_cmdstan simulates the posterior of cost and
# QALYs based on a fitted cmdstan model
# @microsimulation_posterior simulates the posterior based on
# a stanreg obejct
# @cmdstan_epred generates posterior predictive distribution
# for the conditional mean given a fitted stan model

# cmdstan_epred <- function(fit, newdata, fitted_basis, trt){
#
#   draws_df <- fit$draws(format = "df")
#
#   index_s0 <- grep("^s_0\\[", names(draws_df))
#   index_s1 <- grep("^s_1\\[", names(draws_df))
#
#   s0_df <- draws_df[, index_s0] %>% as.matrix()
#   s1_df <- draws_df[, index_s1] %>% as.matrix()
#   b0 <- draws_df$Intercept
#   b_trt <- draws_df$b
#
#
#   new_basis <- predict(fitted_basis, newdata)
#
#
#
#   y_pred <- arm::invlogit(b0 + trt*b_trt + s0_df %*% t(new_basis) + trt * s1_df %*% t(new_basis))
#   attributes(y_pred)$dimnames[[2]] <- NULL
#   return(y_pred)
# }

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

microsim_post_cmdstan <- function(
  fit,
  target_age,
  fitted_basis,
  trt,
  n_simulations,
  n_patients,
  n_cycles,
  n_states,
  treatment,
  treatment_cost,
  utility_healthy,
  utility_diseased,
  cost_healthy,
  cost_diseased,
  sep = FALSE
) {
  QALYs <- numeric(n_simulations)
  costs <- numeric(n_simulations)

  current_state_matrix <- matrix(1, nrow = n_simulations, ncol = n_patients)

  for (cycle in 1:n_cycles) {
    start_state_counts <- table(as.vector(current_state_matrix))
    cat("State counts at start of cycle", cycle, ":")
    print(start_state_counts)

    if (sep) {
      posterior_probs <- cmdstan_epred_sep(fit, fitted_basis, target_age)
    } else {
      posterior_probs <- cmdstan_epred(
        fit,
        target_age,
        fitted_basis,
        target_age
      )
    }

    # Update age
    target_age <- target_age + 1

    # Generate all random numbers
    rand_num_matrix <- matrix(
      runif(n_simulations * n_patients, 0, 1),
      nrow = n_simulations,
      ncol = n_patients
    )

    # Transition 1 to 2
    transition_1_to_2 <- current_state_matrix == 1 &
      rand_num_matrix < posterior_probs
    current_state_matrix[transition_1_to_2] <- 2

    # Transition 2 to 3, but only consider those that did not transition from 1 to 2 in this cycle
    eligible_for_2_to_3 <- !transition_1_to_2 & current_state_matrix == 2
    transition_2_to_3 <- eligible_for_2_to_3 & rand_num_matrix < 0.1
    current_state_matrix[transition_2_to_3] <- 3

    # Calculate QALYs and costs
    QALYs <- QALYs +
      utility_healthy * rowSums(current_state_matrix == 1) +
      utility_diseased * rowSums(current_state_matrix == 2)
    costs <- costs +
      (cost_healthy + treatment_cost * treatment) *
        rowSums(current_state_matrix == 1) +
      cost_diseased * rowSums(current_state_matrix == 2)

    cat(
      "Transition counts after cycle",
      cycle,
      ": Transition 1 to 2:",
      sum(transition_1_to_2),
      ", Transition 2 to 3:",
      sum(transition_2_to_3),
      "\n"
    )

    end_state_counts <- table(as.vector(current_state_matrix))
    cat("State counts after cycle", cycle, ":")
    print(end_state_counts)
  }

  return(list(QALYs = QALYs, costs = costs))
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

# posterior_long <- function(pred, treatment_label, age_vector) {
#   pred <- as.data.frame(t(pred)) %>%
#     mutate(age = age_vector) %>%
#     pivot_longer(cols = -age, names_to = "sample", values_to = "probability") %>%
#     mutate(treatment = treatment_label)
# }

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
  # ... (include the outcome_data, posterior_long, and plot_outcome_surface functions here) ...

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

plot_proportions <- function(control_input, treatment_input, n_cycle, n_state) {
  if (identical(control_input, treatment_input)) {
    control_matrix <- control_input$trace_ctrl
    treatment_matrix <- control_input$trace_trt
  } else {
    control_matrix <- control_input$aggregate_state_distribution
    treatment_matrix <- treatment_input$aggregate_state_distribution
  }

  # Add column names and row names to matrices
  colnames(control_matrix) <- c("Healthy", "Diseased", "Death")
  colnames(treatment_matrix) <- c("Healthy", "Diseased", "Death")

  # Find the cycle where Diseased overtake Healthy for both groups
  # control_overtake_cycle <- which.max(control_matrix[, "Diseased"] > control_matrix[, "Healthy"])
  # treatment_overtake_cycle <- which.max(treatment_matrix[, "Diseased"] > treatment_matrix[, "Healthy"])
  #
  # # Create a data frame for vertical lines
  # vline_data <- data.frame(
  #   xintercept = c(control_overtake_cycle, treatment_overtake_cycle),
  #   group = c("Control", "Treatment")
  # )
  #

  # Convert matrices to data frames and add cycle and group information
  control_df <- as.data.frame(control_matrix) %>%
    mutate(cycle = 0:n_cycle) %>%
    pivot_longer(
      cols = -cycle,
      names_to = "state",
      values_to = "proportion"
    ) %>%
    mutate(group = "Control")

  treatment_df <- as.data.frame(treatment_matrix) %>%
    mutate(cycle = 0:n_cycle) %>%
    pivot_longer(
      cols = -cycle,
      names_to = "state",
      values_to = "proportion"
    ) %>%
    mutate(group = "Treatment")

  # Combine both data frames
  combined_df <- bind_rows(control_df, treatment_df)

  # Create the ggplot
  p <- ggplot(combined_df, aes(x = cycle, y = proportion, color = state)) +
    geom_line() +
    # geom_vline(data = data.frame(xintercept = control_overtake_cycle, group = "Control"), aes(xintercept = xintercept), color = "black", linetype = "dashed") +
    # geom_text(data = data.frame(xintercept = control_overtake_cycle, group = "Control"), aes(x = xintercept, y = Inf, label = round(xintercept, 1)), vjust = 2, hjust = 0.5, inherit.aes = FALSE) +
    # geom_vline(data = data.frame(xintercept = treatment_overtake_cycle, group = "Treatment"), aes(xintercept = xintercept), color = "black", linetype = "dashed") +
    # geom_text(data = data.frame(xintercept = treatment_overtake_cycle, group = "Treatment"), aes(x = xintercept, y = Inf, label = round(xintercept, 1)), vjust = 2, hjust = 0.5, inherit.aes = FALSE) +
    facet_grid(. ~ group) +
    labs(title = "State Transition Over Time", x = "Cycle", y = "") +
    theme_minimal()

  return(p)
}

clean_method_names <- function(method_name) {
  gsub("cohort_", "", method_name) %>%
    gsub("unadjusted", "Constant", .) %>%
    gsub("adjusted", "Covariate adjusted", .) %>%
    gsub("linear", "Linear interaction", .) %>%
    gsub("bspline", "B-Spline RW", .) %>%
    gsub("spline", "Monotonic spline", .) %>%
    gsub("_ex", "", .)
}


# Function to generate combined ggplot
plot_delta_trace <- function(method_list, true_trace, y_axis_label = "Delta") {
  process_method <- function(df, method_name) {
    df %>%
      as.data.frame() %>%
      setNames(seq_len(ncol(df))) %>%
      pivot_longer(everything(), names_to = "cycle", values_to = "delta") %>%
      mutate(
        cycle = as.numeric(cycle),
        method = clean_method_names(method_name),
        scenario = ifelse(grepl("_ex", method_name), "Extended", "Limited")
      )
  }

  df_list <- lapply(names(method_list), function(m) {
    process_method(do.call(rbind, method_list[[m]]), m)
  }) %>%
    bind_rows()

  # Create a data frame for the true trace
  true_df <- data.frame(cycle = seq_along(true_trace), delta = true_trace)

  # Generate ggplot
  p <- ggplot(
    df_list,
    aes(x = cycle, y = delta, fill = method, colour = method)
  ) +
    stat_lineribbon(aes(fill_ramp = after_stat(level)), alpha = 0.4) %>%
      partition(vars(method)) %>%
      blend("multiply") +
    geom_line(
      data = true_df,
      aes(x = cycle, y = delta),
      color = "red",
      inherit.aes = FALSE,
      show.legend = F
    ) %>%
      copy_under(color = "white", linewidth = 2.5) + # Add true trace
    facet_grid(
      scenario ~ .,
      scales = "free_y",
      space = "free_y",
      as.table = FALSE
    ) +
    scale_fill_brewer(palette = "Set2") +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "Cycle", y = y_axis_label)

  return(p)
}

plot_INMB <- function(method_list, res_true, k) {
  # Calculate INMB for each data frame in the list
  res_INMB <- lapply(names(method_list), function(method_name) {
    df <- method_list[[method_name]]
    INMB <- k * df$delta_QALYs - df$delta_costs
    INMB_true <- k * res_true$qalys - res_true$costs
    data.frame(
      method = clean_method_names(method_name),
      scenarios = ifelse(grepl("_ex", method_name), "Extended", "Limited"),
      INMB = INMB,
      INMB_true = INMB_true
    )
  })

  # Bind the data frames together
  df_combined <- bind_rows(res_INMB)

  # Create the plot
  plot <- df_combined %>%
    ggplot(aes(y = method, x = INMB)) +
    stat_slab(
      aes(
        fill = after_stat(
          x >= min(df_combined$INMB_true) & x <= max(df_combined$INMB_true)
        )
      ),
      show.legend = F
    ) +
    stat_pointinterval() +
    scale_fill_brewer() +
    labs(
      title = paste0("Incremental Net Monetary Benefit,\n k = ", k)
    ) +
    facet_grid(
      scenarios ~ .,
      scales = "free_y",
      space = "free_y",
      as.table = FALSE
    )

  return(plot)
}
