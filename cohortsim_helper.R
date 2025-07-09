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

run_multi_cohort_parallel <- function(
  start_ages,
  final_tp_ctrl,
  final_tp_trt,
  n_cycles,
  utilities,
  cost_ctrl,
  cost_trt
) {
  final_tp_ctrl <<- final_tp_ctrl
  final_tp_trt <<- final_tp_trt

  print("done super assignment")

  cl <- makeCluster(5)

  # Exporting all necessary variables and functions to the cluster
  clusterExport(
    cl,
    c(
      "cohort_state_transition",
      "n_cycles",
      "utilities",
      "compute_delta",
      "cost_ctrl",
      "cost_trt",
      "extended_age_grid",
      "n_target"
    )
  )

  multi_cohort_stm <- parLapply(cl, start_ages, function(start_age) {
    cohort_list <- list()
    age_idx_range <- which(
      extended_age_grid >= start_age & extended_age_grid < (start_age + 40)
    )
    print("done calculate idx")
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

  stopCluster(cl)

  return(multi_cohort_stm)
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

read_file_avg <- function(file_names, scenario_name, n_target, weights) {
  # Find the index of the file name that contains the scenario name
  file_index <- grep(scenario_name, file_names)

  # Read the file
  dt <- readRDS(file_names[[file_index]])

  dt_collapsed <- multi_cohort_average(dt, weights)

  # Calculate the averages and store them in a list
  avg_list <- lapply(dt_collapsed, function(method_result) {
    lapply(method_result, function(y) y / n_target)
  })

  dt_collapsed <- map(avg_list, bind_rows)

  # Return the averages
  return(dt_collapsed)
}

import_plots <- function(plots_names, scenario_name) {
  # Detect the plots based on the scenario_name and model_name
  get_y_limits <- function(plot) {
    data_list <- ggplot_build(plot)$data
    y_range_list <- lapply(data_list, function(layer_data) {
      range(layer_data$y, na.rm = TRUE)
    })
    do.call("range", y_range_list)
  }

  linear_plots <- plots_names[str_detect(
    plots_names,
    paste0(scenario_name, ".*", "linear")
  )] %>%
    readRDS(.)
  unrestricted_spline_plots <- plots_names[str_detect(
    plots_names,
    paste0(scenario_name, ".*", "unrestricted_spline")
  )] %>%
    readRDS(.)
  mono_spline_plots <- plots_names[str_detect(
    plots_names,
    paste0(scenario_name, ".*(?<!unrestricted_)spline")
  )] %>%
    readRDS(.)
  adjusted_plots <- plots_names[str_detect(
    plots_names,
    paste0(scenario_name, ".*", "adjusted")
  )] %>%
    readRDS(.)
  # Assuming selected_plots is a data frame, extract the outcomes and effects

  # Combine y-limits to find the global range
  all_y_limits <- c(
    sapply(linear_plots, get_y_limits),
    sapply(unrestricted_spline_plots, get_y_limits),
    sapply(mono_spline_plots, get_y_limits),
    sapply(adjusted_plots, get_y_limits)
  )

  y_max <- max(all_y_limits, na.rm = TRUE) + 0.2
  y_min <- min(all_y_limits, na.rm = TRUE)
  y_min <- ifelse(y_min - 0.2 < 0, 0, y_min - 0.2)

  global_y_range <- c(y_min, y_max)

  # Function to set the global y-axis limits
  set_global_ylimits <- function(plot) {
    plot + scale_y_continuous(limits = global_y_range)
  }
  # Apply the new y-axis limits to all plots
  linear_plots <- lapply(linear_plots, set_global_ylimits)
  unrestricted_spline_plots <- lapply(
    unrestricted_spline_plots,
    set_global_ylimits
  )
  mono_spline_plots <- lapply(mono_spline_plots, set_global_ylimits)
  adjusted_plots <- lapply(adjusted_plots, set_global_ylimits)

  outcome_plots <- ggarrange(
    plotlist = c(
      adjusted_plots,
      linear_plots,
      unrestricted_spline_plots,
      mono_spline_plots
    ),
    nrow = 4,
    ncol = 2,
    common.legend = TRUE,
    legend = "bottom"
  )
  outcome_title <- paste0(
    "outcome surface under limited/extended data",
    "\n",
    "Method from top to bottom: adjusted regressions,linear interaction, unrestricted spline, monotonic spline"
  )

  outcome_plots <- annotate_figure(
    outcome_plots,
    top = text_grob(outcome_title, hjust = 0, x = 0, face = "bold")
  )

  return(list(outcome_plots = outcome_plots))
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

plot_INMB_multi <- function(method_list, res_true, k) {
  # Calculate INMB for each data frame in the list
  res_INMB <- lapply(names(method_list), function(method_name) {
    df <- method_list[[method_name]]
    INMB <- k * df$delta_QALYs - df$delta_costs
    INMB_true <- k * res_true$delta_QALYs - res_true$delta_costs
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
    ggplot(aes(y = method, x = INMB, fill = method)) +
    stat_slab() +
    geom_vline(xintercept = df_combined$INMB_true, linetype = "dashed") +
    stat_pointinterval(alpha = 0.6) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = paste0("Incremental Net Monetary Benefit, k = ", k),
      ylab = "Method",
    ) +
    facet_grid(
      scenarios ~ .,
      scales = "free_y",
      space = "free_y",
      as.table = FALSE
    ) +
    theme(
      plot.title = element_text(face = "bold"), # Bold title
      axis.text.y = element_text(face = "bold") # Bold y-axis labels
    )

  return(plot)
}
