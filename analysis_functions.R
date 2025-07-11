compute_INMB <- function(sublist, wtp) {
  INMB = (wtp * sublist$delta_QALYs - sublist$delta_costs) / n_target

  return(INMB)
}

log_emp_density <- function(post_values, INMB_true) {
  dens <- density(post_values)
  dens_fun <- approxfun(dens$x, dens$y, rule = 2)
  log_dens_true <- log(dens_fun(INMB_true))
  return(log_dens_true)
}

combine_df_lpd <- function(list_all_mean, list_all_se) {
  df_all_limited_se <- list_all_se %>%
    bind_rows() %>%
    mutate_if(is.numeric, round, digits = 2)

  df_all_limited <- list_all_mean %>%
    bind_rows() %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    mutate(across(everything(), as.character))

  df_combined <- mapply(
    function(mean, se) {
      paste0(mean, " (", se, ")")
    },
    df_all_limited,
    df_all_limited_se,
    SIMPLIFY = F
  ) %>%
    as_tibble() %>%
    add_column(scenario = combinations$scenario_name, .before = 1)

  return(df_combined)
}

combine_df_pp <- function(list_all_mean, inmb_true_list) {
  elpd_all_limited_df <- list_all_mean %>%
    bind_rows() %>%
    mutate_if(is.numeric, round, digits = 2) %>%
    mutate(across(everything(), as.character))

  df_combined <- elpd_all_limited_df %>%
    as_tibble() %>%
    add_column(scenario = combinations$scenario_name, .before = 1) %>%
    add_column(INMB_true = round(unlist(inmb_true_list), 0), .before = 2)

  return(df_combined)
}

extract_info <- function(nested_list, scenario_name) {
  map2_dfr(
    nested_list,
    scenario_name,
    ~ imap_dfr(.x, function(method_list, method_name) {
      tibble(
        scenario = .y,
        Methods = method_name,
        lower_CI = method_list[["LL"]],
        upper_CI = method_list[["UL"]]
      )
    })
  )
}

plot_avg_uncertainty <- function(df, plot_title, x_label, true_df) {
  plot_pointrange <- function(df, plot_title, x_label, true_df) {
    # Assign a unique numeric ID to each scenario in df
    df$scenario_id <- as.numeric(factor(
      df$scenario,
      levels = unique(df$scenario)
    )) *
      0.5

    # Create a new column for y-axis position with a small offset for each method in df
    method_offsets <- seq(
      from = -0.1,
      to = 0.1,
      length.out = length(unique(df$Methods))
    )
    names(method_offsets) <- unique(df$Methods)
    df$y_position <- df$scenario_id + method_offsets[df$Methods]

    p <- ggplot(
      df,
      aes(x = value, y = y_position, shape = Methods, color = Methods)
    ) +
      geom_pointrange(
        aes(xmin = lower_CI, xmax = upper_CI, group = Methods),
        size = 0.5
      ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")
      ) +
      scale_y_continuous(
        labels = unique(df$scenario),
        breaks = unique(df$scenario_id),
        expand = c(0, 0.1)
      ) +
      labs(x = x_label, y = "") +
      ggtitle(plot_title)

    # Map the scenarios in true_df to the corresponding numeric IDs used in df
    true_df$scenario <- str_wrap(true_df$scenario, width = 25)

    # Filter true_df to include only the scenarios present in df
    true_df_filtered <- true_df[true_df$scenario %in% df$scenario, ]
    true_df_filtered$scenario_id <- match(
      true_df_filtered$scenario,
      unique(df$scenario)
    ) *
      0.5
    true_df_filtered$y_position <- true_df_filtered$scenario_id

    p <- p +
      geom_point(
        data = true_df_filtered,
        aes(x = INMB, y = y_position),
        color = "red",
        size = 3,
        shape = 17
      )

    return(p)
  }

  df1 <- df %>% filter(scenario %in% combinations$scenario_name[1:6])
  df1$scenario <- str_wrap(df1$scenario, width = 25)
  df2 <- df %>% filter(scenario %in% combinations$scenario_name[7:12])
  df2$scenario <- str_wrap(df2$scenario, width = 25)

  p1 <- plot_pointrange(df1, "", x_label, true_df)
  p2 <- plot_pointrange(df2, "", x_label, true_df)

  comb_plot <- ggarrange(
    p1,
    p2,
    ncol = 2,
    common.legend = TRUE,
    legend = "bottom"
  )

  comb_plot <- annotate_figure(
    comb_plot,
    top = text_grob(
      "Posterior mean INMB across all Monte Carlo replications",
      face = "bold",
      size = 14
    )
  )

  return(comb_plot)
}

generate_se_plots <- function(data, metric_col, plot_title, x_label, y_label) {
  # Validate input data structure
  stopifnot(ncol(data) >= 6, any(colnames(data) == "scenario"))

  # Helper function to create plot
  create_plot <- function(
    df,
    metric_col,
    SE,
    add_true,
    plot_title,
    x_label,
    y_label
  ) {
    df_long <- melt(df, id.vars = "scenario")
    colnames(df_long)[colnames(df_long) == "variable"] <- "Methods"
    df_long$scenario <- str_wrap(df_long$scenario, width = 25)
    df_long[[metric_col]] <- as.numeric(sub("\\(.*\\)", "", df_long$value))
    df_long$SE <- as.numeric(sub(".*\\(", "", sub("\\).*", "", df_long$value)))

    # Assign a unique numeric ID to each scenario in df
    df_long$scenario_id <- as.numeric(factor(
      df_long$scenario,
      levels = unique(df_long$scenario)
    ))

    # Create a new column for y-axis position with a small offset for each method in df
    method_offsets <- seq(
      from = -0.1,
      to = 0.1,
      length.out = length(unique(df_long$Methods))
    )
    names(method_offsets) <- unique(df_long$Methods)
    df_long$y_position <- df_long$scenario_id + method_offsets[df_long$Methods]

    p <- ggplot(
      df_long,
      aes(x = get(metric_col), y = y_position, shape = Methods, color = Methods)
    ) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")
      ) +
      scale_y_continuous(
        labels = unique(df_long$scenario),
        breaks = unique(df_long$scenario_id)
      ) +
      labs(x = x_label, y = y_label) +
      ggtitle(plot_title)

    p <- p +
      geom_pointrange(
        aes(
          xmin = get(metric_col) - 1.96 * SE,
          xmax = get(metric_col) + 1.96 * SE,
          group = Methods
        ),
        size = 0.5
      )

    return(p)
  }

  # Split into two datasets
  df1 <- data[1:6, ]
  df2 <- data[7:12, ]

  # Create plots
  p1 <- create_plot(
    df1,
    metric_col,
    SE,
    add_true,
    paste0(plot_title),
    x_label,
    y_label
  )
  p2 <- create_plot(
    df2,
    metric_col,
    SE,
    add_true,
    paste0(plot_title),
    x_label,
    y_label
  )

  # Combine plots
  comb_plot <- ggarrange(
    p1,
    p2,
    ncol = 2,
    common.legend = TRUE,
    legend = "bottom"
  )

  return(comb_plot)
}

generate_dot_plots <- function(
  data,
  metric_col,
  plot_title,
  x_label,
  y_label,
  jitter_width = 0.01
) {
  # Validate input data structure
  stopifnot(ncol(data) >= 6, any(colnames(data) == "scenario"))

  # Helper function to create plot
  create_plot <- function(df, metric_col, plot_title, x_label, y_label) {
    df_long <- melt(df, id.vars = "scenario")
    colnames(df_long)[colnames(df_long) == "variable"] <- "Methods"
    df_long$scenario <- str_wrap(df_long$scenario, width = 25)
    df_long[[metric_col]] <- as.numeric(sub("\\(.*\\)", "", df_long$value))

    p <- ggplot(
      df_long,
      aes(x = get(metric_col), y = scenario, shape = Methods, color = Methods)
    ) +
      geom_point(size = 4, position = position_dodge(width = jitter_width)) +
      theme_minimal() +
      theme(
        plot.title = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold")
      ) +
      labs(x = x_label, y = y_label) +
      ggtitle(plot_title)

    return(p)
  }

  # Split into two datasets
  df1 <- data[1:6, ]
  df2 <- data[7:12, ]

  # Create plots
  p1 <- create_plot(df1, metric_col, paste0(plot_title), x_label, y_label)
  p2 <- create_plot(df2, metric_col, paste0(plot_title), x_label, y_label)

  # Combine plots
  comb_plot <- ggarrange(
    p1,
    p2,
    ncol = 2,
    common.legend = TRUE,
    legend = "bottom"
  )

  return(comb_plot)
}

pp_approval <- function(filepath, scenario, wtp) {
  dat <- readRDS(paste0(filepath, scenario, ".RDs"))
  INMB_list <- map(dat, ~ compute_INMB(., wtp))

  pp_list <- map(INMB_list, ~ mean(.x > 0))
  avg_pp <- mean(unlist(pp_list))
  CI <- round(quantile(unlist(pp_list), c(0.025, 0.975)), 2)
  CI_pp <- paste0(CI[1], ",", CI[2])

  return(list(avg_pp = avg_pp, CI_pp = CI_pp))
}

elpd_emp_INMB <- function(filepath, scenario, wtp, true) {
  dat <- readRDS(paste0(filepath, scenario, ".RDs"))
  INMB_list <- map(dat, ~ compute_INMB(., wtp))
  lpd_INMB <- map(INMB_list, ~ log_emp_density(., true))
  elpd_INMB <- mean(unlist(lpd_INMB))
  se_lpd <- sd(unlist(lpd_INMB))

  return(list(elpd_INMB = elpd_INMB, se_lpd = se_lpd))
}

INMB_avg <- function(filepath, scneario, wtp) {
  dat <- readRDS(paste0(filepath, scneario, ".RDs"))
  INMB_list <- map(dat, ~ compute_INMB(., wtp))

  avg_list <- lapply(INMB_list, mean)
  modExp_INMB <- mean(unlist(avg_list))
  ll_INMB <- quantile(unlist(avg_list), 0.025)
  ul_INMB <- quantile(unlist(avg_list), .975)

  return(list(average = modExp_INMB, ll_INMB = ll_INMB, ul_INMB = ul_INMB))
}

model_scenario_pp <- function(scenario_name, model_type, suffix, wtp) {
  path <- file.path(here("results", model_type), scenario_name)
  return(pp_approval(path, suffix, wtp)) # pass wtp
}

model_scenario_elpd <- function(scenario_name, model_type, suffix, wtp, true) {
  path <- file.path(here("results", model_type), scenario_name)
  return(elpd_emp_INMB(path, suffix, wtp, true)) # pass wtp and true to elpd_emp_INMB
}

model_scenario_avg <- function(scenario_name, model_type, suffix, wtp) {
  path <- file.path(here("results", model_type), scenario_name)
  return(INMB_avg(path, suffix, wtp)) # pass wtp
}