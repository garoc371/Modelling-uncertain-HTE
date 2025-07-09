library(tidyverse)
library(stringr)
library(patchwork)
library(ggpubr)
# library(Rcpp)
source("microsim_spline/cohortsim_helper.R")
source("microsim_spline/microsim_helper.R")
library(ggdist)
library(reshape2)



weight_names <- list.files(path = "./microsim_spline/MC_sim/results2/true", 
                           pattern = "*weights_true.RDs", 
                           full.names = TRUE)


cohort_true_names <- list.files(path = "./microsim_spline/MC_sim/results2/true", 
                                pattern = "*cohort_true.RDs", 
                                full.names = TRUE)

n_target <- 20000
wtp_1 <- 15000
wtp_2 <- 17500


start_ages <- c(40,50,60,70,80)

# function to calculate INMB for each MC replication

compute_INMB <-function(sublist, wtp){
  INMB = (wtp * sublist$delta_QALYs - sublist$delta_costs)/n_target
  
  return(INMB)
}

log_emp_density <- function(post_values, INMB_true){
  dens <- density(post_values)
  dens_fun <- approxfun(dens$x, dens$y, rule = 2)
  log_dens_true <- log(dens_fun(INMB_true))
  return(log_dens_true)
}


# wrapper for reading and calculating elpd for INMB
# @ scenario is either limited or extended

elpd_emp_INMB <- function(filepath, scenario, wtp, true){
  dat <- readRDS(paste0(filepath, scenario, ".RDs"))
  INMB_list <- map(dat, ~compute_INMB(., wtp))
  lpd_INMB <- map(INMB_list, ~log_emp_density(., true))
  elpd_INMB <- mean(unlist(lpd_INMB))
  se_lpd <- sd(unlist(lpd_INMB))
  
  return(list(elpd_INMB = elpd_INMB, se_lpd = se_lpd))
}

model_scenario_elpd <- function(scenario_name, model_type, suffix, wtp, true) {
  path <- paste0("microsim_spline/MC_sim/results/", model_type, "/", scenario_name)
  return(elpd_emp_INMB(path, suffix, wtp, true))  # pass wtp and true to elpd_emp_INMB
}

# calculate the average posterior probability of approval
# across all MC simulations and the MC standard errors

pp_approval <- function(filepath, scenario, wtp){
  dat <- readRDS(paste0(filepath, scenario, ".RDs"))
  INMB_list <- map(dat, ~compute_INMB(., wtp))
  
  pp_list <- map(INMB_list, ~mean(.x>0))
  avg_pp <- mean(unlist(pp_list))
  CI <- round(quantile(unlist(pp_list), c(0.025, 0.975)),2)
  CI_pp <- paste0(CI[1], ",", CI[2])
  
  return(list(avg_pp = avg_pp, CI_pp = CI_pp))
}


INMB_avg <- function(filepath, scneario, wtp){
  dat <- readRDS(paste0(filepath, scneario, ".RDs"))
  INMB_list <- map(dat, ~compute_INMB(., wtp))
  
  avg_list <- lapply(INMB_list, mean)
  modExp_INMB <- mean(unlist(avg_list))
  ll_INMB <- quantile(unlist(avg_list), 0.025)
  ul_INMB <- quantile(unlist(avg_list), .975)
  
  
  return(list(average = modExp_INMB, ll_INMB = ll_INMB, ul_INMB = ul_INMB))
}

model_scenario_avg <- function(scenario_name, model_type, suffix, wtp) {
  path <- paste0("microsim_spline/MC_sim/results/", model_type, "/", scenario_name)
  return(INMB_avg(path, suffix, wtp))  # pass wtp
}



model_scenario_pp <- function(scenario_name, model_type, suffix, wtp) {
  path <- paste0("microsim_spline/MC_sim/results/", model_type, "/", scenario_name)
  return(pp_approval(path, suffix, wtp))  # pass wtp
}



combinations <- expand.grid(
  control_type = c("non_linear_mon_increase", "non_linear_mon_decrease", "non_linear_non_mon"),
  treatment_type = c("constant", 
                     "non_linear_mon_increase", "non_linear_mon_decrease", "non_linear_non_mon")
)

# for wtp = 20k
pp_all_limited <- vector(mode = "list", length = length(combinations))
pp_all_extended <- vector(mode = "list", length = length(combinations))
elpd_all_limited <- vector(mode = "list", length = length(combinations))
elpd_all_extended <- vector(mode = "list", length = length(combinations))
elpd_all_limited_se <- vector(mode = "list", length = length(combinations))
elpd_all_extended_se <- vector(mode = "list", length = length(combinations))

avg_INMB_all_limited <- vector(mode = "list", length = length(combinations))
avg_INMB_all_extended <- vector(mode = "list", length = length(combinations))
CI_INMB_all_limited <- vector(mode = "list", length = length(combinations))
CI_INMB_all_extended <- vector(mode = "list", length = length(combinations))



# for wtp = 30k
pp_all_limited_30k <- vector(mode = "list", length = length(combinations))
pp_all_extended_30k <- vector(mode = "list", length = length(combinations))
elpd_all_limited_30k <- vector(mode = "list", length = length(combinations))
elpd_all_extended_30k <- vector(mode = "list", length = length(combinations))
elpd_all_limited_30k_se <- vector(mode = "list", length = length(combinations))
elpd_all_extended_30k_se <- vector(mode = "list", length = length(combinations))


model_types <- c("unadjusted", "adjusted", "linear_interaction", "unrestricted_spline", "monotonic_spline")
INMB_true_20k_all <- INMB_true_30k_all <- list()

for (i in 1:nrow(combinations)) {
  control_type <- combinations$control_type[i]
  treatment_type <- combinations$treatment_type[i]
  
  pp_limited <- pp_extended <- pp_limited_30k <- pp_extended_30k  <-
    setNames(vector("list", length(model_types)), model_types)
  
  avg_INMB_limited <- avg_INMB_extended <- 
    CI_INMB_limited <- CI_INMB_extended <-
    setNames(vector("list", length(model_types)), model_types)
  
  
  elpd_limited <- elpd_extended <- elpd_limited_se <- elpd_extended_se <-
    elpd_limited_30k <- elpd_extended_30k <- elpd_limited_30k_se <- elpd_extended_30k_se <-
    setNames(vector("list", length(model_types)), model_types)
  
                                         
  
  # Shorten names
  control_short <- gsub("non_linear", "non", gsub("linear", "lin", as.character(control_type)))
  treatment_short <- gsub("non_linear", "non", gsub("linear", "lin", as.character(treatment_type)))
  
  scenario_name <- paste0(control_short, "_ctrl_", treatment_short, "_trt")
  
  true_weights_current <- readRDS(weight_names[[grep(scenario_name, weight_names)]])
  weights_current_g5 <- aggregate_weights(true_weights_current, start_ages)
  
  res_true_current <- readRDS(cohort_true_names[[grep(scenario_name, cohort_true_names)]])
  INMB_true_20k <- (wtp_1 * res_true_current$delta_QALYs - res_true_current$delta_costs)/n_target
  INMB_true_30k <- (wtp_2 * res_true_current$delta_QALYs - res_true_current$delta_costs)/n_target
  
  INMB_true_20k_all[[i]] <- INMB_true_20k; INMB_true_30k_all[[i]] <- INMB_true_30k
  

  
  for (model_type in model_types){
    elpd_temp_limited <- model_scenario_elpd(scenario_name, model_type, "_limited", wtp_1, INMB_true_20k)
    elpd_temp_limited_30k <- model_scenario_elpd(scenario_name, model_type, "_limited", wtp_2, INMB_true_30k)
    elpd_temp_extended <- model_scenario_elpd(scenario_name, model_type, "_extended", wtp_1, INMB_true_20k)
    elpd_temp_extended_30k <- model_scenario_elpd(scenario_name, model_type, "_extended", wtp_2, INMB_true_30k)
    
    elpd_limited[[model_type]] <- elpd_temp_limited$elpd_INMB
    elpd_limited_se[[model_type]] <- elpd_temp_limited$se_lpd
    elpd_limited_30k[[model_type]] <- elpd_temp_limited_30k$elpd_INMB
    elpd_limited_30k_se[[model_type]] <- elpd_temp_limited_30k$se_lpd
    
    elpd_extended[[model_type]] <- elpd_temp_extended$elpd_INMB
    elpd_extended_se[[model_type]] <- elpd_temp_extended$se_lpd
    elpd_extended_30k[[model_type]] <- elpd_temp_extended_30k$elpd_INMB
    elpd_extended_30k_se[[model_type]] <- elpd_temp_extended_30k$se_lpd
    
  }
  
  elpd_all_limited[[i]] <- elpd_limited
  elpd_all_extended[[i]] <- elpd_extended
  elpd_all_limited_30k[[i]] <- elpd_limited_30k
  elpd_all_extended_30k[[i]] <- elpd_extended_30k
  
  elpd_all_limited_se[[i]] <- elpd_limited_se
  elpd_all_extended_se[[i]] <- elpd_extended_se
  elpd_all_limited_30k_se[[i]] <- elpd_limited_30k_se
  elpd_all_extended_30k_se[[i]] <- elpd_extended_30k_se

  for (model_type in model_types){
    pp_temp_limited <- model_scenario_pp(scenario_name, model_type, "_limited", wtp_1)
    pp_temp_limited_30k <- model_scenario_pp(scenario_name, model_type, "_limited", wtp_2)
    pp_temp_extended <- model_scenario_pp(scenario_name, model_type, "_extended", wtp_1)
    pp_temp_extended_30k <- model_scenario_pp(scenario_name, model_type, "_extended", wtp_2)

    pp_limited[[model_type]] <- pp_temp_limited$avg_pp
    pp_limited_30k[[model_type]] <- pp_temp_limited_30k$avg_pp

    pp_extended[[model_type]] <- pp_temp_extended$avg_pp
    pp_extended_30k[[model_type]] <- pp_temp_extended_30k$avg_pp

  }

  pp_all_limited[[i]] <- pp_limited
  pp_all_extended[[i]] <- pp_extended
  pp_all_limited_30k[[i]] <- pp_limited_30k
  pp_all_extended_30k[[i]] <- pp_extended_30k
  
  
  for (model_type in model_types){
    temp_INMB_limited <- model_scenario_avg(scenario_name, model_type, "_limited", wtp_1)
    temp_INMB_extended <- model_scenario_avg(scenario_name, model_type, "_extended", wtp_1)
    
    avg_INMB_limited[[model_type]] <- temp_INMB_limited$average
    CI_INMB_limited[[model_type]] <- list(LL = temp_INMB_limited$ll_INMB, UL = temp_INMB_limited$ul_INMB)
    avg_INMB_extended[[model_type]] <- temp_INMB_extended$average
    CI_INMB_extended[[model_type]] <- list(LL = temp_INMB_extended$ll_INMB, UL = temp_INMB_extended$ul_INMB)
    
  }

  avg_INMB_all_limited[[i]] <- avg_INMB_limited
  avg_INMB_all_extended[[i]] <- avg_INMB_extended
  CI_INMB_all_limited[[i]] <- CI_INMB_limited
  CI_INMB_all_extended[[i]] <- CI_INMB_extended
}




combinations$scenario_name <- paste0(gsub("_", " ", combinations$control_type),
                                     " control ",
                                     gsub("_", " ", combinations$treatment_type),
                                     " treatment")

combine_df_lpd <- function(list_all_mean, list_all_se){
  
  
  df_all_limited_se <- list_all_se %>% 
    bind_rows() %>% 
    mutate_if(is.numeric, round, digits = 2) 
  
  df_all_limited <- list_all_mean %>%
    bind_rows() %>%
    mutate_if(is.numeric, round, digits = 2)  %>% 
    mutate(across(everything(), as.character))
  
  df_combined <- mapply(function(mean, se){
    paste0(mean, " (", se, ")")
  },df_all_limited, df_all_limited_se, SIMPLIFY = F) %>% 
    as_tibble() %>% 
    add_column(scenario = combinations$scenario_name, .before = 1)
  
  return(df_combined)
}

elpd_limited_15k <- combine_df_lpd(elpd_all_limited, elpd_all_limited_se)
elpd_limited_17k <- combine_df_lpd(elpd_all_limited_30k, elpd_all_limited_30k_se)

# avg_all_limited <- combine_df_lpd(list_all_mean = avg_INMB_all_limited, empse_INMB_all_limited) 

true_df <- tibble(INMB = unlist(INMB_true_20k_all),
                  scenario = combinations$scenario_name)


avg_all_limited_with_true <- left_join(avg_all_limited,true_df, by = "scenario")

# kable(limited_15k, "latex", booktabs = TRUE, escape = FALSE) %>%
#   kable_styling(latex_options = c("striped", "hold_position"), full_width = F) %>%
#   column_spec(2, width = "5cm")


combine_df_pp <- function(list_all_mean, inmb_true_list){


  elpd_all_limited_df <- list_all_mean %>%
    bind_rows() %>%
    mutate_if(is.numeric, round, digits = 2)  %>%
    mutate(across(everything(), as.character))

  df_combined <- elpd_all_limited_df %>%
    as_tibble() %>%
    add_column(scenario = combinations$scenario_name, .before = 1) %>%
    add_column(INMB_true = round(unlist(inmb_true_list), 0), .before = 2)

  return(df_combined)
}
#
pp_limited_15k <- combine_df_pp(pp_all_limited, INMB_true_20k_all)

# Function to extract information and return a tibble
extract_info <- function(nested_list, scenario_name) {
  map2_dfr(nested_list, scenario_name, ~imap_dfr(.x, function(method_list, method_name) {
    tibble(
      scenario = .y,
      Methods = method_name,
      lower_CI = method_list[["LL"]],
      upper_CI = method_list[["UL"]]
    )
  }))
}

CI_all_limited <- extract_info(CI_INMB_all_limited, combinations$scenario_name)
avg_all_limited <- avg_INMB_all_limited %>% bind_rows() %>% 
  add_column(scenario = combinations$scenario_name, .before = 1) %>% 
  melt(., id.vars = "scenario")


colnames(avg_all_limited)[colnames(avg_all_limited) == "variable"] <- "Methods"
avg_all_limited_full <- left_join(avg_all_limited, CI_all_limited, by = c("scenario", "Methods")) %>% 
  mutate_at(vars(last_col():(last_col()-2)), round, digits = 1)


plot_avg_uncertainty <- function(df, plot_title, x_label, true_df){
  
  
  plot_pointrange <- function(df, plot_title, x_label, true_df){
    
    # Assign a unique numeric ID to each scenario in df
    df$scenario_id <- as.numeric(factor(df$scenario, levels = unique(df$scenario)))*0.5
    
    # Create a new column for y-axis position with a small offset for each method in df
    method_offsets <- seq(from = -0.1, to = 0.1, length.out = length(unique(df$Methods)))
    names(method_offsets) <- unique(df$Methods)
    df$y_position <- df$scenario_id + method_offsets[df$Methods]
    
    
    p <- ggplot(df, aes(x = value, y = y_position, shape = Methods, color = Methods)) +
      geom_pointrange(aes(xmin = lower_CI, xmax = upper_CI, group = Methods),
                       size = 0.5) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            axis.text.y = element_text(face = "bold")) +
      scale_y_continuous(labels = unique(df$scenario), breaks = unique(df$scenario_id),
                         expand = c(0,0.1)) +
      labs(x = x_label, y = "") +
      ggtitle(plot_title)
    
    # Map the scenarios in true_df to the corresponding numeric IDs used in df
    true_df$scenario <- str_wrap(true_df$scenario, width = 25)
    
    
    # Filter true_df to include only the scenarios present in df
    true_df_filtered <- true_df[true_df$scenario %in% df$scenario, ]
    true_df_filtered$scenario_id <- match(true_df_filtered$scenario, unique(df$scenario))*0.5
    true_df_filtered$y_position <- true_df_filtered$scenario_id
    
    p <- p + geom_point(data = true_df_filtered, aes(x = INMB, y = y_position),
                        color = "red", size = 3, shape = 17)
    
    return(p)
  }
  
  df1 <- df %>% filter(scenario %in% combinations$scenario_name[1:6])
  df1$scenario <- str_wrap(df1$scenario, width = 25)
  df2 <- df %>% filter(scenario %in% combinations$scenario_name[7:12])
  df2$scenario <- str_wrap(df2$scenario, width = 25)
  
  p1 <- plot_pointrange(df1, "", x_label, true_df)
  p2 <- plot_pointrange(df2, "", x_label, true_df)
  
  comb_plot <- ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom")
  
  comb_plot <- annotate_figure(comb_plot,
                               top = text_grob("Posterior mean INMB across all Monte Carlo replications", 
                                               face = "bold", size = 14))
  
  return(comb_plot)
}

p_avg1 <- plot_avg_uncertainty(avg_all_limited_full, "Posterior mean INMB",
                               x_label = "Posterior mean INMB", true_df = true_df)

p_avg1

# ggsave("microsim_spline/MC_sim/avg_limited.png", width = 45, height = 45, unit = "cm")






pp_extended_15k <- combine_df_pp(pp_all_extended, INMB_true_20k_all)
# extended_17k <- combine_df_lpd(pp_all_extended_30k, INMB_true_30k_all)
# # # 
# # # 
# #
# 
library(kableExtra)
# library(webshot)
kable(pp_limited_15k, "latex", booktabs = TRUE, escape = FALSE) %>%
  kable_styling( full_width = F) %>%
  column_spec(2, extra_css = "font-size: 80%")

kable(pp_extended_15k, "latex", booktabs = TRUE, escape = FALSE) %>%
  kable_styling( full_width = F) %>%
  column_spec(2, extra_css = "font-size: 80%")
#   save_kable(file = "microsim_spline/MC_sim/limited_15k.html")
# webshot("microsim_spline/MC_sim/limited_15k.html", file = "microsim_spline/MC_sim/limited_15k.png")
# 
# #
# kable(limited_17k, "html", booktabs = TRUE, escape = FALSE) %>%
#   kable_styling(latex_options = "striped", full_width = F) %>%
#   column_spec(2, extra_css = "font-size: 80%")
# #
# kable(extended_15k, "html", booktabs = TRUE, escape = FALSE) %>%
#   kable_styling(latex_options = "striped", full_width = F) %>%
#   column_spec(2, extra_css = "font-size: 80%")
# 
# kable(extended_17k, "html", booktabs = TRUE, escape = FALSE) %>%
#   kable_styling(latex_options = "striped", full_width = F) %>%
#   column_spec(2, extra_css = "font-size: 80%")
# 
# 
# kable(limited_15k, "latex", booktabs = TRUE, escape = FALSE) %>%
#   kable_styling(latex_options = c("striped", "hold_position"), full_width = F) %>%
#   column_spec(2, width = "5cm")
# # 
# kable(extended_15k, "latex", booktabs = TRUE, escape = FALSE) %>%
#   kable_styling(latex_options = c("hold_position"), full_width = F) %>%
#   column_spec(2, width = "5cm")
# 
kable(elpd_limited_15k, "latex", booktabs = TRUE, escape = FALSE) %>%
  kable_styling(full_width = F) %>%
  column_spec(2, extra_css = "font-size: 80%")



generate_se_plots <- function(data, metric_col, plot_title, x_label, y_label) {
  # Validate input data structure
  stopifnot(ncol(data) >= 6, any(colnames(data) == "scenario"))
  
  # Helper function to create plot
  create_plot <- function(df, metric_col, SE, add_true, plot_title, x_label, y_label) {
    df_long <- melt(df, id.vars = "scenario")
    colnames(df_long)[colnames(df_long) == "variable"] <- "Methods"
    df_long$scenario <- str_wrap(df_long$scenario, width = 25)
    df_long[[metric_col]] <- as.numeric(sub("\\(.*\\)", "", df_long$value))
    df_long$SE <- as.numeric(sub(".*\\(", "", sub("\\).*", "", df_long$value)))
    
    # Assign a unique numeric ID to each scenario in df
    df_long$scenario_id <- as.numeric(factor(df_long$scenario, levels = unique(df_long$scenario)))
    
    # Create a new column for y-axis position with a small offset for each method in df
    method_offsets <- seq(from = -0.1, to = 0.1, length.out = length(unique(df_long$Methods)))
    names(method_offsets) <- unique(df_long$Methods)
    df_long$y_position <- df_long$scenario_id + method_offsets[df_long$Methods]
    

    p <- ggplot(df_long, aes(x = get(metric_col), y = y_position, shape = Methods, color = Methods)) +
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            axis.text.y = element_text(face = "bold")) +
      scale_y_continuous(labels = unique(df_long$scenario), breaks = unique(df_long$scenario_id)) +
      labs(x = x_label, y = y_label) +
      ggtitle(plot_title)
    
    
      p <- p + geom_pointrange(aes(xmin = get(metric_col) - 1.96*SE, xmax = get(metric_col) + 1.96*SE, group = Methods),
                               size = 0.5)
    
    return(p)
  }
  
  # Split into two datasets
  df1 <- data[1:6, ]
  df2 <- data[7:12, ]
  
  # Create plots
  p1 <- create_plot(df1, metric_col, SE, add_true, paste0(plot_title), x_label, y_label)
  p2 <- create_plot(df2, metric_col, SE, add_true, paste0(plot_title), x_label, y_label)
  
  # Combine plots
  comb_plot <- ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom")
  
  return(comb_plot)
}


lpd_plot <- generate_dot_plots(elpd_limited_15k, "LPD", "Log predictive density of true INMB", 
                               x_label = "Log Predictive Density (LPD)", y_label = "")
lpd_plot
# ggsave("microsim_spline/MC_sim/lpd_plot.png", width = 30, height = 25, unit = "cm")

generate_dot_plots <- function(data, metric_col, plot_title, x_label, y_label, jitter_width = 0.01) {
  # Validate input data structure
  stopifnot(ncol(data) >= 6, any(colnames(data) == "scenario"))
  
  # Helper function to create plot
  create_plot <- function(df, metric_col, plot_title, x_label, y_label) {
    df_long <- melt(df, id.vars = "scenario")
    colnames(df_long)[colnames(df_long) == "variable"] <- "Methods"
    df_long$scenario <- str_wrap(df_long$scenario, width = 25)
    df_long[[metric_col]] <- as.numeric(sub("\\(.*\\)", "", df_long$value))

    
    
    p <- ggplot(df_long, aes(x = get(metric_col), y = scenario, shape = Methods, color = Methods)) +
      geom_point(size = 4, position = position_dodge(width = jitter_width))+
      theme_minimal() +
      theme(plot.title = element_text(face = "bold"),
            axis.text.y = element_text(face = "bold")) +
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
  comb_plot <- ggarrange(p1, p2, ncol = 2, common.legend = TRUE, legend = "bottom")
  
  return(comb_plot)
}



pp_limited_15k_2 <- pp_limited_15k %>% dplyr::select(-INMB_true)
pp_plot <- generate_dot_plots(pp_limited_15k_2, "P_ce", "Probability of Cost-effectiveness, limited data",
                              x_label = expression("Probability of Cost-effectiveness, " * P[CE]), y_label = "")
pp_plot
# ggsave("microsim_spline/MC_sim/pp_limited.png", width = 40, height = 30, unit = "cm")


pp_extended_15k_2 <- pp_extended_15k %>% dplyr::select(-INMB_true)
pp_extended_plot <- generate_dot_plots(pp_extended_15k_2, "P_ce", "Probability of Cost-effectiveness, extended data",
                              x_label = expression("Probability of Cost-effectiveness, " * P[CE]), y_label = "",
                              jitter_width = 0.01)
pp_extended_plot
# ggsave("microsim_spline/MC_sim/pp_extended.png", width = 40, height = 30, unit = "cm")


