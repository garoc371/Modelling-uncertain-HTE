library(tidyverse)
library(stringr)
library(patchwork)
library(ggpubr)
library(Rcpp)
source("microsim_spline/cohortsim_helper.R")
source("microsim_spline/microsim_helper.R")
library(ggdist)

# library(ggblend)
library(Cairo)
theme_set(theme_ggdist())

weight_names <- list.files(path = "./microsim_spline/res_1123", 
                           pattern = "*weights_true.RDs", 
                           full.names = TRUE)


cohort_true_names <- list.files(path = "./microsim_spline/res_1123", 
                         pattern = "*cohort_true.RDs", 
                         full.names = TRUE)

res_g5_names <- list.files(path = "./microsim_spline/res_1123", 
                           pattern = "*cohort_all_g5.RDs", 
                           full.names = TRUE)

# res_g10_names <- list.files(path = "./microsim_spline/res_1024", 
#                             pattern = "*cohort_all_g10.RDs", 
#                             full.names = TRUE)

# res_full_names <- list.files(path = "./microsim_spline/res_104", 
#                              pattern = "*cohort_all_full.RDs", 
#                              full.names = TRUE)

plots_names <- list.files(path = "./microsim_spline/plot_1123", 
                          pattern = "*linear.RDs|*adjusted.RDs|*spline.RDs|", 
                          full.names = TRUE)


n_target <- 20000

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
time_horizon <-40


n_simulations <- 4000
n_patients <- n_target
n_cycles <- 40

baseline_probs = 0.3

utilities <- c(utility_healthy, utility_diseased, 0)
cost_ctrl <- c(cost_healthy, cost_diseased,0)
cost_trt <- c(cost_healthy+treatment_cost_per_cycle , cost_diseased, 0)

start_age1 <- seq(40,80, by = 10)
start_age2 <- seq(40, 80, by = 5)

combinations <- expand.grid(
  control_type = c("linear_increasing", "linear_decreasing", 
                   "non_linear_mon_increase", "non_linear_mon_decrease", "non_linear_non_mon"),
  treatment_type = c("constant", "linear_increasing", "linear_decreasing", 
                     "non_linear_mon_increase", "non_linear_mon_decrease", "non_linear_non_mon")
)

# weight_list <- multi_file_import(weight_names)
# res_true_list <- multi_file_import(res_true_names)
# res_g5_list <- multi_file_import(res_g5_names)
# res_g10_list <- multi_file_import(res_g10_names)
# res_full_list <- multi_file_import(res_full_names)
# 
# res_true_avg <- map(res_true_list, ~map(.x, ~.x/n_target))
# cohort_g5_avg <- map(res_g5_list, ~map_depth(.x, 3, ~.x/n_target))
# cohort_g10_avg <- map(res_g10_list, ~map_depth(.x, 3, ~.x/n_target))
# cohort_full_avg <- map(res_full_list, ~map_depth(.x, 3, ~.x/n_target))


# starting age for five cohort simulations
scenario_cohort <- c("unadjusted_limited","unadjusted_extended",
                     "standard_limited", "standard_extended",
                     "linear_limited", "linear_extended",
                     "unrestricted_spline_limited", "unrestricted_spline_extended",
                     "Monotonic_spline_limited", "Monotonic_spline_extended" )

plot_cohort <- function(res_cohort, res_true){
  res_cohort <- res_cohort %>% 
    mutate(source = rep("cohort_sims", length(delta_QALYs)))
  
  cep_data_true <- data.frame(
    delta_QALYs = res_true$delta_QALYs,
    delta_costs = res_true$delta_costs,
    source = rep("true_value", length(res_true$delta_QALYs))
  )
  
  cep_data <- rbind(res_cohort, cep_data_true)
  
  # Generate plot
  ce_plot <- plot_ce(cep_data)
  
  return(list(plot = ce_plot))
  
}

for (i in 1:nrow(combinations)) {
  control_type <- combinations$control_type[i]
  treatment_type <- combinations$treatment_type[i]
  
  # Shorten names
  control_short <- gsub("non_linear", "non", gsub("linear", "lin", as.character(control_type)))
  treatment_short <- gsub("non_linear", "non", gsub("linear", "lin", as.character(treatment_type)))
  
  control_description <- gsub("_", " ", control_type)
  treatment_description <- gsub("_", " ", treatment_type)
  
  annotation_title <- paste0(
    control_description, " control surface, \n", treatment_description, " treatment surface"
  )
  
  # Scenario name
  scenario_name <- paste0(control_short, "_ctrl_", treatment_short, "_trt")
  
  true_weights_current <- readRDS(weight_names[[grep(scenario_name, weight_names)]])
  weights_current_g5 <- aggregate_weights(true_weights_current, start_age1)
  # weights_current_g10 <- aggregate_weights(true_weights_current, start_age2)
  
  res_true_current <- readRDS(cohort_true_names[[grep(scenario_name, cohort_true_names)]])
  res_true_avg <- lapply(res_true_current, function(y) y/n_target)
  res_g5_current <- read_file_avg(res_g5_names, scenario_name, n_target, weights_current_g5)
  # res_g10_current <- read_file_avg(res_g10_names, scenario_name, n_target, weights_current_g10)
  # res_full_current <- read_file_avg(res_full_names, scenario_name, n_target, true_weights_current)
  
  CE_g5_list <- map(res_g5_current, ~plot_cohort(.x, res_true_avg))
  plot_cohort_list <- map(CE_g5_list, extract_plot_data)
  all_delta_qalys <- unlist(map(plot_cohort_list, "delta_QALYs"))
  all_delta_costs <- unlist(map(plot_cohort_list, "delta_costs"))
  rm(plot_cohort_list)
  range_delta_qalys <- range(all_delta_qalys)
  range_delta_costs <- range(all_delta_costs)
  rm(all_delta_costs, all_delta_qalys)
  for (j in 1:length(scenario_cohort)) {
      CE_g5_list[[j]] <- CE_g5_list[[j]]$plot +
      scale_x_continuous(limits = range_delta_qalys) +
      scale_y_continuous(limits = range_delta_costs) +
      ggtitle(scenario_cohort[[j]])
  }

  CE_g5_plots <- ggarrange(plotlist = CE_g5_list, nrow = 5, ncol = 2, common.legend = T,
                           legend = "bottom")
  CE_g5_plots <-annotate_figure(
    CE_g5_plots,
    top = text_grob(paste0(annotation_title, "\n5 age groups"),x = 0, hjust = 0, face = "bold")
  )

  outcome_plot <- import_plots(plots_names = plots_names, scenario_name)[["outcome_plots"]]

  CE_g5_plots | outcome_plot
  ggsave(paste0("microsim_spline/plot_ce1124/CE_g5_outcome_",scenario_name, ".png"),
         width = 45, height = 30, units = "cm")
  # rm(CE_g5_list, CE_g5_plots)
  
  INMB_20k <- plot_INMB_multi(res_g5_current, res_true_avg, 20000)
  INMB_30k <- plot_INMB_multi(res_g5_current, res_true_avg, 30000)
  
  combined_INMB <- ggarrange(INMB_20k, INMB_30k, nrow = 1, ncol = 2, common.legend = T,
                             legend = "bottom")
  combined_INMB <-annotate_figure(
    combined_INMB,
    top = text_grob(annotation_title,x = 0, hjust = 0, face = "bold")
  )
  combined_INMB <- combined_INMB + theme(
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key = element_rect(fill = "transparent", colour = NA)
  )
  
  
  combined_INMB  
  ggsave(paste0("microsim_spline/plot_ce1124/INMB_",scenario_name, ".png"),
         width = 45, height = 30, units = "cm", device = "png", bg = "white")
  
  
  # CE_g10_list <- map(res_g10_current, ~plot_cohort(.x, res_true_avg))
  # plot_cohort_list <- map(CE_g10_list, extract_plot_data)
  # all_delta_qalys <- unlist(map(plot_cohort_list, "delta_QALYs"))
  # all_delta_costs <- unlist(map(plot_cohort_list, "delta_costs"))
  # rm(plot_cohort_list)
  # range_delta_qalys <- range(all_delta_qalys)
  # range_delta_costs <- range(all_delta_costs)
  # rm(all_delta_costs, all_delta_qalys)
  # 
  # for (j in 1:length(scenario_cohort)) {
  #   CE_g10_list[[j]] <- CE_g10_list[[j]]$plot +
  #     scale_x_continuous(limits = range_delta_qalys) +
  #     scale_y_continuous(limits = range_delta_costs) +
  #     ggtitle(scenario_cohort[[j]])
  #   
  # }
  
  # CE_g10_plots <- ggarrange(plotlist = CE_g10_list, nrow = 4, ncol = 2, common.legend = T,
  #                          legend = "bottom")
  # CE_g10_plots <-annotate_figure(
  #   CE_g10_plots,
  #   top = text_grob(paste0(annotation_title, "\n10 age groups"),x = 0, hjust = 0, face = "bold")
  # )
  
  # CE_g10_plots | outcome_plot
  # ggsave(paste0("microsim_spline/plot_ce1024/CE_g10_outcome_",scenario_name, ".png"),
  #        width = 45, height = 30, units = "cm")
  # rm(CE_g10_list, CE_g10_plots)
  
  
  # CE_full_list <- map(res_full_current, ~plot_cohort(.x, res_true_avg))
  # plot_cohort_list <- map(CE_full_list, extract_plot_data)
  # all_delta_qalys <- unlist(map(plot_cohort_list, "delta_QALYs"))
  # all_delta_costs <- unlist(map(plot_cohort_list, "delta_costs"))
  # rm(plot_cohort_list)
  # range_delta_qalys <- range(all_delta_qalys)
  # range_delta_costs <- range(all_delta_costs)
  # rm(all_delta_costs, all_delta_qalys)
  # 
  # for (j in 1:length(scenario_cohort)) {
  #   CE_full_list[[j]] <- CE_full_list[[j]]$plot +
  #     scale_x_continuous(limits = range_delta_qalys) +
  #     scale_y_continuous(limits = range_delta_costs) +
  #     ggtitle(scenario_cohort[[j]])
  # }
  # 
  # CE_full_plots <- ggarrange(plotlist = CE_full_list, nrow = 4, ncol = 2, common.legend = T,
  #                          legend = "bottom")
  # CE_full_plots <-annotate_figure(
  #   CE_full_plots,
  #   top = text_grob(paste0(annotation_title, "\nFull age groups"),x = 0, hjust = 0, face = "bold")
  # )
  # 
  # 
  # CE_full_plots | outcome_plot
  # ggsave(paste0("microsim_spline/plot_ce104/CE_full_outcome_",scenario_name, ".png"),
  #        width = 45, height = 30, units = "cm")
  # rm(CE_full_list, CE_full_plots)
  # 
  
}  
