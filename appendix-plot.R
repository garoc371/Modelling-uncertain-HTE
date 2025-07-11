# NOTICE: ggplot2 version 3.4.x, ggpubr 0.6.0,
# cowplot 1.1.3, and ggdist 3.3.1 are required
# to reproduce the plot. ggplot2 3.5.x introduces
# a computed mapping layer that is incompatible
# with the stored ggplot object

library(ggpubr)
library(here)
library(Cairo)
library(tidyverse)
library(ggdist)
library(patchwork)
source(here("simulation_functions.R"))

# relevant parameters
n_target <- 2e4
start_age <- seq(40, 80, by = 10)

# read in files
pop_weights <- readRDS(paste0(
  "appendix_data/",
  list.files(here("data"), pattern = "*weights")
))
cohort_true <- readRDS(paste0(
  "appendix_data/",
  list.files(here("data"), pattern = "*cohort_true")
))
res_g5 <- readRDS(paste0(
  "appendix_data/",
  list.files(here("data"), pattern = "*cohort_all")
))
linear_plots <- readRDS(paste0(
  "appendix_data/",
  list.files(here("data"), pattern = "*linear")
))
adjusted_plots <- readRDS(paste0(
  "appendix_data/",
  list.files(here("data"), pattern = "*adjusted")
))
unrestricted_splines <- readRDS(paste0(
  "appendix_data/",
  list.files(here("data"), pattern = "*unrestricted_spline")
))
monotonic_splines <- readRDS(paste0(
  "appendix_data/",
  list.files(here("data"), "monotonic_spline")
))


# relevant helper functions
plot_cohort <- function(res_cohort, res_true) {
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

get_y_limits <- function(plot) {
  data_list <- ggplot_build(plot)$data
  y_range_list <- lapply(data_list, function(layer_data) {
    range(layer_data$y, na.rm = TRUE)
  })
  do.call("range", y_range_list)
}

set_global_ylimits <- function(plot) {
  plot + scale_y_continuous(limits = global_y_range)
}

###################Cost-effectiveness plane################
# relevant plot descriptions

scenario_cohort <- c(
  "unadjusted_limited",
  "unadjusted_extended",
  "standard_limited",
  "standard_extended",
  "linear_limited",
  "linear_extended",
  "unrestricted_spline_limited",
  "unrestricted_spline_extended",
  "Monotonic_spline_limited",
  "Monotonic_spline_extended"
)

control_description <- "non linear mon decrease control surface"
treatment_description <- "non linear mon decrease treatment surface"

annotion_title <- paste0(control_description, "\n", treatment_description)

cohort_true_avg <- lapply(cohort_true, function(y) y / n_target)
pop_weights_g5 <- aggregate_weights(pop_weights, start_age)

cohort_g5_avg <- multi_cohort_average(res_g5, pop_weights_g5) %>%
  lapply(., function(method_results) {
    lapply(method_results, function(y) y / n_target)
  }) %>%
  map(., bind_rows)

CE_g5_list <- map(cohort_g5_avg, ~ plot_cohort(.x, cohort_true_avg))
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

CE_g5_plots <- ggarrange(
  plotlist = CE_g5_list,
  nrow = 5,
  ncol = 2,
  common.legend = T,
  legend = "bottom"
)

CE_g5_plots <- annotate_figure(
  CE_g5_plots,
  top = text_grob(
    paste0(annotion_title, "\n5age groups"),
    x = 0,
    hjust = 0,
    face = "bold"
  )
)

#########################Outcome surface plot######################

all_y_limits <- c(
  sapply(linear_plots, get_y_limits),
  sapply(unrestricted_splines, get_y_limits),
  sapply(monotonic_splines, get_y_limits),
  sapply(adjusted_plots, get_y_limits)
)

y_max <- max(all_y_limits, na.rm = TRUE) + 0.2
y_min <- min(all_y_limits, na.rm = TRUE)
y_min <- ifelse(y_min - 0.2 < 0, 0, y_min - 0.2)

global_y_range <- c(y_min, y_max)

linear_plots <- lapply(linear_plots, set_global_ylimits)
unrestricted_splines <- lapply(
  unrestricted_splines,
  set_global_ylimits
)
monotonic_splines <- lapply(monotonic_splines, set_global_ylimits)
adjusted_plots <- lapply(adjusted_plots, set_global_ylimits)

outcome_plots <- ggarrange(
  plotlist = c(
    adjusted_plots,
    linear_plots,
    unrestricted_splines,
    monotonic_splines
  ),
  nrow = 4,
  ncol = 2,
  common.legend = TRUE,
  legend = "bottom"
)
outcome_title <- paste0(
  "Outcome surface under limited/extended data",
  "\n",
  "Method from top to bottom: adjusted regressions,linear interaction, unrestricted spline, monotonic spline"
)

outcome_plots <- annotate_figure(
  outcome_plots,
  top = text_grob(outcome_title, hjust = 0, x = 0, face = "bold")
)

CE_g5_plots | outcome_plots
ggsave(
  "figures/CE_outcome_plots.png",
  width = 40,
  height = 30,
  units = "cm",
  dpi = 600
)
# for vertical plots, simply adjust the layout in `ggarrange` call
