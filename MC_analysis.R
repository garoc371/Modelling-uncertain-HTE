library(tidyverse)
library(stringr)
library(patchwork)
library(ggpubr)
library(here)
# library(Rcpp)
source(here("simulation_functions.R"))
source(here("analysis_functions.R"))
library(ggdist)
library(reshape2)


weight_names <- list.files(
  path = here("results", "true"),
  pattern = "*weights_true.RDs",
  full.names = TRUE
)


cohort_true_names <- list.files(
  path = here("results", "true"),
  pattern = "*cohort_true.RDs",
  full.names = TRUE
)

n_target <- 20000
wtp_1 <- 15000


start_ages <- c(40, 50, 60, 70, 80)

# Functions moved to analysis_functions.R

combinations <- expand.grid(
  control_type = c(
    "non_linear_mon_increase",
    "non_linear_mon_decrease",
    "non_linear_non_mon"
  ),
  treatment_type = c(
    "constant",
    "non_linear_mon_increase",
    "non_linear_mon_decrease",
    "non_linear_non_mon"
  )
)

# for wtp = 15k
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


model_types <- c(
  "unadjusted",
  "adjusted",
  "linear_interaction",
  "unrestricted_spline",
  "monotonic_spline"
)
INMB_true_15k_all <- list()

for (i in 1:nrow(combinations)) {
  control_type <- combinations$control_type[i]
  treatment_type <- combinations$treatment_type[i]

  pp_limited <- pp_extended <-
    setNames(vector("list", length(model_types)), model_types)

  avg_INMB_limited <- avg_INMB_extended <-
    CI_INMB_limited <- CI_INMB_extended <-
      setNames(vector("list", length(model_types)), model_types)

  elpd_limited <- elpd_extended <- elpd_limited_se <- elpd_extended_se <-
    setNames(vector("list", length(model_types)), model_types)

  # Shorten names
  control_short <- gsub(
    "non_linear",
    "non",
    gsub("linear", "lin", as.character(control_type))
  )
  treatment_short <- gsub(
    "non_linear",
    "non",
    gsub("linear", "lin", as.character(treatment_type))
  )

  scenario_name <- paste0(control_short, "_ctrl_", treatment_short, "_trt")

  true_weights_current <- readRDS(weight_names[[grep(
    scenario_name,
    weight_names
  )]])
  weights_current_g5 <- aggregate_weights(true_weights_current, start_ages)

  res_true_current <- readRDS(cohort_true_names[[grep(
    scenario_name,
    cohort_true_names
  )]])
  INMB_true_15k <- (wtp_1 *
    res_true_current$delta_QALYs -
    res_true_current$delta_costs) /
    n_target

  INMB_true_15k_all[[i]] <- INMB_true_15k

  for (model_type in model_types) {
    elpd_temp_limited <- model_scenario_elpd(
      scenario_name,
      model_type,
      "_limited",
      wtp_1,
      INMB_true_15k
    )
    elpd_temp_extended <- model_scenario_elpd(
      scenario_name,
      model_type,
      "_extended",
      wtp_1,
      INMB_true_15k
    )

    elpd_limited[[model_type]] <- elpd_temp_limited$elpd_INMB
    elpd_limited_se[[model_type]] <- elpd_temp_limited$se_lpd

    elpd_extended[[model_type]] <- elpd_temp_extended$elpd_INMB
    elpd_extended_se[[model_type]] <- elpd_temp_extended$se_lpd
  }

  elpd_all_limited[[i]] <- elpd_limited
  elpd_all_extended[[i]] <- elpd_extended

  elpd_all_limited_se[[i]] <- elpd_limited_se
  elpd_all_extended_se[[i]] <- elpd_extended_se

  for (model_type in model_types) {
    pp_temp_limited <- model_scenario_pp(
      scenario_name,
      model_type,
      "_limited",
      wtp_1
    )
    pp_temp_extended <- model_scenario_pp(
      scenario_name,
      model_type,
      "_extended",
      wtp_1
    )

    pp_limited[[model_type]] <- pp_temp_limited$avg_pp

    pp_extended[[model_type]] <- pp_temp_extended$avg_pp
  }

  pp_all_limited[[i]] <- pp_limited
  pp_all_extended[[i]] <- pp_extended

  for (model_type in model_types) {
    temp_INMB_limited <- model_scenario_avg(
      scenario_name,
      model_type,
      "_limited",
      wtp_1
    )
    temp_INMB_extended <- model_scenario_avg(
      scenario_name,
      model_type,
      "_extended",
      wtp_1
    )

    avg_INMB_limited[[model_type]] <- temp_INMB_limited$average
    CI_INMB_limited[[model_type]] <- list(
      LL = temp_INMB_limited$ll_INMB,
      UL = temp_INMB_limited$ul_INMB
    )
    avg_INMB_extended[[model_type]] <- temp_INMB_extended$average
    CI_INMB_extended[[model_type]] <- list(
      LL = temp_INMB_extended$ll_INMB,
      UL = temp_INMB_extended$ul_INMB
    )
  }

  avg_INMB_all_limited[[i]] <- avg_INMB_limited
  avg_INMB_all_extended[[i]] <- avg_INMB_extended
  CI_INMB_all_limited[[i]] <- CI_INMB_limited
  CI_INMB_all_extended[[i]] <- CI_INMB_extended
}


combinations$scenario_name <- paste0(
  gsub("_", " ", combinations$control_type),
  " control ",
  gsub("_", " ", combinations$treatment_type),
  " treatment"
)


elpd_limited_15k <- combine_df_lpd(elpd_all_limited, elpd_all_limited_se)


true_df <- tibble(
  INMB = unlist(INMB_true_15k_all),
  scenario = combinations$scenario_name
)


# kable(limited_15k, "latex", booktabs = TRUE, escape = FALSE) %>%
#   kable_styling(latex_options = c("striped", "hold_position"), full_width = F) %>%
#   column_spec(2, width = "5cm")

#
pp_limited_15k <- combine_df_pp(pp_all_limited, INMB_true_15k_all)

CI_all_limited <- extract_info(CI_INMB_all_limited, combinations$scenario_name)
avg_all_limited <- avg_INMB_all_limited %>%
  bind_rows() %>%
  add_column(scenario = combinations$scenario_name, .before = 1) %>%
  melt(., id.vars = "scenario")


colnames(avg_all_limited)[colnames(avg_all_limited) == "variable"] <- "Methods"
avg_all_limited_full <- left_join(
  avg_all_limited,
  CI_all_limited,
  by = c("scenario", "Methods")
) %>%
  mutate_at(vars(last_col():(last_col() - 2)), round, digits = 1)


p_avg1 <- plot_avg_uncertainty(
  avg_all_limited_full,
  "Posterior mean INMB",
  x_label = "Posterior mean INMB",
  true_df = true_df
)

p_avg1

# ggsave(here("figures", "avg_limited.png"), width = 45, height = 45, unit = "cm")

pp_extended_15k <- combine_df_pp(pp_all_extended, INMB_true_15k_all)

#
library(kableExtra)
# library(webshot)
kable(pp_limited_15k, "latex", booktabs = TRUE, escape = FALSE) %>%
  kable_styling(full_width = F) %>%
  column_spec(2, extra_css = "font-size: 80%")

kable(pp_extended_15k, "latex", booktabs = TRUE, escape = FALSE) %>%
  kable_styling(full_width = F) %>%
  column_spec(2, extra_css = "font-size: 80%")

kable(elpd_limited_15k, "latex", booktabs = TRUE, escape = FALSE) %>%
  kable_styling(full_width = F) %>%
  column_spec(2, extra_css = "font-size: 80%")

lpd_plot <- generate_se_plots(
  elpd_limited_15k,
  "LPD",
  "Log predictive density of true INMB",
  x_label = "Log Predictive Density (LPD)",
  y_label = ""
)
lpd_plot
# ggsave(here("figures", "lpd_plot.png"), width = 12, height = 9, dpi = 600)

pp_limited_15k_2 <- pp_limited_15k %>% dplyr::select(-INMB_true)
pp_plot <- generate_dot_plots(
  pp_limited_15k_2,
  "P_ce",
  "Probability of Cost-effectiveness, limited data",
  x_label = expression("Probability of Cost-effectiveness, " * P[CE]),
  y_label = ""
)
pp_plot
# ggsave(here("figures", "pp_limited.png"), width = 10, height = 15, dpi = 600)

pp_extended_15k_2 <- pp_extended_15k %>% dplyr::select(-INMB_true)
pp_extended_plot <- generate_dot_plots(
  pp_extended_15k_2,
  "P_ce",
  "Probability of Cost-effectiveness, extended data",
  x_label = expression("Probability of Cost-effectiveness, " * P[CE]),
  y_label = "",
  jitter_width = 0.01
)
pp_extended_plot
# ggsave(here("figures", "pp_extended.png"), width = 10, height = 15, dpi = 600)
