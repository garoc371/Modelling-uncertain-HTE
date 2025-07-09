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
wtp <- seq(0, 30000, by = 1000)

start_ages <- c(40,50,60,70,80)

# function to calculate INMB for each MC replication

compute_INMB <-function(sublist, wtp){
  INMB = (wtp * sublist$delta_QALYs - sublist$delta_costs)/n_target
  
  return(INMB)
}

pp_approval <- function(filepath, scenario, wtp){
  dat <- readRDS(paste0(filepath, scenario, ".RDs"))
  INMB_list <- map(dat, ~compute_INMB(., wtp))
  
  pp_list <- map(INMB_list, ~mean(.x>0))
  avg_pp <- mean(unlist(pp_list))
  
  
  return(avg_pp = avg_pp)
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

pp_all_limited <- vector(mode = "list", length = length(combinations))
pp_all_extended <- vector(mode = "list", length = length(combinations))
model_types <- c("unadjusted", "adjusted", "linear_interaction", "unrestricted_spline", "monotonic_spline")
INMB_list_all <- vector(mode = "list", length = length(combinations))



for (i in 1:nrow(combinations)) {
    control_type <- combinations$control_type[i]
    treatment_type <- combinations$treatment_type[i]
    
  
    # Shorten names
    control_short <- gsub("non_linear", "non", gsub("linear", "lin", as.character(control_type)))
    treatment_short <- gsub("non_linear", "non", gsub("linear", "lin", as.character(treatment_type)))
    
    scenario_name <- paste0(control_short, "_ctrl_", treatment_short, "_trt")
    
    true_weights_current <- readRDS(weight_names[[grep(scenario_name, weight_names)]])
    weights_current_g5 <- aggregate_weights(true_weights_current, start_ages)
    
    res_true_current <- readRDS(cohort_true_names[[grep(scenario_name, cohort_true_names)]])
    INMB_all_k <- vector(mode = "numeric", length = length(wtp))
    for (j in seq_along(wtp)) {
      INMB_true <- (wtp[j] * res_true_current$delta_QALYs - res_true_current$delta_costs)/n_target
      INMB_all_k[j] <- INMB_true
    }
    INMB_list_all[[i]] <- INMB_all_k
    pp_all_limited[[i]] <- list(scenario = scenario_name)
    pp_all_limited[[i]]$models <- setNames(vector("list", length(model_types)), model_types)
    
    pp_all_extended[[i]] <- list(scenario = scenario_name)
    pp_all_extended[[i]]$models <- setNames(vector("list", length(model_types)), model_types)
    
    for (model_type in model_types) {
      pp_temp_limited <- vector(mode = "numeric", length = length(wtp))
      pp_temp_extended <- vector(mode = "numeric", length = length(wtp))
      
      for (j in seq_along(wtp)) {
        pp_temp_limited[j] <- model_scenario_pp(scenario_name, model_type, "_limited", wtp[j])
        pp_temp_extended[j] <- model_scenario_pp(scenario_name, model_type, "_extended", wtp[j])
      }
      
      pp_all_limited[[i]]$models[[model_type]] <- pp_temp_limited
      pp_all_extended[[i]]$models[[model_type]] <- pp_temp_extended
      
    }
    
  }

combinations$scenario_name <- paste0(
  gsub("_", " ", combinations$control_type), " control, ",
  gsub("_", " ", combinations$treatment_type), " treatment"
)

combinations$scenario_abbr <- paste0(
  toupper(substr(gsub("_", "", combinations$control_type), 1, 2)),
  "_",
  toupper(substr(gsub("_", "", combinations$treatment_type), 1, 2))
)


ceac_data_limited <- vector(mode = "list", length = length(pp_all_limited))
ceac_data_extended <- vector(mode = "list", length = length(pp_all_extended))


for (i in seq_along(pp_all_limited)) {
  scenario_name <- pp_all_limited[[i]]$scenario
  model_data_limited <- pp_all_limited[[i]]$models
  model_data_extended <- pp_all_extended[[i]]$models
  
  
  ceac_data_limited[[i]] <- data.frame(
    wtp = rep(wtp, times = length(model_types)),
    probability = unlist(model_data_limited),
    model = rep(names(model_data_limited), each = length(wtp)),
    scenario = combinations$scenario_name[i] 
  ) 
  
  ceac_data_extended[[i]] <- data.frame(
    wtp = rep(wtp, times = length(model_types)),
    probability = unlist(model_data_extended),
    model = rep(names(model_data_extended), each = length(wtp)),
    scenario = combinations$scenario_name[i] 
  ) 
  
}

ceac_data_limited <- bind_rows(ceac_data_limited, .id = "scenario_index")
ceac_data_extended <- bind_rows(ceac_data_extended, .id = "scenario_index")

inmb_positive_wtp <- tibble(
  scenario_index = as.character(seq_along(INMB_list_all)),
  wtp_positive = map_dbl(INMB_list_all, function(INMB_vec) wtp[{which(INMB_vec>0)[1]}])
)

ceac_data_limited <- ceac_data_limited %>% 
  left_join(inmb_positive_wtp, by = "scenario_index")

ceac_data_extended <- ceac_data_extended %>% 
  left_join(inmb_positive_wtp, by = "scenario_index")


ggplot(ceac_data_limited, aes(x = wtp, y = probability)) +
  geom_line(aes(color = model)) +
  geom_vline(aes(xintercept = wtp_positive), linetype = "dashed", color = "red") +
  facet_wrap(~ scenario, ncol = 4, labeller = labeller(scenario = label_wrap_gen(width = 40)))   +
  labs(x = "Willingness-to-Pay", y = "Probability of Cost-Effectiveness",
       title = "Cost-Effectiveness Acceptability Curves") +
  theme_bw() +
  guides(color = guide_legend(title = "Model Type")) 

ggsave("microsim_spline/MC_sim/ceac_all_limited2.png", height = 12, width = 9)

ggplot(ceac_data_extended, aes(x = wtp, y = probability)) +
  geom_line(aes(color = model)) +
  geom_vline(aes(xintercept = wtp_positive), linetype = "dashed", color = "red") +
  facet_wrap(~ scenario, ncol = 4, labeller = labeller(scenario = label_wrap_gen(width = 40)))   +
  labs(x = "Willingness-to-Pay", y = "Probability of Cost-Effectiveness",
       title = "Cost-Effectiveness Acceptability Curves") +
  theme_bw() +
  guides(color = guide_legend(title = "Model Type")) 
ggsave("microsim_spline/MC_sim/ceac_all_extended2.png", height = 6, width = 10, dpi = 600)


  
