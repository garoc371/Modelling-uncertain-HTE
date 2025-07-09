
gen_trial_data <- function(n_trial, n_extra, probs_control, probs_treatment, 
                           probs_control_extra, probs_treatment_extra, ages_control,
                           ages_treatment, extra_control, extra_trt){
  age_grid <- 40:80
  
  control_data <- tibble(probs = probs_control,
                         age = ages_control,
                         trt = rep(0, n_trial))
  control_data$outcome <- rbinom(n_trial,1, control_data$probs)
  
  treatment_data <- tibble(probs = probs_treatment,
                           age = ages_treatment,
                           trt = rep(1, n_trial))
  treatment_data$outcome <- rbinom(n_trial,1, treatment_data$probs)
  
  trial_data <- rbind(control_data, treatment_data)
  trial_data$age_std <- as.vector(scale(trial_data$age))
  
  control_data_extra <- tibble(probs = probs_control_extra,
                               age = extra_control,
                               trt = rep(0, n_extra))
  control_data_extra$outcome <- rbinom(n_extra,1, control_data_extra$probs)
  control_data_extra <- rbind(control_data,control_data_extra)
  
  
  treatment_data_extra <- tibble(probs = probs_treatment_extra,
                                 age = extra_trt,
                                 trt = rep(1, n_extra))
  treatment_data_extra$outcome <- rbinom(n_extra,1, treatment_data_extra$probs)
  treatment_data_extra <- rbind(treatment_data, treatment_data_extra)
  
  trial_data_extra <- rbind(control_data_extra, treatment_data_extra)
  trial_data_extra$age_std <- as.vector(scale(trial_data_extra$age))
  
  return(list(trial_data = trial_data, trial_data_extra = trial_data_extra))

}  

adjusted_model <- function(trial_data, trial_data_extra, filename, plot_outcome = F){
    
  adjusted_model <- stan_glm(
    outcome ~  trt + age_std,
    data = trial_data,
    family = binomial(link = "logit"),
    prior = normal(0, 2.5, autoscale = T),
    prior_intercept = normal(0, 2.5, autoscale = T),
    chains = 4,
    iter = 2000
  )
  
  adjusted_model_ex <- stan_glm(
    outcome ~  trt + age_std,
    data = trial_data_extra,
    family = binomial(link = "logit"),
    prior = normal(0, 2.5, autoscale = T),
    prior_intercept = normal(0, 2.5, autoscale = T),
    chains = 4,
    iter = 2000
  )
  
  if (plot_outcome){
    
    n_samples = 50
    
    true_outcome_control <- outcome_data(target_probs_control, "Control", age_target) %>% mutate(sample = "True")
    true_outcome_trt <- outcome_data(target_probs_treatment, "Treatment", age_target) %>% mutate(sample = "True")
    
    adjusted_plots <- vis_outcome_effect(true_outcome_control, true_outcome_trt,
                                         adjusted_model, adjusted_model_ex, age_grid)
    
    saveRDS(adjusted_plots, paste0(filename, "_adjusted.RDs"))
    rm(adjusted_plots)
    gc()
  }

 
  
  
  adjusted_tp_ctrl <- posterior_transition_probs(adjusted_model, trt_status = 0, age_varying = T, standardize_age = T)
  adjusted_tp_trt <- posterior_transition_probs(adjusted_model, trt_status = 1, age_varying = T, standardize_age = T)
  
  res_adjusted_multi_g5 <- run_multi_cohort_sequential(start_ages = start_ages[[1]], adjusted_tp_ctrl,
                                                     adjusted_tp_trt, n_cycles = n_cycles, utilities = utilities,
                                                     cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  adjusted_tp_ctrl_ex <- posterior_transition_probs(adjusted_model_ex, trt_status = 0, age_varying = T, standardize_age = T)
  adjusted_tp_trt_ex <- posterior_transition_probs(adjusted_model_ex, trt_status = 1, age_varying = T, standardize_age = T)
  
  res_adjusted_multi_ex_g5 <- run_multi_cohort_sequential(start_ages = start_ages[[1]], adjusted_tp_ctrl_ex,
                                                        adjusted_tp_trt_ex, n_cycles = n_cycles, utilities = utilities,
                                                        cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  return(list(adjusted_limited = res_adjusted_multi_g5, adjusted_extended = res_adjusted_multi_ex_g5))
  
}  
    
linear_interaction_model <- function(trial_data, trial_data_extra, filename, plot_outcome = F){

  bayesian_model <- stan_glm(
    outcome ~  trt*age_std,
    data = trial_data,
    family = binomial(link = "logit"),
    prior = normal(0, 2.5, autoscale = T),
    prior_intercept = normal(0, 2.5, autoscale = T),
    chains = 4,
    iter = 2000
  )
  
  bayesian_model_ex <- stan_glm(
    outcome ~  trt*age_std,
    data = trial_data_extra,
    family = binomial(link = "logit"),
    prior = normal(0, 2.5, autoscale = T),
    prior_intercept = normal(0, 2.5, autoscale = T),
    chains = 4,
    iter = 2000
  )
  
  if (plot_outcome){
    
    n_samples = 50
    
    true_outcome_control <- outcome_data(target_probs_control, "Control", age_target) %>% mutate(sample = "True")
    true_outcome_trt <- outcome_data(target_probs_treatment, "Treatment", age_target) %>% mutate(sample = "True")
    
    linear_plots <- vis_outcome_effect(true_outcome_control, true_outcome_trt,
                                       bayesian_model, bayesian_model_ex, age_grid)
    saveRDS(linear_plots, paste0(filename, "_linear.RDs"))
    rm(linear_plots)
    gc()
  }

  linear_tp_ctrl <- posterior_transition_probs(bayesian_model, trt_status = 0, age_varying = T, standardize_age = T)
  linear_tp_trt <- posterior_transition_probs(bayesian_model, trt_status = 1, age_varying = T, standardize_age = T)
  
  res_linear_multi_g5 <- run_multi_cohort_sequential(start_ages = start_ages[[1]], linear_tp_ctrl,
                                                   linear_tp_trt, n_cycles = n_cycles, utilities = utilities,
                                                   cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  linear_tp_ctrl_ex <- posterior_transition_probs(bayesian_model_ex, trt_status = 0, age_varying = T, standardize_age = T)
  linear_tp_trt_ex <- posterior_transition_probs(bayesian_model_ex, trt_status = 1, age_varying = T, standardize_age = T)
  
  res_linear_multi_ex_g5 <- run_multi_cohort_sequential(start_ages = start_ages[[1]], linear_tp_ctrl_ex,
                                                      linear_tp_trt_ex, n_cycles = n_cycles, utilities = utilities,
                                                      cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  return(list(linear_limited = res_linear_multi_g5, linear_extended = res_linear_multi_ex_g5))
  
} 
  
  
  
  
  
unrestricted_spline_model <- function(stan_model, trial_data, trial_data_extra, filename, plot_outcome = F){  
  
  file_name <- paste0(filename, "_unrestricted_spline.RDs")
  
  control_data <- trial_data %>% filter(trt == 0)
  treatment_data <- trial_data %>% filter(trt == 1)
  control_data_extra <- trial_data_extra %>% filter(trt == 0)
  treatment_data_extra <- trial_data_extra %>% filter(trt == 1)
  bspline_rw <- cmdstan_fit_plot(mod =stan_model, increment = 0.15, target_age = age_grid, degree = 3,
                                 ctrl_dat = control_data, ctrl_dat_ex = control_data_extra,
                                 trt_dat = treatment_data, trt_dat_ex = treatment_data_extra,
                                 target_probs_control = target_probs_control, target_probs_treatment = target_probs_treatment,
                                 ispline = F, plot_outcome = plot_outcome, file_name = file_name)

  
  bspline_tp_ctrl <- posterior_spline_matrix(bspline_rw$fit_ctrl,bspline_rw$ctrl_spline_basis, newx = age_grid) 
  bspline_tp_ctrl_ex <- posterior_spline_matrix(bspline_rw$fit_ctrl_ex, bspline_rw$ex_ctrl_spline_basis, newx = age_grid)
  bspline_tp_trt <- posterior_spline_matrix(bspline_rw$fit_trt, bspline_rw$trt_spline_basis, newx = age_grid)
  bspline_tp_trt_ex <- posterior_spline_matrix(bspline_rw$fit_trt_ex, bspline_rw$ex_trt_spline_basis, newx = age_grid)
  
  res_bspline_multi_g5 <- run_multi_cohort_sequential(start_ages = start_ages[[1]], bspline_tp_ctrl,
                                                    bspline_tp_trt, n_cycles = n_cycles, utilities = utilities,
                                                    cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  
  
  res_bspline_multi_ex_g5 <- run_multi_cohort_sequential(start_ages = start_ages[[1]], bspline_tp_ctrl_ex,
                                                       bspline_tp_trt_ex, n_cycles = n_cycles, utilities = utilities,
                                                       cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  return(list(unres_spline_limited = res_bspline_multi_g5, unres_spline_extended = res_bspline_multi_ex_g5))
  
}


unadjusted_model <- function(trial_data, trial_data_extra){  
  unadjusted_model <- stan_glm(
    outcome ~  trt,
    data = trial_data,
    family = binomial(link = "logit"),
    prior = normal(0, 2.5, autoscale = T),
    prior_intercept = normal(0, 2.5, autoscale = T),
    chains = 4,
    iter = 2000
  )
  
  unadjusted_model_ex <- stan_glm(
    outcome ~  trt,
    data = trial_data_extra,
    family = binomial(link = "logit"),
    prior = normal(0, 2.5, autoscale = T),
    prior_intercept = normal(0, 2.5, autoscale = T),
    chains = 4,
    iter = 2000
  )
  
# calculate the transition probability for unadjusted model
# and run the mutli-cohort simulation loop
  
  unadjusted_tp_ctrl <- posterior_transition_probs(unadjusted_model, trt_status = 0, age_varying = F)
  unadjusted_tp_trt <- posterior_transition_probs(unadjusted_model, trt_status = 1, age_varying = F)
  
  res_unadjusted_multi_g5 <- run_multi_cohort_sequential(start_ages = start_ages[[1]], unadjusted_tp_ctrl,
                                                    unadjusted_tp_trt, n_cycles = n_cycles, utilities = utilities,
                                                    cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  
  
  # res_unadjusted_multi_g10 <- run_multi_cohort_parallel(start_ages = start_ages[[2]], unadjusted_tp_ctrl,
  #                                                      unadjusted_tp_trt, n_cycles = n_cycles, utilities = utilities,
  #                                                      cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  gc()
  
  # res_unadjusted_multi_all <- run_multi_cohort_parallel(start_ages = start_ages[[3]], unadjusted_tp_ctrl,
  #                                                      unadjusted_tp_trt, n_cycles = n_cycles, utilities = utilities,
  #                                                      cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  # 
  # gc()
  
  unadjusted_tp_ctrl_ex <- posterior_transition_probs(unadjusted_model_ex, trt_status = 0, age_varying = F)
  unadjusted_tp_trt_ex <- posterior_transition_probs(unadjusted_model_ex, trt_status = 1, age_varying = F)
  
  res_unadjusted_multi_ex_g5 <- run_multi_cohort_sequential(start_ages = start_ages[[1]], unadjusted_tp_ctrl_ex,
                                                    unadjusted_tp_trt_ex, n_cycles = n_cycles, utilities = utilities,
                                                    cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  gc()
  
  return(list(unadjusted_limited = res_unadjusted_multi_g5, unadjusted_extended = res_unadjusted_multi_ex_g5))
}  
  # res_unadjusted_multi_ex_g10 <- run_multi_cohort_parallel(start_ages = start_ages[[2]], unadjusted_tp_ctrl_ex,
  #                                                         unadjusted_tp_trt_ex, n_cycles = n_cycles, utilities = utilities,
  #                                                         cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  
  
  # res_unadjusted_multi_ex_all <- run_multi_cohort_parallel(start_ages = start_ages[[3]], unadjusted_tp_ctrl_ex,
  #                                                         unadjusted_tp_trt_ex, n_cycles = n_cycles, utilities = utilities,
  #                                                         cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  # 
  # gc()
  
  
  
  
  
  # res_adjusted_multi_g10 <- run_multi_cohort_parallel(start_ages = start_ages[[2]], adjusted_tp_ctrl,
  #                                                    adjusted_tp_trt, n_cycles = n_cycles, utilities = utilities,
  #                                                    cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  # gc()
  
  # res_adjusted_multi_all <- run_multi_cohort_parallel(start_ages = start_ages[[3]], adjusted_tp_ctrl,
  #                                                    adjusted_tp_trt, n_cycles = n_cycles, utilities = utilities,
  #                                                    cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  # 
  # gc()
  
  
  

  # res_adjusted_multi_ex_g10 <- run_multi_cohort_parallel(start_ages = start_ages[[2]], adjusted_tp_ctrl_ex,
  #                                                    adjusted_tp_trt_ex, n_cycles = n_cycles, utilities = utilities,
  #                                                    cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  # gc()
  
  # res_adjusted_multi_ex_all <- run_multi_cohort_parallel(start_ages = start_ages[[3]], adjusted_tp_ctrl_ex,
  #                                                    adjusted_tp_trt_ex, n_cycles = n_cycles, utilities = utilities,
  #                                                    cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  # 
  # gc()
  
 
  
  
  
  # res_linear_multi_g10 <- run_multi_cohort_parallel(start_ages = start_ages[[2]], linear_tp_ctrl,
  #                                               linear_tp_trt, n_cycles = n_cycles, utilities = utilities,
  #                                               cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  
  
  # res_linear_multi_all <- run_multi_cohort_parallel(start_ages = start_ages[[3]], linear_tp_ctrl,
  #                                               linear_tp_trt, n_cycles = n_cycles, utilities = utilities,
  #                                               cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  # 
  # gc()
  
  
  
  
  # res_linear_multi_ex_g10 <- run_multi_cohort_parallel(start_ages = start_ages[[2]], linear_tp_ctrl_ex,
  #                                                  linear_tp_trt_ex, n_cycles = n_cycles, utilities = utilities,
  #                                                  cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  # gc()
  
  # res_linear_multi_ex_all <- run_multi_cohort_parallel(start_ages = start_ages[[3]], linear_tp_ctrl_ex,
  #                                                  linear_tp_trt_ex, n_cycles = n_cycles, utilities = utilities,
  #                                                  cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  # 
  # gc()
  
  
  
  
  
  
  
  
monotonic_spline_model <- function(stan_model, trial_data, trial_data_extra, filename, plot_outcome = F) {
  
  filename <- paste0(filename, "_monotonic_spline.RDs")
  
  control_data <- trial_data %>% filter(trt == 0)
  treatment_data <- trial_data %>% filter(trt == 1)
  control_data_extra <- trial_data_extra %>% filter(trt == 0)
  treatment_data_extra <- trial_data_extra %>% filter(trt == 1)
  
  smol_d3 <- cmdstan_fit_plot(mod = stan_model, increment = 0.1, target_age = age_grid, degree = 3,
                              ctrl_dat = control_data, ctrl_dat_ex = control_data_extra,
                              trt_dat = treatment_data, trt_dat_ex = treatment_data_extra,
                              target_probs_control = target_probs_control, target_probs_treatment = target_probs_treatment,
                              file_name = filename, plot_outcome = plot_outcome)
  
 
  
  
  spline_tp_ctrl <- posterior_spline_matrix(smol_d3$fit_ctrl, smol_d3$ctrl_spline_basis, newx = age_grid) 
  spline_tp_ctrl_ex <- posterior_spline_matrix(smol_d3$fit_ctrl_ex, smol_d3$ex_ctrl_spline_basis, newx = age_grid)
  spline_tp_trt <- posterior_spline_matrix(smol_d3$fit_trt, smol_d3$trt_spline_basis, newx = age_grid)
  spline_tp_trt_ex <- posterior_spline_matrix(smol_d3$fit_trt_ex, smol_d3$ex_trt_spline_basis, newx = age_grid)
  
  res_spline_multi_g5 <- run_multi_cohort_sequential(start_ages = start_ages[[1]], spline_tp_ctrl, spline_tp_trt,
                                                n_cycles = n_cycles, utilities = utilities,
                                                cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  
  # res_spline_multi_g10 <- run_multi_cohort_parallel(start_ages = start_ages[[2]], spline_tp_ctrl, spline_tp_trt,
  #                                               n_cycles = n_cycles, utilities = utilities,
  #                                               cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  
  
  # res_spline_multi_all <- run_multi_cohort_parallel(start_ages = start_ages[[3]], spline_tp_ctrl, spline_tp_trt,
  #                                               n_cycles = n_cycles, utilities = utilities,
  #                                               cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  # gc()
  
  res_spline_multi_ex_g5 <- run_multi_cohort_sequential(start_ages = start_ages[[1]], spline_tp_ctrl_ex, spline_tp_trt_ex,
                                                n_cycles = n_cycles, utilities = utilities,
                                                cost_ctrl = cost_ctrl, cost_trt = cost_trt)
 
  
  # res_spline_multi_ex_g10 <- run_multi_cohort_parallel(start_ages = start_ages[[2]], spline_tp_ctrl_ex, spline_tp_trt_ex,
  #                                                  n_cycles = n_cycles, utilities = utilities,
  #                                                  cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  gc()
  
  # res_spline_multi_ex_all <- run_multi_cohort_parallel(start_ages = start_ages[[3]], spline_tp_ctrl_ex, spline_tp_trt_ex,
  #                                                  n_cycles = n_cycles, utilities = utilities,
  #                                                  cost_ctrl = cost_ctrl, cost_trt = cost_trt)
  # 
  # gc()
  
  return(list(mono_spline_limited = res_spline_multi_g5, mono_spline_extended = res_spline_multi_ex_g5))
  
}
