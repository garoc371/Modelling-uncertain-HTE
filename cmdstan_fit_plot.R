# wrapper function for running bayesian spline models by arm,while
# producing plots of outcome surface at the same time

cmdstan_fit_plot <- function(mod, increment, target_age, degree,
                             ctrl_dat, ctrl_dat_ex, trt_dat, trt_dat_ex,
                             target_probs_control, target_probs_treatment, ispline = TRUE, plot_outcome = F, file_name){
  
  knots1 <- knots_quantile(ctrl_dat$age, increment)
  knots2 <- knots_quantile(trt_dat$age,increment)
  knots3 <- knots_quantile(ctrl_dat_ex$age, increment)
  knots4 <- knots_quantile(trt_dat_ex$age, increment)
  
  mspline_degree <- degree
  
  
  pred_age_target <- target_age
  bknots <- range(pred_age_target)
  bknots[2] <- bknots[2] + 0.5
  
  if (ispline){
    ctrl_spline_basis <- iSpline(ctrl_dat$age, knots = knots1, Boundary.knots = bknots, degree = mspline_degree, intercept = T)
    trt_spline_basis <- iSpline(trt_dat$age, knots = knots2, Boundary.knots = bknots, degree = mspline_degree, intercept = T)
    
    ex_ctrl_spline_basis <- iSpline(ctrl_dat_ex$age, knots = knots3, Boundary.knots = bknots, degree = mspline_degree, intercept = T)
    ex_trt_spline_basis <- iSpline(trt_dat_ex$age, knots = knots4, Boundary.knots = bknots, degree = mspline_degree, intercept = T)
  }else{
    ctrl_spline_basis <- bSpline(ctrl_dat$age, knots = knots1[-1], Boundary.knots = bknots, degree = mspline_degree, intercept = T)
    trt_spline_basis <- bSpline(trt_dat$age, knots = knots2[-1], Boundary.knots = bknots, degree = mspline_degree, intercept = T)
    
    ex_ctrl_spline_basis <- bSpline(ctrl_dat_ex$age, knots = knots3[-1], Boundary.knots = bknots, degree = mspline_degree, intercept = T)
    ex_trt_spline_basis <- bSpline(trt_dat_ex$age, knots = knots4[-1], Boundary.knots = bknots, degree = mspline_degree, intercept = T)
  }
  nbasis <- dim(ctrl_spline_basis)[2]
  # mod <- cmdstan_model("rw_spline.stan")
  ctrl_data_list <- list(N = nrow(ctrl_dat),
                         Y = ctrl_dat$outcome,
                         m = nbasis,
                         S = ctrl_spline_basis)
  
  ex_ctrl_data_list <- list(N = nrow(ctrl_dat_ex),
                            Y = ctrl_dat_ex$outcome,
                            m = nbasis,
                            S = ex_ctrl_spline_basis)
  
  trt_data_list <- list(N = nrow(trt_dat),
                        Y = trt_dat$outcome,
                        m = nbasis,
                        S = trt_spline_basis)
  
  ex_trt_data_list <- list(N = nrow(trt_dat_ex),
                           Y = trt_dat_ex$outcome,
                           m = nbasis,
                           S = ex_trt_spline_basis)
  
  
  fit_ctrl <- mod$sample(
    data = ctrl_data_list,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.95
  )
  
  fit_ctrl_ex <- mod$sample(
    data = ex_ctrl_data_list,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.95
  )
  
  fit_trt <- mod$sample(
    data = trt_data_list,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.95
  )
  
  fit_trt_ex <- mod$sample(
    data = ex_trt_data_list,
    seed = 1234,
    chains = 4,
    parallel_chains = 4,
    refresh = 500,
    adapt_delta = 0.95
  )
  
  ctrl_pred_basis <- predict(ctrl_spline_basis, newx = pred_age_target)
  ctrl_pred_basis_ex <- predict(ex_ctrl_spline_basis, newx = pred_age_target)
  trt_pred_basis <- predict(trt_spline_basis, newx = pred_age_target)
  trt_pred_basis_ex <- predict(ex_trt_spline_basis, newx = pred_age_target)
  

  if (plot_outcome){
    pred_y0s <- cmdstan_epred_sep(fit_ctrl, ctrl_spline_basis, pred_age_target)
    pred_y1s <- cmdstan_epred_sep(fit_trt, trt_spline_basis, pred_age_target)
    
    pred_diff <- pred_y1s - pred_y0s

    posterior_control_sep <- posterior_summary(pred_y0s, "Control", pred_age_target)
    posterior_trt_sep <- posterior_summary(pred_y1s, "Treatment", pred_age_target)
    posterior_effect <- posterior_summary(pred_diff, "Treatment_effect", pred_age_target)
  
    
    rm(pred_y0s, pred_y1s)
    gc()
    
    pred_y0s_ex <- cmdstan_epred_sep(fit_ctrl_ex, ex_ctrl_spline_basis, pred_age_target)
    pred_y1s_ex <- cmdstan_epred_sep(fit_trt_ex, ex_trt_spline_basis, pred_age_target)
    pred_diff_ex <- pred_y1s_ex - pred_y0s_ex

    posterior_control_exsep <- posterior_summary(pred_y0s_ex, "Control", pred_age_target)
    posterior_trt_exsep <- posterior_summary(pred_y1s_ex, "Treatment", pred_age_target)
    posterior_effect_ex <- posterior_summary(pred_diff_ex, "Treatment_effect", pred_age_target)
    
   
    rm(pred_y0s_ex, pred_y1s_ex)
    
    true_outcome_control <- outcome_data(target_probs_control, "Control", age_target) %>% mutate(sample = "True")
    true_outcome_trt <- outcome_data(target_probs_treatment, "Treatment", age_target) %>% mutate(sample = "True")
    
    true_effect <- target_probs_treatment - target_probs_control
    true_effect <- outcome_data(true_effect, "Treatment_effect", age_target) %>% mutate(sample = "True")
    # true_effect_logit <- logit(target_probs_treatment)-logit(target_probs_control)
    # true_effect_logit <- outcome_data(true_effect_logit, "Treatment_effect(logit)", age_target) %>% mutate(sample = "True")
    
    sep_spline_combined_data <- bind_rows(
      posterior_control_sep,
      posterior_trt_sep
    ) 
    
    true_data <- bind_rows(true_outcome_control,
                           true_outcome_trt)
    
    s1 <- plot_outcome_surface(sep_spline_combined_data, true_data = true_data)
    rm(sep_spline_combined_data)
    gc()
    
    ex_sep_ispline_combined_data <- bind_rows(
      posterior_control_exsep,
      posterior_trt_exsep
    )
    s2 <- plot_outcome_surface(ex_sep_ispline_combined_data, true_data = true_data)
    rm(ex_sep_ispline_combined_data, true_data)
    gc()
    
    s3 <- plot_outcome_surface(posterior_effect, true_data = true_effect)
    s4 <- plot_outcome_surface(posterior_effect_ex, true_data = true_effect)
    
    
    plots_all <- list(limited_outcome = s1, extended_outcome = s2,
                      limited_effect = s3, extended_effect = s4)
    
    
   
    saveRDS(plots_all, file_name)
    
    
    
  }
  
  
  
  return(list(fit_ctrl = fit_ctrl, fit_trt = fit_trt, fit_ctrl_ex = fit_ctrl_ex,
              fit_trt_ex = fit_trt_ex, 
              ctrl_spline_basis = ctrl_spline_basis, trt_spline_basis = trt_spline_basis,
              ex_ctrl_spline_basis = ex_ctrl_spline_basis, ex_trt_spline_basis = ex_trt_spline_basis))
  
}  
