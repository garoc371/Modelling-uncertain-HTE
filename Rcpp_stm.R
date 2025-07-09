# cppFunction('
#   List state_transition_both_rcpp(int n, NumericMatrix control_probs_matrix, NumericMatrix treatment_probs_matrix, NumericVector utility_values, NumericVector cost_values, int time_horizon, double treatment_cost = 0) {
#     NumericMatrix qalys(n, 2);  // 2 columns for control and treatment
#     NumericMatrix costs(n, 2);  // 2 columns for control and treatment
#     
#     // Initialize discount rate
#     double discount_rate = 0.03;
#     
#     for (int i = 0; i < n; ++i) {
#       std::string state_control = "healthy";
#       std::string state_treatment = "healthy";
#       double cum_qalys_control = 0;
#       double cum_qalys_treatment = 0;
#       double cum_costs_control = 0;
#       double cum_costs_treatment = 0;
#       double discount_factor = 1 / (1 + discount_rate);  // Initialize discount factor for the first year
#       
#       for (int t = 0; t < time_horizon; ++t) {
#         
#         double rand_num1 = R::runif(0, 1);  // Use the same random number for both arms
#         double rand_num2 = R::runif(0, 1);  // Use the same random number for both arms
# 
#         
#         if (state_control == "healthy") {
#           if (rand_num1 <= 0.5*control_probs_matrix(i, t)) {
#             state_control = "diseased";
#           }
#         } else if (state_control == "diseased") {
#           if (rand_num2 <= control_probs_matrix(i, t)) {
#             state_control = "dead";
#           }
#         }
#         
#         // Transition for treatment
#         if (state_treatment == "healthy") {
#           if (rand_num1 <= 0.5*treatment_probs_matrix(i, t)) {
#             state_treatment = "diseased";
#           }
#         } else if (state_treatment == "diseased") {
#           if (rand_num2 <= treatment_probs_matrix(i, t)) {
#             state_treatment = "dead";
#           }
#         }
#         
#         // Apply discount factor for the current year
#         
#         // accumulate costs and qalys after transition for control
#          if (state_control == "healthy") {
#            cum_qalys_control += utility_values[0]* discount_factor;
#            cum_costs_control += cost_values[0]* discount_factor;
#          } else if (state_control == "diseased") {
#            cum_qalys_control += utility_values[1]* discount_factor;
#            cum_costs_control += (cost_values[1])* discount_factor;
#          }
#         
#         // accumulate costs and qalys after transition for treatment
#          if (state_treatment == "healthy") {
#            cum_qalys_treatment += utility_values[0]* discount_factor;
#            cum_costs_treatment += cost_values[0]* discount_factor;
#          } else if (state_treatment == "diseased") {
#            cum_qalys_treatment += utility_values[1]* discount_factor;
#            cum_costs_treatment += (cost_values[1]+treatment_cost)* discount_factor;
#          }
#         
#         discount_factor = discount_factor * (1 / (1 + discount_rate));
#       }
#       
#       qalys(i, 0) = cum_qalys_control;
#       qalys(i, 1) = cum_qalys_treatment;
#       costs(i, 0) = cum_costs_control;
#       costs(i, 1) = cum_costs_treatment;
#     }
#     
#     return List::create(Named("qalys") = qalys, Named("costs") = costs);
#   }
# ')
cppFunction('
  List state_transition_both_rcpp(int n, NumericVector control_probs, NumericVector treatment_probs, NumericVector utility_values, NumericVector cost_values, NumericVector age, int time_horizon, double treatment_cost = 0) {
    NumericMatrix qalys(n, 2);  // 2 columns for control and treatment
    NumericMatrix costs(n, 2);  // 2 columns for control and treatment
    
    // Initialize discount rate
    double discount_rate = 0.03;
    
    for (int i = 0; i < n; ++i) {
      int age_idx = age[i] - 40;  // Assuming age 40 corresponds to the 0-based index
      std::string state_control = "healthy";
      std::string state_treatment = "healthy";
      double cum_qalys_control = 0;
      double cum_qalys_treatment = 0;
      double cum_costs_control = 0;
      double cum_costs_treatment = 0;
      double discount_factor = 1 / (1 + discount_rate);  // Initialize discount factor for the first year
      
      for (int t = 0; t < time_horizon; ++t) {
        
        double rand_num1 = R::runif(0, 1);  // Use the same random number for both arms
        double rand_num2 = R::runif(0, 1);  // Use the same random number for both arms

        int current_age_idx = age_idx + t;
        
        // Transition for control
        if (state_control == "healthy") {
          if (rand_num1 <= 0.5 * control_probs[current_age_idx]) {
            state_control = "diseased";
          }
        } else if (state_control == "diseased") {
          if (rand_num2 <= control_probs[current_age_idx]) {
            state_control = "dead";
          }
        }
        
        // Transition for treatment
        if (state_treatment == "healthy") {
          if (rand_num1 <= 0.5 * treatment_probs[current_age_idx]) {
            state_treatment = "diseased";
          }
        } else if (state_treatment == "diseased") {
          if (rand_num2 <= treatment_probs[current_age_idx]) {
            state_treatment = "dead";
          }
        }
        
         // Apply discount factor for the current year
         
         // accumulate costs and qalys after transition for control
          if (state_control == "healthy") {
            cum_qalys_control += utility_values[0]* discount_factor;
            cum_costs_control += cost_values[0]* discount_factor;
          } else if (state_control == "diseased") {
            cum_qalys_control += utility_values[1]* discount_factor;
            cum_costs_control += (cost_values[1])* discount_factor;
          }
         
         // accumulate costs and qalys after transition for treatment
          if (state_treatment == "healthy") {
            cum_qalys_treatment += utility_values[0]* discount_factor;
            cum_costs_treatment += (cost_values[0]+treatment_cost)* discount_factor;
          } else if (state_treatment == "diseased") {
            cum_qalys_treatment += utility_values[1]* discount_factor;
            cum_costs_treatment += (cost_values[1])* discount_factor;
          }
         
         discount_factor = discount_factor * (1 / (1 + discount_rate));
       }
       
       qalys(i, 0) = cum_qalys_control;
       qalys(i, 1) = cum_qalys_treatment;
       costs(i, 0) = cum_costs_control;
       costs(i, 1) = cum_costs_treatment;
     }
    
    return List::create(Named("qalys") = qalys, Named("costs") = costs);
  }
')


# sim_true <- function(control_type, treatment_type) {
#   # sim_age_target <- sample(40:60, size = n_target, replace = TRUE)
#   
#   sim_age_target <- rep(40, n_target) %>% as.integer()
#   
#   utility_values <- c(healthy = utility_healthy, diseased = utility_diseased, dead = utility_dead)
#   cost_values <- c(healthy = cost_healthy, diseased = cost_diseased, dead = cost_dead)
#   
#   control_grid <- treatment_grid <- expand.grid(age = sim_age_target, cycle = 1:n_cycles) %>%
#     arrange(cycle) %>%
#     mutate(age = age + cycle - 1)
#   
#   ctrl_probs_logit <- control_outcome_logit(baseline_probs, control_grid$age, control_type)
#   trt_effect_logit <- treatment_effect_logit(treatment_grid$age, treatment_type)
#   
#   control_grid$outcome_prob <- plogis(control_outcome_logit(baseline_probs, control_grid$age, control_type))
#   treatment_grid$outcome_prob <- generate_probabilities(treatment_grid$age, baseline_probs, control_type = control_type, treatment_type = treatment_type)
#   
#   control_matrix <- matrix(control_grid$outcome_prob, nrow = length(sim_age_target), ncol = n_cycles, byrow = FALSE)
#   treatment_matrix <- matrix(treatment_grid$outcome_prob, nrow = length(sim_age_target), ncol = n_cycles, byrow = FALSE)
#   
#   # Perform microsimulations
#   microsim_true <- replicate(n = 1000, state_transition_both_rcpp(n_target, control_matrix, treatment_matrix,
#                                                                   utility_values, cost_values, time_horizon, 
#                                                                   treatment_cost_per_cycle), simplify = F)
#   
#   microsim_true_simplified <- lapply(microsim_true, function(sim_results){
#     delta_qaly = sim_results$qalys[,2] - sim_results$qalys[,1]
#     delta_cost <- sim_results$costs[,2] - sim_results$costs[,1]
#     
#     list(delta_qalys = sum(delta_qaly), delta_costs = sum(delta_cost))
#   })
#   
#   microsim_res <- bind_rows(microsim_true_simplified)
#   
#   
#   
#   true_cohort_ctrl <- plogis(control_outcome_logit(baseline_probs = baseline_probs, age = age_grid, type = control_type))
#   true_cohort_trt <- generate_probabilities(age_grid, baseline_probs, control_type = control_type, treatment_type = treatment_type)
#   
#   utilities <- c(utility_healthy, utility_diseased, 0)
#   cost_ctrl <- c(cost_healthy, cost_diseased,0)
#   cost_trt <- c(cost_healthy , cost_diseased+treatment_cost_per_cycle, 0)
#   
#   res_true_cohort <- cohort_state_transition(true_cohort_ctrl, true_cohort_trt, utilities = utilities, cost_ctrl = cost_ctrl, cost_trt = cost_trt,
#                                               time_horizon = n_cycles)
#   
#   
#   return(list(qalys = microsim_res$delta_qalys, costs = microsim_res$delta_costs, res_cohort = res_true_cohort))
# }

sim_true_microsim <- function(control_type, treatment_type) {
  
  age_grid <- 40:80
  limited_probs <- uniform_probs(age_grid)
  extra_probs <- half_normal_probs(age_grid, sd = 5)
  combined_probs <- c(0.6*limited_probs[1:21], 0.4*extra_probs[22:41])
  
  sim_age_target <- sample(40:80, size = n_target, prob = combined_probs, replace = T)
  
  target_probs <- target_prob_vector(control_type, treatment_type, max_age = 120)
  control_probs <- target_probs$control
  treatment_probs <- target_probs$treatment
  
  utility_values <- c(healthy = utility_healthy, diseased = utility_diseased, dead = utility_dead)
  cost_values <- c(healthy = cost_healthy, diseased = cost_diseased, dead = cost_dead)
  
  # Perform microsimulations
  microsim_true <- replicate(n = 1000, state_transition_both_rcpp(n_target, control_probs, treatment_probs,
                                                                  utility_values, cost_values, age = sim_age_target, 
                                                                  time_horizon,treatment_cost_per_cycle), simplify = F)
  
  microsim_true_simplified <- lapply(microsim_true, function(sim_results){
    delta_qaly = sim_results$qalys[,2] - sim_results$qalys[,1]
    delta_cost <- sim_results$costs[,2] - sim_results$costs[,1]
    
    list(delta_qalys = sum(delta_qaly), delta_costs = sum(delta_cost))
  })
  
  microsim_res <- bind_rows(microsim_true_simplified)
  return(list(qalys = microsim_res$delta_qalys, costs = microsim_res$delta_costs, 
              target_prop = combined_probs))
}

sim_true_cohort <- function(control_type, treatment_type) {
  
  age_grid <- 40:80
  limited_probs <- uniform_probs(age_grid)
  extra_probs <- half_normal_probs(age_grid, sd = 5)
  combined_probs <- c(0.6*limited_probs[1:21], 0.4*extra_probs[22:41])
  
  true_cohort <- target_prob_vector(control_type, treatment_type, max_age = 120)
  true_cohort_ctrl <- matrix(true_cohort$control, nrow = 1)
  true_cohort_trt <- matrix(true_cohort$treatment, nrow =1)
  
  utilities <- c(utility_healthy, utility_diseased, 0)
  cost_ctrl <- c(cost_healthy, cost_diseased, 0)
  cost_trt <- c(cost_healthy+treatment_cost_per_cycle , cost_diseased , 0)
  
  ncores <- 32
  cl <- makeCluster(ncores, type = "SOCK")
  clusterExport(cl, c("cohort_state_transition", "n_cycles", "n_target", "compute_delta"))
  start_ages <- 40:80
  extended_age_grid <- 40:120
  res_true_cohort_list <- parLapply(cl,start_ages, function(start_age) {
    
    age_idx_range = which(extended_age_grid >= start_age & extended_age_grid < (start_age + 40))
    cohort_state_transition(true_cohort_ctrl, true_cohort_trt, utilities, cost_ctrl, cost_trt, n_cycles, age_idx_range)
  })
  
  stopCluster(cl)
  
  cohort_true_collapsed <- bind_rows(res_true_cohort_list)
  res_true_cohort <- cohort_true_collapsed %>% 
    summarise(delta_QALYs = weighted.mean(cohort_true_collapsed$delta_QALYs, combined_probs),
              delta_costs = weighted.mean(cohort_true_collapsed$delta_costs, combined_probs))
  
  return(list(res_true_cohort = res_true_cohort, target_prop = combined_probs))
}



