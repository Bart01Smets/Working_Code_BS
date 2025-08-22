#
# Workbench script for the EPI-ECON modelling work
#
# - Deterministic and stochastic modelling
# - Still includes option for both binomial transitions and stochastic beta
######################################################### #
getwd()


setwd("C:/Users/Bart Smets/OneDrive/Documenten/GitHub/Working_Code_BS")

# clear workbench
rm(list=ls())

# load functions
source("scenarios/epi_econ_lib_policy.R")

# SETUP   ####
####################

# Parameters
parameters <- list(
  gamma = 1/7,             # Recovery rate (γ) 1/7, alternatively 1/4
  beta = 3/10+1/7,       # Transmission rate (β) 3/10
  rho = 0.72 / 365,        # Discounting rate (ρ)
  pi = 0.0062,             # Infection fatality rate (π)
  kappa = 197,             # Expected cost of infection (κ)
  ni0 = 0.0005,         # Initial infected population (Ni0)
  ns0 = 0.9995,     # Initial susceptible population (Ns0)
  nr0 = 0,                 # Initial recovered population (Nr0)
  nd0 = 0, # Initial dead population (Nd0)
  v = 31755,
  fx = 123,                # Exchange rate multiplier for USD conversion
  time_horizon = 1500,      # Time of shock
  rng_seed = 150,
  R0= (3/10+1/7)/(1/7),
  pop_size = 1e6,
  infect_thres = 1,
  utility_type = "Log",
  sigma = 0.19,
  bool_regular_sird = FALSE,  # NEW FLAG
  bool_daily_cost_minimizing = FALSE,
  # --- Policy knobs (default-off) ---
  # 1) Activity cap / lockdown triggers
  a_cap = 0.6,                     # <=1; e.g. 0.6 for partial lockdown
  lock_trigger_Ni = 5,            # integer infections to activate cap; NA disables
  lock_release_Ni = 5,            # integer infections to release; NA disables
  
  # 2) Pigouvian tax on activity (can be prevalence-linked)
  tau0 = 0,                      # baseline tax in “utility units” per unit a
  tau_per_Ni = 0,                # extra tax × prevalence (per capita)
  
  # 3) Masks / distancing mandate
  mask_effect = 0,               # 0.35 = 35% per-contact reduction, example: 0.4
  mask_compliance = 0,           # share complying, example: 0.7
  
  # 4) Testing & isolation (minimal: increases effective removal)
  test_rate = 0,                 # fraction of I tested per day
  isolation_effect = 0,          # 0..1: how much testing adds to removal
  
  # 5) Hospital capacity (smooth excess mortality)
  healthcare_capacity = NULL,      # e.g. 800; NULL = disabled
  excess_mortality_multiplier = 1.0, # >=1
  capacity_slope = 0.0,            # >0 for smooth ramp; 0 for step
  
  # 6) Vaccination rollout
  nu = 0.0,                        # daily S->R rate;0.01
  start_vax_day = Inf,             # day to start vaccinating 90
  vax_efficacy_inf = 1,          # if you later add a V state; here kept simple;1
  
  # 7) Border controls / importations
  lambda_import = 0.0,             # mean imported cases per day (Poisson)
  
  # 8) Sectoral NPIs (optional β decomposition; leave zeros to ignore)
  beta_home   = 0,
  beta_work   = 0,
  beta_school = 0,
  beta_other  = 0,
  school_open = TRUE,
  
  
  # 9) Information campaign (behavioral risk multiplier used in a_t choice)
  risk_multiplier = 1,#example value: 1.5
  campaign_start = Inf,           # inclusive example value: 40
  campaign_end   = Inf          # inclusive example value:1000
  
)

# parameters$beta_home   <- parameters$beta * 0.30 * 1.05
# parameters$beta_work   <- parameters$beta * 0.35 * 0.60
# parameters$beta_school <- parameters$beta * 0.20 * 1.00
# parameters$beta_other  <- parameters$beta * 0.15 * 0.60
# parameters$school_open <- TRUE
# 
# # Assumptions:
# hosp_rate <- 0.02          # 2% of infectives need ICU-level care
# icu_beds_per_100k <- 20    # ICU beds effectively available for this disease
# 
# parameters$healthcare_capacity <- round((parameters$pop_size / 1e5) * icu_beds_per_100k / hosp_rate)
# parameters$excess_mortality_multiplier <- 2.0   # IFR can double when ICU saturates
# parameters$capacity_slope <- 6.0               # steeper ramp around the threshold


parameters$pop_size<-1e4
# define number of stochastic runs
num_experiments <- 100

# define fadeout threshold
fadeout_threshold = 100

# define default parameters (for development)
plot_tag <- 'dev'
update_function <- get_transitions_stochastic

# RUN DETERMINISTIC - ODE solver   ####
########################################
# Initial state variables including separate costs
initial_state <- c(Ns = parameters$ns0, 
                   Ni = parameters$ni0, 
                   Nr = parameters$nr0, 
                   Nd = parameters$nd0, 
                   HealthCost = 0, 
                   SocialActivityCost = 0, 
                   TotalCost = 0,
                   a_t = NA, 
                   u_t = NA, 
                   Rt= NA) # add a_t and u_t to keep track of this over time

# Time sequence for pre-shock
times <- seq(0, parameters$time_horizon, by = 1)

# RUN STOCHASTIC BETA REALISATIONS  ####
########################################
# note: make sure the sampled beta is also used for a_t, u_t, Lambda_s and Lambda_i

# get reference: deterministic model
output_sim_deterministic <- run_sir_binomial(initial_state = initial_state, 
                                             times = times, 
                                             parameters = parameters,
                                             bool_stochastic_beta = FALSE,
                                             update_function = get_transitions_deterministic)
# Print deterministic peak infection time
det_peak_time <- which.max(output_sim_deterministic$Ni) - 1
cat(sprintf("Deterministic Peak Infection Time: %d\n", det_peak_time))


output_experiments <- run_experiments(initial_state = initial_state, 
                                      times = times, 
                                      parameters = parameters, 
                                      bool_stochastic_beta = TRUE, 
                                      update_function = get_transitions_deterministic, 
                                      num_experiments= num_experiments)

# inspect results
compare_sim_output(output_experiments, output_sim_deterministic, plot_tag='stoch beta')




# RUN STOCHASTIC BINIOMIAL MODEL REALISATIONS   ####
####################################################

# get reference: deterministic model
output_sim_deterministic <- run_sir_binomial(initial_state = initial_state, 
                                             times = times, 
                                             parameters = parameters,
                                             update_function = get_transitions_deterministic)
# Print deterministic peak infection time
det_peak_time <- which.max(output_sim_deterministic$Ni) - 1
cat(sprintf("Deterministic Peak Infection Time: %d\n", det_peak_time))

head(output_sim_deterministic)
output_experiments <- run_experiments(initial_state = initial_state, 
                                      times = times, 
                                      parameters = parameters,
                                      bool_stochastic_beta = FALSE,
                                      update_function = get_transitions_stochastic, 
                                      num_experiments)

# inspect all results
compare_sim_output(output_experiments, output_sim_deterministic, plot_tag='binomial')

# inspect results excl fadeout
compare_sim_output(output_experiments, output_sim_deterministic, plot_tag='binomial',
                   fadeout_threshold = fadeout_threshold)


# # ==============================
# # MULTI-PARAMETER SENSITIVITY ANALYSIS
# # ==============================
# 
# library(ggplot2)
# 
# # Define parameters and their sensitivity ranges
# param_ranges <- list(
#   v = seq(10000, 50000, length.out = 10),
#   beta = seq(0.2, 0.6, length.out = 10),
#   pi = seq(0.001, 0.01, length.out = 10),
#   gamma = seq(1/14, 1/3, length.out = 10),
#   ni0 = seq(0.0005, 0.01, length.out = 10),
#   ns0 = seq(0.85, 0.999, length.out = 10),
#   pop_size = round(seq(10000, 100000, length.out = 10)),
#   time_horizon = seq(250, 2500, length.out = 10),
#   rho = seq(0, 100/365, length.out = 10)
# )
# 
# # Open a multipage PDF
# pdf("sensitivity_all.pdf", width = 8, height = 6)
# 
# # Loop over each parameter
# for (sensitivity_target in names(param_ranges)) {
#   cat("Running sensitivity for:", sensitivity_target, "\n")
#   
#   param_grid <- data.frame(value = param_ranges[[sensitivity_target]])
#   colnames(param_grid) <- sensitivity_target
#   baseline_value <- parameters[[sensitivity_target]]
#   sensitivity_results <- data.frame()
#   
#   for (i in 1:nrow(param_grid)) {
#     cat("  Scenario", i, "of", nrow(param_grid), "\n")
#     
#     # Update parameter
#     parameters[[sensitivity_target]] <- param_grid[[sensitivity_target]][i]
#     
#     # Recalculate R0 if needed
#     if (sensitivity_target == "beta" || sensitivity_target == "gamma") {
#       parameters$R0 <- parameters$beta / parameters$gamma
#     }
#     
#     # Initial state
#     initial_state <- c(Ns = parameters$ns0,
#                        Ni = parameters$ni0,
#                        Nr = parameters$nr0,
#                        Nd = parameters$nd0,
#                        HealthCost = 0,
#                        SocialActivityCost = 0,
#                        TotalCost = 0,
#                        a_t = NA, u_t = NA, Rt = NA)
#     
#     times <- seq(0, parameters$time_horizon, by = 1)
#     
#     # Deterministic run
#     output_det <- run_sir_binomial(
#       initial_state = initial_state,
#       times = times,
#       parameters = parameters,
#       update_function = get_transitions_deterministic
#     )
#     total_cost_det <- output_det[nrow(output_det), "TotalCost"]
#     
#     # Stochastic run
#     output_stoch <- run_experiments(
#       initial_state = initial_state,
#       times = times,
#       parameters = parameters,
#       bool_stochastic_beta = FALSE,
#       update_function = get_transitions_stochastic,
#       num_experiments = num_experiments
#     )
#     total_costs_stoch <- sapply(1:num_experiments, function(j) {
#       output_stoch$output_all[j, length(times), "TotalCost"]
#     })
#     
#     # Statistics
#     mean_cost_stoch <- mean(total_costs_stoch)
#     q_lo <- quantile(total_costs_stoch, 0.025)
#     q_hi <- quantile(total_costs_stoch, 0.975)
#     se <- sd(total_costs_stoch) / sqrt(num_experiments)
#     ci_lower <- mean_cost_stoch - 1.96 * se
#     ci_upper <- mean_cost_stoch + 1.96 * se
#     val <- param_grid[[sensitivity_target]][i]
#     
#     # Store results
#     sensitivity_results <- rbind(
#       sensitivity_results,
#       data.frame(Parameter = sensitivity_target, Value = val, Type = "Deterministic",
#                  TotalCost = total_cost_det, Q025 = NA, Q975 = NA, CI025 = NA, CI975 = NA),
#       data.frame(Parameter = sensitivity_target, Value = val, Type = "Stochastic",
#                  TotalCost = mean_cost_stoch, Q025 = q_lo, Q975 = q_hi,
#                  CI025 = ci_lower, CI975 = ci_upper)
#     )
#   }
#   
#   # Save individual CSV
#   write.csv(sensitivity_results, paste0("sensitivity_", sensitivity_target, ".csv"), row.names = FALSE)
#   
#   # Plot for current parameter
#   p <- ggplot(sensitivity_results, aes(x = Value, y = TotalCost, color = Type)) +
#     geom_ribbon(data = subset(sensitivity_results, Type == "Stochastic"),
#                 aes(ymin = Q025, ymax = Q975, fill = "95% Credible Interval"),
#                 alpha = 0.2, inherit.aes = TRUE) +
#     geom_ribbon(data = subset(sensitivity_results, Type == "Stochastic"),
#                 aes(ymin = CI025, ymax = CI975, fill = "95% Confidence Interval"),
#                 alpha = 0.3, inherit.aes = TRUE) +
#     geom_line(data = subset(sensitivity_results, Type == "Stochastic"),
#               aes(y = CI025), color = "red", linetype = "solid", size = 0.7) +
#     geom_line(data = subset(sensitivity_results, Type == "Stochastic"),
#               aes(y = CI975), color = "red", linetype = "solid", size = 0.7) +
#     geom_vline(xintercept = baseline_value, linetype = "dashed", color = "blue", size = 1, alpha = 0.8) +
#     geom_line() +
#     geom_point(size = 2) +
#     labs(title = paste("Sensitivity of Total Cost to", sensitivity_target),
#          x = paste("Value of", sensitivity_target),
#          y = "Total Cost (per capita)") +
#     theme_minimal() +
#     scale_color_manual(name = "Transitions",
#                        values = c("Deterministic" = "black", "Stochastic" = "steelblue")) +
#     scale_fill_manual(name = "Uncertainty",
#                       values = c(
#                         "95% Credible Interval" = "steelblue",
#                         "95% Confidence Interval" = "red"
#                       ))
#   
#   print(p)
# }
# 
# # Close the multipage PDF
# dev.off()
# 
# # ==============================
# # END MULTI-PARAMETER SENSITIVITY
# # ==============================