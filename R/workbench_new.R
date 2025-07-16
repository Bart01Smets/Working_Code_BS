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
source("R/epi_econ_lib_new.R")

# SETUP   ####
####################

# Parameters
parameters <- list(
  gamma = 1/7,             # Recovery rate (γ) 1/7, alternatively 1/3
  beta = 3/10+1/7,       # Transmission rate (β) 3/10
  delta = 0.67 / 365,      # Arrival rate of cure (δ)
  rho = 0.05 / 365,        # Discounting rate (ρ)
  pi = 0.0062,             # Infection fatality rate (π)
  sigma = 0.10,            # Determines stochastic beta_t =  N(beta,sigma)
  kappa = 197,             # Expected cost of infection (κ)
  ni0 = 0.000527,         # Initial infected population (Ni0)297
  ns0 = 0.999447004,     # Initial susceptible population (Ns0)527
  nr0 = 0.000024840,                 # Initial recovered population (Nr0)
  nd0 = 0.000000156, # Initial dead population (Nd0)
  v = 31755,
  alpha = 1,               # Altruism parameter (default set to no altruism)
  fx = 123,                # Exchange rate multiplier for USD conversion
  time_horizon = 1500,      # Time of shock
  utility_type = "Log",    # Utility type: "Log" or "Quadratic"
  #rng_seed = 123,
  R0= (3/10+1/7)/(1/7),
  pop_size = 1e6,
  infect_thres = 1
)

# define population size
parameters$pop_size <- 1e4

# define number of stochastic runs
num_experiments <- 100

# define fadeout threshold
fadeout_threshold = 100

# define default parameters (for development)
plot_tag <- 'dev'
bool_stochastic_beta <- FALSE
update_function <- get_transitions_stochastic

# RUN DETERMINISTIC - ODE solver   ####
########################################

# Setting tolerance values
parameters$nifinal <- 10^-9
parameters$tolerance <- 10^-6

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
                                      num_experiments)

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
# # SENSITIVITY ANALYSIS STARTS HERE
# # ==============================
# 
# # Choose the parameter for sensitivity
# sensitivity_target <- "pop_size"  # options: "v", "beta", "pi", "gamma", "ni0", "ns0", "pop_size", "time_horizon"
# 
# # Define ranges for each parameter
# param_ranges <- list(
#   v = seq(10000, 50000, length.out = 10),
#   beta = seq(0.2, 0.6, length.out = 10),
#   pi = seq(0.001, 0.01, length.out = 10),
#   gamma = seq(1/14, 1/3, length.out = 10),
#   ni0 = seq(0.0005, 0.01, length.out = 10),
#   ns0 = seq(0.85, 0.999, length.out = 10),
#   pop_size = round(seq(10000, 100000, length.out = 10)),
#   time_horizon = seq(250, 2500, length.out = 10)
# )
# 
# # Extract grid for selected parameter
# param_grid <- data.frame(value = param_ranges[[sensitivity_target]])
# colnames(param_grid) <- sensitivity_target
# 
# sensitivity_results <- data.frame()
# 
# for (i in 1:nrow(param_grid)) {
#   cat("Running scenario", i, "of", nrow(param_grid), "for", sensitivity_target, "\n")
# 
#   # Update selected parameter
#   parameters[[sensitivity_target]] <- param_grid[[sensitivity_target]][i]
# 
#   # Recompute derived quantities if needed
#   if (sensitivity_target == "beta" || sensitivity_target == "gamma") {
#     parameters$R0 <- parameters$beta / parameters$gamma
#   }
# 
#   # Initial state
#   initial_state <- c(Ns = parameters$ns0,
#                      Ni = parameters$ni0,
#                      Nr = parameters$nr0,
#                      Nd = parameters$nd0,
#                      HealthCost = 0,
#                      SocialActivityCost = 0,
#                      TotalCost = 0,
#                      a_t = NA, u_t = NA, Rt = NA)
# 
#   times <- seq(0, parameters$time_horizon, by = 1)
# 
#   # Deterministic run
#   output_det <- run_sir_binomial(
#     initial_state = initial_state,
#     times = times,
#     parameters = parameters,
#     update_function = get_transitions_deterministic
#   )
#   total_cost_det <- output_det[nrow(output_det), "TotalCost"]
# 
#   # Stochastic run
#   output_stoch <- run_experiments(
#     initial_state = initial_state,
#     times = times,
#     parameters = parameters,
#     bool_stochastic_beta = FALSE,
#     update_function = get_transitions_stochastic,
#     num_experiments = num_experiments
#   )
#   total_costs_stoch <- sapply(1:num_experiments, function(j) {
#     output_stoch$output_all[j, length(times), "TotalCost"]
#   })
# # Compute 95% quantile interval for stochastic costs
# q_lo <- quantile(total_costs_stoch, 0.025)
# q_hi <- quantile(total_costs_stoch, 0.975)
# 
# 
#   mean_cost_stoch <- mean(total_costs_stoch)
# 
#   # Store results
#   val <- param_grid[[sensitivity_target]][i]
# 
#   sensitivity_results <- rbind(
#     sensitivity_results,
#     data.frame(Parameter = sensitivity_target, Value = val, Type = "Deterministic",
#                TotalCost = total_cost_det, Q025 = as.numeric(NA), Q975 = as.numeric(NA)),
#     data.frame(Parameter = sensitivity_target, Value = val, Type = "Stochastic",
#                TotalCost = mean_cost_stoch, Q025 = q_lo, Q975 = q_hi)
#   )
# }
# 
# 
# 
# # # Save results
#  write.csv(sensitivity_results, paste0("sensitivity_", sensitivity_target, ".csv"), row.names = FALSE)
# 
# # # Plot
# # library(ggplot2)
# # ggplot(sensitivity_results, aes(x = Value, y = TotalCost, color = Type)) +
# #   geom_line() +
# #   geom_point(size = 2) +
# #   labs(title = paste("Sensitivity of Total Cost to", sensitivity_target),
# #        x = paste("Value of", sensitivity_target),
# #        y = "Total Cost (per capita)") +
# #   theme_minimal() +
# #   scale_color_manual(values = c("Deterministic" = "black", "Stochastic" = "steelblue"))
# library(ggplot2)
# 
# ggplot(sensitivity_results, aes(x = Value, y = TotalCost, color = Type)) +
#   # 95% quantile ribbon for stochastic only
#   geom_ribbon(data = subset(sensitivity_results, Type == "Stochastic"),
#               aes(x = Value, ymin = Q025, ymax = Q975, fill ="95% Quantile Band"),
#               alpha = 0.2, inherit.aes = FALSE) +
# 
#   # Lines and points
#   geom_line() +
#   geom_point(size = 2) +
# 
#   labs(title = paste("Sensitivity of Total Cost to", sensitivity_target),
#        x = paste("Value of", sensitivity_target),
#        y = "Total Cost (per capita)") +
# 
#   theme_minimal() +
#   scale_color_manual(name="Transitions", values = c("Deterministic" = "black", "Stochastic" = "steelblue")) +
#   scale_fill_manual(name="Uncertainty", values = c("95% Quantile Band" = "steelblue"))
# 
# # ==============================
# # SENSITIVITY ANALYSIS ENDS HERE
# # ==============================


