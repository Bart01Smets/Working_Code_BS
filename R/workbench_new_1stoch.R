#
# Workbench script for the EPI-ECON modelling work
# Single binomial-draw stochastic run vs deterministic baseline
######################################################### #

getwd()
setwd("C:/Users/Bart Smets/OneDrive/Documenten/GitHub/Working_Code_BS")

# clear workbench
rm(list = ls())

# load functions
source("R/epi_econ_lib_new_1stoch.R")

# ========== SETUP ==========
parameters <- list(
  gamma = 1/7,                      # Recovery rate (γ)
  beta = 3/10 + 1/7,               # Transmission rate (β)
  rho = 0.72 / 365,                 # Discount rate (ρ)
  pi = 0.0062,                      # Infection fatality rate (π)
  kappa = 197,                      # Infection cost (κ) - used in a_function
  ni0 = 0.0005,                     # Initial infected share
  ns0 = 0.9995,                     # Initial susceptible share
  nr0 = 0,                          # Initial recovered share
  nd0 = 0,                          # Initial dead share
  v = 31755,                        # Value of life (per capita)
  fx = 123,                         # FX multiplier (to USD)
  time_horizon = 1500,              # Simulation length (days)
  rng_seed = 1500,
  R0 = (3/10 + 1/7)/(1/7),
  pop_size = 1e6,                   # Will reset just below
  infect_thres = 1,                 # Fadeout threshold in persons
  utility_type = "Log",
  sigma = 0.19,                     # Only used if stochastic β is enabled
  bool_regular_sird = FALSE,
  bool_daily_cost_minimizing = FALSE
  # healthcare_capacity = NA,
  # excess_mortality_multiplier = 1
)

# scale population (dev default)
parameters$pop_size <- 1e4

# Initial state (fractions; kernel converts to counts)
initial_state <- c(
  Ns = parameters$ns0,
  Ni = parameters$ni0,
  Nr = parameters$nr0,
  Nd = parameters$nd0,
  HealthCost = 0,
  SocialActivityCost = 0,
  TotalCost = 0,
  a_t = NA,
  u_t = NA,
  Rt = NA
)

# Time grid
times <- seq(0, parameters$time_horizon, by = 1)

# ========== DETERMINISTIC BASELINE ==========
output_sim_deterministic <- run_sir_binomial(
  initial_state = initial_state,
  times = times,
  parameters = parameters,
  bool_stochastic_beta = FALSE,
  update_function = get_transitions_deterministic
)

det_peak_time <- which.max(output_sim_deterministic$Ni) - 1
cat(sprintf("Deterministic Peak Infection Time: %d\n", det_peak_time))

# ========== ONE BINOMIAL STOCHASTIC RUN ==========
set.seed(parameters$rng_seed)

output_stoch_binom_one <- run_one_simulation(
  initial_state = initial_state,
  times = times,
  parameters = parameters,
  bool_stochastic_beta = FALSE,                  # fixed β
  update_function = get_transitions_stochastic   # binomial draws for transitions
)

# Overlay plots + console summary
plot_sim_vs_det(
  output_stochastic   = output_stoch_binom_one,
  output_deterministic = output_sim_deterministic,
  plot_tag = "single_binomial"
)

# Resulting PDF: ./figures/single_binomial.pdf

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