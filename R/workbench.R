#
# Workbench script for the EPI-ECON modelling work
#
# - Deterministic and stochastic modelling
#
######################################################### #

# clear workbench
rm(list=ls())

# load functions
source("R/epi_econ_lib.R")

# SETUP   ####
####################

# Parameters
parameters <- list(
  gamma = 1/7,             # Recovery rate (γ)
  beta = 3/10 + 1/7,       # Transmission rate (β)
  delta = 0.67 / 365,      # Arrival rate of cure (δ)
  rho = 0.05 / 365,        # Discounting rate (ρ)
  pi = 0.0062,             # Infection fatality rate (π)
  kappa = 197,             # Expected cost of infection (κ)
  ni0 = 0.0000527,         # Initial infected population (Ni0)
  ns0 = 1 - 0.0000527,     # Initial susceptible population (Ns0)
  nr0 = 0,                 # Initial recovered population (Nr0)
  nd0 = 0,                 # Initial dead population (Nd0)
  alpha = 0,               # Altruism parameter (default set to no altruism)
  fx = 123,                # Exchange rate multiplier for USD conversion
  T1 = 500,                # Time of shock
  T2 = 1,                  # Duration after shock
  z = 1,                   # Shock multiplier for β
  export = FALSE,          # Export results flag
  utility_type = "Log"     # Utility type: "Log" or "Quadratic"
)

# RUN DETERMINISTIC - ODE solver   ####
########################################

# Setting tolerance values
parameters$nifinal <- 10^-9
parameters$tolerance <- 10^-6

# Initial state variables including separate costs
initial_state <- c(Ns = parameters$ns0, Ni = parameters$ni0, Nr = parameters$nr0, Nd = parameters$nd0, 
                   Lambda_s = 0, Lambda_i = 0, HealthCost = 0, SocialActivityCost = 0, TotalCost = 0)

# Time sequence for pre-shock
time_pre_shock <- seq(0, parameters$T1, by = 1)

# Solve the model
output_pre_shock <- ode(y = initial_state, times = time_pre_shock, func = sir_costate_model, parms = parameters)
output_pre_shock_df <- as.data.frame(output_pre_shock)

# Final Health Cost, Social Activity Cost, and Total Cost
final_health_cost <- tail(output_pre_shock_df$HealthCost, 1)
final_social_activity_cost <- tail(output_pre_shock_df$SocialActivityCost, 1)
final_total_cost <- tail(output_pre_shock_df$TotalCost, 1)

# Print final costs
cat("Final Health Cost (per capita):", final_health_cost, "\n")
cat("Final Social Activity Cost (per capita):", final_social_activity_cost, "\n")
cat("Final Total Cost (per capita):", final_total_cost, "\n")

# RUN DETERMINISTIC - FOR loop by day   ####
########################################

# Solve the model
output_loop <- run_sir_costate_model(y = initial_state, times = time_pre_shock, func = sir_costate_model, parms = parameters)
output_loop_df <- data.frame(output_loop)

# Final Health Cost, Social Activity Cost, and Total Cost
final_health_cost <- tail(output_loop_df$HealthCost, 1)
final_social_activity_cost <- tail(output_loop_df$SocialActivityCost, 1)
final_total_cost <- tail(output_loop_df$TotalCost, 1)

# Print final costs
cat("Final Health Cost (per capita):", final_health_cost, "\n")
cat("Final Social Activity Cost (per capita):", final_social_activity_cost, "\n")
cat("Final Total Cost (per capita):", final_total_cost, "\n")

# COMPARE ODE VS LOOP SOLVERS   ####
########################################

head(output_pre_shock)
head(output_loop)

par(mfrow=c(2,2))
plot(output_pre_shock_df$Ns)
lines(output_loop$Ns,col=2)
plot(output_pre_shock_df$Ns-output_loop$Ns)
plot(output_pre_shock_df$Ni-output_loop$Ni)
plot(output_pre_shock_df$TotalCost-output_loop$TotalCost)

