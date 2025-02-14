#
# Script for the EPI-ECON modelling work: SIR dynamics
#
# - Reproduction of Farboodi et al (2021) code: "SIR"
# - Disease dynamics without intervention
# - Mathematica code shows Health Cost and Total Cost: 172.541 utils
# - Manuscript shows Health Cost and Total Cost in dollars of $21,223

#
######################################################### #

library(deSolve)

# Parameters
parameters <- list(
  gamma = 1/7,             # Recovery rate (γ)
  beta = 3/10 + 1/7,       # Transmission rate (β)
  delta = 0.67 / 365,      # Arrival rate of cure (δ)
  rho = 0.05 / 365,        # Discounting rate (ρ)
  kappa = 197,             # Expected cost of infection (κ)
  ni0 = 0.0000527,         # Initial infected population (Ni0)
  ns0 = 0.9999223,         # Initial susceptible population (Ns0)  (Table 1, page 15)
  fx = 123                 # To convert utils to US dollars 
)

# Basic reproduction number R0 = beta / gamma
# Table 1, page 15: R0 = 3.1
parameters$beta / parameters$gamma

# boundaries
nifinal = 10^-9;

# Define SIR model with cost tracking
sir_model <- function(time, state, parameters) { 
  with(as.list(c(state, parameters)), {
    
    # Extract state variables
    Ns <- state["Ns"]        
    Ni <- state["Ni"]        
    
    # stopping condition
    if(Ni < nifinal){
      return(list(rep(0,3)))
    }
    
    # SIR dynamics
    dNs <- -beta * Ns * Ni
    dNi <- beta * Ns * Ni - gamma * Ni
    
    # HealthCost equation
    dHealthCost <- exp(-(rho + delta) * time) * gamma * kappa * Ni

    # Return rate of change for each state variable
    return(list(c(dNs, dNi, dHealthCost)))
  })
}

# Initial state
initial_state <- c(Ns = parameters$ns0, Ni = parameters$ni0, HealthCost = 0)

# Time sequence for integration
time_sequence <- seq(0, 200, by = 1)

# Solve the ODE system with stopping condition
output <- ode(
  y = initial_state, 
  times = time_sequence, 
  func = sir_model, 
  parms = parameters,
  method = "lsoda"
)

# Convert output to data frame
output_df <- as.data.frame(output)

# Get final total discounted health cost
total_discounted_health_cost <- tail(output_df$HealthCost, 1)

# Print results
cat("Total Discounted Health Cost:", total_discounted_health_cost, " (utils)\n")
cat("Total Discounted Health Cost:", parameters$fx * total_discounted_health_cost, " (dollar)\n")

