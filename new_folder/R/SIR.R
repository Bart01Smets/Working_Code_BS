library(deSolve)

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
  export = FALSE,          # Export results flag
  utility_type = "Log"     # Utility type: "Log" or "Quadratic"
)

# Setting tolerance values
nifinal <- 10^-9
tolerance <- 10^-6

# Function to compute nr0 using uniroot
find_nr0 <- function(ni0, beta, gamma) {
  root_function <- function(nr0) {
    ni0 - (1 - nr0 - exp(-beta / gamma * nr0))
  }
  result <- uniroot(root_function, lower = 0, upper = 1, tol = 1e-9)
  return(result$root)
}

# Compute nr0 using the uniroot method
parameters$nr0 <- find_nr0(parameters$ni0, parameters$beta, parameters$gamma)
parameters$ns0 <- 1 - parameters$ni0 - parameters$nr0

# Print values to verify
print(parameters$nr0)
print(parameters$ns0)

# Define SIR model with cost tracking
sir_costate_model <- function(time, state, parameters) { 
  with(as.list(c(state, parameters)), {
    
    # Extract state variables
    Ns <- state["Ns"]        
    Ni <- state["Ni"]        
    
    # SIR dynamics
    dNs <- -beta * Ns * Ni
    dNi <- beta * Ns * Ni - gamma * Ni
    
    # Corrected HealthCost equation
    dHealthCost <- (rho + delta) * HealthCost - gamma * kappa * Ni

    # Return rate of change for each state variable
    return(list(c(dNs, dNi, dHealthCost)))
  })
}

# Initial state
initial_state <- c(Ns = parameters$ns0, Ni = parameters$ni0, HealthCost = 0)

# Time sequence for integration
time_sequence <- seq(0, 200, by = 1)

# Stopping condition function for `ode()`
root_function <- function(time, state, parameters) {
  Ni <- state["Ni"]
  return(Ni - nifinal)  # Stops when Ni reaches `nifinal`
}

event_function <- function(time, state, parameters) {
  if (state["Ni"] < nifinal) 
    return(c(0, 0, 0))  # Stop condition met
  else
    return(c(1, 1, 1))  # Continue integration
}

# Solve the ODE system with stopping condition
output <- ode(
  y = initial_state, 
  times = time_sequence, 
  func = sir_costate_model, 
  parms = parameters,
  method = "lsoda",  # Adaptive ODE solver
  rootfun = root_function,
  events = list(func = event_function, root = TRUE)
)

# Convert output to data frame
output_df <- as.data.frame(output)

# Get final total discounted health cost
total_discounted_health_cost <- tail(output_df$HealthCost, 1)

# Print results
cat("Total Discounted Health Cost:", total_discounted_health_cost, "\n")

# Print the last 10 values of Infected and HealthCost
print(tail(output_df[, c("time", "Ni", "HealthCost")], 10))
