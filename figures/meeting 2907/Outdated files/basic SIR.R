library(deSolve)
#install.packages("lamW")
library(lamW)

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

# Function to solve for nr0 using the equilibrium condition
find_nr0 <- function(ni0, beta, gamma) {
  root_function <- function(nr0) {
    ni0 - (1 - nr0 - exp(-beta/gamma * nr0))
  }
  result <- uniroot(root_function, lower = 0, upper = 1, tol = 1e-9)
  return(result$root)
}

#find_nr0_lambert <- function(ni0, beta, gamma) {
 # lambert_term <- - (1 - ni0) * (beta / gamma)
#  nr0 <- -lambertW0(lambert_term) / (beta / gamma)
#  return(nr0)
#}



# Compute initial nr0 and ns0 based on equilibrium conditions
parameters$nr0 <- find_nr0_lambert(parameters$ni0, parameters$beta, parameters$gamma)
parameters$ns0 <- 1 - parameters$ni0 - parameters$nr0

sir_costate_model <- function(time, state, parameters) { 
  with(as.list(c(state, parameters)), {
    
    # Extract the state variables
    Ns <- state["Ns"]        
    Ni <- state["Ni"]        
    #  Nr <- state["Nr"]        
    #    Nd <- state["Nd"]        
   # Lambda_s <- state["Lambda_s"]
  #  Lambda_i <- state["Lambda_i"]
  #  HealthCost <- state["HealthCost"]
  #  SocialActivityCost <- state["SocialActivityCost"]
  #  TotalCost <- state["TotalCost"]

    # Define the SIR model ODE system
    #sir_model <- function(time, state, parameters) {
    #  with(as.list(state), {
        dNs <- -beta * Ns * Ni
        dNi <- beta * Ns * Ni - gamma * Ni
     #   return(list(c(dNs, dNi)))
    #  })
    #}
    
    dHealthCost <- fx * exp(-(rho + delta) * time) * gamma * kappa * (-Ni)
  #  dSocialActivityCost <- fx * exp(-(rho + delta) * time) * (Ns + Ni) * (u_t) #abs
  #  dTotalCost <- dHealthCost + dSocialActivityCost
    
    
    # Return the rate of change for each state variable
    return(list(c(dNs, dNi, dHealthCost)))#dNr, dNd, , dSocialActivityCost, dTotalCost
  })
}


# Initial state
initial_state <- c(Ns = parameters$ns0, Ni = parameters$ni0, HealthCost=0)

# Time sequence for integration
time_sequence <- seq(0, 200, by = 1)

# Stopping condition function for `ode()`
root_function <- function(time, state, parameters) {
  Ni <- state["Ni"]
  return(Ni - nifinal)  # Stops when Ni reaches `nifinal`
}

# Custom function to stop integration when Ni < nifinal
#event_function <- function(time, state, parameters) {
 # if (state["Ni"] < nifinal) return(0)  # Stop condition met
#  else return(1)  # Continue integration
#}

# Solve the ODE system with stopping condition
output <- ode(
  y = initial_state, 
  times = time_sequence, 
  func = sir_costate_model, 
  parms = parameters,
  #atol = 1e-6, 
  #rtol = 1e-6, 
  method = "lsoda",  # Adaptive ODE solver
  rootfun= root_function
 # events = list(func = event_function, root=TRUE)
)

# Convert output to data frame
output_df <- as.data.frame(output)

stop_time <- min(output_df$time[output_df$Ni < nifinal], na.rm = TRUE)
# Interpolate Ni(t) for cost integration
ni_interp <- approxfun(output_df$time, output_df$Ni, rule = 2)

# Compute total discounted health cost via numerical integration
discounted_health_cost <- integrate(
  function(t) exp(-(parameters$rho + parameters$delta) * t) * 
    parameters$gamma * parameters$kappa * ni_interp(t),
  lower = 0,
  upper = stop_time
)$value

# Extract health cost at t = 200
health_cost_at_200 <- output_df$HealthCost[output_df$time == 200]

# Print results
cat("Health Cost at t = 200:", health_cost_at_200, "\n")
cat("Total Discounted Health Cost:", discounted_health_cost, "\n")