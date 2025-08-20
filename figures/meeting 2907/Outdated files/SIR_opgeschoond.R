library(deSolve)
#install.packages("lamW")
library(lamW)
#install.packages("nleqslv")
#install.packages("pracma")

library(nleqslv)
library(pracma)

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

#find_nr0 <- function(ni0, beta, gamma) {
#  root_function <- function(nr0) {
#    ni0 - (1 - nr0 - exp(-beta / gamma * nr0))
#  }
  
  # Initial guess
 # initial_guess <- ni0 / (beta / gamma - 0.99)
  
  # Solve using uniroot()
#  result <- uniroot(root_function, lower = 0, upper = 1, tol = 1e-9)
#  return(result$root)
#}

# Compute nr0 using the uniroot method
#parameters$nr0 <- find_nr0(parameters$ni0, parameters$beta, parameters$gamma)
#parameters$ns0 <- 1 - parameters$ni0 - parameters$nr0

# Print values to verify
print(parameters$nr0)
print(parameters$ns0)


#find_nr0 <- function(ni0, beta, gamma) {
 # f <- function(nr0) {
#    ni0 - (1 - nr0 - exp(-beta / gamma * nr0))
 # }
  
  # Use nleqslv() for root solving
#  solution <- nleqslv(0.5, f)  # Initial guess = 0.5 (can be adjusted)
  
#  if (solution$termcd == 1) {
 #   return(solution$x)
#  } else {
#    stop("Root finding failed.")
 # }
#}

#find_nr0_lambert <- function(ni0, beta, gamma) {
 # lambert_term <- - (1 - ni0) * (beta / gamma)
#  nr0 <- -lambertW0(lambert_term) / (beta / gamma)
#  return(nr0)
#}

# Compute nr0 using the corrected method
#parameters$nr0 <- find_nr0_lambert(parameters$ni0, parameters$beta, parameters$gamma)
#parameters$ns0 <- 1 - parameters$ni0 - parameters$nr0

# Compute initial nr0 and ns0 based on equilibrium conditions
#parameters$nr0 <- find_nr0(parameters$ni0, parameters$beta, parameters$gamma)
#parameters$ns0 <- 1 - parameters$ni0 - parameters$nr0

#print(parameters$nr0)


sir_costate_model <- function(time, state, parameters) { 
  with(as.list(c(state, parameters)), {
    
    # Extract the state variables
    Ns <- state["Ns"]        
    Ni <- state["Ni"]        
    dNs <- -beta * Ns * Ni
    dNi <- beta * Ns * Ni - gamma * Ni
    dHealthCost <- (rho + delta) * HealthCost - gamma * kappa * Ni

      
     # fx * exp(-(rho + delta) * time) * gamma * kappa * (-Ni) Dit is niet-afgeleide kost
    
    # Return the rate of change for each state variable
    return(list(c(dNs, dNi, dHealthCost)))
  })
}

# Initial state
initial_state <- c(Ns = parameters$ns0, Ni = parameters$ni0, HealthCost=0)

# Time sequence for integration
time_sequence <- seq(0, 200, by = 1)

# Stopping condition function for `ode()`
#root_function <- function(time, state, parameters) {
 # Ni <- state["Ni"]
#  return(Ni - nifinal)  # Stops when Ni reaches `nifinal`
#}

#event_function <- function(time, state, parameters) {
#  if (state["Ni"] < nifinal) 
#  return(c(0,0,0))  # Stop condition met
 # else
 #   return(c(1,1,1))  # Continue integration
#}

print(parameters$nr0)
# Solve the ODE system with stopping condition
output <- ode(
  y = initial_state, 
  times = time_sequence, 
  func = sir_costate_model, 
  parms = parameters,
  #atol = 1e-6, 
  #rtol = 1e-6, 
  method = "lsoda",  # Adaptive ODE solver
#  rootfun= root_function,
 # events = list(func = event_function, root=TRUE)
)

# Convert output to data frame
output_df <- as.data.frame(output)

stop_time <- min(output_df$time[output_df$Ni < nifinal], na.rm = TRUE)
# Interpolate Ni(t) for cost integration
#ni_interp <- approxfun(output_df$time, output_df$Ni, rule = 2)

# Function to interpolate Ni(t) for cost integration
ni_interp <- approxfun(output_df$time, output_df$Ni, method = "linear", rule = 2)

# Compute total discounted health cost
discounted_health_cost <- integral(
  function(t) exp(-(parameters$rho + parameters$delta) * t) * 
    parameters$gamma * parameters$kappa * ni_interp(t),
  0, stop_time
)
print(tail(output_df[, c("time", "Ni", "HealthCost")], 200))

# Print the last 10 values of Infected and Recovered
print(tail(output_df[, c("time", "Ni")], 10))

# Compute Recovered as 1 - Susceptible - Infected
output_df$Nr <- 1 - output_df$Ns - output_df$Ni

# Print the last 10 values of Infected and Recovered
print(tail(output_df[, c("time", "Ni", "Nr")], 10))

# Compute Recovered as 1 - Susceptible - Infected
output_df$Nr <- 1 - output_df$Ns - output_df$Ni

# Print the full table of time, infected, and recovered
print(output_df[, c("time", "Ni", "Nr")])



# Extract health cost at t = 200
health_cost_at_200 <- output_df$HealthCost[output_df$time == 200]

# Print results
cat("Health Cost at t = 200:", health_cost_at_200, "\n")
cat("Total Discounted Health Cost:", discounted_health_cost, "\n")