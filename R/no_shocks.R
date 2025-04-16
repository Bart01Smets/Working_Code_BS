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
  utility_type = "Log"     # Utility type: "Log" or "Quadratic"
)

# Setting tolerance values
nifinal <- 10^-9
tolerance <- 10^-6

# Function to calculate the optimal action (a_t)
a_function <- function(alpha, beta, Ns, Ni, Lambda_s, Lambda_i, utility_type) {
  #denom <- max(abs(Lambda_s - Lambda_i), tolerance)  # Avoid division by zero
  denom <- max((Lambda_s - Lambda_i), tolerance)  # Avoid division by zero
 # if (utility_type == "Log") {
  #  sqrt_arg <- (Ns + Ni) + 4 * (1 + alpha) * beta * Ni * Ns * denom
  ##  return( (-(Ns + Ni) + sqrt(sqrt_arg)) / (2 * (1 + alpha) * beta * Ni * Ns * denom) )
  ##} else if (utility_type == "Quadratic") {
  ##  return( (Ns + Ni) / (Ns + Ni + (1 + alpha) * beta * Ni * Ns * denom) )
  ##}
#}


if (utility_type == "Log") {
  sqrt_arg <- (Ns + Ni) + 4 * (1 + alpha) * beta * Ni * Ns * denom
  return( (-(Ns + Ni) + sqrt(Ns + Ni) * sqrt(sqrt_arg)) / (2 * (1 + alpha) * beta * Ni * Ns * denom) )
  
} else if (utility_type == "Quadratic") {
  return( (Ns + Ni) / (Ns + Ni + (1 + alpha) * beta * Ni * Ns * denom) )
}
}
#   if (abs(denom) < tolerance)
#    {denom <- sign(denom) * tolerance } # Preserve sign, avoid zero

# Utility function
utility_function <- function(a_t, utility_type) {
  if (utility_type == "Log") {
    return(log(a_t) - a_t + 1)
  } else {
    return(-1/2 * (1 - a_t)^2)
  }
}


# Define the SIR model function including costates and costs
sir_costate_model <- function(time, state, parameters) { 
  with(as.list(c(state, parameters)), {
    
    # Extract the state variables
    Ns <- state["Ns"]        
    Ni <- state["Ni"]        
    Lambda_s <- state["Lambda_s"]
    Lambda_i <- state["Lambda_i"]
    HealthCost <- state["HealthCost"]
    SocialActivityCost <- state["SocialActivityCost"]
    TotalCost <- state["TotalCost"]
    
    # Calculate optimal action based on utility type
    a_t <- a_function(alpha, beta, Ns, Ni, Lambda_s, Lambda_i, utility_type)
    
    # Calculate utility of action
    u_t <- utility_function(a_t, utility_type)
    
    # Differential equations for SIR model with costates
    dNs <- -beta * a_t^2 * Ns * Ni
    dNi <- beta * a_t^2 * Ns * Ni - gamma * Ni
    dLambda_s <- (rho + delta) * Lambda_s - (u_t + (Lambda_i - Lambda_s) * beta * a_t^2 * Ni)
    dLambda_i <- (rho + delta) * Lambda_i - (u_t + alpha * (Lambda_i - Lambda_s) * beta * a_t^2 * Ns) + gamma * (kappa + Lambda_i)
    
    # Cost differential equations (per capita)
  #  dHealthCost <- (rho + delta) * HealthCost - gamma * kappa * Ni
 #   dSocialActivityCost <- (rho + delta) * SocialActivityCost + (Ns + Ni) * abs(u_t)
    #dTotalCost <- dHealthCost + dSocialActivityCost
    
    # Adjusted cost differential equations (per capita)
    #dHealthCost <- fx*exp(-(rho + delta) * time) * gamma * kappa * Ni#(Ns + Ni + Nr + Nd) fx
    #dSocialActivityCost <- fx*exp(-(rho + delta) * time) * (Ns + Ni) * abs(u_t)#/ (Ns + Ni + Nr + Nd)fx
    
    dHealthCost <- (rho + delta) * HealthCost - gamma * kappa * Ni
    dSocialActivityCost <- (rho + delta) * SocialActivityCost + (Ns + Ni) * abs(u_t)
    
    
    dTotalCost <- dHealthCost + dSocialActivityCost
    
    
    # Return the rate of change for each state variable
    return(list(c(dNs, dNi, dLambda_s, dLambda_i, dHealthCost, dSocialActivityCost, dTotalCost)))
  })
}

# Initial state variables including costs
initial_state <- c(
  Ns = parameters$ns0, 
  Ni = parameters$ni0, 
  Lambda_s = 0, 
  Lambda_i = 0, 
  HealthCost = 0, 
  SocialActivityCost = 0, 
  TotalCost = 0
)

# Time sequence for the simulation
time_sequence <- seq(0, 500, by = 1)  # Simulate for 500 days

# Solve the model
output <- ode(
  y = initial_state, 
  times = time_sequence, 
  func = sir_costate_model, 
  parms = parameters, 
  rtol = 1e-6, 
  atol = 1e-6
)

# Convert output to a data frame
output_df <- as.data.frame(output)

# Final costs
final_health_cost <- tail(output_df$HealthCost, 1)
final_social_activity_cost <- tail(output_df$SocialActivityCost, 1)
final_total_cost <- tail(output_df$TotalCost, 1)

# Print final costs
cat("Final Health Cost (per capita):", final_health_cost, "\n")
cat("Final Social Activity Cost (per capita):", final_social_activity_cost, "\n")
cat("Final Total Cost (per capita):", final_total_cost, "\n")
