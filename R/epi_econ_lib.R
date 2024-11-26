#
# Help function for the EPI-ECON modelling work
#
# - Deterministic and stochastic transitions
################################################################ #

library(deSolve)

# Function to calculate the optimal action (a_t)
a_function <- function(alpha, beta, Ns, Ni, Lambda_s, Lambda_i, utility_type, tolerance) {
  denom <- max(abs(Lambda_s - Lambda_i), tolerance)  # Avoid division by zero
  
  if (utility_type == "Log") {
    sqrt_arg <- (Ns + Ni) + 4 * (1 + alpha) * beta * Ni * Ns * denom
    return( (-(Ns + Ni) + sqrt(sqrt_arg)) / (2 * (1 + alpha) * beta * Ni * Ns * denom) )
  } else if (utility_type == "Quadratic") {
    return( (Ns + Ni) / (Ns + Ni + (1 + alpha) * beta * Ni * Ns * denom) )
  }
}

# Utility function
utility_function <- function(a_t, utility_type) {
  if (utility_type == "Log") {
    return(log(a_t) - a_t + 1)
  } else {
    return(-1/2 * (1 - a_t)^2)
  }
}

# Define the SIR model function including costates, recovered, dead, and separate costs
sir_costate_model <- function(time, state, parameters) { 
  with(as.list(c(state, parameters)), {
    
    # Extract the state variables
    Ns <- state["Ns"]        
    Ni <- state["Ni"]        
    Nr <- state["Nr"]        
    Nd <- state["Nd"]        
    Lambda_s <- state["Lambda_s"]
    Lambda_i <- state["Lambda_i"]
    HealthCost <- state["HealthCost"]
    SocialActivityCost <- state["SocialActivityCost"]
    TotalCost <- state["TotalCost"]
    
    # Calculate optimal action based on utility type
    a_t <- a_function(alpha, beta, Ns, Ni, Lambda_s, Lambda_i, utility_type, tolerance)
    
    # Calculate utility of action
    u_t <- utility_function(a_t, utility_type)
    
    # Differential equations for SIR model with costates
    dNs <- -beta * a_t^2 * Ns * Ni
    dNi <- beta * a_t^2 * Ns * Ni - gamma * Ni - pi * gamma * Ni
    dNr <- gamma * Ni * (1 - pi)
    dNd <- gamma * Ni * pi
    dLambda_s <- (rho + delta) * Lambda_s - (u_t + (Lambda_i - Lambda_s) * beta * a_t^2 * Ni)
    dLambda_i <- (rho + delta) * Lambda_i - (u_t + alpha * (Lambda_i - Lambda_s) * beta * a_t^2 * Ns) - gamma * (kappa + Lambda_i)
    
    # Adjusted cost differential equations (per capita)
    dHealthCost <- (fx / (Ns + Ni + Nr + Nd)) * exp(-(rho + delta) * time) * gamma * kappa * Ni
    dSocialActivityCost <- (fx / (Ns + Ni + Nr + Nd)) * exp(-(rho + delta) * time) * (Ns + Ni) * abs(u_t)
    dTotalCost <- dHealthCost + dSocialActivityCost
    
    # Return the rate of change for each state variable
    return(list(c(dNs, dNi, dNr, dNd, dLambda_s, dLambda_i, dHealthCost, dSocialActivityCost, dTotalCost)))
  })
}

