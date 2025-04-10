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
  ns0 = 0.9999223,     # Initial susceptible population (Ns0)
  nr0 = 0.00002484,                 # Initial recovered population (Nr0)
  nd0 = 0.000000156,                 # Initial dead population (Nd0)
  alpha = 0,               # Altruism parameter (default set to no altruism)
  fx = 123,                # Exchange rate multiplier for USD conversion
  T1 = 500,                # Time of shock
  T2 = 1,                  # Duration after shock
  z = 1,                   # Shock multiplier for β
  export = FALSE,          # Export results flag
  utility_type = "Log"     # Utility type: "Log" or "Quadratic"
)

# Setting tolerance values
nifinal <- 10^-9
tolerance <- 10^-6

# Function to calculate the optimal action (a_t)

a_function <- function(alpha, beta, Ns, Ni, Lambda_s, Lambda_i, utility_type) {
  denom <- max(abs(Lambda_s - Lambda_i), tolerance)  # Avoid division by zero
  #  denom <- Lambda_s - Lambda_i
  
  
  # if (abs(denom) < tolerance)
  #    {denom <- sign(denom) * tolerance } 
  
  
  if (utility_type == "Log") {
    sqrt_arg <- (Ns + Ni) + 4 * (1 + alpha) * beta * Ni * Ns * denom
    return( (-(Ns + Ni) + sqrt(Ns + Ni) * sqrt(sqrt_arg)) / (2 * (1 + alpha) * beta * Ni * Ns * denom) )
    
    
    
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
    #  Nr <- state["Nr"]        
    #    Nd <- state["Nd"]        
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
    dNi <- beta * a_t^2 * Ns * Ni  - gamma * Ni #- pi * gamma * Ni
    #dNr <- gamma * Ni * (1 - pi)
    #   dNd <- gamma * Ni * pi
    dLambda_s <- (rho + delta) * Lambda_s - (u_t + (Lambda_i - Lambda_s) * beta * a_t^2 * Ni)
    
    dLambda_i <- (rho + delta) * Lambda_i - (u_t + alpha * (Lambda_i - Lambda_s) * beta * a_t^2 * Ns) + gamma * (kappa + Lambda_i)
    
    # Instantaneous cost calculations (per capita)
  #  dHealthCost <- fx * exp(-(rho + delta) * time) * gamma * kappa * (-Ni)
   # dSocialActivityCost <- fx * exp(-(rho + delta) * time) * (Ns + Ni) * (u_t) #abs
    
    dHealthCost <- fx*((rho + delta) * HealthCost - gamma * kappa * (-Ni))
    dSocialActivityCost <- fx*((rho + delta) * SocialActivityCost + (Ns + Ni) * u_t)
    
    # Total cost is the sum of cumulative costs
    dTotalCost <- dHealthCost + dSocialActivityCost
    #final_health_cost <- sum(exp(-(rho + delta) * output_combined_df$time) * 
    #  gamma * kappa * output_combined_df$Ni)
    
    
    # Return the rate of change for each state variable
    return(list(c(dNs, dNi, dLambda_s, dLambda_i, dHealthCost, dSocialActivityCost, dTotalCost)))#dNr, dNd
  })
}

# Initial state variables including separate costs
initial_state <- c(Ns = parameters$ns0, Ni = parameters$ni0, 
                   Lambda_s = 0, Lambda_i = 0, HealthCost = 0, SocialActivityCost = 0, TotalCost = 0)#Nr = parameters$nr0, Nd = parameters$nd0, 

# Time sequence for pre-shock
time_pre_shock <- seq(0, parameters$T1, by = 1)

# Solve the model before the shock
output_pre_shock <- ode(y = initial_state, times = time_pre_shock, func = sir_costate_model, parms = parameters)
output_pre_shock_df <- as.data.frame(output_pre_shock)

# Adjust the transmission rate by applying the shock
parameters$beta <- parameters$z * parameters$beta

# Extract state at T1 for post-shock initial conditions including costs
Ns_T1 <- tail(output_pre_shock_df$Ns, 1)
Ni_T1 <- tail(output_pre_shock_df$Ni, 1)
#Nr_T1 <- tail(output_pre_shock_df$Nr, 1)
#Nd_T1 <- tail(output_pre_shock_df$Nd, 1)
Lambda_s_T1 <- tail(output_pre_shock_df$Lambda_s, 1)
Lambda_i_T1 <- tail(output_pre_shock_df$Lambda_i, 1)
HealthCost_T1 <- tail(output_pre_shock_df$HealthCost, 1)
SocialActivityCost_T1 <- tail(output_pre_shock_df$SocialActivityCost, 1)
TotalCost_T1 <- tail(output_pre_shock_df$TotalCost, 1)

# Initial state for post-shock period including costs
initial_state_shock <- c(Ns = Ns_T1, Ni = Ni_T1,  
                         Lambda_s = Lambda_s_T1, Lambda_i = Lambda_i_T1, 
                         HealthCost = HealthCost_T1, SocialActivityCost = SocialActivityCost_T1, TotalCost = TotalCost_T1) #Nr = Nr_T1, Nd = Nd_T1,

# Time sequence for post-shock period
time_post_shock <- seq(parameters$T1, parameters$T1 + parameters$T2, by = 1)

# Solve the model after the shock
output_post_shock <- ode(y = initial_state_shock, times = time_post_shock, func = sir_costate_model, parms = parameters, rtol=1e-4, atol=1e-4)
output_post_shock_df <- as.data.frame(output_post_shock)

# Combine results of pre-shock and post-shock phases
output_combined_df <- rbind(output_pre_shock_df, output_post_shock_df[-1, ])  # Remove overlapping initial state row

# Apply terminal conditions for costates
Ns_T <- tail(output_combined_df$Ns, 1)
Ni_T <- tail(output_combined_df$Ni, 1)
Lambda_i_T <- tail(output_combined_df$Lambda_i, 1)

Lambda_s_T <- (parameters$beta * Lambda_i_T * Ni_T) / (parameters$rho + parameters$delta + parameters$beta * Ni_T)
Lambda_i_T <- - (parameters$gamma * parameters$kappa * (parameters$delta + parameters$rho + parameters$beta * Ni_T)) /
  (parameters$alpha * parameters$beta^2 * Ni_T * Ns_T + 
  (parameters$delta + parameters$rho + parameters$beta * Ni_T) * 
  (parameters$rho + parameters$delta + parameters$gamma - parameters$alpha * parameters$beta * Ns_T))

# Append terminal conditions
output_combined_df$Lambda_s[nrow(output_combined_df)] <- Lambda_s_T
output_combined_df$Lambda_i[nrow(output_combined_df)] <- Lambda_i_T

# Print final costs and costates
cat("Final Lambda_s:", Lambda_s_T, "\n")
cat("Final Lambda_i:", Lambda_i_T, "\n")
cat("Final Health Cost (per capita):", tail(output_combined_df$HealthCost, 1), "\n")
cat("Final Social Activity Cost (per capita):", tail(output_combined_df$SocialActivityCost, 1), "\n")
cat("Final Total Cost (per capita):", tail(output_combined_df$TotalCost, 1), "\n")
