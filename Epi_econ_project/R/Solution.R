#
# Script for the EPI-ECON modelling work: SIR with intervention (optimal policy)
#
# - Reproduction of Farboodi et al (2021) code: "Solution" 
# - Disease dynamics with intervention
# - Laissez-Faire: Mathematica code with alpha=1 shows Health Cost: 11374, and Activity Cost: 1277
# - optimal policy: 
#      - Mathematica code with alpha=0 shows Health Cost: 860, and Activity Cost: 7227
#      - manuscript shows a(t)=0.6
#
######################################################### #

rm(list=ls())

library(deSolve)

# Plot results
library(ggplot2)
library(gridExtra)

# Parameters
parameters <- list(
  gamma = 1/7,             # Recovery rate (γ)
  beta = 3/10 + 1/7,       # Transmission rate (β)
  delta = 0.67 / 365,      # Arrival rate of cure (δ)
  rho = 0.05 / 365,        # Discounting rate (ρ)
  pi = 0.0062,             # Infection fatality rate (π)
  kappa = 197,             # Expected cost of infection (κ)
  ni0 = 0.0000527,         # Initial infected population (Ni0)
  ns0 = 0.9999223,         # Initial susceptible population (Ns0)
  nr0 = 0.00002484,        # Initial recovered population (Nr0)
  nd0 = 0.000000156,       # Initial dead population (Nd0)
  alpha = 1,               # Altruism parameter (0 means no altruism)
  fx = 123                 # Exchange rate multiplier for USD conversion
)

# tolerance values
nifinal <- 10^-9
tolerance <- 10^-6

# Optimal activity a(t) based on utility function
a_function <- function(alpha, beta, ns, ni, Lambda_s, Lambda_i) {
  #return(1)
  return(0.60)
}

# Utility function
# utility_type = "log"
utility_function <- function(a_t) {
    return(log(a_t) - a_t + 1)
}

# Define the SIR model function including costates, recovered, dead, and separate costs
sir_costate_model <- function(time, state, parameters) { 
  with(as.list(c(state, parameters)), {
    
    # Calculate optimal action based on utility type
    a_t <- a_function(alpha, beta, Ns, Ni, Lambda_s, Lambda_i)
 
    # Calculate utility of action
    u_a <- utility_function(a_t)
  
    # Differential equations for SIR model with costates
    dNs <- -beta * a_t^2 * Ns * Ni
    dNi <- beta * a_t^2 * Ns * Ni  - gamma * Ni 
    dNr <- gamma * Ni * (1 - pi)
    dNd <- gamma * Ni * pi
    
    # # Costate equations
    dLambda_s <- (rho + delta) * Lambda_s - u_a - (Lambda_i - Lambda_s) * beta * a_t^2 * Ni
    dLambda_i <- (rho + delta) * Lambda_i - u_a - alpha * (Lambda_i - Lambda_s) * beta * a_t^2 * Ns + gamma * (kappa + Lambda_i)
    
    # cost calculations
    dHealthCost <- fx * (exp(-(rho + delta) * time) * gamma * kappa * Ni)
    dSocialActivityCost <- -fx * (exp(-(rho + delta) * time) * (Ns + Ni) * u_a)
    
    # Total cost is the sum of cumulative costs
    dTotalCost <- dHealthCost + dSocialActivityCost

    # Return the rate of change for each state variable
    return(list(c(dNs, dNi, dNr, dNd, dLambda_s, dLambda_i, dHealthCost, dSocialActivityCost, dTotalCost, a_t, u_t)))
  })
}

# Compute terminal conditions dynamically
compute_terminal_conditions <- function(ns_T, ni_T, parameters) {
  Lambda_i_T <- - (parameters$gamma * parameters$kappa * (parameters$delta + parameters$rho + parameters$beta * ni_T)) /
    (parameters$alpha * parameters$beta^2 * ni_T * ns_T + (parameters$delta + parameters$rho + parameters$beta * ni_T) * 
       (parameters$rho + parameters$delta + parameters$gamma - parameters$alpha * parameters$beta * ns_T))
  
  Lambda_s_T <- (parameters$beta * Lambda_i_T * ni_T) / (parameters$rho + parameters$delta + parameters$beta * ni_T)
  
  return(c(lambda_s = Lambda_s_T, lambda_i = Lambda_i_T))
}
Lambda_values <- compute_terminal_conditions(ns_T = parameters$ns0, ni_T = parameters$ni0, parameters = parameters)
Lambda_values

# Initial state variables including separate costs
initial_state <- c(Ns = parameters$ns0, 
                   Ni = parameters$ni0, 
                   Nr = parameters$nr0, 
                   Nd = parameters$nd0,
                   Lambda_s = Lambda_values[1], 
                   Lambda_i = Lambda_values[2], 
                   HealthCost = 0, 
                   SocialActivityCost = 0, 
                   TotalCost = 0,
                   a_t = 1, 
                   u_t = 0) 

# Time sequence for pre-shock
time_pre_shock <- seq(0, 610, by = 1)

# Solve the model
output_intervention <- ode(y = initial_state, 
                           times = time_pre_shock, 
                           func  = sir_costate_model, 
                           parms = parameters)
output_intervention_df <- as.data.frame(output_intervention)

# Main SIR plot
p1 <- ggplot(output_intervention_df, aes(x = time)) +
  geom_line(aes(y = Ns, color = "Susceptible")) +
  geom_line(aes(y = Ni, color = "Infected")) +
  geom_line(aes(y = Nr, color = "Recovered")) +
  labs(y = "Proportion", title = "SIR Model with Optimal Control") +
  theme_minimal()

# Lambda values plot
p2 <- ggplot(output_intervention_df, aes(x = time)) +
  geom_line(aes(y = Lambda_s, color = "Lambda_s")) +
  geom_line(aes(y = Lambda_i, color = "Lambda_i")) +
  labs(y = "Lambda Values", title = "Costate Variables") +
  theme_minimal()

#Cost values plot
p3 <- ggplot(output_intervention_df, aes(x = time)) +
  geom_line(aes(y = HealthCost, color = "HealthCost")) +
  geom_line(aes(y = SocialActivityCost, color = "SocialActivityCost")) +
  geom_line(aes(y = TotalCost, color = "TotalCost")) +
  labs(y = "Cost", title = "Cost results") +
  theme_minimal()

output_intervention_df$a_t <- c(output_intervention_df$a_t[1],diff(output_intervention_df$a_t))
output_intervention_df$u_t <- c(output_intervention_df$u_t[1],diff(output_intervention_df$u_t))

#utility and activity values plot
p4 <- ggplot(output_intervention_df, aes(x = time)) +
  geom_line(aes(y = a_t, color = "a_t")) +
  geom_line(aes(y = u_t, color = "u_T")) +
  labs(y = "a_t and u_t", title = "Activity and Utility loss") +
  theme_minimal()

# Arrange the plots
grid.arrange(p1, p2, p3, p4, ncol = 2)

# Print final costs and costates
output_final <- tail(output_intervention_df, 1)
cat("Final Lambda_s:", output_final$Lambda_s, "\n")
cat("Final Lambda_i:", output_final$Lambda_i, "\n")
cat("Final Health Cost (per capita):", output_final$HealthCost, "\n")
cat("Final Social Activity Cost (per capita):", output_final$SocialActivityCost, "\n")
cat("Final Total Cost (per capita):", output_final$TotalCost, "\n")
