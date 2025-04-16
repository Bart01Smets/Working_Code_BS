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
library(ggplot2)
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
  ni0 = 0.0000527025093365949,         # Initial infected population (Ni0)
  ns0 = 0.999922203320773,         # Initial susceptible population (Ns0)
  nr0 = 0.00002484,        # Initial recovered population (Nr0)
  nd0 = 0.000000156,       # Initial dead population (Nd0)
  alpha =0,               # Altruism parameter (0 means optimal policy)
  fx = 123
 # scenario="optimal-policy" # Exchange rate multiplier for USD conversion
)

# tolerance values
nifinal <- 10^-9    # not used yet
tolerance <- 10^-6  # not used yet

# Optimal activity a(t) based on utility function
a_function <- function(alpha, beta, ns, ni, lambda_s, lambda_i) {

  diff_lambda <- lambda_s - lambda_i
  sqrt_term <- ns + ni + 4 * (1 + alpha) * beta * ni * ns * diff_lambda
  
  print(paste("diff_lambda =", diff_lambda))
  print(paste("sqrt_term =", sqrt_term))

  # defensive check
  if (sqrt_term < 0 ||diff_lambda == 0) {
    return(1)  # Prevent complex values or division by zero
  }

  num <- -(ns + ni) + sqrt(ns + ni) * sqrt(sqrt_term)
  denom <- 2 * (1 + alpha) * beta * ns * ni * diff_lambda
  return(num/denom)  # Ensure activity remains between [0,1] return(max(min(num / denom, 1), 0))
}


# Utility function
# utility_type = "log"
utility_function <- function(a_t) {
    return(log(a_t) - a_t + 1)
}

# explore activity-utility relationship
a_t <- seq(0,1,0.01)
plot(a_t,utility_function(seq(0,1,0.01)),ylab='utility')

## for debugging
# attach(parameters)
# attach(as.list(initial_state))

# Define the SIR model function including costates, recovered, dead, and separate costs
sir_costate_model <- function(time, state, parameters) { 
  with(as.list(c(state, parameters)), { # this statements enables the use of the parameter and state names
     
    # Debugging print statements
    print(paste("t =", time))
    print(paste("Ns =", Ns))
    print(paste("Ni =", Ni))
    print(paste("Lambda_s =", Lambda_s))
   print(paste("Lambda_i =", Lambda_i))
   
    # Calculate optimal action based on utility type
    a_t <- a_function(alpha, beta, Ns, Ni, Lambda_s, Lambda_i)

    
    # stopping condition
    if(Ni < nifinal){
      return(list(rep(0,12)))
    }
 
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
    dSocialActivityCost <- -fx * (exp(-(rho + delta) * time) * (Ns + Ni) * u_a)  # note: negative
    
    # Total cost is the sum of cumulative costs
    dTotalCost <- dHealthCost + dSocialActivityCost

    # Reproduction number
    R0 = beta / gamma
    R_t = R0 * a_t^2 * Ns
    
    
    print(paste("a_t =", a_t))
    # Return the rate of change for each state variable
    return(list(c(dNs, dNi, dNr, dNd, dLambda_s, dLambda_i, dHealthCost, dSocialActivityCost, dTotalCost, a_t, u_t , R_t)))
  })
}

#a_t2<-a_function(parameters$alpha, parameters$beta, parameters$ns0, parameters$ni0, -102.859, -194.343)
u_t0<-utility_function(0.9978751)

# Initial state variables including separate costs
initial_state <- c(Ns = parameters$ns0, 
                   Ni = parameters$ni0, 
                   Nr = parameters$nr0, 
                   Nd = parameters$nd0,
                   Lambda_s = -102.859768242236,#Lambda_values[1], -102.859	-194.343, -164.924, -29822.368
                   Lambda_i = -194.343405878049,#Lambda_values[2], -165.651	-21985.5 #-165.618 -22490.5
                   HealthCost = 0, 
                   SocialActivityCost = 0, 
                   TotalCost = 0,
                   a_t = 0.997873982271116, #0.61904
                   u_t = u_t0,#-0101193
                   R_t = 0) # Lambda_s = -165,#Lambda_values[1], -102.859	-194.343 Lambda_i = -21985.5,#Lambda_values[2], -165.651	-21985.5
#initial_state <- as.numeric(initial_state)
#names(initial_state) <- c("Ns", "Ni", "Nr", "Nd", "Lambda_s", "Lambda_i", 
             #             "HealthCost", "SocialActivityCost", "TotalCost", "a_t", "u_t", "R_t")

# Adjust values dynamically based on `alpha`
#if (parameters$alpha == 1) {
#  initial_state$Lambda_s <- -165.651
#  initial_state$Lambda_i <- -21958.5
#}

# Convert list to named numeric vector for `ode()`
#initial_state <- unlist(initial_state)
# Time sequence for pre-shock
time_pre_shock <- seq(0, 20, by = 1)

# Solve the model
#output_intervention <- ode(y = initial_state, 
       #                    times = time_pre_shock, 
         #                  func  = sir_costate_model, 
         #                  parms = parameters,
                           #method= "lsoda",#method="rk4" "bdf" "euler"
         #                  method="lsoda"
                           output_intervention <- ode(
                             y = initial_state, 
                             times = time_pre_shock, 
                             func = sir_costate_model, 
                             parms = parameters,
                             method = "lsoda"   # lsoda automatically handles stiff/non-stiff transitions
                                    # tighten absolute tolerance
                                   # tighten relative tolerance
                                    # limit maximum step size
                                 # allow plenty of steps to avoid premature termination
                           )
#) #tcrit = time_pre_shock)
print(output_intervention)
output_intervention_df <- as.data.frame(output_intervention)
print(output_intervention_df$a_t)


#output_intervention_df$a_t <- c(output_intervention_df$a_t[-1], output_intervention_df$a_t[nrow(output_intervention_df)])


#output_intervention_df <- output_intervention_df[output_intervention_df$time %% 1 == 0, ]

print(unique(output_intervention_df$time))
output_intervention_df$a_t <- c(output_intervention_df$a_t[1],diff(output_intervention_df$a_t))
print(output_intervention_df$a_t)
output_intervention_df$u_t <- c(output_intervention_df$u_t[1],diff(output_intervention_df$u_t))
print(output_intervention_df$u_t)



#output_intervention_df$R_t <- c(output_intervention_df$R_t[1],diff(output_intervention_df$R_t))


print(output_intervention_df[, c("time", "Ns", "Ni", "a_t", "u_t", "Lambda_s", "Lambda_i")])
# Install and load the writexl package if not already installed

library(writexl)

# Prepare the data frame with key variables: time, Ns, Ni, activity (a_t), utility (u_t), Lambda_s, and Lambda_i
export_data <- output_intervention_df[, c("time", "Ns", "Ni", "a_t", "u_t", "Lambda_s", "Lambda_i")]
# Define the file path
file_path <- "C:/Users/Bart Smets/OneDrive/Documenten/GitHub/github map/SIR_OutputLF.xlsx"
# Export the data frame to an Excel file
write_xlsx(export_data, file_path)
print(export_data)


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

#utility and activity values plot
p4 <- ggplot(output_intervention_df, aes(x = time)) +
  geom_line(aes(y = a_t, color = "a_t")) +
  geom_line(aes(y = u_t, color = "u_T")) +
  labs(y = "a_t and u_t", title = "Activity and Utility loss") +
  theme_minimal()

#utility and activity values plot
p5 <- ggplot(output_intervention_df, aes(x = time)) +
  geom_line(aes(y = R_t, color = "R_t")) 
  labs(y = "R_t", title = "Reproduction number") +
  theme_minimal()

 
  output_intervention_df_limited<-subset(output_intervention_df,time<=100)
  
  # Plot Lambda_s and Lambda_i over time
  ggplot(output_intervention_df_limited, aes(x = time)) +
    geom_line(aes(y = Lambda_s, color = "Lambda_s")) +
    geom_line(aes(y = Lambda_i, color = "Lambda_i")) +
    labs(y = "Lambda values", title = "Costate Variables Over Time", x = "Time (days)") +
    scale_color_manual(name = "Legend", values = c("Lambda_s" = "blue", "Lambda_i" = "red")) +
    theme_minimal()
  
# Arrange the plots
grid.arrange(p1, p2, p3, p4, p5, ncol = 2)

# Print final costs and costates
output_final <- tail(output_intervention_df, 1)
cat("Final Lambda_s:", output_final$Lambda_s, "\n")
cat("Final Lambda_i:", output_final$Lambda_i, "\n")
cat("Final Health Cost (per capita):", output_final$HealthCost, "\n")
cat("Final Social Activity Cost (per capita):", output_final$SocialActivityCost, "\n")
cat("Final Total Cost (per capita):", output_final$TotalCost, "\n")




