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
    sigma=1.689,            # Determines size of drift term ()
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
  utility_type = "Log",    # Utility type: "Log" or "Quadratic"
  rng_seed = 123
)

# RUN DETERMINISTIC - ODE solver   ####
########################################

# Setting tolerance values
parameters$nifinal <- 10^-9
parameters$tolerance <- 10^-6

# Initial state variables including separate costs
initial_state <- c(Ns = parameters$ns0, Ni = parameters$ni0, Nr = parameters$nr0, Nd = parameters$nd0, 
                   Lambda_s = 0, Lambda_i = 0, HealthCost = 0, SocialActivityCost = 0, TotalCost = 0,
                   a_t = NA, u_t = NA) # add a_t and u_t to keep track of this over time

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
cat("*** ODE SOLVER *** ","\n")
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
cat("*** LOOP-BASED ODE *** ","\n")
cat("Final Health Cost (per capita):", final_health_cost, "\n")
cat("Final Social Activity Cost (per capita):", final_social_activity_cost, "\n")
cat("Final Total Cost (per capita):", final_total_cost, "\n")

# COMPARE ODE VS LOOP SOLVERS   ####
########################################

head(output_pre_shock)
head(output_loop)

par(mfrow=c(2,2)) # get 4 sub-panels
plot(output_pre_shock_df$Ns,ylab='Ns')
lines(output_loop$Ns,col=2)
plot(output_pre_shock_df$Ns-output_loop$Ns,ylab='diff Ns')
plot(output_pre_shock_df$Ni-output_loop$Ni,ylab='diff Ni')
plot(output_pre_shock_df$TotalCost-output_loop$TotalCost,ylab='diff TotalCost')


# ALTERNATIVE FOR LOOP-BASED SOLVER   ####
########################################

# call function
output_sim <- run_sir_update(initial_state = initial_state, 
                             times = time_pre_shock, 
                             parameters = parameters)
# Print final costs
cat("*** USER-DEFINED *** ","\n")
cat("Final Health Cost (per capita):", output_sim$HealthCost[nrow(output_sim)], "\n")
cat("Final Social Activity Cost (per capita):", output_sim$SocialActivityCost[nrow(output_sim)], "\n")
cat("Final Total Cost (per capita):", output_sim$TotalCost[nrow(output_sim)], "\n")


# RUN MULTIPLE MODEL REALISATIONS   ####
########################################

num_experiments <- 10

# set random number generator seed
# redundant for deterministic modelling, but to be complete...
set.seed(parameters$rng_seed)

# init result matrix
output_experiments <- data.frame(matrix(NA,nrow=num_experiments,ncol=length(initial_state)))
names(output_experiments) <- names(initial_state)
dim(output_experiments)

# run multiple model realisations and store results
for(i_exp in 1:num_experiments){
  
  # get run-specific parameters
  parameters_exp <- parameters
  
  # example with slightly different beta per run
 # parameters_exp$beta <- parameters_exp$beta + rnorm(1,mean=0,sd=1e-3)
#   parameters_exp$beta <- max(0, parameters_exp$beta + parameters$sigma * rnorm(1, mean = 0, sd = 1))
  # run model
  output_sim <- run_sir_update(initial_state = initial_state, 
                               times = time_pre_shock, 
                               parameters = parameters_exp)
  # store results
  output_experiments[i_exp,] <- output_sim[nrow(output_sim),]
}

# inspect results
output_experiments

# Calculate  average of each column
column_means <- colMeans(output_experiments)


# Format the output as a data frame for a cleaner tabular display
column_means_df <- data.frame(Column = names(column_means), Mean = column_means)

# Print the formatted data frame
print(column_means_df)

# some graphical exploration
par(mfrow=c(1,1)) # reset sub-panels
boxplot(output_experiments$SocialActivityCost,main='SocialActivityCost')

# RUN STOCHASTIC MODEL REALISATIONS   ####
########################################

# get reference: deterministic model
output_sim_deterministic <- run_sir_update(initial_state = initial_state, 
                                       times = time_pre_shock, 
                                       parameters = parameters,
                                       infections_function = get_infections_determistic)
# define number of stochastic runs
num_experiments <- 10

# set random number generator seed
# redundant for deterministic modelling, but to be complete...
set.seed(parameters$rng_seed)

# init result matrix
output_experiments <- data.frame(matrix(NA,nrow=num_experiments,ncol=length(initial_state)))
names(output_experiments) <- names(initial_state)
dim(output_experiments)

# run multiple model realisations and store results
for(i_exp in 1:num_experiments){
  
  # run model
  output_sim <- run_sir_update(initial_state = initial_state, 
                               times = time_pre_shock, 
                               parameters = parameters,
                               infections_function = get_infections_stochastic_beta)
  # store results
  output_experiments[i_exp,] <- output_sim[nrow(output_sim),]
}

# inspect results
output_experiments

# Calculate  average of each column
column_means <- colMeans(output_experiments)


# Format the output as a data frame for a cleaner tabular display
column_means_df <- data.frame(Column = names(column_means), Mean = column_means)

# Print the formatted data frame
print(column_means_df)

# some graphical exploration
par(mfrow=c(2,2)) # reset sub-panels
xx <- boxplot(output_experiments$SocialActivityCost,main='SocialActivityCost')
abline(h=output_sim_deterministic$SocialActivityCost[nrow(output_sim_deterministic)],col=4)
text(x=0.6,y=output_sim_deterministic$SocialActivityCost[nrow(output_sim_deterministic)],
     labels=('deterministic'), pos=3,col=4)

# explore last simulation: infections
plot(output_sim$Ni,col=2,main='Infections',type='l',lwd=2) # infections, stochastic
lines(output_sim_deterministic$Ni,col=1) # infections, stochastic

# explore last simulation: all health states
plot(output_sim$Ns,main='Health states',ylim=0:1,type='l')
lines(output_sim$Ni,col=2)
lines(output_sim$Nr,col=3)
lines(output_sim$Nd,col=4)
legend('topright',
       c('Ns','Ni','Nr','Nd'),
       col = 1:4,
       lwd=2,
       ncol = 4)

