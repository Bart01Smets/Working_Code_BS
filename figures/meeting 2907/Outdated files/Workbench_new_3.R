#
# Workbench script for the EPI-ECON modelling work
#
# - Deterministic and stochastic modelling
#
######################################################### #
getwd()



setwd("C:/Users/Bart Smets/OneDrive/Documenten/GitHub/Working_Code_BS")


# clear workbench
rm(list=ls())

# load functions
source("R/epi_econ_lib_new3.R")

# SETUP   ####
####################

# Parameters
parameters <- list(
  gamma = 1/7,             # Recovery rate (γ)
  beta = 3/10 + 1/7,       # Transmission rate (β)
  delta = 0.67 / 365,      # Arrival rate of cure (δ)
  rho = 0.05 / 365,        # Discounting rate (ρ)
  pi = 0.0062,             # Infection fatality rate (π)
  sigma = 0.10,            # Determines stochastic beta_t =  N(beta,sigma)
  kappa = 197,             # Expected cost of infection (κ)
  ni0 = 0.000527,         # Initial infected population (Ni0)
  ns0 = 1 - 0.000527,     # Initial susceptible population (Ns0)
  nr0 = 0,                 # Initial recovered population (Nr0)
  nd0 = 0, # Initial dead population (Nd0)
  v = 31755,
  alpha = 1,               # Altruism parameter (default set to no altruism)
  fx = 123,                # Exchange rate multiplier for USD conversion
  time_horizon = 500,      # Time of shock
  utility_type = "Log",    # Utility type: "Log" or "Quadratic"
  rng_seed = 123,
  scenario = "laissez-faire",
  R0= (3/10+1/7)/(1/7),
  infect_thres = 1
)
print(parameters$v)

# define population size
parameters$pop_size <- 1e4

# define number of stochastic runs
num_experiments <- 10

# define fadeout threshold
fadeout_threshold = 100

# define default parameters (for development)
plot_tag <- 'dev'
bool_stochastic_beta <- FALSE
update_function <- get_transitions_stochastic

# RUN DETERMINISTIC - ODE solver   ####
########################################

# Setting tolerance values
parameters$nifinal <- 10^-9
parameters$tolerance <- 10^-6

# Initial state variables including separate costs
initial_state <- c(Ns = parameters$ns0, 
                   Ni = parameters$ni0, 
                   Nr = parameters$nr0, 
                   Nd = parameters$nd0, 
                   #      lambda_s = 0, 
                   #     lambda_i = 0, 
                   HealthCost = 0, 
                   SocialActivityCost = 0, 
                   TotalCost = 0,
                   a_t = NA, 
                   u_t = NA, 
                   Rt= NA) # add a_t and u_t to keep track of this over time

# Time sequence for pre-shock
times <- seq(0, parameters$time_horizon, by = 1)

# RUN STOCHASTIC BETA REALISATIONS  ####
########################################
# note: make sure the sampled beta is also used for a_t, u_t, Lambda_s and Lambda_i

# get reference: deterministic model
output_sim_deterministic <- run_sir_binomial(initial_state = initial_state, 
                                             times = times, 
                                             parameters = parameters,
                                             bool_stochastic_beta = FALSE,
                                             update_function = get_transitions_deterministic)
# Compute and Track Effective Reproduction Number R(t)



output_experiments <- run_experiments(initial_state = initial_state, 
                                      times = times, 
                                      parameters = parameters, 
                                      bool_stochastic_beta = TRUE, 
                                      update_function = get_transitions_deterministic, 
                                      num_experiments)

# inspect results
compare_sim_output(output_experiments, output_sim_deterministic, plot_tag='stoch beta')


# RUN STOCHASTIC BINIOMIAL MODEL REALISATIONS   ####
####################################################

# get reference: deterministic model
output_sim_deterministic <- run_sir_binomial(initial_state = initial_state, 
                                             times = times, 
                                             parameters = parameters,
                                             update_function = get_transitions_deterministic)

head(output_sim_deterministic)
output_experiments <- run_experiments(initial_state = initial_state, 
                                      times = times, 
                                      parameters = parameters, 
                                      bool_stochastic_beta = FALSE, 
                                      update_function = get_transitions_stochastic, 
                                      num_experiments)

# inspect all results
compare_sim_output(output_experiments, output_sim_deterministic, plot_tag='binomial')

# inspect results excl fadeout
compare_sim_output(output_experiments, output_sim_deterministic, plot_tag='binomial',
                   fadeout_threshold = fadeout_threshold)


# RUN MODEL CHECK   ####
####################################################

# load reference input and output
parameters <- readRDS('test/parameters.rds')
parameters$scenario <- "laissez-faire"
parameters$infect_thres = 0
initial_state <- readRDS('test/initial_state.rds')
num_experiments <- 4
output_reference <- readRDS('test/output_test.rds')

# run stochastic experiments
output_test <- run_experiments(initial_state = initial_state, 
                               times = times, 
                               parameters = parameters, 
                               bool_stochastic_beta = FALSE, 
                               update_function = get_transitions_stochastic, 
                               num_experiments)




# Update reference test files to match latest model structure
update_reference_values <- function(){
  saveRDS(parameters, 'test/parameters.rds')
  saveRDS(initial_state, 'test/initial_state.rds')
  saveRDS(output_test, 'test/output_test.rds')
}

update_reference_values()  # Run it immediately


# compare each element
bool_output_unchanged <- TRUE
for(i_out in 1:length(output_test)){
  output_name <- names(output_test)[i_out]
  if(any(dim(output_test[[i_out]]) != dim(output_reference[[i_out]])) |
     any(output_test[[i_out]] != output_reference[[i_out]],na.rm=TRUE)){
    # print warning
    warning(paste0('Model output "',output_name,'" changed'))
    bool_output_unchanged <- FALSE
  } 
}
# print statement if output did not change
if(bool_output_unchanged){
  print('MODEL OUTPUT DID NOT CHANGE')
}

# Rebuild .rds reference files to match current model version
# ------------------------------------------------------------

# Define parameters (current version)
parameters <- list(
  gamma = 1/7,
  beta = 3/10 + 1/7,
  delta = 0.67 / 365,
  rho = 0.05 / 365,
  pi = 0.0062,
  sigma = 0.10,
  kappa = 197,
  ni0 = 0.000527,
  ns0 = 1 - 0.000527,
  nr0 = 0,
  nd0 = 0,
  v = 31755,
  alpha = 1,
  fx = 123,
  time_horizon = 500,
  utility_type = "Log",
  rng_seed = 123,
  scenario = "laissez-faire",
  R0 = (3/10+1/7)/(1/7),
  infect_thres = 0,
  pop_size = 1e4
)

# Define initial state (without lambda variables)
initial_state <- c(
  Ns = parameters$ns0,
  Ni = parameters$ni0,
  Nr = parameters$nr0,
  Nd = parameters$nd0,
  HealthCost = 0,
  SocialActivityCost = 0,
  TotalCost = 0,
  a_t = NA,
  u_t = NA,
  Rt = NA
)

#test

# Time sequence
times <- seq(0, parameters$time_horizon, by = 1)

# Run reference simulation using updated setup
output_test <- run_experiments(initial_state = initial_state,
                               times = times,
                               parameters = parameters,
                               bool_stochastic_beta = FALSE,
                               update_function = get_transitions_stochastic,
                               num_experiments = 4)

# Save the new .rds files (overwrites old ones)
saveRDS(parameters, 'test/parameters.rds')
saveRDS(initial_state, 'test/initial_state.rds')
saveRDS(output_test, 'test/output_test.rds')


# function to update test input/output
update_reference_values <- function(){
  saveRDS(parameters,'test/parameters.rds')
  saveRDS(initial_state,'test/initial_state.rds')
  saveRDS(output_test,'test/output_test.rds')
}
update_reference_values()