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
source("R/epi_econ_lib_new.R")

# SETUP   ####
####################

# Parameters
parameters <- list(
  gamma = 1/7,             # Recovery rate (γ) 1/7, alternatively 1/3
  beta = 3/10+1/7,       # Transmission rate (β) 3/10
  delta = 0.67 / 365,      # Arrival rate of cure (δ)
  rho = 0.05 / 365,        # Discounting rate (ρ)
  pi = 0.0062,             # Infection fatality rate (π)
  sigma = 0.10,            # Determines stochastic beta_t =  N(beta,sigma)
  kappa = 197,             # Expected cost of infection (κ)
  ni0 = 0.000527,         # Initial infected population (Ni0)297
  ns0 = 1 - 0.000527,     # Initial susceptible population (Ns0)527
  nr0 = 0.000024840,                 # Initial recovered population (Nr0)
  nd0 = 0.000000156, # Initial dead population (Nd0)
  v = 31755,
  alpha = 0,               # Altruism parameter (default set to no altruism)
  fx = 123,                # Exchange rate multiplier for USD conversion
  time_horizon = 1500,      # Time of shock
  utility_type = "Log",    # Utility type: "Log" or "Quadratic"
  #rng_seed = 123,
  R0= (3/10+1/7)/(1/7),
  pop_size = 1e6,
  infect_thres = 1
)


# define population size
parameters$pop_size <- 10000#1e4

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
