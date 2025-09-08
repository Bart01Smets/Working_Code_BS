#
# Workbench script for the EPI-ECON modelling work
#
# - Deterministic and stochastic modelling
######################################################### #
getwd()


setwd("C:/Users/Bart Smets/OneDrive/Documenten/GitHub/Working_Code_BS")

# clear workbench
rm(list=ls())
library(ggplot2)
# load functions
source("R/epi_econ_lib_new.R")

# SETUP   ####
####################

# Parameters
parameters <- list(
  gamma = 1/7,             # Recovery rate (γ) 1/7, alternatively 1/4
  beta = 3/10+1/7,       # Transmission rate (β) 3/10
  rho = 0.05 / 365,        # Discounting rate (ρ)
  pi = 0.0062,             # Infection fatality rate (π)
  kappa = 197,             # Expected cost of infection (κ)
  ni0 = 0.0005,         # Initial infected population (Ni0)
  ns0 = 0.9995,     # Initial susceptible population (Ns0)
  nr0 = 0,                 # Initial recovered population (Nr0)
  nd0 = 0, # Initial dead population (Nd0)
  v = 31755,
  fx = 123,                # Exchange rate multiplier for USD conversion
  time_horizon = 1500,      # Time of shock
  rng_seed = 150,
  R0= (3/10+1/7)/(1/7),
  pop_size = 1e4,
  infect_thres = 1
)

parameters$beta <- 0.3 + parameters$gamma

# define number of stochastic runs
num_experiments <- 300

# define fadeout threshold
fadeout_threshold = 100

# define default parameters (for development)
plot_tag <- 'dev'
update_function <- get_transitions_stochastic

# RUN DETERMINISTIC - ODE solver   ####
########################################
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

# RUN STOCHASTIC BINIOMIAL MODEL REALISATIONS   ####
####################################################

# get reference: deterministic model
output_sim_deterministic <- run_sir_binomial(initial_state = initial_state, 
                                             times = times, 
                                             parameters = parameters,
                                             update_function = get_transitions_deterministic)
# Print deterministic peak infection time
det_peak_time <- which.max(output_sim_deterministic$Ni) - 1
cat(sprintf("Deterministic Peak Infection Time: %d\n", det_peak_time))

head(output_sim_deterministic)
output_experiments <- run_experiments(initial_state = initial_state, 
                                      times = times, 
                                      parameters = parameters,
                                      update_function = get_transitions_stochastic, 
                                      num_experiments)

# inspect all results
compare_sim_output(output_experiments, output_sim_deterministic, plot_tag='binomial')

# inspect results excl fadeout
compare_sim_output(output_experiments, output_sim_deterministic, plot_tag='binomial',
                   fadeout_threshold = fadeout_threshold)


# ==============================
# MULTI-PARAMETER SENSITIVITY ANALYSIS
# ==============================

baseline_params <- parameters

# convenient helpers for baseline values
baseline_beta  <- baseline_params$beta          # ~0.442857
baseline_gamma <- baseline_params$gamma         # ~0.142857
#-
# labels for plots
param_labels <- list(
  time_horizon = "time horizon",
  v            = "value of life (v)",
  beta         = expression(beta),
  pi           = expression(pi),
  gamma        = expression(gamma),
  ni0          = expression(N[i0]),
  pop_size     = "population size",
  rho          = expression(rho)
)

make_labels <- function(key) {
  lab <- param_labels[[key]]
  if (is.null(lab)) {
    fallback <- gsub("_", " ", key)
    return(list(
      title = paste("Sensitivity of total cost to", fallback),
      x     = paste("Value of", fallback)
    ))
  }
  if (is.expression(lab)) {
    sym <- lab[[1]]
    return(list(
      title = bquote("Sensitivity of total cost to " ~ .(sym)),
      x     = bquote("Value of " ~ .(sym))
    ))
  } else {
    return(list(
      title = paste("Sensitivity of total cost to", lab),
      x     = paste("Value of", lab)
    ))
  }
}

# ranges: β centered on baseline ±0.3; γ from 1/14 .. 1/3
beta_halfwidth <- (0.7 - 0.1) / 2  # = 0.3
param_ranges <- list(
  time_horizon = seq(200, 2000, length.out = 10),
  v            = seq(10000, 50000, length.out = 10),
  beta         = seq(baseline_beta - beta_halfwidth,
                     baseline_beta + beta_halfwidth,
                     length.out = 10),
  pi           = seq(0.001, 0.01, length.out = 10),
  gamma        = seq(1/14, 1/3, length.out = 10),
  ni0          = seq(0.0005, 0.05, length.out = 10),
  pop_size     = round(seq(10000, 100000, length.out = 10)),
  rho          = seq(0, 0.2/365, length.out = 10)
)

pdf("sensitivity_all.pdf", width = 8, height = 6)

for (sensitivity_target in names(param_ranges)) {
  cat("Running sensitivity for:", sensitivity_target, "\n")
  
  # reset to baseline each target
  parameters <- baseline_params
  
  # build grid for the current target
  param_grid <- data.frame(value = param_ranges[[sensitivity_target]])
  colnames(param_grid) <- sensitivity_target
  
  # dashed line at baseline of the target parameter
  baseline_value <- if (sensitivity_target == "beta") {
    baseline_params$beta
  } else {
    baseline_params[[sensitivity_target]]
  }
  
  sensitivity_results <- data.frame()
  
  for (i in 1:nrow(param_grid)) {
    cat("  Scenario", i, "of", nrow(param_grid), "\n")
    parameters <- baseline_params  # fresh copy each scenario
    
    if (sensitivity_target == "beta") {
      # Vary β only; γ stays at baseline
      parameters$beta <- param_grid[[sensitivity_target]][i]
      parameters$R0   <- parameters$beta / parameters$gamma
      val <- parameters$beta
      
    } else if (sensitivity_target == "gamma") {
      # Vary γ only; β stays at baseline
      parameters$gamma <- param_grid[[sensitivity_target]][i]
      parameters$R0    <- parameters$beta / parameters$gamma
      val <- parameters$gamma
      
    } else {
      # Other parameters vary normally (β and γ untouched)
      parameters[[sensitivity_target]] <- param_grid[[sensitivity_target]][i]
      # keep R0 consistent with whatever β,γ currently are
      parameters$R0 <- parameters$beta / parameters$gamma
      val <- param_grid[[sensitivity_target]][i]
    }
    
    # Initial state and runs
    initial_state <- c(
      Ns = parameters$ns0, Ni = parameters$ni0, Nr = parameters$nr0, Nd = parameters$nd0,
      HealthCost = 0, SocialActivityCost = 0, TotalCost = 0,
      a_t = NA, u_t = NA, Rt = NA
    )
    times <- seq(0, parameters$time_horizon, by = 1)
    
    output_det <- run_sir_binomial(
      initial_state, times, parameters,
      update_function = get_transitions_deterministic
    )
    total_cost_det <- output_det[nrow(output_det), "TotalCost"]
    
    output_stoch <- run_experiments(
      initial_state, times, parameters,
      get_transitions_stochastic, num_experiments
    )
    total_costs_stoch <- sapply(
      1:num_experiments,
      function(j) output_stoch$output_all[j, length(times), "TotalCost"]
    )
    
    mean_cost_stoch <- mean(total_costs_stoch)
    q_lo <- quantile(total_costs_stoch, 0.025); q_hi <- quantile(total_costs_stoch, 0.975)
    se <- sd(total_costs_stoch) / sqrt(num_experiments)
    ci_lower <- mean_cost_stoch - 1.96 * se; ci_upper <- mean_cost_stoch + 1.96 * se
    
    sensitivity_results <- rbind(
      sensitivity_results,
      data.frame(Parameter = sensitivity_target, Value = val, Type = "Deterministic",
                 TotalCost = total_cost_det, Q025 = NA, Q975 = NA, CI025 = NA, CI975 = NA),
      data.frame(Parameter = sensitivity_target, Value = val, Type = "Stochastic",
                 TotalCost = mean_cost_stoch, Q025 = q_lo, Q975 = q_hi, CI025 = ci_lower, CI975 = ci_upper)
    )
  }
  
  sensitivity_results <- sensitivity_results[order(sensitivity_results$Type, sensitivity_results$Value), ]
  
  labs_ <- make_labels(sensitivity_target)
  out_csv <- paste0("sensitivity_", sensitivity_target, ".csv")
  write.csv(sensitivity_results, out_csv, row.names = FALSE)
  
  p <- ggplot(sensitivity_results, aes(x = Value, y = TotalCost, color = Type, group = Type)) +
    geom_ribbon(data = subset(sensitivity_results, Type == "Stochastic"),
                aes(ymin = Q025, ymax = Q975, fill = "95% Credible Interval"),
                alpha = 0.2, inherit.aes = TRUE) +
    geom_ribbon(data = subset(sensitivity_results, Type == "Stochastic"),
                aes(ymin = CI025, ymax = CI975, fill = "95% Confidence Interval"),
                alpha = 0.3, inherit.aes = TRUE) +
    geom_line(data = subset(sensitivity_results, Type == "Stochastic"),
              aes(y = CI025), linetype = "solid", size = 0.7) +
    geom_line(data = subset(sensitivity_results, Type == "Stochastic"),
              aes(y = CI975), linetype = "solid", size = 0.7) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    geom_vline(xintercept = baseline_value, linetype = "dashed", size = 1, alpha = 0.8) +
    labs(title = labs_$title, x = labs_$x, y = "Total cost (per capita)") +
    theme_minimal() +
    scale_color_manual(values = c("Deterministic" = "black", "Stochastic" = "steelblue"), guide = "none") +
    scale_fill_manual(values = c("95% Credible Interval" = "steelblue",
                                 "95% Confidence Interval" = "red"),
                      guide = "none") +
    theme(legend.position = "none")
  
  print(p)
}

dev.off()

