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
  gamma = 1/7,             # Recovery rate (γ)
  beta = 3/10+1/7,       # Transmission rate (β)
  rho = 0.05 / 365,        # Discounting rate (ρ)
  pi = 0.0062,             # Infection fatality probability (π)
  ni0 = 0.0005,         # Initial infected population (ni0)
  ns0 = 0.9995,     # Initial susceptible population (ns0)
  nr0 = 0,                 # Initial recovered population (nr0)
  nd0 = 0,             # Initial dead population (nd0)
  v = 31755,
  fx = 123,                # Exchange rate multiplier for USD conversion
  time_horizon = 1500,
  rng_seed = 150,
  kappa= 197,
  R0= (3/10+1/7)/(1/7),
  pop_size = 1e4,
  infect_thres = 1,
  bool_regular_sird = FALSE,  # TRUE => regular SIRD (a=1, no ActivityCost)
  
  integer_with_carry = FALSE,   # Use discrete transitions in deterministic case option
  costs_from = "realized"      # "realized" => Deterministic integer costs
  )

# define number of stochastic runs
num_experiments <- 500

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

## ==============================
## FIGURES for Sensitivity Analysis (drop-in)
## - Δ Total Cost vs # of simulations
## - Univariate sensitivity: Δ Total Cost vs parameter
## ==============================

if (!dir.exists("figures")) dir.create("figures")

suppressPackageStartupMessages({
  library(ggplot2)
  library(confintr)
})

# ## --------------------------------------------------------------
# ## FIGURE: Δ Total Cost (Stochastic − Deterministic) vs # Sims
# ## --------------------------------------------------------------
# n_grid <- seq(100, 1000, length.out = 10)
# diff_df <- compute_delta_vs_num_sims(parameters, n_grid, fadeout_threshold_for_diff = 100)
# 
# p_diff_nsims <- ggplot(diff_df, aes(x = num_sims, y = diff_mean)) +
#   geom_ribbon(aes(ymin = diff_ci_lo, ymax = diff_ci_hi), fill = "red", alpha = 0.25) +
#   geom_line(size = 1) +
#   geom_point(size = 1.8) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   scale_x_continuous(breaks = n_grid) +
#   labs(
#     title = "Impact of simulations on difference in total cost",
#     x = "Number of stochastic simulations",
#     y = "Cost difference"
#   ) +
#   theme_minimal(base_size = 12)
# 
# print(p_diff_nsims)
# ggsave("figures/diff_totalcost_vs_num_sims.pdf", p_diff_nsims, width = 8, height = 5)
# ggsave("figures/diff_totalcost_vs_num_sims.png", p_diff_nsims, width = 8, height = 5, dpi = 200)
# print(diff_df)

# ## --------------------------------------------------------------
# ## SENSITIVITY: Δ Total Cost (Stochastic − Deterministic)
# ## --------------------------------------------------------------
# 
# # Math-safe axis labels / titles 
# x_label_expr <- function(param_name) {
#   switch(param_name,
#          "beta"         = expression(paste("Value of ", beta, " (transmission rate)")),
#          "gamma"        = expression(paste("Value of ", gamma, " (recovery rate)")),
#          "rho"          = expression(paste("Value of ", rho, " (discount rate)")),
#          "kappa"        = expression(paste("Value of ", kappa, " (infection cost parameter)")),
#          "Ni0"          = expression(paste("Value of ", N[i](0), " (initial infected)")),
#          "pop_size"     = expression(paste("Value of ", N, " (population size)")),
#          "time_horizon" = expression(paste("Value of ", T[max], " (time horizon)")),
#          as.expression(bquote(.(paste0("Value of ", gsub("_", " ", param_name)))))
#   )
# }
# title_expr <- function(param_name) {
#   switch(param_name,
#          "beta"         = expression(paste("Total cost difference for different values of ", beta)),
#          "gamma"        = expression(paste("Total cost difference for different values of ", gamma)),
#          "rho"          = expression(paste("Total cost difference for different values of ", rho)),
#          "kappa"        = expression(paste("Total cost difference for different values of ", kappa)),
#          "Ni0"          = expression(paste("Total cost difference for different values of ", N[i](0))),
#          "pop_size"     = expression(paste("Total cost difference for different values of ", N)),
#          "time_horizon" = expression(paste("Total cost difference for different values of ", T[max])),
#          as.expression(bquote(.(paste("Total cost difference for different values of", param_name))))
#   )
# }
# 
# param_ranges <- list(
#   time_horizon = seq(200, 2000, length.out = 10),
#   kappa        = seq(50, 500, length.out = 10),
#   beta         = seq(0.15, 0.65, length.out = 10),
#   gamma        = seq(1/14, 1/3, length.out = 10),
#   Ni0          = round(seq(parameters$ni0 * parameters$pop_size,
#                            0.05 * parameters$pop_size,
#                            length.out = 10)),
#   pop_size     = round(seq(10000, 100000, length.out = 10)),
#   rho          = seq(0, 0.2/365, length.out = 10)
# )
# 
# out_pdf_all <- "figures/sensitivity_diff_ALL.pdf"
# pdf(out_pdf_all, width = 8, height = 5)
# 
# all_results <- list()
# for (param_name in names(param_ranges)) {
#   results <- compute_sensitivity_grid(parameters, param_name, param_ranges[[param_name]],
#                                       num_experiments_sens = 500,
#                                       fadeout_threshold_for_diff = 100)
#   
#   # Baseline vline (same logic as manual block)
#   baseline_val <- switch(param_name,
#                          "Ni0"          = parameters$ni0 * parameters$pop_size,
#                          "kappa"        = parameters$kappa,
#                          "beta"         = parameters$beta,
#                          "gamma"        = parameters$gamma,
#                          "rho"          = parameters$rho,
#                          "pop_size"     = parameters$pop_size,
#                          "time_horizon" = parameters$time_horizon,
#                          parameters[[param_name]])
#   if (param_name %in% c("Ni0", "pop_size")) baseline_val <- as.integer(round(baseline_val))
#   add_baseline <- (baseline_val >= min(results$Value, na.rm = TRUE)) &&
#     (baseline_val <= max(results$Value, na.rm = TRUE))
#   
#   p <- ggplot(results, aes(x = Value, y = Diff_Mean)) +
#     geom_ribbon(aes(ymin = Diff_Q_Lo, ymax = Diff_Q_Hi,
#                     fill = "95% Quantile Interval"), alpha = 0.20) +
#     geom_ribbon(aes(ymin = Diff_CI_Lo, ymax = Diff_CI_Hi,
#                     fill = "95% Confidence Interval"), alpha = 0.30) +
#     geom_line(size = 1, color = "black") +
#     geom_point(size = 2, color = "black") +
#     geom_hline(yintercept = 0, linetype = "dashed") +
#     { if (add_baseline) geom_vline(xintercept = baseline_val,
#                                    linetype = "dotted", color = "blue") else NULL } +
#     labs(
#       title = title_expr(param_name),
#       x     = x_label_expr(param_name),
#       y     = "Cost difference"
#     ) +
#     theme_minimal(base_size = 12) +
#     scale_fill_manual(values = c("95% Quantile Interval" = "steelblue",
#                                  "95% Confidence Interval" = "red"),
#                       guide = "none")
#   
#   print(p)
#   all_results[[param_name]] <- results
# }
# dev.off()
# 
# sensitivity_results_all <- do.call(rbind, all_results)
# write.csv(sensitivity_results_all, "sensitivity_diff_ALL.csv", row.names = FALSE)
# cat("Saved PDF to:", out_pdf_all, "\n")
# cat("Saved CSV to: sensitivity_diff_ALL.csv\n")



