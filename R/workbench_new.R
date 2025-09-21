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
  ni0 = 0.0005,         # Initial infected population (Ni0)
  ns0 = 0.9995,     # Initial susceptible population (Ns0)
  nr0 = 0,                 # Initial recovered population (Nr0)
  nd0 = 0, # Initial dead population (Nd0)
  v = 31755,
  fx = 123,                # Exchange rate multiplier for USD conversion
  time_horizon = 1500,      # Time of shock
  rng_seed = 150,
  kappa= 197,
  R0= (3/10+1/7)/(1/7),
  pop_size = 1e4,
  infect_thres = 1,
  bool_regular_sird = FALSE  # TRUE => regular SIRD (a=1, no SocialActivityCost)
)

parameters$beta <- 0.3 + parameters$gamma

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
#
# # ==============================
# # FIGURE: Δ Total Cost (Stochastic − Deterministic) vs Number of Simulations
# # ==============================
# 
# # Grid of simulation counts
# n_grid <- seq(100, 1000, length.out = 10)
# fadeout_threshold_for_diff <- 100  # set 0 to include all; uses the same logic as elsewhere
# 
# # Initial state and time grid (use your current time_horizon)
# init_state <- c(Ns = parameters$ns0, Ni = parameters$ni0, Nr = parameters$nr0, Nd = parameters$nd0,
#                 HealthCost = 0, SocialActivityCost = 0, TotalCost = 0,
#                 a_t = NA, u_t = NA, Rt = NA)
# times <- 0:parameters$time_horizon
# 
# # Deterministic baseline (run once)
# det_out <- run_sir_binomial(
#   initial_state = init_state, times = times, parameters = parameters,
#   update_function = get_transitions_deterministic
# )
# det_total <- det_out[nrow(det_out), "TotalCost"]
# 
# # Storage
# diff_df <- data.frame(
#   num_sims   = integer(),
#   det_total  = numeric(),
#   stoch_mean = numeric(),
#   stoch_ci_lo = numeric(),
#   stoch_ci_hi = numeric(),
#   diff_mean  = numeric(),   # (stoch_mean - det_total)
#   diff_ci_lo = numeric(),   # (stoch_ci_lo - det_total)
#   diff_ci_hi = numeric()    # (stoch_ci_hi - det_total)
# )
# 
# for (N in n_grid) {
#   # Run N stochastic experiments
#   sims <- run_experiments(
#     initial_state = init_state, times = times, parameters = parameters,
#     update_function = get_transitions_stochastic,
#     num_experiments = N
#   )
# 
#   # Optionally exclude fadeouts
#   summ <- sims$output_summary
#   if (fadeout_threshold_for_diff > 0) {
#     keep <- (summ$Nr + summ$Ni) >= fadeout_threshold_for_diff
#     if (!any(keep)) {
#       message(sprintf("All %d sims filtered out at N=%d; skipping this point.", nrow(summ), N))
#       next
#     }
#     summ <- summ[keep, , drop = FALSE]
#   }
# 
#   # 95% CI for the mean stochastic total cost (uses {confintr})
#   stoch_ci <- ci_mean(summ$TotalCost)
#   stoch_mean <- as.numeric(stoch_ci$estimate)
#   stoch_lo   <- as.numeric(stoch_ci$interval[1])
#   stoch_hi   <- as.numeric(stoch_ci$interval[2])
# 
# 
#   # Differences vs deterministic
#   diff_mean <- stoch_mean - det_total
#   diff_lo   <- stoch_lo   - det_total
#   diff_hi   <- stoch_hi   - det_total
# 
#   diff_df <- rbind(diff_df, data.frame(
#     num_sims = N,
#     det_total = det_total,
#     stoch_mean = stoch_mean,
#     stoch_ci_lo = stoch_lo,
#     stoch_ci_hi = stoch_hi,
#     diff_mean = diff_mean,
#     diff_ci_lo = diff_lo,
#     diff_ci_hi = diff_hi
#   ))
# }
# 
# # Ensure output dir
# if (!dir.exists("figures")) dir.create("figures")
# 
# # Plot
# library(ggplot2)
# p_diff_nsims <- ggplot(diff_df, aes(x = num_sims, y = diff_mean)) +
#   geom_ribbon(aes(ymin = diff_ci_lo, ymax = diff_ci_hi), fill="red", alpha = 0.25) +
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
# 
# # Optional: inspect the table
# print(diff_df)


# # ==============================
# # SENSITIVITY: Δ Total Cost (Stochastic − Deterministic)
# # ==============================
# 
# # --- Config ---
# baseline_params <- parameters
# num_experiments_sens <- 500       # tweak as needed
# fadeout_threshold_for_diff <- 100    # set 0 to include all
# if (!dir.exists("figures")) dir.create("figures")
# 
# suppressPackageStartupMessages({
#   library(ggplot2)
#   library(confintr)
# })
# 
# param_ranges <- list(
#   time_horizon = seq(200, 2000, length.out = 10),
#   kappa        = seq(50, 500, length.out = 10),   # instead of pi and v
#   beta         = seq(0.15, 0.65, length.out = 10),
#   gamma        = seq(1/14, 1/3, length.out = 10),
#   Ni0          = round(seq(baseline_params$ni0 * baseline_params$pop_size,
#                            0.05 * baseline_params$pop_size,
#                            length.out = 10)),
#   pop_size     = round(seq(10000, 100000, length.out = 10)),
#   rho          = seq(0, 0.2/365, length.out = 10)
# )
# 
# # --- Math-safe axis labels (plotmath expressions) ---
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
# 
# # --- Titles that correspond to the variable name (matching x-axis) ---
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
# run_sensitivity_diff <- function(param_name, values) {
#   cat("Running Δ-cost sensitivity for:", param_name, "\n")
#   results <- data.frame()
# 
#   for (i in seq_along(values)) {
#     val <- values[i]
#     cat(sprintf("  %s = %s (%d/%d)\n", param_name, format(val, digits = 6), i, length(values)))
# 
#     # Copy + set parameter under test
#     params <- baseline_params
# if (param_name == "Ni0") {
#   params$ni0 <- val / params$pop_size
# } else {
#   params[[param_name]] <- val   # κ is just another primitive here
# }
#     params$R0 <- params$beta / params$gamma
# 
#     # Initial state and horizon
#     init_state <- c(
#       Ns = params$ns0, Ni = params$ni0, Nr = params$nr0, Nd = params$nd0,
#       HealthCost = 0, SocialActivityCost = 0, TotalCost = 0,
#       a_t = NA, u_t = NA, Rt = NA
#     )
#     times <- 0:params$time_horizon
# 
#     # Deterministic run
#     det_out <- run_sir_binomial(
#       initial_state = init_state, times = times, parameters = params,
#       update_function = get_transitions_deterministic
#     )
#     det_total <- det_out[nrow(det_out), "TotalCost"]
# 
#     # Stochastic experiments
#     sims <- run_experiments(
#       initial_state = init_state, times = times, parameters = params,
#       update_function = get_transitions_stochastic,
#       num_experiments = num_experiments_sens
#     )
#     summ <- sims$output_summary
# 
#     # Optional filtering for fadeout (use Nr + Ni as observed size)
#     if (fadeout_threshold_for_diff > 0) {
#       keep <- (summ$Nr + summ$Ni) >= fadeout_threshold_for_diff
#       if (!any(keep)) {
#         warning(sprintf("All sims filtered out for %s=%s; skipping point.", param_name, as.character(val)))
#         next
#       }
#       summ <- summ[keep, , drop = FALSE]
#     }
# 
#     # CI/quantiles on stochastic totals
#     stoch_totals <- summ$TotalCost
#     ci <- ci_mean(stoch_totals)
#     mean_stoch <- as.numeric(ci$estimate)
#     ci_lo <- as.numeric(ci$interval[1]); ci_hi <- as.numeric(ci$interval[2])
#     q <- quantile(stoch_totals, c(0.025, 0.975), na.rm = TRUE)
#     q_lo <- as.numeric(q[1]); q_hi <- as.numeric(q[2])
# 
#     # Differences (stochastic − deterministic)
#     diff_mean <- mean_stoch - det_total
#     diff_ci_lo <- ci_lo - det_total; diff_ci_hi <- ci_hi - det_total
#     diff_q_lo  <- q_lo  - det_total; diff_q_hi  <- q_hi  - det_total
# 
#     results <- rbind(results, data.frame(
#       Parameter = param_name, Value = val,
#       Det_Total = det_total, Stoch_Mean = mean_stoch,
#       CI_Lo = ci_lo, CI_Hi = ci_hi, Q_Lo = q_lo, Q_Hi = q_hi,
#       Diff_Mean = diff_mean, Diff_CI_Lo = diff_ci_lo, Diff_CI_Hi = diff_ci_hi,
#       Diff_Q_Lo  = diff_q_lo,  Diff_Q_Hi  = diff_q_hi
#     ))
#   }
# 
#   # --- Baseline value for vertical line ---
#   baseline_val <- switch(param_name,
#                          "Ni0"      = baseline_params$ni0 * baseline_params$pop_size,  # back to counts
#                          "kappa"    = baseline_params$kappa,  # you set this above as pi*v
#                          "beta"     = baseline_params$beta,
#                          "gamma"    = baseline_params$gamma,
#                          "rho"      = baseline_params$rho,
#                          "pop_size" = baseline_params$pop_size,
#                          "time_horizon" = baseline_params$time_horizon,
#                          baseline_params[[param_name]]
#   )
#   # If the x grid is integer (e.g., Ni0, pop_size), match the scale for aesthetics
#   if (param_name %in% c("Ni0", "pop_size")) baseline_val <- as.integer(round(baseline_val))
# 
#   # Only draw the vline if baseline is within plotting range
#   add_baseline <- (baseline_val >= min(results$Value, na.rm = TRUE)) &&
#     (baseline_val <= max(results$Value, na.rm = TRUE))
# 
# 
#   # Plot with math-safe labels/titles
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
# 
# 
#   print(p)
#   return(results)
# }
# 
# # --- Run for all parameters in param_ranges, save to ONE pdf ---
# out_pdf_all <- "figures/sensitivity_diff_ALL.pdf"
# 
# # If fonts render oddly, use CairoPDF() instead of pdf():
# # Cairo::CairoPDF(out_pdf_all, width = 8, height = 5)
# pdf(out_pdf_all, width = 8, height = 5)
# 
# all_results <- list()
# for (sensitivity_target in names(param_ranges)) {
#   results <- run_sensitivity_diff(sensitivity_target, param_ranges[[sensitivity_target]])
#   all_results[[sensitivity_target]] <- results
# }
# dev.off()
# 
# # Optional: combined CSV with all results
# sensitivity_results_all <- do.call(rbind, all_results)
# write.csv(sensitivity_results_all, "sensitivity_diff_ALL.csv", row.names = FALSE)
# 
# cat("Saved PDF to:", out_pdf_all, "\n")
# cat("Saved CSV to: sensitivity_diff_ALL.csv\n")
# 
