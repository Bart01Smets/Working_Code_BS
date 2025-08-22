#
# Workbench script for the EPI-ECON modelling work
#
# - Deterministic and stochastic modelling
# - Side-by-side baselines for (A) stochastic β and (B) binomial transitions
# - Dual-uncertainty PSA summaries + bands
#########################################################

# Optional: make paths robust (comment out if you prefer setwd)
# install.packages("here") # once
# setwd(here::here())


getwd()
setwd("C:/Users/Bart Smets/OneDrive/Documenten/GitHub/Working_Code_BS")

rm(list = ls())

source("scenarios/epi_econ_lib_par_uncert2.R")

# =====================
# SETUP
# =====================

parameters <- list(
  gamma = 1/7,                   # Recovery rate (γ)
  beta = 3/10 + 1/7,             # Transmission rate (β)
  rho = 0.72 / 365,              # Discounting rate (ρ)
  pi = 0.0062,                   # Infection fatality rate (π)
  kappa = 197,                   # Expected cost of infection (κ)  [not used in current a()]
  ni0 = 0.0005,                  # Initial infected proportion
  ns0 = 0.9995,                  # Initial susceptible proportion
  nr0 = 0,                       # Initial recovered
  nd0 = 0,                       # Initial dead
  v = 31755,
  fx = 123,                      # USD conversion multiplier
  time_horizon = 1500,
  rng_seed = 150,
  R0 = (3/10 + 1/7) / (1/7),
  pop_size = 1e6,
  infect_thres = 1,
  utility_type = "Log",
  sigma = 0.19,                  # SD for stochastic β
  bool_regular_sird = FALSE,
  bool_daily_cost_minimizing = FALSE
)

# Workbench tweaks
parameters$pop_size <- 1e4
num_experiments     <- 100
fadeout_threshold   <- 100

# Initial state (proportions; kernel converts to counts internally)
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

times <- seq(0, parameters$time_horizon, by = 1)

# =====================
# DETERMINISTIC REFERENCE
# =====================

output_det <- run_sir_binomial(
  initial_state = initial_state,
  times         = times,
  parameters    = parameters,
  bool_stochastic_beta = FALSE,
  update_function      = get_transitions_deterministic
)

det_peak_time <- which.max(output_det$Ni) - 1
cat(sprintf("[Deterministic] Peak Infection Time: %d\n", det_peak_time))

# =====================
# (A) STOCHASTIC β (no binomial; deterministic transitions with random β_t)
# =====================

output_stochbeta <- run_experiments(
  initial_state        = initial_state,
  times                = times,
  parameters           = parameters,
  bool_stochastic_beta = TRUE,
  update_function      = get_transitions_deterministic,
  num_experiments      = num_experiments
)

compare_sim_output(
  output_experiments   = output_stochbeta,
  output_deterministic = output_det,
  plot_tag             = "stoch_beta (all)"
)

compare_sim_output(
  output_experiments   = output_stochbeta,
  output_deterministic = output_det,
  plot_tag             = "stoch_beta (excl. fadeout)",
  fadeout_threshold    = fadeout_threshold
)

# =====================
# (B) BINOMIAL STOCHASTICITY (fixed β; binomial transitions)
# =====================

output_binom <- run_experiments(
  initial_state        = initial_state,
  times                = times,
  parameters           = parameters,
  bool_stochastic_beta = FALSE,
  update_function      = get_transitions_stochastic,
  num_experiments      = num_experiments
)

compare_sim_output(
  output_experiments   = output_binom,
  output_deterministic = output_det,
  plot_tag             = "binomial (all)"
)

compare_sim_output(
  output_experiments   = output_binom,
  output_deterministic = output_det,
  plot_tag             = "binomial (excl. fadeout)",
  fadeout_threshold    = fadeout_threshold
)

# =====================
# DUAL-UNCERTAINTY PSA (parameter uncertainty + stochastic)
# We generate TWO reports: one using the β-stochastic baseline, one using the binomial baseline
# =====================

if (!dir.exists("figures")) dir.create("figures")
library(ggplot2)

# --- PSA engine outputs
psa_output_all <- run_psa_full(
  initial_state        = initial_state,
  times                = times,
  base_parameters      = parameters,
  n_param_draws        = 10,     # increase later
  n_stoch              = 5,
  bool_stochastic_beta = TRUE,   # TRUE = include β-stochasticity in PSA paths
  update_function      = get_transitions_stochastic
)

psa <- run_psa(
  initial_state        = initial_state,
  times                = times,
  base_parameters      = parameters,
  n_param_draws        = 30,     # increase later
  n_stoch              = 10,
  bool_stochastic_beta = TRUE,
  update_function      = get_transitions_stochastic
)

# Save numeric PSA summaries
write.csv(psa$psa_means, "figures/psa_means.csv", row.names = FALSE)
write.csv(psa$psa_bands, "figures/psa_bands.csv", row.names = FALSE)

# =====================
# PSA REPORT — baseline = stochastic β
# =====================

pdf("figures/uncertainty_report_stochbeta.pdf", width = 8, height = 5)

# Page 1: Total cost vs R0 (per-draw mean)
print(
  ggplot(psa$psa_means, aes(x = R0, y = TotalCost)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = "PSA: Total Cost vs R0 (per-draw means)",
         x = "Sampled R0", y = "Total Cost (per capita)") +
    theme_minimal()
)

# Pages 2–5: Dual-uncertainty bands with stochastic-β baseline
for (vv in c("TotalCost","HealthCost","SocialActivityCost","Ni")) {
  ttl <- switch(vv,
                "TotalCost" = "Dual Uncertainty Bands — Total Cost",
                "HealthCost" = "Dual Uncertainty Bands — Health Cost",
                "SocialActivityCost" = "Dual Uncertainty Bands — Social Activity Cost",
                "Ni" = "Dual Uncertainty Bands — Infectives (Ni)",
                vv)
  print(
    plot_dual_uncertainty_bands(
      baseline_output_all = output_stochbeta$output_all,
      psa_output_all      = psa_output_all,
      deterministic_df    = output_det,
      var                 = vv,
      title               = ttl
    )
  )
}
dev.off()

# =====================
# PSA REPORT — baseline = binomial
# (same PSA pool; only the baseline band changes)
# =====================

pdf("figures/uncertainty_report_binomial.pdf", width = 8, height = 5)

# Page 1: Total cost vs R0
print(
  ggplot(psa$psa_means, aes(x = R0, y = TotalCost)) +
    geom_point(alpha = 0.6) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(title = "PSA: Total Cost vs R0 (per-draw means)",
         x = "Sampled R0", y = "Total Cost (per capita)") +
    theme_minimal()
)

# Pages 2–5: Dual-uncertainty bands with binomial baseline
for (vv in c("TotalCost","HealthCost","SocialActivityCost","Ni")) {
  ttl <- switch(vv,
                "TotalCost" = "Dual Uncertainty Bands — Total Cost",
                "HealthCost" = "Dual Uncertainty Bands — Health Cost",
                "SocialActivityCost" = "Dual Uncertainty Bands — Social Activity Cost",
                "Ni" = "Dual Uncertainty Bands — Infectives (Ni)",
                vv)
  print(
    plot_dual_uncertainty_bands(
      baseline_output_all = output_binom$output_all,
      psa_output_all      = psa_output_all,
      deterministic_df    = output_det,
      var                 = vv,
      title               = ttl
    )
  )
}
dev.off()

# =====================
# OPTIONAL: MULTI-PARAMETER SENSITIVITY (unchanged; keep commented until needed)
# =====================

# library(ggplot2)
# param_ranges <- list(
#   v = seq(10000, 50000, length.out = 10),
#   beta = seq(0.1, 0.7, length.out = 10),
#   pi = seq(0.001, 0.01, length.out = 10),
#   gamma = seq(1/14, 1/3, length.out = 10),
#   ni0 = seq(0.0001, 0.05, length.out = 10),
#   pop_size = round(seq(10000, 100000, length.out = 10)),
#   time_horizon = seq(200, 2000, length.out = 10),
#   rho = seq(0, 100/365, length.out = 10)
# )
# pdf("sensitivity_all.pdf", width = 8, height = 6)
# ... (your loop block unchanged) ...
# dev.off()
