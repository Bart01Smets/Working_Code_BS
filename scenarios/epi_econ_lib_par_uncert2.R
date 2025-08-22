##
# Helper functions for the EPI-ECON modelling work
#
# - Deterministic and stochastic transitions
################################################################
library(confintr)
library(scales)

# ---------------------------
# Activity rule (with β override)
# ---------------------------
a_function <- function(Ni, Ns, parameters, beta_override = NULL) {
  if (isTRUE(parameters$bool_regular_sird)) return(1)  # Always full activity
  
  # Convert to proportions
  Ni_prop <- Ni / parameters$pop_size
  Ns_prop <- Ns / parameters$pop_size
  
  # Defensive early exit if Ni or Ns is ~0
  if (Ni_prop < 1e-10 || Ns_prop < 1e-10) return(1)
  
  # Use step-specific beta if provided
  beta_use <- if (is.null(beta_override)) parameters$beta else beta_override
  
  # Quadratic terms (your chosen myopic/planner rule)
  multiplier <- beta_use * parameters$pi * parameters$v * Ns_prop * Ni_prop
  num_sum    <- Ns_prop + Ni_prop
  discrim    <- num_sum^2 + 4 * multiplier * num_sum
  denom      <- 2 * multiplier
  
  if (denom <= 0 || discrim < 0) return(1)
  
  a_t <- (-num_sum + sqrt(discrim)) / denom
  max(0, min(1, a_t))  # bound to [0,1]
}

# (Optional) daily cost-minimizing variant (left commented)
# a_function_cost_minimizing <- function(Ns, Ni, t, parameters) { ... }

# ---------------------------
# Utility
# ---------------------------
utility_function <- function(a_t, parameters) {
  if (isTRUE(parameters$bool_regular_sird)) return(0)
  # Log utility by default
  log(a_t) - a_t + 1
  # For quadratic: return(-0.5 * (1 - a_t)^2)
}

# ---------------------------
# Effective reproduction number
# ---------------------------
calculate_Rt <- function(R0, a_t, Ns_prop, Ni) {
  if (Ni <= 0) return(0)
  R0 * a_t^2 * Ns_prop
}

# ---------------------------
# Transition helpers
# ---------------------------
get_transitions_stochastic <- function(n, prob) {
  if (is.na(n) || is.na(prob) || !is.numeric(n) || !is.numeric(prob) ||
      n < 0 || prob < 0 || prob > 1) {
    warning("Invalid rbinom() inputs: returning 0")
    return(0)
  }
  sum(rbinom(n = n, size = 1, prob = prob))
}

get_transitions_deterministic <- function(n, prob){
  if (length(prob) == n) sum(prob) else n * prob
}

# ---------------------------
# Core SIRD loop (binomial kernel)
# ---------------------------
run_sir_binomial <- function(initial_state,
                             times,
                             parameters,
                             bool_stochastic_beta = FALSE,
                             update_function = get_transitions_stochastic){
  
  # Copy & convert initial states (proportions -> counts)
  states <- data.frame(t(initial_state))
  states[grepl('N.', names(states))] <-
    round(states[grepl('N.', names(states))] * parameters$pop_size)
  
  # Fix total pop by adjusting Ns if needed
  states[grepl('Ns', names(states))] <-
    parameters$pop_size - sum(states[grepl('N.', names(states)) & !grepl('Ns', names(states))])
  
  # Output matrix
  states_out <- matrix(NA_real_, nrow = length(times), ncol = length(states))
  colnames(states_out) <- names(states)
  
  # Unpack initial states
  Ns <- states$Ns; Ni <- states$Ni; Nr <- states$Nr; Nd <- states$Nd
  HealthCost <- 0;   SocialActivityCost <- 0
  fx_per_capita <- parameters$fx / parameters$pop_size
  
  # Write initial row (time 0)
  states_out[1, ] <- c(Ns, Ni, Nr, Nd,
                       HealthCost, SocialActivityCost,
                       HealthCost + SocialActivityCost,
                       a_t = NA, u_t = NA, Rt = NA)
  
  # Time loop (skip t=0 which is written)
  for (i_day in times[-1]) {
    
    # Possibly stochastic β_t
    beta_t <- if (bool_stochastic_beta) {
      max(0.0001, rnorm(1, mean = parameters$beta, sd = parameters$sigma))
    } else {
      parameters$beta
    }
    
    # Activity (pass β_t)
    a_t <- if (isTRUE(parameters$bool_daily_cost_minimizing)) {
      a_function_cost_minimizing(Ns, Ni, i_day, parameters)
    } else {
      a_function(Ni, Ns, parameters, beta_override = beta_t)
    }
    
    # Utility (pass parameters explicitly)
    u_t <- utility_function(a_t, parameters)
    
    # Capacity-adjusted IFR
    pi_effective <- if (!is.null(parameters$healthcare_capacity) &&
                        Ni > parameters$healthcare_capacity) {
      parameters$pi * parameters$excess_mortality_multiplier
    } else {
      parameters$pi
    }
    
    # Transition probabilities
    p_infect  <- 1 - exp(- beta_t * a_t^2 * Ni / parameters$pop_size)
    p_recover <- 1 - exp(- parameters$gamma)
    p_death   <- 1 - exp(- pi_effective)
    
    # Draw/compute transitions
    new_infections <- update_function(Ns, prob = p_infect)
    new_recoveries <- update_function(Ni, prob = p_recover)
    
    # Fadeout rule may fully clear infections
    if ((Ni - new_recoveries) < parameters$infect_thres) {
      new_recoveries <- Ni
      new_infections <- 0
    }
    
    # Deaths must be based on the *final* number of recoveries
    new_death <- update_function(new_recoveries, prob = p_death)
    
    # State deltas
    dNs <- -new_infections
    dNi <-  new_infections - new_recoveries
    dNr <-  new_recoveries - new_death
    dNd <-  new_death
    
    # Costs (per capita, discounted)
    if (!isTRUE(parameters$bool_regular_sird)) {
      disc <- exp(-parameters$rho * i_day)
      HealthCost <- HealthCost +
        fx_per_capita * disc * parameters$v * new_death
      SocialActivityCost <- SocialActivityCost +
        fx_per_capita * disc * (Ns + Ni) * abs(u_t)
    }
    
    # Rt at current step (using pre-update Ns)
    Rt <- calculate_Rt(parameters$R0, a_t, Ns / parameters$pop_size, Ni)
    
    # Apply transitions
    Ns <- Ns + dNs
    Ni <- Ni + dNi
    Nr <- Nr + dNr
    Nd <- Nd + dNd
    
    # Record row (time index = i_day + 1)
    states_out[i_day + 1, ] <- c(Ns, Ni, Nr, Nd,
                                 HealthCost, SocialActivityCost,
                                 HealthCost + SocialActivityCost,
                                 a_t, u_t, Rt)
  }
  
  data.frame(states_out)
}

# ---------------------------
# Parameter sampling (no I/O)
# ---------------------------
sample_parameters <- function(base) {
  draw <- base
  draw$beta  <- max(1e-4, rnorm(1, mean = base$beta,  sd = 0.05))
  draw$gamma <- max(1e-4, rnorm(1, mean = base$gamma, sd = 0.02))
  draw$pi    <- pmax(1e-6, rnorm(1, mean = base$pi,   sd = 0.002))
  draw$v     <- pmax(1e3,  rlnorm(1, meanlog = log(base$v), sdlog = 0.3))
  draw$rho   <- pmax(0,    rnorm(1, mean = base$rho,  sd = 0.2/365))
  draw$ni0   <- pmax(1e-6, rnorm(1, mean = base$ni0,  sd = base$ni0 * 0.30))
  total_fixed <- draw$ni0 + draw$nr0 + draw$nd0
  draw$ns0 <- max(1e-9, 1 - total_fixed)
  draw$R0 <- draw$beta / draw$gamma
  draw
}

# ---------------------------
# PSA (summary tables)
# ---------------------------
run_psa <- function(initial_state, times, base_parameters,
                    n_param_draws = 200, n_stoch = 20,
                    bool_stochastic_beta = TRUE,
                    update_function = get_transitions_stochastic) {
  
  set.seed(base_parameters$rng_seed)
  
  keep_params <- c("beta","gamma","pi","v","rho","ni0","ns0","R0")
  metric_cols <- c("HealthCost","SocialActivityCost","TotalCost","PeakTime")
  
  out_params  <- data.frame(matrix(NA, nrow = n_param_draws, ncol = length(keep_params)))
  names(out_params) <- keep_params
  
  psa_means <- data.frame(matrix(NA, nrow = n_param_draws, ncol = length(metric_cols)))
  names(psa_means) <- metric_cols
  
  band_cols <- c("Health_mean","Health_q025","Health_q975",
                 "Social_mean","Social_q025","Social_q975",
                 "Total_mean","Total_q025","Total_q975",
                 "Peak_mean","Peak_q025","Peak_q975")
  psa_bands <- data.frame(matrix(NA, nrow = n_param_draws, ncol = length(band_cols)))
  names(psa_bands) <- band_cols
  
  for (d in 1:n_param_draws) {
    pars_d <- sample_parameters(base_parameters)
    init_d <- c(Ns = pars_d$ns0, Ni = pars_d$ni0, Nr = pars_d$nr0, Nd = pars_d$nd0,
                HealthCost = 0, SocialActivityCost = 0, TotalCost = 0,
                a_t = NA, u_t = NA, Rt = NA)
    pars_d$rng_seed <- base_parameters$rng_seed + d
    
    sims <- run_experiments(initial_state = init_d,
                            times = times,
                            parameters = pars_d,
                            bool_stochastic_beta = bool_stochastic_beta,
                            update_function = update_function,
                            num_experiments = n_stoch)
    
    summ <- sims$output_summary
    hc <- summ$HealthCost; sac <- summ$SocialActivityCost
    tot <- summ$TotalCost;  pt  <- summ$PeakTime
    
    psa_means[d, ] <- c(mean(hc,  na.rm = TRUE),
                        mean(sac, na.rm = TRUE),
                        mean(tot, na.rm = TRUE),
                        mean(pt,  na.rm = TRUE))
    
    psa_bands[d, ] <- c(mean(hc,  na.rm = TRUE), quantile(hc,  c(.025,.975), na.rm = TRUE),
                        mean(sac, na.rm = TRUE), quantile(sac, c(.025,.975), na.rm = TRUE),
                        mean(tot, na.rm = TRUE), quantile(tot, c(.025,.975), na.rm = TRUE),
                        mean(pt,  na.rm = TRUE), quantile(pt,  c(.025,.975), na.rm = TRUE))
    
    out_params[d, ] <- unlist(pars_d[keep_params])
  }
  
  list(psa_means = cbind(out_params, psa_means),
       psa_bands = cbind(out_params, psa_bands))
}

# ---------------------------
# PSA (full pooled runs)
# ---------------------------
run_psa_full <- function(initial_state, times, base_parameters,
                         n_param_draws = 60, n_stoch = 10,
                         bool_stochastic_beta = TRUE,
                         update_function = get_transitions_stochastic) {
  set.seed(base_parameters$rng_seed)
  
  tmp <- run_sir_binomial(initial_state, times, base_parameters,
                          bool_stochastic_beta, update_function)
  var_names <- colnames(tmp)
  Tlen <- nrow(tmp)
  total_runs <- n_param_draws * n_stoch
  
  output_all_psa <- array(NA_real_, dim = c(total_runs, Tlen, length(var_names)),
                          dimnames = list(NULL, NULL, var_names))
  
  run_idx <- 0
  for (d in 1:n_param_draws) {
    pars_d <- sample_parameters(base_parameters)
    init_d <- c(Ns = pars_d$ns0, Ni = pars_d$ni0, Nr = pars_d$nr0, Nd = pars_d$nd0,
                HealthCost = 0, SocialActivityCost = 0, TotalCost = 0,
                a_t = NA, u_t = NA, Rt = NA)
    pars_d$rng_seed <- base_parameters$rng_seed + d
    
    for (s in 1:n_stoch) {
      run_idx <- run_idx + 1
      sim <- run_sir_binomial(initial_state = init_d,
                              times = times,
                              parameters = pars_d,
                              bool_stochastic_beta = bool_stochastic_beta,
                              update_function = update_function)
      output_all_psa[run_idx, , ] <- as.matrix(sim)
    }
  }
  output_all_psa
}

# ---------------------------
# Dual-uncertainty bands plot (ggplot2)
# ---------------------------
# Requires: library(ggplot2)
plot_dual_uncertainty_bands <- function(baseline_output_all,   # 3D array from run_experiments()$output_all
                                        psa_output_all,        # 3D array from run_psa_full()
                                        deterministic_df,      # data.frame from deterministic run
                                        var = "Ni",
                                        title = NULL) {
  
  base_mat <- baseline_output_all[, , var, drop = TRUE]
  psa_mat  <- psa_output_all[, , var, drop = TRUE]
  
  time_vec <- 0:(ncol(base_mat) - 1)
  
  ribbon_summ <- function(M) {
    data.frame(
      time = time_vec,
      q025 = apply(M, 2, function(x) quantile(x, 0.025, na.rm = TRUE)),
      mean = apply(M, 2, function(x) mean(x, na.rm = TRUE)),
      q975 = apply(M, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
    )
  }
  
  base_sum <- ribbon_summ(base_mat)  # stochastic-only
  psa_sum  <- ribbon_summ(psa_mat)   # parameter + stochastic
  
  det_df <- data.frame(time = 0:(nrow(deterministic_df) - 1),
                       value = deterministic_df[[var]])
  
  ylab_txt <- switch(var,
                     "Ni" = "Infectives",
                     "Ns" = "Susceptibles",
                     "Nr" = "Recovered",
                     "Nd" = "Deceased",
                     "a_t" = "Activity",
                     "u_t" = "Utility",
                     "Rt" = "Reproduction Number",
                     "HealthCost" = "Health Cost",
                     "SocialActivityCost" = "Social Activity Cost",
                     "TotalCost" = "Total Cost",
                     var)
  
  if (is.null(title)) title <- paste("Dual Uncertainty Bands for", ylab_txt)
  
  ggplot() +
    geom_ribbon(data = psa_sum,
                aes(x = time, ymin = q025, ymax = q975),
                alpha = 0.25) +
    geom_ribbon(data = base_sum,
                aes(x = time, ymin = q025, ymax = q975),
                alpha = 0.40) +
    geom_line(data = psa_sum, aes(x = time, y = mean)) +
    geom_line(data = base_sum, aes(x = time, y = mean), linetype = "dashed") +
    geom_line(data = det_df, aes(x = time, y = value), linetype = "dotted") +
    labs(title = title, x = "Time", y = ylab_txt,
         subtitle = "Outer band: parameter + stochastic | Inner band: stochastic-only | Dotted: deterministic") +
    theme_minimal()
}

# ---------------------------
# Summary helpers
# ---------------------------
get_summary_stats <- function(sim_output){
  c(get_mean_ci_text(sim_output$HealthCost),
    get_mean_ci_text(sim_output$SocialActivityCost),
    get_mean_ci_text(sim_output$TotalCost))
}

get_mean_ci_text <- function(vect){
  if (length(vect) == 1) return(vect)
  summary_value <- ci_mean(vect)
  paste0(round(summary_value$estimate), ' [',
         paste(round(summary_value$interval), collapse = ','), ']')
}

# ---------------------------
# Compare stochastic vs deterministic (self-contained)
# ---------------------------
compare_sim_output <- function(output_experiments, output_deterministic,
                               plot_tag, fadeout_threshold = 0){
  
  output_summary <- output_experiments$output_summary
  output_all     <- output_experiments$output_all
  
  # Deterministic final row
  output_summary_deterministic <- output_deterministic[nrow(output_deterministic), ]
  
  # Optional fadeout filter
  if (fadeout_threshold > 0) {
    exp_fade_out   <- (output_summary$Nr + output_summary$Ni) < fadeout_threshold
    output_summary <- output_summary[!exp_fade_out, ]
    output_all     <- output_all[!exp_fade_out, , ]
    plot_tag <- paste(plot_tag, '(excl. fadeout)')
  } else {
    plot_tag <- paste(plot_tag, '(all)')
  }
  
  if (!dir.exists("figures")) dir.create("figures")
  pdf(file = paste0("figures/", gsub("[^[:alnum:]]", "_", plot_tag), ".pdf"))
  
  # --- Boxplot of costs
  output_cost <- output_summary[, grepl('Cost', names(output_summary))]
  y_lim <- range(output_cost, na.rm = TRUE)
  boxplot(output_cost, las = 2, ylim = y_lim, main = plot_tag)
  points(colMeans(output_cost), pch = 8) # stochastic mean
  points((1:3) + 0.2, output_summary_deterministic[names(output_cost)], col = 4, pch = 8, lwd = 3)
  
  # --- Time-series bands (explicit selection; no globals)
  sel_states <- c("Ni","Nr","Nd","a_t","u_t","Rt","HealthCost","SocialActivityCost","TotalCost")
  
  pretty_labels <- c(
    "a_t" = "Activity",
    "Ni" = "Infectives", "Nr" = "Recovered", "Nd" = "Deceased",
    "u_t" = "Utility", "Rt" = "Reproduction Number",
    "HealthCost" = "Health Cost", "SocialActivityCost" = "Social Activity Cost",
    "TotalCost" = "Total Cost"
  )
  
  for (i_state in sel_states) {
    if (!i_state %in% dimnames(output_all)[[3]]) next
    
    plot_label <- ifelse(i_state %in% names(pretty_labels), pretty_labels[[i_state]], i_state)
    var_matrix <- output_all[, , i_state, drop = TRUE]
    time_vec <- 0:(ncol(var_matrix) - 1)
    
    # Bands
    ci_lower <- apply(var_matrix, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
    ci_upper <- apply(var_matrix, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
    n_sim <- nrow(var_matrix)
    ci_mean <- apply(var_matrix, 2, mean, na.rm = TRUE)
    ci_sd <- apply(var_matrix, 2, sd, na.rm = TRUE)
    ci_se <- ci_sd / sqrt(n_sim)
    ci_lo_meanband <- ci_mean - 1.96 * ci_se
    ci_hi_meanband <- ci_mean + 1.96 * ci_se
    
    plot(NULL,
         xlim = range(time_vec),
         ylim = range(ci_lower, ci_upper, ci_lo_meanband, ci_hi_meanband,
                      output_deterministic[[i_state]], na.rm = TRUE),
         xlab = 'time', ylab = plot_label, main = plot_label)
    
    polygon(c(time_vec, rev(time_vec)), c(ci_upper, rev(ci_lower)),
            col = alpha("lightblue", 0.4), border = NA)
    
    polygon(c(time_vec, rev(time_vec)), c(ci_hi_meanband, rev(ci_lo_meanband)),
            col = alpha("red", 0.4), border = NA)
    
    lines(time_vec, ci_mean, col = "blue", lwd = 2)
    lines(0:(nrow(output_deterministic) - 1), output_deterministic[[i_state]],
          col = "black", lwd = 2, lty = 2)
  }
  
  # --- Console summaries
  mean_health <- mean(output_summary$HealthCost, na.rm = TRUE)
  ci_health <- ci_mean(output_summary$HealthCost)$interval
  
  mean_social <- mean(output_summary$SocialActivityCost, na.rm = TRUE)
  ci_social <- ci_mean(output_summary$SocialActivityCost)$interval
  
  cat("Stochastic Summary (from boxplot data):\n")
  cat(sprintf("Health Cost: %.0f [%.0f, %.0f]\n", mean_health, ci_health[1], ci_health[2]))
  cat(sprintf("Social Activity Cost: %.0f [%.0f, %.0f]\n", mean_social, ci_social[1], ci_social[2]))
  
  # % differences
  health_quant <- quantile(output_summary$HealthCost, probs = c(0.025, 0.975), na.rm = TRUE)
  pct_health_upper <- 100 * (health_quant[2] - mean_health) / mean_health
  pct_health_lower <- 100 * (mean_health - health_quant[1]) / mean_health
  
  mean_activity_cost <- mean(output_summary$SocialActivityCost, na.rm = TRUE)
  activity_quant <- quantile(output_summary$SocialActivityCost, probs = c(0.025, 0.975), na.rm = TRUE)
  pct_activity_upper <- 100 * (activity_quant[2] - mean_activity_cost) / mean_activity_cost
  pct_activity_lower <- 100 * (mean_activity_cost - activity_quant[1]) / mean_activity_cost
  
  health_det <- output_deterministic[nrow(output_deterministic), "HealthCost"]
  activity_det <- output_deterministic[nrow(output_deterministic), "SocialActivityCost"]
  pct_health_mean_diff <- 100 * (health_det - mean_health) / health_det
  pct_activity_mean_diff <- 100 * (activity_det - mean_activity_cost) / activity_det
  
  # Peak-time metrics
  t_peak <- which.max(output_deterministic$Ni)
  a_vals <- output_all[, t_peak, "a_t"]
  ni_vals <- output_all[, t_peak, "Ni"]
  
  a_mean <- mean(a_vals, na.rm = TRUE); ni_mean <- mean(ni_vals, na.rm = TRUE)
  a_q <- quantile(a_vals, probs = c(0.025, 0.975), na.rm = TRUE)
  ni_q <- quantile(ni_vals, probs = c(0.025, 0.975), na.rm = TRUE)
  
  a_det <- output_deterministic[t_peak, "a_t"]
  ni_det <- output_deterministic[t_peak, "Ni"]
  
  pct_a_upper <- 100 * (a_q[2] - a_mean) / a_mean
  pct_a_lower <- 100 * (a_mean - a_q[1]) / a_mean
  pct_ni_upper <- 100 * (ni_q[2] - ni_mean) / ni_mean
  pct_ni_lower <- 100 * (ni_mean - ni_q[1]) / ni_mean
  
  pct_a_mean_diff <- 100 * (a_det - a_mean) / a_det
  pct_ni_mean_diff <- 100 * (ni_det - ni_mean) / ni_det
  
  cat("\nCOST DIFFERENCES (BOXPLOT METRICS)\n")
  cat(sprintf("Health Cost – Credible Interval: +%.2f%% / -%.2f%%\n", pct_health_upper, pct_health_lower))
  cat(sprintf("Activity Cost – Credible Interval: +%.2f%% / -%.2f%%\n", pct_activity_upper, pct_activity_lower))
  cat(sprintf("Health Cost – Stochastic vs Det: %.2f%%\n", pct_health_mean_diff))
  cat(sprintf("Activity Cost – Stochastic vs Det: %.2f%%\n", pct_activity_mean_diff))
  
  cat("\nPEAK-TIME DIFFERENCES (Ni and a_t at peak Ni)\n")
  cat(sprintf("Activity a(t) – Credible Interval: +%.2f%% / -%.2f%%\n", pct_a_upper, pct_a_lower))
  cat(sprintf("Infectives Ni – Credible Interval: +%.2f%% / -%.2f%%\n", pct_ni_upper, pct_ni_lower))
  cat(sprintf("Activity a(t) – Stochastic vs Det: %.2f%%\n", pct_a_mean_diff))
  cat(sprintf("Infectives Ni – Stochastic vs Det: %.2f%%\n", pct_ni_mean_diff))
  
  # Cost boxplots with CI overlays + deterministic points
  boxplot(output_summary$HealthCost,
          main = "Health Cost", ylab = "Health Cost",
          ylim = range(c(output_summary$HealthCost, output_summary_deterministic$HealthCost), na.rm = TRUE),
          col = "gray90")
  points(1, mean(output_summary$HealthCost, na.rm = TRUE), pch = 8)
  ci_health <- ci_mean(output_summary$HealthCost)
  segments(1, ci_health$interval[1], 1, ci_health$interval[2], col = "forestgreen", lwd = 2)
  quant_health <- quantile(output_summary$HealthCost, probs = c(0.025, 0.975), na.rm = TRUE)
  segments(1.1, quant_health[1], 1.1, quant_health[2], col = "blue", lty = 2, lwd = 2)
  points(1.2, output_summary_deterministic$HealthCost, col = "red", pch = 18, cex = 1.5)
  
  boxplot(output_summary$SocialActivityCost,
          main = "Social Activity Cost", ylab = "Social Activity Cost",
          ylim = range(c(output_summary$SocialActivityCost, output_summary_deterministic$SocialActivityCost), na.rm = TRUE),
          col = "gray90")
  points(1, mean(output_summary$SocialActivityCost, na.rm = TRUE), pch = 8)
  ci_social <- ci_mean(output_summary$SocialActivityCost)
  segments(1, ci_social$interval[1], 1, ci_social$interval[2], col = "forestgreen", lwd = 2)
  quant_social <- quantile(output_summary$SocialActivityCost, probs = c(0.025, 0.975), na.rm = TRUE)
  segments(1.1, quant_social[1], 1.1, quant_social[2], col = "blue", lty = 2, lwd = 2)
  points(1.2, output_summary_deterministic$SocialActivityCost, col = "red", pch = 18, cex = 1.5)
  
  # Printable summary table
  print_out <- data.frame(
    output = c('Health Cost (per capita)', 'Social Activity Cost (per capita)', 'Total Cost (per capita)'),
    deterministic = get_summary_stats(output_deterministic[nrow(output_deterministic), ]),
    stochastic    = get_summary_stats(output_summary)
  )
  names(print_out)[1] <- plot_tag
  print(paste('SUMMARY STATISTICS:', plot_tag))
  print(print_out)
  
  # Peak time summary
  peak_times <- output_summary$PeakTime
  peak_ci <- ci_mean(peak_times)
  cat("\nPeak Infection Timing Summary:\n")
  cat(sprintf("Mean Peak Time: %.1f [%.1f, %.1f]  (mean ± 1.96 * SE)\n",
              peak_ci$estimate, peak_ci$interval[1], peak_ci$interval[2]))
  cat(sprintf("Quantiles: 2.5%% = %.1f, 50%% (median) = %.1f, 97.5%% = %.1f\n",
              quantile(peak_times, 0.025, na.rm = TRUE),
              quantile(peak_times, 0.50, na.rm = TRUE),
              quantile(peak_times, 0.975, na.rm = TRUE)))
  
  dev.off()
}

# ---------------------------
# Experiment runner
# ---------------------------
run_experiments <- function(initial_state, times, parameters,
                            bool_stochastic_beta, update_function,
                            num_experiments)
{
  set.seed(parameters$rng_seed)
  
  # Learn structure
  temp_output <- run_sir_binomial(initial_state, times, parameters,
                                  bool_stochastic_beta, update_function)
  
  # Output summary + PeakTime
  output_summary <- data.frame(matrix(NA_real_, nrow = num_experiments,
                                      ncol = ncol(temp_output) + 1))
  names(output_summary) <- c(names(temp_output), "PeakTime")
  
  # Full time-series array
  output_all <- array(NA_real_,
                      dim = c(num_experiments, length(times), ncol(temp_output)),
                      dimnames = list(NULL, NULL, names(temp_output)))
  
  for (i_exp in 1:num_experiments) {
    # message(i_exp)
    output_sim <- run_sir_binomial(initial_state = initial_state,
                                   times = times,
                                   parameters = parameters,
                                   bool_stochastic_beta = bool_stochastic_beta,
                                   update_function = update_function)
    
    output_summary[i_exp, 1:ncol(temp_output)] <- output_sim[nrow(output_sim), ]
    output_summary$PeakTime[i_exp] <- which.max(output_sim$Ni) - 1
    output_all[i_exp, , ] <- as.matrix(output_sim)
  }
  
  list(output_summary = output_summary, output_all = output_all)
}
