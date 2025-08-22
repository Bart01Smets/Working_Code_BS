##
# EPI-ECON helper (NO HOSPITALIZATION) — imperfect recovery awareness (zeta)
# Drop-in replacement for epi_econ_lib_hosprec.R
# - No Nh state or hospital admissions/resolution
# - People exiting I resolve directly to R or D
# - Activity rule uses imperfect awareness (zeta) and quadratic closed-form
# - Costs and plotting utilities compatible with the workbench
# ------------------------------------------------------------------------

# Dependencies
suppressWarnings({
  if (!requireNamespace("confintr", quietly = TRUE)) {
    stop("Package 'confintr' is required. Please install.packages('confintr').")
  }
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required. Please install.packages('scales').")
  }
})
library(confintr)
library(scales)

# ----------------------------
# Activity rule: imperfect info (no hospital)
# ----------------------------
a_function <- function(Ni, Ns, parameters) {
  if (isTRUE(parameters$bool_regular_sird)) return(1)

  Ni_prop <- Ni / parameters$pop_size
  Ns_prop <- Ns / parameters$pop_size
  if (!is.finite(Ni_prop) || !is.finite(Ns_prop) || Ni_prop <= 0 || Ns_prop <= 0) return(1)

  zeta <- if (!is.null(parameters$zeta)) parameters$zeta else 0
  zeta <- max(0, min(1, zeta))  # clamp

  # Mass choosing a(t) and effective per-infection cost under imperfect info
  W <- zeta + (1 - zeta) * (Ns_prop + Ni_prop)             # chooser mass
  kappa_tilde <- (1 - zeta) * parameters$pi * parameters$v # perceived per-inf cost

  # Quadratic: D*A^2 + W*A - W = 0  with D = beta * kappa_tilde * Ns * Ni
  D <- parameters$beta * kappa_tilde * Ns_prop * Ni_prop
  if (!is.finite(W) || !is.finite(D) || D <= 0 || W <= 0) return(1)

  disc <- W * W + 4 * D * W
  if (!is.finite(disc) || disc < 0) return(1)

  a_t <- (-W + sqrt(disc)) / (2 * D)
  a_t <- max(0, min(1, a_t))
  if (!is.finite(a_t)) a_t <- 1
  return(a_t)
}

# ----------------------------
# Utility and Rt
# ----------------------------
utility_function <- function(a_t) {
  # Returns 0 when using regular SIRD (no econ)
  if (exists("parameters", inherits = TRUE) && isTRUE(get("parameters")$bool_regular_sird)) return(0)
  return(log(a_t) - a_t + 1)
}

calculate_Rt <- function(R0, a_t, Ns_prop, Ni) {
  if (Ni <= 0) return(0)
  return(R0 * a_t^2 * Ns_prop)
}

# ----------------------------
# Transition helpers
# ----------------------------
get_transitions_stochastic <- function(n, prob) {
  if (is.na(n) || is.na(prob) || !is.numeric(n) || !is.numeric(prob) || n < 0 || prob < 0 || prob > 1) {
    warning("Invalid rbinom() inputs: returning 0")
    return(0)
  }
  return(sum(rbinom(n = n, size = 1, prob = prob)))
}

get_transitions_deterministic <- function(n, prob){
  if (length(prob) == n) {
    transitions <- sum(prob)
  } else{
    transitions <- (n * prob)
  }
  return(transitions)
}

# ----------------------------
# Main kernel (no hospitalization)
# ----------------------------
run_sir_binomial <- function(initial_state,
                             times,
                             parameters,
                             bool_stochastic_beta = FALSE,
                             update_function = get_transitions_stochastic){

  # copy & convert initial fractions to counts
  states <- data.frame(t(initial_state))
  states[grepl('N.', names(states))] <- round(states[grepl('N.', names(states))] * parameters$pop_size)
  # keep total pop fixed by adjusting Ns
  states[grepl('Ns', names(states))] <- parameters$pop_size -
    sum(states[grepl('N.', names(states)) & !grepl('Ns', names(states))])

  # output holder
  states_out <- matrix(NA, nrow = length(times), ncol = length(states))
  colnames(states_out) <- names(states)

  # local copies
  Ns <- states$Ns; Ni <- states$Ni; Nr <- states$Nr; Nd <- states$Nd
  HealthCost <- 0; SocialActivityCost <- 0
  fx_per_capita <- parameters$fx / parameters$pop_size

  # write initial row
  states_out[1,] <- c(Ns, Ni, Nr, Nd,
                      HealthCost, SocialActivityCost, HealthCost + SocialActivityCost,
                      NA, NA, NA)

  for (i_day in times[-1]) {
    # stochastic beta (optional)
    beta_t <- if (bool_stochastic_beta) max(1e-4, rnorm(1, mean = parameters$beta, sd = parameters$sigma)) else parameters$beta

    # activity & utility
    a_t <- a_function(Ni, Ns, parameters)
    if (isTRUE(parameters$bool_daily_cost_minimizing) && exists("a_function_cost_minimizing")) {
      a_t <- a_function_cost_minimizing(Ns, Ni, i_day, parameters)
    }
    u_t <- utility_function(a_t)

    # capacity-adjusted mortality (optional)
    if (!is.null(parameters$healthcare_capacity) && Ni > parameters$healthcare_capacity) {
      pi_effective <- parameters$pi * parameters$excess_mortality_multiplier
    } else {
      pi_effective <- parameters$pi
    }

    # infections
    p_infect <- 1 - exp(- beta_t * a_t^2 * Ni / parameters$pop_size)
    new_infections <- update_function(Ns, prob = p_infect)

    # exits from I (direct resolution to R or D)
    p_exit_I <- 1 - exp(-parameters$gamma)
    exits_I  <- update_function(Ni, prob = p_exit_I)

    new_death      <- if (exits_I > 0) update_function(exits_I, prob = pi_effective) else 0
    new_recoveries <- exits_I - new_death

    # state deltas
    dNs <- -new_infections
    dNi <-  new_infections - exits_I
    dNr <-  new_recoveries
    dNd <-  new_death

    # fadeout guard
    if ((Ni + dNi) < parameters$infect_thres) {
      rem <- max(0, Ni + dNi)
      dNr <- dNr + rem
      dNi <- -Ni
    }

    # discounted costs (per capita)
    if (!isTRUE(parameters$bool_regular_sird)) {
      zeta <- if (!is.null(parameters$zeta)) parameters$zeta else 0
      W_mass <- zeta * parameters$pop_size + (1 - zeta) * (Ns + Ni)
      HealthCost         <- HealthCost + fx_per_capita * exp(-parameters$rho * i_day) * parameters$v * new_death
      SocialActivityCost <- SocialActivityCost + fx_per_capita * exp(-parameters$rho * i_day) * W_mass * abs(u_t)
    }

    # update states
    Ns <- Ns + dNs; Ni <- Ni + dNi; Nr <- Nr + dNr; Nd <- Nd + dNd
    Rt <- calculate_Rt(parameters$R0, a_t, Ns / parameters$pop_size, Ni)

    states_out[i_day + 1, ] <- c(Ns, Ni, Nr, Nd,
                                 HealthCost, SocialActivityCost, HealthCost + SocialActivityCost,
                                 a_t, u_t, Rt)
  }

  data.frame(states_out)
}

# ----------------------------
# Summary helpers
# ----------------------------
get_mean_ci_text <- function(vect){
  if(length(vect) == 1){
    return(vect)
  } else{
    summary_value <- confintr::ci_mean(vect)
    paste0(round(summary_value$estimate),' [',
           paste(round(summary_value$interval),collapse=','),']')
  }
}

get_summary_stats <- function(sim_output){
  c(get_mean_ci_text(sim_output$HealthCost),
    get_mean_ci_text(sim_output$SocialActivityCost),
    get_mean_ci_text(sim_output$TotalCost))
}

# ----------------------------
# Comparison figure & console summary
# ----------------------------
compare_sim_output <- function(output_experiments, output_deterministic,
                               plot_tag, fadeout_threshold = 0){

  output_summary <- output_experiments$output_summary
  output_all     <- output_experiments$output_all

  # option to exclude stochastic fadeout
  if(fadeout_threshold > 0){
    exp_fade_out   <- (output_summary$Nr + output_summary$Ni) < fadeout_threshold
    output_summary <- output_summary[!exp_fade_out, , drop = FALSE]
    output_all     <- output_all[!exp_fade_out, , , drop = FALSE]
    plot_tag <- paste(plot_tag,'(excl. fadeout)')
  } else {
    plot_tag <- paste(plot_tag,'(all)')
  }

  # ensure figures dir
  if (!dir.exists("figures")) dir.create("figures")
  pdf(file = paste0("figures/", gsub("[^[:alnum:]]", "_", plot_tag), ".pdf"))

  # COST BOXES
  output_cost <- output_summary[, grepl('Cost', names(output_summary)), drop = FALSE]
  y_lim <- range(output_cost, na.rm = TRUE)
  boxplot(output_cost, las = 2, ylim = y_lim, main = plot_tag)
  points(colMeans(output_cost, na.rm = TRUE), pch = 8) # mean
  points((1:NCOL(output_cost)) + 0.2,
         as.numeric(output_deterministic[nrow(output_deterministic), names(output_cost)]),
         col = 4, pch = 8, lwd = 3)

  # Which variables to time-plot if present
  preferred <- c("a_t","Ni","Nr","Nd","u_t","Rt","HealthCost","SocialActivityCost","TotalCost")
  present <- intersect(preferred, dimnames(output_all)[[3]])

  time_vec <- 0:(dim(output_all)[2]-1)

  for (i_state in present){
    plot_label <- switch(i_state,
                         a_t = "Activity", Ni = "Infectives", Nr = "Recovered",
                         Nd = "Deceased", u_t = "Utility", Rt = "Reproduction Number",
                         HealthCost = "Health Cost", SocialActivityCost = "Social Activity Cost",
                         TotalCost = "Total Cost", i_state)

    var_matrix <- output_all[,,i_state, drop = TRUE]

    ci_lower <- apply(var_matrix, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
    ci_upper <- apply(var_matrix, 2, function(x) quantile(x, 0.975, na.rm = TRUE))

    n_sim <- nrow(var_matrix)
    ci_mean <- apply(var_matrix, 2, mean, na.rm = TRUE)
    ci_sd <- apply(var_matrix, 2, sd, na.rm = TRUE)
    ci_se <- ci_sd / sqrt(n_sim)
    ci_lo_meanband <- ci_mean - 1.96 * ci_se
    ci_hi_meanband <- ci_mean + 1.96 * ci_se

    ylim_all <- range(ci_lower, ci_upper, ci_lo_meanband, ci_hi_meanband,
                      output_deterministic[, i_state], na.rm = TRUE)

    plot(NULL, xlim = range(time_vec), ylim = ylim_all,
         xlab = 'time', ylab = plot_label, main = plot_label)
    polygon(c(time_vec, rev(time_vec)), c(ci_upper, rev(ci_lower)), col = scales::alpha("lightblue", 0.4), border = NA)
    polygon(c(time_vec, rev(time_vec)), c(ci_hi_meanband, rev(ci_lo_meanband)), col = scales::alpha("red", 0.4), border = NA)
    lines(time_vec, ci_mean, col = "blue", lwd = 2)
    lines(output_deterministic[, i_state], col = "black", lwd = 2, lty = 2)
  }

  # Console summaries
  mean_health <- mean(output_summary$HealthCost, na.rm = TRUE)
  ci_health <- confintr::ci_mean(output_summary$HealthCost)$interval
  mean_social <- mean(output_summary$SocialActivityCost, na.rm = TRUE)
  ci_social <- confintr::ci_mean(output_summary$SocialActivityCost)$interval

  cat("Stochastic Summary:\n")
  cat(sprintf("Health Cost: %.0f [%.0f, %.0f]\n", mean_health, ci_health[1], ci_health[2]))
  cat(sprintf("Social Activity Cost: %.0f [%.0f, %.0f]\n", mean_social, ci_social[1], ci_social[2]))

  # Credible interval % widths
  health_quant <- quantile(output_summary$HealthCost, probs = c(0.025, 0.975), na.rm = TRUE)
  pct_health_upper <- 100 * (health_quant[2] - mean_health) / mean_health
  pct_health_lower <- 100 * (mean_health - health_quant[1]) / mean_health

  activity_quant <- quantile(output_summary$SocialActivityCost, probs = c(0.025, 0.975), na.rm = TRUE)
  pct_activity_upper <- 100 * (activity_quant[2] - mean_social) / mean_social
  pct_activity_lower <- 100 * (mean_social - activity_quant[1]) / mean_social

  # Mean vs deterministic (final)
  health_det <- output_deterministic[nrow(output_deterministic), "HealthCost"]
  activity_det <- output_deterministic[nrow(output_deterministic), "SocialActivityCost"]
  pct_health_mean_diff <- 100 * (health_det - mean_health) / health_det
  pct_activity_mean_diff <- 100 * (activity_det - mean_social) / activity_det

  # Peak-time differences at deterministic peak
  t_peak <- which.max(output_deterministic$Ni)
  a_vals <- output_all[, t_peak, "a_t", drop = TRUE]
  ni_vals <- output_all[, t_peak, "Ni", drop = TRUE]
  a_mean <- mean(a_vals, na.rm = TRUE); ni_mean <- mean(ni_vals, na.rm = TRUE)
  a_q <- quantile(a_vals, probs = c(0.025, 0.975), na.rm = TRUE)
  ni_q <- quantile(ni_vals, probs = c(0.025, 0.975), na.rm = TRUE)
  a_det <- output_deterministic[t_peak, "a_t"]; ni_det <- output_deterministic[t_peak, "Ni"]

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

  # Additional cost boxplots with CI/quantile bars
  boxplot(output_summary$HealthCost,
          main = "Health Cost", ylab = "Health Cost",
          ylim = range(c(output_summary$HealthCost, output_deterministic$HealthCost), na.rm = TRUE),
          col = "gray90")
  points(1, mean(output_summary$HealthCost, na.rm = TRUE), pch = 8)
  seg_h <- confintr::ci_mean(output_summary$HealthCost)
  segments(x0 = 1, x1 = 1, y0 = seg_h$interval[1], y1 = seg_h$interval[2], col = "forestgreen", lwd = 2)
  q_h <- quantile(output_summary$HealthCost, probs = c(0.025, 0.975), na.rm = TRUE)
  segments(x0 = 1.1, x1 = 1.1, y0 = q_h[1], y1 = q_h[2], col = "blue", lwd = 2, lty = 2)
  points(1.2, output_deterministic[nrow(output_deterministic), "HealthCost"], col = "red", pch = 18, cex = 1.5)

  boxplot(output_summary$SocialActivityCost,
          main = "Social Activity Cost", ylab = "Social Activity Cost",
          ylim = range(c(output_summary$SocialActivityCost, output_deterministic$SocialActivityCost), na.rm = TRUE),
          col = "gray90")
  points(1, mean(output_summary$SocialActivityCost, na.rm = TRUE), pch = 8)
  seg_s <- confintr::ci_mean(output_summary$SocialActivityCost)
  segments(x0 = 1, x1 = 1, y0 = seg_s$interval[1], y1 = seg_s$interval[2], col = "forestgreen", lwd = 2)
  q_s <- quantile(output_summary$SocialActivityCost, probs = c(0.025, 0.975), na.rm = TRUE)
  segments(x0 = 1.1, x1 = 1.1, y0 = q_s[1], y1 = q_s[2], col = "blue", lwd = 2, lty = 2)
  points(1.2, output_deterministic[nrow(output_deterministic), "SocialActivityCost"], col = "red", pch = 18, cex = 1.5)

  # Print summary table
  print_out <- data.frame(
    output = c('Health Cost (per capita)','Social Activity Cost (per capita)','Total Cost (per capita)'),
    deterministic = c(
      get_mean_ci_text(output_deterministic[nrow(output_deterministic), "HealthCost"]),
      get_mean_ci_text(output_deterministic[nrow(output_deterministic), "SocialActivityCost"]),
      get_mean_ci_text(output_deterministic[nrow(output_deterministic), "TotalCost"])
    ),
    stochastic = c(
      get_mean_ci_text(output_summary$HealthCost),
      get_mean_ci_text(output_summary$SocialActivityCost),
      get_mean_ci_text(output_summary$TotalCost)
    )
  )
  names(print_out)[1] <- plot_tag
  print(paste('SUMMARY STATISTICS:', plot_tag))
  print(print_out)

  # Peak infection timing summary across sims
  if ("PeakTime" %in% names(output_summary)) {
    peak_ci <- confintr::ci_mean(output_summary$PeakTime)
    cat("\nPeak Infection Timing Summary:\n")
    cat(sprintf("Mean Peak Time: %.1f [%.1f, %.1f]  (mean ± 1.96 * SE)\n",
                peak_ci$estimate, peak_ci$interval[1], peak_ci$interval[2]))
    cat(sprintf("Quantiles: 2.5%% = %.1f, 50%% (median) = %.1f, 97.5%% = %.1f\n",
                quantile(output_summary$PeakTime, 0.025, na.rm = TRUE),
                quantile(output_summary$PeakTime, 0.50, na.rm = TRUE),
                quantile(output_summary$PeakTime, 0.975, na.rm = TRUE)))
  }

  dev.off()
}

# ----------------------------
# Experiment harness
# ----------------------------
run_experiments <- function(initial_state, times, parameters,
                            bool_stochastic_beta, update_function,
                            num_experiments)
{
  set.seed(parameters$rng_seed)

  temp_output <- run_sir_binomial(initial_state, times, parameters, bool_stochastic_beta, update_function)

  output_summary <- data.frame(matrix(NA, nrow = num_experiments, ncol = ncol(temp_output) + 1))
  names(output_summary) <- c(names(temp_output), "PeakTime")

  output_all <- array(NA, dim = c(num_experiments, length(times), ncol(temp_output)),
                      dimnames = list(NULL, NULL, names(temp_output)))

  for(i_exp in 1:num_experiments){
    message(sprintf("Experiment %d/%d", i_exp, num_experiments))
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
