
##
# Help function for the EPI-ECON modelling work
#
# - Deterministic and stochastic transitions
################################################################ #
#setwd("C:/Users/Bart Smets/OneDrive/Documenten/GitHub/Working_Code_BS/R")
library(confintr)
library(scales)

a_function <- function(Ni, Ns, parameters) {
  if (isTRUE(parameters$bool_regular_sird)) return(1)  # Always full activity
  # Convert to proportions
  Ni_prop <- Ni / parameters$pop_size
  Ns_prop <- Ns / parameters$pop_size
  
  # Defensive early exit if Ni or Ns is ~0
  if (Ni_prop < 1e-10 || Ns_prop < 1e-10) return(1)
  
  # Define terms for quadratic
  multiplier <- parameters$beta * parameters$pi*parameters$v * Ns_prop*Ni_prop
  num_sum <- Ns_prop+Ni_prop
  discrim <- num_sum^2 + 4 * multiplier * num_sum
  denom <- 2 * multiplier
  
  if (denom <= 0 || discrim < 0) return(1)
  
  a_t <- (-num_sum + sqrt(discrim)) / denom
  
  # Bound a_t to [0,1]
  return(max(0, min(1, a_t)))
}

# Activity rule with imperfect recovery information (zeta)
# Uses: W = zeta + (1 - zeta) * (n_s + n_i)
#       kappa_tilde = (1 - zeta) * pi * v
a_function <- function(Ni, Ns, parameters) {
  if (isTRUE(parameters$bool_regular_sird)) return(1)
  
  # proportions
  Ni_prop <- Ni / parameters$pop_size
  Ns_prop <- Ns / parameters$pop_size
  
  # defensive exits
  if (!is.finite(Ni_prop) || !is.finite(Ns_prop) ||
      Ni_prop <= 0 || Ns_prop <= 0) return(1)
  
  # imperfect information about recovery
  zeta <- if (!is.null(parameters$zeta)) parameters$zeta else 0
  
  # weights and effective per-infection cost
  W <- zeta + (1 - zeta) * (Ns_prop + Ni_prop)                 # mass choosing a(t)
  kappa_tilde <- (1 - zeta) * parameters$pi * parameters$v     # κ̃
  
  # quadratic coefficients: D * A^2 + W * A - W = 0
  D <- parameters$beta * kappa_tilde * Ns_prop * Ni_prop
  
  # guard rails
  if (!is.finite(W) || !is.finite(D) || D <= 0 || W <= 0) return(1)
  
  discrim <- W * W + 4 * D * W
  if (!is.finite(discrim) || discrim < 0) return(1)
  
  a_t <- (-W + sqrt(discrim)) / (2 * D)
  
  # bound to [0, 1]
  a_t <- max(0, min(1, a_t))
  if (!is.finite(a_t)) a_t <- 1
  return(a_t)
}
#Hospitalization
a_function <- function(Ni, Ns, parameters) {
  if (isTRUE(parameters$bool_regular_sird)) return(1)

  Ni_prop <- Ni / parameters$pop_size
  Ns_prop <- Ns / parameters$pop_size
  if (Ni_prop <= 0 || Ns_prop <= 0) return(1)

  zeta <- parameters$zeta
  rho  <- parameters$rho
  xi   <- parameters$xi

  # Mass choosing a(t):
  W <- zeta + (1 - zeta) * (Ns_prop + Ni_prop)

  # Hospitalization-adjusted present-value infection cost κ_h
  kappa_h <- (1 - zeta) * parameters$pi * xi * parameters$v / (rho + xi)  # C.3
  D <- parameters$beta * kappa_h * Ns_prop * Ni_prop

  if (!is.finite(W) || !is.finite(D) || D <= 0 || W <= 0) return(1)
  disc <- W*W + 4*D*W
  if (disc < 0) return(1)

  a_t <- (-W + sqrt(disc)) / (2 * D)
  return(max(0, min(1, a_t)))
}


# Utility function
utility_function <- function(a_t) {
  #if (utility_type == "Log") {
  if (isTRUE(parameters$bool_regular_sird)) return(0)
  return(log(a_t) - a_t + 1)
  #} else {
  #  return(-1/2 * (1 - a_t)^2)
}


#Calculation of effective reproduction number
calculate_Rt <- function(R0, a_t, Ns_prop, Ni) {
  if (Ni <= 0) return(0)
  return(R0 * a_t^2 * Ns_prop)
}

get_transitions_stochastic <- function(n, prob) {
  #cat("DEBUG — rbinom() inputs:\n")
  #cat("  n =", n, "\n")
  #cat("  prob =", prob, "\n")
  
  if (is.na(n) || is.na(prob) || !is.numeric(n) || !is.numeric(prob) || n < 0 || prob < 0 || prob > 1) {
    warning("Invalid rbinom() inputs: returning 0")
    return(0)
  }
  
  return(sum(rbinom(n = n, size = 1, prob = prob)))
}


get_transitions_deterministic <- function(n, prob){
  if(length(prob)==n){
    transitions=sum(prob)
  } else{
    transitions<-(n*prob)
  }
  return(transitions)
}

# to run the SIRD kernel with a for loop and binomial transitions
run_sir_binomial <- function(initial_state,
                             times,
                             parameters,
                             bool_stochastic_beta = FALSE,
                             update_function = get_transitions_stochastic){
  
  # copy initial states
  states <- data.frame(t(initial_state))
  
  # convert population fractions into numbers
  states[grepl('N.',names(states))] <- round(states[grepl('N.',names(states))] * parameters$pop_size)
  
  # make sure the total population size is fixed (by changing Ns if needed)
  states[grepl('Ns',names(states))] <- parameters$pop_size - sum(states[grepl('N.',names(states)) & !grepl('Ns',names(states))])
  
  # set summary data.table
  states_out        <- matrix(NA, nrow=length(times), ncol=length(states))
  colnames(states_out) <- names(states)
  dim(states_out)
  
  # start with initial states.
  # note: use single values as parameter is much faster in the iterative process
  names(states)
  Ns <- states$Ns
  Ni <- states$Ni
  Nr <- states$Nr
  Nd <- states$Nd
  Nh <- if ("Nh" %in% names(states)) states$Nh else 0   # <<< ADD THIS
  HealthCost <- 0
  SocialActivityCost <- 0
  
  # compute some results only once
  fx_per_capita <- parameters$fx / parameters$pop_size
  
  # run over times (excl the first one)
  for(i_day in times[-1]){
    
    #  beta_t <- parameters$beta
    if(bool_stochastic_beta){
      beta_t <- max(0.0001, rnorm(1, mean = parameters$beta, sd = parameters$sigma))
    } else {
      beta_t <- parameters$beta
    }
    
    
    Ns_prop <- Ns / parameters$pop_size
    # a_t <- a_function(Ni, Ns, parameters)
    if (isTRUE(parameters$bool_daily_cost_minimizing)) {
      a_t <- a_function_cost_minimizing(Ns, Ni, i_day, parameters)
    } else {
      a_t <- a_function(Ni, Ns, parameters)
    }
    
    # Calculate utility of action
    u_t <- utility_function(a_t)
    # Adjust death rate if over healthcare capacity
    if (!is.null(parameters$healthcare_capacity) && Ni > parameters$healthcare_capacity) {
      pi_effective <- parameters$pi * parameters$excess_mortality_multiplier
    } else {
      pi_effective <- parameters$pi
    }
    
    
    
    # --- replace the current recover/death block with this ---
    xi   <- parameters$xi
    zeta <- parameters$zeta
    
    # 1) Infection probability (unchanged)
    p_infect <- 1 - exp(- beta_t * a_t^2 * Ni / parameters$pop_size)
    new_infections <- update_function(Ns, prob = p_infect)
    
    # 2) Exits from I this day
    p_exit_I <- 1 - exp(-parameters$gamma)
    exits_I  <- update_function(Ni, prob = p_exit_I)
    
    #   Split exits: to hospital vs unaware recovery (joins R but still counted in W)
    p_enter_H <- (1 - zeta)
    new_admit_H <- if (exits_I > 0) update_function(exits_I, prob = p_enter_H) else 0
    new_unaware_rec <- exits_I - new_admit_H
    
    # 3) Resolve hospital stock
    p_resolveH    <- 1 - exp(-xi)
    resolves_H    <- if (Nh > 0) update_function(Nh, prob = p_resolveH) else 0
    p_die_from_H  <- pi_effective
    new_death     <- if (resolves_H > 0) update_function(resolves_H, prob = p_die_from_H) else 0
    new_rec_from_H <- resolves_H - new_death
    
    # 4) State deltas
    dNs <- -new_infections
    dNi <-  new_infections - exits_I
    dNh <-  new_admit_H - resolves_H
    dNr <-  new_unaware_rec + new_rec_from_H
    dNd <-  new_death
    # --- end replacement ---
    
    
    if((Ni + dNi) < parameters$infect_thres){
      # new_recoveries = Ni
      # new_infections = 0
      
      rem <- max(0, Ni + dNi)
      dNr <- dNr + rem
      dNi <- -Ni
    }
    
    # get current costs (per capita)
    if (!isTRUE(parameters$bool_regular_sird)) {
      W_mass <- parameters$zeta * parameters$pop_size + (1 - parameters$zeta) * (Ns + Ni)
      HealthCost <-  HealthCost + fx_per_capita * exp(-parameters$rho * i_day) * parameters$v * new_death
      SocialActivityCost <- SocialActivityCost + fx_per_capita * exp(-parameters$rho * i_day) * W_mass * abs(u_t)
    }
    
    # HealthCost <-  HealthCost+  fx_per_capita * exp(-parameters$rho * i_day) *parameters$v*new_death
    # SocialActivityCost <- SocialActivityCost+ fx_per_capita * exp(-parameters$rho * i_day) * (Ns + Ni) * abs(u_t)
    Rt <- calculate_Rt(parameters$R0, a_t, Ns/parameters$pop_size, Ni)
    
    # Update states
    
    Ns <- Ns + dNs
    Ni <- Ni + dNi
    Nr <- Nr + dNr
    Nd <- Nd + dNd
    Nh <- Nh + dNh
    # keep track of the states
    states_out[i_day+1,] = c(Ns, Ni, Nr, Nd, Nh,
                             HealthCost, SocialActivityCost,
                             HealthCost + SocialActivityCost,
                             a_t, u_t, Rt)
  }
  
  # return as data.frame
  return(data.frame(states_out))
}

# get a vector with the economic summary stats
get_summary_stats <- function(sim_output){
  return(c(get_mean_ci_text(sim_output$HealthCost),
           get_mean_ci_text(sim_output$SocialActivityCost),
           get_mean_ci_text(sim_output$TotalCost)))
}

###Computing 95% confidence and 'quantile' bands'

# get the mean and CI in text format
get_mean_ci_text <- function(vect){
  
  # if only one value, return this value
  if(length(vect)==1){
    return(vect)
  } else{ # else, calculate mean and 95% confidence intervals
    summary_value <- ci_mean(vect)
    summary_text  <- paste0(round(summary_value$estimate),' [',
                            paste(round(summary_value$interval),collapse=','),']')
    return(summary_text)
  }
}

# function to compare stochastic and deterministic output
compare_sim_output <- function(output_experiments, output_deterministic,
                               plot_tag, fadeout_threshold = 0){
  
  # split output
  output_summary <- output_experiments$output_summary
  output_all     <- output_experiments$output_all
  
  # get final states of deterministic model
  output_summary_deterministic <- output_sim_deterministic[nrow(output_sim_deterministic),]
  
  # option to identify simulations without stochastic fade out
  if(fadeout_threshold > 0){
    exp_fade_out   <- (output_summary$Nr+output_summary$Ni) < fadeout_threshold
    output_summary <- output_summary[!exp_fade_out,]
    output_all <- output_all[!exp_fade_out,,]
    plot_tag <- paste(plot_tag,'(excl. fadeout)')
  } else {
    plot_tag <- paste(plot_tag,'(all)')
  }
  
  if (!dir.exists("figures")) dir.create("figures")
  
  pdf(file = paste0("figures/", gsub("[^[:alnum:]]", "_", plot_tag), ".pdf"))
  
  # set figure sub panels
  #par(mfrow=c(3,4))
  
  # explore costs
  output_cost <- output_summary[,grepl('Cost',names(output_summary))]
  y_lim <- range(output_cost)
  
  
  boxplot(output_cost,las=2,ylim=y_lim,main=plot_tag)
  points(colMeans(output_cost),pch=8) # mean
  points((1:3)+0.2,output_summary_deterministic[names(output_cost)],col=4,pch=8,lwd=3) # mean
  
  # explore states
  sel_states <- names(initial_state)[!grepl('Total',names(initial_state)) & !grepl('Ns',names(initial_state))]
  
  # Label mapping for selected variables
  pretty_labels <- c(
    "a_t" = "Activity",
    "Ni" = "Infectives",
    "Ns" = "Susceptibles",
    "Nr" = "Recovered",
    "Nd" = "Deceased",
    "Nh" = "Hospitalized",
    "u_t" = "Utility",
    "Rt" = "Reproduction Number",
    "HealthCost" = "Health Cost",
    "SocialActivityCost" = "Social Activity Cost",
    "TotalCost" = "Total Cost"
  )
  
  for(i_state in sel_states){
    plot_label <- ifelse(i_state %in% names(pretty_labels), pretty_labels[[i_state]], i_state)
    
    # Extract simulation matrix: [experiment, time]
    var_matrix <- output_all[,,i_state]
    
    time_vec <- 0:(ncol(var_matrix)-1)
    #
    
    # Quantile-based band 
    ci_lower <- apply(var_matrix, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
    ci_upper <- apply(var_matrix, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
    
    # Compute mean-based 95% CI band 
    n_sim <- nrow(var_matrix)
    ci_mean <- apply(var_matrix, 2, mean, na.rm = TRUE)
    ci_sd <- apply(var_matrix, 2, sd, na.rm = TRUE)
    ci_se <- ci_sd / sqrt(n_sim)
    ci_lo_meanband <- ci_mean - 1.96 * ci_se
    ci_hi_meanband <- ci_mean + 1.96 * ci_se
    
    # Plot both bands
    plot(NULL, xlim=range(time_vec), ylim=range(ci_lower, ci_upper, ci_lo_meanband, ci_hi_meanband, output_sim_deterministic[,i_state], na.rm=TRUE),
         xlab='time', ylab=plot_label, main=plot_label)
    
    # Plot quantile-based (light blue)
    polygon(c(time_vec, rev(time_vec)), c(ci_upper, rev(ci_lower)), col=alpha("lightblue", 0.4), border=NA)
    
    # Plot mean-based CI (light green)
    polygon(c(time_vec, rev(time_vec)), c(ci_hi_meanband, rev(ci_lo_meanband)), col=alpha("red", 0.4), border=NA)
    
    # Plot mean line
    lines(time_vec, ci_mean, col="blue", lwd=2)
    
    # Plot deterministic path
    lines(output_sim_deterministic[,i_state], col="black", lwd=2, lty=2)
    
  }
  # Compute mean and 95% CI for health cost
  mean_health <- mean(output_summary$HealthCost, na.rm = TRUE)
  ci_health <- ci_mean(output_summary$HealthCost)$interval
  
  # Compute mean and 95% CI for social activity cost
  mean_social <- mean(output_summary$SocialActivityCost, na.rm = TRUE)
  ci_social <- ci_mean(output_summary$SocialActivityCost)$interval
  
  # Print summary (formatted like compare_sim_output)
  cat("Stochastic Summary (from boxplot data):\n")
  cat(sprintf("Health Cost: %.0f [%.0f, %.0f]\n", mean_health, ci_health[1], ci_health[2]))
  cat(sprintf("Social Activity Cost: %.0f [%.0f, %.0f]\n", mean_social, ci_social[1], ci_social[2]))
  
  
  
  ### Credible and confidence interval values calculation ###
  # 1. Health Cost Credible Interval % differences
  health_quant <- quantile(output_summary$HealthCost, probs = c(0.025, 0.975), na.rm = TRUE)
  pct_health_upper = 100 * (health_quant[2] - mean_health) / mean_health
  pct_health_lower = 100 * (mean_health - health_quant[1]) / mean_health
  
  # 2. Activity Cost Credible Interval % differences
  mean_activity_cost <- mean(output_summary$SocialActivityCost, na.rm = TRUE)
  activity_quant <- quantile(output_summary$SocialActivityCost, probs = c(0.025, 0.975), na.rm = TRUE)
  pct_activity_upper = 100 * (activity_quant[2] - mean_activity_cost) / mean_activity_cost
  pct_activity_lower = 100 * (mean_activity_cost - activity_quant[1]) / mean_activity_cost
  
  # 3. Health & Activity cost mean vs deterministic
  health_det <- output_deterministic[nrow(output_deterministic), "HealthCost"]
  activity_det <- output_deterministic[nrow(output_deterministic), "SocialActivityCost"]
  pct_health_mean_diff = 100 * (health_det - mean_health) / health_det
  pct_activity_mean_diff = 100 * (activity_det - mean_activity_cost) / activity_det
  
  
  # Activity & Infectives at Peak Infection #
  
  # Peak infection time (in deterministic)
  t_peak <- which.max(output_deterministic$Ni)
  
  # From stochastic output: activity and Ni across sims at t_peak
  a_vals <- output_all[, t_peak, "a_t"]
  ni_vals <- output_all[, t_peak, "Ni"]
  
  # Means
  a_mean <- mean(a_vals, na.rm = TRUE)
  ni_mean <- mean(ni_vals, na.rm = TRUE)
  
  # Credible intervals
  a_q <- quantile(a_vals, probs = c(0.025, 0.975), na.rm = TRUE)
  ni_q <- quantile(ni_vals, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Deterministic values
  a_det <- output_deterministic[t_peak, "a_t"]
  ni_det <- output_deterministic[t_peak, "Ni"]
  
  # Credible interval width (% from mean)
  pct_a_upper = 100 * (a_q[2] - a_mean) / a_mean
  pct_a_lower = 100 * (a_mean - a_q[1]) / a_mean
  pct_ni_upper = 100 * (ni_q[2] - ni_mean) / ni_mean
  pct_ni_lower = 100 * (ni_mean - ni_q[1]) / ni_mean
  
  # Mean vs deterministic
  pct_a_mean_diff = 100 * (a_det - a_mean) / a_det
  pct_ni_mean_diff = 100 * (ni_det - ni_mean) / ni_det
  
  # ==== OUTPUT ====
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
  
  
  # Health Cost boxplot with mean, 95% CI, and deterministic point
  boxplot(output_summary$HealthCost,
          main = "Health Cost",
          ylab = "Health Cost",
          ylim = range(c(output_summary$HealthCost, output_summary_deterministic$HealthCost), na.rm = TRUE),
          col = "gray90")
  
  # Add stochastic mean
  points(1, mean(output_summary$HealthCost, na.rm = TRUE), pch = 8)
  
  # Add 95% confidence interval
  
  ci_health <- ci_mean(output_summary$HealthCost)
  
  segments(x0 = 1, x1 = 1, y0 = ci_health$interval[1], y1 = ci_health$interval[2], col = "forestgreen", lwd = 2)
  
  # Add 95% quantile interval (e.g., 2.5% - 97.5%)
  quant_health <- quantile(output_summary$HealthCost, probs = c(0.025, 0.975), na.rm = TRUE)
  segments(x0 = 1.1, x1 = 1.1, y0 = quant_health[1], y1 = quant_health[2], col = "blue", lwd = 2, lty = 2)
  
  
  # Add deterministic outcome
  points(1.2, output_summary_deterministic$HealthCost, col = "red", pch = 18, cex = 1.5)
  
  
  # Social Activity Cost boxplot with mean, 95% CI, and deterministic point
  boxplot(output_summary$SocialActivityCost,
          main = "Social Activity Cost",
          ylab = "Social Activity Cost",
          ylim = range(c(output_summary$SocialActivityCost, output_summary_deterministic$SocialActivityCost), na.rm = TRUE),
          col = "gray90")
  
  # Add stochastic mean
  points(1, mean(output_summary$SocialActivityCost, na.rm = TRUE), pch = 8)
  
  # Add 95% confidence interval
  ci_social <- ci_mean(output_summary$SocialActivityCost)
  segments(x0 = 1, x1 = 1, y0 = ci_social$interval[1], y1 = ci_social$interval[2], col = "forestgreen", lwd = 2)
  
  quant_social <- quantile(output_summary$SocialActivityCost, probs = c(0.025, 0.975), na.rm = TRUE)
  segments(x0 = 1.1, x1 = 1.1, y0 = quant_social[1], y1 = quant_social[2], col = "blue", lwd = 2, lty = 2)
  
  
  # Add deterministic outcome
  points(1.2, output_summary_deterministic$SocialActivityCost, col = "red", pch = 18, cex = 1.5)
  
  # prepare summary statistics
  print_out <- data.frame(output=c('Health Cost (per capita)','Social Activity Cost (per capita)','Total Cost (per capita)'),
                          deterministic = get_summary_stats(output_sim_deterministic[nrow(output_sim_deterministic),]),
                          stochastic = get_summary_stats(output_summary))
  
  names(print_out)[1] <- plot_tag
  
  # print summary statistics to console
  print(paste('SUMMARY STATISTICS:',plot_tag))
  print(print_out)
  # Peak Infection Timing Summary 
  peak_times <- output_summary$PeakTime
  
  # Compute mean ± 1.96 × SE using ci_mean()
  peak_ci <- ci_mean(peak_times)
  
  # Print mean CI
  cat("\nPeak Infection Timing Summary:\n")
  cat(sprintf("Mean Peak Time: %.1f [%.1f, %.1f]  (mean ± 1.96 * SE)\n",
              peak_ci$estimate, peak_ci$interval[1], peak_ci$interval[2]))
  
  # Print quantiles
  cat(sprintf("Quantiles: 2.5%% = %.1f, 50%% (median) = %.1f, 97.5%% = %.1f\n",
              quantile(peak_times, 0.025, na.rm = TRUE),
              quantile(peak_times, 0.50, na.rm = TRUE),
              quantile(peak_times, 0.975, na.rm = TRUE)))
  dev.off()
}

run_experiments <- function(initial_state, times, parameters,
                            bool_stochastic_beta, update_function,
                            num_experiments)
{
  
  # Set random number generator seed for reproducibility
  set.seed(parameters$rng_seed)
  
  # Run a single simulation to get the structure of the output
  temp_output <- run_sir_binomial(initial_state, times, parameters, bool_stochastic_beta, update_function)
  
  # Create output_summary with an extra column for PeakTime
  output_summary <- data.frame(matrix(NA, nrow = num_experiments, ncol = ncol(temp_output) + 1))
  names(output_summary) <- c(names(temp_output), "PeakTime")
  
  # Create 3D array to store full time-series output for all simulations
  output_all <- array(NA, dim = c(num_experiments, length(times), ncol(temp_output)),
                      dimnames = list(NULL, NULL, names(temp_output)))
  
  # Run simulations
  for(i_exp in 1:num_experiments){
    print(i_exp)
    
    # Run the simulation
    output_sim <- run_sir_binomial(initial_state = initial_state,
                                   times = times,
                                   parameters = parameters,
                                   bool_stochastic_beta = bool_stochastic_beta,
                                   update_function = update_function)
    
    # Store final row of simulation in output_summary
    output_summary[i_exp, 1:ncol(temp_output)] <- output_sim[nrow(output_sim), ]
    
    # Store PeakTime (subtract 1 to match time indexing)
    output_summary$PeakTime[i_exp] <- which.max(output_sim$Ni) - 1
    
    # Store full simulation output
    output_all[i_exp, , ] <- as.matrix(output_sim)
  }
  return(list(output_summary = output_summary,
              output_all = output_all))
}