
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

a_function <- function(Ni, Ns, parameters, i_day = NULL) {
  if (isTRUE(parameters$bool_regular_sird)) return(1)
  
  # proportions
  pop <- parameters$pop_size
  Ni_prop <- Ni / pop
  Ns_prop <- Ns / pop
  if (Ni_prop < 1e-10 || Ns_prop < 1e-10) return(1)
  
  # Info campaign: inflate perceived risk during campaign window
  in_campaign <- !is.null(i_day) &&
    i_day >= parameters$campaign_start &&
    i_day <= parameters$campaign_end
  risk_mult <- if (in_campaign) parameters$risk_multiplier else 1.0
  
  # Pigouvian tax (baseline + prevalence-linked component)
  tau_t <- parameters$tau0 + parameters$tau_per_Ni * Ni_prop
  
  # Original closed-form structure, augmented:
  # (i) risk term scaled by risk_mult
  # (ii) linear tax term folded in as an added penalty ~ tau_t * (Ns+Ni)
  base_mult <- parameters$beta * parameters$pi * parameters$v * Ns_prop * Ni_prop
  multiplier <- risk_mult * base_mult + tau_t * (Ns_prop + Ni_prop)
  
  num_sum <- Ns_prop + Ni_prop
  discrim <- num_sum^2 + 4 * multiplier * num_sum
  denom <- 2 * multiplier
  
  if (denom <= 0 || discrim < 0) return(1)
  a_t <- (-num_sum + sqrt(discrim)) / denom
  
  return(max(0, min(1, a_t)))
}


# a_function_cost_minimizing <- function(Ns, Ni, t, parameters) {
#   pop_size <- parameters$pop_size
#   Ns_prop <- Ns / pop_size
#   Ni_prop <- Ni / pop_size
#   
#   if (Ni <= 1e-12 || Ns <= 1e-12) return(1)
#   
#   beta <- parameters$beta
#   gamma <- parameters$gamma
#   pi <- parameters$pi
#   v <- parameters$v
#   rho <- parameters$rho
#   
#   fx_per_capita <- parameters$fx / pop_size
#   discount <- exp(-rho * t)
#   
#   utility_function <- function(a) {
#     if (parameters$utility_type == "Log") return(log(a) - a + 1)
#     else return(-0.5 * (1 - a)^2)
#   }
#   
#   total_cost_fn <- function(a) {
#     # Compute probabilities
#     p_infect <- 1 - exp(- beta * a^2 * Ni / pop_size)
#     p_recover <- 1 - exp(-gamma)
#     p_death <- 1 - exp(-pi)
#     
#     # Expected transitions
#     E_deaths <- Ns * p_infect * p_recover * p_death
#     
#     # Cost components
#     HC <- fx_per_capita * discount * v * E_deaths
#     SAC <- fx_per_capita * discount * abs(utility_function(a)) * (Ns + Ni)
#     
#     return(HC + SAC)
#   }
#   
#   res <- optimize(total_cost_fn, lower = 0.01, upper = 1, tol = 1e-4)
#   return(res$minimum)
# }


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
beta_with_policies <- function(parameters, beta_t) {
  # Masks/distancing
  beta_masked <- beta_t * (1 - parameters$mask_effect * parameters$mask_compliance)
  
  # Sectoral override (if any sector betas are set > 0, use their sum)
  sector_sum <- parameters$beta_home + parameters$beta_work +
    (if (isTRUE(parameters$school_open)) parameters$beta_school else 0) +
    parameters$beta_other
  
  if (sector_sum > 0) return(sector_sum)
  return(beta_masked)
}

run_sir_binomial <- function(initial_state,
                             times,
                             parameters,
                             bool_stochastic_beta = FALSE,
                             update_function = get_transitions_stochastic){
  
  # copy & convert states
  states <- data.frame(t(initial_state))
  states[grepl('N.', names(states))] <- round(states[grepl('N.', names(states))] * parameters$pop_size)
  states[grepl('Ns', names(states))] <- parameters$pop_size -
    sum(states[grepl('N.', names(states)) & !grepl('Ns', names(states))])
  
  states_out <- matrix(NA, nrow = length(times), ncol = length(states))
  colnames(states_out) <- names(states)
  
  # init scalars
  Ns <- states$Ns; Ni <- states$Ni; Nr <- states$Nr; Nd <- states$Nd
  HealthCost <- 0; SocialActivityCost <- 0
  fx_per_capita <- parameters$fx / parameters$pop_size
  
  # lockdown state with hysteresis (activity cap)
  lock_active <- FALSE
  
  for (i_day in times[-1]) {
    
    # 0) stochastic β
    beta_t <- if (bool_stochastic_beta) {
      max(1e-4, rnorm(1, mean = parameters$beta, sd = parameters$sigma))
    } else parameters$beta
    
    # 1) masks + sectoral NPIs
    beta_policy <- beta_with_policies(parameters, beta_t)
    
    # 2) choose activity a_t (with Pigouvian tax + info campaign)
    a_t <- if (isTRUE(parameters$bool_daily_cost_minimizing)) {
      a_function_cost_minimizing(Ns, Ni, i_day, parameters)  # (if you later re-enable)
    } else {
      a_function(Ni, Ns, parameters, i_day)
    }
    
    # 3) lockdown / activity cap trigger logic
    if (!is.na(parameters$lock_trigger_Ni) && Ni >= parameters$lock_trigger_Ni) lock_active <- TRUE
    if (!is.na(parameters$lock_release_Ni) && Ni <= parameters$lock_release_Ni) lock_active <- FALSE
    if (isTRUE(lock_active)) a_t <- min(a_t, parameters$a_cap)
    
    # 4) probabilities with current policies
    p_infect  <- 1 - exp(- beta_policy * a_t^2 * Ni / parameters$pop_size)
    
    # testing/isolation increases effective removal
    gamma_eff <- parameters$gamma + parameters$test_rate * parameters$isolation_effect
    p_recover <- 1 - exp(- gamma_eff)
    
    # hospital capacity: smooth excess IFR
    if (!is.null(parameters$healthcare_capacity) && parameters$healthcare_capacity > 0 &&
        parameters$excess_mortality_multiplier > 1 && parameters$capacity_slope > 0) {
      overload <- max(0, (Ni - parameters$healthcare_capacity)) / max(1, parameters$healthcare_capacity)
      mult <- 1 + (parameters$excess_mortality_multiplier - 1) * (1 - exp(-parameters$capacity_slope * overload))
      pi_effective <- parameters$pi * mult
    } else if (!is.null(parameters$healthcare_capacity) && Ni > parameters$healthcare_capacity &&
               parameters$excess_mortality_multiplier > 1) {
      pi_effective <- parameters$pi * parameters$excess_mortality_multiplier
    } else {
      pi_effective <- parameters$pi
    }
    p_death <- 1 - exp(- pi_effective)
    
    # 5) vaccination (S -> R) before infections that day
    new_vax <- 0
    if (i_day >= parameters$start_vax_day && parameters$nu > 0) {
      p_vax <- 1 - exp(-parameters$nu)
      new_vax <- update_function(Ns, prob = p_vax)
    }
    Ns_after_vax <- max(0, Ns - new_vax)
    
    # 6) infections, recoveries, deaths
    new_infections <- update_function(Ns_after_vax, prob = p_infect)
    new_recoveries <- update_function(Ni, prob = p_recover)
    new_death      <- update_function(new_recoveries, prob = p_death)
    
    # 7) imported cases (add to infections; still limited by remaining S)
    imports <- if (parameters$lambda_import > 0) rpois(1, parameters$lambda_import) else 0
    can_infect <- max(0, Ns_after_vax - new_infections)
    add_imports <- min(imports, can_infect)
    new_infections <- new_infections + add_imports
    
    # fadeout safeguard (your existing logic)
    if ((Ni - new_recoveries) < parameters$infect_thres) {
      new_recoveries <- Ni
      new_infections <- max(0, new_infections - min(new_infections, new_recoveries))  # keep S nonnegative
    }
    
    # 8) state deltas
    dNs <- -(new_vax + new_infections)
    dNi <- new_infections - new_recoveries
    dNr <- new_recoveries - new_death + new_vax   # vaccinated go to R directly
    dNd <- new_death
    
    # 9) costs
    u_t <- utility_function(a_t)
    if (!isTRUE(parameters$bool_regular_sird)) {
      HealthCost        <- HealthCost        + fx_per_capita * exp(-parameters$rho * i_day) * parameters$v * new_death
      SocialActivityCost<- SocialActivityCost+ fx_per_capita * exp(-parameters$rho * i_day) * (Ns + Ni) * abs(u_t)
    }
    
    # 10) update states
    Ns <- Ns + dNs; Ni <- Ni + dNi; Nr <- Nr + dNr; Nd <- Nd + dNd
    Rt <- calculate_Rt(parameters$R0, a_t, Ns/parameters$pop_size, Ni)
    
    states_out[i_day + 1, ] <- c(Ns, Ni, Nr, Nd,
                                 HealthCost, SocialActivityCost,
                                 HealthCost + SocialActivityCost,
                                 a_t, u_t, Rt)
  }
  
  data.frame(states_out)
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
  output_summary_deterministic <- output_deterministic[nrow(output_deterministic),]
  
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
    plot(NULL, xlim=range(time_vec), ylim=range(ci_lower, ci_upper, ci_lo_meanband, ci_hi_meanband, output_deterministic[,i_state], na.rm=TRUE),
         xlab='time', ylab=plot_label, main=plot_label)
    
    # Plot quantile-based (light blue)
    polygon(c(time_vec, rev(time_vec)), c(ci_upper, rev(ci_lower)), col=alpha("lightblue", 0.4), border=NA)
    
    # Plot mean-based CI (light green)
    polygon(c(time_vec, rev(time_vec)), c(ci_hi_meanband, rev(ci_lo_meanband)), col=alpha("red", 0.4), border=NA)
    
    # Plot mean line
    lines(time_vec, ci_mean, col="blue", lwd=2)
    
    # Plot deterministic path
    lines(output_deterministic[,i_state], col="black", lwd=2, lty=2)
    
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
                          deterministic = get_summary_stats(output_deterministic[nrow(output_deterministic),]),
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