
##
# Help function for the EPI-ECON modelling work
#
# - Deterministic and stochastic transitions
################################################################ #
#setwd("C:/Users/Bart Smets/OneDrive/Documenten/GitHub/Working_Code_BS/R")
library(confintr)
library(scales)

a_function <- function(Ni, Ns, parameters) {
  if (!is.null(parameters$bool_regular_sird) && isTRUE(parameters$bool_regular_sird)) {
    return(1)
  }
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

# Utility function
utility_function <- function(a_t, parameters) {
  if (!is.null(parameters$bool_regular_sird) && isTRUE(parameters$bool_regular_sird)) {
    return(0)
  }
  return(log(a_t) - a_t + 1)
}

#Calculation of effective reproduction number
calculate_Rt <- function(R0, a_t, Ns_prop, Ni) {
  if (Ni <= 0) return(0)
  return(R0 * a_t^2 * Ns_prop)
}

get_transitions_stochastic <- function(n, prob) {
  
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




run_sir_binomial <- function(initial_state,
                             times,
                             parameters,
                             update_function = get_transitions_stochastic) {
  # copy initial states
  states <- data.frame(t(initial_state))

  # convert population fractions into numbers
  states[grepl('N.', names(states))] <- round(states[grepl('N.', names(states))] * parameters$pop_size)

  # keep total population fixed (adjust Ns if needed)
    states[grepl('Ns', names(states))] <- parameters$pop_size -
    sum(states[grepl('N.', names(states))] & !grepl('Ns', names(states)))


  # output matrix
  states_out <- matrix(NA, nrow = length(times), ncol = length(states))
  colnames(states_out) <- names(states)

  # init stocks and costs
  Ns <- states$Ns; Ni <- states$Ni; Nr <- states$Nr; Nd <- states$Nd
  HealthCost <- 0; SocialActivityCost <- 0

  fx_per_capita <- parameters$fx / parameters$pop_size
  regular_mode  <- !is.null(parameters$bool_regular_sird) && isTRUE(parameters$bool_regular_sird)

  for (i_day in times[-1]) {
    beta_t  <- parameters$beta
    Ns_prop <- Ns / parameters$pop_size

    if (regular_mode) {
      # ----- Regular SIRD: fixed activity, no activity cost -----
      a_t <- 1
      u_t <- 0

      # NOTE: no a_t^2 here in regular mode
      p_infect  <- 1 - exp(- beta_t * (Ni / parameters$pop_size))
      p_recover <- 1 - exp(- (1 - parameters$pi) * parameters$gamma)
      p_death   <- 1 - exp(- parameters$pi * parameters$gamma)

      new_infections <- update_function(Ns, prob = p_infect)
      new_recoveries <- update_function(Ni, prob = p_recover)
      new_death      <- update_function(Ni, prob = p_death)

      # fade-out short-circuit
      if ((Ni - new_recoveries) < parameters$infect_thres) {
        new_recoveries <- Ni
        new_infections <- 0
      }

      dNs <- -new_infections
      dNi <-  new_infections - new_recoveries - new_death
      dNr <-  new_recoveries
      dNd <-  new_death

      # costs: health only
      HealthCost <- HealthCost + fx_per_capita * exp(-parameters$rho * i_day) * parameters$v * new_death

      Rt <- calculate_Rt(parameters$R0, a_t, Ns_prop, Ni)

    } else {
      # ----- Behavioral epi-econ branch -----
      a_t <- a_function(Ni, Ns, parameters)
      u_t <- utility_function(a_t, parameters)

      p_infect  <- 1 - exp(- beta_t * a_t^2 * (Ni / parameters$pop_size))
      p_recover <- 1 - exp(- (1 - parameters$pi) * parameters$gamma)
      p_death   <- 1 - exp(- parameters$pi * parameters$gamma)

      new_infections <- update_function(Ns, prob = p_infect)
      new_recoveries <- update_function(Ni, prob = p_recover)
      new_death      <- update_function(Ni, prob = p_death)

      # fade-out short-circuit
      if ((Ni - new_recoveries) < parameters$infect_thres) {
        new_recoveries <- Ni
        new_infections <- 0
      }

      dNs <- -new_infections
      dNi <-  new_infections - new_recoveries - new_death
      dNr <-  new_recoveries
      dNd <-  new_death

      HealthCost <- HealthCost +
        fx_per_capita * exp(-parameters$rho * i_day) * parameters$v * new_death

      SocialActivityCost <- SocialActivityCost +
        fx_per_capita * exp(-parameters$rho * i_day) * (Ns + Ni) * abs(u_t)

      Rt <- calculate_Rt(parameters$R0, a_t, Ns_prop, Ni)
    }

    # update stocks
    Ns <- Ns + dNs; Ni <- Ni + dNi; Nr <- Nr + dNr; Nd <- Nd + dNd

    # record
    states_out[i_day + 1, ] <- c(Ns, Ni, Nr, Nd,
                                 HealthCost, SocialActivityCost,
                                 HealthCost + SocialActivityCost,
                                 a_t, u_t, Rt)
  }

  data.frame(states_out)
}
#run_sir_binomial function with carry/accumulator mechanism for discrete integer transitions, comment out in normal mode.
# run_sir_binomial <- function(initial_state,
#                              times,
#                              parameters,
#                              update_function = get_transitions_stochastic){
# 
#   # copy initial states
#   states <- data.frame(t(initial_state))
# 
#   # convert population fractions into numbers
#   states[grepl('N.',names(states))] <- round(states[grepl('N.',names(states))] * parameters$pop_size)
# 
#   # keep total population fixed (adjust Ns if needed)
#   states[grepl('Ns',names(states))] <- parameters$pop_size -
#     sum(states[grepl('N.',names(states)) & !grepl('Ns',names(states))])
# 
#   # output matrix (your code assumes these names exist in initial_state; keeping same behavior)
#   states_out <- matrix(NA, nrow = length(times), ncol = length(states))
#   colnames(states_out) <- names(states)
# 
#   # init stocks and costs
#   Ns <- states$Ns; Ni <- states$Ni; Nr <- states$Nr; Nd <- states$Nd
#   HealthCost <- 0; SocialActivityCost <- 0
# 
#   fx_per_capita <- parameters$fx / parameters$pop_size
#   regular_mode  <- !is.null(parameters$bool_regular_sird) && isTRUE(parameters$bool_regular_sird)
# 
#   # --- NEW: carry/accumulator switch & helper ---
#  # use_integer_with_carry <- isTRUE(parameters$integer_with_carry)
#   # BEFORE
#   use_integer_with_carry <- isTRUE(parameters$integer_with_carry)
#   
#   # AFTER  (only carry when deterministic updater is selected)
#   use_integer_with_carry <- isTRUE(parameters$integer_with_carry) &&
#     identical(update_function, get_transitions_deterministic)
#   
#   # local carry buffers (only used when use_integer_with_carry == TRUE)
#   carry <- list(inf = 0, rec = 0, dea = 0)
# 
#   take_with_carry <- function(expected, key, cap = Inf) {
#     # Only used in deterministic integer path; guards + caps to avoid negatives / overflows
#     if (!is.finite(expected) || expected < 0) expected <- 0
#     carry[[key]] <<- carry[[key]] + expected
#     realized <- floor(carry[[key]])
#     if (realized > cap) realized <- cap
#     if (realized < 0)  realized <- 0
#     carry[[key]] <<- carry[[key]] - realized
#     return(as.numeric(realized))
#   }
# 
#   # record initial row (time 0)
#   states_out[1,] <- c(Ns, Ni, Nr, Nd,
#                       HealthCost, SocialActivityCost,
#                       HealthCost + SocialActivityCost,
#                       a_t = if (regular_mode) 1 else a_function(Ni, Ns, parameters),
#                       u_t = if (regular_mode) 0 else utility_function(if (regular_mode) 1 else a_function(Ni, Ns, parameters), parameters),
#                       Rt  = calculate_Rt(parameters$R0, if (regular_mode) 1 else a_function(Ni, Ns, parameters), Ns/parameters$pop_size, Ni))
# 
#   for(i_day in times[-1]){
#     beta_t   <- parameters$beta
#     Ns_prop  <- Ns / parameters$pop_size
# 
#     if (regular_mode) {
#       # ----- Regular SIRD: fixed activity, no activity cost -----
#       a_t <- 1
#       u_t <- 0
# 
#       # NOTE: no a_t^2 here in regular mode
#       p_infect  <- 1 - exp(- beta_t * (Ni / parameters$pop_size))
#       p_recover <- 1 - exp(- (1 - parameters$pi) * parameters$gamma)
#       p_death   <- 1 - exp(- parameters$pi * parameters$gamma)
# 
#       if (use_integer_with_carry) {
#         # expected flows
#         exp_infections <- Ns * p_infect
#         exp_recoveries <- Ni * p_recover
#         exp_deaths     <- Ni * p_death
# 
#         # realize integers with carry (order: infections, deaths, recoveries)
#         new_infections <- take_with_carry(exp_infections, "inf", cap = Ns)
#         new_death      <- take_with_carry(exp_deaths,     "dea", cap = Ni)
#         new_recoveries <- take_with_carry(exp_recoveries, "rec", cap = max(Ni - new_death, 0))
# 
#       } else {
#         # original path: delegate to update_function (stoch or deterministic)
#         new_infections <- update_function(Ns, prob = p_infect)
#         new_recoveries <- update_function(Ni, prob = p_recover)
#         new_death      <- update_function(Ni, prob = p_death)
#       }
# 
#       # fade-out short-circuit (apply after computing flows)
#       if ((Ni - new_recoveries) < parameters$infect_thres) {
#         new_recoveries <- Ni
#         new_infections <- 0
#         # optional: clear carries so they don't drip post-fadeout
#         carry$rec <- 0; carry$inf <- 0; carry$dea <- 0
#       }
# 
#       dNs <- -new_infections
#       dNi <-  new_infections - new_recoveries - new_death
#       dNr <-  new_recoveries
#       dNd <-  new_death
# 
#       # costs: choose expected or realized
#       if (identical(parameters$costs_from, "realized")) {
#         HealthCost <- HealthCost + fx_per_capita * exp(-parameters$rho * i_day) * parameters$v * new_death
#       } else {
#         # default: expected (smooth, avoids "collapse")
#         exp_deaths_today <- Ni * p_death
#         HealthCost <- HealthCost + fx_per_capita * exp(-parameters$rho * i_day) * parameters$v * exp_deaths_today
#       }
# 
#       Rt <- calculate_Rt(parameters$R0, a_t, Ns_prop, Ni)
# 
#     } else {
#       # ----- Behavioral epi-econ branch -----
#       a_t <- a_function(Ni, Ns, parameters)
#       u_t <- utility_function(a_t, parameters)
# 
#       p_infect  <- 1 - exp(- beta_t * a_t^2 * (Ni / parameters$pop_size))
#       p_recover <- 1 - exp(- (1 - parameters$pi) * parameters$gamma)
#       p_death   <- 1 - exp(- parameters$pi * parameters$gamma)
# 
#       if (use_integer_with_carry) {
#         # expected flows
#         exp_infections <- Ns * p_infect
#         exp_recoveries <- Ni * p_recover
#         exp_deaths     <- Ni * p_death
# 
#         # realize integers with carry (order: infections, deaths, recoveries)
#         new_infections <- take_with_carry(exp_infections, "inf", cap = Ns)
#         new_death      <- take_with_carry(exp_deaths,     "dea", cap = Ni)
#         new_recoveries <- take_with_carry(exp_recoveries, "rec", cap = max(Ni - new_death, 0))
# 
#       } else {
#         # original path: delegate to update_function (stoch or deterministic)
#         new_infections <- update_function(Ns, prob = p_infect)
#         new_recoveries <- update_function(Ni, prob = p_recover)
#         new_death      <- update_function(Ni, prob = p_death)
#       }
# 
#       # fade-out short-circuit
#       if ((Ni - new_recoveries) < parameters$infect_thres) {
#         new_recoveries <- Ni
#         new_infections <- 0
#         carry$rec <- 0; carry$inf <- 0; carry$dea <- 0
#       }
# 
#       dNs <- -new_infections
#       dNi <-  new_infections - new_recoveries - new_death
#       dNr <-  new_recoveries
#       dNd <-  new_death
# 
#       # Health cost: expected vs realized
#       if (identical(parameters$costs_from, "realized")) {
#         HealthCost <-  HealthCost + fx_per_capita * exp(-parameters$rho * i_day) * parameters$v * new_death
#       } else {
#         exp_deaths_today <- Ni * p_death
#         HealthCost <-  HealthCost + fx_per_capita * exp(-parameters$rho * i_day) * parameters$v * exp_deaths_today
#       }
# 
#       SocialActivityCost <- SocialActivityCost +
#         fx_per_capita * exp(-parameters$rho * i_day) * (Ns + Ni) * abs(u_t)
# 
#       Rt <- calculate_Rt(parameters$R0, a_t, Ns_prop, Ni)
#     }
# 
#     # update stocks
#     Ns <- Ns + dNs; Ni <- Ni + dNi; Nr <- Nr + dNr; Nd <- Nd + dNd
# 
#     # record
#     states_out[i_day+1,] <- c(Ns, Ni, Nr, Nd,
#                               HealthCost, SocialActivityCost,
#                               HealthCost + SocialActivityCost,
#                               a_t, u_t, Rt)
#   }
# 
#   data.frame(states_out)
# }

# get a vector with the economic summary stats
get_summary_stats <- function(sim_output){
  return(c(get_mean_ci_text(sim_output$HealthCost),
           get_mean_ci_text(sim_output$SocialActivityCost),
           get_mean_ci_text(sim_output$TotalCost)))
}

# get the mean and CI in text format
get_mean_ci_text <- function(vect){
  if(length(vect)==1){
    return(vect)
  } else{
    summary_value <- ci_mean(vect)
    summary_text  <- paste0(round(summary_value$estimate),' [',
                            paste(round(summary_value$interval),collapse=','),']')
    return(summary_text)
  }
}

# function to compare stochastic and deterministic output
compare_sim_output <- function(output_experiments, output_deterministic,
                               plot_tag, fadeout_threshold = 100){
  
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
  
  # explore costs
  output_cost <- output_summary[,grepl('Cost',names(output_summary))]
  y_lim <- range(output_cost)
  
  
  boxplot(output_cost,las=2,ylim=y_lim,main=plot_tag)
  points(colMeans(output_cost),pch=8) # mean
  points((1:3)+0.2,output_summary_deterministic[names(output_cost)],col=4,pch=8,lwd=3) # mean
  
  # explore states
  sel_states <- names(initial_state)[!grepl('Total',names(initial_state)) & !grepl('Ns',names(initial_state))]
  
  # Label mapping for selected variables
  labels <- c(
    "a_t" = "Activity",
    "Ni" = "Infectives",
    "Ns" = "Susceptibles",
    "Nr" = "Recovered",
    "Nd" = "Deceased",
    "u_t" = "Utility",
    "Rt" = "Reproduction Number",
    "HealthCost" = "Health Cost",
    "SocialActivityCost" = "Activity Cost",
    "TotalCost" = "Total Cost"
  )
  
  for(i_state in sel_states){
    plot_label <- ifelse(i_state %in% names(labels), labels[[i_state]], i_state)
    
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
    
    # choose y-lims that include the lower bands (can be negative for Utility)
    y_min <- min(ci_lower, ci_lo_meanband, output_deterministic[, i_state], na.rm = TRUE)
    y_max <- max(ci_upper, ci_hi_meanband, output_deterministic[, i_state], na.rm = TRUE)
    
    plot(NULL,
         xlim = c(min(time_vec), max(time_vec)),
         ylim = c(y_min, y_max),
         xlab = "time", ylab = plot_label, main = plot_label)
    abline(h = 0, lty = 3, col = "grey60")
    
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
  
  
  
  ### Quantile and confidence interval values calculation ###
  # 1. Health Cost Quantile Interval % differences
  health_quant <- quantile(output_summary$HealthCost, probs = c(0.025, 0.975), na.rm = TRUE)
  pct_health_upper = 100 * (health_quant[2] - mean_health) / mean_health
  pct_health_lower = 100 * (mean_health - health_quant[1]) / mean_health
  
  # 2. Activity Cost Quantile Interval % differences
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
  
  # Quantile intervals
  a_q <- quantile(a_vals, probs = c(0.025, 0.975), na.rm = TRUE)
  ni_q <- quantile(ni_vals, probs = c(0.025, 0.975), na.rm = TRUE)
  
  # Deterministic values
  a_det <- output_deterministic[t_peak, "a_t"]
  ni_det <- output_deterministic[t_peak, "Ni"]
  
  # Quantile interval width (% from mean)
  pct_a_upper = 100 * (a_q[2] - a_mean) / a_mean
  pct_a_lower = 100 * (a_mean - a_q[1]) / a_mean
  pct_ni_upper = 100 * (ni_q[2] - ni_mean) / ni_mean
  pct_ni_lower = 100 * (ni_mean - ni_q[1]) / ni_mean
  
  # Mean vs deterministic
  pct_a_mean_diff = 100 * (a_det - a_mean) / a_det
  pct_ni_mean_diff = 100 * (ni_det - ni_mean) / ni_det
  
  # ==== OUTPUT ====
  cat("\nCOST DIFFERENCES (BOXPLOT METRICS)\n")
  cat(sprintf("Health Cost – Quantile Interval: +%.2f%% / -%.2f%%\n", pct_health_upper, pct_health_lower))
  cat(sprintf("Activity Cost – Quantile Interval: +%.2f%% / -%.2f%%\n", pct_activity_upper, pct_activity_lower))
  cat(sprintf("Health Cost – Stochastic vs Det: %.2f%%\n", pct_health_mean_diff))
  cat(sprintf("Activity Cost – Stochastic vs Det: %.2f%%\n", pct_activity_mean_diff))
  
  cat("\nPEAK-TIME DIFFERENCES (Ni and a_t at peak Ni)\n")
  cat(sprintf("Activity a(t) – Quantile Interval: +%.2f%% / -%.2f%%\n", pct_a_upper, pct_a_lower))
  cat(sprintf("Infectives Ni – Quantile Interval: +%.2f%% / -%.2f%%\n", pct_ni_upper, pct_ni_lower))
  cat(sprintf("Activity a(t) – Stochastic vs Det: %.2f%%\n", pct_a_mean_diff))
  cat(sprintf("Infectives Ni – Stochastic vs Det: %.2f%%\n", pct_ni_mean_diff))
  
  
  # Health Cost boxplot with mean, 95% CI, and deterministic point
  boxplot(output_summary$HealthCost,
          main = "Health Cost",
          ylab = "Health Cost",
          ylim = range(c(output_summary$HealthCost, output_summary_deterministic$HealthCost), na.rm = TRUE),
          col = "gray90")
  
  # Add stochastic mean
 # points(1, mean(output_summary$HealthCost, na.rm = TRUE), pch = 8)
  
  # Add 95% confidence interval
  
  ci_health <- ci_mean(output_summary$HealthCost)
  
  segments(x0 = 1, x1 = 1, y0 = ci_health$interval[1], y1 = ci_health$interval[2], col = "forestgreen", lwd = 2)
  
  # Add 95% quantile interval (e.g., 2.5% - 97.5%)
  quant_health <- quantile(output_summary$HealthCost, probs = c(0.025, 0.975), na.rm = TRUE)
  segments(x0 = 1.1, x1 = 1.1, y0 = quant_health[1], y1 = quant_health[2], col = "blue", lwd = 2, lty = 2)
  
  
  
  # Deterministic and stochastic outcome lines for healthcost
  abline(h = output_summary_deterministic$HealthCost, col = "red", lwd = 2)
  abline(h = mean(output_summary$HealthCost, na.rm = TRUE), col = "purple", lwd = 2)
  
  # Social Activity Cost boxplot with mean, 95% CI, and deterministic point
  boxplot(output_summary$SocialActivityCost,
          main = "Activity Cost",
          ylab = "Activity Cost",
          ylim = range(c(output_summary$SocialActivityCost, output_summary_deterministic$SocialActivityCost), na.rm = TRUE),
          col = "gray90")
  
  # Add these two lines 
  abline(h = mean(output_summary$SocialActivityCost, na.rm = TRUE), col = "purple", lwd = 2)
  abline(h = output_summary_deterministic$SocialActivityCost, col = "red", lwd = 2)
  

  # Add 95% confidence interval
  ci_social <- ci_mean(output_summary$SocialActivityCost)
  segments(x0 = 1, x1 = 1, y0 = ci_social$interval[1], y1 = ci_social$interval[2], col = "forestgreen", lwd = 2)
  
  quant_social <- quantile(output_summary$SocialActivityCost, probs = c(0.025, 0.975), na.rm = TRUE)
  segments(x0 = 1.1, x1 = 1.1, y0 = quant_social[1], y1 = quant_social[2], col = "blue", lwd = 2, lty = 2)
  
  
  # prepare summary statistics
  print_out <- data.frame(output=c('Health Cost (per capita)','Activity Cost (per capita)','Total Cost (per capita)'),
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

run_experiments <- function(initial_state, times, parameters, update_function,
                            num_experiments){
  
  # Set random number generator seed for reproducibility
  set.seed(parameters$rng_seed)
  
  # Run a single simulation to get the structure of the output
  temp_output <- run_sir_binomial(initial_state, times, parameters, update_function)
  
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
