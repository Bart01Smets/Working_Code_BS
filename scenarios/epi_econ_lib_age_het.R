
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

# --- NEW: tiny helpers ---

# Split an integer 'total' across groups using shares (robust to rounding)
split_counts <- function(total, shares) {
  shares <- shares / sum(shares)
  v <- round(total * shares)
  diff <- total - sum(v)
  if (diff != 0) v[length(v)] <- v[length(v)] + diff
  v
}

# Compute a by-age activity vector using your existing a_function
# by temporarily plugging in pi/v for each group.
a_function_by_age <- function(Ni_vec, Ns_vec, parameters) {
  G <- length(parameters$age_groups)
  a_vec <- numeric(G)
  for (g in seq_len(G)) {
    tmp <- parameters
    if (!is.null(parameters$pi_by_age)) tmp$pi <- parameters$pi_by_age[g]
    if (!is.null(parameters$v_by_age))  tmp$v  <- parameters$v_by_age[g]
    a_vec[g] <- a_function(Ni_vec[g], Ns_vec[g], tmp)
  }
  pmax(0, pmin(1, a_vec))
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
  Ns <- states$Ns; Ni <- states$Ni; Nr <- states$Nr; Nd <- states$Nd
  HealthCost <- 0; SocialActivityCost <- 0
  
  fx_per_capita <- parameters$fx / parameters$pop_size
  
  # --- NEW: initialize by-age state if heterogeneity is ON ---
  if (isTRUE(parameters$bool_age_heterogeneity)) {
    G <- length(parameters$age_groups)
    shares <- parameters$age_shares
    shares <- shares / sum(shares)
    
    Ns_g <- split_counts(Ns, shares)
    Ni_g <- split_counts(Ni, if (is.null(parameters$infect_init_shares)) shares else parameters$infect_init_shares)
    Nr_g <- split_counts(Nr, shares)
    Nd_g <- split_counts(Nd, shares)
    
    # Fallbacks if user didn't set vectors/matrix
    pi_vec <- if (!is.null(parameters$pi_by_age)) parameters$pi_by_age else rep(parameters$pi, G)
    v_vec  <- if (!is.null(parameters$v_by_age))  parameters$v_by_age  else rep(parameters$v,  G)
    
    C <- parameters$contact_matrix
    if (is.null(dim(C)) || any(dim(C) != c(G,G))) {
      C <- diag(G) # identity as safe default
    } else {
      # gently normalize rows to sum to 1
      C <- sweep(C, 1, rowSums(C), FUN = "/")
      C[!is.finite(C)] <- 0
    }
  }
  # --- END NEW ---
  # --- Age-dependent beta/gamma setup ---
  G <- length(parameters$age_groups)
  
  # γ by age (vector length G); fallback to scalar gamma
  gamma_vec <- if (!is.null(parameters$gamma_by_age)) {
    parameters$gamma_by_age
  } else {
    rep(parameters$gamma, G)
  }
  
  # Build β_{g,h} matrix
  if (!is.null(parameters$beta_matrix)) {
    beta_mat <- parameters$beta_matrix
    # sanity: coerce to GxG and nonnegative
    beta_mat <- matrix(beta_mat, nrow = G, ncol = G, byrow = FALSE)
    beta_mat[beta_mat < 0] <- 0
  } else {
    susc <- if (!is.null(parameters$beta_susc_by_age)) parameters$beta_susc_by_age else rep(1, G)
    inf  <- if (!is.null(parameters$beta_inf_by_age))  parameters$beta_inf_by_age  else rep(1, G)
    base <- parameters$beta
    beta_mat <- outer(susc, inf, FUN = "*") * base
  }
  
  # Combine β with mixing C: element-wise (keeps your interpretation of C as contact weights)
  beta_contact_mat <- beta_mat * C
  # --- end age-dependent beta/gamma setup ---
  
  for(i_day in times[-1]){
    
    beta_t <- if (bool_stochastic_beta) {
      max(0.0001, rnorm(1, mean = parameters$beta, sd = parameters$sigma))
    } else {
      parameters$beta
    }
    
    # === BRANCH: Heterogeneous vs Homogeneous ===
    if (!isTRUE(parameters$bool_age_heterogeneity)) {
      # ---- ORIGINAL PATH (unchanged logic) ----
      Ns_prop <- Ns / parameters$pop_size
      
      a_t <- if (isTRUE(parameters$bool_daily_cost_minimizing)) {
        a_function_cost_minimizing(Ns, Ni, i_day, parameters)
      } else {
        a_function(Ni, Ns, parameters)
      }
      
      u_t <- utility_function(a_t)
      
      if (!is.null(parameters$healthcare_capacity) && Ni > parameters$healthcare_capacity) {
        pi_effective <- parameters$pi * parameters$excess_mortality_multiplier
      } else {
        pi_effective <- parameters$pi
      }
      
      p_infect  <- 1 - exp(- beta_t * a_t^2 * Ni / parameters$pop_size)
      p_recover <- 1 - exp(- parameters$gamma)
      p_death   <- 1 - exp(- pi_effective)
      
      new_infections <- update_function(Ns, prob = p_infect)
      new_recoveries <- update_function(Ni, prob = p_recover)
      new_death      <- update_function(new_recoveries, prob = p_death)
      
      dNs <- -new_infections
      dNi <-  new_infections - new_recoveries
      dNr <-  new_recoveries - new_death
      dNd <-  new_death
      
      if((Ni - new_recoveries) < parameters$infect_thres){
        new_recoveries <- Ni
        new_infections <- 0
        dNs <- 0; dNi <- -Ni; dNr <- Ni; dNd <- 0
      }
      
      if (!isTRUE(parameters$bool_regular_sird)) {
        HealthCost        <- HealthCost + fx_per_capita * exp(-parameters$rho * i_day) * parameters$v * new_death
        SocialActivityCost<- SocialActivityCost + fx_per_capita * exp(-parameters$rho * i_day) * (Ns + Ni) * abs(u_t)
      }
      
      Rt <- calculate_Rt(parameters$R0, a_t, Ns/parameters$pop_size, Ni)
      
      Ns <- Ns + dNs; Ni <- Ni + dNi; Nr <- Nr + dNr; Nd <- Nd + dNd
      
    } else {
      # ---- NEW: AGE-HETEROGENEOUS PATH ----
      # ---- NEW: AGE-HETEROGENEOUS PATH ----
      
      G <- length(parameters$age_groups)
      
      # Activity per age
      if (isTRUE(parameters$bool_age_specific_activity)) {
        a_vec <- a_function_by_age(Ni_g, Ns_g, parameters)
      } else {
        a_common <- if (isTRUE(parameters$bool_daily_cost_minimizing)) {
          a_function_cost_minimizing(sum(Ns_g), sum(Ni_g), i_day, parameters)
        } else {
          a_function(sum(Ni_g), sum(Ns_g), parameters)
        }
        a_vec <- rep(a_common, G)
      }
      
      # Utility per age
      u_vec <- vapply(a_vec, utility_function, numeric(1))
      
      # Capacity bump for IFR (already vector by age)
      if (!is.null(parameters$healthcare_capacity) && sum(Ni_g) > parameters$healthcare_capacity) {
        pi_eff <- pi_vec * parameters$excess_mortality_multiplier
      } else {
        pi_eff <- pi_vec
      }
      
      # Stochastic β scaling: scale the whole β matrix (keeps relative age structure)
      beta_scale_t <- if (bool_stochastic_beta) {
        # convert absolute sigma on scalar beta into multiplicative noise on the matrix
        sd_mult <- parameters$sigma / max(1e-8, parameters$beta)
        max(1e-4, rnorm(1, mean = 1, sd = sd_mult))
      } else {
        1
      }
      B <- beta_scale_t * beta_contact_mat  # GxG
      
      # Force of infection for each susceptible age g:
      # λ_g = a_g * sum_h B[g,h] * a_h * I_h / Pop
      Ia <- a_vec * Ni_g                    # elementwise
      lambda_g <- a_vec * as.numeric(B %*% Ia) / parameters$pop_size
      
      # Transition probabilities (clamped into [0,1])
      p_infect_g  <- pmin(pmax(1 - exp(-lambda_g), 0), 1)
      p_recover_g <- pmin(pmax(1 - exp(-gamma_vec), 0), 1)
      p_death_g   <- pmin(pmax(1 - exp(-pi_eff),    0), 1)
      
      # Draw (or compute) transitions
      new_inf_g <- mapply(update_function, n = Ns_g, prob = p_infect_g)
      new_rec_g <- mapply(update_function, n = Ni_g, prob = p_recover_g)
      new_dth_g <- mapply(update_function, n = new_rec_g, prob = p_death_g)
      
      # Fadeout threshold on totals
      if ((sum(Ni_g) - sum(new_rec_g)) < parameters$infect_thres) {
        new_rec_g <- Ni_g
        new_inf_g <- rep(0, G)
        new_dth_g <- rep(0, G)
      }
      
      # Update by-age states
      Ns_g <- Ns_g - new_inf_g
      Ni_g <- Ni_g + new_inf_g - new_rec_g
      Nr_g <- Nr_g + new_rec_g - new_dth_g
      Nd_g <- Nd_g + new_dth_g
      
      # Aggregate
      Ns <- sum(Ns_g); Ni <- sum(Ni_g); Nr <- sum(Nr_g); Nd <- sum(Nd_g)
      
      # Costs (per capita)
      if (!isTRUE(parameters$bool_regular_sird)) {
        HealthCost <- HealthCost + fx_per_capita * exp(-parameters$rho * i_day) * sum(v_vec * new_dth_g)
        SocialActivityCost <- SocialActivityCost + fx_per_capita * exp(-parameters$rho * i_day) *
          sum( (Ns_g + Ni_g) * abs(u_vec) )
      }
      
      # Population-weighted outputs
      w <- (Ns_g + Ni_g); if (sum(w) == 0) w[] <- 1
      a_t <- sum(a_vec * w) / sum(w)
      u_t <- sum(u_vec * w) / sum(w)
      
      # Approximate Rt using per-age λ_g and γ_g
      Rt <- sum( (Ns_g / parameters$pop_size) * (lambda_g / pmax(1e-12, gamma_vec)) )
      
    }
    
    # keep track of the states (aggregate outputs, as before)
    states_out[i_day+1,] = c(Ns, Ni, Nr, Nd,
                             HealthCost, SocialActivityCost,
                             HealthCost + SocialActivityCost,
                             a_t, u_t, Rt)
  }
  
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