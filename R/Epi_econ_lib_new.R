
##
# Help function for the EPI-ECON modelling work
#
# - Deterministic and stochastic transitions
################################################################ #
#setwd("C:/Users/Bart Smets/OneDrive/Documenten/GitHub/Working_Code_BS/R")
library(confintr)
library(scales)
#Activity function calculation
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
  # setup & sanity 
  stopifnot(is.list(parameters), length(times) >= 1)
  
  # copy initial states
  states <- data.frame(t(initial_state))
  
  # convert population fractions into counts
  isN <- grepl("^N", names(states))
  states[isN] <- round(states[isN] * parameters$pop_size)
  
  # keep total population fixed: Ns = pop - (Ni+Nr+Nd)
  Ns_idx <- which(names(states) == "Ns")
  other_idx <- which(isN & names(states) != "Ns")
  states[[Ns_idx]] <- parameters$pop_size - sum(states[other_idx])
  
  # output matrix
  states_out <- matrix(NA_real_, nrow = length(times), ncol = length(states))
  colnames(states_out) <- names(states)
  
  # init stocks and costs
  Ns <- states$Ns; Ni <- states$Ni; Nr <- states$Nr; Nd <- states$Nd
  HealthCost <- 0; SocialActivityCost <- 0
  
  fx_per_capita <- parameters$fx / parameters$pop_size
  regular_mode  <- !is.null(parameters$bool_regular_sird) && isTRUE(parameters$bool_regular_sird)
  
  # ---- carry/accumulator switch ----
  use_integer_with_carry <- isTRUE(parameters$integer_with_carry)
  only_with_det_updater <- TRUE
  if (use_integer_with_carry && only_with_det_updater) {
    use_integer_with_carry <- identical(update_function, get_transitions_deterministic)
  }
  
  # local carry buffers (used only if use_integer_with_carry == TRUE)
  carry <- list(inf = 0, rec = 0, dea = 0)
  
  take_with_carry <- function(expected, key, cap = Inf) {
    if (!is.finite(expected) || expected < 0) expected <- 0
    carry[[key]] <<- carry[[key]] + expected
    realized <- floor(carry[[key]])
    if (realized > cap) realized <- cap
    if (realized < 0)  realized <- 0
    carry[[key]] <<- carry[[key]] - realized
    return(as.numeric(realized))
  }
  
  # helper to compute flows either via carry or updater
  realize_flows <- function(Ns_cur, Ni_cur, p_inf, p_rec, p_dea) {
    if (use_integer_with_carry) {
      exp_inf <- Ns_cur * p_inf
      exp_rec <- Ni_cur * p_rec
      exp_dea <- Ni_cur * p_dea
      new_inf <- take_with_carry(exp_inf, "inf", cap = Ns_cur)
      new_dea <- take_with_carry(exp_dea, "dea", cap = Ni_cur)
      new_rec <- take_with_carry(exp_rec, "rec", cap = max(Ni_cur - new_dea, 0))
      list(new_inf = new_inf, new_rec = new_rec, new_dea = new_dea)
    } else {
      new_inf <- update_function(Ns_cur, prob = p_inf)
      new_rec <- update_function(Ni_cur, prob = p_rec)
      new_dea <- update_function(Ni_cur, prob = p_dea)
      list(new_inf = new_inf, new_rec = new_rec, new_dea = new_dea)
    }
  }
  
  # record t = 0 row
  if (regular_mode) {
    a0 <- 1; u0 <- 0
  } else {
    a0 <- a_function(Ni, Ns, parameters)
    u0 <- utility_function(a0, parameters)
  }
  Rt0 <- calculate_Rt(parameters$R0, a0, Ns / parameters$pop_size, Ni)
  states_out[1, ] <- c(Ns, Ni, Nr, Nd, HealthCost, SocialActivityCost,
                       HealthCost + SocialActivityCost, a0, u0, Rt0)
  
  # ---- main loop ----
  for (i_day in times[-1]) {
    beta_t  <- parameters$beta
    Ns_prop <- Ns / parameters$pop_size
    
    if (regular_mode) {
      a_t <- 1
      u_t <- 0
      
      p_infect  <- 1 - exp(- beta_t * (Ni / parameters$pop_size))
      p_recover <- 1 - exp(- (1 - parameters$pi) * parameters$gamma)
      p_death   <- 1 - exp(- parameters$pi * parameters$gamma)
      
      flows <- realize_flows(Ns, Ni, p_infect, p_recover, p_death)
      new_infections <- flows$new_inf
      new_recoveries <- flows$new_rec
      new_death      <- flows$new_dea
      
      if ((Ni - new_recoveries) < parameters$infect_thres) {
        new_recoveries <- Ni
        new_infections <- 0
        carry$rec <- 0; carry$inf <- 0; carry$dea <- 0
      }
      
      dNs <- -new_infections
      dNi <-  new_infections - new_recoveries - new_death
      dNr <-  new_recoveries
      dNd <-  new_death
      
      # always realized health cost only
      HealthCost <- HealthCost + fx_per_capita * exp(-parameters$rho * i_day) * 
        parameters$v * new_death
      
      Rt <- calculate_Rt(parameters$R0, a_t, Ns_prop, Ni)
      
    } else {
      a_t <- a_function(Ni, Ns, parameters)
      u_t <- utility_function(a_t, parameters)
      
      p_infect  <- 1 - exp(- beta_t * a_t^2 * (Ni / parameters$pop_size))
      p_recover <- 1 - exp(- (1 - parameters$pi) * parameters$gamma)
      p_death   <- 1 - exp(- parameters$pi * parameters$gamma)
      
      flows <- realize_flows(Ns, Ni, p_infect, p_recover, p_death)
      new_infections <- flows$new_inf
      new_recoveries <- flows$new_rec
      new_death      <- flows$new_dea
      
      if ((Ni - new_recoveries) < parameters$infect_thres) {
        new_recoveries <- Ni
        new_infections <- 0
        carry$rec <- 0; carry$inf <- 0; carry$dea <- 0
      }
      
      dNs <- -new_infections
      dNi <-  new_infections - new_recoveries - new_death
      dNr <-  new_recoveries
      dNd <-  new_death
      
      # always realized health cost only
      HealthCost <- HealthCost + fx_per_capita * exp(-parameters$rho * i_day) * 
        parameters$v * new_death
      
      SocialActivityCost <- SocialActivityCost +
        fx_per_capita * exp(-parameters$rho * i_day) * (Ns + Ni) * abs(u_t)
      
      Rt <- calculate_Rt(parameters$R0, a_t, Ns_prop, Ni)
    }
    
    Ns <- Ns + dNs; Ni <- Ni + dNi; Nr <- Nr + dNr; Nd <- Nd + dNd
    
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

compare_sim_output <- function(output_experiments, output_deterministic,
                               plot_tag, fadeout_threshold = 100,
                               compare_carry = FALSE,
                               output_det_nocarry = NULL,
                               output_det_carry   = NULL) {
  if (!dir.exists("figures")) dir.create("figures")
  pdf(file = paste0("figures/", gsub("[^[:alnum:]]", "_", plot_tag), ".pdf"))
  op <- par(cex.lab = 1.4, cex.axis = 1.2)
  on.exit({ par(op); dev.off() }, add = TRUE)
  
  # -------------------------
  # Labels for plotted series
  # -------------------------
  labels <- c(
    "a_t" = "Activity",
    "Ni"  = "Number of infectious individuals",
    "Ns"  = "Susceptibles",
    "Nr"  = "Recovered",
    "Nd"  = "Number of deaths",
    "u_t" = "Utility",
    "Rt"  = "Reproduction Number",
    "HealthCost" = "Health cost",
    "SocialActivityCost" = "Activity cost",
    "TotalCost" = "Health cost"
  )
  
  # =======================================================================================
  # MODE A: Deterministic relative vs deterministic absolute (no stochastic distributions)
  # =======================================================================================
  if (isTRUE(compare_carry)) {
    stopifnot(!is.null(output_det_nocarry), !is.null(output_det_carry))
    
    rel_name <- "Relative frequencies"
    abs_name <- "Absolute frequencies"
    
    # 1) Final costs barplot
    last_rel <- output_det_nocarry[nrow(output_det_nocarry), c("HealthCost","SocialActivityCost","TotalCost")]
    last_abs <- output_det_carry  [nrow(output_det_carry  ), c("HealthCost","SocialActivityCost","TotalCost")]
    
    cost_mat <- rbind(
      `Deterministic relative frequencies (continuous compartments)` = as.numeric(last_rel),
      `Deterministic absolute frequencies (discrete integer idea)`   = as.numeric(last_abs)
    )
    colnames(cost_mat) <- c("HealthCost","ActivityCost","TotalCost")
    
    y_lim <- range(cost_mat, na.rm = TRUE)
    barplot(cost_mat, beside = TRUE, ylim = y_lim, las = 2,
            ylab = "Per-capita cost",
            main = paste(plot_tag, "(final costs)"),
            col = gray(c(0.8, 0.5)))
   # legend("topleft", fill = gray(c(0.8, 0.5)),
       #    legend = c(rel_name, abs_name), bty = "n")
    
    # Console summary table
    print_out <- data.frame(
      output = c('Health Cost (per capita)','Activity Cost (per capita)','Total cost (per capita)'),
      rel    = c(last_rel$HealthCost, last_rel$SocialActivityCost, last_rel$TotalCost),
      abs    = c(last_abs$HealthCost, last_abs$SocialActivityCost, last_abs$TotalCost)
    )
    names(print_out)[2:3] <- c(rel_name, abs_name)
    names(print_out)[1]   <- plot_tag
    print(paste('SUMMARY STATISTICS:', plot_tag))
    print(print_out)
    
    # 2) Time series: relative vs absolute
    all_cols   <- intersect(colnames(output_det_nocarry), colnames(output_det_carry))
    sel_states <- setdiff(all_cols, c("TotalCost","Ns"))
    time_vec <- seq_len(nrow(output_det_nocarry)) - 1
    
    for (i_state in sel_states) {
      plot_label <- ifelse(i_state %in% names(labels), labels[[i_state]], i_state)
      y_min <- min(output_det_nocarry[, i_state], output_det_carry[, i_state], na.rm = TRUE)
      y_max <- max(output_det_nocarry[, i_state], output_det_carry[, i_state], na.rm = TRUE)
      
      # plot(NULL, xlim = c(min(time_vec), max(time_vec)),
      #      ylim = c(y_min, y_max), xlab = "Days", ylab = plot_label)
      
      if (i_state == "Ni") {
        # Force both axes to start at 0
        plot(NULL,
             xlim = c(0, max(time_vec, na.rm = TRUE)),
             ylim = c(0, max(y_max, 0, na.rm = TRUE)),
             xlab = "Days",
             ylab = "Number of infectious individuals",
             xaxs = "i", yaxs = "i")
        abline(h = 0, lty = 3, col = "grey60")   # ← drop v = 0
        
      } else if (i_state == "a_t") {
        # Force only x-axis to start at 0
        plot(NULL,
             xlim = c(0, max(time_vec, na.rm = TRUE)),
             ylim = c(y_min, y_max),
             xlab = "Days",
             ylab = plot_label,
             xaxs = "i")
        abline(h = 0, lty = 3, col = "grey60")
        
        
      } else {
        # Default for all other variables
        plot(NULL,
             xlim = c(min(time_vec, na.rm = TRUE), max(time_vec, na.rm = TRUE)),
             ylim = c(y_min, y_max),
             xlab = "Days",
             ylab = plot_label)
      }
      
      # Add faint reference lines for clarity
   #   abline(h = 0, v = 0, lty = 3, col = "grey60")
  
      # make a semi-transparent yellow once
      trans_yellow <- grDevices::adjustcolor("yellow", alpha.f = 0.35)
      # Relative (continuous) = dashed black; Absolute (discrete) = solid blue
      lines(time_vec, output_det_nocarry[, i_state], col = "black", lwd = 2, lty = 2)
      lines(time_vec, output_det_carry  [, i_state], col = trans_yellow,  lwd = 2)
      
    #  legend("topleft",
       #      legend = c(rel_name, abs_name),
        #     col = c("black",trans_yellow), lwd = 2, lty = c(2,1), bty = "n")
    }
    
    return(invisible(NULL))
  }
  
  # =========================================================================
  # MODE B (default): stochastic vs deterministic  — original behavior
  # =========================================================================
  output_summary <- output_experiments$output_summary
  output_all     <- output_experiments$output_all
  
  output_summary_deterministic <- output_deterministic[nrow(output_deterministic), ]
  
  if (fadeout_threshold > 0) {
    exp_fade_out   <- (output_summary$Nr + output_summary$Ni) < fadeout_threshold
    output_summary <- output_summary[!exp_fade_out, , drop = FALSE]
    output_all     <- output_all[!exp_fade_out, , , drop = FALSE]
    plot_tag       <- paste(plot_tag, '(excl. fadeout)')
  } else {
    plot_tag <- paste(plot_tag, '(all)')
  }
  
  if (nrow(output_summary) == 0) {
    plot.new(); title("No simulations left after fadeout filter")
    return(invisible(NULL))
  }
  
  # Boxplot for costs (stochastic) + deterministic markers
  output_cost <- output_summary[, grepl('Cost', names(output_summary)), drop = FALSE]
  y_lim <- range(output_cost, na.rm = TRUE)
  boxplot(output_cost, las = 2, ylim = y_lim)
  points(colMeans(output_cost), pch = 8)
  points((1:3) + 0.2, output_summary_deterministic[names(output_cost)], col = 4, pch = 8, lwd = 3)
  
  # Time series with quantile & CI bands
  all_cols   <- colnames(output_deterministic)
  sel_states <- setdiff(all_cols, c("TotalCost","Ns"))
  
  for (i_state in sel_states) {
    plot_label <- ifelse(i_state %in% names(labels), labels[[i_state]], i_state)
    var_matrix <- output_all[, , i_state]
    time_vec   <- 0:(ncol(var_matrix) - 1)
    
    ci_lower <- apply(var_matrix, 2, function(x) quantile(x, 0.025, na.rm = TRUE))
    ci_upper <- apply(var_matrix, 2, function(x) quantile(x, 0.975, na.rm = TRUE))
    
    n_sim   <- nrow(var_matrix)
    ci_mean <- apply(var_matrix, 2, mean, na.rm = TRUE)
    ci_sd   <- apply(var_matrix, 2, sd, na.rm = TRUE)
    ci_se   <- ci_sd / sqrt(n_sim)
    ci_lo_meanband <- ci_mean - 1.96 * ci_se
    ci_hi_meanband <- ci_mean + 1.96 * ci_se
    
    y_min <- min(ci_lower, ci_lo_meanband, output_deterministic[, i_state], na.rm = TRUE)
    y_max <- max(ci_upper, ci_hi_meanband, output_deterministic[, i_state], na.rm = TRUE)
    
    # plot(NULL, xlim = c(min(time_vec), max(time_vec)),
    #      ylim = c(y_min, y_max), xlab = "Days", ylab = plot_label)
    # 
    if (i_state == "Ni") {
      # Force both axes to start at 0
      plot(NULL,
           xlim = c(0, max(time_vec, na.rm = TRUE)),
           ylim = c(0, max(y_max, 0, na.rm = TRUE)),
           xlab = "Days",
           ylab = "Number of infectious individuals",
           xaxs = "i", yaxs = "i")
    #  abline(h = 0, lty = 3, col = "grey60")   # ← drop v = 0
      
    } else if (i_state == "a_t") {
      # Force only x-axis to start at 0
      plot(NULL,
           xlim = c(0, max(time_vec, na.rm = TRUE)),
           ylim = c(y_min, y_max),
           xlab = "Days",
           ylab = plot_label,
           xaxs = "i")
      abline(h = 0, lty = 3, col = "grey60")
      
    } else {
      # Default for all other variables
      plot(NULL,
           xlim = c(min(time_vec, na.rm = TRUE), max(time_vec, na.rm = TRUE)),
           ylim = c(y_min, y_max),
           xlab = "Days",
           ylab = plot_label)
    }
    
    # Add faint reference lines for clarity
    abline(h = 0, v = 0, lty = 3, col = "grey60")
    
    
    #abline(h = 0, lty = 3, col = "grey60")
    
    polygon(c(time_vec, rev(time_vec)), c(ci_upper, rev(ci_lower)),
            col = scales::alpha("lightblue", 0.4), border = NA)
    polygon(c(time_vec, rev(time_vec)), c(ci_hi_meanband, rev(ci_lo_meanband)),
            col = scales::alpha("red", 0.4), border = NA)
    
    lines(time_vec, ci_mean, col = "blue", lwd = 2)
    lines(output_deterministic[, i_state], col = "black", lwd = 2, lty = 2)
  }
  
  # === Original console summaries ===
  
  # Mean & 95% CI for costs
  mean_health <- mean(output_summary$HealthCost, na.rm = TRUE)
  ci_health <- ci_mean(output_summary$HealthCost)$interval
  
  mean_social <- mean(output_summary$SocialActivityCost, na.rm = TRUE)
  ci_social <- ci_mean(output_summary$SocialActivityCost)$interval
  
  cat("Stochastic Summary (from boxplot data):\n")
  cat(sprintf("Health Cost: %.0f [%.0f, %.0f]\n", mean_health, ci_health[1], ci_health[2]))
  cat(sprintf("Social Activity Cost: %.0f [%.0f, %.0f]\n", mean_social, ci_social[1], ci_social[2]))
  
  # Quantile interval % widths
  health_quant <- quantile(output_summary$HealthCost, probs = c(0.025, 0.975), na.rm = TRUE)
  pct_health_upper <- 100 * (health_quant[2] - mean_health) / mean_health
  pct_health_lower <- 100 * (mean_health - health_quant[1]) / mean_health
  
  mean_activity_cost <- mean(output_summary$SocialActivityCost, na.rm = TRUE)
  activity_quant <- quantile(output_summary$SocialActivityCost, probs = c(0.025, 0.975), na.rm = TRUE)
  pct_activity_upper <- 100 * (activity_quant[2] - mean_activity_cost) / mean_activity_cost
  pct_activity_lower <- 100 * (mean_activity_cost - activity_quant[1]) / mean_activity_cost
  
  # Mean vs deterministic (%)
  health_det <- output_deterministic[nrow(output_deterministic), "HealthCost"]
  activity_det <- output_deterministic[nrow(output_deterministic), "SocialActivityCost"]
  pct_health_mean_diff <- 100 * (health_det - mean_health) / health_det
  pct_activity_mean_diff <- 100 * (activity_det - mean_activity_cost) / activity_det
  
  # Peak-time metrics
  t_peak <- which.max(output_deterministic$Ni)
  a_vals <- output_all[, t_peak, "a_t"]
  ni_vals <- output_all[, t_peak, "Ni"]
  
  a_mean <- mean(a_vals, na.rm = TRUE)
  ni_mean <- mean(ni_vals, na.rm = TRUE)
  
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
  cat(sprintf("Health Cost – Quantile Interval: +%.2f%% / -%.2f%%\n", pct_health_upper, pct_health_lower))
  cat(sprintf("Activity Cost – Quantile Interval: +%.2f%% / -%.2f%%\n", pct_activity_upper, pct_activity_lower))
  cat(sprintf("Health Cost – Stochastic vs Det: %.2f%%\n", pct_health_mean_diff))
  cat(sprintf("Activity Cost – Stochastic vs Det: %.2f%%\n", pct_activity_mean_diff))
  
  cat("\nPEAK-TIME DIFFERENCES (Ni and a_t at peak Ni)\n")
  cat(sprintf("Activity a(t) – Quantile Interval: +%.2f%% / -%.2f%%\n", pct_a_upper, pct_a_lower))
  cat(sprintf("Infectives Ni – Quantile Interval: +%.2f%% / -%.2f%%\n", pct_ni_upper, pct_ni_lower))
  cat(sprintf("Activity a(t) – Stochastic vs Det: %.2f%%\n", pct_a_mean_diff))
  cat(sprintf("Infectives Ni – Stochastic vs Det: %.2f%%\n", pct_ni_mean_diff))
  
  # Health Cost box + CIs & quantiles
  boxplot(output_summary$HealthCost,
          ylab = "Health cost",
          ylim = range(c(output_summary$HealthCost, output_summary_deterministic$HealthCost), na.rm = TRUE),
          col = "gray90")
  ci_health_full <- ci_mean(output_summary$HealthCost)
  segments(x0 = 1, x1 = 1, y0 = ci_health_full$interval[1], y1 = ci_health_full$interval[2],
           col = "forestgreen", lwd = 2)
  quant_health <- quantile(output_summary$HealthCost, probs = c(0.025, 0.975), na.rm = TRUE)
  segments(x0 = 1.1, x1 = 1.1, y0 = quant_health[1], y1 = quant_health[2],
           col = "blue", lwd = 2, lty = 2)
  abline(h = output_summary_deterministic$HealthCost, col = "red", lwd = 2)
  abline(h = mean(output_summary$HealthCost, na.rm = TRUE), col = "purple", lwd = 2)
  
  # Activity Cost box + CIs & quantiles
  boxplot(output_summary$SocialActivityCost,
          ylab = "Activity cost",
          ylim = range(c(output_summary$SocialActivityCost, output_summary_deterministic$SocialActivityCost), na.rm = TRUE),
          col = "gray90")
  abline(h = mean(output_summary$SocialActivityCost, na.rm = TRUE), col = "purple", lwd = 2)
  abline(h = output_summary_deterministic$SocialActivityCost, col = "red", lwd = 2)
  ci_social_full <- ci_mean(output_summary$SocialActivityCost)
  segments(x0 = 1, x1 = 1, y0 = ci_social_full$interval[1], y1 = ci_social_full$interval[2],
           col = "forestgreen", lwd = 2)
  quant_social <- quantile(output_summary$SocialActivityCost, probs = c(0.025, 0.975), na.rm = TRUE)
  segments(x0 = 1.1, x1 = 1.1, y0 = quant_social[1], y1 = quant_social[2],
           col = "blue", lwd = 2, lty = 2)
  
  # Summary table
  print_out <- data.frame(
    output = c('Health Cost (per capita)','Activity Cost (per capita)','Total Cost (per capita)'),
    deterministic = get_summary_stats(output_deterministic[nrow(output_deterministic),]),
    stochastic    = get_summary_stats(output_summary)
  )
  names(print_out)[1] <- plot_tag
  print(paste('SUMMARY STATISTICS:', plot_tag))
  print(print_out)
  
  # Peak Infection Timing Summary
  peak_times <- output_summary$PeakTime
  peak_ci <- ci_mean(peak_times)
  cat("\nPeak Infection Timing Summary:\n")
  cat(sprintf("Mean Peak Time: %.1f [%.1f, %.1f]  (mean ± 1.96 * SE)\n",
              peak_ci$estimate, peak_ci$interval[1], peak_ci$interval[2]))
  cat(sprintf("Quantiles: 2.5%% = %.1f, 50%% (median) = %.1f, 97.5%% = %.1f\n",
              quantile(peak_times, 0.025, na.rm = TRUE),
              quantile(peak_times, 0.50, na.rm = TRUE),
              quantile(peak_times, 0.975, na.rm = TRUE)))
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

################################################################
# Sensitivity analysis helpers  (no kappa)
# - Parameter sweeps (compute_sensitivity_grid)
# - Number-of-simulations analysis (compute_delta_vs_num_sims)
################################################################

# Core engine to compute one Δ-cost point for a given parameter value
# Notes:
# - For param_name == "Ni0", 'value' is an ABSOLUTE COUNT; we convert to proportion.
# - For pi and v (and any other numeric param), we set them directly.
# - We recompute R0 = beta/gamma after any change to beta or gamma.
compute_delta_cost_point <- function(baseline_params, param_name, value,
                                     num_experiments_sens = 300,
                                     fadeout_threshold_for_diff = 100,
                                     update_function_stoch = get_transitions_stochastic,
                                     update_function_det   = get_transitions_deterministic) {
  # clone + set
  params <- baseline_params
  
  if (param_name == "Ni0") {
    # 'value' is count of initially infected persons
    params$ni0 <- value / params$pop_size
  } else {
    # any direct numeric parameter (e.g., pi, v, beta, gamma, rho, pop_size, time_horizon)
    params[[param_name]] <- value
  }
  
  # keep R0 consistent with beta/gamma if either changed
  params$R0 <- params$beta / params$gamma
  
  # init state & horizon
  init_state <- c(
    Ns = params$ns0, Ni = params$ni0, Nr = params$nr0, Nd = params$nd0,
    HealthCost = 0, SocialActivityCost = 0, TotalCost = 0,
    a_t = NA, u_t = NA, Rt = NA
  )
  times <- 0:params$time_horizon
  
  # deterministic baseline
  det_out   <- run_sir_binomial(init_state, times, params, update_function_det)
  det_total <- det_out[nrow(det_out), "TotalCost"]
  
  # stochastic experiments
  sims <- run_experiments(init_state, times, params, update_function_stoch,
                          num_experiments = num_experiments_sens)
  summ <- sims$output_summary
  
  # optionally drop fadeouts (very small outbreaks)
  if (fadeout_threshold_for_diff > 0) {
    keep <- (summ$Nr + summ$Ni) >= fadeout_threshold_for_diff
    if (!any(keep)) return(NULL)
    summ <- summ[keep, , drop = FALSE]
  }
  
  stoch_totals <- summ$TotalCost
  ci <- confintr::ci_mean(stoch_totals)
  q  <- stats::quantile(stoch_totals, c(0.025, 0.975), na.rm = TRUE)
  
  data.frame(
    Parameter   = param_name,
    Value       = value,
    Det_Total   = det_total,
    Stoch_Mean  = as.numeric(ci$estimate),
    CI_Lo       = as.numeric(ci$interval[1]),
    CI_Hi       = as.numeric(ci$interval[2]),
    Q_Lo        = as.numeric(q[1]),
    Q_Hi        = as.numeric(q[2]),
    Diff_Mean   = as.numeric(ci$estimate) - det_total,
    Diff_CI_Lo  = as.numeric(ci$interval[1]) - det_total,
    Diff_CI_Hi  = as.numeric(ci$interval[2]) - det_total,
    Diff_Q_Lo   = as.numeric(q[1]) - det_total,
    Diff_Q_Hi   = as.numeric(q[2]) - det_total
  )
}

# Vectorized sweep over a grid for one parameter
compute_sensitivity_grid <- function(baseline_params, param_name, values,
                                     num_experiments_sens = 300,
                                     fadeout_threshold_for_diff = 100) {
  out <- lapply(values, function(v)
    compute_delta_cost_point(
      baseline_params, param_name, v,
      num_experiments_sens, fadeout_threshold_for_diff
    )
  )
  do.call(rbind, Filter(Negate(is.null), out))
}

# Δ-cost vs number of simulations (matches manual spec setup)
compute_delta_vs_num_sims <- function(baseline_params, n_grid,
                                      fadeout_threshold_for_diff = 100) {
  # local copy and recompute R0 exactly as elsewhere
  params <- baseline_params
  params$R0 <- params$beta / params$gamma
  
  init_state <- c(
    Ns = params$ns0, Ni = params$ni0, Nr = params$nr0, Nd = params$nd0,
    HealthCost = 0, SocialActivityCost = 0, TotalCost = 0,
    a_t = NA, u_t = NA, Rt = NA
  )
  times <- 0:params$time_horizon
  
  det_out   <- run_sir_binomial(init_state, times, params, get_transitions_deterministic)
  det_total <- det_out[nrow(det_out), "TotalCost"]
  
  rows <- lapply(n_grid, function(N) {
    sims <- run_experiments(init_state, times, params, get_transitions_stochastic,
                            num_experiments = N)
    summ <- sims$output_summary
    if (fadeout_threshold_for_diff > 0) {
      keep <- (summ$Nr + summ$Ni) >= fadeout_threshold_for_diff
      if (!any(keep)) return(NULL)
      summ <- summ[keep, , drop = FALSE]
    }
    ci <- confintr::ci_mean(summ$TotalCost)
    data.frame(
      num_sims    = N,
      det_total   = det_total,
      stoch_mean  = as.numeric(ci$estimate),
      stoch_ci_lo = as.numeric(ci$interval[1]),
      stoch_ci_hi = as.numeric(ci$interval[2]),
      diff_mean   = as.numeric(ci$estimate) - det_total,
      diff_ci_lo  = as.numeric(ci$interval[1]) - det_total,
      diff_ci_hi  = as.numeric(ci$interval[2]) - det_total
    )
  })
  
  do.call(rbind, Filter(Negate(is.null), rows))
}
