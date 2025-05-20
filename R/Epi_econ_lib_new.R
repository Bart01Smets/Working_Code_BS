
#
# Help function for the EPI-ECON modelling work
#
# - Deterministic and stochastic transitions
################################################################ #
#setwd("C:/Users/Bart Smets/OneDrive/Documenten/GitHub/Working_Code_BS/R")
library(confintr)
library(scales)

# Functions to calculate the optimal action (a_t)

# #1: Original function
# 
# a_function <- function(alpha, beta, Ns, Ni, lambda_s, lambda_i, utility_type, scenario, tolerance) {
#   denom <- max(abs(lambda_s - lambda_i), tolerance)  # Avoid division by zero#Lambda_s -Lambda_i
# 
# #test
# 
#   if (utility_type == "Log") {
#     sqrt_arg <- (Ns + Ni)^2 + 8 * beta * Ns * Ni * denom * (Ns + Ni)
#    return((- (Ns + Ni) + sqrt(sqrt_arg)) / (4 * beta * Ns * Ni * denom))
# 
#   if  sqrt_arg <- (Ns + Ni) + 4 * (1 + alpha) * beta * Ni * Ns * denom
#   return( (-(Ns + Ni) + sqrt(sqrt_arg)) / (2 * (1 + alpha) * beta * Ni * Ns * denom) )



# #2: Optimal activity keeping the repreoduction number below 1.
# 
# a_function <- function(R0, Ns_prop) {
#   if (is.na(R0) || is.na(Ns_prop) || R0 <= 0 || Ns_prop <= 0) {
#     warning("Invalid R0 or Ns_prop")
#     return(0)
#   }
# 
#   a_t <- sqrt(1 / (R0 * Ns_prop))
# 
#   return(max(0, min(1, a_t)))  # Ensure stays between 0 and 1
# }


# #3: Simple infections rule
# 
# a_function <- function(Ni, parameters) {
#   Ni_prop <- Ni / parameters$pop_size
#   sensitivity <- 50  # Tune this number to control steepness
#   a_t <- 1 / (1 + sensitivity * Ni_prop)#Ni_prop
#   if (is.na(Ni) || is.null(parameters$pop_size) || is.na(parameters$pop_size)) {
#     warning("Invalid inputs to a_function")
#     return(0)
#   }
# 
#  return(max(0, min(1, a_t)))
# }

# #4: Simple rule based on previous activity.
# 
# #Updated a_function with memory of previous a_t
# a_function <- function(Rt, a_t_previous, sensitivity_up = 0.05, sensitivity_down = 0.05) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
# 
#   if (Rt > 1) {
#     # If Rt > 1, decrease a_t (scaled by sensitivity factor)
#     delta <- -sensitivity_down #* (Rt - 1)
#   } else {
#     # If Rt < 1, increase a_t (scaled by sensitivity factor)
#     delta <- sensitivity_up #* (1 - Rt)
#   }
# 
#   # Update a_t
#   a_t_new <- a_t_previous*(1+ delta)
# 
#   # Clamp between 0 and 1
#   a_t_new <- max(0, min(1, a_t_new))
# 
#   return(a_t_new)
# }


#5: Myopic Laissez-faire
# Myopic altruistic activity function (laissez-faire equilibrium)
a_function <- function(Ns, Ni, parameters) {#
  Ni_prop <- Ni / parameters$pop_size
  Ns_prop <-  Ns/parameters$pop_size#parameters$ns0  # Assume constant or use from state if dynamic
  epsilon <- 1e-12
  Ni_prop <- max(Ni_prop, epsilon)

  # Compute multiplier
  multiplier <- parameters$beta * Ni_prop * parameters$kappa * (1 + parameters$alpha*Ns_prop)#Ns_prop

  sqrt_term <- sqrt(1 + 4 * multiplier)
  a_t <- (-1 + sqrt_term) / (2 * multiplier)

 #  Clamp a_t between [0,1]
 return (a_t)

 #(max(0, min(1, a_t)))
 }
# a_function <- function(Ni, Ns, parameters) {
#   Ni_prop <- Ni / parameters$pop_size
#   Ns_prop <- Ns / parameters$pop_size
#   multiplier <- parameters$beta * Ni_prop * Ns_prop * parameters$kappa * (1 + parameters$alpha * Ns_prop)
#   sqrt_term <- sqrt(1 + 8 * multiplier)
#   a_t <- (-1 + sqrt_term) / (4 * multiplier)
#   return(max(0, min(1, a_t)))
# }

# # McAdams–Day-style Laissez-Faire Activity Rule
# a_function <- function(Ni, Ns, parameters) {
#   Ni_prop <- Ni / parameters$pop_size
#   Ns_prop <- Ns / parameters$pop_size
# 
#   b1   <- parameters$b1
#   b2   <- parameters$b2
#   beta <- parameters$beta
#   H    <- parameters$H  # expected harm from infection
# 
#   denom <- 2 * beta * Ns_prop * Ni_prop * H - 2 * b2
#   if (denom <= 0) return(1)  # no infection risk → full activity
# 
#   a_star <- b1 / denom
#   return(max(0, min(1, a_star)))  # clamp to [0, 1]
# }
# #McAdams–Day-style Optimal Policy Activity Rule (Social Planner)
# a_function <- function(Ni, Ns, parameters) {
#   Ni_prop <- Ni / parameters$pop_size
#   Ns_prop <- Ns / parameters$pop_size
# 
#   b1   <- parameters$b1
#   b2   <- parameters$b2
#   beta <- parameters$beta
#   H    <- parameters$H
# 
#   total_prop <- Ns_prop + Ni_prop
#   denom <- 2 * beta * Ns_prop * Ni_prop * H - 2 * total_prop * b2
#   if (denom <= 0) return(1)
# 
#   numer <- total_prop * b1
#   a_star <- numer / denom
#   return(max(0, min(1, a_star)))
# }

# # Activity function: Laissez-Faire with Keppo-style myopic utility
# a_function <- function(Ni, Ns, parameters) {
#   pi <- Ni / parameters$pop_size
#   beta <- parameters$beta
#   gamma <- parameters$gamma
#   L <- parameters$L
#   sigma <- parameters$sigma
#   q_pi <- sigma / (sigma + pi)
# 
#   denom <- 2 * gamma * beta * pi * q_pi * L
#   if (denom <= 0) return(1)
# 
#   a_star <- denom^(-1 / (2 * gamma - 1))
#   return(max(0, min(1, a_star)))
# }
# 
# #Bethune, Korinek
# 
# a_function <- function(Ni, Ns, parameters) {
#   # Parameters needed
#   phi <- parameters$phi       # weight on utility
#   beta <- parameters$beta     # transmission rate
#   kappa <- parameters$kappa   # expected cost of infection
#   A <- 1                      # average activity level (externalized in laissez-faire)
# 
#   # Defensive programming
#   if (Ni <= 0 || is.na(Ni)) return(1)  # if no infections, full activity
# 
#   # Core expression
#   numerator <- -phi + sqrt(phi^2 + 8 * beta * A * Ni * kappa * phi)
#   denominator <- 4 * beta * A * Ni * kappa
# 
#   a_t <- numerator / denominator
# 
#   # Clamp a_t between 0 and 1
#   return(max(0, min(1, a_t)))
# }


# if (nrow(output_summary) == 0) {
#   warning("All simulations faded out — no data to plot.")
#   return(NULL)
# }







# # Activity function: Optimal Policy (Planner) with Keppo-style welfare
# a_function <- function(Ni, Ns, parameters) {
#   S <- Ns / parameters$pop_size
#   I <- Ni / parameters$pop_size
# 
#   beta <- parameters$beta
#   gamma <- parameters$gamma
#   L <- parameters$L
# 
#   numerator <- S + I
#   denominator <- 2 * gamma * beta * S * I * L
# 
#   if (denominator <= 0) return(1)
# 
#   a_star <- (numerator / denominator)^(1 / (2 * gamma - 1))
#   return(max(0, min(1, a_star)))
# }


# #6 a_t for optimal policy in the myopic scenario
# # Optimal policy activity function (derived from planner's FOC)
# a_function <- function(Ni, Ns, parameters) {
#   # Defensive checks
#     Ni_prop <- Ni / parameters$pop_size
#      Ns_prop <- Ns / parameters$pop_size
#      epsilon <- 1e-12  # for numerical safety
#      Ni_prop <- max(Ni_prop, epsilon)
#      Ns_prop <- max(Ns_prop, epsilon)
#   
#   
#   if (is.na(Ns) || is.na(Ni) || Ns <= 0 || Ni <= 0) {
#     warning("Invalid Ns or Ni in a_function()")
#     return(0)
#   }
# 
#   beta <- parameters$beta
#   kappa <- parameters$kappa
#   denom <- beta * Ns_prop * Ni * kappa / (Ns_prop + Ni_prop)
#   sqrt_term <- sqrt(1 + 8 * denom)
# 
#   a_t <- (-1 + sqrt_term) / (4 * denom)
# 
#   return(max(0, min(1, a_t)))  # Clamp to [0, 1]
# }

# #7: Optimal policy for behavioral rule
# 
# # # Optimal policy activity function: conservative and persistent distancing
# a_function <- function(Ni, parameters) {
#   Ni_prop <- Ni / parameters$pop_size
#   base_reduction <- 0.4    # baseline reduction even at low infection
#   max_suppression <- 0.9   # maximum reduction at high prevalence
#   sensitivity <- 300       # higher = more reactive to small infections
# 
#   a_t <- 1 - (base_reduction + (1 - base_reduction) * (1 - exp(-sensitivity * Ni_prop)))
# 
#   # Bound the activity level between 0 and 1
#  # a_t <- max(0, min(1, a_t))
# 
#   if (is.na(a_t)) {
#     warning("Invalid inputs to a_function")
#     return(0)
#   }
# 
#   return(a_t)
# }
# #8
# # Myopic Optimal Activity Function (Social Planner)
# a_function <- function(Ni, Ns, parameters) {
#   Ni_prop <- Ni / parameters$pop_size
#   Ns_prop <- Ns / parameters$pop_size
#   epsilon <- 1e-12  # for numerical safety
#   Ni_prop <- max(Ni_prop, epsilon)
#   Ns_prop <- max(Ns_prop, epsilon)
#   
#   # Compute M = β * Ni * Ns * κ * (1 + α)
#   multiplier <- parameters$beta * Ni_prop * Ns_prop * parameters$kappa * (1 + parameters$alpha)
#   
#   # Avoid division by zero
#   if (multiplier <= 0) return(1)
#   
#   # Compute a(t)
#   sqrt_term <- sqrt(1 + 8 * multiplier)
#   a_t <- (-1 + sqrt_term) / (4 * multiplier)
#   
#   # Clamp to [0, 1] for realism
#   return(max(0, min(1, a_t)))
# }
# #9: rule of thumb social activity based on Ni for optimal policy
# # Rule-of-thumb social activity function based on infection level Ni
# a_function <- function(Ni) {
#   # Hardcoded parameters
#   vigilance_threshold <- 0.001   # threshold infection rate (e.g. 0.1%)
#   vigilance_sensitivity <- 5000  # steepness of response
#   
#   # Edge case: if no infections, full activity
#   if (Ni <= 0) {
#     return(1)
#   }
#   
#   # Vigilance-based suppression (inverse logistic-like decay)
#   a_t <- 1 / (1 + vigilance_sensitivity * (Ni / vigilance_threshold)^2)
#   
#   # Ensure output is in [0,1]
#   return(max(min(a_t, 1), 0))
# }

# a_function <- function(Ni, Ns, Rt, a_t_prev, parameters) {
#   scenario <- parameters$scenario
# 
#   if (scenario == "log-optimal") {
#   #   Implement closed-form from FOC (lambda-based)
#   #   [Insert your original utility-based solution]
#   } else if (scenario == "Rt-threshold") {
#     # a_t decreases when Rt > 1
#     if (Rt > 1) {
#       return(max(0, a_t_prev * 0.9))
#     } else {
#       return(min(1, a_t_prev * 1.01))
#     }
#   } else if (scenario == "infection-vigilance") {
#     # Simple infection-based suppression (your current rule)
#     Ni_prop <- Ni / parameters$pop_size
#     threshold <- 0.001
#     sensitivity <- 5000
#     a_t <- 1 / (1 + sensitivity * (Ni_prop / threshold)^2)
#     return(max(min(a_t, 1), 0))
#   } else {
#     warning("Unrecognized scenario in a_function")
#     return(1)
#   }
# }

# a_function <- function(Ni, parameters) {
#   Ni_prop <- Ni / parameters$pop_size  # Convert to proportion
# 
#   vigilance_threshold <- 0.002   # e.g. 0.2% threshold
#   vigilance_sensitivity <- 10000  # reduced slope for smoother decay
# 
#   if (Ni_prop <= 0) return(1)
# 
#   a_t <- 1 / (1 + vigilance_sensitivity * (Ni_prop / vigilance_threshold)^2)
#   return(max(min(a_t, 1), 0))
# }






# Utility function
utility_function <- function(a_t, utility_type) {
  if (utility_type == "Log") {
    return(log(a_t) - a_t + 1)
  } else {
    return(-1/2 * (1 - a_t)^2)
  }
}

calculate_Rt <- function(R0, a_t, Ns_prop) {
  return(R0 * a_t^2 * Ns_prop)
}

# calculate the number of infections or recoveries
#get_transitions_stochastic <- function(n, prob){
#  return(sum(rbinom(n = n,size = 1, prob = prob)))
 # cat("Calling rbinom with n =", n, ", prob =", prob, "\n")
 # 
#}

get_transitions_stochastic <- function(n, prob) {
  cat("DEBUG — rbinom() inputs:\n")
  cat("  n =", n, "\n")
  cat("  prob =", prob, "\n")
  
  if (is.na(n) || is.na(prob) || !is.numeric(n) || !is.numeric(prob) || n < 0 || prob < 0 || prob > 1) {
    warning("Invalid rbinom() inputs: returning 0")
    return(0)
  }
  
  return(sum(rbinom(n = n, size = 1, prob = prob)))
}



get_transitions_deterministic_farboodi <- function(n, prob){
  transitions <- (n * prob)
  # if (transitions < 0.5) {
  #  return(0)
  #  } else {
  return(transitions)
  # }
}

get_transitions_deterministic <- function(n, prob){
 # transitions <- (n * prob)
  # if ((n - transitions) < 1) {
  #  transitions = n
  #  }
  if(length(prob)==n){
    transitions=sum(prob)
  } else{
    transitions<-(n*prob)
  }
  return(transitions)
}



# to run the SIRD kernel with a for loop and binomial transistions
#y = initial_state; times = time_pre_shock; func = sir_costate_model; parms = parameters
run_sir_binomial <- function(initial_state, 
                             times, 
                             parameters,
                             bool_stochastic_beta = FALSE,
                             update_function = get_transitions_stochastic){ # note: the default is get_infections_determistic()
  
  # copy initial states
  states <- data.frame(t(initial_state))
  
  # convert population fractions into numbers
  states[grepl('N.',names(states))] <- round(states[grepl('N.',names(states))] * parameters$pop_size)
  
  # make sure the total population size is fixed (by changing Ns if needed)
  states[grepl('Ns',names(states))] <- parameters$pop_size - sum(states[grepl('N.',names(states)) & !grepl('Ns',names(states))])
  
  # set summary data.table
  states_out        <- matrix(NA, nrow=length(times), ncol=length(states))
  #colnames(states_out) <- names(states)
  colnames(states_out) <- c("Ns", "Ni", "Nr", "Nd", 
                            "HealthCost", "SocialActivityCost", "TotalCost",
                            "a_t", "u_t", "Rt")

  
  
  dim(states_out)
  
  # start with initial states.  
  # note: use single values as parameter is much faster in the iterative process
  names(states)
  Ns <- states$Ns
  Ni <- states$Ni
  Nr <- states$Nr
  Nd <- states$Nd
 # lambda_s <- 0
#  lambda_i <- 0
  HealthCost <- 0
  SocialActivityCost <- 0
  
  
  
  # compute some results only once
  fx_per_capita <- parameters$fx / parameters$pop_size
  rho_plus_delta <- parameters$rho + parameters$delta
  a_t_prev <- 1
  # run over times (excl the first one)
  for(i_day in times[-1]){
    
    if(bool_stochastic_beta){
      # random noise on beta
      beta_t <- max(0.0001,  rnorm(1, mean = parameters$beta, sd =parameters$sigma))#1
    } else{
      beta_t <- parameters$beta
    }
    
    # Calculate optimal action based on utility type
 #   a_t <- a_function(Ni, parameters$alpha, beta_t, 
    #                  Ns/parameters$pop_size, Ni/parameters$pop_size, 
     #                 lambda_s, lambda_i, 
      #                parameters$utility_type, parameters$scenario, parameters$tolerance)#Lambda_s, Lambda_i
  #  a_t <- a_function(Rt, parameters)
    
    # Initialize a_t_prev = 1 before the loop (outside)
     # if (i_day == 1) {
     #   a_t_prev <- 1
     # }
    
    Rt <- calculate_Rt(parameters$R0, a_t_prev, Ns / parameters$pop_size)
    Ns_prop <- Ns / parameters$pop_size
    a_t <- a_function(Ni, Ns, parameters)
#Rt parameters Ni Rt a_t_prev
    a_t_prev <- a_t  # update for next time step
    
    
  #  cat("DEBUG -- Ni =", Ni, "Ns =", Ns, "a_t =", a_t, "beta_t =", beta_t, "\n")
  #  cat("Calculated prob =", beta_t * a_t^2 * Ni / parameters$pop_size, "\n")
    
    # Calculate utility of action
    u_t <- utility_function(a_t, parameters$utility_type)
    
    # Calculate transitions
    new_infections <- update_function(Ns, prob = beta_t * a_t^2 * Ni/parameters$pop_size)
    new_recoveries <- update_function(Ni, prob = parameters$gamma)
    new_death      <- update_function(new_recoveries, prob = parameters$pi)
    
    if((Ni - new_recoveries) < parameters$infect_thres){
      new_recoveries = Ni
      new_infections = 0
    }
    
    # get health transitions
    dNs <- -new_infections
    dNi <- new_infections - new_recoveries
    dNr <- new_recoveries - new_death
    dNd <- new_death
    
    # calculate change in lambda 
    #  dLambda_s <- rho_plus_delta * Lambda_s - 
    #   (u_t + (Lambda_i - Lambda_s) *  beta_t * a_t^2 * Ni/parameters$pop_size)
    # dLambda_i <- rho_plus_delta * Lambda_i - 
    #  (u_t + parameters$alpha * (Lambda_i - Lambda_s) *  beta_t * a_t^2 * Ns/parameters$pop_size) - 
    #  parameters$gamma * (parameters$kappa + Lambda_i)
  
      
  #  } else if (parameters$scenario == "optimal-policy") {
   #   d_lambda_s <- rho_plus_delta * lambda_s - 
  #      u_t - (lambda_i - lambda_s) * beta_t * a_t^2 * Ni / parameters$pop_size
  #    d_lambda_i <- rho_plus_delta * lambda_i - 
   #     u_t - (lambda_i - lambda_s) * beta_t * a_t^2 * Ns / parameters$pop_size - 
    #    parameters$gamma * (parameters$kappa + lambda_i)
    #}
    
    
    # get current costs (per capita)
    HealthCost <- HealthCost + fx_per_capita * exp(-rho_plus_delta * i_day)  *  Ni*parameters$kappa/parameters$pop_size
    SocialActivityCost <- SocialActivityCost + fx_per_capita * exp(-rho_plus_delta * i_day) * (Ns + Ni) * abs(u_t)/parameters$pop_size
    
    Rt <- calculate_Rt(parameters$R0, a_t, Ns/parameters$pop_size) 
 
    
    
    # Update states
    # note: this needs to be done at the end, otherwise it affects the Lambda calculation
    Ns <- Ns + dNs
    Ni <- Ni + dNi
    Nr <- Nr + dNr
    Nd <- Nd + dNd
    #  Lambda_s <- Lambda_s + dLambda_s
    #  Lambda_i <- Lambda_i + dLambda_i
    #lambda_s <- lambda_s+d_lambda_s
    #lambda_i <- lambda_i+d_lambda_i
    # keep track of the states
    states_out[i_day+1,] = c(Ns, Ni, Nr, Nd,
                             #lambda_s, lambda_i,
                             HealthCost, SocialActivityCost, 
                             HealthCost + SocialActivityCost,
                             a_t, u_t, Rt)#Lambda_s, Lamba_i
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
    exp_fade_out   <- (output_summary$Nr+output_summary$Nr) < fadeout_threshold
    output_summary <- output_summary[!exp_fade_out,]
    output_all <- output_all[!exp_fade_out,,]
    plot_tag <- paste(plot_tag,'(excl. fadeout)')
  } else {
    plot_tag <- paste(plot_tag,'(all)')
  }
  # Plot R(t) in place of lambda_s
  #plot(times, R_t, type = "l", col = "blue", lwd = 2,
  #    xlab = "Time (days)", ylab = "Effective Reproduction Number R(t)",
  #    main = "Effective Reproduction Number Over Time")
  #abline(h = 1, col = "red", lty = 2)  # Threshold of R(t) = 1
  
  # set figure sub panels
  par(mfrow=c(3,4))
  
  # explore costs
  output_cost <- output_summary[,grepl('Cost',names(output_summary))]
  y_lim <- range(output_cost)
  boxplot(output_cost,las=2,ylim=y_lim,main=plot_tag)
  points(colMeans(output_cost),pch=8) # mean
  points((1:3)+0.2,output_summary_deterministic[names(output_cost)],col=4,pch=8,lwd=3) # mean
  
  # plot activity vs health cost 
  plot(output_summary$SocialActivityCost,
       output_summary$HealthCost,
       xlab='SocialActivityCost',
       ylab='HealthCost',
       main='health vs activity')
  points(output_summary_deterministic$SocialActivityCost,
         output_summary_deterministic$HealthCost,
         col=4,pch=8,lwd=3) # mean
  
  
  
  # explore states
  sel_states <- names(initial_state)[!grepl('Total',names(initial_state)) & !grepl('Ns',names(initial_state))]
  for(i_state in sel_states){
    y_lim <- range(0,output_all[,,i_state],output_sim_deterministic[,i_state],na.rm = TRUE)
    plot(output_sim_deterministic[,i_state], col=1, main=i_state, ylab = i_state, xlab = 'time', ylim = y_lim, type='l', lwd=2) # infections, deterministic
    for(i_exp in 1:dim(output_all)[1]){
      output_sim <- data.frame(output_all[i_exp,,i_state])
      lines(output_sim,col=alpha(2,0.3)) # infections, stochastic
    }
    lines(output_sim_deterministic[,i_state],col=1,lwd=2) # add line again
  }
  
  # prepare summary statistics
  print_out <- data.frame(output=c('Health Cost (per capita)','Social Activity Cost (per capita)','Total Cost (per capita)'),
                          deterministic = get_summary_stats(output_sim_deterministic[nrow(output_sim_deterministic),]),
                          stochastic = get_summary_stats(output_summary))
  names(print_out)[1] <- plot_tag
  
  # print summary statistics to console
  print(paste('SUMMARY STATISTICS:',plot_tag))
  print(print_out)
  
} # end function



# to run multiple model realisations
run_experiments <- function(initial_state, times, parameters, 
                            bool_stochastic_beta, update_function, 
                            num_experiments){
  
  # set random number generator seed
  set.seed(parameters$rng_seed)
  
  # init result matrix
  output_summary <- data.frame(matrix(NA,nrow=num_experiments,ncol=length(initial_state)))
  names(output_summary) <- names(initial_state)
  dim(output_summary)
  
  # init array to capture all results
  output_all <- array(NA,dim=c(num_experiments,length(times),length(initial_state)),dimnames=list(NULL,NULL,names(initial_state)))
  dim(output_all)
  
  # run multiple model realisations and store results
  for(i_exp in 1:num_experiments){
    print(i_exp)
    # run model
    output_sim <- run_sir_binomial(initial_state = initial_state, 
                                   times = times, 
                                   parameters = parameters,
                                   bool_stochastic_beta = bool_stochastic_beta,
                                   update_function = update_function)
    # store results
    output_summary[i_exp,] <- output_sim[nrow(output_sim),]
    output_all[i_exp,,] <- as.matrix(output_sim)
  }
  
  return(list(output_summary=output_summary,
              output_all=output_all))
}

#Alternative optimal activity calculations


# a_function <- function(Rt) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
# 
#   x <- abs(Rt - 1)  # quadratic distance from 1
# 
#   if (Rt > 1) {
#     a_t <- 1 - x
#   } else {
#     a_t <- 1 + x
#   }
# 
#   return(max(0, min(1, a_t)))  # Clamp to [0, 1]
# }





# a_function <- function(Rt) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
#   
#   x <- (Rt - 1)^2  # quadratic distance from 1
#   
#   if (Rt > 1) {
#     a_t <- 1 - x
#   } else {
#     a_t <- 1 + x
#   }
#   
#   return(max(0, min(1, a_t)))  # Clamp to [0, 1]
# }

# a_function <- function(Rt) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
#   
#   deviation <- Rt - 1
#   sensitivity <- 3  # you can tune 2, 3, 5, depending on how strict you want
#   
#   adjusted_deviation <- deviation^2 * sign(deviation)  # quadratic + directional
#   adjustment <- tanh(sensitivity * adjusted_deviation)
#   
#   a_t <- 1 - 0.5 * adjustment
#   
#   return(a_t)
# }

# a_function <- function(Rt) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
#   
#   deviation <- Rt - 1
#   sensitivity <- 3  # sensitivity tuning
#   
#   # Quadratic adjustment
#   adjusted_deviation <- deviation^2 * sign(deviation)
#   adjustment <- tanh(sensitivity * adjusted_deviation)
#   
#   if (Rt <= 1) {
#     a_t <- 1  # Full normal activity if Rt <= 1
#   } else {
#     a_t <- 1 - 0.5 * adjustment  # Reduce activity if Rt > 1
#   }
#   
#   return(a_t)  # Always between 0 and 1 naturally
# }
# a_function <- function(Rt) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
#   
#   sensitivity <- 2.5  # tuned softer
#   threshold <- 1
#   
#   if (Rt <= threshold) {
#     a_t <- 1
#   } else {
#     penalty <- 1 / (1 + exp(-sensitivity * (threshold - Rt)))
#     a_t <- 1 - 0.4 * penalty  # smaller multiplier
#   }
#   
#   return(a_t)
# }

# a_function <- function(Rt) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
#   
#   deviation <- Rt - 1
#   sensitivity <- 1.5  # important: soft
#   max_penalty <- 0.5  # max 50% activity reduction
#   
#   adjusted_deviation <- deviation^2 * sign(deviation)
#   penalty <- tanh(sensitivity * adjusted_deviation)
#   
#   if (Rt <= 1) {
#     a_t <- 1  # Full activity if Rt <= 1
#   } else {
#     a_t <- 1 - max_penalty * penalty
#   }
#   
#   return(a_t)
# }
######################################################

# 
# a_function <- function(Ni, parameters) {
#   if (is.na(Ni) || is.na(parameters$pop_size)) {
#     warning("Invalid Ni or pop_size")
#     return(0)
#   }
#   Ni_prop <- Ni / parameters$pop_size  # Proportion of infectives
#   threshold <- 0.01                    # 1% prevalence triggers concern
#   sensitivity <- 200                  # Controls steepness of response
# 
#   # Logistic-style vigilance
#   a_t <- 1 / (1 + sensitivity * max(0, Ni_prop - threshold))
# 
#   return(max(0, min(1, a_t)))
# }

#

# a_function <- function(Rt, parameters) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
#   
#   # Customizable parameters
#   midpoint <- 1             # Where change is neutral
#   base_level <- 0.5         # Activity at Rt = 1
#   max_change <- 0.5         # How far activity can shift from base_level
#   steepness <- 3            # Controls how sharp the curve is around Rt = 1
#   
#   # Sigmoid-like function centered at Rt = 1
#   # Positive shift below 1, negative shift above 1
#   delta <- Rt - midpoint
#   direction <- -tanh(steepness * delta)  # tanh gives smooth transition [-1, 1]
#   
#   # Compute activity
#   a_t <- base_level + max_change * direction
#   
#   # Clamp to [0, 1]
#   return(max(0, min(1, a_t)))
# }






# a_function <- function(Rt) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
#   
#   deviation <- Rt - 1
#   sensitivity <- 3  # sensitivity tuning
#   
#   # Quadratic adjustment
#   adjusted_deviation <- deviation^2 * sign(deviation)
#   adjustment <- tanh(sensitivity * adjusted_deviation)
#   
#   if (Rt <= 1) {
#     a_t <- 1  # Full normal activity if Rt <= 1
#   } else {
#     a_t <- 1 - 0.5 * adjustment  # Reduce activity if Rt > 1
#   }
#   
#   return(a_t)  # Always between 0 and 1 naturally
# }



# a_function <- function(R0, Ns_prop) {
#  if (is.na(R0) || is.na(Ns_prop) || R0 <= 0 || Ns_prop <= 0) {
#    warning("Invalid R0 or Ns_prop")
#    return(0)
#  }
#  
#  a_t <- sqrt(1 / (R0 * Ns_prop))
#  
#  return(max(0, min(1, a_t)))  # Ensure stays between 0 and 1
#}









# a_function <- function(Rt, parameters) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
# 
#   # Parameters for logistic shape
#   midpoint <- 1       # Where activity drops to 50%
#   steepness <- 5      # Controls how quickly activity falls
# 
#   a_t <- 1 / (1 + exp(steepness * (Rt - midpoint)))
# 
#   return(max(0, min(1, a_t)))  # Clamp to [0, 1]
# }

# a_function <- function(Rt, parameters) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
#   
#   Rt_min <- 0.75   # below this, full activity
#   Rt_max <- 1.25   # above this, minimum activity
#   
#   # Linear interpolation
#   if (Rt <= Rt_min) {
#     a_t <- 1
#   } else if (Rt >= Rt_max) {
#     a_t <- 0
#   } else {
#     a_t <- 1 - (Rt - Rt_min) / (Rt_max - Rt_min)
#   }
#   
#   return(max(0, min(1, a_t)))  # Clamp to [0, 1]
# }


# a_function <- function(Rt, parameters) {
#  # if (is.na(Rt)) {
#   #  warning("Invalid Rt")
#  #   return(0)
#  # }
# 
#   sensitivity <- 0.6  # similar to infections 50
#   a_t <- 1 / (1 + sensitivity * Rt)
# 
#   return(max(0, min(1, a_t)))
# }

# a_function <- function(parameters) {
#   R0 <- parameters$R0
#   if (is.na(R0)) {
#     warning("Invalid R0")
#     return(0)
#   }
# 
#   k <- 0.5  # Caution multiplier
#   a_t <- 1 / (1 + k * max(0, R0 - 1))
# 
#   return(max(0, min(1, a_t)))
# }

# a_function <- function(Rt, parameters) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
# 
#   threshold <- 1
#   sensitivity <- 2.5  # between 1–5
# 
#   if (Rt <= threshold) {
#     a_t <- 1
#   } else {
#     a_t <- 1 / (1 + sensitivity * (Rt - threshold))
#   }
# 
#   return(max(0, min(1, a_t)))  # Clamp to [0, 1]
# }





#a_function <- function(Rt, parameters) {
# sensitivity <- 0.6  # tuned to mirror the old Ni-based sensitivity of 50
#  a_t <- 1 / (1 + sensitivity * Rt)
# return(a_t)
#}

# a_function <- function(Rt) {
#   if (is.na(Rt)) {
#     warning("Invalid Rt")
#     return(0)
#   }
#   
#   sensitivity <- 1.5   # much softer now
#   threshold <- 1
#   
#   if (Rt <= threshold) {
#     a_t <- 1
#   } else {
#     penalty <- 1 / (1 + exp(-sensitivity * (threshold - Rt)))
#     a_t <- 1 - 0.2 * penalty  # SMALLER max penalty
#   }
#   
#   return(a_t)
# }
# 



# # Myopic altruistic activity function (laissez-faire equilibrium)
# a_function <- function(Ni, parameters) {#
#   Ni_prop <- Ni / parameters$pop_size
#   Ns_prop <- parameters$ns0  # Assume constant or use from state if dynamic
#   epsilon <- 1e-12
#   Ni_prop <- max(Ni_prop, epsilon)
# 
#   # Compute multiplier
#   multiplier <- parameters$beta * Ni_prop * parameters$kappa * (1 + parameters$alpha * Ns_prop)
# 
#   sqrt_term <- sqrt(1 + 4 * multiplier)
#   a_t <- (-1 + sqrt_term) / (2 * multiplier)
# 
#  #  Clamp a_t between [0,1]
#  return (a_t)
# 
#   #(max(0, min(1, a_t)))
# }



# if (utility_type == "Log") {
#  if (scenario == "laissez-faire") {
#  sqrt_arg <- (Ns + Ni) + 4 * (1 + alpha) * beta * Ni * Ns * denom
#  return((-(Ns + Ni) + sqrt(sqrt_arg)) / (2 * (1 + alpha) * beta * Ni * Ns * denom))
# } else if (scenario == "optimal-policy") {
#  sqrt_arg <- (Ns + Ni)^2 + 8 * beta * Ns * Ni * denom * (Ns + Ni)
#    return((-(Ns + Ni) + sqrt(sqrt_arg)) / (4 * beta * Ns * Ni * denom))
#   }


# } else if (utility_type == "Quadratic") {
#    return( (Ns + Ni) / (Ns + Ni + (1 + alpha) * beta * Ni * Ns * denom) )
# }
#}

#dlambda's



# Derivative of mu_s 
#  d_mu_s <- rho_plus_delta * mu_s - u_t - (mu_i - mu_s) * beta_t * a_t^2 * Ni/parameters$pop_size
#  R0 <- parameters$beta / parameters$gamma

# Function to calculate R(t)
#  calculate_Rt <- function(R0, a_t, Ns) {
#    return(R0 * a_t^2 * Ns)
# }
# Derivative of mu_i 
# d_mu_i <- rho_plus_delta * mu_i - u_t - (mu_i - mu_s) * beta_t * a_t^2 * Ns/parameters$pop_size + parameters$gamma * (parameters$kappa + mu_i)
# if (parameters$scenario == "laissez-faire") {
#  d_lambda_s <- rho_plus_delta * lambda_s - 
#   (u_t + (lambda_i - lambda_s) * beta_t * a_t^2 * Ni / parameters$pop_size)
#  d_lambda_i <- rho_plus_delta * lambda_i - 
#   (u_t + parameters$alpha * (lambda_i - lambda_s) * beta_t * a_t^2 * Ns / parameters$pop_size) + 
#    parameters$gamma * (parameters$kappa + lambda_i)

#d_lambda_s <- rho_plus_delta * lambda_s - 
#  u_t - (1 + parameters$alpha) * (lambda_i - lambda_s) * beta_t * a_t^2 * Ni / parameters$pop_size
#d_lambda_i <- rho_plus_delta * lambda_i - 
#  u_t - (1 + parameters$alpha) * (lambda_i - lambda_s) * beta_t * a_t^2 * Ns / parameters$pop_size + 
#   parameters$gamma * (parameters$kappa + lambda_i)

# Corrected equation for d_lambda_s
#   d_lambda_s <- rho_plus_delta * lambda_s - 
#      u_t - parameters$alpha * (lambda_i - lambda_s) * beta_t * a_t^2 * Ni / parameters$pop_size

# Corrected equation for d_lambda_i
#   d_lambda_i <- rho_plus_delta * lambda_i - 
#    u_t - parameters$alpha * (lambda_i - lambda_s) * beta_t * a_t^2 * Ns / parameters$pop_size + 
#   parameters$gamma * (parameters$kappa + lambda_i)


