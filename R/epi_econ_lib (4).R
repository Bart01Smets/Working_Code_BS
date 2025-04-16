
#
# Help function for the EPI-ECON modelling work
#
# - Deterministic and stochastic transitions
################################################################ #
#setwd("C:/Users/Bart Smets/OneDrive/Documenten/GitHub/Working_Code_BS/R")
library(confintr)
library(scales)

# Function to calculate the optimal action (a_t)
a_function <- function(alpha, beta, Ns, Ni, lambda_s, lambda_i, utility_type, scenario, tolerance) {
  denom <- max(abs(lambda_s - lambda_i), tolerance)  # Avoid division by zero#Lambda_s -Lambda_i
  
  if(Ni == 0){
    return(1)
  }

  # if (utility_type == "Log") {
  #  sqrt_arg <- (Ns + Ni)^2 + 8 * beta * Ns * Ni * denom * (Ns + Ni)
  # return((- (Ns + Ni) + sqrt(sqrt_arg)) / (4 * beta * Ns * Ni * denom))
  
  #if  #sqrt_arg <- (Ns + Ni) + 4 * (1 + alpha) * beta * Ni * Ns * denom
  #return( (-(Ns + Ni) + sqrt(sqrt_arg)) / (2 * (1 + alpha) * beta * Ni * Ns * denom) )
  
  if (utility_type == "Log") {
    if (scenario == "laissez-faire") {
      sqrt_arg <- (Ns + Ni) + 4 * (1 + alpha) * beta * Ni * Ns * denom
      return((-(Ns + Ni) + sqrt(sqrt_arg)) / (2 * (1 + alpha) * beta * Ni * Ns * denom))
    } else if (scenario == "optimal-policy") {
      sqrt_arg <- (Ns + Ni)^2 + 8 * beta * Ns * Ni * denom * (Ns + Ni)
      return((-(Ns + Ni) + sqrt(sqrt_arg)) / (4 * beta * Ns * Ni * denom))
    }
    
    
  } else if (utility_type == "Quadratic") {
    return( (Ns + Ni) / (Ns + Ni + (1 + alpha) * beta * Ni * Ns * denom) )
  }
}

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
get_transitions_stochastic <- function(n, prob){
  return(sum(rbinom(n = n,size = 1, prob = prob)))
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
  transitions <- (n * prob)
  # if ((n - transitions) < 1) {
   #  transitions = n
  #  } 
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
  colnames(states_out) <- names(states)
  dim(states_out)
  
  # start with initial states.  
  # note: use single values as parameter is much faster in the iterative process
  names(states)
  Ns <- states$Ns
  Ni <- states$Ni
  Nr <- states$Nr
  Nd <- states$Nd
  lambda_s <- 0
  lambda_i <- 0
  HealthCost <- 0
  SocialActivityCost <- 0
  
  
  
  # compute some results only once
  fx_per_capita <- parameters$fx / parameters$pop_size
  rho_plus_delta <- parameters$rho + parameters$delta
  
  # run over times (excl the first one)
  for(i_day in times[-1]){
    
    if(bool_stochastic_beta){
      # random noise on beta
      beta_t <- max(0.0001,  rnorm(1, mean = parameters$beta, sd =parameters$sigma))#1
    } else{
      beta_t <- parameters$beta
    }
    
    # Calculate optimal action based on utility type
    a_t <- a_function(parameters$alpha, beta_t, 
                      Ns/parameters$pop_size, Ni/parameters$pop_size, 
                      lambda_s, lambda_i, 
                      parameters$utility_type, parameters$scenario, parameters$tolerance)#Lambda_s, Lambda_i
    
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
    
    
    # Derivative of mu_s 
    #  d_mu_s <- rho_plus_delta * mu_s - u_t - (mu_i - mu_s) * beta_t * a_t^2 * Ni/parameters$pop_size
    #  R0 <- parameters$beta / parameters$gamma
    
    # Function to calculate R(t)
    #  calculate_Rt <- function(R0, a_t, Ns) {
    #    return(R0 * a_t^2 * Ns)
    # }
    # Derivative of mu_i 
    # d_mu_i <- rho_plus_delta * mu_i - u_t - (mu_i - mu_s) * beta_t * a_t^2 * Ns/parameters$pop_size + parameters$gamma * (parameters$kappa + mu_i)
    if (parameters$scenario == "laissez-faire") {
      d_lambda_s <- rho_plus_delta * lambda_s - 
        (u_t + (lambda_i - lambda_s) * beta_t * a_t^2 * Ni / parameters$pop_size)
      d_lambda_i <- rho_plus_delta * lambda_i - 
        (u_t + parameters$alpha * (lambda_i - lambda_s) * beta_t * a_t^2 * Ns / parameters$pop_size) - 
        parameters$gamma * (parameters$kappa + lambda_i)
      
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
      
      
    } else if (parameters$scenario == "optimal-policy") {
      d_lambda_s <- rho_plus_delta * lambda_s - 
        u_t - (lambda_i - lambda_s) * beta_t * a_t^2 * Ni / parameters$pop_size
      d_lambda_i <- rho_plus_delta * lambda_i - 
        u_t - (lambda_i - lambda_s) * beta_t * a_t^2 * Ns / parameters$pop_size + 
        parameters$gamma * (parameters$kappa + lambda_i)
    }
    
    
    # get current costs (per capita)
    HealthCost <- HealthCost + fx_per_capita * exp(-rho_plus_delta * i_day) * parameters$gamma * parameters$kappa * Ni
    SocialActivityCost <- SocialActivityCost + fx_per_capita * exp(-rho_plus_delta * i_day) * (Ns + Ni) * abs(u_t)
    
    Rt <- calculate_Rt(parameters$R0, a_t, Ns/parameters$pop_size) 
    
    
    # Update states
    # note: this needs to be done at the end, otherwise it affects the Lambda calculation
    Ns <- Ns + dNs
    Ni <- Ni + dNi
    Nr <- Nr + dNr
    Nd <- Nd + dNd
    #  Lambda_s <- Lambda_s + dLambda_s
    #  Lambda_i <- Lambda_i + dLambda_i
    lambda_s <- lambda_s+d_lambda_s
    lambda_i <- lambda_i+d_lambda_i
    # keep track of the states
    states_out[i_day+1,] = c(Ns, Ni, Nr, Nd,
                             lambda_s, lambda_i,
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
