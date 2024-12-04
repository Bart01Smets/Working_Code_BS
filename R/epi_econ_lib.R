#
# Help function for the EPI-ECON modelling work
#
# - Deterministic and stochastic transitions
################################################################ #

library(confintr)
library(scales)

# Function to calculate the optimal action (a_t)
a_function <- function(alpha, beta, Ns, Ni, Lambda_s, Lambda_i, utility_type, tolerance) {
  denom <- max(abs(Lambda_s - Lambda_i), tolerance)  # Avoid division by zero
  
  if(Ni == 0){
    return(1)
  }
  
  if (utility_type == "Log") {
    sqrt_arg <- (Ns + Ni) + 4 * (1 + alpha) * beta * Ni * Ns * denom
    return( (-(Ns + Ni) + sqrt(sqrt_arg)) / (2 * (1 + alpha) * beta * Ni * Ns * denom) )
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

# calculate the number of infections or recoveries
get_transitions_stochastic <- function(n, prob){
  return(sum(rbinom(n = n,size = 1, prob = prob)))
}

get_transitions_deterministic <- function(n, prob){
  return((n * prob))
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
  states_out        <- data.frame(matrix(NA, nrow=length(times), ncol=length(states)))
  names(states_out) <- names(states)
  dim(states_out)
  
  # start with initial states
  states_out[1,]    <- states
  dim(states_out)
  
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
                      states$Ns/parameters$pop_size, states$Ni/parameters$pop_size, 
                      states$Lambda_s, states$Lambda_i, 
                      parameters$utility_type, parameters$tolerance)
    
    # Calculate utility of action
    u_t <- utility_function(a_t, parameters$utility_type)
    
    # Calculate transitions
    new_infections <- update_function(states$Ns, prob = beta_t * a_t^2 * states$Ni/parameters$pop_size)
    new_recoveries <- update_function(states$Ni, prob = parameters$gamma)
    new_death      <- update_function(new_recoveries, prob = parameters$pi)
    
    # get health transitions
    dNs <- -new_infections
    dNi <- new_infections - new_recoveries
    dNr <- new_recoveries - new_death
    dNd <- new_death
    
    # calculate change in lambda 
    dLambda_s <- (parameters$rho + parameters$delta) * states$Lambda_s - 
      (u_t + (states$Lambda_i - states$Lambda_s) *  beta_t * a_t^2 * states$Ni/parameters$pop_size)
    dLambda_i <- (parameters$rho + parameters$delta) * states$Lambda_i - 
      (u_t + parameters$alpha * (states$Lambda_i - states$Lambda_s) *  beta_t * a_t^2 * states$Ns/parameters$pop_size) - 
      parameters$gamma * (parameters$kappa + states$Lambda_i)
    
    # get current costs (per capita)
    states$HealthCost <- states$HealthCost + (parameters$fx / (states$Ns + states$Ni + states$Nr + states$Nd)) * exp(-(parameters$rho + parameters$delta) * i_day) * parameters$gamma * parameters$kappa * states$Ni
    states$SocialActivityCost <- states$SocialActivityCost + (parameters$fx / (states$Ns + states$Ni + states$Nr + states$Nd)) * exp(-(parameters$rho + parameters$delta) * i_day) * (states$Ns + states$Ni) * abs(u_t)
    states$TotalCost <- states$HealthCost + states$SocialActivityCost
    
    # Update states
    # note: this needs to be done at the end, otherwise it affects the Lambda calculation
    states$Ns <- states$Ns + dNs
    states$Ni <- states$Ni + dNi
    states$Nr <- states$Nr + dNr
    states$Nd <- states$Nd + dNd
    states$Lambda_s <- states$Lambda_s + dLambda_s
    states$Lambda_i <- states$Lambda_i + dLambda_i
    
    states$a_t <- a_t
    states$u_t <- u_t
    
    # keep track of the states
    states_out[i_day+1,] = states
  }
  
  return(states_out)
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
                               plot_tag, bool_excl_fadeout = FALSE){
  
  # split output
  output_summary <- output_experiments$output_summary
  output_all     <- output_experiments$output_all
  
  # get final states of deterministic model
  output_summary_deterministic <- output_sim_deterministic[nrow(output_sim_deterministic),]
  
  # identify simulations with stochastic fade out 
  if(bool_excl_fadeout){
    exp_fade_out <- output_summary$Ni == 0
    output_summary <- output_summary[!exp_fade_out,]
    output_all <- output_all[!exp_fade_out,,]
    plot_tag <- paste(plot_tag,'(excl. fadeout)')
  } else {
    plot_tag <- paste(plot_tag,'(all)')
  }
  
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
  sel_states <- names(initial_state)[!grepl('Total',names(initial_state))]
  for(i_state in sel_states){
    y_lim <- range(0,output_all[,,i_state],na.rm = TRUE)
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


