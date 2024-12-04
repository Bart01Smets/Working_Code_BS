#
# Help function for the EPI-ECON modelling work
#
# - Deterministic and stochastic transitions
################################################################ #

library(confintr)

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
compare_sim_output <- function(output_sim, output_experiments, output_deterministic, plot_tag){
  
  # identify simulations with stochastic fade out 
  output_experiments$bool_fade_out <- output_experiments$Ni == 0
  
  # some graphical exploration
  par(mfrow=c(2,2)) # reset sub-panels
  boxplot(output_experiments$SocialActivityCost, ylab='SocialActivityCost', main=paste(plot_tag,'(all)'))
  points(1,mean(output_experiments$SocialActivityCost),pch=8) # mean
  abline(h=output_sim_deterministic$SocialActivityCost[nrow(output_sim_deterministic)],col=4)
  text(x=0.6,y=output_sim_deterministic$SocialActivityCost[nrow(output_sim_deterministic)],
       labels=('deterministic'), pos=3,col=4)
  
  # some graphical exploration (excl stochastic fade out)
  boxplot(output_experiments$SocialActivityCost[!output_experiments$bool_fade_out], ylab='SocialActivityCost', main=paste(plot_tag,'(excl fadeout)'))
  points(1,mean(output_experiments$SocialActivityCost[!output_experiments$bool_fade_out]),pch=8) # mean
  abline(h=output_sim_deterministic$SocialActivityCost[nrow(output_sim_deterministic)],col=4)
  text(x=0.6,y=output_sim_deterministic$SocialActivityCost[nrow(output_sim_deterministic)],
       labels=('deterministic'), pos=3,col=4)
  
  # explore last simulation: infections
  plot(output_sim$Ni,col=2,main='Infections',type='l',lwd=2) # infections, stochastic
  lines(output_sim_deterministic$Ni,col=1) # infections, stochastic
  
  # explore last simulation: all health states
  plot(output_sim$Ns,main='Health states',ylim=c(0,parameters$pop_size),type='l')
  lines(output_sim$Ni,col=2)
  lines(output_sim$Nr,col=3)
  lines(output_sim$Nd,col=4)
  legend('topright',
         c('Ns','Ni','Nr','Nd'),
         col = 1:4,
         lwd = 2,
         ncol = 4,
         cex=0.5)
  
  # Print final stats
  print(paste('SUMMARY STATISTICS:',plot_tag))
  print(data.frame(output=c('Health Cost (per capita)','Social Activity Cost (per capita)','Total Cost (per capita)'),
             determ=get_summary_stats(output_sim_deterministic[nrow(output_sim_deterministic),]),
             stochastic_all = get_summary_stats(output_experiments),
             stochstic_select = get_summary_stats(output_experiments[!output_experiments$bool_fade_out,])))
  
  }

