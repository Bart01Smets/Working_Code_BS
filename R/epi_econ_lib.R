#
# Help function for the EPI-ECON modelling work
#
# - Deterministic and stochastic transitions
################################################################ #

library(deSolve)

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

# Define the SIR model function including costates, recovered, dead, and separate costs
sir_costate_model <- function(time, state, parameters) { 
  with(as.list(c(state, parameters)), {
    
    # Extract the state variables
    Ns <- state["Ns"]        
    Ni <- state["Ni"]        
    Nr <- state["Nr"]        
    Nd <- state["Nd"]        
    Lambda_s <- state["Lambda_s"]
    Lambda_i <- state["Lambda_i"]
    HealthCost <- state["HealthCost"]
    SocialActivityCost <- state["SocialActivityCost"]
    TotalCost <- state["TotalCost"]
    
    # Calculate optimal action based on utility type
    a_t <- a_function(alpha, beta, Ns, Ni, Lambda_s, Lambda_i, utility_type, tolerance)
    
    # Calculate utility of action
    u_t <- utility_function(a_t, utility_type)
    
    # Differential equations for SIR model with costates
    dNs <- -beta * a_t^2 * Ns * Ni
    dNi <- beta * a_t^2 * Ns * Ni - gamma * Ni
    dNr <- gamma * Ni * (1 - pi)
    dNd <- gamma * Ni * pi
    dLambda_s <- (rho + delta) * Lambda_s - (u_t + (Lambda_i - Lambda_s) * beta * a_t^2 * Ni)
    dLambda_i <- (rho + delta) * Lambda_i - (u_t + alpha * (Lambda_i - Lambda_s) * beta * a_t^2 * Ns) - gamma * (kappa + Lambda_i)
    
    # Adjusted cost differential equations (per capita)
    dHealthCost <- (fx / (Ns + Ni + Nr + Nd)) * exp(-(rho + delta) * time) * gamma * kappa * Ni
    dSocialActivityCost <- (fx / (Ns + Ni + Nr + Nd)) * exp(-(rho + delta) * time) * (Ns + Ni) * abs(u_t)
    dTotalCost <- dHealthCost + dSocialActivityCost
    
    # Return the rate of change for each state variable
    return(list(c(dNs, dNi, dNr, dNd, dLambda_s, dLambda_i, dHealthCost, dSocialActivityCost, dTotalCost, a_t, u_t)))
  })
}

# Run the SIRD kernel with a for loop (hence, not ODE solver)
# note: 'list' and 'unlist' operations are time consuming
run_sir_costate_model <- function(y , times, func, parms){
  
  states_out        <- data.frame(matrix(NA, ncol=length(y), nrow=length(times)))
  names(states_out) <- names(y)
  states_out[1,]    <- y
  dim(states_out)
  
  for(i in times[-1]){
    states_out[i+1,] <- states_out[i,] + unlist(func(time = i, state = states_out[i,], parameters = parameters))
  }
  
  return(states_out)
}

# calculate the number of infections: deterministic
get_infections_determistic <- function(beta,a_t,Ns,Ni){
  return(beta * a_t * Ns * Ni)
}

# calculate the number of infections: stochastic beta
get_infections_stochastic_beta <- function(beta,a_t,Ns,Ni){
  
  # random noise on beta
  beta_sample <- max(0.0001, beta + parameters$sigma * rnorm(1, mean = 0, sd =beta/2))#1
  
  # return new infections
  return(beta_sample * a_t * Ns * Ni)
}

# alternative to run the SIRD kernel with a for loop (hence, not ODE solver)
#y = initial_state; times = time_pre_shock; func = sir_costate_model; parms = parameters
run_sir_update <- function(initial_state, times, parameters, 
                           infections_function = get_infections_determistic){ # note: the default is get_infections_determistic()

  # copy initial states
  states <- data.frame(t(initial_state))
  
  # set summary data.table
  states_out        <- data.frame(matrix(NA, nrow=length(times), ncol=length(states)))
  names(states_out) <- names(states)
  dim(states_out)

  # start with initial states
  states_out[1,]    <- states
  dim(states_out)
  
  # run over times (excl the first one)
  for(i_day in times[-1]){
      
      # Calculate optimal action based on utility type
      a_t <- a_function(parameters$alpha, parameters$beta, 
                        states$Ns, states$Ni, 
                        states$Lambda_s, states$Lambda_i, 
                        parameters$utility_type, parameters$tolerance)
      
      # Calculate utility of action
      u_t <- utility_function(a_t, parameters$utility_type)
      
      # Calculate transitions
      new_infections <- infections_function(parameters$beta, a_t^2, states$Ns, states$Ni)
      new_recoveries <- parameters$gamma * states$Ni
      
      # get health transitions
      dNs <- -new_infections
      dNi <- new_infections - new_recoveries
      dNr <- new_recoveries * (1 - parameters$pi)
      dNd <- new_recoveries * parameters$pi
      
      # calculate change in lambda 
      dLambda_s <- (parameters$rho + parameters$delta) * states$Lambda_s - 
                      (u_t + (states$Lambda_i - states$Lambda_s) *  parameters$beta * a_t^2 * states$Ni)
      dLambda_i <- (parameters$rho + parameters$delta) * states$Lambda_i - 
                    (u_t + parameters$alpha * (states$Lambda_i - states$Lambda_s) *  parameters$beta * a_t^2 * states$Ns) - 
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

# calculate the number of infections or recovereis
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
    
    # Calculate optimal action based on utility type
    a_t <- a_function(parameters$alpha, parameters$beta, 
                      states$Ns/parameters$pop_size, states$Ni/parameters$pop_size, 
                      states$Lambda_s, states$Lambda_i, 
                      parameters$utility_type, parameters$tolerance)
    
    # Calculate utility of action
    u_t <- utility_function(a_t, parameters$utility_type)
    
    # Calculate transitions
    new_infections <- update_function(states$Ns, prob = parameters$beta * a_t^2 * states$Ni/parameters$pop_size)
    new_recoveries <- update_function(states$Ni, prob = parameters$gamma)
    new_death      <- update_function(new_recoveries, prob = parameters$pi)
    
    # get health transitions
    dNs <- -new_infections
    dNi <- new_infections - new_recoveries
    dNr <- new_recoveries - new_death
    dNd <- new_death
    
    # calculate change in lambda 
    dLambda_s <- (parameters$rho + parameters$delta) * states$Lambda_s - 
      (u_t + (states$Lambda_i - states$Lambda_s) *  parameters$beta * a_t^2 * states$Ni/parameters$pop_size)
    dLambda_i <- (parameters$rho + parameters$delta) * states$Lambda_i - 
      (u_t + parameters$alpha * (states$Lambda_i - states$Lambda_s) *  parameters$beta * a_t^2 * states$Ns/parameters$pop_size) - 
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