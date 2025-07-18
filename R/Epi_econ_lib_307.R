
#
# Help function for the EPI-ECON modelling work
#
# - Deterministic and stochastic transitions
################################################################ #

library(confintr)
library(scales)


a_function <- function(Ni, Ns, parameters) {
  Ni_prop <- Ni / parameters$pop_size
  Ns_prop <- Ns / parameters$pop_size
 
  multiplier <- parameters$beta *Ni_prop*Ns_prop*parameters$kappa#* (1 + parameters$alpha* Ns_prop)#*Ns_prop Ni_prop parameters$kappa#Ndparameters$v* Ni_prop
  sqrt_term <- sqrt(1 + 8 * multiplier) #_prop
  a_t <- (-1 + sqrt_term) / (4 * multiplier)
  
 
  if (Ni < 1e-10 || Ns < 1e-10) return(1)
  
  return(max(0, min(1, a_t)))
}


# Utility function
utility_function <- function(a_t, utility_type) {
  if (utility_type == "Log") {
    return(log(a_t) - a_t + 1)
  } else {
    return(-1/2 * (1 - a_t)^2)
  }
}

calculate_Rt <- function(R0, a_t, Ns_prop, Ni) {
  if (Ni <= 0) return(0)
  return(R0 * a_t^2 * Ns_prop)
}

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
  
  return(transitions)
  # }
}

get_transitions_deterministic <- function(n, prob){
 
  if(length(prob)==n){
    transitions=sum(prob)
  } else{
    transitions<-(n*prob)
  }
  return(transitions)
}

# to run the SIRD kernel with a for loop and binomial transistions

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
    Rt <- calculate_Rt(parameters$R0, a_t_prev, Ns / parameters$pop_size,  Ni)
    Ns_prop <- Ns / parameters$pop_size
    a_t <- a_function(Ni, Ns, parameters)

    a_t_prev <- a_t  # update for next time step
    
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
    
    # get current costs (per capita)
    scale_factor <- 2132
    HealthCost <-  HealthCost+ scale_factor* fx_per_capita * exp(-rho_plus_delta * i_day) *  parameters$kappa  * Ni/parameters$pop_size#/parameters$pop_size#HealthCost + parameters$kappaNd/parameters$pop_size*parameters$v
    SocialActivityCost <- SocialActivityCost+ scale_factor* fx_per_capita * exp(-rho_plus_delta * i_day) * (Ns + Ni) * abs(u_t)/parameters$pop_size#/parameters$pop_size#SocialActivityCost +
  
    
    Rt <- calculate_Rt(parameters$R0, a_t, Ns/parameters$pop_size, Ni) 
 
    
    
    # Update states
    # note: this needs to be done at the end, otherwise it affects the Lambda calculation
    Ns <- Ns + dNs
    Ni <- Ni + dNi
    Ni <- max(Ni, 0)
    Nr <- Nr + dNr
    Nd <- Nd + dNd
  
    # keep track of the states
    states_out[i_day+1,] = c(Ns, Ni, Nr, Nd,
                             
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
 
  # Define custom plot labels for selected variables
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
   
    y_vals <- c(0, output_all[,,i_state], output_sim_deterministic[,i_state])
    y_vals <- y_vals[is.finite(y_vals)]  # remove NA, NaN, Inf
    if (length(y_vals) == 0) {
      y_lim <- c(0, 1)  # fallback if all values invalid
    } else {
      y_lim <- range(y_vals)
    }
    
    
    plot_label <- ifelse(i_state %in% names(pretty_labels), pretty_labels[[i_state]], i_state)
    
    plot(output_sim_deterministic[,i_state], 
         col=1, 
         main=plot_label, 
         ylab=plot_label, 
         xlab='time', 
         ylim=y_lim, 
         type='l', 
         lwd=2)
   
    for(i_exp in 1:dim(output_all)[1]){
      output_sim <- data.frame(output_all[i_exp,,i_state])
      lines(output_sim,col=alpha(2,0.3)) # infections, stochastic
    }

    
    lines(output_sim_deterministic[,i_state],col=1,lwd=2) # add line again
  }
 
  boxplot(output_summary$HealthCost,
          main = "Health Cost",
          ylab = "Health Cost",
          ylim = range(c(output_summary$HealthCost, output_summary_deterministic$HealthCost)),
          col = "gray90")  # optional for better contrast
  points(1, mean(output_summary$HealthCost), pch = 8)  # mean (black star)
  points(1.2, output_summary_deterministic$HealthCost, col = "red", pch = 18, cex = 1.5)  # deterministic (blue diamond)
  
  # Social Activity Cost boxplot
  boxplot(output_summary$SocialActivityCost,
          main = "Social Activity Cost",
          ylab = "Social Activity Cost",
          ylim = range(c(output_summary$SocialActivityCost, output_summary_deterministic$SocialActivityCost)),
          col = "gray90")
  points(1, mean(output_summary$SocialActivityCost), pch = 8)
  points(1.2, output_summary_deterministic$SocialActivityCost, col = "red", pch = 18, cex = 1.5)
  
  
  
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
  
  # Run one simulation to get the actual variable names
  temp_output <- run_sir_binomial(initial_state, times, parameters,
                                  bool_stochastic_beta, update_function)

  
  output_all <- array(NA, dim=c(num_experiments, length(times), ncol(temp_output)),
                      dimnames=list(NULL, NULL, names(temp_output)))
  
  # Update summary matrix as well
  output_summary <- data.frame(matrix(NA, nrow=num_experiments, ncol=ncol(temp_output)))
  names(output_summary) <- names(temp_output)
  
  
  
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

