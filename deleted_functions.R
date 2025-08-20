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
# # Myopic altruistic activity function (laissez-faire equilibrium)
# a_function <- function(Ni, Ns, parameters) {#
#   Ni_prop <- Ni / parameters$pop_size
#   Ns_prop <-  Ns/parameters$pop_size#parameters$ns0  # Assume constant or use from state if dynamic
#   epsilon <- 1e-12
#   Ni_prop <- max(Ni_prop, epsilon)
# 
#   # Compute multiplier
#   multiplier <- parameters$beta * Ni_prop * parameters$kappa * (1 + parameters$alpha*Ns_prop)#Ns_prop
# 
#   sqrt_term <- sqrt(1 + 4 * multiplier)
#   a_t <- (-1 + sqrt_term) / (2 * multiplier)
# 
#  #  Clamp a_t between [0,1]
#  return (a_t)
# 
#  #(max(0, min(1, a_t)))
# }
# #nd, v


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

#   HealthCostDaily = 0,
#  SocialActivityCostDaily = 0,
#      lambda_s = 0, 
#     lambda_i = 0, 
# beta = 0.3,
# H = 333.333,     # expected cost of infection
# b1 = 1,      # private benefit
# b2 = 0.5,     # complementarity
# # === Disease and belief dynamics ===
# beta = 0.3,             # biological transmission rate
# gamma = 1.5,            # curvature parameter in filter function: f(v) = (1 + v)^(-gamma)
# sigma = 0.25,           # probability a random contact is with a susceptible (used in q(pi))
# 
# # === Economic damage of infection ===
# L = 200,
# phi=1# expected cost of infection (moderate scenario, paper uses 100–500 range)

# # RUN MODEL CHECK   ####
# ####################################################
# 
# # load reference input and output
# parameters <- readRDS('test/parameters.rds')
# parameters$scenario <- "laissez-faire"
# parameters$infect_thres = 0
# initial_state <- readRDS('test/initial_state.rds')
# num_experiments <- 4
# output_reference <- readRDS('test/output_test.rds')
# 
# # run stochastic experiments
# output_test <- run_experiments(initial_state = initial_state, 
#                                times = times, 
#                                parameters = parameters, 
#                                bool_stochastic_beta = FALSE, 
#                                update_function = get_transitions_stochastic, 
#                                num_experiments)
# 
# 
# 
# 
# 
# # compare each element
# bool_output_unchanged <- TRUE
# for(i_out in 1:length(output_test)){
#   output_name <- names(output_test)[i_out]
#   if(any(dim(output_test[[i_out]]) != dim(output_reference[[i_out]])) |
#      any(output_test[[i_out]] != output_reference[[i_out]],na.rm=TRUE)){
#     # print warning
#     warning(paste0('Model output "',output_name,'" changed'))
#     bool_output_unchanged <- FALSE
#   } 
# }
# # print statement if output did not change
# if(bool_output_unchanged){
#   print('MODEL OUTPUT DID NOT CHANGE')
# }
# 
# # Update reference test files to match latest model structure
# update_reference_values <- function(){
#   saveRDS(parameters, 'test/parameters.rds')
#   saveRDS(initial_state, 'test/initial_state.rds')
#   saveRDS(output_test, 'test/output_test.rds')
# }
# 
# update_reference_values()  # Run it immediately
# 

# calculate the number of infections or recoveries
#get_transitions_stochastic <- function(n, prob){
#  return(sum(rbinom(n = n,size = 1, prob = prob)))
# cat("Calling rbinom with n =", n, ", prob =", prob, "\n")
# 
#}
# if (transitions < 0.5) {
#  return(0)
#  } else {
# transitions <- (n * prob)
# if ((n - transitions) < 1) {
#  transitions = n
#  }
#y = initial_state; times = time_pre_shock; func = sir_costate_model; parms = parameters
#colnames(states_out) <- c("Ns", "Ni", "Nr", "Nd", 
#                          "HealthCost", "SocialActivityCost", "TotalCost",
#                          "a_t", "u_t", "Rt", "HealthCostDaily", "SocialActivityCostDaily")

# HealthCostDaily <- 0
#  SocialActivityCostDaily <- 0
# lambda_s <- 0
#  lambda_i <- 0
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
#Rt parameters Ni Rt a_t_prev
#Nd, v
#  cat("DEBUG -- Ni =", Ni, "Ns =", Ns, "a_t =", a_t, "beta_t =", beta_t, "\n")
#  cat("Calculated prob =", beta_t * a_t^2 * Ni / parameters$pop_size, "\n")
#/parameters$pop_size

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
# HealthCostDaily <- scale_factor * fx_per_capita * exp(-rho_plus_delta * i_day) * Ni * parameters$kappa / parameters$pop_size
# SocialActivityCostDaily <- scale_factor * fx_per_capita * exp(-rho_plus_delta * i_day) * (Ns + Ni) * abs(u_t) / parameters$pop_size

# HealthCostDaily <- scale_factor * fx_per_capita * exp(-rho_plus_delta * i_day) * Ni * parameters$kappa / parameters$pop_size
# SocialActivityCostDaily <- scale_factor * fx_per_capita * exp(-rho_plus_delta * i_day) * (Ns + Ni) * abs(u_t) / parameters$pop_size
# 
# Still track cumulative for consistency if you like
#  HealthCost <- HealthCost + HealthCostDaily
# SocialActivityCost <- SocialActivityCost + SocialActivityCostDaily
#  Lambda_s <- Lambda_s + dLambda_s
#  Lambda_i <- Lambda_i + dLambda_i
#lambda_s <- lambda_s+d_lambda_s
#lambda_i <- lambda_i+d_lambda_i
#lambda_s, lambda_i,
# HealthCostDaily, SocialActivityCostDaily)#Lambda_s, Lamba_i
# Plot R(t) in place of lambda_s
#plot(times, R_t, type = "l", col = "blue", lwd = 2,
#    xlab = "Time (days)", ylab = "Effective Reproduction Number R(t)",
#    main = "Effective Reproduction Number Over Time")
#abline(h = 1, col = "red", lty = 2)  # Threshold of R(t) = 1
# sel_states <- #c("HealthCostDaily", "SocialActivityCostDaily",
#    names(initial_state)[!grepl('Total',names(initial_state)) & !grepl('Ns',names(initial_state))])


# sel_states <- c("Ni", "Nr", "Nd", "a_t", "u_t", "Rt", "HealthCostDaily", "SocialActivityCostDaily")
# y_lim <- range(0,output_all[,,i_state],output_sim_deterministic[,i_state],na.rm = TRUE)
# plot_label <- ifelse(i_state %in% names(pretty_labels), pretty_labels[[i_state]], i_state)
# plot(output_sim_deterministic[,i_state], col=1, main=i_state, ylab = i_state, xlab = 'time', ylim = y_lim, type='l', lwd=2) # infections, deterministic
# for(i_exp in 1:dim(output_all)[1]){
#   col_index <- which(dimnames(output_all)[[3]] == i_state)
#   lines(output_all[i_exp, , col_index], col = alpha('grey30', 0.5))
# }
# # Individual boxplot for HealthCost
# boxplot(output_summary$HealthCost,
#         main = "Health Cost",
#         ylab = "Health Cost",
#         ylim = range(c(output_summary$HealthCost, output_summary_deterministic$HealthCost)))
# points(1, mean(output_summary$HealthCost), pch = 8)
# points(1.2, output_summary_deterministic$HealthCost, col = "blue", pch = 18, cex = 1.5)
# 
# # Individual boxplot for SocialActivityCost
# boxplot(output_summary$SocialActivityCost,
#         main = "Social Activity Cost",
#         ylab = "Social Activity Cost",
#         ylim = range(c(output_summary$SocialActivityCost, output_summary_deterministic$SocialActivityCost)))
# points(1, mean(output_summary$SocialActivityCost), pch = 8)
# points(1.2, output_summary_deterministic$SocialActivityCost, col = "blue", pch = 18, cex = 1.5)
# 
# Health Cost boxplot (stochastic + deterministic overlay)
# temp_output <- run_sir_binomial(
#    initial_state = initial_state,
#    times = times,
#    parameters = parameters,
#    bool_stochastic_beta = bool_stochastic_beta,
#    update_function = update_function
#  )
# cat("Preview of temp_output columns:\n")
# print(colnames(temp_output))
# #fadeout_threshold=fadeout_threshold
# Compute and Track Effective Reproduction Number R(t)
#if (is.nan(a_t) || is.infinite(a_t) || is.na(a_t)) return(1)
# Nd_prop <- Nd / parameters$pop_size
#setwd("C:/Users/Bart Smets/OneDrive/Documenten/GitHub/Working_Code_BS/R")
