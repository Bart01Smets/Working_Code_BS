
#   # Expression 1: Altruism case (Equation 30)
#   if (version == "altruism") {
#     denom <- (1 + alpha) * mult
#     a_t <- sqrt(num / denom)
#   }
# 
#   # Expression 2: Optimal policy (Equation 21)
#   else if (version == "optimal") {
#     denom <- 2 * mult
#     a_t <- sqrt(num / denom)
#   }
# 
#   # Expression 3: Quadratic root (Equation 11)
#   else if (version == "quadratic") {
#     sqrt_term <- 1 + 8 * mult * num
#     if (sqrt_term < 0 || mult == 0) return(1)
#     a_t <- (-1 + sqrt(sqrt_term)) / (4 * mult)
#   }
#   else if (version == "laissez-faire") {
#     # Derived: (-num + sqrt(num^2 + 8 * mult * num)) / (4 * mult)
#     sqrt_term <- num^2 + 8 * mult * num
#     if (sqrt_term < 0 || mult == 0) return(1)
#     a_t <- (-num + sqrt(sqrt_term)) / (4 * mult)
#   }
#   else if (version == "altruistic_myopic") {
#     denom <- 2 * (1 + alpha) * mult
#     if (denom == 0) return(1)
#     a_t <- sqrt(num / denom)
#   }#ok results
#   else if (version == "altr_myo") {
#     # Myopic altruism where external harm is proxied by kappa
#     harm_term <- kappa + alpha * kappa  # i.e., (1 + alpha) * kappa
#     denom <- 2 * beta * Ns_prop * Ni_prop * harm_term
#     if (denom <= 0) return(1)
#     a_t <- sqrt(num / denom)}
#     
#     else if (version == "myopic_altruistic") {
#       # Full expression: a* = (-1 + sqrt(1 + 8βκnsni(1 + αβns))) / [4βκnsni(1 + αβns)]
#       denom <- 4 * beta * kappa * Ns_prop * Ni_prop * (1 + alpha * beta * Ns_prop)
#       sqrt_term <- 1 + 8 * beta * kappa * Ns_prop * Ni_prop * (1 + alpha * beta * Ns_prop)
#       if (denom <= 0 || sqrt_term < 0) return(1)
#       a_t <- (-1 + sqrt(sqrt_term)) / denom
#   }
# 

}
#   
#   # Fallback if invalid version is passed
#   else {
#     warning("Invalid 'version' passed to a_function(). Returning 1.")
#     return(1)
#   }
#   
#   return(max(0, min(1, a_t)))
# }
# Expression 4: Laissez-faire (Myopic, α = 0)




# 
# a_function <- function(Ni, Ns, parameters) {
#   if (Ni < 1e-10 || Ns < 1e-10) return(1)
# 
#   # Compute proportions
#   Ni_prop <- Ni / parameters$pop_size
#   Ns_prop <- Ns / parameters$pop_size
#   denom_pop <- Ns_prop + Ni_prop
# 
#   # Guard against division by zero
#   if (denom_pop < 1e-10) return(1)
# 
#   # Define M = β * κ * (n_s * n_i) / (n_s + n_i)
#   num <- Ns_prop * Ni_prop
#   M <- parameters$beta * parameters$kappa * num / denom_pop
# 
#   sqrt_term <- 1 + 8 * M
#   if (M == 0 || sqrt_term < 0) return(1)
# 
#   # Optimal activity
#   a_t <- (-1 + sqrt(sqrt_term)) / (4 * M)
# 
#   return(max(0, min(1, a_t)))
# }