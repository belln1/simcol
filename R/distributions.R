# Pre-defined Value Distributions

# Mean of Profit in Log-Normal Distribution
get_mu <- function(meanlog, sdlog, min, max) {
  integrand <- function(x) {
    # PDF evaluated at (x - pftmin) = (x - 1)
    x * dlnorm(x - min, meanlog = meanlog, sdlog = sdlog)
  }
  result <- integrate(integrand, lower = min, upper = max)
  return(result$value)
}

# Stochastic Profit in every Time Period
get_pi_t <- function(n, meanlog, sdlog, min, max){
  x <- EnvStats::rlnormTrunc(n, meanlog=meanlog, sdlog=sdlog, min=min, max=max)
  return(x + min)
}

# Industry-specific Probability of Opportunity to Cartelize
get_kappa <- function(n, min, max) {
  runif(n, min=min, max=max)
}  

# Industry-specific Proneness to Collude
get_eta <- function(n, meanlog, sdlog, min, max) {
  EnvStats::rlnormTrunc(n, meanlog=meanlog, sdlog=sdlog, min=min, max=max)
}

# ICC with Leniency if theta < sigma, Model 2 (constant sigma)
get_ICC_leniency <- function(input_ind, input_enf){
  icc <- (input_ind$delta * (1-input_enf$sigma) * (1-input_ind$K) * input_ind$MU ) / 
    ((input_ind$eta - 1) * (1 - (input_ind$delta * (1-input_ind$K))) ) - 
    ((input_enf$sigma - input_enf$theta) * input_enf$gamma * input_ind$MU) /
    (((input_ind$eta) - 1) )
  return(icc)
}

# ICC with Leniency if theta < sigma, Model 3 (time-variant sigma)
get_ICC_leniency_time <- function(input_ind, time_enf, input_enf){
  icc <- (input_ind$delta * (1-time_enf$sigma) * (1-input_ind$K) * input_ind$MU ) / 
    ((input_ind$eta - 1) * (1 - (input_ind$delta * (1-input_ind$K))) ) -
    ((time_enf$sigma - input_enf$theta) * input_enf$gamma * input_ind$MU) /
    (((input_ind$eta) - 1) )
  return(icc)
}

# ICC without Leniency if sigma < theta, all Models
get_ICC_no_leniency <- function(input_ind, time_enf){
  icc <- (input_ind$delta * (1-time_enf$sigma) * (1-input_ind$K) * input_ind$MU ) / 
    ((input_ind$eta - 1) * (1 - (input_ind$delta * (1-input_ind$K))) )
  return(icc)
}

# Functions Model 2B
# Probability of Leniency Application (given theta < sigma)
get_prob_leniency <- function(K){
  (1 - (10 * K))/2
}

# Functions Model 3
# Probability sigma of Detection and Penalty
get_sigma <- function(q, r, s){
  q * r * s
}
# Probability s of CA Success
get_s <- function(lambda, sum_L, sum_R, tau=0.8, xi=1, v=500, rho=1.5) {
  tau / (xi + (v * (((lambda * sum_L) + sum_R)^rho)))
}
#assert_that(get_s(lambda = 0.5, sum_L= 0, sum_R= 0, tau = 0.8, xi =1, v = 500, rho = 1.5)==0.8)
