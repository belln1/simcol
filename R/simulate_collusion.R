is_IC <- function(pi, icc){
  pi < icc
}

get_bernoulli <- function(n, prob) {
  as.numeric(runif(n) <= prob)
}

enter_cartel <- function(n, k, pi, icc){
  get_bernoulli(n, k) * is_IC(pi, icc)
}

stay_cartel <- function(pi, icc){
  is_IC(pi, icc)
}

get_rate <- function(sums, n_ind) {
  sums/n_ind
}

#' Simulate Collusion
#'
#' Simulation of industries with collusive and non-collusive periods over time, detected and undetected, based on Harrington/Chang (2015)
#' @param model Modelname (1, 2, 2B, 3)
#' @param periods Number of time periods
#' @param n_ind Number of industries
#'
#' @return List with panel data cartels, detected cartels, leniency cases, profit over time, time varying parameters, ICC over time
#' @import dplyr
#' @import EnvStats
#' @examples
#' sim_list <- sim_col();
#' sim_list <- sim_col(model = "2", periods = 1000, n_ind = 1000, theta = 0.035);
#' sim_list <- sim_col(model = "2B", periods = 100, n_ind = 100, theta = 0.035);
#' sim_list <- sim_col(model = "3", periods = 10, n_ind = 10);
#' @export
sim_col <- function(model="1", periods=100, n_ind=100,
                     delta=0.85,
                     pi_mean = 0, pi_sd = 1.5, pi_min = 1, pi_max = 100,
                     k_min = 0, k_max = 0.1,
                     eta_mean = 1, eta_sd = 1.5, eta_min = 1.1, eta_max = 10,
                     S = 0.8, Q = 0.2, R = 0.3, gamma = 2,
                     theta = 1, prob_leniency = 1,
                     lambda = 0.5, tau=0.8, xi=1, v=500, rho=1.5
                     ){
  MU <- get_mu(meanlog = pi_mean, sdlog = pi_sd, min = pi_min, max = Inf)
  K <- get_kappa(n_ind, min=k_min, max=k_max)
  eta <- get_eta(n_ind, meanlog=eta_mean, sdlog=eta_sd, min=eta_min, max=eta_max)
  sigma <- get_sigma(Q, R, S)

  industry <- 1:n_ind
  input_ind <- data.frame(cbind(industry, delta, pi_mean, pi_sd, MU, pi_min, pi_max, K, eta, prob_leniency))
  input_enforce <- data.frame(cbind(Q, R, S, sigma, gamma, theta, lambda, tau, xi, v, rho))

  if (model == 3) {
    s_t <- rep(0, periods)
    s_t[1] <- get_s(lambda, 0, 0, tau, xi, v, rho)
    s_t[1]
    sigma_t <- rep(0, periods)
    time_enforce <- data.frame(cbind(s_t, sigma_t))
    sim_list <- simulate_collusion_M3(periods, n_ind, input_ind, input_enforce, time_enforce)

  } else if (model == "2B") {
      input_ind$prob_leniency <- get_prob_leniency(input_ind$K)
      sim_list <- simulate_collusion_M1_M2(periods, n_ind, input_ind, input_enforce)

  } else {
    sim_list <- simulate_collusion_M1_M2(periods, n_ind, input_ind, input_enforce)
  }
  return(sim_list)
}

simulate_collusion_M3 <- function(periods, n_ind, input_ind, input_enforce, time_enforce) {
  # Initialize Matrices
  in_cartel <- matrix(0, nrow = periods, ncol = n_ind)
  in_q_r <- matrix(0, nrow = periods, ncol = n_ind)
  detection <- matrix(0, nrow = periods, ncol = n_ind)
  leniency <- matrix(0, nrow = periods, ncol = n_ind)
  pi <- matrix(0, nrow = periods, ncol = n_ind)
  icc <- matrix(0, nrow = periods, ncol = n_ind)

  # First Period
  time_enforce$sigma_t[1] <- get_sigma(input_enforce$Q, input_enforce$R, time_enforce$s_t[1])

  # ICC
  if (time_enforce$sigma_t[1] <= input_enforce$theta){
    icc[1,] <- get_ICC_no_leniency(input_ind, time_enforce[1,])
  } else {
    icc[1,] <- get_ICC_leniency_time(input_ind, time_enforce[1,], input_enforce)
  }

  pi[1,] <- get_pi_t(n_ind, mean=input_ind$pi_mean, sd=input_ind$pi_sd, min=input_ind$pi_min, max=input_ind$pi_max)
  in_cartel[1,] <- enter_cartel(n_ind, input_ind$K, pi[1,], icc[1])
  in_q_r[1,] <- in_cartel[1,] * get_bernoulli(n_ind, input_enforce$Q*input_enforce$R)
  detection[1,] <- in_q_r[1,] * get_bernoulli(n_ind, time_enforce$s_t[1])
#  assert_that(min(as.numeric(in_cartel - detection >= 0)) > 0, msg = "cartel >= detection")
  sum_L <- 0
  sum_R <- get_rate(sum(in_q_r[1,]), n_ind)

  # Set time-variant values Second Period
  time_enforce$s_t[2] <- get_s(input_enforce$lambda, sum_L, sum_R, input_enforce$tau, input_enforce$xi, input_enforce$v, input_enforce$rho)
  time_enforce$sigma_t[2] <- get_sigma(input_enforce$Q, input_enforce$R, time_enforce$s_t[2])
  if (time_enforce$sigma_t[2] <= input_enforce$theta){
    icc[2,] <- get_ICC_no_leniency(input_ind, time_enforce[2,])
  } else {
    icc[2,] <- get_ICC_leniency_time(input_ind, time_enforce[2,], input_enforce)
  }

  # All Other Periods
  for (i in 2:periods) {
    pi[i,] <- get_pi_t(n_ind, mean=input_ind$pi_mean, sd=input_ind$pi_sd, min=input_ind$pi_min, max=input_ind$pi_max)
    previous_cartel <- in_cartel[i-1,]
    previous_in_q_r <- in_q_r[i-1,]
    previous_detection <- detection[i-1,]
    previous_leniency <- leniency[i-1,]
    # entry: no cartel previous and ICC fit
    # stay: cartel previous and ICC fit
    in_cartel[i,] <- (1 - previous_detection) * (1 - previous_leniency) * (
      previous_cartel * stay_cartel(pi[i,], icc[i,]) +
        (1 - previous_cartel) * enter_cartel(n_ind, input_ind$K, pi[i,], icc[i,])
    )
#    assert_that(min(in_cartel[i,])>=0, msg = "cartel >=0")
    # detection today: in cartel and sigma_t
    in_q_r[i,] <- in_cartel[i,] * get_bernoulli(n_ind, input_enforce$Q*input_enforce$R)
    detection[i,] <- in_q_r[i,] * get_bernoulli(n_ind, time_enforce$s_t[i])
#    assert_that(min(detection[i,])>=0, msg="detection>=0")
    # leniency yesterday: undetected cartel yesterday and ICC do not fit today and theta < sigma_t
    leniency[i-1,] <- (input_enforce$theta < time_enforce$sigma_t[i]) * (1 - in_cartel[i,]) * (previous_cartel) * (1 - previous_detection)
#    assert_that(min(leniency[i-1,])>=0, msg="leniency>=0")
    # detection yesterday: might be already detected. or not detected, but ICC do not fit today and no leniency and sigma_t
    # previous_in_q_r: filled in last period. in_q_r[i-1,]: filled in this period.
    in_q_r[i-1,] <- previous_in_q_r +
      (1 - in_cartel[i,]) * previous_cartel * (1 - previous_in_q_r) * (1 - leniency[i-1,]) * get_bernoulli(n_ind, input_enforce$Q*input_enforce$R)
    detection[i-1,] <- previous_detection +
      (1 - previous_in_q_r) * in_q_r[i-1,] * get_bernoulli(n_ind, time_enforce$s_t[i])

#    assert_that(min(detection[i-1,])>=0, msg="detection(i-1)>=0")
    sum_L <- get_rate(sum(leniency[i-1,]), n_ind)
    sum_R <- get_rate(sum(in_q_r[i-1,]), n_ind)
    # assert_that(sum_L >= 0, msg="L>=0")
    # assert_that(sum_R >= 0, msg="R>=0")

    if (i < periods) {
      time_enforce$s_t[i+1] <- get_s(input_enforce$lambda, sum_L, sum_R, input_enforce$tau, input_enforce$xi, input_enforce$v, input_enforce$rho)
      time_enforce$sigma_t[i+1] <- get_sigma(input_enforce$Q, input_enforce$R, time_enforce$s_t[i+1])

      if (time_enforce$sigma_t[i+1] <= input_enforce$theta){
        icc[i+1,] <- get_ICC_no_leniency(input_ind, time_enforce[i+1,])
      } else {
        icc[i+1,] <- get_ICC_leniency_time(input_ind, time_enforce[i+1,], input_enforce)
      }
    }
  }
  return(list(cartels = in_cartel, detection = detection, leniency = leniency, profit = pi, time_enforce = time_enforce, icc = icc, in_q_r=in_q_r, input_ind=input_ind))
}

simulate_collusion_M1_M2 <- function(periods, n_ind, input_ind, input_enforce) {
  # Initialize Matrices
  in_cartel <- matrix(0, nrow = periods, ncol = n_ind)
  in_q_r <- matrix(0, nrow = periods, ncol = n_ind)
  detection <- matrix(0, nrow = periods, ncol = n_ind)
  leniency <- matrix(0, nrow = periods, ncol = n_ind)
  pi <- matrix(0, nrow = periods, ncol = n_ind)

  # ICC
  if (input_enforce$sigma <= input_enforce$theta){
    ICC <- get_ICC_no_leniency(input_ind, input_enforce)
  } else {
    ICC <- get_ICC_leniency(input_ind, input_enforce)
  }

  # First Period
  pi[1,] <- get_pi_t(n_ind, mean=input_ind$pi_mean, sd=input_ind$pi_sd, min=input_ind$pi_min, max=input_ind$pi_max)
  in_cartel[1,] <- enter_cartel(n_ind, input_ind$K, pi[1,], ICC)
  in_q_r[1,] <- in_cartel[1,] * get_bernoulli(n_ind, input_enforce$Q*input_enforce$R)
  detection[1,] <- in_cartel[1,] * get_bernoulli(n_ind, input_enforce$S)
  # assert_that(min(as.numeric(in_cartel - detection >= 0)) > 0)

  # All Other Periods
  for (i in 2:periods) {
    pi[i,] <- get_pi_t(n_ind, mean=input_ind$pi_mean, sd=input_ind$pi_sd, min=input_ind$pi_min, max=input_ind$pi_max)
    previous_cartel <- in_cartel[i-1,]
    previous_in_q_r <- in_q_r[i-1,]
    previous_detection <- detection[i-1,]
    previous_leniency <- leniency[i-1,]
    # entry: not cartel previous and ICC fit
    # stay: cartel previous and ICC fit
    in_cartel[i,] <- (1 - previous_detection) * (1 - previous_leniency) * (
      previous_cartel * stay_cartel(pi[i,], ICC) +
        (1 - previous_cartel) * enter_cartel(n_ind, input_ind$K, pi[i,], ICC)
    )
    # detection today: in cartel and sigma
    in_q_r[i,] <- in_cartel[i,] * get_bernoulli(n_ind, input_enforce$Q*input_enforce$R)
    detection[i,] <- in_q_r[i,] * get_bernoulli(n_ind, input_enforce$S)
    # leniency yesterday: undetected cartel yesterday and ICC do not fit today and theta < sigma
    leniency[i-1,] <- (input_enforce$theta < input_enforce$sigma) * get_bernoulli(n_ind, input_ind$prob_leniency) * (1 - in_cartel[i,]) * (previous_cartel) * (1 - previous_detection)
    # detection yesterday: might be already detected. or not detected, but ICC do not fit today and no leniency and sigma
    in_q_r[i-1,] <- previous_in_q_r +
      (1 - in_cartel[i,]) * previous_cartel * (1 - previous_in_q_r) * (1 - leniency[i-1,]) * get_bernoulli(n_ind, input_enforce$Q*input_enforce$R)
    detection[i-1,] <- previous_detection +
      (1 - previous_in_q_r) * in_q_r[i-1,] * get_bernoulli(n_ind, input_enforce$S)
  }
  return(list(cartels = in_cartel, detection = detection, leniency = leniency, profit = pi, input_enforce = input_enforce, icc = ICC, in_q_r=in_q_r, input_ind=input_ind))
}



