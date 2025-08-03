rm(list = ls())
library(simcol)
# library(assertthat)

# Set Seed for Reproducibility
sim_seed <- 123
set.seed(sim_seed)

# Model 1: Baseline Model, sigma < theta, therefore no leniency applications
sim_list <- sim_col()
cartels_duration <- get_cartel_duration(sim_list$cartels, sim_list$detection, sim_list$leniency)
# assert_that(sum(cartels_duration$leniency) == 0)
# assert_that(nrow(cartels_duration) == 254)

################################################################################################

# Model 2A: sigma > theta, all firms apply for leniency when cartel collapses, no inactive undetected cartels
sim_list <- sim_col(model = "2", periods = 100, n_ind = 100, theta = 0.035)
cartels_duration <- get_cartel_duration(sim_list$cartels, sim_list$detection, sim_list$leniency)
# assert_that(nrow(filter(cartels_duration, in_sample==0 & end < 10)) == 0)
# assert_that(nrow(cartels_duration) == 1)

################################################################################################

# Model 2B: cartels apply for leniency with probability prob_leniency
sim_list <- sim_col(model = "2B", periods = 100, n_ind = 100, theta = 0.035)
cartels_duration <- get_cartel_duration(sim_list$cartels, sim_list$detection, sim_list$leniency)
input_ind <- sim_list$input_ind
cartels_duration <- left_join(cartels_duration, input_ind, by="industry")
# assert_that(nrow(cartels_duration) == 314)
################################################################################################

# Model 3: Probability of CA succes (s) and therefore sigma depend on CA caseload
sim_list <- sim_col(model = "3", periods = 10, n_ind = 10)
cartels_duration <- get_cartel_duration(sim_list$cartels, sim_list$detection, sim_list$leniency)
input_ind <- sim_list$input_ind
cartels_duration <- left_join(cartels_duration, input_ind, by="industry")
# assert_that(nrow(cartels_duration) == 1)
