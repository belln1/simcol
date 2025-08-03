# Functions to calculate cartel durations and data analytics examples  --------------------------------------------------------------

get_starttimes <- function(cartels) {
  d <- dplyr::lag(cartels)
  d[is.na(d)] = 0
  e <- as.numeric((cartels - d) == 1)
  starttimes <- matrix(e, ncol = ncol(cartels))
}

get_endtimes <- function(cartels) {
  d <- dplyr::lead(cartels)
  d[is.na(d)] = 0
  e <- as.numeric((cartels - d) == 1)
  endtimes <- matrix(e, ncol = ncol(cartels))
}

#' Cartel Duration
#'
#' Function to combine paneldata cartels, detected cartels, leniency cases returned by function sim_col into a cross-sectional dataset with duration
#' @param cartel_dates binary collusive state
#' @param detect_dates binary collusive state that eventually got detected
#' @param leniency_dates binary collusive state that eventually got found due to leniency application
#'
#' @return cross-sectional dataset with duration
#' @examples
#' sim_list <- sim_col();
#' df <- get_cartel_duration(sim_list$cartels, sim_list$detection, sim_list$leniency);
#' @export
get_cartel_duration <- function(cartel_dates, detect_dates, leniency_dates) {
  starttimes <- get_starttimes(cartel_dates)
  endtimes <- get_endtimes(cartel_dates)
  start <- as_tibble(which(starttimes==1, arr.ind = TRUE))
  start <- dplyr::rename(start, "start" = row)
  end <- as_tibble(which(endtimes==1, arr.ind = TRUE))
  end <- dplyr::rename(end, "end" = row)
  detect <- as_tibble(which(detect_dates==1, arr.ind = TRUE))
  detect <- dplyr::rename(detect, "detect_date" = row)
  detect$detected = 1
  lenienc <- as_tibble(which(leniency_dates==1, arr.ind = TRUE))
  lenienc <- dplyr::rename(lenienc, "leniency_date" = row)
  lenienc$leniency = 1
  df <- cbind(start, end=end$end)
  df <- df %>%
    dplyr::rename(industry = col) %>%
    dplyr::mutate(duration = end - start + 1,
                  cartel = 1) %>%
    left_join(detect, by = join_by(industry==col, end==detect_date)) %>%
    mutate(detected = ifelse(is.na(detected), 0, detected)) %>%
    left_join(lenienc, by = join_by(industry==col, end==leniency_date)) %>%
    mutate(leniency = ifelse(is.na(leniency), 0, leniency),
           in_sample = ifelse((detected==1) | (leniency==1), 1, 0)) %>%
    group_by(industry) %>%
    arrange(industry, start) %>%  # order rows
    dplyr::mutate(cartel = cumsum(cartel),
                  nTc = cumsum(in_sample),
                  rep_off = ifelse(nTc > 1, 1, 0)) %>%
    ungroup %>%
    relocate(industry, start, end, duration)  # order columns
  return(df)
}

get_sample_panel <- function(sample_duration, PERIODS, N_INDUSTRIES) {
  grid_sample <- expand.grid(
    period = 1:PERIODS,
    industry = 1:N_INDUSTRIES
  )
  grid_sample$in_cartel <- 0
  for (i in 1:nrow(sample_duration)) {
    grid_sample$in_cartel <- pmax(
      grid_sample$in_cartel,
      as.integer(
        grid_sample$industry == sample_duration$industry[i] &
          grid_sample$period >= sample_duration$start[i] &
          grid_sample$period <= sample_duration$end[i]
      )
    )
  }
  output <- grid_sample %>%
    spread(industry, in_cartel) %>%
    select(-period)
  return(output)
}

# Convert to Paneldata
get_endtimes_vector <- function(in_cartel) {
  d <- lead(in_cartel)
  d[is.na(d)] <- 0
  e <- as.numeric((in_cartel - d) == 1)
}

get_rep_off <- function(in_cartel) {
  e <- get_endtimes_vector(in_cartel)
  ec <- cumsum(e)
  lec <- lag(ec)
  lec[is.na(lec)] <- 0
  rep <- ifelse(lec>0, 1, 0)
}

get_paneldata <- function(cartels_population, cartels_duration, periods, num_industries) {
  sample_duration <- filter(cartels_duration, detected==1)
  cartels_sample <- get_sample_panel(sample_duration, periods, num_industries)
  leniency_duration <- filter(cartels_duration, leniency==1)
  cartels_sample$period <- 1:periods
  cartels_sample1 <- cartels_sample %>%
    gather(key = "industry", value = "detected", 1:num_industries)
  mean(cartels_sample1$detected)

  colnames(cartels_population) <- 1:num_industries
  cartels_population$period <- 1:periods
  cartels_population <- cartels_population %>%
    gather(key = "industry", value = "in_cartel", 1:num_industries)
  mean(cartels_population$in_cartel)

  cartels_sample2 <- cartels_sample1 %>%
    mutate(in_cartel = cartels_population$in_cartel,
           industry = as.numeric(industry))

  if (nrow(leniency_duration)>0){
    cartels_leniency <- get_sample_panel(leniency_duration, periods, num_industries)
    cartels_leniency$period <- 1:periods
    cartels_leniency <- cartels_leniency %>%
      gather(key = "industry", value = "leniency", 1:num_industries)
    mean(cartels_leniency$leniency)
    cartels_sample2 <- cartels_sample2 %>%
      mutate(leniency = cartels_leniency$leniency)
  } else {
    cartels_sample2 <- cartels_sample2 %>%
      mutate(leniency = 0)
  }

  df <- cartels_sample2 %>%
    mutate(in_sample = detected + leniency) %>%
    group_by(industry) %>%
    mutate(R = get_rep_off(in_sample)) %>%
    ungroup() %>%
    relocate(industry, period, in_cartel, detected, leniency, in_sample, R)

  return(df)
}

