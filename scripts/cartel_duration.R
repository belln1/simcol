# Functions to calculate cartel durations  --------------------------------------------------------------

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

get_cartel_duration <- function(cartels, detect_dates, leniency_dates) {
  starttimes <- get_starttimes(cartels)
  endtimes <- get_endtimes(cartels)
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

