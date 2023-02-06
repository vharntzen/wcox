#' Prepare outcome-dependently sampled time-to-event data for weight calculation.
#'
#' @description This function prepares a data set to use for the calculation of
#' weights that correct for ascertainment bias, @seealso [Calculate_weights()].
#' The specific ascertainment bias addressed is introduced by outcome-dependent
#' sampling of clustered time-to-events, for example by inclusion of high-risk
#' families in genetic epidemiological studies in cancer research.
#'
#' @details In the corresponding paper the sample interval-based survival (S_k.)
#' is distinguished from the population equivalent with a hyperscript 0 instead
#'  of a full stop (.). The same holds for the sample interval-based hazard
#'  (p_k.).
#'
#' @export

Prepare_data <- function(dat, population_incidence, breaks){

  ### Input:

#' @param dat Data.frame with one row per individual which at least includes
#' a column **d** with event indicator (1 for event, 0 for censored), a column
#' **y** with event/censoring time.
#' @param population_incidence A vector (in combination with breaks) or
#' a data.frame (columns 1) 'start age group', 2) 'end age group', 3)'S_pop') with
#' population incidence per 100,000 per interval k.
#' @param breaks Cut-points for the (age/time) groups. Only needed when population_incidence is a vector.

  ### Output:

#' @return Data.frame one row per individual and a.o. columns *d*
  #' non-censoring indicator; **k** interval of (age) group; **S_k**
  #' population interval-based proportion of individuals experiencing the
  #' event in intervals later than k; **p_k** population proportion of
  #' individuals experiencing the event within interval k; **S_k.** sample
  #' proportion of individuals experiencing the event in intervals later
  #' than k; **p_k.** sample proportion of individuals experiencing the event
  #' within interval k.

  require(survival)
  require(dplyr)
  require(tidyr)

  ### From population incidence rate to p_k and S_k:

  # Organize external information (population incidence rate)
  if(is.vector(population_incidence)){
  if(is.numeric(breaks)){ # Option 1: two vectors population_incidence and breaks.
    if(length(breaks)!=(length(population_incidence)+1) ) warning("Number of breaks should equal the number of groups plus one.")
    n_agegroups <- length(population_incidence)
    from <- breaks[-(n_agegroups + 1)]  # Remove last
    to <- breaks[-1] # Remove first
    dat_inc <- data.frame(start = from, end = to, incidence_rate = population_incidence/100000)
  } else if(breaks == "5 years"){ # Option 2: 5 years age groups.
    n_agegroups <- length(population_incidence)
    breaks <- seq(0, 5*(n_agegroups), by = 5)
    from <- breaks[-(n_agegroups + 1)]  # Remove last
    to <- breaks[-1] # Remove first
    dat_inc <- data.frame(start = from, end = to, incidence_rate = population_incidence/100000)
  } else if(breaks == "10 years"){ # Option 3: 10 years age groups.
    n_agegroups <- length(population_incidence)
    breaks <- seq(0, 10*(n_agegroups), by = 10)
    from <- breaks[-(n_agegroups + 1)]  # Remove last
    to <- breaks[-1] # Remove first
    dat_inc <- data.frame(start = from, end = to, incidence_rate = population_incidence/100000)
  }
  } else { # If population incidence is given as a data.frame:

    dat_inc <- data.frame(start = population_incidence[,1],
                          end = population_incidence[,2],
                          incidence_rate = population_incidence[,3]/100000)
    breaks <- population_incidence[,1]
  }

  dat_inc <- dat_inc %>% mutate( k = paste("(", start, ",", end, "]", sep=""),
                                 t = end - start)

  # Calculate S_k and p_k of population from population incidence.
  Surv_k <- exp(-cumsum(dat_inc$t * dat_inc$incidence_rate)) # Survival up to end of interval.
  # P(ti > end of interval k), i.e. probability of experiencing event after (survival S_k)
  p_k <- c() # Make container for p_k.
  p_k[1] <- 1 - Surv_k[1]
  for (i in 2:n_agegroups){
    p_k[i] <- Surv_k[i]/Surv_k[i-1]
  }
  S_k <- 1 - p_k

  # Subset the required information.
  external <- dat_inc %>% mutate(S_k = S_k, p_k = p_k) %>%
    select(k, incidence_rate, S_k, p_k)


  ### Transform data to long format:
  dat_long <- survSplit(data = dat, cut = breaks, end = "y", event = "d", episode = "interval")

  dat_long$k <- cut(dat_long$y, breaks = breaks)
   index <- which(is.na(dat_long$k)==T)
  if(length(index)>0){
  dat_long <- dat_long[-index,]}  # Include only those that fall within an age category.

  ### Calculate p_k' and S_k' from sample:
  # Note that in the code the hyperscript 0 is replaced by '.'
  # Store number of cases and unaffected individuals per age group
  agegroups_info <- aggregate(id ~ k + d, data = dat_long, length)

  agegroups_info <- agegroups_info %>% spread(d, -k)
  colnames(agegroups_info) <- c("k","censored","event")

  # Define p_k. and S_k.:
  agegroups_info$p_k. <- agegroups_info$event/(agegroups_info$event + agegroups_info$censored)
  agegroups_info$S_k. <- 1 - agegroups_info$p_k.

  # Add S_k' and p_k' to the data.
  dat$k <- cut(dat$y, breaks = breaks)
  dat$tstart <- breaks[as.numeric(dat$k)]
  dat$p_k <- external$p_k[as.numeric(dat$k)]
  dat$S_k <- external$S_k[as.numeric(dat$k)]
  dat$p_k. <- agegroups_info$p_k[as.numeric(dat$k)]
  dat$S_k. <- agegroups_info$S_k[as.numeric(dat$k)]

  ### Output the data.frame:
  dat

}
