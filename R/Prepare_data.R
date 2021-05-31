#' Prepare outcome-dependently sampled time-to-event data for weight calculation.
#'
#' @description This function prepares a data set to use for the calculation of
#' weights that correct for ascertainment bias, @seealso [Calculate_weights()].
#' The specific ascertainment bias addressed is introduced by outcome-dependent
#' sampling of clustered time-to-events, for example by inclusion of high-risk
#' families in genetic epidemiological studies in cancer research.
#'
#' @details In the corresponding paper the sample interval-based survival (S_k.)
#' is distinguished from the population equivalent with an apostroph (') instead
#'  of a full stop (.). The same holds for the sample interval-based hazard
#'  (p_k.).
#'
#' @export

Prepare_data <- function(dat, population_incidence, breaks) {

  ### Input:

#' @param dat Data.frame with one row per individual which at least includes
#' a column **d** with event indicator (1 for event, 0 for censored), a column
#' **y** with event/censoring time.
#' @param population_incidence  Either a vector (in combination with breaks) or
#' a data.frame (columns start age group, end age group, incidence rate).
#' @param breaks Cut-points for the (age/time) groups.

  ### Output:

#' @return Data.frame in long format (one interval of an individual per row)
#' with a.o. the columns **S_k** population interval-based survival; **p_k**
#' population interval-based hazard; **S_k.** sample interval-based survival;
#' **p_k.** sample interval-based hazard.

  require(survival)
  require(dplyr)
  require(tidyr)

  ### From population incidence rate to p_k and S_k:

  # Organize external information (population incidence rate)
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
  } else if(breaks == "10 years"){ # Option 2: 5 years age groups.
    n_agegroups <- length(population_incidence)
    breaks <- seq(0, 10*(n_agegroups), by = 10)
    from <- breaks[-(n_agegroups + 1)]  # Remove last
    to <- breaks[-1] # Remove first
    dat_inc <- data.frame(start = from, end = to, incidence_rate = population_incidence/100000)
  }

  dat_inc <- dat_inc %>% mutate( y_cat = paste("(", start, ",", end, "]", sep=""),
                                 t = end - start)

  # Calculate S_k and p_k of population from population incidence.
  S_k <- exp(-cumsum(dat_inc$t * dat_inc$incidence_rate))
  # P(ti > end of interval k), i.e. probability of experiencing event after (survival S_k)
  p_k <- c() # Make container for p_k.
  p_k[1] <- 1-S_k[1]
  for (i in 2:n_agegroups){
    p_k[i] <- S_k[i-1]-S_k[i]
  }

  # Subset the required information.
  external <- dat_inc %>% mutate(S_k = S_k, p_k = p_k) %>%
    select(y_cat, incidence_rate, S_k, p_k)


  ### Transform data to long format:

  dat_long <- survSplit(data = dat, cut = breaks, end = "y", event = "d", episode = "interval")
#  dat_long$k <- cut(dat_long$y, breaks = breaks)
  dat_long$y_cat <- cut(dat_long$y, breaks = breaks)
  dat_long$p_k <- external$p_k[as.numeric(dat_long$y_cat)]
  dat_long$S_k <- external$S_k[as.numeric(dat_long$y_cat)]
  index <- which(is.na(dat_long$y_cat)==T)
  if(length(index)>0){
  dat_long <- dat_long[-index,]}  # Include only those that fall within an age category.

  ### Calculate p_k' and S_k' from sample:
  # Note that in the code the apostroph is replaced by ......
  # Store number of cases and unaffected individuals per age group
  agegroups_info <- aggregate(id ~ y_cat + d, data = dat_long, length)

  agegroups_info <- agegroups_info %>% spread(d, -y_cat)
  colnames(agegroups_info) <- c("k","s_k","r_k")

  # Calculate person years per outcome and age group.
  personyears <- aggregate(y ~ y_cat + d,
                           data = dat_long[,c("d","y_cat", "y"),], sum) %>%
    spread(d, -y_cat)
  colnames(personyears) <- c("k","q_k","p_k")

  # Other preparations.
  agegroups_info <- merge(agegroups_info, personyears, by="k")  # Merge.
  labs <- levels(as.factor(dat_long$k))  # Store names of age groups.
  agegroups_info$t <- dat_inc$t # Add age group length in years.

  # Calculate sample incidence (mu).
  for (k in 1:n_agegroups){

    if(k!=n_agegroups){calcsum <- sum(agegroups_info$r_k[(k+1):n_agegroups], agegroups_info$s_k[(k+1):n_agegroups])} else {calcsum <- 0}

    agegroups_info$mu_k[k] <- agegroups_info$r_k[k]/
      (agegroups_info$p_k[k] + agegroups_info$q_k[k] + agegroups_info$t[k] * calcsum)

  }

  # Calculate S_k' and p_k' from sample incidence rate.
  # In the code, ' is replaced by .
  agegroups_info<- agegroups_info %>% select(k, t, mu_k)
  agegroups_info$S_k. <- exp(-cumsum(agegroups_info$t * agegroups_info$mu_k))

  agegroups_info$p_k.[1] <- 1 - agegroups_info$S_k.[1]
  for (i in 2:n_agegroups){
    agegroups_info$p_k.[i]<-agegroups_info$S_k.[i-1] - agegroups_info$S_k.[i]
  }

  sample_info <- agegroups_info %>% select(k, p_k., S_k.)

  # Add S_k' and p_k' to the data in long format.
  df_long_out <- merge(dat_long, sample_info, by.x = "y_cat", by.y = "k")
  df_long_out$y_cat <- NULL # Remove y_cat (double).

  ### Output the data.frame:
  df_long_out

}
