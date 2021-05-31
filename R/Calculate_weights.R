#' Calculate time-dependent weights.
#'
#' @description
#' This function calculates weights to correct for ascertainment bias in
#' time-to-event data where clusters are outcome-dependently sampled,
#' for example high-risk families in genetic epidemiological studies in
#' cancer research.
#'
#' @details
#' Weights are based on a comparison between the survival and incidence
#' rate (hazard) between sample and population. Therefore, besides the
#' sample data, the population incidence rate (per 100 000) is needed as
#' input, as well as the cut-offs of the (age/time-to-event) groups for which
#' this is available. The function provides two options for the latter:
#' cut-offs can be provided manually or using the standard 5- or 10-years (age)
#' categories (0-4, 5-9, ... or 0-9, 10-14, ...). Note that resulting intervals
#' are of the form [xx, xx).
#'
#' @export


Calculate_weights <- function(dat_long) {

### Input:

#' @param dat_long Data.frame in long format (i.e. one row per time interval)
#' with columns **k** interval of (age) group; **S_k** population interval-based
#' survival; **p_k** population interval-based hazard; **S_k.** sample interval-
#' based survival; **p_k.** sample interval-based hazard.

### Output:

#' @return Vector with weights.

# Extract variables from input data.frame:
  k <- dat_long$k
  S_k <- dat_long$S_k
  p_k <- dat_long$p_k
  S_k. <- dat_long$S_k.
  p_k. <- dat_long$p_k.

# Create empty containers:
  v <- w <- rep(NA, nrow(dat_long))

# Calculate the weights:
  for (n in 1:nrow(dat_long)){
    v[n] = ( (S_k[n]) / p_k[n] ) * (p_k.[n] / (S_k.[n]))  # Weights for unaffected
    w[n] = 1  # Weights for cases
  }

  merged <- cbind(dat_long, w, v)
  merged$weight <- NA

# Assign w for uncensored, v for censored (interval-wise):
  merged$weight[which(merged$d == 1)] <- merged$w[which(merged$d == 1)]
  merged$weight[which(merged$d == 0)] <- merged$v[which(merged$d == 0)]

  merged$w <- merged$v <- NULL  # Remove old variable

# Collect output:
  vec_weights <- merged$weight
  vec_weights[which(vec_weights==0)] <- 0.0000001  # If weight is 0, add very small value (coxph does not accept weights of 0).

# Print warning if weights are invalid:
  ifelse((sum(vec_weights<0)>0), print("Invalid (negative) weights!"), print("No negative weights"))  # If there are weights of zero, skip this iteration

# Return output:
  vec_weights

}
