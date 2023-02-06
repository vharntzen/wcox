## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----Loading, message = F-----------------------------------------------------

library(twcox)

# Also needed, amongst others for %>%.
require(dplyr)
require(survival)

## ----Load data set, message = F-----------------------------------------------
# Load toy data.
data("fam_dat")

## ----Show first lines data----------------------------------------------------
# Show the first few lines of the toy data.
head(fam_dat)

## ----Distribution of family size----------------------------------------------
# Show the number of families and top rows.
cat( "There are", unique(fam_dat$family_id) %>% length() , "families in data set. ")

# Examine family sizes.
fam_sizes <- fam_dat %>% group_by(family_id) %>% count()

cat( "Family size is on average ", 
     fam_sizes$n %>% mean() %>% round(1), 
     " and ranges from ", 
     fam_sizes$n %>% min(), " to ", 
     fam_sizes$n %>% max(), ".", sep = "")

## ----External information-----------------------------------------------------

# Enter the incidence by age group per 100 000, starting with the youngest age group.
incidence_rate <- c(2882, 1766, 1367, 1155, 987, 845, 775, 798, 636, 650)

# Define the age group cutoffs.
breaks_manually <- c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
breaks_10yrs <- "10 years"

## ----Rename variables---------------------------------------------------------
# Rename variables.
my_dat <- rename(fam_dat, d = event_indicator, y = age)

## ----Prepare data, message = F------------------------------------------------

# Using option 1 (manual cut-offs):
my_dat_out <- Prepare_data(dat = my_dat, population_incidence = incidence_rate, 
                           breaks = breaks_manually)

# Unhash to use option 2 (pre-set cut-offs):
# my_dat_out <- Prepare_data(dat = my_dat, population_incidence = incidence_rate, 
# breaks = "10 years")

## ----Look at prepared data----------------------------------------------------

# Select the newly add columns.
my_dat_out %>% select(id, tstart, d, S_k, p_k, S_k., p_k.) %>% arrange(id) %>% head(8)

## ----Calculate weights--------------------------------------------------------
# Calculate weights using the object that is output from the Prepare_data() function.

w <- Calculate_weights(dat_long = my_dat_out)

## ----Show weights-------------------------------------------------------------
# Show the weights.
my_dat_out %>% mutate(weight = w) %>% 
  select(id, d, weight) %>% 
  arrange(id) %>% 
  filter(id %in% c(1,3))

## ----Fit model using weights, message = F-------------------------------------
# Fit the model.
fit <- coxph(Surv(y, d==1) ~ x + cluster(family_id), weights = w, data =  my_dat_out)

fit

## ----Examine estimates of fitted model----------------------------------------
# Extract estimates.
fit.summary <- summary(fit)

# Summarize findings.
cat("Covariate effect (95% CI): ",
exp(fit.summary$coefficients[1]), " (",
exp(fit.summary$coefficients[1] - 1.96*fit.summary$coefficients[4]), ", ",
exp(fit.summary$coefficients[1] + 1.96*fit.summary$coefficients[4]), ").", sep = "")

