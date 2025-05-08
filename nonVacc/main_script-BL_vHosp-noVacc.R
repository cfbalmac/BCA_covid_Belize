# Benefit-Cost Analysis to assess a vaccination strategy for Belize
# Main script
# Author: Carlos Balmaceda

library(readr)
library(dplyr)
library(tidyverse)
library(kableExtra)
library(knitr)
library(ggplot2)
library(foreach)
library(doMC)
library(lubridate)
library(magrittr)
library(coda)
library(tidyverse)
library(rootSolve)
library(mgcv)
library(doParallel)
library(imputeTS)
library(rjags)
library(bayescount)
library(signal)

#rm(list=ls(all=TRUE))


# Load datasets, functions and parameters ----------------------------------------------
# Read the database
#owid_covid_data <- read_csv("data/owid-covid-data_v2-onl.csv") # Data from owid for Belize
covid_data_MoH <- read_csv("data/covid-data-MoH.csv")   # Data from the MoH of Belize - from 22/03/2020 to 30/06/2023

start_date <- min(covid_data_MoH$date) # first case
end_date <- max(covid_data_MoH$date) 
date_range <- seq(start_date,end_date,1)

covid_data_MoH1 <- tibble(date = date_range) %>%
  left_join(covid_data_MoH, by = "date") %>%
  mutate(new_cases_miss = ifelse(is.na(new_cases), NA, new_cases))

covid_data_MoH2 <- covid_data_MoH1 %>% 
  dplyr::filter(date <= as.Date("2021-03-01"))

# Imputation using Kalman Filter algorithm
covid_data_MoH2$new_cases <- na_kalman(covid_data_MoH2$new_cases_miss)
covid_data_MoH2$new_cases <- round(covid_data_MoH2$new_cases)
covid_data_MoH2 <- covid_data_MoH2 %>%
  select(date, new_cases)

# --------------------------------------------------------------
daily_cases <- covid_data_MoH2 %>%
  mutate(total_cases= cumsum(new_cases))

# --------------------------------------------------------------
# Creating smooth data using a Bayesian poisson-gamma model

# --------------------------------------------------------------
# Define the Bayesian Poisson Smoothing function
#bayesian_poisson_smoothing <- function(data, alpha = 1, beta = 1) {
#  smoothed_data <- numeric(length(data))

#  for (i in 1:length(data)) {
# Apply the Gamma-Poisson smoothing
#    smoothed_data[i] <- (alpha + data[i]) / (beta + 1)
#  }

# Round the smoothed values to integers since it's count data
#  return(round(smoothed_data))
#}
#smooth_cases <- bayesian_poisson_smoothing(daily_cases$new_cases, alpha = 10, beta = 5)

#plot(date_range, smooth_cases)

# --------------------------------------------------------------
# Creating smooth data using a Spline interpolation (with rounding)
#spline_fit <- spline(x = 1:length(daily_cases$date), y = daily_cases$new_cases, n = nrow(daily_cases))
#smooth_cases.sp <- round(spline_fit$y)
#plot(date_range, smooth_cases.sp)

# --------------------------------------------------------------
# Creating smooth data using median filtering
#smooth_cases.mf <- runmed(daily_cases$new_cases, 7)  # 7 is the window size
#plot(date_range, smooth_cases.mf)

# --------------------------------------------------------------
#daily_cases <- daily_cases %>%
#  mutate(new_cases=smooth_cases.mf)

# write.csv(daily_cases,"daily_cases.csv")

# --------------------------------------------------------------
# Only for owid dataset 
# --------------------------------------------------------------

# Filter the dababase to Belize only (owid dataset)
#bl_df <- owid_covid_data %>%
#  filter(location=="Belize")

#summary(bl_df)  # 1 NA (obs 827)
#bl_df$new_cases <- ifelse(is.na(bl_df$new_cases),0,bl_df$new_cases)
#bl_df$new_cases_smoothed <- ifelse(is.na(bl_df$new_cases_smoothed),0,bl_df$new_cases_smoothed)

#df_cases_death <- bl_df %>%
#  select(date,total_cases,new_cases,new_cases_smoothed,total_deaths,new_deaths,new_deaths_smoothed) %>%
#  filter(date >= as.Date("2020-03-29"))

# dataset with daily cases
#daily_cases <- owid_covid_data %>%
#    filter(date >= as.Date("2020-08-03"))   # 139.04718 total cases
#daily_cases$new_cases <- round(daily_cases$new_cases)
#daily_cases$total_cases <- round(daily_cases$total_cases)

# --------------------------------------------------------------
# --------------------------------------------------------------
# Calculate daily log returns of new cases
daily_cases2 <- daily_cases %>%
  mutate(log_return = log(new_cases / lag(new_cases)))

# Calculate volatility (standard deviation of log returns)
volatility <- sd(daily_cases2$log_return, na.rm = TRUE)
#volatility <- 0.3

# --------------------------------------------------------------
# Load Hospitalisation data
# --------------------------------------------------------------
hosp_data_IHME <- read_csv("data/hosp-ihme-bl.csv")   # Data from IHME - from 22/03/2020 to 07/12/2022
start_date2 <- min(hosp_data_IHME$date) # first case
end_date2 <- max(hosp_data_IHME$date) 
date_range2 <- seq(start_date2,end_date2,1)

hosp_data_IHME <- tibble(date = date_range2) %>%
  left_join(hosp_data_IHME, by = "date") %>%
  mutate(hosp = ifelse(is.na(value), NA, value))

# Imputation using Kalman Filter algorithm
hosp_data_IHME$hosp <- na_kalman(hosp_data_IHME$hosp)
hosp_data_IHME$hosp <- round(hosp_data_IHME$hosp)
hosp_data_IHME <- hosp_data_IHME %>%
  select(date, hosp) %>%
  dplyr::filter(date<= as.Date("2021-03-01"))

# --------------------------------------------------------------
# Estimation of hospitalisation rate
covid_data_BL <- inner_join(daily_cases,hosp_data_IHME, by="date")
covid_data_BL <- covid_data_BL %>%
  mutate(hosp2 = ifelse(hosp > new_cases, 1, hosp),
         p.hosp = hosp2 / new_cases)


# --------------------------------------------------------------

# --------------------------------------------------------------
t_step <- 0.25

# --------------------------------------------------------------
# Load model and plotting functions
source("nonVacc/R-nv/main_functionsBLvHosp-noVacc.R")
source("nonVacc/R-nv/plotting_functions-BLvHosp-noVacc.R")
source("nonVacc/R-nv/load_timeseries_data-BLvHosp-noVacc.R")

# --------------------------------------------------------------

# - - -
# Load model parameters
#thetaR_IC <- read_csv("inputs/theta_initial_conditions.csv")
#theta <- list(
#  r0=2.5, # note this is only IC - SMC estimates this
#  beta=NA,
#  betavol=volatility,
#  incubation = 1/6.4,
#  report = 1/6.1,
#  recover = 1/2.9, #2.9
#  init_cases=1,
#  pop=391483,   # https://sib.org.bz/census/2022-census/ - 405272 (World Bank 2023)
#  rep_prop=1, # propn reported - 0.0066   |0.6
# onset_prop=1, # propn onsets known  **1
#  confirmed_prop=1, # propn confirmed reported   (owid positivity rate)
# vacc_frac=1/0.583, # Estimate fraction that received a vaccine
#  r0_decline =1, # decline in R0 for scenario analysis
#  rep_var=1, # dispersion in local reporting confirmed cases - 6.1,
#  recover_hosp= NA, # rate of hospitalisation - 1/time of hospitalisation,
#  hosp_prop= 0.01,#sum(covid_data_BL$hosp)/sum(covid_data_BL$total_cases),#mean(covid_data_BL$p.hosp),  # propn hospitalised by infection
  #hosp_prop= mean(covid_data_BL$p.hosp),  # propn hospitalised by infection
#  death_hosp= 0.06
#)

#theta[["beta"]] <- theta[["r0"]]*(theta[["recover"]]) # Scale initial value of R0

#theta[["OR_vacc"]] <- 1/0.19

source("nonVacc/theta_nvacc.R")

theta_initNames <- c("sus","exp1","exp2","inf1","inf2","hosp","rec","waiting","cases","reports") # also defines groups to use in model


# - - -
# Load timeseries -  specify travel data being used
# NOTE: USES REPORTING DELAY AS INPUT



# Run set up check --------------------------------------------------------------

clusterExport(cl, varlist = c("smc_model", "theta","t_period","theta_initNames","vacc_time","process_model",
                              "AssignWeights","data_list"))

# - - -
# Run SMC and check likelihood
output_smc <- smc_model(theta,
                        nn=1e3, # number of particles
                        dt=t_step
)


output_smc$lik

# Run main outputs --------------------------------------------------------------

source("nonVacc/R-nv/outputs_main-BLvHosp-noVacc.R")















