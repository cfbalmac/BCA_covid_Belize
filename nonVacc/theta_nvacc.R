# parameters non-vaccinated scenario

theta <- list(
  r0=2.5, # note this is only IC - SMC estimates this
  beta=NA,
  betavol=volatility,
  incubation = 1/5.2,
  report = 1/6.1,
  recover = 1/2.9, #2.9
  init_cases=1,
  pop=391483,   # https://sib.org.bz/census/2022-census/ - 405272 (World Bank 2023)
  rep_prop=1, # propn reported - 0.0066   |0.6
  onset_prop=1, # propn onsets known  **1
  confirmed_prop=1, # propn confirmed reported   (owid positivity rate)
  vacc_frac=1/0.583, # Estimate fraction that received a vaccine
  r0_decline =1, # decline in R0 for scenario analysis
  rep_var=1, # dispersion in local reporting confirmed cases - 6.1,
  recover_hosp= NA, # rate of hospitalisation - 1/time of hospitalisation,
  hosp_prop= 0.01,#sum(covid_data_BL$hosp)/sum(covid_data_BL$total_cases),#mean(covid_data_BL$p.hosp),  # propn hospitalised by infection
  #hosp_prop= mean(covid_data_BL$p.hosp),  # propn hospitalised by infection
  death_hosp= 0.06
)

theta[["beta"]] <- theta[["r0"]]*(theta[["recover"]]) # Scale initial value of R0

theta[["OR_vacc"]] <- 1/0.19