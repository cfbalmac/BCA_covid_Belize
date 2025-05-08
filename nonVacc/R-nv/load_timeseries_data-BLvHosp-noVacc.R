 # Define values
pre_peak <- 3 # -1 is 2 before peak, 2 is 2 after
omit_recent <- 0
omit_conf <- 0

# Set up start/end dates
new_date_bl <- as.Date("2023-06-30")
start_date <- as.Date("2020-03-22") # first case

#end_date <- max(case_data_in$date) # omit recent day?
end_date <- new_date_bl #as.Date("2020-03-03") # period to forecast ahead
date_range <- seq(start_date,end_date,1)

# When restrictions started
Bl_vacc <- as.Date("2021-03-01")
vacc_time <- as.numeric(Bl_vacc - start_date + 1)

fix_r0_tt <- as.numeric(new_date_bl - start_date + 1) #as.Date("2020-03-18") 

###########################################################

t_period <- as.numeric(end_date-start_date)+1


################################################################################

# Load Belize onset data --------------------------------------------

case_data_bl <- daily_cases
cutoff_time_bl <- max(case_data_bl$date) - omit_recent # omit final days of time points
case_data_bl[case_data_bl$date>cutoff_time_bl,"new_cases"] <- NA
case_data_bl_time <- rep(0,length(date_range))

# ensure final points are omitted
for(ii in 1:length(date_range)){
  case_data_bl_time[ii] <- sum(case_data_bl[case_data_bl$date==date_range[ii],]$new_cases)
}
case_data_bl_time[date_range>cutoff_time_bl] <- NA # omit final points


# tally cases
case_time <- rep(0,length(date_range))

for(ii in 1:length(date_range)){
  case_time[ii] = sum(case_data_bl[case_data_bl$date==date_range[ii],]$new_cases)
}

################################################################################
################################################################################

################################################################################
################################################################################

################################################################################

################################################################################

# Extract Belize data --------------------------------------------------

cases_BL <- daily_cases
cases_BL <-  cases_BL %>% mutate(new_cases = NA)
cases_BL$new_cases <- cases_BL$total_cases - c(NA,head(cases_BL$total_cases,-1)) 

case_data_bl_conf_time <- rep(NA,length(date_range))
cutoff_time_bl <- max(cases_BL$date)
cutoff_min_bl <- as.Date("2020-03-22")

for(ii in 1:length(date_range)){
  case_data_bl_conf_time[ii] = sum(cases_BL[cases_BL$date==date_range[ii],]$new_cases)
}

case_data_bl_conf_time[date_range>cutoff_time_bl | date_range<cutoff_min_bl] <- NA # omit all old and recent points

# Remove outlier:
case_data_bl_conf_time[case_data_bl_conf_time>1e4] <- NA

# NEED REPORTING PARAMETER

# Load Belize hospitalisation data --------------------------------------------

hosp_data <- hosp_data_IHME
cutoff_time_hosp <- max(hosp_data$date) - omit_recent # omit final days of time points
hosp_data[hosp_data$date>cutoff_time_hosp,"new_cases"] <- NA
hosp_data_time <- rep(0,length(date_range))

# ensure final points are omitted
for(ii in 1:length(date_range)){
  hosp_data_time[ii] <- sum(hosp_data[hosp_data$date==date_range[ii],]$hosp)
}
hosp_data_time[date_range>cutoff_time_hosp] <- NA # omit final points


# tally cases
case_time <- rep(0,length(date_range))

for(ii in 1:length(date_range)){
  case_time[ii] = sum(case_data_bl[case_data_bl$date==date_range[ii],]$new_cases)
}


# Compile list of data to use:
data_list <- list(case_data_onset = case_data_bl_time,
                  case_data_conf = case_data_bl_conf_time,
                  hospital_data = hosp_data_time)



# Quick plot
#plot(case_data_belize$date,case_data_belize$new_cases,xlim=as.Date(c("2020-03-22","2023-06-30")),ylim=c(0,500)); points(case_data_belize$date,case_data_belize$new_cases,col="blue")


