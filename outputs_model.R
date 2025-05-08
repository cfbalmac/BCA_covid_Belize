# Benefit-Cost Analysis to assess a vaccination strategy for Belize
# Output script
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
library(parallel)

# --------------------------------------------------------------
# --------------------------------------------------------------
library(lubridate)
library(FinancialMath)
library(ggplot2)

rm(list=ls(all=TRUE))
# --------------------------------------------------------------
set.seed(1234)

num_cores <- detectCores() - 1  # Use all but one core to keep the system responsive
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, 123)  # Set seed for parallel processing

registerDoParallel(cl)

# --------------------------------------------------------------
# Parameters Costs and Benefits
# --------------------------------------------------------------
date_range <- seq(as.Date("2020-03-22"),as.Date("2023-06-30"),1)
pop <- 391483   # https://sib.org.bz/census/2022-census/ - 405272 (World Bank 2023) # Belize population

td.c <- 0.03
td.b <- 0.03

# Costs
# --------------------------------------------------------------

c.hosp.MV <- 41769/12.1 # USA - Di Fusco et al. (2021) from Gholipour 2023
c.hosp.UCI <- 10348.36/7  # IRAN - Ghaffari Darab et al (2021) from Gholipour 2023
c.hosp <- 2323.62/7 # IRAN - Ghaffari Darab et al (2021) from Gholipour 2023
c.oop <- 439*0.07  # IRAN - Ebrahimipour et al (2022) from Gholipour 2023 - 7% of population incur in oop expense
m.salary.average <- 1427*0.49  # BZ$1,427 (Apr-2024) 1BLZdollar = 0.49USD https://sib.org.bz/wp-content/uploads/LabourForce_2024-04.pdf
d.salary.average <- m.salary.average/20   # Divided by 20 because are the monthly days worked

gob.vacc_programme_2020 <- 1384622/2
gob.vacc_programme_2021 <- 16107526.58/2
gob.vacc_programme_2022 <- 16990686.70/2
gob.vacc_programme_2023 <- 5740053.06/2

c.gob_exp <- gob.vacc_programme_2020+gob.vacc_programme_2021+gob.vacc_programme_2022+gob.vacc_programme_2023

# Benefits
# --------------------------------------------------------------
th <- 365.25 # Time horizon in days (3.28 years).

u.overall <- 0.83/th # Utility overall population - Mao et al (2024)
u.Inf <- 0.79/th
u.Hosp <- 0.74/th
u.ICU <- 0.68/th
u.Rec <- 0.78/th # Assumption from Mao et al (2024). Not all the recovered patients will back to their normal life. Some of them has consequences. Therefore, we assume the utility of the overall population.

p.Hosp <- 0.73
p.UCI <- 0.27
p.MV <- p.UCI*0.13

n.days_hosp <- 8
n.days_hosp_nv <- 24.4

ndays_medleave <- 14+n.days_hosp
ndays_medleave_nv <- 14+n.days_hosp_nv


# Intervention scenario (Vaccination)
# --------------------------------------------------------------

source("main_script-BL_vHosp.R")
#source("R/main_functionsBL-Hosp.R")
#source("R/plotting_functions-BL-Hosp.R")
#source("R/load_timeseries_data-BL-Hosp.R")


filename="1"
load(paste0("bootstrap_fit_",filename,".RData"))

rep_plot_v <- 100; nn <- 2e3; cut_off=0; dt=0.25

out_rep_v <- foreach(kk = 1:rep_plot_v) %dopar% {
  output_smc <- smc_model(theta,nn,dt)
  output_smc
}

S_plot_v = matrix(NA,ncol=rep_plot_v,nrow=t_period)
I_plot_v = matrix(NA,ncol=rep_plot_v,nrow=t_period)
H_plot_v = matrix(NA,ncol=rep_plot_v,nrow=t_period)
R_plot_v = matrix(NA,ncol=rep_plot_v,nrow=t_period)
C_plot_v = matrix(NA,ncol=rep_plot_v,nrow=t_period)
Rep_plot_v = matrix(NA,ncol=rep_plot_v,nrow=t_period)
R0_plot_v = matrix(NA,ncol=rep_plot_v,nrow=t_period)
Rec_plot_v = matrix(NA,ncol=rep_plot_v,nrow=t_period)
CC_plot_v = matrix(NA,ncol=rep_plot_v,nrow=t_period)

for(kk in 1:rep_plot_v){
  output_smc <- out_rep_v[[kk]]
  if(output_smc$lik != - Inf){
    I_plot_v[,kk] <- output_smc$E_trace + output_smc$I_trace*(1-theta[["confirmed_prop"]])
    H_plot_v[,kk] <- output_smc$H_trace
    R_plot_v[,kk] <- output_smc$R_trace
    S_plot_v[,kk] <- output_smc$S_trace
    case_pos_v <- theta[["confirmed_prop"]]*theta[["rep_prop"]]*(output_smc$C_trace - c(0,head(output_smc$C_trace,-1)))
    C_plot_v[,kk] <- rpois(length(case_pos_v),lambda=case_pos_v)
    rep_pos_v <- theta[["confirmed_prop"]]*(output_smc$Rep_trace - c(0,head(output_smc$Rep_trace,-1))) # ALL CASES: theta[["rep_prop"]]*
    Rep_plot_v[,kk] <- rpois(length(rep_pos_v),lambda=rep_pos_v)
    R0_plot_v[,kk] <- output_smc$beta_trace/(theta[["recover"]])
    recov_case_v <- output_smc$R_trace - c(0,head(output_smc$R_trace,-1))
    Rec_plot_v[,kk] <- rpois(length(recov_case_v), lambda=recov_case_v)
    cum_case_pos_v <- theta[["confirmed_prop"]]*theta[["rep_prop"]]*(output_smc$C_trace)
    CC_plot_v[,kk] <- rpois(length(cum_case_pos_v),lambda=cum_case_pos_v)
  }
}


# Remove NA fits
S_plot_v = S_plot_v[,!is.na(S_plot_v[t_period,])]
I_plot_v = I_plot_v[,!is.na(I_plot_v[t_period,])]
H_plot_v = H_plot_v[,!is.na(H_plot_v[t_period,])]
R_plot_v = R_plot_v[,!is.na(R_plot_v[t_period,])]
C_plot_v = C_plot_v[,!is.na(C_plot_v[t_period,])]
Rep_plot_v = Rep_plot_v[,!is.na(Rep_plot_v[t_period,])]
R0_plot_v = R0_plot_v[,!is.na(R0_plot_v[t_period,])]
Rec_plot_v = Rec_plot_v[,!is.na(Rec_plot_v[t_period,])]
CC_plot_v = CC_plot_v[,!is.na(CC_plot_v[t_period,])]

# Calculate quantiles
S_quantile_v <- apply(S_plot_v,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)}) # number
S_quantile_n_average_v <- apply(S_plot_v,1,mean) # number
Inf_quantile_v <- apply(C_plot_v,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)}) # number
Inf_quantile_n_average_v <- apply(theta[["onset_prop"]]*C_plot_v,1,mean) # number
Inf_quantile_cum_v <- cumsum(Inf_quantile_n_average_v)
H_quantile_v <- apply(H_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)})
H_quantile_n_average_v <- apply(H_plot_v,1,mean)
H_quantile_cum_v <- cumsum(H_quantile_n_average_v)
R_quantile_v <- apply(Rec_plot_v,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)})
R_quantile_n_average_v <- apply(Rec_plot_v,1,mean)
R_quantile_cum_v <- cumsum(R_quantile_n_average_v)

trace_Vacc <- tibble(date= date_range,
                     Sv= S_quantile_n_average_v,
                     Iv = Inf_quantile_n_average_v,
                     Hv = H_quantile_n_average_v,
                     Rv = R_quantile_cum_v,
                     Dv = R_quantile_n_average_v*0.12,#theta[["death_hosp"]],
                     SuSv = pop-(Iv+Hv+Rv),
                     Totalv = SuSv+Iv+Hv+Rv,
                     C.hosp = Hv*(c.hosp*p.Hosp+c.hosp.UCI*p.UCI+c.hosp.MV*p.MV)*n.days_hosp,
                     C.oop = Iv*c.oop,  # out-of-pocket costs
                     C.exp_gob = ifelse(date <= as.Date("2020-12-31"), gob.vacc_programme_2020/as.numeric(as.Date("2020-12-31")-min(as.Date(date_range))),
                                        ifelse((date >= as.Date("2021-01-01") & date <= as.Date("2021-12-31")), gob.vacc_programme_2021/as.numeric(as.Date("2021-12-31")-as.Date("2021-01-01")),
                                               ifelse((date >= as.Date("2022-01-01") & date <= as.Date("2022-12-31")), gob.vacc_programme_2022/as.numeric(as.Date("2022-12-31")-as.Date("2022-01-01")),
                                                      gob.vacc_programme_2023/as.numeric(max(as.Date(date_range))-as.Date("2023-01-01"))))),
                     C.total = C.hosp+C.oop,#+C.exp_gob,
                     Year = as.numeric(format(date, "%Y")),
                     Year_Index = Year - min(Year),  # Compute year index relative to the minimum year in the dataset
                     C.total_disc = C.total / ((1 + td.c) ^ Year_Index),
                     Q.total = SuSv*u.overall+Iv*u.Inf+Hv*(p.Hosp*u.Hosp+p.UCI*u.ICU)+Rv*u.Rec,
                     Q.total_disc = Q.total / ((1 + td.b) ^ Year_Index),
                     QoLpc = Q.total/Totalv,
                     C.absv = d.salary.average*(Hv*ndays_medleave+Iv*14)
                     )

save(
  trace_Vacc,
  S_plot_v,
  I_plot_v,
  H_plot_v,
  R_plot_v,
  C_plot_v,
  Rep_plot_v,
  R0_plot_v,
  Rec_plot_v,
  CC_plot_v,
  S_quantile_v,
  S_quantile_n_average_v,
  Inf_quantile_v,
  Inf_quantile_n_average_v,
  Inf_quantile_cum_v,
  H_quantile_v,
  H_quantile_n_average_v,
  H_quantile_cum_v,
  R_quantile_v,
  R_quantile_n_average_v,
  R_quantile_cum_v,
  file=paste0("bootstrap_fit_vacc",filename,".RData")) 


ggplot(trace_Vacc)+
  #geom_line(aes(date,Sv))+
  geom_line(aes(date,Iv))+
  geom_line(aes(date,Hv))+
  #geom_line(aes(date,Rv))+
  geom_line(aes(date,Dv))

ggplot(trace_Vacc)+
  geom_line(aes(date,Iv))

# --------------------------------------------------------------
# --------------------------------------------------------------

# Contrafactual scenario (Without Vaccination)
# --------------------------------------------------------------

source("CBA-param.R")
source("nonVacc/main_script-BL_vHosp-noVacc.R")
#source("nonVacc/R-nv/main_functionsBLvHosp-noVacc.R")
#source("nonVacc/R-nv/plotting_functions-BLvHosp-noVacc.R")
#source("nonVacc/R-nv/load_timeseries_data-BLvHosp-noVacc.R")

filename="1"
load(paste0("nonVacc/outputs-nv/bootstrap_fit_",filename,".RData"))

#rep_plot_nv <- 100; nn <- 2e4; cut_off=0; dt=0.25

#out_rep_nv <- foreach(kk = 1:rep_plot_nv) %dopar% {
#  output_smc <- smc_model(theta,nn,dt)
#  output_smc
#}

S_plot_nv = matrix(NA,ncol=rep_plot_nv,nrow=t_period)
I_plot_nv = matrix(NA,ncol=rep_plot_nv,nrow=t_period)
H_plot_nv = matrix(NA,ncol=rep_plot_nv,nrow=t_period)
R_plot_nv = matrix(NA,ncol=rep_plot_nv,nrow=t_period)
C_plot_nv = matrix(NA,ncol=rep_plot_nv,nrow=t_period)
Rep_plot_nv = matrix(NA,ncol=rep_plot_nv,nrow=t_period)
R0_plot_nv = matrix(NA,ncol=rep_plot_nv,nrow=t_period)
Rec_plot_nv = matrix(NA,ncol=rep_plot_nv,nrow=t_period)
CC_plot_nv = matrix(NA,ncol=rep_plot_nv,nrow=t_period)

for(kk in 1:rep_plot_nv){
  output_smc <- out_rep_nv[[kk]]
  if(output_smc$lik != - Inf){
    I_plot_nv[,kk] <- output_smc$E_trace + output_smc$I_trace*(1-theta[["confirmed_prop"]])
    H_plot_nv[,kk] <- output_smc$H_trace
    R_plot_nv[,kk] <- output_smc$R_trace
    S_plot_nv[,kk] <- output_smc$S_trace
    case_pos_nv <- theta[["confirmed_prop"]]*theta[["rep_prop"]]*(output_smc$C_trace - c(0,head(output_smc$C_trace,-1)))
    C_plot_nv[,kk] <- rpois(length(case_pos_nv),lambda=case_pos_nv)
    rep_pos_nv <- theta[["confirmed_prop"]]*(output_smc$Rep_trace - c(0,head(output_smc$Rep_trace,-1))) # ALL CASES: theta[["rep_prop"]]*
    Rep_plot_nv[,kk] <- rpois(length(rep_pos_nv),lambda=rep_pos_nv)
    R0_plot_nv[,kk] <- output_smc$beta_trace/(theta[["recover"]])
    recov_case_nv <- output_smc$R_trace - c(0,head(output_smc$R_trace,-1))
    Rec_plot_nv[,kk] <- rpois(length(recov_case_nv), lambda=recov_case_nv)
    cum_case_pos_nv <- theta[["confirmed_prop"]]*theta[["rep_prop"]]*(output_smc$C_trace)
    CC_plot_nv[,kk] <- rpois(length(cum_case_pos_nv),lambda=cum_case_pos_nv)
  }
}


# Remove NA fits
S_plot_nv = S_plot_nv[,!is.na(S_plot_nv[t_period,])]
I_plot_nv = I_plot_nv[,!is.na(I_plot_nv[t_period,])]
H_plot_nv = H_plot_nv[,!is.na(H_plot_nv[t_period,])]
R_plot_nv = R_plot_nv[,!is.na(R_plot_nv[t_period,])]
C_plot_nv = C_plot_nv[,!is.na(C_plot_nv[t_period,])]
Rep_plot_nv = Rep_plot_nv[,!is.na(Rep_plot_nv[t_period,])]
R0_plot_nv = R0_plot_nv[,!is.na(R0_plot_nv[t_period,])]
Rec_plot_nv = Rec_plot_nv[,!is.na(Rec_plot_nv[t_period,])]
CC_plot_nv = CC_plot_nv[,!is.na(CC_plot_nv[t_period,])]

# Calculate quantiles
S_quantile_nv <- apply(S_plot_nv,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)}) # number
S_quantile_n_average_nv <- apply(S_plot_nv,1,function(x){mean(x,na.rm=TRUE)}) # number
Inf_quantile_nv <- apply(C_plot_nv,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)}) # number
Inf_quantile_n_average_nv <- apply(theta[["onset_prop"]]*C_plot_nv,1,function(x){mean(x,na.rm=TRUE)}) # number
Inf_quantile_cum_nv <- cumsum(Inf_quantile_n_average_nv)
H_quantile_nv <- apply(H_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)})
H_quantile_n_average_nv <- apply(H_plot_nv,1,function(x){mean(x,na.rm=TRUE)})
H_quantile_cum_nv <- cumsum(H_quantile_n_average_nv)
R_quantile_nv <- apply(Rec_plot_nv,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)})
R_quantile_n_average_nv <- apply(Rec_plot_nv,1,function(x){mean(x,na.rm=TRUE)})
R_quantile_cum_nv <- cumsum(R_quantile_n_average_nv)

trace_nVacc <- tibble(date= date_range,
                     Snv= S_quantile_n_average_nv,
                     Inv = Inf_quantile_n_average_nv,
                     Hnv = H_quantile_n_average_nv,
                     Rnv = R_quantile_cum_nv,
                     Dnv = R_quantile_n_average_nv*0.12*(1/0.86),#theta[["death_hosp"]],
                     SuSnv = pop-(Inv+Hnv+Rnv),
                     Totalnv = SuSnv+Inv+Hnv+Rnv,
                     C.hosp = Hnv*(c.hosp*p.Hosp+c.hosp.UCI*p.UCI+c.hosp.MV*p.MV)*n.days_hosp_nv,
                     C.oop = Inv*c.oop,  # out-of-pocket costs
                     C.exp_gob = 0,
                     C.total = C.hosp+C.oop+C.exp_gob,
                     Year = as.numeric(format(date, "%Y")),
                     Year_Index = Year - min(Year),  # Compute year index relative to the minimum year in the dataset
                     C.total_disc = C.total / ((1 + td.c) ^ Year_Index),
                     Q.total = SuSnv*u.overall+Inv*u.Inf+Hnv*(p.Hosp*u.Hosp+p.UCI*u.ICU)+Rnv*u.Rec,
                     Q.total_disc = Q.total / ((1 + td.b) ^ Year_Index),
                     QoLpc = Q.total/Totalnv,
                     C.absnv = d.salary.average*(Hnv*ndays_medleave_nv+Inv*14)
                     )


save(
  trace_nVacc,
  S_plot_nv,
  I_plot_nv,
  H_plot_nv,
  R_plot_nv,
  C_plot_nv,
  Rep_plot_nv,
  R0_plot_nv,
  Rec_plot_nv,
  CC_plot_nv,
  S_quantile_nv,
  S_quantile_n_average_nv,
  Inf_quantile_nv,
  Inf_quantile_n_average_nv,
  Inf_quantile_cum_nv,
  H_quantile_nv,
  H_quantile_n_average_nv,
  H_quantile_cum_nv,
  R_quantile_nv,
  R_quantile_n_average_nv,
  R_quantile_cum_nv,
  file=paste0("bootstrap_fit_nvacc",filename,".RData")) 

ggplot(trace_nVacc)+
  #geom_line(aes(date,Snv))+
  geom_line(aes(date,Inv))+
  geom_line(aes(date,Hnv))+
  #geom_line(aes(date,Rnv))+
  geom_line(aes(date,Dnv))


# --------------------------------------------------------------
# --------------------------------------------------------------
# Gathering datasets

# --------------------------------------------------------------

filename="1"
load(paste0("bootstrap_fit_vacc",filename,".RData"))
load(paste0("bootstrap_fit_nvacc",filename,".RData"))

times <- seq(1:nrow(trace_nVacc))
vsl.us <- 9400000   # VSL for US - World Bank 2024 (https://openknowledge.worldbank.org/entities/publication/e37f9361-51be-4603-8bfc-feb74b449b77)
GNI.us <- 57900   # World Bank 2023
GNI.BL <- 7190  # World Bank 2023
elasticity <- 1.2

gdp.pc_BL <- 6984.22    # World Bank 2022

vsl.bl <- vsl.us*((GNI.BL/GNI.us)^elasticity)
#iwtp.vrr <- sum((trace_nVacc$C.hosp.i-trace_Vacc$C.hosp.i),na.rm = TRUE)
iwtp.vrr <- sum(trace_nVacc$C.hosp-trace_Vacc$C.hosp)+sum(trace_nVacc$C.oop-trace_Vacc$C.oop)+sum(trace_nVacc$C.absnv-trace_Vacc$C.absv)
iwtp.vrr <- iwtp.vrr/pop

ce.threshold <- median(c(0.1,0.64))*gdp.pc_BL   # Median value for % range of the threshold proposed by Ochalek et al (iDSI 2015)
#ce.threshold <- 0.49*gdp.pc_BL   #% range of the threshold proposed by Pichon-Rivieri 2023

trace_diff <- tibble(date= date_range,
                     Year = as.numeric(format(date, "%Y")),
                     Year_Index = Year - min(Year),  # Compute year index relative to the minimum year in the dataset
                     SuSdiff = trace_Vacc$SuSv-trace_nVacc$SuSnv,
                     Idiff = trace_Vacc$Iv-trace_nVacc$Inv,
                     Hdiff = trace_Vacc$Hv-trace_nVacc$Hnv,
                     Rdiff = trace_Vacc$Rv-trace_nVacc$Rnv,
                     Ddiff = trace_Vacc$Dv-trace_nVacc$Dnv,
                     Cdiff.exp_gob = trace_Vacc$C.exp_gob-trace_nVacc$C.exp_gob, # Incremental Government Expenses - Programme cost
                     Cdiff.hosp = trace_Vacc$C.hosp-trace_nVacc$C.hosp, # Incremental Hospitalisation Cost 
                     Cdiff.oop = trace_Vacc$C.oop-trace_nVacc$C.oop,# Incremental Out-of-Pocket Cost 
                     C.abs = trace_Vacc$C.absv-trace_nVacc$C.absnv, # Incremental Absenteeism Cost
                     Cdiff.total = (trace_Vacc$C.total - trace_nVacc$C.total) + C.abs, # Undiscounted Incremental Total Cost
                     Cdiff.total_disc = Cdiff.total/(1+td.c)^Year_Index, # Discounted Incremental Total Cost
                     QoLdiff = trace_Vacc$Q.total-trace_nVacc$Q.total, # Undiscounted Incremental Total QALYs
                     QoLdiff.pc = QoLdiff/pop, # Undiscounted Incremental Total QALYs per capita
                     QoL = QoLdiff*ce.threshold, # Undiscounted Incremental Total QALYs Monetised
                     B.premD = -Ddiff*m.salary.average, # Premature Death Costs - negative because we need to use the averted deaths (positive)
                     VRR = -Hdiff*iwtp.vrr,  # Value per Nonfatal Health Risk Reductions - negative because we need to use the averted hospitalisations (positive)
                     VSL = -Ddiff*vsl.bl, # Value per Statistical Life - negative because we need to use the averted deaths (positive)
                     Qdiff.total = VRR + VSL + B.premD, # Undiscounted Total Monetised Benefits - for base case use only
                     Qdiff.total_disc = Qdiff.total/(1+td.b)^Year_Index, # Discounted Total Monetised Benefits - for base case use only
                     B.premD_disc = B.premD/(1+td.b)^Year_Index, # Discounted Premature Death Costs - for scenario analysis use only
                     QoL_disc = QoL/(1+td.b)^Year_Index, # Discounted Total Monetised QALYs - for scenario analysis use only,
                     VSL_disc = VSL/(1+td.b)^Year_Index # Discounted Value per Statistical Life - for scenario analysis use only,
                     )


#ggplot(trace_diff)+
#  geom_line(aes(date,Idiff))

#ggplot(trace_diff)+
#  geom_line(aes(date,Ddiff))


# --------------------------------------------------------------
# BCA Results - Benefits as averted events
# --------------------------------------------------------------
# --------------------------------------------------------------

benefits <- sum(trace_diff$Qdiff.total_disc,na.rm = TRUE)
costs <- sum(trace_diff$Cdiff.total_disc,na.rm = TRUE)+gob.vacc_programme_2020/(1+td.c)^0+
  gob.vacc_programme_2021/(1+td.c)^1+
  gob.vacc_programme_2022/(1+td.c)^2+
  gob.vacc_programme_2023/(1+td.c)^3

#results
npv <- benefits - costs
bc.ratio <- benefits/costs

cash_flows <- trace_diff$Qdiff.total_disc - trace_diff$Cdiff.total_disc
cash_flows_yearly <- unname(tapply(cash_flows, format(date_range,"%Y"),sum)) # From 2020 to 2023
cash_flows_yearly <- as.vector(cash_flows_yearly)
cash_flows_yearly_mirr <- c(-c.gob_exp,cash_flows_yearly)
#cash_flows <- cash_flows[!is.na(cash_flows)]
#cash_flows2 <- c(-c.gob_exp,trace_diff$Qdiff.total_disc - (trace_diff$Cdiff.total_disc))

#results_irr <- tibble(date=date_range,
#                  cash_flows = cash_flows)
#write.csv(results_irr,"results_irr.csv")

# --------------------------------------------------------------
# IRR estimation

irr.bc <- IRR(cf0=-c.gob_exp,cf=cash_flows_yearly,times = seq(1:4))
irr.bc

# --------------------------------------------------------------
# modified IRR

# Define the rates
finance_rate <- 0.04  # Financing rate
reinvest_rate <- 0.05  # Reinvestment rate

# Number of periods
n <- length(cash_flows_yearly_mirr) - 1

# Calculate the present value of negative cash flows
PV_negatives <- sum(cash_flows_yearly_mirr[cash_flows_yearly_mirr < 0] / (1 + finance_rate) ** (seq_along(cash_flows_yearly_mirr[cash_flows_yearly_mirr < 0]) - 1))

# Calculate the future value of positive cash flows
FV_positives <- sum(cash_flows_yearly_mirr[cash_flows_yearly_mirr > 0] * (1 + reinvest_rate) ** (n - seq_along(cash_flows_yearly_mirr[cash_flows_yearly_mirr > 0]) + 1))

# Calculate the MIRR
mirr_bc <- (FV_positives / abs(PV_negatives)) ** (1 / n) - 1
mirr_bc*100


# --------------------------------------------------------------
# Overall results

sum(trace_diff$Idiff)   # Averted Infections
sum(trace_diff$Hdiff)   # Averted Hospitalisations
sum(trace_diff$Ddiff)   # Averted Deaths

sum(trace_diff$B.premD)   # Averted premature death costs

sum(trace_diff$Cdiff.hosp)   # Averted hospitalisation costs
sum(trace_diff$Cdiff.oop)   # Averted out-of-cost costs

sum(trace_diff$C.abs)   # Averted absenteeism costs

ggplot()+
  geom_line(aes(date_range,cash_flows),colour="blue",linewidth=1.5)+
  ylim(-90000000,250000000)+
  geom_hline(yintercept=0, colour="red",linetype="dashed")+
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y")+
  theme_test() +
  ggtitle("NPV cash flows")

yearly_time <- c(2020,2021,2022,2023)

ggplot()+
  geom_col(aes(yearly_time,cash_flows_yearly),colour="blue", fill="blue",width = 0.25)+
  ylim(-2600000000,17000000000)+
  geom_hline(yintercept=0, colour="red",linetype="dashed")+
  theme_test() +
  ggtitle("NPV cash flows")

# --------------------------------------------------------------
# Scenario Analysis (Benefits using QALYs)
# --------------------------------------------------------------
# --------------------------------------------------------------

benefits_sa1 <- sum(trace_diff$QoL_disc,na.rm = TRUE)+sum(trace_diff$B.premD_disc)+sum(trace_diff$VSL_disc)
costs_sa1 <- sum(trace_diff$Cdiff.total_disc,na.rm = TRUE)+gob.vacc_programme_2020/(1+td.c)^0+
  gob.vacc_programme_2021/(1+td.c)^1+
  gob.vacc_programme_2022/(1+td.c)^2+
  gob.vacc_programme_2023/(1+td.c)^3

# Results
npv_sa1 <- benefits_sa1 - costs_sa1
br.ratio_sa1 <- benefits_sa1/costs_sa1


cash_flows_sa1 <- (trace_diff$QoL_disc+trace_diff$B.premD_disc) - (trace_diff$Cdiff.total_disc)
cash_flows_yearly_sa1 <- unname(tapply(cash_flows_sa1, format(date_range,"%Y"),sum)) # From 2020 to 2023
cash_flows_yearly_sa1 <- as.vector(cash_flows_yearly_sa1)
cash_flows_yearly_mirr_sa1 <- c(-c.gob_exp,cash_flows_yearly_sa1)


#write.csv(results_irr,"results_irr.csv")

# --------------------------------------------------------------
# IRR estimation
irr.sa1 <- IRR(cf0=-c.gob_exp,cf=cash_flows_yearly_sa1,times = seq(1:4))
irr.sa1

# --------------------------------------------------------------
# modified IRR

# Calculate the present value of negative cash flows
PV_negatives_sa <- sum(cash_flows_yearly_mirr_sa1[cash_flows_yearly_mirr_sa1 < 0] / (1 + finance_rate) ** (seq_along(cash_flows_yearly_mirr_sa1[cash_flows_yearly_mirr_sa1 < 0]) - 1))

# Calculate the future value of positive cash flows
FV_positives_sa <- sum(cash_flows_yearly_mirr_sa1[cash_flows_yearly_mirr_sa1 > 0] * (1 + reinvest_rate) ** (n - seq_along(cash_flows_yearly_mirr_sa1[cash_flows_yearly_mirr_sa1 > 0]) + 1))

# Calculate the MIRR
mirr_sa <- (FV_positives_sa / abs(PV_negatives_sa)) ** (1 / n) - 1
mirr_sa*100


# --------------------------------------------------------------
# Overall results

sum(trace_diff$QoLdiff)   # Incremental Total QoL
sum(trace_diff$QoLdiff.pc)  # Incremental QoL per capita

sum(trace_diff$QoL_disc)   # Incremental Total QoL
sum(trace_diff$QoLdiff.pc)*ce.threshold

ggplot()+
  geom_line(aes(date_range,cash_flows_sa1),colour="blue",linewidth=1.5)+
  #ylim(-90000000,200000000)+
  geom_hline(yintercept=0, colour="red",linetype="dashed")+
  scale_x_date(date_breaks = "3 months" , date_labels = "%b-%y")+
  theme_test() +
  ggtitle("NPV cash flows")

ggplot()+
  geom_col(aes(yearly_time,cash_flows_yearly),colour="blue", fill="blue",width = 0.25)+
  ylim(-2600000000,17000000000)+
  geom_hline(yintercept=0, colour="red",linetype="dashed")+
  theme_test() +
  ggtitle("NPV cash flows")



