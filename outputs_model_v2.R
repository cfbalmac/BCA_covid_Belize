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

# --------------------------------------------------------------
# --------------------------------------------------------------
library(lubridate)
library(FinancialMath)
library(ggplot2)

rm(list=ls(all=TRUE))
# --------------------------------------------------------------
#set.seed(1234)

num_cores <- detectCores() - 1  # Use all but one core to keep the system responsive
cl <- makeCluster(num_cores)
clusterSetRNGStream(cl, 123)  # Set seed for parallel processing

registerDoParallel(cl)


#stopCluster(cl)


# Intervention scenario (Vaccination)
# --------------------------------------------------------------

source("CBA-param.R")
source("main_script-BL_vHosp.R")
#source("R/main_functionsBL-Hosp.R")
#source("R/plotting_functions-BL-Hosp.R")
#source("R/load_timeseries_data-BL-Hosp.R")

filename="1_vacc"
load(paste0("Vacc/outputs_vacc/results_",filename,".RData"))


trace_Vacc <- tibble(date= date_range,
                     Sv= S_quantile_n_average,
                     Iv = Inf_quantile_n_average,
                     Hv = H_quantile_n_average,
                     Rv = R_quantile_cum,
                     Dv = R_quantile_n_average*0.12,#theta[["death_hosp"]],
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
  S_plot,
  I_plot,
  H_plot,
  R_plot,
  C_plot,
  Rep_plot,
  R0_plot,
  Rec_plot,
  CC_plot,
  S_quantile,
  S_quantile_n_average,
  Inf_quantile,
  Inf_quantile_n_average,
  Inf_quantile_cum,
  H_quantile,
  H_quantile_n_average,
  H_quantile_cum,
  R_quantile,
  R_quantile_n_average,
  R_quantile_cum,
  file=paste0("Vacc/outputs_vacc/bootstrap_fit_vacc",filename,".RData")) 


# --------------------------------------------------------------
# --------------------------------------------------------------

# Contrafactual scenario (Without Vaccination)
# --------------------------------------------------------------

source("CBA-param.R")
source("nonVacc/main_script-BL_vHosp-noVacc.R")
#source("nonVacc/R-nv/main_functionsBLvHosp-noVacc.R")
#source("nonVacc/R-nv/plotting_functions-BLvHosp-noVacc.R")
#source("nonVacc/R-nv/load_timeseries_data-BLvHosp-noVacc.R")

filename="1_nvacc"
load(paste0("nonVacc/outputs-nv/results_",filename,".RData"))

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
  file=paste0("nonVacc/outputs-nv/bootstrap_fit_nvacc",filename,".RData")) 


# --------------------------------------------------------------
# --------------------------------------------------------------
# Gathering datasets

# --------------------------------------------------------------

filename="1_vacc"
filename1="1_nvacc"
load(paste0("Vacc/outputs_vacc/bootstrap_fit_vacc",filename,".RData"))
load(paste0("nonVacc/outputs-nv/bootstrap_fit_nvacc",filename1,".RData"))

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
  ylim(-3000000000,20000000000)+
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
  ylim(-3000000000,20000000000)+
  geom_hline(yintercept=0, colour="red",linetype="dashed")+
  theme_test() +
  ggtitle("NPV cash flows")


stopCluster(cl)

