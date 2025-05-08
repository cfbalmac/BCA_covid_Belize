

#rm(list=ls(all=TRUE))
# --------------------------------------------------------------
# --------------------------------------------------------------
# --------------------------------------------------------------
# Parameters Costs and Benefits
# --------------------------------------------------------------
date_range <- seq(as.Date("2020-03-22"),as.Date("2023-06-30"),1)
pop <- 391483   # https://sib.org.bz/census/2022-census/ - 405272 (World Bank 2023) # Belize population

td.c <- 0.05
td.b <- 0.05

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

p.Hosp <- 0.79 # 1-p.UCI
p.UCI <- 0.21 # Chang et al 2021 COVID-19 ICU and mechanical ventilation patient characteristics and outcomes — A systematic review and meta-analysis
p.MV <- p.UCI*0.69 # Chang et al 2021 COVID-19 ICU and mechanical ventilation patient characteristics and outcomes — A systematic review and meta-analysis

n.days_hosp <- 11.91     # Fang et al 2024 - Vaccines reduced hospital length of stay and fraction of inspired oxygen of COVID-19 patients: A retrospective cohort study 
n.days_hosp_nv <- 14.26  # Fang et al 2024 - Vaccines reduced hospital length of stay and fraction of inspired oxygen of COVID-19 patients: A retrospective cohort study 

ndays_medleave <- 14+n.days_hosp
ndays_medleave_nv <- 14+n.days_hosp_nv
