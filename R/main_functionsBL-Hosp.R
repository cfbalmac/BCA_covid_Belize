# Process model for simulation --------------------------------------------


# --------------------------------------------------------------

process_model <- function(t_start,t_end,dt,theta,simTab,simzetaA,vaccF,recover_hosp){
  
  # simTab <- storeL[,tt-1,]; t_start = 1; t_end = 2; dt = 0.1; simzetaA <- simzeta[1,]; vaccF=theta[["vacc_frac"]]
  
  susceptible_t <- simTab[,"sus"] # input function
  exposed_t1 <- simTab[,"exp1"] # input function
  exposed_t2 <- simTab[,"exp2"] # input function

  infectious_t1 <- simTab[,"inf1"] # input function
  infectious_t2 <- simTab[,"inf2"] # input function
  
  hospitalized_t <- simTab[,"hosp"]  # New compartment for hospitalized individuals
  recovered_t <- simTab[,"rec"]      # Adjusting for recovery from hospitalization
  
  waiting_t <- simTab[,"waiting"] # input function
  cases_t <- simTab[,"cases"] # input function
  reports_t <- simTab[,"reports"] # input function
  
  # scale transitions
  inf_rate <- (simzetaA/theta[["pop"]])*dt
  inc_rate <- theta[["incubation"]]*2*dt
  rec_rate <- theta[["recover"]]*dt
  rec_rate_hosp <- 1/recover_hosp * dt    # Recovery rate for hospitalized
  hosp_rate <- theta[["hosp_prop"]] * dt         # Hospitalization rate
  death_rate_hosp <- theta[["death_hosp"]] * dt    # Death rate from hospitalization
  
  rep_rate <- theta[["report"]]*dt
  vacc_frac <- vaccF
  prob_rep <- exp(-theta[["report"]]*theta[["recover"]]) # probability case is reported rather than recovers

  for(ii in seq((t_start+dt),t_end,dt) ){
    
    # transitions
    S_to_E1 <- susceptible_t*(infectious_t1+infectious_t2)*inf_rate # stochastic transmission
    
    # Delay until symptoms
    E1_to_E2 <- exposed_t1*inc_rate # as two compartments
    E2_to_I1 <- exposed_t2*inc_rate
    
    # Infectious to Hospitalized and Recovered
    I1_to_I2 <- infectious_t1 * rec_rate
    I2_to_Rec <- infectious_t2 * rec_rate 
    I2_to_Hosp <- infectious_t2 * hosp_rate
    Hosp_to_Rec <- hospitalized_t * rec_rate_hosp
    Hosp_to_Death <- hospitalized_t * death_rate_hosp
    
    # Delay until recovery
    #I1_to_I2 <- infectious_t1*rec_rate
    #I2_to_R <- infectious_t2*rec_rate
    
    # Delay until reported
    W_to_Rep <- waiting_t*rep_rate
    
    # Process model for SEIR
    susceptible_t <- susceptible_t - S_to_E1
    exposed_t1 <- exposed_t1 + S_to_E1*(1-vacc_frac) - E1_to_E2
    exposed_t2 <- exposed_t2 + E1_to_E2 - E2_to_I1

    infectious_t1 <- infectious_t1 + E2_to_I1 - I1_to_I2
    infectious_t2 <- infectious_t2 + I1_to_I2 - I2_to_Hosp - I2_to_Rec
    
    hospitalized_t <- hospitalized_t + I2_to_Hosp - Hosp_to_Rec - Hosp_to_Death
    recovered_t <- recovered_t + Hosp_to_Rec + Hosp_to_Death + I2_to_Rec
    
    # Case tracking - including removal of cases within Q compartment
    waiting_t <- waiting_t + E2_to_I1*prob_rep - W_to_Rep
    cases_t <- cases_t + E2_to_I1*prob_rep
    reports_t <- reports_t + W_to_Rep
  }
  
  simTab[,"sus"] <- susceptible_t # output
  simTab[,"exp1"] <- exposed_t1 # output
  simTab[,"exp2"] <- exposed_t2 # output
  simTab[,"inf1"] <- infectious_t1 # output
  simTab[,"inf2"] <- infectious_t2 # output
  simTab[,"hosp"] <- hospitalized_t
  simTab[,"rec"] <- recovered_t
  
  simTab[,"waiting"] <- waiting_t # output
  simTab[,"cases"] <- cases_t # output
  simTab[,"reports"] <- reports_t # output

  simTab
  
}

# --------------------------------------------------------------

# SMC function -------------------------------------------------
# --------------------------------------------------------------

smc_model <- function(theta,nn,dt=1,seed=NULL){
  if (!is.null(seed)) set.seed(seed)
  
  # nn = 100;   dt <- 0.25
  
  # Assumptions - using daily growth rate
  ttotal <- t_period
  t_length <- ttotal
  
  storeL <- array(0,dim=c(nn,t_length, length(theta_initNames)),dimnames = list(NULL,NULL,theta_initNames))
  
  # Add initial condition
  #storeL[,1,"exp1"] <- theta[["init_cases"]]
  #storeL[,1,"exp2"] <- theta[["init_cases"]]
  storeL[,1,"inf1"] <- theta[["init_cases"]]/2
  storeL[,1,"inf2"] <- theta[["init_cases"]]/2
  storeL[,1,"sus"] <- theta[["pop"]] - theta[["init_cases"]]
  
  #simzeta <- matrix(rlnorm(nn*t_length, mean = -theta[["betavol"]]^2/2, sd = theta[["betavol"]]),ncol=ttotal)
  simzeta <- matrix(rnorm(nn*t_length, mean = 0, sd = theta[["betavol"]]),nrow=ttotal)
  simzeta[1,] <- exp(simzeta[1,])*theta[["beta"]] # define IC
  
  # Fix R for forward simulation
  #simzeta[fix_r0_tt:ttotal,] <- log(theta[["r0_decline"]])
  
  # Latent variables
  S_traj = matrix(NA,ncol=1,nrow=ttotal)
  C_traj = matrix(NA,ncol=1,nrow=ttotal)
  Rep_traj = matrix(NA,ncol=1,nrow=ttotal)
  E_traj = matrix(NA,ncol=1,nrow=ttotal)
  I_traj = matrix(NA,ncol=1,nrow=ttotal)
  H_traj = matrix(NA,ncol=1,nrow=ttotal)
  R_traj = matrix(NA,ncol=1,nrow=ttotal)
  beta_traj = matrix(NA,ncol=1,nrow=ttotal);
  w <- matrix(NA,nrow=nn,ncol=ttotal); w[,1] <- 1  # weights
  W <- matrix(NA,nrow=nn,ncol=ttotal)
  A <- matrix(NA,nrow=nn,ncol=ttotal) # particle parent matrix
  l_sample <- rep(NA,ttotal)
  lik_values <- rep(NA,ttotal)
  
  # Iterate through steps
  
  for(tt in 2:ttotal){
    
    # DEBUG  tt=2
    
    # Add random walk on transmission ?
    simzeta[tt,] <- simzeta[tt-1,]*exp(simzeta[tt,])
    
    # vaccination?
    if(tt<vacc_time){vaccF <- theta[["vacc_frac"]]}else{vaccF <- 0}

    if(tt<vacc_time){recover_hosp <- 14.26}else{recover_hosp <- 11.91}
    
    
    # run process model
    storeL[,tt,] <- process_model(tt-1,tt,dt,theta,storeL[,tt-1,],simzeta[tt,],vaccF,recover_hosp)
    
    # calculate weights
    w[,tt] <- AssignWeights(data_list,storeL,nn,theta,tt)
    
    # check likelihood isn't NA
    if(is.na(max(w[1:nn,tt])) | max(w[1:nn,tt]) == 0){
      likelihood0 = -Inf
      return(list(S_trace=S_traj,
                  C_trace=C_traj,
                  Rep_trace=Rep_traj,
                  I_trace=I_traj,
                  H_trace=H_traj,
                  R_trace=R_traj,
                  beta_trace=beta_traj,
                  lik=likelihood0 ))
    }
    
    #c_local_val,c_val,rep_val,local_case_data_tt,case_data_tt,rep_data_tt,nn){
    
    # normalise particle weights
    sum_weights <- sum(w[1:nn,tt])
    W[1:nn,tt] <- w[1:nn,tt]/sum_weights
    
    # resample particles by sampling parent particles according to weights:
    A[, tt] <- sample(1:nn,prob = W[1:nn,tt],replace = T)
    
    # DEPRECATED
    # for (j in 1:nn){
    #   locs <- pickA[cumsum_W >= rand_vals[j]]; A[j, tt] <- locs[1]
    # }
    
    # Resample particles for corresponding variables
    storeL[,tt,] <- storeL[ A[, tt] ,tt,]
    simzeta[tt,] <- simzeta[tt, A[, tt]] #- needed for random walk on beta
    
    
  } # END PARTICLE LOOP
  
  
  # Estimate likelihood:
  for(tt in 1:ttotal){
    lik_values[tt] = log(sum(w[1:nn,tt])) # log-likelihoods
  }
  
  likelihood0 = -ttotal*log(nn)+ sum(lik_values) # add full averaged log-likelihoods
  
  # Sample latent variables:
  locs <-  sample(1:nn,prob = W[1:nn,tt],replace = T)
  l_sample[ttotal] <- locs[1]
  C_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"cases"]
  Rep_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"reports"]
  S_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"sus"]
  E_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"exp2"] +storeL[l_sample[ttotal],ttotal,"exp1"]
  I_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"inf1"]+storeL[l_sample[ttotal],ttotal,"inf2"]
  H_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"hosp"]
  R_traj[ttotal,] <- storeL[l_sample[ttotal],ttotal,"rec"]
  beta_traj[ttotal,] <- simzeta[ttotal,l_sample[ttotal]]
  
  for(ii in seq(ttotal,2,-1)){
    l_sample[ii-1] <- A[l_sample[ii],ii] # have updated indexing
    C_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"cases"]
    Rep_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"reports"]
    S_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"sus"]
    E_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"exp2"]+ storeL[l_sample[ii-1],ii-1,"exp1"]
    I_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"inf1"]+ storeL[l_sample[ii-1],ii-1,"inf2"]
    H_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"hosp"]
    R_traj[ii-1,] <- storeL[l_sample[ii-1],ii-1,"rec"]
    beta_traj[ii-1,] <- simzeta[ii-1,l_sample[ii-1]]
  }
  
  # DEBUG  plot(Rep_traj[,1]-C_traj[,1])
  
  
  return(list(S_trace=S_traj,
              C_trace=C_traj,
              Rep_trace=Rep_traj,
              E_trace=E_traj,
              I_trace=I_traj,
              H_trace=H_traj,
              R_trace=R_traj,
              beta_trace=beta_traj,
              lik=likelihood0 ))
  
  
}

# --------------------------------------------------------------

# Likelihood calc for SMC --------------------------------------------
# --------------------------------------------------------------

AssignWeights <- function(data_list,storeL,nn,theta,tt){
  
  # c_local_val=case_localDiff; c_val=caseDiff; rep_val=repDiff;
  
  # Gather data
  local_case_data_tt <- data_list$case_data_onset[tt]
  local_case_data_conf_tt <- data_list$case_data_conf[tt]
  
  local_case_hosp <- data_list$hospital_data[tt]

  # Scale for reporting lag
  local_case_data_tt_scale <- 1#data_list$local_case_data_onset_scale[tt] # deprecated
  
  # Gather variables
  hosp_t <- storeL[, tt, "hosp"]
  rep <- storeL[,tt,"reports"]
  caseDiff <- storeL[,tt,"cases"] - storeL[,tt-1,"cases"]
  repDiff <- storeL[,tt,"reports"] - storeL[,tt-1,"reports"]
  hospDiff <- storeL[,tt,"hosp"] - storeL[,tt-1,"hosp"]
  
  # Prevalence - scale by asymptomatics - second half only // storeL[,tt,"exp1"] + 
  inf_prev <- storeL[,tt,"exp1"] + storeL[,tt,"exp2"] + (storeL[,tt,"inf1"] + storeL[,tt,"inf2"])*(1-theta[["confirmed_prop"]])
  
  c_val <- pmax(0,caseDiff)
  rep_val <- pmax(0,repDiff) # NOTE CHECK FOR POSITIVITY
  r_rep_cum <- rep
  hosp_val <- pmax(0,hospDiff)
  
  # Local confirmed cases (by onset)
  
  if(!is.na(local_case_data_tt)){
    expected_val <- c_val* theta[["confirmed_prop"]]*theta[["onset_prop"]]*theta[["rep_prop"]]*local_case_data_tt_scale # scale by reporting proportion and known onsets
    loglikSum_local_onset <- dnbinom(x=local_case_data_tt, mu=expected_val,size = 1/theta[["rep_var"]],log = T)
    #loglikSum_local_onset <- dpois(local_case_data_tt,lambda = expected_val,log=T)
  }else{
    loglikSum_local_onset <- 0
  }
  
  # Local confirmed cases (by confirmation) -- HOLD OUT FOR NOW AS LIKELIHOOD LOW
  
  if(!is.na(local_case_data_conf_tt)){
    expected_val <- rep_val * theta[["confirmed_prop"]] * theta[["rep_prop"]] # scale by reporting proportion and known onsets
    #loglikSum_local_conf <- dpois(local_case_data_conf_tt,lambda = expected_val,log=T)
    
    loglikSum_conf <- dnbinom(x=local_case_data_conf_tt,mu=expected_val,size=1/theta[["rep_var"]],log=T)
    
    #print(loglikSum_conf[1])
    
  }else{
    loglikSum_conf <- 0
  }
  
  # Calculate log-likelihood of observed vs. expected hospitalizations
  if (!is.na(local_case_hosp)) {
    expected_hospitalizations <- hosp_val * theta[["hosp_prop"]]
    
    loglik_hosp <- dnbinom(x=local_case_hosp, mu = expected_hospitalizations,size = 1/theta[["rep_var"]], log = TRUE)
  } else {
    loglik_hosp <- 0
  }
  
  
  # - - -
  # Tally up likelihoods
  # loglikSum_local_conf   
  loglikSum <- loglikSum_local_onset + loglikSum_conf #+ loglik_hosp
  exp(loglikSum) # convert to normal probability
  # w[,tt] <- exp(loglikSum)   DEBUG
}

# --------------------------------------------------------------

# Simple simulation function calc (DEPRECATED) --------------------------------------------
# --------------------------------------------------------------

simple_sim <- function(){
  
  
  # Assumptions - using daily growth rate
  ttotal <- t_period
  nn <- 100
  dt <- 1
  t_length <- ttotal/dt
  
  storeL <- array(0,dim=c(nn,t_length, length(theta_initNames)),dimnames = list(NULL,NULL,theta_initNames))
  
  # Initial condition
  storeL[,1,"inf"] <- theta[["init_cases"]]
  simzeta <- matrix(rlnorm(nn*t_length, mean = -theta[["betavol"]]^2/2, sd = theta[["betavol"]]),ncol=ttotal)
  
  for(ii in 2:ttotal){
    storeL[,ii,] <- process_model(ii-1,ii,dt,theta,storeL[,ii-1,],simzeta[,ii-1])
  }
  
  # Calculate likelihood
  log_lik <- apply(storeL[,,"cases"],1,function(x){likelihood_calc(x,case_data_matrix)})
  
  #Calculate relative probability of curves
  relative_prob <- exp(log_lik)/max(exp(log_lik))
  
  par(mfrow=c(3,1),mar=c(2,3,1,1),mgp=c(2,0.7,0))
  
  # Plot outputs
  
  plot(date_range,storeL[1,,1],col="white",ylim=c(0,1e5),xlab="",ylab="cases in Belize")
  for(ii in 1:nn){
    lines(date_range,storeL[ii,,"inf"],col=rgb(0,0,1,relative_prob[ii]))
  }
  
  # Plot daily growth rate
  plot(date_range,simzeta[1,],col="white",ylim=c(0,3),xlab="",ylab="daily growth rate")
  for(ii in 1:nn){
    lines(date_range,simzeta[ii,],col=rgb(0,0,1,relative_prob[ii]))
  }
  
  dev.copy(png,paste("plots/case_inference.png",sep=""),units="cm",width=10,height=15,res=150)
  dev.off()
  
  
}

# MLE grid search -  1D ---------------------------------------------------------

MLE_check <- function(p_name = "rep_prop", theta_tab,nn=1e3){
  
  # theta_tab <- seq(0.001,0.01,0.001)
  store_lik <- NULL
  
  for(ii in 1:length(theta_tab)){
    
    theta[[p_name]] <- theta_tab[ii]
    
    # Run SMC and output likelihood
    output_smc <- smc_model(theta,
                            nn=1e3 # number of particles
    )
    store_lik <- rbind(store_lik,c(theta_tab[ii],output_smc$lik))
    
  }
  
  colnames(store_lik) <- c("param","lik")
  store_lik <- as_tibble(store_lik)
  
}

# MLE grid search -  2D ---------------------------------------------------------

MLE_check_2D <- function(p1_name = "rep_prop", p2_name = "confirmed_prop", 
                         theta_tab1, theta_tab2,nn=1e3,
                         filename = 1){
  
  # p1_name = "rep_prop"; p2_name = "confirmed_prop"; theta_tab1 = seq(0.01,0.05,0.01); theta_tab2 = seq(0.3,1,0.1)
  
  store_lik <- NULL
  
  for(ii in 1:length(theta_tab1)){
    
    for(jj in 1:length(theta_tab2)){
      
      theta[[p1_name]] <- theta_tab1[ii]
      theta[[p2_name]] <- theta_tab2[jj]
      
      # Run SMC and output likelihooda
      output_smc <- smc_model(theta,
                              nn=1e3 # number of particles
      )
      store_lik <- rbind(store_lik,c(theta_tab1[ii],theta_tab2[jj],output_smc$lik))
      
    } 
  }
  
  colnames(store_lik) <- c("param1","param2","lik")
  store_lik <- as_tibble(store_lik)
  
  write_csv(store_lik,paste0("outputs/param_search_",filename,".csv"))
  
}

# MLE grid search -  2D ---------------------------------------------------------

MLE_check_3D <- function(p1_name = "rep_prop", 
                         p2_name = "confirmed_prop", 
                         p3_name = "betavol", 
                         theta_tab1, 
                         theta_tab2,
                         theta_tab3,
                         nn=1e3,
                         filename = 1){
  
  # p1_name = "rep_prop"; p2_name = "confirmed_prop"; p3_name = "betavol"; theta_tab1 = seq(0.01,0.05,0.02); theta_tab2 = seq(0.6,1,0.2); theta_tab3 = seq(0.1,0.3,0.1)
  
  store_lik <- NULL
  
  out_fit <- foreach(ii = 1:length(theta_tab1)) %dopar% {
    
    #for(ii in 1:length(theta_tab1)){
    
    for(jj in 1:length(theta_tab2)){
      
      for(kk in 1:length(theta_tab3)){
        
        theta[[p1_name]] <- theta_tab1[ii]
        theta[[p2_name]] <- theta_tab2[jj]
        theta[[p3_name]] <- theta_tab3[kk]
        
        # Run SMC and output likelihooda
        output_smc <- smc_model(theta,
                                nn=1e3 # number of particles
        )
        store_lik <- rbind(store_lik,c(theta_tab1[ii],theta_tab2[jj],theta_tab3[kk],output_smc$lik))
        
      }
    } 
    store_lik
  }
  
  # Collate results
  
  store_lik <- NULL
  for(ii in 1:length(theta_tab1)){
    
    store_lik <- rbind(store_lik,out_fit[[ii]])
    
  }
  
  
  colnames(store_lik) <- c("param1","param2","param3","lik")
  store_lik <- as_tibble(store_lik)
  
  write_csv(store_lik,paste0("outputs/param_search_",filename,".csv"))
  
}

# Compute acceptance probability ------------------------------------------


# Run MCMC loop ------------------------------------------