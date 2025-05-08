# Helper functions --------------------------------------------

c.text<-function(x,sigF=3){
  bp1=signif(c(median(x),quantile(x,0.025),quantile(x,0.975,na.rm = TRUE)),sigF)
  paste(bp1[1]," (",bp1[2],"-",bp1[3],")",sep="")
}

c.nume<-function(x){
  bp1=c(median(x),quantile(x,0.025),quantile(x,0.975,na.rm = TRUE))
  as.numeric(bp1)
}

bin_conf <- function(x,n){
  htest <- binom.test(x,n)
  h_out <- c(x/n, htest$conf.int[1], htest$conf.int[2])
  h_out
}

# Run SMC to get bootstrap estimates --------------------------------------------

run_fits <- function(rep_plot,nn,cut_off,dt,filename="1_vacc"){
  
  # rep_plot <- 100; nn <- 2e3; cut_off=0; dt=0.25
  
  out_rep <- foreach(kk = 1:rep_plot, .options.RNG = 123) %dopar% {
    output_smc <- smc_model(theta,nn,dt,seed=123+kk)
    output_smc
  }
  
  S_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  I_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  H_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  R_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  C_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  Rep_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  R0_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  Rec_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  CC_plot = matrix(NA,ncol=rep_plot,nrow=t_period)
  
  for(kk in 1:rep_plot){
    output_smc <- out_rep[[kk]]
    if(output_smc$lik != - Inf){
      I_plot[,kk] <- output_smc$E_trace + output_smc$I_trace*(1-theta[["confirmed_prop"]])
      H_plot[,kk] <- output_smc$H_trace
      R_plot[,kk] <- output_smc$R_trace
      S_plot[,kk] <- output_smc$S_trace
      case_pos <- theta[["confirmed_prop"]]*theta[["rep_prop"]]*(output_smc$C_trace - c(0,head(output_smc$C_trace,-1)))
      C_plot[,kk] <- rpois(length(case_pos),lambda=case_pos)
      rep_pos <- theta[["confirmed_prop"]]*(output_smc$Rep_trace - c(0,head(output_smc$Rep_trace,-1))) # ALL CASES: theta[["rep_prop"]]*
      Rep_plot[,kk] <- rpois(length(rep_pos),lambda=rep_pos)
      R0_plot[,kk] <- output_smc$beta_trace/(theta[["recover"]])
      recov_case <- theta[["confirmed_prop"]]*theta[["rep_prop"]]*(output_smc$R_trace - c(0,head(output_smc$R_trace,-1)))
      Rec_plot[,kk] <- rpois(length(recov_case), lambda=recov_case)
      cum_case_pos <- theta[["confirmed_prop"]]*theta[["rep_prop"]]*(output_smc$C_trace)
      CC_plot[,kk] <- rpois(length(cum_case_pos),lambda=cum_case_pos)
    }
  }
  
  save(
    S_plot,
    I_plot,
    H_plot,
    R_plot,
    C_plot,
    Rep_plot,
    R0_plot,
    Rec_plot,
    CC_plot,
    file=paste0("Vacc/outputs_vacc/bootstrap_fit_",filename,".RData")) 
  
}

# Plot outputs from SMC --------------------------------------------


plot_outputs <- function(filename="1_vacc"){
  
  # filename="1"
  
  cut_off <- 0 #end_date - as.Date("2020-01-23")
  
  load(paste0("Vacc/outputs_vacc/bootstrap_fit_",filename,".RData"))
  
  # Remove NA fits
  S_plot = S_plot[,!is.na(S_plot[t_period,])]
  I_plot = I_plot[,!is.na(I_plot[t_period,])]
  H_plot = H_plot[,!is.na(H_plot[t_period,])]
  R_plot = R_plot[,!is.na(R_plot[t_period,])]
  C_plot = C_plot[,!is.na(C_plot[t_period,])]
  Rep_plot = Rep_plot[,!is.na(Rep_plot[t_period,])]
  R0_plot = R0_plot[,!is.na(R0_plot[t_period,])]
  Rec_plot = Rec_plot[,!is.na(Rec_plot[t_period,])]
  CC_plot = CC_plot[,!is.na(CC_plot[t_period,])]
  
  # Calculate quantiles
  S_quantile <- apply(S_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)}) # number
  S_quantile_n_average <- apply(S_plot,1,function(x){mean(x,na.rm=TRUE)}) # number
  Inf_quantile <- apply(C_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)}) # number
  Inf_quantile_n_average <- apply(theta[["onset_prop"]]*C_plot,1,function(x){mean(x,na.rm=TRUE)}) # number
  Inf_quantile_cum <- cumsum(Inf_quantile_n_average)
  H_quantile <- apply(H_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)})
  H_quantile_n_average <- apply(H_plot,1,function(x){mean(x,na.rm=TRUE)})
  H_quantile_cum <- cumsum(H_quantile_n_average)
  R_quantile <- apply(Rec_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)})
  R_quantile_n_average <- apply(Rec_plot,1,function(x){mean(x,na.rm=TRUE)})
  R_quantile_cum <- cumsum(R_quantile_n_average)
  
  R0_quantile <- apply(R0_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)}) #*S_plot/theta[["pop"]]
  R0_average <- apply(R0_plot,1,mean)
  
  Case_quantile <- apply(C_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)}) 
  Case_quantile_onset <- theta[["onset_prop"]]*Case_quantile
  Case_average_onset <- theta[["onset_prop"]]*apply(C_plot,1,mean)
  
  CCase_quantile <- apply(CC_plot,1,function(x){quantile(x,c(0.025,0.25,0.5,0.75,0.975),na.rm = TRUE)}) 
  CCase_quantile_onset <- theta[["onset_prop"]]*CCase_quantile
  CCase_average_onset <- theta[["onset_prop"]]*apply(CC_plot/10,1,mean)
  
  
  # Remove final few points (as estimation less reliable)
  S_quantileA <- S_quantile[,1:(ncol(R0_quantile)-cut_off)]
  Inf_quantileA <- Inf_quantile[,1:(ncol(R0_quantile)-cut_off)]
  H_quantileA <- H_quantile[,1:(ncol(R0_quantile)-cut_off)]
  R0_quantileA <- R0_quantile[,1:(ncol(R0_quantile)-cut_off)]
  date_rangeA <- date_range[1:(length(date_range)-cut_off)]
  Case_quantile_onsetA <- Case_quantile_onset[,1:(length(date_range)-cut_off)]
  CCase_quantile <- CCase_quantile[,1:(length(date_range)-cut_off)]
  
  save(
    R0_average,
    R0_quantile,
    Case_quantile,
    Case_quantile_onset,
    Case_average_onset,
    CCase_quantile,
    CCase_quantile_onset,
    CCase_average_onset,
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
    file=paste0("Vacc/outputs_vacc/results_",filename,".RData")) 
  
  
  
  # forecast window
  date_rangeF <- date_range[date_range>max(case_data_bl$date)]
  yyF <- rep(0,length(date_rangeF))
  
  # - - - - - - - 
  # Calculate daily incidence
  #Case_diff_quantile <- Case_quantile[,1:ncol(Case_quantile)] - cbind(c(0,0,0),Case_quantile[,1:(ncol(Case_quantile)-1)])
  
  #par(mar=c(2,3,1,1),mgp=c(2,0.55,0)) #mfrow=c(4,2),
  #layout(matrix(c(1,1,1,2,3,4,5,6,7), 3, 3, byrow = TRUE))
  
  # Plot outputs
  a_col <- 0.4 # alpha
  xMin1 <- as.Date("2020-03-22") #min(as.Date("2019-12-01")) 
  xMin <- xMin1 #min(as.Date("2020-01-01"))
  xMax <- end_date-1 #max(date_range)
  
  
  # Plot reproduction number
  xMaxR <- as.Date("2023-6-30")
  date_rangeB <- date_rangeA#[date_rangeA>as.Date("2020-03-03")]
  R0_quantileB <- R0_quantileA#[,date_rangeA>as.Date("2020-03-03")]
  xMax1 <- xMax #as.Date("2020-03-31") #xMax #xMax #
  
  plot(date_rangeB,R0_quantileB[1,],col="white",ylim=c(0,15),xlim=c(xMin1,xMaxR),xlab="",ylab=expression(paste(R[t])))
  polygon(c(date_rangeF,rev(date_rangeF)),c(yyF,rev(yyF+1e5)),lty=0,col=rgb(0.9,0.9,0.9))
  
  polygon(c(date_rangeB,rev(date_rangeB)),c(R0_quantileB[2,],rev(R0_quantileB[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeB,rev(date_rangeB)),c(R0_quantileB[1,],rev(R0_quantileB[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_rangeB,R0_quantileB[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  lines(date_rangeB,1+0*R0_quantileB[3,],lty=2)
  
  #text(labels="model",x=xMin1,y=9,adj=0,col="blue")
  
  lines(c(Bl_vacc,Bl_vacc),c(0,20),col="red")
  
  title(LETTERS[1],adj=0); letR = 2
  
  dev.copy(png,paste("Vacc/plots/reproduction_number.png",sep=""),units="cm",width=10,height=15,res=150)
  dev.off()
  
  
  # Plot local case onsets
  ym1 <- 1500
  plot(date_rangeA,Case_quantile_onsetA[1,],col="white",ylim=c(0,ym1),xlim=c(xMin1,xMax),xlab="",ylab="new onsets in Belize")

  polygon(c(date_rangeF,rev(date_rangeF)),c(yyF,rev(yyF+1e5)),lty=0,col=rgb(0.9,0.9,0.9))
  
  polygon(c(date_rangeA,rev(date_rangeA)),c(Case_quantile_onsetA[2,],rev(Case_quantile_onsetA[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeA,rev(date_rangeA)),c(Case_quantile_onsetA[1,],rev(Case_quantile_onsetA[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  points(daily_cases$date,daily_cases$new_cases,pch=1,cex=1, col="darkgrey")
  #lines(date_rangeA,Case_quantile_onsetA[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  lines(date_rangeA,Case_average_onset,type="l",col="darkgreen",xaxt="n",yaxt="n",xlab="",ylab="",lty=1,lwd=3)
  lines(date_rangeA,H_quantile_n_average,type="l",col="blue",xaxt="n",yaxt="n",xlab="",ylab="",lty=1,lwd=3)
  
  #text(labels="vaccination",x=bl_vacc+0.5,y=0.9*ym1,adj=0,col="red")
  
  #text(labels="model (median)",x=xMin1,y=0.9*ym1,adj=0,col="blue")
  text(labels="daily cases model (mean)",x=xMin1,y=0.8*ym1,adj=0,col="darkgreen")
  text(labels="Hospitalisation",x=xMin1,y=0.6*ym1,adj=0,col="blue")
  text(labels="data daily cases",x=xMin1,y=0.7*ym1,adj=0,col="black")
  
 
  #title(LETTERS[letR],adj=0); letR = letR + 1
  
  dev.copy(png,paste("Vacc/plots/new_cases.png",sep=""),units="cm",width=10,height=15,res=150)
  dev.off()
  
  
  # - - -
  
  # Plot susceptibles
  # 
  plot(date_rangeA,S_quantileA[1,],col="white",ylim=c(0,1),xlim=c(xMin1,xMax),xlab="",ylab="propn prevalence in Belize")
  polygon(c(date_rangeA,rev(date_rangeA)),c(S_quantileA[2,],rev(S_quantileA[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_rangeA,rev(date_rangeA)),c(S_quantileA[1,],rev(S_quantileA[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_rangeA,S_quantileA[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  #   
  text(labels="model estimate",x=min(date_rangeA),y=9,adj=0,col="blue")
  #   
  lines(c(Bl_vacc,Bl_vacc),c(0,10),col="red")
  title(LETTERS[letR],adj=0); letR = letR + 1
  
  dev.copy(png,paste("Vacc/plots/susceptibles.png",sep=""),units="cm",width=10,height=15,res=150)
  dev.off()
  

  # Plot case total predictions
  ym1 <- 100000
  plot(date_range,CCase_quantile_onset[1,],col="white",ylim=c(0,ym1),xlim=c(xMin1,xMax),xlab="",ylab="cases in Belize")
  
  polygon(c(date_rangeF,rev(date_rangeF)),c(yyF,rev(yyF+1e5)),lty=0,col=rgb(0.9,0.9,0.9))
  
  polygon(c(date_range,rev(date_range)),c(CCase_quantile_onset[2,],rev(CCase_quantile_onset[4,])),lty=0,col=rgb(0,0.3,1,0.35))
  polygon(c(date_range,rev(date_range)),c(CCase_quantile_onset[1,],rev(CCase_quantile_onset[5,])),lty=0,col=rgb(0,0.3,1,0.2))
  lines(date_range,CCase_quantile_onset[3,],type="l",col=rgb(0,0,1),xaxt="n",yaxt="n",xlab="",ylab="")
  #lines(date_range,CCase_average_onset,type="l",col="green",xaxt="n",yaxt="n",xlab="",ylab="")
  lines(c(Bl_vacc,Bl_vacc),c(0,1e6),col="red")

  text(labels="model (median)",x=xMin1,y=0.9*ym1,adj=0,col="blue")
  #text(labels="model (mean)",x=xMin1,y=0.8*ym1,adj=0,col="blue")
  text(labels="data",x=xMin1,y=0.7*ym1,adj=0,col="black")
  
  title(LETTERS[letR],adj=0); letR = letR + 1
  
  points(daily_cases$date,daily_cases$total_cases,pch=18,cex=1.1)
  
  dev.copy(png,paste("Vacc/plots/total_infections.png",sep=""),units="cm",width=10,height=15,res=150)
  dev.off()
  
  
  # output figure
  #dev.copy(png,paste("plots/cases_inference.png",sep=""),units="cm",width=18,height=18,res=150)
  #dev.copy(pdf,paste("plots/Figure_2.pdf",sep=""),width=7.5,height=8)
  #dev.off()
  
  
}



# Output R0 estimates over time --------------------------------------------

r0_value_output <- function(filename="1"){
  
  # filename="1"
  
  load(paste0("bootstrap_fit_",filename,".RData"))
  
  # Extract R0 estimates
  
  period_interest <- as.Date("2021-03-01")
  
  R0_plot = R0_plot[,!is.na(R0_plot[t_period,])]
  
  med_R0 <- apply(R0_plot,1,function(x){quantile(x,c(0.5),na.rm = TRUE)})
  c.text(med_R0[match(period_interest,date_range)])
  
  #file_out <- as_tibble(cbind(date_range,c.text(t(med_R0))))
  
  # Median R0 after before closure
  
  period_interest <- as.Date(c("2020-03-22","2023-06-30"))
  index_pick <- match(period_interest,date_range)
  R0_all <- R0_plot[index_pick[1]:index_pick[2],]
  
  med_R03 <- apply(R0_all,1,function(x){quantile(x,c(0.5),na.rm = TRUE)})
  
  # Median R0 range before closure
  period_interest <- as.Date(c("2020-03-22","2021-03-30"))
  index_pick <- match(period_interest,date_range)
  R0_all <- R0_plot[index_pick[1]:index_pick[2],]
  
  med_R02 <- apply(R0_all,1,function(x){quantile(x,c(0.5),na.rm = TRUE)})
  med_R02 <- R0_all; dim(med_R02) <- NULL
  
  
  # Median R0 after before closure
  
  period_interest <- as.Date(c("2020-03-22","2021-02-28"))
  index_pick <- match(period_interest,date_range)
  R0_all <- R0_plot[index_pick[1]:index_pick[2],]
  
  #med_R0 <- apply(R0_all,1,function(x){quantile(x,c(0.5))})
  med_R0 <- R0_all; dim(med_R0) <- NULL
  
  out_r <- cbind(c.text(med_R02),c.text(med_R0),c.text(med_R03))
  out_r <- as_tibble(out_r); names(out_r) <- c("before","after")
  
  write_csv(out_r,"Vacc/outputs_vacc/before_after_R.csv")
  
  # - - - 
  # symptomatic cases in Belize
  
  period_interest <- as.Date("2023-06-30")
  
  c.text(I_plot[date_range==period_interest,])
  
  
}


# Plot distributions ------------------------------------------------------

plot_distn <- function(){
  
  xx <- seq(0,20,0.1)
  yy_recover <- dgamma(xx,shape=2,rate=2/(1/theta[["recover"]]))
  yy_incubation <- dgamma(xx,shape=2,rate=2/(1/theta[["incubation"]]))
  yy_report <- dexp(xx,rate=theta[["report"]])
  
  par(mfrow=c(1,3),mar=c(3,3,1,1),mgp=c(2,0.7,0)) #mfrow=c(4,2),
  plot(xx,yy_incubation,xlab="incubation period",ylab="probability density",type="l",yaxs="i",lwd=2)
  title(LETTERS[1],adj=0)
  plot(xx,yy_recover,xlab="infectious period",ylab="probability density",type="l",yaxs="i",lwd=2)
  title(LETTERS[2],adj=0)
  plot(xx,yy_report,xlab="delay onset-to-confirmation",ylab="probability density",type="l",yaxs="i",lwd=2)
  title(LETTERS[3],adj=0)
  
  dev.copy(pdf,paste("Vacc/plots/Figure_S1.pdf",sep=""),width=8,height=3)
  dev.off()
  
}



