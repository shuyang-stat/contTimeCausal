#' Continuous-time Structural Failure Time Model (ctSFTM)
#'
#' The function estimates the effect of treatment effect for a survival outcome under a SFTM
#' with time-varying treatment and confounding in the presence of dependent censoring.
#'
#' @param V the time to treatment discontinuation or failure or censoring (n x 1)
#' @param deltaV the binary indicator of treatment discontinuation at time V (n x 1)
#' @param U the time to failure or censoring (n x 1)
#' @param deltaD the binary indicator of failure at time U (n x 1)
#' @param Lti 1-dimensional baseline covariate (n x 1)
#' @param Ltd4Vtime a matrix  consisting of the time-dependent covariate (n x ltimeV)
#'
#' ltimeV is the length of uniquely observed treatment discontinuation times (called V times)
#'
#' one row represents one individual's time-dependent covariates
#'
#' columns represent ordered V times
#' @param Ltd4Utime a matrix  consisting of the time-dependent covariate (n x ltimeU)
#'
#' ltimeU is the length of uniquely observed failure times (called U times)
#'
#' one row represents one individual's time-dependent covariates
#'
#' columns represent ordered U times
#'
#' @return
#'
#' \code{est}: estimate of the SFTM parameter
#'
#'
#' @import survival MASS zoo
#'
#' @details
#' The SFTM assumes that the potential failure time \code{U} had the individual never received treatment and the observed failure time \code{T} follow
#'\deqn{U ~ \int_0^T e^{\psi A_u}d u, }
#'where \code{~} means "has the same distribution as", and \eqn{A_u} is the treatment indicator at time \eqn{u}.
#'We assume that the individual continuously received treatment until time \eqn{V}.
#' The observed failure time can be censored assuming the censoring time is independent of the failure time given the treatment and covariate history (the so-called ignorable censoring).
#'The current function provides a template to handle one-dimensional baseline covariate and one-dimensional time-dependent covariate;
#'extension to handling multiple baseline  and  time-dependent covariates is possible.
#'Variance estimate should be implemented by delete-one-group jackknifing and recalling ctSFTM.
#'
#' @seealso \code{\link{ctCoxMSM}}
#'
#'@references Yang, S., K. Pieper, and F. Cools. (2019) Semiparametric estimation of structural failure time model in continuous-time processes.
#'   \url{https://arxiv.org/abs/1808.06408}
#'
#' @examples
#'
#' library("survival")
#' library("MASS")
#' library("zoo")
#'
#'  set.seed(seed=11)
#'  n=1000
#'
#'  ## generate time-indept covariate
#'
#'  Lti<-rbinom(n,1,0.55)
#'  Lti<-Lti-mean(Lti)
#'
#'  ## generate time-dept covariate
#'
#'  Sigma<-matrix(0,3,3)
#'  for(i in 1:3){
#'    for(j in 1:3){
#'      Sigma[i,j]<-0.7^(abs(i-j))
#'    }
#'  }
#'
#'  ## Vtd represents the values of covariate at times t1=0, t2=5, and t3=10.
#'  ## We assume that the time-dependent variable remains constant between measurements.
#'
#'  Vtdtemp<-mvrnorm(n = n, rep(0,3), Sigma)
#'  Vtd<-Vtdtemp
#'  t<-c(0,5,10,100)
#'  colnames(Vtd)<-paste("t",1:3,sep="")
#'
#'  ## generate time-to-events
#'  ## D =time to death if never stop treatment (time-indep Cox)
#'  ## V =time to discontinuity (time-dep Cox)
#'  ## avoiding the same time points for V and U
#'
#'  ## generate D according to an exp distribution
#'
#'  D<-rexp(n=n,0.2)
#'
#'  Vtd<-Vtdtemp+ matrix((D-20)/5,n,3,byrow=FALSE)
#'  colnames(Vtd)<-paste("t",1:3,sep="")
#'
#'  ## generate V according to a tme-dept Cox using Bender et al (2005)
#'
#'  lambdaV <- 0.15;  betaV <- c(0.15,0.15)
#'
#'  v  <- runif(n=n)
#'  temp1 <- (- log(1-v) / (lambdaV * exp(cbind(Lti,Vtd[,1]) %*% betaV)))
#'  v  <- runif(n=n)
#'  temp2 <- (- log(1-v) / (lambdaV * exp(cbind(Lti,Vtd[,2]) %*% betaV)))
#'  v  <- runif(n=n)
#'  temp3 <- (- log(1-v) / (lambdaV * exp(cbind(Lti,Vtd[,3]) %*% betaV)))
#'  id1<-(temp1 < t[2])
#'  id2<-(temp2 < (t[3]-t[2]))
#'  id3<-(temp3 < (t[4]-t[3]))
#'  V2<- id1*temp1 + (1-id1)*id2*(temp2+t[2]) + (1-id1)*(1-id2)*(temp3+t[3])
#'
#'  ## generate Tv according to a SFTM
#'  psi<- 0
#'  true<-exp(psi)
#'
#'  id1<-D<=V2
#'  T.temp11<-D*exp(-psi[1])
#'
#'  id11<-T.temp11<=V2
#'  id12<-T.temp11>V2
#'  T.temp12<-D + V2-exp(psi[1])*V2
#'
#'  id2<-D>V2
#'  T.temp2<-D + V2-exp(psi[1])*V2
#'
#'  Tv<-id11*T.temp11+id12*T.temp12
#'
#'  ## generate censoring according to time-dept Cox
#'  ## nu=time to censoring
#'
#'  lambdaC <- 0.025; betaC <- c(0.15,0.15)
#'
#'  v  <- runif(n=n)
#'  temp3 <- (- log(1-v) / (lambdaC * exp(cbind(Lti,1) %*% betaC)))
#'  v  <- runif(n=n)
#'  temp4 <- (- log(1-v) / (lambdaC * exp(cbind(Lti,0) %*% betaC)))
#'  id3<-(temp3 < V2)
#'  nu<- id3*temp3 + (1-id3)*(V2+temp4)
#'
#'  check1<-sort( c(V2, apply(cbind(Tv,nu),1,min)))
#'  check2<-c(check1,9999)-c(0,check1)
#'
#'  if(min(check2)<10^-6){
#'    print("Please re-generate the data in order to avoid the same time points for V and U")
#'  }
#'
#'  U<-apply( cbind(Tv,nu) ,1,min)
#'  deltaD <- ( U<nu )
#'  deltaV<-(V2<U)&(V2<nu)
#'  V<-apply(cbind(V2,U,nu),1,min)
#'
#'  ## time-dependent covariate
#'  ## Ltd4Vtime is a n x ltimeV matrix consisting of the time-dependent cov
#'  ## each row represents each indiviudal
#'  ## columns represent ordered V times (the realized treatment discontinuation times)
#'
#'  data1<-list(time=V,status=deltaV)
#'  fit<-coxph(Surv(time, status) ~ . , data1)
#'  ss<-survfit(fit)
#'  obsV.times<-ss$time
#'  ltime<-length(obsV.times)
#'
#'  id1<- (obsV.times < t[2])
#'  id2<-((obsV.times < t[3])&(obsV.times > t[2]))
#'  id3<- (obsV.times > t[3])
#'  Ltd4Vtime<-matrix(NA,nrow=n,ncol=ltime)
#'  Ltd4Vtime[,which(id1==1)]<-Vtd[,1]
#'  Ltd4Vtime[,which(id2==1)]<-Vtd[,2]
#'  Ltd4Vtime[,which(id3==1)]<-Vtd[,3]
#'
#'  ## Ltd4Utime is a n x ltimeU matrix consisting of the time-dependent cov
#'  ## each row represents each indiviudal
#'  ## columns represent ordered U times (the realized event times)
#'
#'  data2<-list(time=U,status=1-deltaD)
#'  fit<-coxph(Surv(time, status) ~ . , data2)
#'  ss<-survfit(fit)
#'  obsU.times<-ss$time[ss$n.event==1]
#'  ltimeU<-length(obsU.times)
#'  id1<- (obsU.times < t[2])
#'  id2<-((obsU.times < t[3])&(obsU.times > t[2]))
#'  id3<- (obsU.times > t[3])
#'  Ltd4Utime<-matrix(NA,nrow=n,ncol=ltimeU)
#'  Ltd4Utime[,which(id1==1)]<-Vtd[,1]
#'  Ltd4Utime[,which(id2==1)]<-Vtd[,2]
#'  Ltd4Utime[,which(id3==1)]<-Vtd[,3]
#'
#'  true
#'   contTimeCausal::ctSFTM(V,deltaV,U,deltaD,Lti,Ltd4Vtime,Ltd4Utime)$est
#'
#'
#' @export
#'

ctSFTM<-function(V,deltaV,U,deltaD,Lti,Ltd4Vtime,Ltd4Utime){

  ## Fit a time-dependent Cox proportional hazards model for V:
  ## The baseline hazard function estimated by "survfit.coxph"
  ## is the hazard when all covariates are equal to their sample means.
  ## In order to evaluate the hazard function correctly, we first
  ## need to centerize the covariates by subtracting their means.

  ## This step does not not fit a time-dependent Cox PH model yet
  ## this step only prepares obsV.times for later use

  n<-length(V)
  data1<-list(time=V,status=deltaV)
  fit<-coxph(Surv(time, status) ~ . , data1)
  ss<-survfit(fit)
  obsV.times<-ss$time
  ltime<-length(obsV.times)
  cumu1.hazard<- -log(ss$surv)
  hazard.KM<-cumu1.hazard-c(0,cumu1.hazard[1:(length(cumu1.hazard)-1)])
  ss.KM<-ss$surv

  ## This step fits a time-dependent Cox PH model
  ## creat time-dependent covariate

  cnames<-c("patnum","start","stop","deltaD","deltaV","V","D.time","V.time","Lti","Vtd")
  dataforTDcoxPH<-matrix(0,(ltime)*n,length(cnames))
  colnames(dataforTDcoxPH)<-cnames
  dataforTDcoxPH[,"start"]   <-c(0,obsV.times[1:(ltime-1)])
  dataforTDcoxPH[,"stop"]    <-obsV.times[1:(ltime)]
  dataforTDcoxPH[,"patnum"]  <-rep(1:n   ,each=(ltime))
  dataforTDcoxPH[,"V"]       <-rep(V ,each=(ltime))
  dataforTDcoxPH[,"deltaD"] <-rep(deltaD,each=(ltime))
  dataforTDcoxPH[,"deltaV"] <-rep(deltaV,each=(ltime))
  dataforTDcoxPH[,"Lti"]     <-rep(Lti   ,each=(ltime))
  dataforTDcoxPH[,"D.time"]  <-(dataforTDcoxPH[,"V"]>dataforTDcoxPH[,"start"])*
    (dataforTDcoxPH[,"V"]<=dataforTDcoxPH[,"stop"])*
    (dataforTDcoxPH[,"deltaD"]==1)
  dataforTDcoxPH[,"V.time"]  <-(dataforTDcoxPH[,"V"]>dataforTDcoxPH[,"start"])*
    (dataforTDcoxPH[,"V"]<=dataforTDcoxPH[,"stop"])*
    (dataforTDcoxPH[,"deltaV"]==1)
  dataforTDcoxPH[,"Vtd"]<-as.vector(t(Ltd4Vtime))
  retain<-which( (dataforTDcoxPH[,"V"]>dataforTDcoxPH[,"start"]) )
  dataforTDcoxPH<-dataforTDcoxPH[retain,]
  dataforTDcoxPH.list<-list(start=dataforTDcoxPH[,"start"],
                            stop=dataforTDcoxPH[,"stop"],
                            V.time=dataforTDcoxPH[,"V.time"],
                            Lti=dataforTDcoxPH[,"Lti"]-mean(dataforTDcoxPH[,"Lti"]),
                            Vtd=dataforTDcoxPH[,"Vtd"]-mean(dataforTDcoxPH[,"Vtd"]))
  Lti4Ltime<-Lti
  Lti4Ltime<-Lti4Ltime-mean(dataforTDcoxPH[,"Lti"])
  Ltd4Vtime<-Ltd4Vtime-mean(dataforTDcoxPH[,"Vtd"])

  fit<- coxph(Surv(start, stop, V.time) ~Lti+Vtd,data=dataforTDcoxPH.list)
  gammahat1<-fit$coefficients
  fit$coefficient
  ss<-survfit(fit)
  cumu1.hazard<- -log(ss$surv)
  base.hazard1<-cumu1.hazard-c(0,cumu1.hazard[1:(length(cumu1.hazard)-1)])
  hazard1<-base.hazard1
  obsV.times<-ss$time
  ltime<-length(obsV.times)

  ## overall survival function for U using KM estimator

  data4<-list(time=U,status=deltaD)
  fit<-coxph(Surv(time, status) ~ . , data4)
  ss<-survfit(fit)
  obs2.times<-ss$time[ss$n.event==1]
  ltime2<-length(obs2.times)

  ############################################
  # create dMv
  # for each pat at each v time
  ############################################

  hazard2<-piTtilde<-matrix(hazard1,nrow=n,ncol=ltime,byrow=TRUE)*
    exp(Lti4Ltime*gammahat1[1]+Ltd4Vtime*gammahat1[2])
  hazard2[hazard2>1]<-1
  temp1<-1-piTtilde
  temp1[temp1<0]<-1
  Kt<-t(apply(temp1,1,cumprod))
  c4<-exp(-(Lti4Ltime*gammahat1[1]+Ltd4Vtime*gammahat1[2]))

  dM1u <- piTtilde
  matobs.times<-matrix(obsV.times,nrow=n,ncol=ltime,byrow=TRUE)
  A.mat<- matobs.times<=matrix(V,nrow=n,ncol=ltime,byrow=FALSE)
  V.mat<-matrix(V,nrow=n,ncol=ltime,byrow=FALSE)
  U.mat<-matrix(U,nrow=n,ncol=ltime,byrow=FALSE)
  dNu<-(V.mat==matobs.times)*(matrix(deltaV,nrow=n,ncol=ltime,byrow=FALSE))
  iXs1 <- (matobs.times < V.mat)|(matobs.times == V.mat) # ontrt time points indicator
  iXs2 <- (matobs.times < U.mat)|(matobs.times == U.mat) # survival time points indicator
  dMu<- ( dNu-dM1u *iXs1*iXs2)

  #############################################################
  ## Calculate IPCW
  ## using S data first
  ## then smoothing the weights for the U data
  ##
  ## exp(-integrated hazard) is approx by prod of (1-hazard);
  #############################################################

  data4<-list(time=U,status=1-deltaD)
  fit<-coxph(Surv(time, status) ~ . , data4)
  ss<-survfit(fit)
  cumu3.hazard<- -log(ss$surv)
  K3t.KM<- ss$surv[ss$n.event==1]
  obsU.times<-ss$time[ss$n.event==1]
  ltimeU<-length(obsU.times)

  cnames<-c("patnum","start","stop","deltaD","U","V","S.time","Z","Lti","Vtd")
  dataforTScoxPH<-matrix(0,(ltimeU)*n,length(cnames))
  colnames(dataforTScoxPH)<-cnames
  dataforTScoxPH[,"start"]   <-c(0,obsU.times[1:(ltimeU-1)])
  dataforTScoxPH[,"stop"]    <-obsU.times[1:(ltimeU)]
  dataforTScoxPH[,"patnum"]  <-rep(1:n   ,each=(ltimeU))
  dataforTScoxPH[,"U"]       <-rep(U ,each=(ltimeU))
  dataforTScoxPH[,"V"]       <-rep(V*deltaV+(V+1)*(1-deltaV) ,each=(ltimeU))
  dataforTScoxPH[,"Z"]       <-dataforTScoxPH[,"V"]>dataforTScoxPH[,"stop"]
  dataforTScoxPH[,"deltaD"] <-rep(deltaD,each=(ltimeU))
  dataforTScoxPH[,"Lti"]     <-rep(Lti   ,each=(ltimeU))
  dataforTScoxPH[,"S.time"]<-(dataforTScoxPH[,"U"]>dataforTScoxPH[,"start"])*
    (dataforTScoxPH[,"U"]<=dataforTScoxPH[,"stop"])*
    (dataforTScoxPH[,"deltaD"]==0)
  dataforTScoxPH[,"Vtd"]<-as.vector(t(Ltd4Utime))
  retain<-which( (dataforTScoxPH[,"U"]>dataforTScoxPH[,"start"]) )
  dataforTScoxPH<-dataforTScoxPH[retain,]
  dataforTScoxPH.list<-list(start=dataforTScoxPH[,"start"],
                            stop=dataforTScoxPH[,"stop"],
                            S.time=dataforTScoxPH[,"S.time"],
                            Lti=dataforTScoxPH[,"Lti"]-mean(dataforTScoxPH[,"Lti"]),
                            Vtd=dataforTScoxPH[,"Vtd"]-mean(dataforTScoxPH[,"Vtd"]),
                            Z=dataforTScoxPH[,"Z"])
  Lti4Utime<-Lti
  Lti4Utime<-Lti4Utime-mean(dataforTScoxPH[,"Lti"])
  Ltd4Utime<-Ltd4Utime-mean(dataforTScoxPH[,"Vtd"])

  fit<- coxph(Surv(start, stop, S.time) ~Lti+Z,
              data=dataforTScoxPH.list)
  gammahat3<-fit$coefficients

  ss<-survfit(fit)
  cumu3.hazard<- -log(ss$surv)
  base.hazard3<-cumu3.hazard-c(0,cumu3.hazard[1:(length(cumu3.hazard)-1)])
  hazard3<-base.hazard3
  if((length(hazard3)-ltimeU)>0){
    hazard3<-c(0,base.hazard3)
  }

  matobsU.times<-matrix(obsU.times,nrow=n,ncol=ltimeU,byrow=TRUE)
  Z<- matobsU.times<=matrix(V,nrow=n,ncol=ltimeU,byrow=FALSE)
  V.mat<-matrix(V,nrow=n,ncol=ltimeU,byrow=FALSE)

  obs32.times<-c(obsU.times,obs2.times)
  jjorder<-sort(obs32.times,index.return = TRUE)
  obs.times<-jjorder$x
  obs.timesid<-jjorder$ix
  augN<-length(obs.times)

  hazard3<-matrix(hazard3,nrow=n,ncol=ltimeU,byrow=TRUE)*
    exp(Lti4Utime*gammahat3[1]+Z*gammahat3[2])
  piTtilde<-hazard3

  dLambda3c <- hazard3
  U3.mat<-matrix(U,nrow=n,ncol=ltimeU,byrow=FALSE)
  dN3c<-(U3.mat==matobsU.times)*(matrix(1-deltaD,nrow=n,ncol=ltimeU,byrow=FALSE))
  iXs3 <- (matobsU.times < U3.mat)|(matobsU.times == U3.mat)
  dM3c<- ( dN3c-dLambda3c *iXs3)

  temp1<-1-piTtilde
  temp1[temp1<=0]<-1
  K3t<-t(apply(temp1,1,cumprod))#*id_piTtilde # now can be zero
  K3t.forD<-K3t.forD.KM<-matrix(NA,n,augN)
  K3t.forD[,which((jjorder$ix<(ltimeU+1)))]<-K3t
  K3t.forD.KM[,which((jjorder$ix<(ltimeU+1)))]<-K3t.KM

  eg1 <- cbind(1:augN,t(K3t.forD))
  egz1 <- zoo(eg1)
  index(egz1) <- egz1[,1]
  temp1 <- na.approx(egz1, rule = 2)
  temp1 <- temp1[,-1]
  temp1 <- as.matrix(temp1)
  K3t.forD<-t(temp1)
  K3t<-K3t.forD[,which((jjorder$ix>(ltimeU)))]

  eg1 <- cbind(1:augN,t(K3t.forD.KM))
  egz1 <- zoo(eg1)
  index(egz1) <- egz1[,1]
  temp1 <- na.approx(egz1, rule = 2)
  temp1 <- temp1[,-1]
  temp1 <- as.matrix(temp1)
  K3t.forD.KM<-t(temp1)
  K3t.KM<-K3t.forD.KM[,which((jjorder$ix>(ltimeU)))]
  K3s<-K3t.forD[,which((jjorder$ix<=(ltimeU)))]

  matobs2.times<-matrix(obs2.times,nrow=n,ncol=ltime2,byrow=TRUE)
  matK3t<-K3t
  matK3t.KM<-K3t.KM

  U2.mat<-matrix(U,nrow=n,ncol=ltime2,byrow=FALSE) # U time mat 2
  Delta2.mat<-matrix(deltaD,nrow=n,ncol=ltime2,byrow=FALSE) # U time mat 2
  K3   <-apply( (U2.mat==matobs2.times)*Delta2.mat*matK3t,1,sum)
  K3.KM<-apply( (U2.mat==matobs2.times)*Delta2.mat*matK3t.KM,1,sum)
  IPCW<-rep(0,n)
  IPCW[deltaD>0]<-1/K3[deltaD>0]

  part1<-U*(1-deltaV)+V*(deltaV)
  part2<-(U-V)*(deltaV)

  lm.out<-lm(part1[IPCW>0]~Lti[IPCW>0]+Vtd[IPCW>0,1],weights=IPCW[IPCW>0])
  coeff<-lm.out$coefficients
  Epart1<-coeff[1]+coeff[2]*Lti+coeff[3]*Vtd[,1]

  lm.out<-lm(part2[IPCW>0]~Lti[IPCW>0]+Vtd[IPCW>0,1],weights=IPCW[IPCW>0])
  coeff<-lm.out$coefficients
  Epart2<-coeff[1]+coeff[2]*Lti+coeff[3]*Vtd[,1]

  x0psi.mat1<-matrix(part1,nrow=n,ncol=ltime,byrow=FALSE)
  x0psi.mat2<-matrix(part2,nrow=n,ncol=ltime,byrow=FALSE)
  x0psi.mat1c<-matrix(part1-Epart1,nrow=n,ncol=ltime,byrow=FALSE)
  x0psi.mat2c<-matrix(part2-Epart2,nrow=n,ncol=ltime,byrow=FALSE)

  sum41c<-sum(c4*x0psi.mat1c*dMu*iXs1*iXs2*IPCW)
  sum42c<-sum(c4*x0psi.mat2c*dMu*iXs1*iXs2*IPCW)
  est<- -sum42c/sum41c

  est

  return(list(est=est))
}
