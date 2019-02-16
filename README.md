---
  output: github_document
---

<!-- rmarkdown v1 -->

<!-- README.md is generated from README.Rmd. Please edit that file -->




# Continuous Time Causal Models

![equation](http://www.sciweavers.org/tex2img.php?eq=1%2Bsin%28mc%5E2%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)

### Continuous-time causal models for evaluating the effect of time-varying treatments  

### Visit the [package website](https://github.com/shuyang1987/contTimeCausal)

## Installation with `devtools`:


```r
devtools::install_github("shuyang1987/contTimeCausal")
```

## Description




Structural failure time models (SFTM) and Cox marginal structural models (MSM) are causal models for estimating the effect of time-varying treatments on a survival outcome. 
This package provides a class of continuous-time structural failure time models (ctSFTM) and continuous-time Cox marginal structural models (ctCoxMSM), which respects the continuous time nature of the underlying data processes in practice (i.e. the variables and processes are measured at irregularly spaced time points, which are not necessarily the same for all subjects). 

The ctSFTM assumes that the potential failure time $U$ had the individual never received treatment and the observed failure time $T$ follow
$$U \thicksim \int_0^T e^{\psi\times A_u}d u, $$
where $\thicksim$ means "has the same distribution as", and $A_u$ is the treatment indicator at time $u$.

The ctCoxMSM assumes that the potential failure time 
$T^{\overline{a}}$ under the treatment regime $\overline{a}$ (a hypothetical treatment regime from baseline to the event time) follows a proportional hazards model with the hazard function at $u$ as $$\lambda_0(u)e^{\psi\times a_u}.$$

We assume that the individual continuously received treatment until time $V$. The observed failure time can be censored. We assume an ignorable censoring mechanism in the sense that the censoring time is independent of the failure time given the treatment and covariate history. The stSFTM() returns [a continuous-time g-estimator](https://arxiv.org/abs/1808.06408)
of $\psi$ and the ctCoxMSM() returns [an inverse probability of treatment weighting estimator](https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.12845) of $\psi$.

The current function provides a template to handle one-dimensional baseline covariate and one-dimensional time-dependent covariate; extension to handling multiple baseline  and  time-dependent covariates is possible. Variance estimation should be implemented by delete-one-group jackknifing and recalling ctSFTM or ctCoxMSM.



### Main Papers: Yang et al. (2018) and Yang et al. (2019)

Yang, S., K. Pieper, and F. Cools. (2019) Semiparametric estimation of structural failure time model in continuous-time processes. https://arxiv.org/abs/1808.06408

Yang, S., A. A. Tsiatis, and M. Blazing (2018). Modeling survival distribution as a function of time to treatment discontinuation: a dynamic treatment regime approach, *Biometrics*,
 72, 1055--1065. https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.12845


## Use

ctSFTM()   estimates the effect of treatment effect for a survival outcome under a ctSFTM with time-varying treatment and confounding in the presence of dependent censoring.

ctCoxMSM() estimates the effect of treatment effect for a survival outcome under a ctCoxMSM with time-varying treatment and confounding in the presence of dependent censoring.


### Toy example


```r
library("survival")
library("MASS")
library("zoo")
#> 
#> Attaching package: 'zoo'
#> 
#> The following objects are masked from 'package:base':
#> 
#>     as.Date, as.Date.numeric

 set.seed(seed=11)
 n=1000

 ## generate time-indept covariate

 Lti<-rbinom(n,1,0.55)
 Lti<-Lti-mean(Lti)

 ## generate time-dept covariate

 Sigma<-matrix(0,3,3)
 for(i in 1:3){
   for(j in 1:3){
     Sigma[i,j]<-0.7^(abs(i-j))
   }
 }

 ## Vtd represents the values of covariate at times t1=0, t2=5, and t3=10.
 ## We assume that the time-dependent variable remains constant between measurements.

 Vtdtemp<-mvrnorm(n = n, rep(0,3), Sigma)
 Vtd<-Vtdtemp
 t<-c(0,5,10,100)
 colnames(Vtd)<-paste("t",1:3,sep="")

 ## generate time-to-events
 ## D =time to death if never stop treatment (time-indep Cox)
 ## V =time to discontinuity (time-dep Cox)
 ## avoiding the same time points for V and U

 ## generate D according to an exp distribution

 D<-rexp(n=n,0.2)

 Vtd<-Vtdtemp+ matrix((D-20)/5,n,3,byrow=FALSE)
 colnames(Vtd)<-paste("t",1:3,sep="")

 ## generate V according to a tme-dept Cox using Bender et al (2005)

 lambdaV <- 0.15;  betaV <- c(0.15,0.15)

 v  <- runif(n=n)
 temp1 <- (- log(1-v) / (lambdaV * exp(cbind(Lti,Vtd[,1]) %*% betaV)))
 v  <- runif(n=n)
 temp2 <- (- log(1-v) / (lambdaV * exp(cbind(Lti,Vtd[,2]) %*% betaV)))
 v  <- runif(n=n)
 temp3 <- (- log(1-v) / (lambdaV * exp(cbind(Lti,Vtd[,3]) %*% betaV)))
 id1<-(temp1 < t[2])
 id2<-(temp2 < (t[3]-t[2]))
 id3<-(temp3 < (t[4]-t[3]))
 V2<- id1*temp1 + (1-id1)*id2*(temp2+t[2]) + (1-id1)*(1-id2)*(temp3+t[3])

 ## generate Tv according to a SFTM
 psi<- 0
 true<-exp(psi)

 id1<-D<=V2
 T.temp11<-D*exp(-psi[1])

 id11<-T.temp11<=V2
 id12<-T.temp11>V2
 T.temp12<-D + V2-exp(psi[1])*V2

 id2<-D>V2
 T.temp2<-D + V2-exp(psi[1])*V2

 Tv<-id11*T.temp11+id12*T.temp12

 ## generate censoring according to time-dept Cox
 ## nu=time to censoring

 lambdaC <- 0.025; betaC <- c(0.15,0.15)

 v  <- runif(n=n)
 temp3 <- (- log(1-v) / (lambdaC * exp(cbind(Lti,1) %*% betaC)))
 v  <- runif(n=n)
 temp4 <- (- log(1-v) / (lambdaC * exp(cbind(Lti,0) %*% betaC)))
 id3<-(temp3 < V2)
 nu<- id3*temp3 + (1-id3)*(V2+temp4)

 check1<-sort( c(V2, apply(cbind(Tv,nu),1,min)))
 check2<-c(check1,9999)-c(0,check1)

 if(min(check2)<10^-6){
   print("Please re-generate the data in order to avoid the same time points for V and U")
 }
#> [1] "Please re-generate the data in order to avoid the same time points for V and U"

 U<-apply( cbind(Tv,nu) ,1,min)
 deltaD <- ( U<nu )
 deltaV<-(V2<U)&(V2<nu)
 V<-apply(cbind(V2,U,nu),1,min)

 ## time-dependent covariate
 ## Ltd4Vtime is a n x ltimeV matrix consisting of the time-dependent cov
 ## each row represents each indiviudal
 ## columns represent ordered V times (the realized treatment discontinuation times)

 data1<-list(time=V,status=deltaV)
 fit<-coxph(Surv(time, status) ~ . , data1)
 ss<-survfit(fit)
 obsV.times<-ss$time
 ltime<-length(obsV.times)

 id1<- (obsV.times < t[2])
 id2<-((obsV.times < t[3])&(obsV.times > t[2]))
 id3<- (obsV.times > t[3])
 Ltd4Vtime<-matrix(NA,nrow=n,ncol=ltime)
 Ltd4Vtime[,which(id1==1)]<-Vtd[,1]
 Ltd4Vtime[,which(id2==1)]<-Vtd[,2]
 Ltd4Vtime[,which(id3==1)]<-Vtd[,3]

 ## Ltd4Utime is a n x ltimeU matrix consisting of the time-dependent cov
 ## each row represents each indiviudal
 ## columns represent ordered U times (the realized event times)

 data2<-list(time=U,status=1-deltaD)
 fit<-coxph(Surv(time, status) ~ . , data2)
 ss<-survfit(fit)
 obsU.times<-ss$time[ss$n.event==1]
 ltimeU<-length(obsU.times)
 id1<- (obsU.times < t[2])
 id2<-((obsU.times < t[3])&(obsU.times > t[2]))
 id3<- (obsU.times > t[3])
 Ltd4Utime<-matrix(NA,nrow=n,ncol=ltimeU)
 Ltd4Utime[,which(id1==1)]<-Vtd[,1]
 Ltd4Utime[,which(id2==1)]<-Vtd[,2]
 Ltd4Utime[,which(id3==1)]<-Vtd[,3]

 true
#> [1] 1
 contTimeCausal::ctSFTM(V,deltaV,U,deltaD,Lti,Ltd4Vtime,Ltd4Utime)$est
#> [1] 0.9923524
 contTimeCausal::ctCoxMSM(V,deltaV,U,deltaD,Lti,Ltd4Vtime,Ltd4Utime)$est
#> [1] 1.001188
```


