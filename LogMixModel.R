library(tidyverse)
library(data.table)
library(R2WinBUGS)
library(mcmcplots)
data <- read.csv("project2.csv", sep = " ")
data <- as.data.table(data)
data <- data[,IDgroup:=.GRP, by = subject_id]
data <- data[,':='(subject_id=NULL, ALSFRS_Total= NULL)]
data <- data[,N:=.N,by=IDgroup]
#data <- data[N>10,] #just half obs!!!
data <-  data[, ID := seq(1,5000)]
data %>% summary()
data %>% dim
str(data)
anyNA(data)

# visualization
library(ggplot2)
library(plotly)
ggplotly(ggplot(data[IDgroup %in% base::sample(1:500,40)],aes(x=time,y=bin_score,group= IDgroup,color= as.factor(IDgroup)))+geom_line())
# Logistic model (not Bayesian!)
library(lme4)

#--------------------------------------------------------
#--  just fixed parameters to check if its working   ----
#--------------------------------------------------------

#classical
results <- glm(bin_score~time,data = data, family = binomial())
results # AIC = 2759.984
plot(results) 

#Bayesian


sink("model.txt")        
cat("model FIRST
 {
  # N observations
  for (i in 1:N) {
    bin_score[i] ~ dbern(p[i])
    logit(p[i]) <- b0 + b1*time[i]
  }
    b0 ~ dnorm(0,0.01)
    b1 ~ dnorm(0,0.01)
 }", fill = TRUE)
  sink()

  bin_score <- data$bin_score 
  time <- data$time
  N = length(bin_score)
  
  dataList = list("time","bin_score","N")
  params = c("b0","b1")
  
  inits <- function(){
     list(b0 = rnorm(1,0,1),b1 = rnorm(1,0,1))
     }
  
  
  nc <- 3    #number of MCMC chains to run
  ni <- 6000  #number of samples for each chain     
  nb <- 3000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  bugs.out <-R2WinBUGS::bugs(data=dataList,inits=inits, parameters.to.save=params, model.file="model.txt", 
                   n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE,
                   bugs.directory = "F:\\WinBUGS14", working.directory=getwd())
  print(bugs.out, digits = 3)
  bugs.out$summary
  bugs.out$DIC #5750.64
  
  plot(bugs.out) # gives a summary plot of coefficients and credible intervals
  mcmcplot(bugs.out) 
  #------------------------------------------------------------
  #--  just ranodm effects, 500 group => 1000 parameters   ----
  #--          b0 and b1 normal, uncorrelated              ----
  #------------------------------------------------------------
  
  #classical
  results2 <- glm(bin_score~as.factor(IDgroup)+time*as.factor(IDgroup),data = data, family = binomial())
  results2  #AIC = 2501 
  
  #Bayesian  
  sink("model2.txt")        
  cat("model FIRST
      {
      # N observations
      for (i in 1:N) {
      bin_score[i] ~ dbern(p[i])
      logit(p[i]) <- b0[group[i]] + b1[group[i]]*time[i]
      }
      for (k in 1:K) {
      b0[k] ~ dnorm(0,0.1)
      b1[k] ~ dnorm(0,0.1)
      }

      }", fill = TRUE)
  sink()
  
  bin_score <- data$bin_score 
  time <- data$time
  group <- data$IDgroup
  K = data$IDgroup %>% unique() %>% length()
  N = length(bin_score)
  
  dataList = list("time","bin_score","group","K","N")
  params = c("b0","b1")
  
  inits <- function(){
   list(b0 = rnorm(500,0,1),b1 = rnorm(500,0,1))
  }
  
  
  nc <- 2    #number of MCMC chains to run
  ni <- 6000  #number of samples for each chain     
  nb <- 3000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  bugs.out2 <-R2WinBUGS::bugs(data=dataList,inits=inits, parameters.to.save=params, model.file="model2.txt", 
                             n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE,
                             bugs.directory = "F:\\WinBUGS14", working.directory=getwd())
  print(bugs.out2, digits = 3)
  bugs.out2$summary
  bugs.out2$DIC #2361.07
  
  plot(bugs.out2)
  mcmcplot(bugs.out2) 
  
  #------------------------------------------------------------
  #--  just ranodm effects, 500 group => 1000 parameters   ----
  #--           b0 and b1 normal, correleted               ----
  #------------------------------------------------------------
  
  sink("model3.txt")        
  cat("model FIRST
      {
      for (i in 1:N) {
      bin_score[i] ~ dbern(p[i])
      logit(p[i]) <- b[group[i],1] + b[group[i],2]*time[i]
      }
      
      for (k in 1:K) {
      b[k,1:2] ~ dmnorm(mu[],prec[,])
      }
      
      q ~ dunif(0,1)
      D[1,2] <- q*var
      D[2,1] <- q*var
      D[1,1] <- var
      D[2,2] <- var
      prec[1:2, 1:2] <- inverse(D[,])
      }", fill = TRUE)
  sink()
  
  bin_score <- data$bin_score 
  time <- data$time
  group <- data$IDgroup
  K = data$IDgroup %>% unique() %>% length()
  N = length(bin_score)
  mu = as.vector(c(0,0))
  var = 10
  dataList = list("time","bin_score","group","K","N","mu","var")
  params = c("b","q")
  
  inits=function(){list(b = matrix(rnorm(1000,0,1),ncol=2), q=runif(1,0,1) )}
  
  nc <- 2    #number of MCMC chains to run
  ni <- 6000  #number of samples for each chain     
  nb <- 3000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  bugs.out3 <-R2WinBUGS::bugs(data=dataList,inits=inits, parameters.to.save=params, model.file="model3.txt", 
                              n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE,
                              bugs.directory = "F:\\WinBUGS14", working.directory=getwd())
  print(bugs.out3, digits = 3)
  bugs.out3$summary
  bugs.out3$DIC #2364.3
  #worst than with corr=0
  
  plot(bugs.out3)
  

  #------------------------------------------------------------
  #--  just ranodm effects, 500 group => 1000 parameters   ----
  #--                    b0 and b1 IW                      ----
  #------------------------------------------------------------
  sink("model4.txt")        
  cat("model FIRST
      {
      for (i in 1:N) {
      bin_score[i] ~ dbern(p[i])
      logit(p[i]) <- b[group[i],1] + b[group[i],2]*time[i]
      }
      
      for (k in 1:K) {
      b[k,1:2] ~ dmnorm(mu[],prec[,])
      }
        
      prec[1:2, 1:2] ~ dwish(R[,], 2)
      R[1, 1] <- precision
      R[1, 2] <- 0
      R[2, 1] <- 0
      R[2, 2] <- precision

      }", fill = TRUE)
  sink()
  
  bin_score <- data$bin_score 
  time <- data$time
  group <- data$IDgroup
  K = data$IDgroup %>% unique() %>% length()
  N = length(bin_score)
  mu = as.vector(c(0,0))
  precision <- 0.0001
  dataList = list("time","bin_score","group","K","N","mu","precision")
  params = c("b")
  
  inits=function(){list(b = matrix(rnorm(1000,0,1),ncol=2))}
  
  nc <- 2    #number of MCMC chains to run
  ni <- 6000  #number of samples for each chain     
  nb <- 3000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  bugs.out4 <-R2WinBUGS::bugs(data=dataList,inits=inits, parameters.to.save=params, model.file="model4.txt", 
                              n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE,
                              bugs.directory = "F:\\WinBUGS14", working.directory=getwd())
  print(bugs.out4, digits = 3)
  bugs.out4$summary
  bugs.out4$DIC #5581.03 haven't converged
  
  plot(bugs.out4)
  
  
  #------------------------------------------------------------
  #--  just ranodm INTERCEPT, 500 group => 500 parameters  ----
  #--                      b0 normal                       ----
  #------------------------------------------------------------
  
  #classical
  results2 <- glm(bin_score~as.factor(IDgroup),data = data, family = binomial())
  results2  #AIC = 3364
  
  #Bayesian  
  sink("model5.txt")        
  cat("model FIRST
      {
      # N observations
      for (i in 1:N) {
      bin_score[i] ~ dbern(p[i])
      logit(p[i]) <- b0[group[i]]
      }
      for (k in 1:K) {
      b0[k] ~ dnorm(0,0.1)
      }
      
      }", fill = TRUE)
  sink()
  
  bin_score <- data$bin_score 
  time <- data$time
  group <- data$IDgroup
  K = data$IDgroup %>% unique() %>% length()
  N = length(bin_score)
  
  dataList = list("time","bin_score","group","K","N")
  params = c("b0")
  
  inits <- function(){
    list(b0 = rnorm(500,0,1))
  }
  
  
  nc <- 2    #number of MCMC chains to run
  ni <- 6000  #number of samples for each chain     
  nb <- 3000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  bugs.out5 <-R2WinBUGS::bugs(data=dataList,inits=inits, parameters.to.save=params, model.file="model2.txt", 
                              n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE,
                              bugs.directory = "F:\\WinBUGS14", working.directory=getwd())
  print(bugs.out5, digits = 3)
  bugs.out5$summary
  bugs.out5$DIC #2362.38 it seems that random slope is useless
  
  plot(bugs.out5)
  mcmcplot(bugs.out5) 
  
  #------------------------------------------------------------
  #--  just ranodm slope, 500 group => 500 parameters      ----
  #--                     b1 normal                        ----
  #------------------------------------------------------------
  
  #classical
  results2 <- glm(bin_score~time*as.factor(IDgroup),data = data, family = binomial())
  results2  #AIC = 2501
  
  #Bayesian  
  sink("model6.txt")        
  cat("model FIRST
      {
      # N observations
      for (i in 1:N) {
      bin_score[i] ~ dbern(p[i])
      logit(p[i]) <- b1[group[i]]*time[i]
      }
      for (k in 1:K) {
      b1[k] ~ dnorm(0,0.1)
      }
      
      }", fill = TRUE)
  sink()
  
  bin_score <- data$bin_score 
  time <- data$time
  group <- data$IDgroup
  K = data$IDgroup %>% unique() %>% length()
  N = length(bin_score)
  
  dataList = list("time","bin_score","group","K","N")
  params = c("b1")
  
  inits <- function(){
    list(b1 = rnorm(500,0,1))
  }
  
  
  nc <- 2    #number of MCMC chains to run
  ni <- 6000  #number of samples for each chain     
  nb <- 3000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  bugs.out6 <-R2WinBUGS::bugs(data=dataList,inits=inits, parameters.to.save=params, model.file="model6.txt", 
                              n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE,
                              bugs.directory = "F:\\WinBUGS14", working.directory=getwd())
  print(bugs.out6, digits = 3)
  bugs.out6$summary
  bugs.out2$DIC # 3917.56 conclusion is obvious
  
  plot(bugs.out6)
  mcmcplot(bugs.out6) 