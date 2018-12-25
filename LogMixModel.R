library(tidyverse)
library(data.table)
library(R2WinBUGS)
library(mcmcplots)
data <- read.csv("project2.csv", sep = " ")
data <- as.data.table(data)
data <- data[,IDgroup:=.GRP, by = subject_id]
data <- data[,':='(subject_id=NULL, ALSFRS_Total= NULL)]
data <- data[,N:=.N,by=IDgroup]
#data[N<4,] %>% dim # 77 patients have less than 4 obs
#data[,sum:=sum(bin_score==0),by=IDgroup][sum==0,]$IDgroup %>% unique %>% length 
# 59 patients had never abnormal results
#data <- data[N>10,] #maybe we sould take just patients with more than X obs?
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
  
  #-------------------------------------|
  #   |--| |--| |--| |----|     |       |
  #   |--| |--| |--|   ||       |       | 
  #   |    |  | |  \   ||       |       |               
  #-------------------------------------|
  
  
  
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
  
  dataList = list("bin_score","group","K","N")
  params = c("b0")
  
  inits <- function(){
    list(b0 = rnorm(500,0,1))
  }
  
  
  nc <- 2    #number of MCMC chains to run
  ni <- 6000  #number of samples for each chain     
  nb <- 3000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  bugs.out5 <-R2WinBUGS::bugs(data=dataList,inits=inits, parameters.to.save=params, model.file="model5.txt", 
                              n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE,
                              bugs.directory = "F:\\WinBUGS14", working.directory=getwd())
  print(bugs.out5, digits = 3)
  bugs.out5$summary
  bugs.out5$DIC # 3178.04, random slope is important
  
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
  
  
  #----------------------------------------|
  #   |--| |--| |--| |----|    |  |        |
  #   |--| |--| |--|   ||      |  |        | 
  #   |    |  | |  \   ||      |  |        |               
  #----------------------------------------|
  
names(data)
# we can add AGE, GENDER, TIME and TREATMENT
data <- data[,Age := ((Age-min(Age))/(max(Age)-min(Age)))]

#-------------------------------------------------------------------------
#--  ranodm intercept and slope + fixed age,gender, time and treatment----
#--                         all priors normal                         ----
#-------------------------------------------------------------------------
#classical
results7 <- glm(bin_score~as.factor(IDgroup)+time*as.factor(IDgroup)+time+Age+treatment+gender,data = data, family = binomial())
results7 #haven't converge

#Bayesian  
sink("model7.txt")        
cat("model FIRST
    {
    # N observations
    for (i in 1:N) {
    bin_score[i] ~ dbern(p[i])
    logit(p[i]) <- b0[group[i]] + b1[group[i]]*time[i] + 
      a1*age[i] + a2*treatment[i] + a3*gender[i] + a4*time[i]
    }
    for (k in 1:K) {
    b0[k] ~ dnorm(0,0.1)
    b1[k] ~ dnorm(0,0.1)
    }
    
    a1 ~ dnorm(0,0.001)
    a2 ~ dnorm(0,0.001)
    a3 ~ dnorm(0,0.001)
    a4 ~ dnorm(0,0.001)
    
    }", fill = TRUE)
  sink()
  
  bin_score <- data$bin_score 
  time <- data$time
  group <- data$IDgroup
  
  K = data$IDgroup %>% unique() %>% length()
  N = length(bin_score)
  
  age <- data$Age
  treatment <- data$treatment 
  gender <- data$gender
  
  dataList = list("time","bin_score","group","K","N","age","treatment","gender")
  params = c("b0","b1","a1","a2","a3","a4")
  
  inits7 <- function(){
    list(b0 = rnorm(500,0,1),b1 = rnorm(500,0,1),
         a1 = rnorm(1,0,1),a2 = rnorm(1,0,1),a3 = rnorm(1,0,1),a4 = rnorm(1,0,1))
  }
  
  
  nc <- 2    #number of MCMC chains to run
  ni <- 6000  #number of samples for each chain     
  nb <- 3000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  bugs.out7 <-R2WinBUGS::bugs(data=dataList,inits=inits7, parameters.to.save=params, model.file="model7.txt", 
                              n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE,
                              bugs.directory = "F:\\WinBUGS14", working.directory=getwd())
  print(bugs.out7, digits = 3)
  bugs.out7$summary
  bugs.out7$DIC #2007.74
  
  plot(bugs.out7)
  
  #-------------------------------------------------------------------------
  #--            ranodm intercept and slope + fixed time                ----
  #--                          b0,b1,a4 normal                          ----
  #-------------------------------------------------------------------------
  
  #Bayesian  
  sink("model8.txt")        
  cat("model FIRST
      {
      # N observations
      for (i in 1:N) {
      bin_score[i] ~ dbern(p[i])
      logit(p[i]) <- b0[group[i]] + b1[group[i]]*time[i] + a4*time[i]
      }
      for (k in 1:K) {
      b0[k] ~ dnorm(0,0.1)
      b1[k] ~ dnorm(0,0.1)
      }
      
      a4 ~ dnorm(0,0.001)
      
      }", fill = TRUE)
  sink()
  
  bin_score <- data$bin_score 
  time <- data$time
  group <- data$IDgroup
  
  K = data$IDgroup %>% unique() %>% length()
  N = length(bin_score)
  
  dataList = list("time","bin_score","group","K","N")
  params = c("b0","b1","a4")
  
  inits8 <- function(){
    list(b0 = rnorm(500,0,1),b1 = rnorm(500,0,1),a4 = rnorm(1,0,1))
  }
  
  
  nc <- 2    #number of MCMC chains to run
  ni <- 6000  #number of samples for each chain     
  nb <- 3000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  bugs.out8 <-R2WinBUGS::bugs(data=dataList,inits=inits8, parameters.to.save=params, model.file="model8.txt", 
                              n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE,
                              bugs.directory = "F:\\WinBUGS14", working.directory=getwd())
  print(bugs.out8, digits = 3)
  bugs.out8$summary
  bugs.out8$DIC # 2013.99  - for me it seems to be the best model
  
  plot(bugs.out8)
  
  #posterior predictive checks to evaluate the chosen model 
  
  #no idea how to do it
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  #verify whether the probit or the complementary log-log link provide a better model.
  
  
  #-------------------------------------------------------------------------
  #--            ranodm intercept and slope + fixed time                ----
  #--                          b0,b1,a4 normal                          ----
  #--                       probit link function                        ----
  #-------------------------------------------------------------------------
 
  sink("model9.txt")        
  cat("model FIRST
      {
      # N observations
      for (i in 1:N) {
      bin_score[i] ~ dbern(p[i])
      p[i] <- phi(b0[group[i]] + b1[group[i]]*time[i] + a4*time[i])
      }
      for (k in 1:K) {
      b0[k] ~ dnorm(0,1)
      b1[k] ~ dnorm(0,1)
      }
      
      a4 ~ dnorm(0,1)
      
      }", fill = TRUE)
  sink()
  #doesn't want to work for 'flater' priors :( -  trap error occur
  bin_score <- data$bin_score 
  time <- data$time
  group <- data$IDgroup
  
  K = data$IDgroup %>% unique() %>% length()
  N = length(bin_score)
  
  dataList = list("time","bin_score","group","K","N")
  params = c("b0","b1","a4")
  
  inits9 <- function(){
    list(b0 = rnorm(500,0,1),b1 = rnorm(500,0,1),a4 = rnorm(1,0,1))
  }
  
  
  nc <- 2    #number of MCMC chains to run
  ni <- 6000  #number of samples for each chain     
  nb <- 3000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  bugs.out9 <-R2WinBUGS::bugs(data=dataList,inits=inits9, parameters.to.save=params, model.file="model9.txt", 
                              n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE,
                              bugs.directory = "F:\\WinBUGS14", working.directory=getwd())
  print(bugs.out9, digits = 3)
  bugs.out9$summary
  bugs.out9$DIC # 2439
  
  plot(bugs.out9)
  
  
  #-------------------------------------------------------------------------
  #--            ranodm intercept and slope + fixed time                ----
  #--                          b0,b1,a4 normal                          ----
  #--              complementary log-log link function                  ----
  #-------------------------------------------------------------------------
  sink("model10.txt")        
  cat("model FIRST
      {
      # N observations
      for (i in 1:N) {
      bin_score[i] ~ dbern(p[i])
      cloglog(p[i]) <- b0[group[i]] + b1[group[i]]*time[i] + a4*time[i]
      }
      for (k in 1:K) {
      b0[k] ~ dnorm(0,1)
      b1[k] ~ dnorm(0,1)
      }
      
      a4 ~ dnorm(0,1)
      
      }", fill = TRUE)
  sink()
  
  bin_score <- data$bin_score 
  time <- data$time
  group <- data$IDgroup
  
  K = data$IDgroup %>% unique() %>% length()
  N = length(bin_score)
  
  dataList = list("time","bin_score","group","K","N")
  params = c("b0","b1","a4")
  
  inits10 <- function(){
    list(b0 = rnorm(500,0,1),b1 = rnorm(500,0,1),a4 = rnorm(1,0,1))
  }
  #trap, even with this priors 
  
  nc <- 2    #number of MCMC chains to run
  ni <- 6000  #number of samples for each chain     
  nb <- 3000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  bugs.out10 <-R2WinBUGS::bugs(data=dataList,inits=inits10, parameters.to.save=params, model.file="model10.txt", 
                              n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE,
                              bugs.directory = "F:\\WinBUGS14", working.directory=getwd())
  print(bugs.out10, digits = 3)
  bugs.out10$summary
  bugs.out10$DIC 
  
  plot(bugs.out10)