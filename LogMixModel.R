library(tidyverse)
library(data.table)
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
results <- glmer(bin_score~(1|IDgroup)+(time|IDgroup),data = data,family = binomial(),control=glmerControl(optimizer="bobyqa",optCtrl=list(maxfun=2e5)))
results
plot(results) #AIC = 1629

results2 <- glm(bin_score~as.factor(IDgroup)+time*as.factor(IDgroup),data = data, family = binomial())
results2
plot(results2) #AIC = 1413 

#Bayesian
library(R2OpenBUGS)
library(R2WinBUGS)

sink("model.txt")        
cat("model FIRST
 {
  # N observations
  for (i in 1:N) {
    bin_score[i] ~ dbern(p[i])
    logit(p[i]) <- b0[group[i]] + b1[group[i]]*time[i]
  }
  # Nsubj groups
  for (j in 1:Nsubj) {
    b0[j] ~ dnorm(0,1000)
    b1[j] ~ dnorm(0,1000)
  }
 }", fill = TRUE)
  sink()

  bin_score <- data$bin_score 
  time <- data$time
  group <- data$IDgroup
  Nsubj = length(unique(data$IDgroup))
  N = length(bin_score)
  
  dataList = list("time","bin_score","group","Nsubj","N")
  params = c("b0","b1")
  
  inits <- function(){
     list(b0 = rnorm(Nsubj,0,11),b1 = rnorm(Nsubj,0,0.1))
     }
  
  
  nc <- 3    #number of MCMC chains to run
  ni <- 10000  #number of samples for each chain     
  nb <- 5000   #number of samples to discard as burn-in
  nt <- 1      #thinning rate, increase this to reduce autocorrelation
  
  bugs.out <-R2WinBUGS::bugs(data=dataList,inits=NULL, parameters.to.save=params, model.file="model.txt", 
                   n.chains=nc, n.iter=ni, n.burnin=nb, n.thin=nt, debug=TRUE, DIC=TRUE,
                   bugs.directory = "F:\\WinBUGS14", working.directory=getwd())
 