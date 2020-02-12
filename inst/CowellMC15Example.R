rm(list=ls())
#install.packages("euroformix", repos="http://R-Forge.R-project.org") 
library(euroformix)

#load data:
load(system.file("mcData.Rdata", package = "euroformix"))
samples <- list(MC15=mcdata$samples$MC15) #use only sample MC15

retlist <- Qassignate(samples,mcdata$popFreq,mcdata$refData) #get Q-assignated allele frequiences
popFreqQ <- retlist$popFreq
refDataQ <- retlist$refData

#Model:
threshT = 50 #peak height threshold
condHp <- c(1,2,3) #Conditioned references under hp
condHd <- c(1,2,0) #Conditioned references under hd
nC <- 3 #number of contributors

#NON Q-assignated data:
set.seed(1)
nDone <- 1 #number of random optimization start points
hptime <- system.time( {  hpSmle <- contLikMLE(nC,samples,mcdata$popFreq,mcdata$refData,condOrder=condHp,threshT=threshT,nDone=nDone) } )[3]
print(hpSmle$fit$thetahat2) #MLE under hp
print(hpSmle$fit$thetaSE) #standard error
set.seed(1)
hdtime <- system.time( {  hdSmle <- contLikMLE(nC,samples,mcdata$popFreq,mcdata$refData,condOrder=condHd,threshT=threshT,nDone=nDone) } )[3] #slow for replicate version
print(hdSmle$fit$thetahat2) #MLE under hp
print(hdSmle$fit$thetaSE) #standard error
print(c(hptime,hdtime)) #time usage
LRmle <- exp(hpSmle$fit$loglik-hdSmle$fit$loglik) #MLE optimized 
print(log10(LRmle))

#For Q-assignated data: (much faster than using all data)
set.seed(1)
nDone <- 3 #number of random optimization start points
hptime <- system.time( {  hpSmleQ <- contLikMLE(nC,samples,popFreqQ,refDataQ,condOrder=condHp,threshT=threshT,nDone=nDone) } )[3] #
print(hpSmleQ$fit$thetahat2) #MLE
print(hpSmleQ$fit$thetaSE) #standard error
hdtime <- system.time( {  hdSmleQ <- contLikMLE(nC,samples,popFreqQ,refDataQ,condOrder=condHd,threshT=threshT,nDone=nDone) } )[3] #
print(c(hptime,hdtime))
print(hdSmleQ$fit$thetahat2) #MLE
print(hdSmleQ$fit$thetaSE) #standard error
LRmleQ <- exp(hpSmleQ$fit$loglik-hdSmleQ$fit$loglik) #MLE optimized 
print(log10(LRmleQ))

#Integrated likelihood
reltol <- 0.05 #relative error
low <- rep(0,nC+2) #lower boundary of integral
up <- rep(1,nC+2) #upper boundary of integral
up[nC] <- 10000
up[nC+1] <- 1
hptime <- system.time( {  hpSintQ <- contLikINT(nC,samples,popFreqQ,low,up,refDataQ,condOrder=condHp,threshT=threshT,reltol=reltol) } )[3]
hdtime <- system.time( {  hdSintQ <- contLikINT(nC,samples,popFreqQ,low,up,refDataQ,condOrder=condHd,threshT=threshT,reltol=reltol) } )[3]
LRintQ <- hpSintQ$margL/hdSintQ$margL #MLE optimized 
print(log10(LRintQ))




