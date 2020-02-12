rm(list=ls())
#install.packages("euroformix", repos="http://R-Forge.R-project.org") 
library(euroformix); #sessionInfo()

#model:
nC <- 2 #number of contributors
threshT = 150 #detection threshold (default is 50)

#true parameters to generate data from
xi=0.1 #true stutter ratio in generated samples
mu = 1500 #amount of dna (on peak height scale)
sigma=0.15 #coefficient of variance
nrep=2 #number of replicates (with same assumed parameters)
mx=rev(nC:1)/sum(1:nC) #mixture propotions for each of the contributors

#load allele-frequncy data:
load(system.file("mcData.Rdata", package = "euroformix"))
popFreq <- mcdata$popFreq #get object with allele-frequncies

#Generate samples using popFreq-object and assumed model parameters
set.seed(1) 
gendata <- genDataset(nC,popFreq,mu,sigma,mx=mx,threshT=threshT,nrep=nrep,stutt=xi)
samples = gendata$samples #generated mixture-samples
print( sapply(samples,function(x) sapply(x,function(x) length(x$adata)))) #show number of allele observations for each samples
retlist <- Qassignate(samples,popFreq=popFreq,refData=gendata$refData) #Q-assignation (speeds up calculations)
popFreq <- retlist$popFreq
refData <- retlist$refData 
nR <- length(refData[[1]]) #number of references in refData-object
print(sum(unlist(refData)=="99")) #number of dropout from references
tabref <- matrix(unlist(refData),ncol=2*nC,byrow=TRUE)
rownames(tabref) <- names(popFreq)

#Hypothesis: Hp: ref2+1unknown , Hd:2 unknowns
hpH <- hdH <- rep(0,nR) #NB: this vector must be of size nR!
hpH[1] <- 1 #condition ref2 to be placed as contributor 1 in model

###################################################
#Maximum likelihood estimation of parameter models#
###################################################
nDone=5 #Number of random start points in the maximizer(default is 1)
hdmle <- contLikMLE(nC,samples,popFreq,refData=refData,condOrder=hdH,threshT=threshT,nDone=nDone) 
print(hdmle$fit$loglik) #fitted log-likelihood
hdmlelogliki <- logLiki(hdmle) #get loglikelihood for each locus under hd

#given hdDmle-object: Do deconvolution:
dlisthd <- deconvolve(mlefit=hdmle,alpha=0.99,maxlist=100) #time depends on discrimination of sample
#barplot(dlisthd$pG)
print(tabref)
print(head(dlisthd$table1))

set.seed(1)
hpmle <- contLikMLE(nC,samples,popFreq,refData=refData,condOrder=hpH,threshT=threshT,nDone=nDone)
print(hpmle$fit$loglik) #fitted log-likelihood
hpmlelogliki <- logLiki(hpmle)  #get loglikelihood for each locus under hp

#Weight-of-Evidence:
LRmlei <- exp(hpmlelogliki - hdmlelogliki) #MLE LR for each locus
LRmle <- exp(hpmle$fit$loglik - hdmle$fit$loglik) #MLE LR = prod(LRmlei)

###################################
#MCMC Metropolis Hasting sampling:#
##Used to explore parameter space##
###################################
niter <- 1e4 #number of posterior-samples
#Under Hd
delta <- 15 #try find a good behaved delta such that acceptance rate is around 0.25 and acf low
hdmcmc <- contLikMCMC(hdmle,niter,delta)
print(hdmcmc$accrat) #acceptance rate of sampler
validMCMC(hdmcmc)

#Under Hp
delta <- 15 #try find a good behaved delta such that acceptance rate is around 0.25 and acf low
hpmcmc <- contLikMCMC(hpmle,niter,delta)
print(hpmcmc$accrat) #acceptance rate of sampler
validMCMC(hpmcmc)

LRmcmc <- hpmcmc$margL/hdmcmc$margL #get estimated LR based on MCMC simulations
#dp <- density(hpmcmc$postlogL)
#dd <- density(hdmcmc$postlogL)
#matplot(cbind(dp$x,dd$x),cbind(dp$y,dd$y),ty="l")
#plot(density((hpmcmc$postlogL - hdmcmc$postlogL)/log(10)))

#######################
#Numerical Integration#
#######################
#specify parameter limits of the integral (we consider it after the MCMC exploring)
#From posterior plot we assume: mu inside [400,800], sigma in [0,0.3]
#theta=(mx1,mu,sigma,xi) is the order of the parameters
hdL <- hpL <- c(rep(0,nC-1),0,0,0,0)
hdU <- hpU <- c(rep(1,nC-1),2000,0.5,1,0.4) #assume upper limit of stutter ratio as 0.4
reltol <- 0.1 #required relative error of the integration estimate
hpint <- contLikINT(nC,samples,popFreq,lower=hpL,upper=hpU,refData=refData,condOrder=hpH,threshT=threshT,reltol=reltol)
hdint <- contLikINT(nC,samples,popFreq,lower=hdL,upper=hdU,refData=refData,condOrder=hdH,threshT=threshT,reltol=reltol)
LRint <- hpint$margL/hdint$margL #The LR weight-of-evidence based on integrated likelihood
LRdev <- c(hpint$dev[1]/hdint$dev[2],hpint$dev[2]/hdint$dev[1]) #error-interval of the integrated based LR

tab <- log10(rbind(LRmle,LRint,LRmcmc))
colnames(tab) <- "log10LR"

