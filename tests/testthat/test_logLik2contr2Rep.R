#Testing that the numerical calculation of loglik is correct
#rm(list=ls());library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results
nDone0=2
steptol0=1e-6
s0 = 3 #signif of checking

kit0 = "testkit" #name of selected kit
examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0("test_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("test_evids.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("test_refs.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))#,threshT)
refs = sample_tableToList(tableReader(reffn))

#Set specific settings for each dye:
cols = c("blue","yellow")
ATdye = c(50,70)
pCdye = c(0.0133 , 0.0097)
lamdye = c(0.025, 0.034)
fstdye = c(0.01,0.02)

#Allign settings wrt color-marker info  
kitdyes = getKit(kit0,"COLOR") #get kitinfo from selected kit
indmatch = match(kitdyes$Color,cols)
ATv = ATdye[indmatch] 
pCv = pCdye[indmatch]
lamv = lamdye[indmatch]
fstv = fstdye[indmatch]
names(ATv) <- names(pCv) <- names(lamv) <- names(fstv)  <- toupper(kitdyes$Marker) #don't need to be same order as of popFreq

#plotEPG2(samples,kit0,refs, AT=ATv) #plot data
dat = prepareData(samples,refs,popFreq,threshT=ATv,minF=NULL, normalize = TRUE) #obtain data to use for analysis
NOC = 2

test_that("check maximum likelihood Hp:", {
  cond = c(1,2)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,xi=NULL,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=NULL,seed=seed0,nDone=nDone0,steptol=steptol0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates

  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)

  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, kit0)
  expect(logLikv,logLikv2)

  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(signif(thhat,s0),c( 0.796,0.204,795,0.151,0.747,0.119,0.055) )
  expect(round(mle$fit$loglik,s0),-330.464) 

  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  #valid2 = getValidProbs(mle)
  compareValid(valid$ProbObs,c(0.757,0.216,0.622,0.26,0.936,0.833,0.186,0.021,0.255,0.113,0.064,0.193,0.606,0.373,0.99,0.935,0.728,0.972,0.355,0.096,0.748,0.467,0.396,0.222,0.584,0.096,0.825,0.95,0.36,0.25,0.634,0.317,0.87,0.031,0.082,0.619,0.713,0.016,0.959,0.631,0.509,0.697,0.876,0.654,0.165,0.242,0.742,0.739,0.112,0.767,0.626,0.974,0.907,0.556,0.584,0.343,0.569))
})

test_that("check maximum likelihood Hd (unrelated):", {
  cond = c(1,0)
  knownRef = 2
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = knownRef,xi=NULL,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=NULL, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, kit0)
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,2),c(0.76,0.24,796.62,0.15,0.74,0.12,0.05) ) #DIFFERED FOR 32 and 64 bit!
  expect(round(mle$fit$loglik,s0),-349.573) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  #valid2 = getValidProbs(mle)
  compareValid(valid$ProbObs,c(0.876,0.213,0.694,0.31,0.908,0.77,0.16,0.015,0.279,0.117,0.043,0.169,0.586,0.351,0.99,0.919,0.65,0.957,0.338,0.104,0.813,0.552,0.416,0.103,0.398,0.128,0.871,0.898,0.212,0.285,0.682,0.385,0.909,0.042,0.107,0.541,0.644,0.023,0.974,0.508,0.446,0.681,0.887,0.637,0.215,0.306,0.89,0.991,0.701,0.692,0.528,0.78,0.528,0.535,0.563,0.412,0.642))
  
  #CHECK deconvolution (DC):
  DC = deconvolve(mle)
  expect( round(as.numeric(DC$table2[,5]),s0),c(0.925 ,0.875, 0.543, 1.000 ,0.999 ,0.807, 0.983))
})


test_that("check maximum likelihood Hd (sibling):", {
  ibd0 = c(1/4,1/2,1/4) #assuming the unknown is a sibling of poi
  knownRef = 2
  refRel = 2 #ref index of related individual
  cond = c(1,0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = knownRef,xi=NULL,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=NULL,knownRel = refRel,ibd=ibd0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, kit0, ibd=ibd0, refRel = refRel )
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c(0.782,0.218,794.674,0.152,0.746,0.119,0.051) )
  expect(round(mle$fit$loglik,s0),-337.189) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  #valid2 = getValidProbs(mle)
  compareValid(valid$ProbObs,c(0.812,0.213,0.652,0.287,0.925,0.81,0.178,0.019,0.253,0.109,0.062,0.195,0.606,0.373,0.991,0.939,0.689,0.963,0.322,0.097,0.768,0.495,0.402,0.166,0.503,0.111,0.842,0.929,0.295,0.27,0.659,0.346,0.884,0.036,0.094,0.589,0.686,0.02,0.964,0.575,0.466,0.661,0.886,0.671,0.183,0.264,0.753,0.985,0.366,0.727,0.565,0.94,0.791,0.548,0.576,0.371,0.597))
  
  #CHECK deconvolution (DC):
  DC = deconvolve(mle)
  #paste0( round(as.numeric(DC$table2[,5]),s0) ,collapse = ",")
  expect( round(as.numeric(DC$table2[,5]),s0),c(0.976,0.955,0.967,1,1,0.978,0.553))
})
