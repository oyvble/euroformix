#Testing of whether the numerical calculation of loglik (maximum likelihood approach) is correct
#BOTH QUALITATIVE AND QUANTITATIVE MODEL IS TESTED
#library(euroformix);library(testthat)

kit0 = "ESX17" #name of selected kit
s0 = 3 #signif
#Helpfunction for checking model after model fit
examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0("ESX17_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("ESX17_evids_LTDNA.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("ESX17_refs_LTDNA.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))#,threshT)
samples1 = samples[1] #only first
refData = sample_tableToList(tableReader(reffn))
kitinfo = getKit(kit0) #get kitinfo

ATv = 100
pCv = 0.05
fstv = 0.01
lamv = 0.01

#plotEPG2(samples,kit0,refs, AT=ATv) #plot data
dat = prepareData(samples,refData,popFreq,threshT=ATv,minF=NULL, normalize = TRUE) #obtain data to use for analysis
dat1 = prepareData(samples1,refData,popFreq,threshT=ATv,minF=NULL, normalize = TRUE) #obtain data to use for analysis

#Hypothesis:
NOC = 3
condhp = c(1,2,3)
condhd = c(1,2,0)
knownRef = 3

test_that("Calculating quantitative model (1 rep):", {
#  nC=NOC;samples=dat1$samples;popFreq=dat1$popFreq;refData=dat1$refData;condOrder=condhp;AT=ATv;pC=pCv;fst=fstv;lambda = lamv ;kit=kit0; DEG = TRUE;BWS = FALSE; FWS=FALSE
  mlefit = calcMLE(nC=NOC,samples=dat1$samples,popFreq=dat1$popFreq,refData=dat1$refData,condOrder=condhp,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE)
  #Check DC
  DC = deconvolve(mlefit)
  
  likhp = mlefit$fit$loglik
  likhd = calcMLE(nC=NOC,samples=dat1$samples,popFreq=dat1$popFreq,refData=dat1$refData,condOrder=condhd,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE, knownRef = knownRef)$fit$loglik
  LR = (likhp-likhd)/log(10)
  expect(round(likhp,s0),-245.016)
  expect(round(likhd,s0),-247.941)
  expect(round(LR,s0),1.271)

  likhp2 = calcMLE(nC=NOC,samples=samples1,popFreq=popFreq,refData=refData,condOrder=condhp,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE)$fit$loglik
  likhd2 = calcMLE(nC=NOC,samples=samples1,popFreq=popFreq,refData=refData,condOrder=condhd,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE, knownRef = knownRef)$fit$loglik
  expect(round(likhp,s0),round(likhp2,s0))
  expect(round(likhd,s0),round(likhd2,s0))
  
  likhp3 = contLikMLE(nC=NOC,samples=dat1$samples,popFreq=dat1$popFreq,refData=dat1$refData,condOrder=condhp,threshT=ATv,prC=pCv,fst=fstv,lambda = lamv ,kit=kit0,xi=0,xiFW=0)$fit$loglik
  likhd3 = contLikMLE(nC=NOC,samples=dat1$samples,popFreq=dat1$popFreq,refData=dat1$refData,condOrder=condhd,threshT=ATv,prC=pCv,fst=fstv,lambda = lamv ,kit=kit0,xi=0,xiFW=0, knownRef = knownRef)$fit$loglik
  expect(round(likhp,s0),round(likhp3,s0))
  expect(round(likhd,s0),round(likhd3,s0))
})


test_that("Calculating quantitative model (3 reps):", {
  likhp = calcMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhp,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE)$fit$loglik
  likhd = calcMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhd,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE, knownRef = knownRef)$fit$loglik
  LR = (likhp-likhd)/log(10)
  expect(round(likhp,s0),-792.665)
  expect(round(likhd,s0),-810.009)
  expect(round(LR,s0),7.533)
  
  likhp2 = calcMLE(nC=NOC,samples=samples,popFreq=popFreq,refData=refData,condOrder=condhp,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE)$fit$loglik
  likhd2 = calcMLE(nC=NOC,samples=samples,popFreq=popFreq,refData=refData,condOrder=condhd,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE, knownRef = knownRef)$fit$loglik
  expect(round(likhp,s0),round(likhp2,s0))
  expect(round(likhd,s0),round(likhd2,s0))
  
  likhp3 = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhp,threshT=ATv,prC=pCv,fst=fstv,lambda = lamv ,kit=kit0,xi=0,xiFW=0)$fit$loglik
  likhd3 = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhd,threshT=ATv,prC=pCv,fst=fstv,lambda = lamv ,kit=kit0,xi=0,xiFW=0, knownRef = knownRef)$fit$loglik
  expect(round(likhp,s0),round(likhp3,s0))
  expect(round(likhd,s0),round(likhd3,s0))
})



#fstv = 0
#condhp = c(0,1,2)
#condhd = c(0,1,0)
#knownRef = 3


test_that("Calculating qualitative model (1 rep):", {
  #nC=NOC;samples=dat1$samples;popFreq=dat1$popFreq;refData=dat1$refData;condOrder=condhp;prC=pCv;fst=fstv
  mlehp = calcQualMLE(nC=NOC,samples=dat1$samples,popFreq=dat1$popFreq,refData=dat1$refData,condOrder=condhp,prC=pCv,fst=fstv)
  mlehd = calcQualMLE(nC=NOC,samples=dat1$samples,popFreq=dat1$popFreq,refData=dat1$refData,condOrder=condhd,prC=pCv,fst=fstv, knownRef = knownRef)
  likhp = mlehp$loglik
  likhd = mlehd$loglik
  LR = (likhp-likhd)/log(10)
  expect(round(likhp,s0),-50.975)
  expect(round(likhd,s0),-49.913)
  expect(round(LR,s0),-0.461)
  
  #CHECK MANUAL CALCULATIONS:
  checkQualModel(mlehp) #check per-marker numeric
  checkQualModel(mlehd) #check per-marker numeric
})

test_that("Calculating qualitative model (all reps):", {
  #nC=NOC;samples=dat1$samples;popFreq=dat1$popFreq;refData=dat1$refData;condOrder=condhp;prC=pCv;fst=fstv
  mlehp = calcQualMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhp,prC=pCv,fst=fstv)
  mlehd = calcQualMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhd,prC=pCv,fst=fstv, knownRef = knownRef)
  likhp = mlehp$loglik
  likhd = mlehd$loglik
  LR = (likhp-likhd)/log(10)
  expect(round(likhp,s0),-146.636)
  expect(round(likhd,s0),-153.451)
  expect(round(LR,s0),2.96)
  
  #CHECK MANUAL CALCULATIONS:
  checkQualModel(mlehp) #check per-marker numeric
  checkQualModel(mlehd) #check per-marker numeric
})


