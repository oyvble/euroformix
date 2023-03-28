#Testing of whether the numerical calculation of loglik (maximum likelihood approach) is correct
#BOTH QUALITATIVE AND QUANTITATIVE MODEL IS TESTED

#library(euroformix);library(testthat)
#rm(list=ls());

kit0 = "SGMplus" #name of selected kit
s0 = 3 #signif
#Helpfunction for checking model after model fit
examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0(kit0,"_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0(kit0,"_evids_LTDNA.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0(kit0,"_refs_LTDNA.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))#,threshT)
refData = sample_tableToList(tableReader(reffn))[1] #use only first ref
kitinfo = getKit(kit0) #get kitinfo

ATv = 100
pCv = 0.05
fstv = 0.01
lamv = 0.01
seed0 = 1

#plotEPG2(samples,kit0, AT=ATv) #plot data
dat = prepareData(samples,refData,popFreq,threshT=ATv,minF=NULL, normalize = TRUE) #obtain data to use for analysis

#Hypothesis:
NOC = 2
condhp = c(1,0)
condhd = c(0,0)
knownRef = 1


test_that("Calculating quantitative model (3 reps):", {
  likhp = calcMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhp,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE)$fit$loglik
  likhd = calcMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhd,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE, knownRef = knownRef)$fit$loglik
  LR = (likhp-likhd)/log(10)
  expect(round(likhp,s0),-221.057)
  expect(round(likhd,s0),-240.318)
  expect(round(LR,s0),8.365)
  
  likhp2 = calcMLE(nC=NOC,samples=samples,popFreq=popFreq,refData=refData,condOrder=condhp,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE)$fit$loglik
  likhd2 = calcMLE(nC=NOC,samples=samples,popFreq=popFreq,refData=refData,condOrder=condhd,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE, knownRef = knownRef)$fit$loglik
  expect(round(likhp,s0),round(likhp2,s0))
  expect(round(likhd,s0),round(likhd2,s0))
  
  likhp3 = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhp,threshT=ATv,prC=pCv,fst=fstv,lambda = lamv ,kit=kit0,xi=0,xiFW=0)$fit$loglik
  likhd3 = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhd,threshT=ATv,prC=pCv,fst=fstv,lambda = lamv ,kit=kit0,xi=0,xiFW=0, knownRef = knownRef)$fit$loglik
  expect(round(likhp,s0),round(likhp3,s0))
  expect(round(likhd,s0),round(likhd3,s0))
})


test_that("Calculating valid/DC of model (3 reps):", {
  #Fit mle again
  mleHp = calcMLE(nC=NOC,samples=samples,popFreq=popFreq,refData=refData,condOrder=condhp,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE,seed=seed0)
  mleHd = calcMLE(nC=NOC,samples=samples,popFreq=popFreq,refData=refData,condOrder=condhd,AT=ATv,pC=pCv,fst=fstv,lambda = lamv ,kit=kit0, DEG = TRUE,BWS = FALSE, FWS=FALSE, knownRef = knownRef,seed=seed0)

  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mleHp,kit=kit0,createplot=FALSE,verbose = FALSE)
  compareValid(valid$ProbObs,c(0.085,0.229,0.868,0.335,0.98,0.379,0.72,0.451,0.88,0.056,0.527,0.011,0.035,0.737,0.392,0.526,0.305,0.322,0.568,0.477,0.108,0.073,0.324,0.997,0.296,0.309,0.088,0.295,0.879,0.86,0.506,0.016,0.244,0.373,0.159,0.166,0.515,0.252))
  
  valid = validMLEmodel(mleHd,kit=kit0,createplot=FALSE,verbose = FALSE)
  compareValid(valid$ProbObs,c(0.1,0.253,0.875,0.329,0.975,0.371,0.685,0.418,0.873,0.055,0.517,0.011,0.035,0.731,0.396,0.529,0.308,0.308,0.55,0.46,0.086,0.058,0.409,0.999,0.378,0.297,0.087,0.289,0.867,0.863,0.524,0.018,0.265,0.379,0.162,0.154,0.491,0.202))
  
  #CHECK deconvolution (DC):
  #paste0(round(as.numeric(DC$table2[,2]),s0),collapse = ",")
  DC = deconvolve(mleHp)
  compareValid(as.numeric(DC$table2[,5]),c(0.99,0.408,1,1,0.853,0.742,0.699,0.977,0.562,0.451))
  
  DC = deconvolve(mleHd)
  compareValid(as.numeric(DC$table2[,2]),c(0.999,0.917,1,1,0.987,0.708,0.631,0.999,0.661,0.757))
})


test_that("Calculating qualitative model (all reps):", {
  #nC=NOC;samples=dat1$samples;popFreq=dat1$popFreq;refData=dat1$refData;condOrder=condhp;prC=pCv;fst=fstv
  mlehp = calcQualMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhp,prC=pCv,fst=fstv)
  mlehd = calcQualMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=condhd,prC=pCv,fst=fstv, knownRef = knownRef)
  likhp = mlehp$loglik
  likhd = mlehd$loglik
  LR = (likhp-likhd)/log(10)
  expect(round(likhp,s0),-56.913)
  expect(round(likhd,s0),-70.584)
  expect(round(LR,s0),5.937)
  
  #CHECK MANUAL CALCULATIONS:
  checkQualModel(mlehp) #check per-marker numeric
  checkQualModel(mlehd) #check per-marker numeric
})


