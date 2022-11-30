#Testing that the numerical calculation of loglik is correct
#library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results
nDone0=3
steptol0=1e-6
s0 = 3 #signif of checking

kit0 = "ESX17" #name of selected kit
examples = paste(path.package("euroformix"),"tutorialdata",sep=.Platform$file.sep)
popfn = paste(examples,paste0("ESX17_NorwayAMEL.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("stain.txt"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("refs.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))#,threshT)
refs = sample_tableToList(tableReader(reffn))
kitinfo = getKit(kit0) #get kitinfo

#Allign settings wrt color-marker info  
ATv = 150
pCv = 0.05
lamv = 0.01
fstv = 0.01

#plotEPG2(samples,kit0,refs, AT=ATv) #plot data
dat = prepareData(samples,refs,popFreq,threshT=ATv,minF=NULL, normalize = TRUE) #obtain data to use for analysis

NOC = 2
test_that("check maximum likelihood Hp:", {
  cond = c(1,2)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,xi=NULL,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates

  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c(0.762,0.238,1904.47,0.285,0.686,0.057) )
  expect(round(mle$fit$loglik,s0),-399.678) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  compareValid(valid$ProbObs,c(0.162,0.079,0.198,0.822,0.679,0.561,0.321,0.269,0.681,0.866,0.279,0.381,0.157,0.048,0.186,0.418,0.787,0.691,0.051,0.937,0.425,0.199,0.358,0.296,0.27,0.515,0.542,0.452,0.001,0.012,0.49,0.595,0.232,0.183,0.504,0.562,0.474,0.409,0.3,0.335,0.998,0.756,1,0.301,0.716,0.904,0.174,0.709,0.712,0.773,0.566,0.327,0.469,0.638))
})

test_that("check maximum likelihood Hd:", {
  cond = c(1,0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = 2,xi=NULL,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c(0.768,0.232,1927.452,0.305,0.676,0.061) )
  expect(round(mle$fit$loglik,s0),-426.181) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  compareValid(valid$ProbObs,c(0.241,0.023,0.111,0.814,0.766,0.584,0.328,0.268,0.641,0.842,0.223,0.329,0.383,0.007,0.089,0.412,0.823,0.689,0.045,0.917,0.429,0.219,0.365,0.314,0.281,0.628,0.467,0.419,0.009,0.001,0.425,0.549,0.374,0.193,0.542,0.531,0.472,0.337,0.442,0.299,0.994,0.959,0.999,0.39,0.695,0.879,0.143,0.681,0.654,0.886,0.562,0.347,0.456,0.638))
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle)
  expect_equal(as.numeric(DC$table2[1,5]),c(0.6688))
})

