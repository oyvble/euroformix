#Testing that the numerical calculation of loglik is correct
#library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results
nDone0=3
steptol0=1e-6
s0 = 3 #signif of checking

kit0 = "testkit" #name of selected kit
examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0("test_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("test_evid1.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("test_ref1.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))#,threshT)
refs = sample_tableToList(tableReader(reffn))
kitinfo = getKit(kit0) #get kitinfo

#Set specific settings for each dye:
cols = c("blue","yellow")
ATdye = c(50,70)
pCdye = c(0 , 0.0097)
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


 NOC = 1
test_that("check maximum likelihood Hp:", {
  cond = 1
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,xi=NULL,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=NULL, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates

  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, kit0)
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c( 1,921.894,0.099,0.771,0.104,0.051) )
  expect(round(mle$fit$loglik,s0),-134.043) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  compareValid(valid$ProbObs,c(0.717,0.975,0.667,0.387,0.579,0.611,0.114,0.298,0.085,0.381,0.434,0.045,0.087,0.204,0.272,0.09,0.441,0.953,0.163,0.781,0.218,0.029,0.644,0.989,0.297))
})

test_that("check maximum likelihood Hd (unrelated):", {
  cond = 0
  knownRef = 1 #index of typed known non-contributor
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef=knownRef,xi=NULL,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=NULL, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, kit0,knownRef=knownRef)
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c( 1,921.887,0.099,0.771,0.104,0.051) )
  expect(round(mle$fit$loglik,s0),-155.291) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  compareValid(valid$ProbObs,c(0.717,0.975,0.667,0.387,0.579,0.611,0.114,0.298,0.085,0.381,0.434,0.045,0.087,0.204,0.272,0.09,0.441,0.953,0.163,0.781,0.218,0.029,0.644,0.989,0.297))
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle)
  expect_equal(as.numeric(DC$table2[,2]),c(1,1,1,1,1,1,1))
})


test_that("check maximum likelihood Hd (sibling):", {
  ibd0 = c(1/4,1/2,1/4) #assuming the unknown is a sibling of poi
  knownRel = 1
  cond=0
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond, xi=NULL,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=NULL,knownRel = knownRel,ibd=ibd0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #mle$prepareC$relGind
  #mle$prepareC$ibd
  
  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, kit0, ibd0=ibd0, knownRel=knownRel, knownRef=knownRel)
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c( 1,921.894,0.099,0.771,0.104,0.051) )
  expect(round(mle$fit$loglik,s0), -141.185) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  compareValid(round(valid$ProbObs,s0),c(0.717,0.975,0.667,0.387,0.579,0.611,0.114,0.298,0.085,0.381,0.434,0.045,0.087,0.204,0.272,0.09,0.441,0.953,0.163,0.781,0.218,0.029,0.644,0.989,0.297))
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle)
  expect_equal(as.numeric(DC$table2[,2]),c(1,1,1,1,1,1,1))
})

