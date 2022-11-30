#Testing that the numerical calculation of loglik is correct
#rm(list=ls());
#library(euroformix);library(testthat)
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
samples = sample_tableToList(tableReader(evidfn))[1]#Use only 1st replicate
refs = sample_tableToList(tableReader(reffn))
kitinfo = getKit(kit0) #get kitinfo

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

test_that("check maximum likelihood Hp: Ref1 + Ref2", {
  cond = c(1,2) #known contributors
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,xi=0,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, kit0, modelStutt=FALSE)
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c( 0.724,0.276,688.462,0.203,0.774) )
  expect(round(mle$fit$loglik,s0), -203.634) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  vals=c(0.94,0.716,0.914,0.201,0.372,0.095,0.465,0.967,0.427,0.16,0.072,0.687,0.124,0.213,0.761,0.527,0.438,0.12,0.798,0.069,0.317,0.551,0.855,0.279,0.258,0.843,0.961,0.782,0.457)
  #plot(sort(valid$ProbObs),sort(vals),ty="l")
  #compareValid(valid$ProbObs,vals)
  
  DC = deconvolve(mle)
  expect_equal(round(as.numeric(DC$table2[,5]),s0),rep(1,nrow(DC$table2)))
})


test_that("check maximum likelihood Hd: 2 unknown (unrelated):", {
  cond = c(0,0)
  knownRefs = c(1,2) #known non-contributors
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond, knownRef = knownRefs, xi=0,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond, pCv,ATv,fstv,lamv, kit0, modelStutt=FALSE)
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c( 0.719,0.281,699.019,0.231,0.771) )
  expect(round(mle$fit$loglik,s0), -234.096) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  vals=c(0.749,0.959,0.73,0.02,0.674,0.261,0.267,0.951,0.471,0.245,0.042,0.695,0.196,0.303,0.634,0.206,0.491,0.269,0.334,0.224,0.146,0.332,0.641,0.3,0.898,0.461,0.969,0.473,0.736)
  #plot(sort(valid$ProbObs),sort(vals),ty="l")
  #compareValid(valid$ProbObs,vals)
  
  #Calculate deconvolution (DC):
  #paste0(round(as.numeric(DC$table2[,5]),s0),collapse = ",")
  DC = deconvolve(mle)
  expect_equal(round(as.numeric(DC$table2[,5]),s0),c(0.627,0.52,0.315,0.752,0.783,0.526,0.57))
  
})


test_that("check maximum likelihood Hd (sibling):", {
  ibd0 = c(0.25,0.5,0.25) #Sibling
  knownRefs = 1:2  #typed non-contributor
  refRel = 2 #ref index of related individual
  cond = c(0,0)
  
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = knownRefs, xi=0,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0, xiFW=0, knownRel = refRel,ibd=ibd0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates

  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  expect(round(mle$fit$loglik,2), -226.54)   
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, kit0, ibd=ibd0, refRel = refRel, modelStutt=FALSE )
  expect(logLikv,logLikv2)
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle)
  expect_equal(round(as.numeric(DC$table2[,5]),s0),c(0.929,0.648,0.864,0.94,0.965,0.558,0.594)) #Prob of Sibling (last)
})


test_that("check maximum likelihood Hd (parent/child):", {
  ibd0 = c(0,1,0) #assuming the unknown is a child of Ref1
  knownRefs = 1:2
  refRel = 2 #ref index of related individual
  cond = c(0,0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = knownRefs,xi=0,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0,knownRel = refRel,ibd=ibd0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik) 
  expect(round(mle$fit$loglik,2), -227.44) #32-bit compatible   
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, kit0, ibd=ibd0, refRel = refRel, modelStutt=FALSE )
  expect(logLikv,logLikv2)
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle)
  expect_equal(round(as.numeric(DC$table2[,5]),s0),c(0.825,0.739,0.462,0.892,0.929,0.549,0.6)) #Prob of Sibling (last)
})


