#Testing that the numerical calculation of loglik is correct
#rm(list=ls());library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results
nDone0=3
steptol0=1e-6
s0 = 3 #signif of checking

kit0=NULL
examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0("test_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("test_evid1.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("test_ref1.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))#,threshT)
refs = sample_tableToList(tableReader(reffn))

#Restrict outcome to those observed (no Q-allele)
markerDO = "D19S433"  #"D3S1358" #marker to perform "incident" on
popFreq[[markerDO]] = popFreq[[markerDO]][names( popFreq[[markerDO]])%in%samples[[1]][[markerDO]]$adata]
popFreq[[markerDO]] = popFreq[[markerDO]]/sum(popFreq[[markerDO]])

#Set specific settings for each dye:
AT = 50
pC = 0.05 
lam =0.01 
fst = 0.01 

#plotEPG2(samples,kit0,refs, AT=ATv) #plot data
dat = prepareData(samples,refs,popFreq,threshT=AT,minF=NULL, normalize = TRUE) #obtain data to use for analysis

#Need marker specific settings: even if constant
init = function(x) setNames(rep(x,length(popFreq)),names(popFreq))
ATv = init(AT)
pCv = init(pC)
lamv =  init(lam)
fstv =  init(fst)

NOC = 1
test_that("check maximum likelihood Hp:", {
  cond = 1
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,xi=NULL,prC=pC,threshT=AT,fst=fst,lambda=lam,xiFW=NULL, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)

  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, modelDEG=FALSE)
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c( 1,755.989,0.166,0.11,0.044) )
  expect(round(mle$fit$loglik,s0),  -146.148) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  compareValid(valid$ProbObs,c(0.665,0.957,0.665,0.242,0.253,0.29,0.052,0.039,0.003,0.55,0.857,0.537,0.323,0.463,0.423,0.293,0.472,0.781,0.336,0.484,0.323,0.456,0.759,0.992,0.241))
})

test_that("check maximum likelihood Hd (unrelated):", {
  cond=0
  knownRef = 1
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef=knownRef,xi=NULL,prC=pC, threshT=AT,fst=fst,lambda=lam,xiFW=NULL, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates

  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, modelDEG=FALSE, knownRef=knownRef)
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c( 1,755.989,0.166,0.11,0.044) )
  expect(round(mle$fit$loglik,s0), -167.102) #-167.322 if instead D3S1358 is used
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  compareValid(valid$ProbObs,c(0.665,0.957,0.665,0.242,0.253,0.29,0.051,0.039,0.003,0.55,0.857,0.541,0.323,0.463,0.423,0.293,0.472,0.781,0.336,0.484,0.323,0.457,0.759,0.992,0.241))
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle)
  expect(as.numeric(DC$table2[,2]),c(1,1,1,1,1,1,1))
})


