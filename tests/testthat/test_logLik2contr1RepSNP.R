#Testing that the numerical calculation of loglik is correct (for SNP data)
#rm(list=ls());library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results
nDone0=4
steptol0=1e-6
s0 = 3 #signif of checking

kit0 = NULL #No kit
examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0("SNP_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("SNP_evids.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("SNP_refs.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))[1]
refs = sample_tableToList(tableReader(reffn))
locs = names(popFreq)#[10:20]

#Set common settings for all markers:
AT = 500 #put high to remove alleles in markers
pC = 0.05
lam = 0.01
fst = 0.02

init = function(x) setNames(rep(x,length(popFreq)),names(popFreq))
ATv = init(AT)
pCv = init(pC)
lamv =  init(lam)
fstv =  init(fst)



test_that("Test that wrong data format causes error:", {

#Filter data based on detection threshold
samples2 =samples[[1]] #copy data
for(loc in names(samples2)) {
  keep =  samples2[[loc]]$hdata>=ATv[loc]
  samples2[[loc]]$hdata = samples2[[loc]]$hdata[keep]
  samples2[[loc]]$adata = samples2[[loc]]$adata[keep]
}

errorOccured = FALSE 
tryCatch({
  mle = contLikMLE(nC=2,samples=list(Evid=samples2),popFreq=popFreq,xi=0,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
},error = function(e) errorOccured <<- TRUE) #print(e))
expect_equal(errorOccured,TRUE) #CHECK that function call caused error
})

#plotEPG2(samples,kit0,refs, AT=ATv) #plot data
dat = prepareData(samples,refs,popFreq[locs],threshT=ATv,minF=NULL, normalize = TRUE) #obtain data to use for analysis
#plotMPS2(dat$samples,AT = ATv,refData = dat$refData,locYmax = FALSE)
#sapply(dat$refData,function(x) lapply(x, function(y) length(y)))

NOC = 2

test_that("check maximum likelihood Hp (P1 + 1U):", {
  cond = c(1,0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = 2, xi=0,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, modelStutt=FALSE,  modelDEG=FALSE)
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c(0.801,0.199,730.601,0.632) )
  expect(round(mle$fit$loglik,s0),-1409.52) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  expect(round(valid$ProbObs,s0),c(0.588,0.616,0.25,0.909,0.683,0.84,0.15,0.159,0.029,0.665,0.725,0.135,0.433,0.458,0.051,0.525,0.721,0.746,0.001,0.108,0.753,0.303,0.279,0.132,0.224,0.728,0.781,0.752,0.738,0.45,0.953,0.639,0.153,0.501,0.024,0.623,0.617,0.586,0.385,0.307,0.76,0.712,0.522,0.722,0.041,0.941,0.597,0.588,0.096,0.708,0.18,0.035,0.722,0.056,0.096,0.815,0.771,0.651,0.044,0.659,0.249,0.48,0.022,0.123,0.92,0.311,0.797,0.065,0.641,0.164,0.568,0.162,0.933,0.777,0.239,0.188,0.153,0.252,0.764,0.688,0.288,0.193,0.628,0.381,0.035,0.592,0.213,0.752,0.322,0.338,0.021,0.491,0.901,0.874,0.588,0.284,0.398,0.373,0.834,0.402,0.815,0.737,0.731,0.79,0.013,0.668,0.577,0.531,0.48,0.877,0.809,0.587,0.326,0.492,0.607,0.421,0.018,0.348,0.678,0.096,0.708,0.449,0.825,0.628,0.763,0.82,0.756,0.072,0.299,0.621,0.459,0.194,0.35,0.894,0.122,0.265,0.538,0.321,0.911,0.699,0.085,0.534,0.851,0.712,0.464,0.468,0.389,0.872,0.587,0.776,0.835,0.339,0.431,0.598,0.599,0.064,0.965,0.641,0.651,0.967,0.489,0.662,0.732,0.614,0.162,0.579,0.087,0.344,0.16,0.537,0.744,0.743,0.366,0.674,0.726))
  #valid$Significant
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  #paste0( round(as.numeric(DC$table2[1:20,5]),s0),collapse = ",")
  expect(round(as.numeric(DC$table2[1:20,5]),s0),c(0.518,0.5,0.569,0.618,0.472,0.506,0.653,0.456,0.521,0.486,1,0.506,0.593,0.466,0.503,0.509,0.505,0.516,1,1))
})

test_that("check maximum likelihood Hd: 2 unknown (unrelated):", {
  NOC = 2
  knownRefs = 1:2 #known non-contributors
  cond = c(0,0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,knownRef = knownRefs,xi=0,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, modelStutt=FALSE,  modelDEG=FALSE)
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c(0.5,0.5,744.064,0.67) )
  expect(round(mle$fit$loglik,s0),-1470.835) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  expect(round(valid$ProbObs,s0),c(0.567,0.581,0.366,0.751,0.834,0.858,0.145,0.128,0.061,0.701,0.636,0.148,0.406,0.426,0.11,0.655,0.759,0.804,0.001,0.075,0.693,0.31,0.491,0.057,0.138,0.813,0.845,0.725,0.718,0.475,0.927,0.708,0.169,0.412,0.026,0.723,0.57,0.577,0.353,0.282,0.748,0.676,0.499,0.692,0.086,0.92,0.6,0.412,0.167,0.623,0.194,0.09,0.637,0.06,0.101,0.751,0.844,0.551,0.048,0.526,0.103,0.732,0.039,0.17,0.897,0.305,0.709,0.073,0.594,0.155,0.38,0.268,0.949,0.843,0.335,0.282,0.121,0.249,0.768,0.631,0.26,0.173,0.824,0.187,0.068,0.721,0.128,0.834,0.285,0.321,0.009,0.695,0.879,0.867,0.64,0.265,0.361,0.411,0.76,0.495,0.778,0.733,0.713,0.763,0.008,0.436,0.753,0.334,0.54,0.849,0.804,0.642,0.295,0.466,0.564,0.405,0.045,0.19,0.64,0.087,0.753,0.559,0.855,0.695,0.795,0.846,0.766,0.116,0.417,0.707,0.604,0.326,0.498,0.942,0.066,0.416,0.658,0.337,0.873,0.626,0.087,0.623,0.646,0.862,0.412,0.461,0.429,0.807,0.674,0.731,0.834,0.325,0.387,0.554,0.586,0.118,0.932,0.701,0.731,0.936,0.549,0.622,0.72,0.682,0.247,0.474,0.097,0.239,0.283,0.607,0.729,0.711,0.536,0.619,0.727))

  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  #paste0( round(as.numeric(DC$table2[1:20,5]),s0),collapse = ",")
  expect(round(as.numeric(DC$table2[1:20,5]),s0),c(0.6,0.511,0.58,0.801,0.585,0.515,0.799,0.558,0.599,0.542,1,0.502,0.77,0.658,0.48,0.569,0.587,0.583,1,1))
})

