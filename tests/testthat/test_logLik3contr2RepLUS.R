#Testing that the numerical calculation of loglik is correct (for SNP data)
#library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results
nDone0=2
steptol0=1e-6
s0 = 3 #signif of checking

kit0 = "ForenSeq" #No kit
examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0(kit0,"_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0(kit0,"_evids.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0(kit0,"_refs.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))[1] #use only first replicate #,threshT)
refData = sample_tableToList(tableReader(reffn))

#Set common settings for all markers:
AT = 30 
pC = 0.05
lam = 0.01
fst = 0.01

init = function(x) setNames(rep(x,length(popFreq)),names(popFreq))
ATv = init(AT)
pCv = init(pC)
lamv =  init(lam)
fstv =  init(fst)

#plotMPS2(samples,refData, AT=ATv) #plot data
dat = prepareData(samples,refData,popFreq,threshT=ATv,minF=NULL, normalize = TRUE) #obtain data to use for analysis
#plotMPS2(dat$samples,dat$refData)
NOC = 3

test_that("check maximum likelihood Hp: No Stutter", {
  cond = c(1,2,3,0)
# mle = contLikMLE(nC=NOC,samples=samples,popFreq=popFreq,refData=refData,condOrder=cond,xi=0,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,xi=0,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c(0.695,0.219,0.087,1082.381,0.46,0.781) )
  expect(round(mle$fit$loglik,s0),-946.952) 
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  vals = c(0.165,0.024,0.048,0.298,0.087,0.212,0.532,0.206,0.572,0.901,0.431,0.483,0.629,0.467,0.542,0.298,0.384,0.149,0.106,0.03,0.768,0.46,0.793,0.956,0.959,0.527,0.042,0.223,0.5,0.087,0.699,0.703,0.646,0.677,0.004,0.056,0.052,0.06,0.007,0.117,0.772,0.813,0.098,0.213,0.649,0.289,0.327,0.43,0.194,0.189,0.765,0.786,0.647,0.745,0.316,0.896,0.536,0.698,0.281,0.944,0.343,0.419,0.235,0.246,0.336,0.947,0.998,0.104,0.349,0.068,0.028,0.136,0.589,0.681,0.154,0.554,0.365,0.039,0.527,0.534,0.304,0.403,0.376,0.36,0.199,0.268,0.645,0.633,0.784,0.552,0.55,0.711,0.044,0.229,0.045,0.341,0.423,0.546,0.606,0.402,0.181,0.611,0.188,0.354,0.213,0.422,0.858,0.997,0.997,0.509,0.148,0.268,0.425,0.49,0.348,0.073,0.143,0.667,0.274,0.897,0.895,0.218)
  #compareValid(valid$ProbObs,vals)
  #plot(sort(valid$ProbObs),sort(vals),ty="l")
})

test_that("check maximum likelihood Hd (unrelated): No Stutter", {
  cond = c(1,0,2,0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = 2,xi=0,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c(0.713,0.105,0.182,1107.704,0.462,0.778) )
  expect(round(mle$fit$loglik,s0),-966.185) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  vals=c(0.146,0.013,0.031,0.294,0.106,0.145,0.354,0.185,0.625,0.876,0.699,0.389,0.595,0.617,0.443,0.413,0.353,0.134,0.095,0.025,0.805,0.421,0.762,0.945,0.981,0.456,0.035,0.275,0.561,0.074,0.486,0.799,0.79,0.626,0.009,0.046,0.029,0.028,0.012,0.135,0.806,0.896,0.087,0.158,0.62,0.263,0.193,0.499,0.172,0.142,0.592,0.859,0.77,0.713,0.225,0.876,0.808,0.787,0.197,0.931,0.215,0.364,0.594,0.22,0.191,0.965,0.999,0.092,0.272,0.051,0.037,0.193,0.449,0.798,0.141,0.524,0.407,0.031,0.488,0.497,0.545,0.458,0.345,0.344,0.216,0.291,0.606,0.578,0.754,0.609,0.61,0.477,0.067,0.274,0.04,0.391,0.285,0.702,0.566,0.37,0.127,0.578,0.168,0.834,0.152,0.392,0.636,0.996,0.999,0.42,0.126,0.337,0.494,0.46,0.308,0.116,0.073,0.708,0.205,0.965,0.87,0.196)
#  compareValid(valid$ProbObs,vals)
  #  plot(sort(valid$ProbObs),sort(vals),ty="l")
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle)
#paste0( round(as.numeric(DC$table2[1:20,8]),s0),collapse = ",")
  expect(round(as.numeric(DC$table2[1:20,8]),s0),c(0.46,0.514,0.885,0.246,0.893,0.94,0.421,0.342,0.393,0.916,0.364,0.563,0.41,0.937,0.675,0.501,0.905,0.877,0.912,0.905))
})


test_that("check maximum likelihood Hp: BWS", {
  cond = c(1,2,3,0)
  # mle = contLikMLE(nC=NOC,samples=samples,popFreq=popFreq,refData=refData,condOrder=cond,xi=0,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,xi=NULL,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c(0.744,0.191,0.065,1122.011,0.408,0.787,0.109) )
  expect(round(mle$fit$loglik,s0),-837.321) 
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  vals = c(0.104,0.017,0.031,0.36,0.071,0.187,0.368,0.244,0.647,0.923,0.532,0.473,0.635,0.539,0.527,0.353,0.363,0.178,0.127,0.068,0.729,0.524,0.719,0.969,0.977,0.551,0.044,0.272,0.249,0.105,0.637,0.723,0.626,0.683,0.003,0.07,0.043,0.042,0.006,0.068,0.704,0.85,0.117,0.251,0.652,0.342,0.123,0.437,0.234,0.288,0.653,0.744,0.719,0.763,0.319,0.904,0.612,0.773,0.285,0.96,0.342,0.433,0.2,0.216,0.421,0.925,0.999,0.125,0.224,0.055,0.035,0.163,0.567,0.651,0.164,0.544,0.417,0.058,0.479,0.529,0.385,0.467,0.364,0.327,0.224,0.302,0.559,0.582,0.785,0.625,0.613,0.441,0.041,0.275,0.054,0.4,0.398,0.575,0.536,0.381,0.212,0.604,0.221,0.473,0.247,0.402,0.799,0.997,0.999,0.613,0.186,0.312,0.347,0.474,0.346,0.039,0.111,0.735,0.382,0.849,0.912,0.254)
  #compareValid(valid$ProbObs,vals)
  #plot(sort(valid$ProbObs),sort(vals),ty="l")
})

test_that("check maximum likelihood Hd (unrelated): BWS", {
  cond = c(1,0,2,0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = 2,xi=NULL,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c(0.787,0.071,0.141,1128.41,0.435,0.774,0.127) )
  expect(round(mle$fit$loglik,s0),-891.753) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  vals=c(0.157,0.014,0.027,0.404,0.107,0.149,0.305,0.222,0.68,0.893,0.686,0.426,0.606,0.646,0.448,0.438,0.352,0.168,0.119,0.043,0.82,0.446,0.68,0.952,0.985,0.454,0.028,0.338,0.438,0.09,0.556,0.789,0.738,0.628,0.011,0.077,0.038,0.019,0.012,0.186,0.813,0.867,0.088,0.226,0.63,0.341,0.086,0.509,0.204,0.179,0.605,0.814,0.797,0.721,0.279,0.875,0.771,0.793,0.251,0.942,0.235,0.335,0.577,0.192,0.334,0.918,0.999,0.118,0.18,0.055,0.051,0.221,0.522,0.723,0.154,0.519,0.47,0.05,0.448,0.49,0.564,0.527,0.378,0.347,0.251,0.338,0.479,0.523,0.753,0.692,0.738,0.339,0.06,0.317,0.051,0.477,0.344,0.681,0.442,0.355,0.187,0.566,0.183,0.776,0.189,0.391,0.736,0.994,1,0.526,0.13,0.368,0.563,0.443,0.334,0.079,0.084,0.859,0.247,0.896,0.88,0.209)
  #  compareValid(valid$ProbObs,vals)
  #  plot(sort(valid$ProbObs),sort(vals),ty="l")
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle)
  #paste0( round(as.numeric(DC$table2[1:20,8]),s0),collapse = ",")
  expect(round(as.numeric(DC$table2[1:20,8]),s0),c(0.314,0.433,0.577,0.801,0.123,0.389,0.157,0.355,0.242,0.262,0.37,0.266,0.554,0.253,0.336,0.957,0.338,0.495,0.881,0.659))
})


