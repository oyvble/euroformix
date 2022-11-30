#Testing that the numerical calculation of loglik is correct
#rm(list=ls());library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results
nDone0=4
steptol0=1e-6
s0 = 3 #signif of checking

kit0 = "GlobalFiler" #name of selected kit
examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0("GF_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("GF_evid.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("GF_refs.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))#,threshT)
refs = sample_tableToList(tableReader(reffn))
kitinfo = getKit(kit0) #get kitinfo

#Force the reference to have a an allele not observed AND other observed spans the freq data
refs[["K48"]][["D16S539"]]$adata = c("7","13") #insert another allele to reference

#Set specific settings
AT = 75
pC = 0.05
lam = 0.01
fst = 0.01
#plotEPG2(samples,kit0,refs, AT=AT) #plot data
dat = prepareData(samples,refs,popFreq,threshT=AT,minF=NULL, normalize = TRUE) #obtain data to use for analysis

#Need marker specific settings: even if constant
init = function(x) setNames(rep(x,length(popFreq)),names(popFreq))
ATv = init(AT)
pCv = init(pC)
lamv =  init(lam)
fstv =  init(fst)

#dat$refData[["D16S539"]][["K48"]]
#The new
#dat$popFreq$D16S539  
#all(sapply(dat$popFreq,function(x) "99"%in%names(x)))
NOC = 4

test_that("check maximum likelihood Hp:", {
  cond = 1:NOC #condition on all contributors
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,xi=NULL,prC=pC,threshT=AT,fst=fst,lambda=lam,kit=kit0,xiFW=NULL, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, kit0)
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,2),c( 0.07,0.27,0.31,0.35,17554.68,0.17,0.99,0.08,0.01) )  #problem with 32-bit
  expect(round(mle$fit$loglik,s0), -1220.274) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  vals=c(0.152,0.128,0.226,0.546,0.551,0.742,0.568,0.649,0.966,0.73,0.967,0.577,0.42,0.604,0.53,0.69,0.473,0.259,0.589,0.378,0.908,0.952,0.93,0.01,0.964,0.109,0.335,0.034,0.292,0.026,0.534,0.85,0.834,0.788,0.63,0.834,0.817,0.224,0.777,0.92,0.371,0.245,0.109,0.027,0.063,0.117,0.386,0.824,0.874,0.859,0.883,0.987,0.368,0.589,0.801,0.69,0.68,0.918,0.839,0.585,0.111,0.367,0.263,0.636,0.558,0.559,0.629,0.767,0.267,0.648,0.416,0.503,0.425,0.542,0.173,0.882,0.472,0.498,0.117,0.022,0.023,0.261,0.382,0.445,0.521,0.795,0.701,0.572,0.752,0.088,0.411,0.753,0.461,0.924,0.68,0.254,0.969,0.629,0.841,0.507,0.963,0.497,0.905,0.223,0.417,0.954,0.54,0.654,0.147,0.052,0.341,0.031,0.105,0.319,0.8,0.653,0.706,0.387,0.644,0.504,0.739,0.58,0.727,0.686,0.162,0.556,0.614,0.031,0.013,0.077,0.024,0.01,0.068,0.497,0.203,0.11,0.013,0.301,0.474,0.715,0.373,0.119,0.152)
  #plot(sort(valid$ProbObs),sort(vals),ty="l")
  compareValid(valid$ProbObs,vals)
})

unknownIndHd = 4 #index of unknown

test_that("check maximum likelihood Hd (unrelated):", {
  
  #CALC loglik under Hd:
  cond = 1:NOC
  cond[unknownIndHd] = 0
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = unknownIndHd,xi=NULL,prC=pC,threshT=AT,fst=fst,lambda=lam,kit=kit0,xiFW=NULL, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv2 = getLogLiki(thhat, dat, NOC, cond,pCv,ATv,fstv,lamv, kit0)
  expect(logLikv,logLikv2)
  
  #Check param and loglik values:
  #paste0(round(thhat,s0),collapse = ",")
 # expect(round(thhat,2),c(0.072,0.269,0.312,0.347,17549.972,0.17,0.994,0.081,0.007) ) #problem with 32-bit
  expect(round(thhat,2),c(0.07,0.27,0.31,0.35,17549.97,0.17,0.99,0.08,0.01) ) #problem with 32-bit
  expect(round(mle$fit$loglik,s0), -1286.861) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  vals = c(0.152,0.134,0.27,0.272,0.268,0.627,0.565,0.65,0.996,0.721,0.964,0.574,0.417,0.61,0.49,0.682,0.469,0.26,0.58,0.379,0.965,0.976,0.924,0.003,0.047,0.107,0.76,0.027,0.27,0.024,0.53,0.846,0.933,0.827,0.625,0.83,0.814,0.221,0.772,0.917,0.368,0.246,0.106,0.026,0.067,0.107,0.384,0.82,0.871,0.85,0.876,0.995,0.361,0.586,0.801,0.688,0.676,0.926,0.835,0.585,0.111,0.371,0.259,0.628,0.559,0.622,0.662,0.761,0.271,0.642,0.415,0.505,0.422,0.541,0.176,0.892,0.559,0.493,0.13,0.034,0.001,0.194,0.046,0.443,0.522,0.774,0.683,0.563,0.837,0.086,0.409,0.748,0.724,0.963,0.673,0.248,0.965,0.624,0.826,0.502,0.993,0.491,0.897,0.222,0.77,0.933,0.526,0.952,0.14,0.039,0.21,0.03,0.055,0.31,0.796,0.649,0.718,0.38,0.637,0.5,0.735,0.577,0.726,0.685,0.159,0.55,0.609,0.031,0.002,0.005,0.007,0.006,0.068,0.027,0.028,0.072,0.01,0.299,0.469,0.67,0.382,0.116,0.187)
  #plot(sort(valid$ProbObs),sort(vals),ty="l")
  compareValid(valid$ProbObs,vals)
  
  #Calculate deconvolution (DC):
  #paste0(round( as.numeric(DC$table2[,ncol(DC$table2)-1]),s0),collapse = ",")
  DC = deconvolve(mle)
  expect(round( as.numeric(DC$table2[,ncol(DC$table2)-1]),s0),c(0.965,0.993,0.998,0.979,0.858,1,1,0.997,1,1,1,0.913,0.995,0.993,0.983,0.946,0.995,1,0.994,0.996,0.999))
})


