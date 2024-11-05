#Testing that the numerical calculation of loglik is correct
#library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results
nDone0=5
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
# paste0(round(thhat,s0),collapse = ",")
#  expect(round(thhat,s0),c(0.762,0.238,1904.481,0.285,0.686,0.057) )
  expect(round(thhat,s0),c(0.765,0.235,1870.515,0.288,0.696,0.057) )
  expect(round(mle$fit$loglik,s0),-400.064) 
  
  #CHECK Cumulative probs    
# paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  compareValid(valid$ProbObs,c(0.271,0.13, 0.199,0.837,0.699,0.572,0.331,0.279,0.682,0.867,0.269,0.369,0.184,0.048,0.203,0.432,0.795,0.7,0.051,0.938,0.429,0.2,0.361,0.294,0.271,0.514,0.531,0.442,0.002,0.014,0.504,0.607,0.247,0.186,0.508,0.56,0.472,0.4,0.293,0.335,0.998,0.779,1,0.314,0.725,0.908,0.173,0.716,0.711,0.78,0.564,0.324,0.452,0.621))
})

test_that("check maximum likelihood Hd:", {
  cond = c(1,0)
  knownRef = 2
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef=knownRef,xi=NULL,prC=pCv,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0, seed=seed0,steptol=steptol0,nDone=nDone0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS
  logLikv = logLiki(mle) #obtain per marker resutls
  expect(sum(logLikv),mle$fit$loglik)
  
  #Check param and loglik values:
# paste0(round(thhat,s0),collapse = ",")
  expect(round(thhat,s0),c(0.771,0.229,1874.725,0.312,0.692,0.061))
  expect(round(mle$fit$loglik,s0),-426.993) 
  
  #CHECK Cumulative probs    
  #paste0(round(valid$ProbObs,s0),collapse = ",")
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
  compareValid(valid$ProbObs,c(0.4,0.059,0.113,0.839,0.79,0.597,0.339,0.285,0.644,0.845,0.213,0.315,0.445,0.008,0.114,0.435,0.83,0.699,0.045,0.921,0.432,0.221,0.368,0.313,0.28,0.624,0.454,0.406,0.015,0.002,0.452,0.566,0.401,0.196,0.545,0.53,0.471,0.326,0.429,0.299,0.995,0.963,0.999,0.402,0.71,0.885,0.143,0.693,0.655,0.887,0.56,0.341,0.433,0.615))
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle)
  expect_equal(as.numeric(DC$table2[1,5]),c(0.6499))
})

