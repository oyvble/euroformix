#TESTING IF CALL TO EFMengine is working.
#require(euroformix);library(testthat);sessionInfo()
kit = "testkit" #name of selected kit
examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0("test_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("test_evid1.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("test_ref1.csv"),sep=.Platform$file.sep)


#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))#,threshT)
refData = sample_tableToList(tableReader(reffn))

#settings  
AT = 50; pC=0.05; lambda=0.01; fst=0#.01;
#Hypotheses
condOrder=NULL;knownRef=NULL

#mod=modMat[1,]
calcLogLik = function(mod,DEG=TRUE) {
#  nC=2;BWS=TRUE;FWS=FALSE
#  nC=3;BWS=TRUE;FWS=FALSE
  #  nC=1;BWS=TRUE;FWS=FALSE
  nC=mod$nC
  BWS=mod$BWS
  FWS=mod$FWS
  param = rep(1,nC)
  param = param/sum(param)
  param = c(param, 1000,0.1,1,0.1,0.01)
  if(!BWS) param[nC+4] = 0
  if(!FWS) param[nC+5] = 0
  c = prepareC(nC,samples,popFreq, refData, condOrder, knownRef, kit,BWS,FWS,AT,pC,lambda,fst,NULL,NULL,NULL,TRUE, FALSE)
  obj = prepareCobj(c) #use wrapper function to obtain C++ pointer
  loglik = obj$calcGenoWeightsMax(as.numeric(param) ) #Calculate
  logliki = obj$logliki()
  obj$close() #free memory
  return(list(loglik,logliki))
}

nCvec = 1:4
BWSvec <- FWSvec <- c(FALSE,TRUE)
modMat = expand.grid(nC=nCvec,BWS=BWSvec,FWS=FWSvec)
modMat = modMat[!(modMat$FWS & !modMat$BWS),]
nMods = nrow(modMat)

logLikTrue = c(-353.5401, -373.5067, -370.1198, -350.5569
            , -253.415, -273.3755, -289.5421, -303.317
            , -244.9223, -264.8718, -280.718, -296.1194)


test_that("All sanity examples are correct", {
  for(i in seq_len(nMods)) {
    logLik = calcLogLik(modMat[i,])[[1]]
    expect(round(logLik,4),logLikTrue[i])
  }
})


