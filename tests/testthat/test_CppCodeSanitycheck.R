#TESTING IF CALL TO EFMengine is working.
#require(euroformix);library(testthat);sessionInfo()
samples = list(evid=list(M1=list(adata=c(10,11,9),hdata=c(637,603,67))))
freq = c(0.055,0.300,0.144)
freq = setNames(c(freq,1-sum(freq)),c(samples$evid$M1$adata,"99"))
popFreq = list(M1=freq)

#settings 
AT = 50; pC=0.05; lambda=0.01; fst=0;
#Hypotheses

test_that("Cpp code is OK", {
  #calcLogLik = function(param,nC,BWS,FWS) {
  #Obtain C-object:
  nC=3;BWS=TRUE;FWS=FALSE;
  param = rep(1,nC)
  param = param/sum(param)
  param = c(param, 1000,0.1,1,0.1,0.01)
  if(!BWS) param[nC+4] = 0
  if(!FWS) param[nC+5] = 0
  c = prepareC(nC,samples,popFreq, NULL, NULL, NULL, NULL,BWS,FWS,AT,pC,lambda,fst,NULL,NULL,NULL,TRUE, FALSE)
  obj = prepareCobj(c) #use wrapper function to obtain C++ pointer
  loglik = obj$calcGenoWeightsMax(as.numeric(param) ) #Calculate
  expect(round(loglik,4),-39.3401)
})
  