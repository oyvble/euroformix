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

  #PREPARE FEEDING DATA
  mod = Rcpp::Module( "mod",PACKAGE="euroformix" ) #load module
  obj = methods::new(mod$ExposedClass) #create object of class
  #Step 1: insert data to exposed class (filldata)
  obj$filldata(c$nStutterModels,c$nMarkers,c$nRepMarkers,c$nAlleles,c$startIndMarker_nAlleles,c$startIndMarker_nAllelesReps,c$peaks,c$freqs,c$dropinWeight, c$nTyped, c$maTyped, c$basepair,
               c$BWfrom, c$FWfrom, c$BWto, c$FWto, c$nPotStutters, c$startIndMarker_nAllelesTot, c$QalleleIndex, c$dropinProb, c$fst, c$AT, c$NOK, c$knownGind, c$relGind, c$ibd, as.integer(0)) 
  #Step 2: Indexing large matrix (doIndex)
  obj$prepare(as.integer(nC))
  #sum((nA+nPS)*(nA*(nA+1)/2)^nC) #number of evaluations
  loglik = obj$loglik(as.numeric(param) ) #Calculate
  expect(round(loglik,4),-39.3401)
  #FALSE=-32.99025 (caused by overflow) 
})
  