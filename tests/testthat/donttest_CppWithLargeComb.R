#TESTING IF CALL TO EFMengine is working with a large number of genotype combinations
#NOTE: THIS TEST WILL REQUIRE A COMPUTER WITH 32GB RAM
require(euroformix);library(testthat);sessionInfo()
nC = 5 #large 
nA = 11 #number of alleles
samples = list(evid=list(M1=list(adata=c(1:nA),hdata=rep(1000,nA))))
freq = rep(0.01,nA)
freq = setNames(c(freq,1-sum(freq)),c(1:nA,"99"))
popFreq = list(M1=freq)
AT = 50; pC=0.05; lambda=0.01; fst=0;

BWS=FALSE;FWS=FALSE;DEG=FALSE
param = rep(1,nC)
param = param/sum(param)
param = c(param, 1000,0.1,1,0.1,0.01)
if(!BWS) param[nC+4] = 0
if(!FWS) param[nC+5] = 0
c = prepareC(nC,samples,popFreq, NULL, NULL, NULL, NULL, BWS,FWS,AT,pC,lambda,fst,NULL,NULL,NULL,TRUE, FALSE)
obj = prepareCobj(c) #use wrapper function to obtain C++ pointer
#obj$close() #free memory


#refData=NULL;condOrder=NULL;knownRef=NULL;kit=NULL;knownRel=NULL;ibd=NULL;minF=NULL;normalize=TRUE;adjFragQallele=FALSE

