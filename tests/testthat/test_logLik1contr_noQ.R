#Testing that the numerical calculation of loglik is correct
#rm(list=ls());library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results

expect_approx = function(tol,x,y) { #helpfunction with specified tolerance
  expect_equal(as.numeric(x),as.numeric(y),tolerance = tol)
}
                
examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0("test_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("test_evid1.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("test_ref1.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))#,threshT)
refs = sample_tableToList(tableReader(reffn))

#Restrict outcome to those observed (no Q-allele)
markerDO = "D3S1358" #marker to perform "accident" on
popFreq[[markerDO]] = popFreq[[markerDO]][names( popFreq[[markerDO]])%in%samples[[1]][[markerDO]]$adata]
popFreq[[markerDO]] = popFreq[[markerDO]]/sum(popFreq[[markerDO]])

#Set specific settings for each dye:
AT = 50
pC = 0.05 
lam =0.01 
fst = 0.01 

#plotEPG2(samples,kit0,refs, AT=ATv) #plot data
dat = prepareData(samples,refs,popFreq,threshT=AT,minF=NULL, normalize = TRUE) #obtain data to use for analysis


test_that("check maximum likelihood Hp:", {
  
  #CALC loglik under Hp:
  mle = contLikMLE(nC=1,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=1,xi=NULL,prC=pC,nDone=2,threshT=AT,fst=fst,lambda=lam,xiFW=NULL,seed=seed0)
  #mle = contLikMLE(nC=1,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=1,xi=0,prC=pC,nDone=2,threshT=AT,fst=fst,lambda=lam,xiFW=0,seed=seed0)
  thhat=mle$fit$thetahat #obtain maximum likelihood estimates

  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #exctract params:  
  shape0 = 1/(thhat[2]^2) #get shape param
  scale0 = thhat[1]/shape0 #get scale param
  xiB = thhat[3] #backwards stutter prop
  xiF = thhat[4] #forward stutter prop
    
  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
#    loc = locsCheck[7]
    freq = dat$popFreq[[loc]]
    Aset <- Aset0 <- names(freq) #get allele outcome 
    if( "99"%in%Aset ) Aset0 = head(Aset,-1) #keep only non-dropouts
    Astutt = c(as.numeric(Aset0)-1,as.numeric(Aset0)+1) #get stutter alleles
    Aset = c(Aset,setdiff(Astutt,Aset0)) #insert potential stutters last
    
    #Assuming 1 known contribtion
    poi =  dat$refData[[loc]][[1]] #get alleles of POI
    nAG = matrix(0,nrow=length(Aset),ncol=1,dimnames = list(Aset,"poi")) #create allele counting matrix for unknown contributors
    nAG[,1] = table(factor(poi,level=Aset)) #insert contribution
    
    #Prepare PH data (same order as Aset):
    Y = rep(0,length(Aset))  #preparing PHs after extending allele vector
    Y[match(dat$samples[[1]][[loc]]$adata,Aset)] =dat$samples[[1]][[loc]]$hdata  #insert PHs
    
    mui <-  c(nAG[,1,drop=FALSE]*shape0) #expected contribution

    #Delegate stutterprod
    FWstutt = match(Aset,as.character(as.numeric(Aset)+1)) #alleleinds to stutter from
    BWstutt = match(Aset,as.character(as.numeric(Aset)-1)) #alleleinds to stutter from
    stuttB <- mui[BWstutt]*xiB #backward-stutter parts
    stuttF <- mui[FWstutt]*xiF #forward-stutter parts
    mui = mui*(1- (xiB+xiF))
    indBW = !is.na(stuttB)
    indFW = !is.na(stuttF)
    mui[indBW] = mui[indBW] + stuttB[indBW]
    mui[indFW] = mui[indFW] + stuttF[indFW]
    
    #Divide set into dropin, contr and dropout
    psiDI <- which( Y>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
    psiYmu <- which( Y>0 & mui>0 )  #contributing to model and observed PH
    psiDO <- which( Y==0 & mui>0 )  #dropout elem
    
    vali = 0
    if(length(psiYmu)>0) vali =  vali + sum(dgamma(Y[psiYmu],shape=mui[psiYmu],scale=scale0,log=TRUE))   # CONTRIBUTION TO PH
    if(length(psiDO)>0) vali = vali + sum(pgamma(AT,shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
    if(length(psiDI)>0) vali = vali + sum( dexp( Y[psiDI]-AT, rate=lam, log=TRUE) + log(pC*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
    if(length(psiDI)==0) vali = vali + log(1-pC) #in case of no dropin

    expect_equal(as.numeric(logLikv[loc]),as.numeric(vali))
  }      
  
  #Excpected joint values:
  expect_equal(as.numeric(thhat),c(755.98654446779,0.166446078342523,0.109723574534928,0.0444312063011257)) 
  expect_equal(mle$fit$loglik,-146.1482463) #check
  
  #Check cumulative probabilities (based on C++ code):
  valid = validMLEmodel(mle,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-6,valid$ProbObs,c( 0.665228142424458,0.957023625523975,0.658866144525858,0.232620963006383,0.252608695014166,0.290361189462231,0.0397127450090101,0.039107639177858,0.00325927921019866,0.544098850341052,0.857026750969566,0.536702691948382,0.314044189940371,0.463488753263063,0.420340957087939,0.293306084190538,0.465427584725022,0.781219649504299,0.327321752206193,0.484424213336163,0.314044189940371,0.456499944452378,0.758087854653409,0.992377585697513,0.210369337324636))
})

test_that("check maximum likelihood Hd (unrelated):", {

  #CALC loglik under Hp:
  mle = contLikMLE(nC=1,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=0,knownRef = 1,xi=NULL,prC=pC,nDone=2,threshT=AT,fst=fst,lambda=lam,xiFW=NULL,seed=seed0)
  thhat=mle$fit$thetahat #obtain maximum likelihood estimates

  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #exctract params:  
  shape0 = 1/(thhat[2]^2) #get shape param
  scale0 = thhat[1]/shape0 #get scale param
  xiB = thhat[3] #backwards stutter prop
  xiF = thhat[4] #forward stutter prop
  
  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
    #loc = locsCheck[1]
    
    freq = dat$popFreq[[loc]]
    Aset <- Aset0 <- names(freq) #get allele outcome 
    if( "99"%in%Aset ) Aset0 = head(Aset,-1) #keep only non-dropouts
    Astutt = c(as.numeric(Aset0)-1,as.numeric(Aset0)+1) #get stutter alleles
    Aset = c(Aset,setdiff(Astutt,Aset0)) #insert potential stutters last
    
    #Prepare PH data (same order as Aset):
    Y = rep(0,length(Aset))  #preparing PHs after extending allele vector
    Y[match(dat$samples[[1]][[loc]]$adata,Aset)] =dat$samples[[1]][[loc]]$hdata  #insert PHs
    
    #Get outcome of unknowns
    poi =  dat$refData[[loc]][[1]] #get alleles of POI (known under Hd)
    Glist = calcGjoint(freq=freq,nU=1,fst=fst,refK=poi)
    Gset = Glist$G #get allele out come of unknowns

    pEvid = 0 #obtain probabiliity of evidence (sum over all genotypes)
    for(gind in 1:nrow(Gset)) { #traverse all genorypes
      
      #Assuming 1 known contribtion
      nAG = matrix(0,nrow=length(Aset),ncol=1,dimnames = list(Aset,"unknown")) #create allele counting matrix for unknown contributors
      nAG[,1] = table(factor(Gset[gind,],level=Aset)) #insert contribution
      
      mui <-  c(nAG[,1,drop=FALSE]*shape0) #expected contribution
      
      #Delegate stutterprod
      FWstutt = match(Aset,as.character(as.numeric(Aset)+1)) #alleleinds to stutter from
      BWstutt = match(Aset,as.character(as.numeric(Aset)-1)) #alleleinds to stutter from
      stuttB <- mui[BWstutt]*xiB #backward-stutter parts
      stuttF <- mui[FWstutt]*xiF #forward-stutter parts
      indBW = !is.na(stuttB)
      indFW = !is.na(stuttF)
      indLooseStutt = indBW & indFW #index of alleles which are assumed to not loose stutter product
      mui[indLooseStutt] = mui[indLooseStutt]*(1- (xiB+xiF)) #loose stutter products
      mui[indBW] = mui[indBW] + stuttB[indBW]
      mui[indFW] = mui[indFW] + stuttF[indFW]
      
      #Divide set into dropin, contr and dropout
      psiDI <- which( Y>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
      psiYmu <- which( Y>0 & mui>0 )  #contributing to model and observed PH
      psiDO <- which( Y==0 & mui>0 )  #dropout elem
      
      vali = 0
      if(length(psiYmu)>0) vali =  vali + sum(dgamma(Y[psiYmu],shape=mui[psiYmu],scale=scale0,log=TRUE))   # CONTRIBUTION TO PH
      if(length(psiDO)>0) vali = vali + sum(pgamma(AT,shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
      if(length(psiDI)>0) vali = vali + sum( dexp( Y[psiDI]-AT, rate=lam, log=TRUE) + log(pC*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
      if(length(psiDI)==0) vali = vali + log(1-pC) #in case of no dropin
      pEvid = pEvid + exp(vali)*Glist$Gprob[gind] #ind
    }    
    expect_equal(as.numeric(logLikv[loc]),as.numeric(log(pEvid)))
  }      
  
  #Excpected joint values:
  expect_equal(as.numeric(thhat),c(755.980871446666,0.166443739251391,0.109723100635356,0.0444284773242095 )) #compare with manual derived
  expect_equal(mle$fit$loglik,-167.322130600091) #check
  
  valid = validMLEmodel(mle,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-6,valid$ProbObs,c(0.665241562961352,0.957028848955695,0.658900523025309,0.232626375327165,0.252613203252035,0.290371036705247,0.0392942242470234,0.0390659570221497,0.00325900163703341,0.544109498624105,0.857036616117881,0.541271432546718,0.3140512811963,0.463490038405816,0.420364120711596,0.293309240299272,0.465437286654174,0.781229032778341,0.327329145645313,0.484431326575639,0.314050807520921,0.456668446507917,0.75810952780816,0.992377554611051,0.2103798384693))

  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_equal(as.numeric(DC$table2[,2]),c(1,1,1,1,1,1,1))
  
})


