#Testing that the numerical calculation of loglik is correct (for SNP data)
#rm(list=ls());library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results
kit0 = NULL #No kit

expect_approx = function(tol,x,y) { #helpfunction with specified tolerance
  expect_equal(as.numeric(x),as.numeric(y),tolerance = tol)
}

examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0("SNP_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("SNP_evids.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("SNP_refs.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn)) #use only first replicate #,threshT)
refs = sample_tableToList(tableReader(reffn))

#Set common settings for all markers:
ATv = 60 #put 60 to remove alleles in one replicate but not the other
pCv = 0.05
lamv = 0.01
fstv = 0.02

#plotEPG2(samples,kit0,refs, AT=ATv) #plot data
dat = prepareData(samples,refs,popFreq,threshT=ATv,minF=NULL, normalize = TRUE) #obtain data to use for analysis

test_that("check maximum likelihood Hp:", {
  NOC = 2
  cond = c(1,2)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,xi=0,prC=pCv,nDone=2,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0,seed=seed0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates

  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #exctract params:  
  mx=thhat[1:NOC]
  shape0 = 1/(thhat[NOC+2]^2) #get shape param
  scale0 = thhat[NOC+1]/shape0 #get scale param

  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
#loc=locsCheck[1]
    freq = dat$popFreq[[loc]]
    Aset = names(freq) #get allele outcome 
    
    condRef =  dat$refData[[loc]][cond] #get alleles of conds
    if(any(sapply(condRef,is.null))) next #Skip if any references are missing
    nAG = matrix(0,nrow=length(Aset),ncol=NOC,dimnames = list(Aset)) #create allele counting matrix for unknown contributors
    for(ind in 1:length(condRef)) nAG[,ind] = table(factor(condRef[[ind]],level=Aset)) #insert contribution
    
    #Prepare PH data (same order as Aset):
    Ylist = list() #list per samples
    for(sample in names(dat$samples)) {
      Ylist[[sample]] = rep(0,length(Aset))  #preparing PHs after extending allele vector
      Ylist[[sample]][match(dat$samples[[sample]][[loc]]$adata,Aset)] = dat$samples[[sample]][[loc]]$hdata  #insert PHs
    }
    mui <-  c(nAG[,drop=FALSE]%*%mx)*shape0 #expected contribution

    vali = 0
    for(sample in names(dat$samples)) {
      psiDI <- which( Ylist[[sample]]>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
      psiYmu <- which( Ylist[[sample]]>0 & mui>0 )  #contributing to model and observed PH
      psiDO <- which( Ylist[[sample]]==0 & mui>0 )  #dropout elem
      
      if(length(psiYmu)>0) vali =  vali + sum(dgamma(Ylist[[sample]][psiYmu],shape=mui[psiYmu],scale=scale0,log=TRUE))   # CONTRIBUTION TO PH
      if(length(psiDO)>0) vali = vali + sum(pgamma(ATv,shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
      if(length(psiDI)>0) vali = vali + sum( dexp( Ylist[[sample]][psiDI]-ATv, rate=lamv, log=TRUE) + log(pCv*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
      if(length(psiDI)==0) vali = vali + log(1-pCv) #in case of no dropin
    }
    expect_equal(as.numeric(logLikv[loc]),as.numeric(vali))
  }      
    
  #Excpected joint values:
  expect_approx(1e-3,thhat,c(0.784287340610599,0.215712659389401,797.885130883847,0.590903414911386))
  expect_equal(mle$fit$loglik,-3980.05359609996) #check
    
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE) #check 20 first only
  expect_approx(1e-3,valid$ProbObs[1:20],c(0.549096125474773,0.535092179262372,0.826887844117223,0.882357942413679,0.23164476702976,0.314145305928015,0.352790002895515,0.430217208963195,0.987545264943279,0.98878315417829,0.604512633654525,0.633268630267128,0.839451971275505,0.780315767107799,0.799021640515115,0.799873487917275,0.548019621819006,0.815724995230886,0.267254047924345,0.589649584713638))
  
})

test_that("check maximum likelihood Hd (unrelated):", {
  NOC = 2
  #CALC loglik under Hp:
  cond = c(1,0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = 2,xi=0,prC=pCv,nDone=2,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0,seed=seed0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #exctract params:  
  mx=thhat[1:NOC]
  shape0 = 1/(thhat[NOC+2]^2) #get shape param
  scale0 = thhat[NOC+1]/shape0 #get scale param

  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
#loc=locsCheck[2]
    freq = dat$popFreq[[loc]]
    Aset = names(freq) #get allele outcome 

    #Assuming 1 known contribtion
    condRef =  dat$refData[[loc]][cond] #get alleles of conds
    if(any(sapply(condRef,is.null))) next #Skip if any references are missing
    
    nAG = matrix(0,nrow=length(Aset),ncol=NOC,dimnames = list(Aset)) #create allele counting matrix for unknown contributors
    for(ind in 1:length(condRef)) nAG[,ind] = table(factor(condRef[[ind]],level=Aset)) #insert contribution
    
    #Prepare PH data (same order as Aset):
    Ylist = list() #list per samples
    for(sample in names(dat$samples)) {
      Ylist[[sample]] = rep(0,length(Aset))  #preparing PHs after extending allele vector
      Ylist[[sample]][match(dat$samples[[sample]][[loc]]$adata,Aset)] = dat$samples[[sample]][[loc]]$hdata  #insert PHs
    }

    nAG0 = nAG #temporary store
    Glist = calcGjoint(freq=freq,nU=1,fst=fstv,refK=unlist( dat$refData[[loc]]))
    Gset = Glist$G #get allele out come of unknowns

    pEvid = 0 #obtain probabiliity of evidence (sum over all genotypes)
    for(gind in 1:nrow(Gset)) { #traverse all genorypes
      nAG = nAG0
      nAG[,2] = table(factor(Gset[gind,],level=Aset)) #insert contribution
      mui <-  c(nAG[,drop=FALSE]%*%mx)*shape0 #expected contribution

      vali = 0 
      for(sample in names(dat$samples)) { #traversing all samples
        #Divide set into dropin, contr and dropout
        psiDI <- which( Ylist[[sample]]>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
        psiYmu <- which( Ylist[[sample]]>0 & mui>0 )  #contributing to model and observed PH
        psiDO <- which( Ylist[[sample]]==0 & mui>0 )  #dropout elem
        
        if(length(psiYmu)>0) vali =  vali + sum(dgamma(Ylist[[sample]][psiYmu],shape=mui[psiYmu],scale=scale0,log=TRUE))   # CONTRIBUTION TO PH
        if(length(psiDO)>0) vali = vali + sum(pgamma(ATv,shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
        if(length(psiDI)>0) vali = vali + sum( dexp( Ylist[[sample]][psiDI]-ATv, rate=lamv, log=TRUE) + log(pCv*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
        if(length(psiDI)==0) vali = vali + log(1-pCv) #in case of no dropin
      }
      pEvid = pEvid + exp( vali + log(Glist$Gprob[gind])) #sum up
    } #end for each genotype
    expect_equal(as.numeric(logLikv[loc]),as.numeric(log(pEvid)))
  }      
  
  #Excpected joint values:
  expect_approx(1e-3,thhat,c(0.666100819168908,0.333899180831092,815.87511749016,0.5130965438557))
  expect_equal(mle$fit$loglik,-3807.17711247693) #check
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-3,valid$ProbObs[1:20],c(0.726790292608432,0.712919160966207,0.699109750829867,0.796738750078469,0.242610970624618,0.339553372676808,0.203522487453136,0.272233089671143,0.961761695281427,0.965808886347945,0.787992312637902,0.811710455551996,0.898829992718985,0.847228444081038,0.70813046277954,0.709339415292811,0.391519590292488,0.75919484571863,0.301499156240431,0.688403837316524))

  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_approx(1e-3,DC$table2[1:20,5],c(0.7031,0.9075,0.7596,0.904,0.6718,0.8387,0.887,0.6318,0.7209,0.7599,0.7601,0.8476,0.9067,0.988,0.6098,0.9391,0.5548,0.7928,0.6857,0.8198))
  
})

