#Testing of whether the numerical calculation of loglik (maximum likelihood approach) is correct
#Different kinds of models are tested

#rm(list=ls());library(euroformix);library(testthat)
kit0 = "testkit" #name of selected kit
s0 = 3 #signif
#Helpfunction for checking model after model fit
checkModel0Unknown = function(mle) {
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = mle$logliki
  thhat=mle$pDhatContr #obtain drop-out probs per contributor
  cond = mle$model$condOrder
  NOC = mle$model$nC
  
  #Obtain data  
  popFreq = mle$model$popFreq
  refData = mle$model$refData
  sample = mle$model$samples
  
  expect_equal(sum(logLikv),mle$loglik) #checking joint loglik
  #exctract params:  
  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
    #loc=locsCheck[1]
    freq = popFreq[[loc]]
    Aset = names(freq) #get allele outcome 
    Aset0 = head(Aset,-1) #keep non-dropouts
    
    #Assuming 2 known contribtions
    condRef =  refData[[loc]][cond] #get alleles of conds
    nAG = matrix(0,nrow=length(Aset),ncol=NOC,dimnames = list(Aset)) #create allele counting matrix for unknown contributors
    for(ind in 1:length(condRef)) nAG[,ind] = table(factor(condRef[[ind]],level=Aset)) #insert contribution
    
    #temporary caluclation
    pDprodTmp = 1 #calculte temporary
    for(k in 1:NOC) { #traverse each contributor
      pDprodTmp = pDprodTmp*thhat[k]^nAG[,k]
    } 
    mui = 1-pDprodTmp #rowSums(nAG) #obtain per-allele contribution
    
    #calculate likelihood for observations
    vali = 0 
    for(sample in names(samples)) { #for each replicates
      #Divide set into dropin, contr and dropout
      #      sample=names(samples)[1]
      yvec = rep(0,length(Aset))  #preparing PHs after extending allele vector
      yvec[match(samples[[sample]][[loc]]$adata,Aset)] = samples[[sample]][[loc]]$hdata  #insert PHs
      
      psiDI <- which( yvec>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
      psiYmu <- which( yvec>0 & mui>0 )  #contributing to model and observed PH
      psiDO <- which( yvec==0 & mui>0 )  #dropout elem
      
      if(length(psiYmu)>0) vali =  vali + sum(log(1-pDprodTmp[psiYmu]))   # CONTRIBUTION TO PH
      if(length(psiDO)>0) vali = vali + sum(log(pDprodTmp[psiDO])) #CONTR TO DROPOUT
      if(length(psiDI)>0) {
        vali = vali + sum(log(pCv[loc]*freq[psiDI])) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
      } else {
        vali = vali + log(1-pCv[loc]) #in case of no dropin
      }
    }
    expect_equal(as.numeric(logLikv[loc]),as.numeric(vali))
  }      
}

#Helpfunction for checking model after model fit (iterating over 1 genotype set)
checkModel1Unknown = function(mle,dat) {
 
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = mle$logliki
  thhat = mle$pDhatContr #obtain drop-out probs per contributor
  cond = mle$model$condOrder
  NOC = mle$model$nC
  
  #Obtain data  
  popFreq = mle$model$popFreq
  refData = mle$model$refData
  sample = mle$model$samples
  
  expect_equal(sum(logLikv),mle$loglik) #checking joint loglik
  
  #exctract params:  
  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
    #loc=locsCheck[1]
    freq = popFreq[[loc]]
    Aset = names(freq) #get allele outcome 
    Aset0 = head(Aset,-1) #keep non-dropouts
    
    #Assuming 2 known contribtions
    condRef =  refData[[loc]][cond] #get alleles of conds
    nAG = matrix(0,nrow=length(Aset),ncol=NOC,dimnames = list(Aset)) #create allele counting matrix for unknown contributors
    for(ind in 1:length(condRef)) nAG[,ind] = table(factor(condRef[[ind]],level=Aset)) #insert contribution
    
    #Obtain all combinations for unknown
    nAG0 = nAG #temporary store
    Glist = calcGjoint(freq=freq,nU=1,fst=fstv[loc],refK=unlist( refData[[loc]]))
    Gset = Glist$G #get allele out come of unknowns
    
    #Traversing all combinations for unknown
    pEvid = 0 #obtain probabiliity of evidence (sum over all genotypes)
    for(gind in 1:nrow(Gset)) { #traverse all genorypes
      nAG = nAG0
      nAG[,2] = table(factor(Gset[gind,],level=Aset)) #insert contribution
      
      #temporary caluclation
      pDprodTmp = 1 #calculte temporary
      for(k in 1:NOC) { #traverse each contributor
        pDprodTmp = pDprodTmp*thhat[k]^nAG[,k]
      } 
      mui = 1-pDprodTmp #rowSums(nAG) #obtain per-allele contribution
      
      #calculate likelihood for observations
      vali = 0 
      for(sample in names(samples)) { #for each replicates
        #Divide set into dropin, contr and dropout
        #      sample=names(samples)[1]
        yvec = rep(0,length(Aset))  #preparing PHs after extending allele vector
        yvec[match(samples[[sample]][[loc]]$adata,Aset)] = samples[[sample]][[loc]]$hdata  #insert PHs
        
        psiDI <- which( yvec>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
        psiYmu <- which( yvec>0 & mui>0 )  #contributing to model and observed PH
        psiDO <- which( yvec==0 & mui>0 )  #dropout elem
        
        if(length(psiYmu)>0) vali =  vali + sum(log(1-pDprodTmp[psiYmu]))   # CONTRIBUTION TO PH
        if(length(psiDO)>0) vali = vali + sum(log(pDprodTmp[psiDO])) #CONTR TO DROPOUT
        if(length(psiDI)>0) {
          vali = vali + sum(log(pCv[loc]*freq[psiDI])) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
        } else {
          vali = vali + log(1-pCv[loc]) #in case of no dropin
        }
      }
      pEvid = pEvid + exp( vali + log(Glist$Gprob[gind])) #sum up
    } #end for each genotype
    expect_equal(as.numeric(logLikv[loc]),as.numeric(log(pEvid)))
  }
}

examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0("test_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("test_evids.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("test_refs.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))#,threshT)
refs = sample_tableToList(tableReader(reffn))
kitinfo = getKit(kit0) #get kitinfo

#Set specific settings for each dye:
cols = c("blue","yellow")
ATdye = c(50,70)
pCdye = c(0.0133 , 0.0097)
lamdye = c(0.025, 0.034)
fstdye = c(0.01,0.02)
maxiter=5 #max number of iter

#Allign settings wrt color-marker info  
kitdyes = getKit(kit0,"COLOR") #get kitinfo from selected kit
indmatch = match(kitdyes$Color,cols)
ATv = ATdye[indmatch] 
pCv = pCdye[indmatch]
fstv = fstdye[indmatch]
names(ATv) <- names(pCv) <- names(fstv)  <- toupper(kitdyes$Marker) #don't need to be same order as of popFreq

#plotEPG2(samples,kit0,refs, AT=ATv) #plot data
dat = prepareData(samples,refs,popFreq,threshT=ATv,minF=NULL, normalize = TRUE) #obtain data to use for analysis


test_that("Hyp 1: R1+R2 (common drop-out parameter):", {
  NOC = 2
  cond = c(1,2)
  #steptol=1e-6; prDv0=c(0.1,0.35,0.7); knownRef=NULL; nC=NOC;samples=dat$samples;popFreq=dat$popFreq;refData=dat$refData;condOrder=cond;prC=pCv;fst=fstv;maxIter=maxiter;prDcontr = c(0,NA);prDcommon = NULL
  mle = qualLikMLE(nC=2,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,prC=pCv,fst=fstv,maxIter=maxiter)

  expect(round(mle$loglik,s0),c( -95.117) ) #compare with manual derived
  expect(round(mle$pDhatContr,s0),c( 0.073, 0.073) ) #compare with manual derived
  
  checkModel0Unknown(mle) #check per-marker numeric
})


test_that("Hyp 2: R1+R2 (cond zero drop-out prob for R1):", {
  NOC = 2
  cond = c(1,2)
#steptol=1e-6; prDv0=c(0.1,0.35,0.7); knownRef=NULL; nC=NOC;samples=dat$samples;popFreq=dat$popFreq;refData=dat$refData;condOrder=cond;prC=pCv;fst=fstv;maxIter=maxiter;prDcontr = c(0,NA);prDcommon = NULL
  mle = qualLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,prC=pCv,fst=fstv,maxIter=maxiter,prDcontr = c(0,NA),prDcommon = NULL)

  expect(round(mle$pDhatContr,s0),c( 0.000, 0.146) ) #compare with manual derived
  expect(round(mle$loglik,s0),c( -92.916) ) #compare with manual derived
  checkModel0Unknown(mle) #check per-marker numeric
})


test_that("Hyp 3: R1+R2 (two drop-out prob params):", {
  NOC = 2
  cond = c(1,2)
  mle = qualLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,prC=pCv,fst=fstv,maxIter=100,prDcommon = c(1,2))
  
  expect(round(mle$pDhatContr,s0),c( 0.000, 0.146) ) #compare with manual derived
  expect(round(mle$loglik,s0),c( -92.916) ) #compare with manual derived
  checkModel0Unknown(mle) #check per-marker numeric
})



test_that("Hyp 4: R1+1U (common drop-out parameter):", {
  NOC = 2
  cond = c(1,0)
  mle = qualLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef=2,prC=pCv,fst=fstv,maxIter=maxiter)
  
  expect(round(mle$loglik,s0),c( -77.659) ) #compare with manual derived
  expect(round(mle$pDhatContr,s0) ,c(   0.038, 0.038) ) #compare with manual derived
  checkModel1Unknown(mle) #check per-marker numeric
})


test_that("Hyp 5: R1+1U (cond zero drop-out prob for R1):", {
  NOC = 2
  cond = c(1,0)
  mle = qualLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef=2,prC=pCv,fst=fstv,maxIter=maxiter,prDcontr = c(0,NA))
  
  expect(round(mle$loglik,s0),c( -76.45) ) #compare with manual derived
  expect(round(mle$pDhatContr,s0) ,c(   0, 0.075) ) #compare with manual derived
  checkModel1Unknown(mle) #check per-marker numeric
})

test_that("Hyp 6: R1+1U (two drop-out prob params):", {
  NOC = 2
  cond = c(1,0)
  mle = qualLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef=2,prC=pCv,fst=fstv,maxIter=100,prDcommon = c(1,2))
  
  expect(round(mle$loglik,s0),c( -76.45) ) #compare with manual derived
  expect(round(mle$pDhatContr,s0) ,c(  0, 0.075) ) #compare with manual derived
  checkModel1Unknown(mle) #check per-marker numeric
})




