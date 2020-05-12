#Testing that the numerical calculation of loglik is correct
#rm(list=ls());library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results
kit0 = "testkit" #name of selected kit

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
kitinfo = getKit(kit0) #get kitinfo

#Set specific settings for each dye:
cols = c("blue","yellow")
ATdye = c(50,70)
pCdye = c(0.0133 , 0.0097)
lamdye = c(0.025, 0.034)
fstdye = c(0.01,0.02)

#Allign settings wrt color-marker info  
kitdyes = getKit(kit0,"COLOR") #get kitinfo from selected kit
indmatch = match(kitdyes$Color,cols)
ATv = ATdye[indmatch] 
pCv = pCdye[indmatch]
lamv = lamdye[indmatch]
fstv = fstdye[indmatch]
names(ATv) <- names(pCv) <- names(lamv) <- names(fstv)  <- toupper(kitdyes$Marker) #don't need to be same order as of popFreq

#plotEPG2(samples,kit0,refs, AT=ATv) #plot data
dat = prepareData(samples,refs,popFreq,threshT=ATv,minF=NULL, normalize = TRUE) #obtain data to use for analysis


test_that("check imported data :", {
  
  #Check imported evid sample
  expect_equal(samples[[1]][["VWA"]]$adata, paste(16:18) ) #Check alleles 
  expect_equal(samples[[1]][["VWA"]]$hdata, c(185,1573,93)) #check PHs
  

  #Check imported ref sample
  expect_equal(refs[[1]][["VWA"]]$adata, c("17","17")) #Check ref
  
  #Check freq
  freq = c(0.087,0.077,0.202,0.308,0.214,0.093,0.018,0.001)
  names(freq) = 14:21 
  expect_equal(popFreq$VWA, freq) #Check freq
  
})

test_that("check calcRMPfst:", {

  #Calculate random match prob for profile:
  rmp = calcRMPfst(dat,POIind=1,condInd=NULL,fst=fstv ) #non-scaled values

  #Check overall:
  LRupper = -sum(log10(rmp)) #maximum attainable LR (log10 scale)
  expect_equal(LRupper, 9.22783157308603) 
    
  #COMPARE rmp with MANUAL FORMULA (only few checked)
  locsCheck = names(rmp) #check all
  for(loc in locsCheck) {
    fst0 = fstv[loc]
    poi = refs[[1]][[loc]]$adata #get alleles of person of intererst genotype
    freq = popFreq[[loc]]
    freq0 = freq[ names(freq)%in%poi]
    
    ishom = length(freq0)==1
    if(ishom) {
      p1 = (2*fst0 + (1-fst0)*freq0[1])/(1+1*fst0)
      p2 = (3*fst0 + (1-fst0)*freq0[1])/(1+2*fst0)
      rmp0 = p1*p2
    } else {
      p1 = (1*fst0 + (1-fst0)*freq0[1])/(1+1*fst0)
      p2 = (1*fst0 + (1-fst0)*freq0[2])/(1+2*fst0)
      rmp0 = 2*p1*p2
    }
    expect_equal(as.numeric(rmp0),as.numeric(rmp[loc])) #compare with manual derived

  } #end for each
})


test_that("check maximum likelihood Hp:", {
  
  #CALC loglik under Hp:
  mle = contLikMLE(nC=1,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=1,xi=NULL,prC=pCv,nDone=2,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=NULL,seed=seed0)
  thhat=mle$fit$thetahat #obtain maximum likelihood estimates

  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #exctract params:  
  shape0 = 1/(thhat[2]^2) #get shape param
  scale0 = thhat[1]/shape0 #get scale param
  beta =  thhat[3] #degrad slope param
  xiB = thhat[4] #backwards stutter prop
  xiF = thhat[5] #forward stutter prop
    
  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
    freq = dat$popFreq[[loc]]
    Aset = names(freq) #get allele outcome 
    Aset0 = head(Aset,-1) #keep non-dropouts
    Astutt = c(as.numeric(Aset0)-1,as.numeric(Aset0)+1)
    Aset = c(Aset,setdiff(Astutt,Aset0)) #insert potential stutters last
    
    #Prepare base pair info for kits:
    bpv = rep(NA,length(Aset)) #obtain base pairs for kit
    subkit <- subset(kitinfo,toupper(kitinfo$Marker)==loc) #control on dyer
    for(aa in Aset) { #for each allele
      ind <- which(subkit$Allele==aa)
      if(length(ind)==0) ind <- which.min(abs(as.numeric(subkit$Allele) - as.numeric(aa))) #nearest neighbour (in allele name)-> maximum of the alleles
      bpv[aa==Aset] <- subkit$Size[ind] #insert base pair
    }
    
    #Assuming 1 known contribtion
    poi =  dat$refData[[loc]][[1]] #get alleles of POI
    nAG = matrix(0,nrow=length(Aset),ncol=1,dimnames = list(Aset,"poi")) #create allele counting matrix for unknown contributors
    nAG[,1] = table(factor(poi,level=Aset)) #insert contribution
    
    #Prepare PH data (same order as Aset):
    Y = rep(0,length(Aset))  #preparing PHs after extending allele vector
    Y[match(dat$samples[[1]][[loc]]$adata,Aset)] =dat$samples[[1]][[loc]]$hdata  #insert PHs
    
    mui <-  c(nAG[,1,drop=FALSE]*shape0) #expected contribution
    mui  <- mui*exp( (bpv-125)/100*log(beta) )  #in case of degradation
    
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
    if(length(psiDO)>0) vali = vali + sum(pgamma(ATv[loc],shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
    if(length(psiDI)>0) vali = vali + sum( dexp( Y[psiDI]-ATv[loc], rate=lamv[loc], log=TRUE) + log(pCv[loc]*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
    if(length(psiDI)==0) vali = vali + log(1-pCv[loc]) #in case of no dropin

    expect_equal(as.numeric(logLikv[loc]),as.numeric(vali))
  }      
    
  #Excpected joint values:
  expect_equal(as.numeric(thhat),c( 921.88130259,0.09885783,0.77117149,0.10356833,0.05099810)) #compare with manual derived
  expect_equal(mle$fit$loglik,-134.0967859) #check
    
  #Check cumulative probabilities (based on C++ code):
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-6,valid$ProbObs,c( 0.716982372630061,0.975447353703181,0.663136550668606,0.374704093478249,0.578846339734187,0.611381628955929,0.0881914359751329,0.297805547030208,0.0846878917395518,0.37187441289122,0.434020852720626,0.0454482079155214,0.0671375172510371,0.203773923345375,0.269376359633339,0.0900625186686192,0.423169613781737,0.95295264062683,0.133733257003317,0.781074688027564,0.214306247406648,0.0294977051814099,0.643651443145833,0.988745051039168,0.262148532739819))
})

test_that("check maximum likelihood Hd (unrelated):", {

  #CALC loglik under Hp:
  mle = contLikMLE(nC=1,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=0,knownRef = 1,xi=NULL,prC=pCv,nDone=2,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=NULL,seed=seed0)
  thhat=mle$fit$thetahat #obtain maximum likelihood estimates

  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #exctract params:  
  shape0 = 1/(thhat[2]^2) #get shape param
  scale0 = thhat[1]/shape0 #get scale param
  beta =  thhat[3] #degrad slope param
  xiB = thhat[4] #backwards stutter prop
  xiF = thhat[5] #forward stutter prop
  
  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
    #loc = locsCheck[1]
    
    freq = dat$popFreq[[loc]]
    Aset = names(freq) #get allele outcome 
    Aset0 = head(Aset,-1) #keep non-dropouts
    Astutt = c(as.numeric(Aset0)-1,as.numeric(Aset0)+1)
    Aset = c(Aset,setdiff(Astutt,Aset0)) #insert potential stutters last
    
    #Prepare base pair info for kits:
    bpv = rep(NA,length(Aset)) #obtain base pairs for kit
    subkit <- subset(kitinfo,toupper(kitinfo$Marker)==loc) #control on dyer
    for(aa in Aset) { #for each allele
      ind <- which(subkit$Allele==aa)
      if(length(ind)==0) ind <- which.min(abs(as.numeric(subkit$Allele) - as.numeric(aa))) #nearest neighbour (in allele name)-> maximum of the alleles
      bpv[aa==Aset] <- subkit$Size[ind] #insert base pair
    }
    
    #Prepare PH data (same order as Aset):
    Y = rep(0,length(Aset))  #preparing PHs after extending allele vector
    Y[match(dat$samples[[1]][[loc]]$adata,Aset)] =dat$samples[[1]][[loc]]$hdata  #insert PHs
    
    #Get outcome of unknowns
    poi =  dat$refData[[loc]][[1]] #get alleles of POI (known under Hd)
    Glist = calcGjoint(freq=freq,nU=1,fst=fstv[loc],refK=poi)
    Gset = Glist$G #get allele out come of unknowns

    pEvid = 0 #obtain probabiliity of evidence (sum over all genotypes)
    for(gind in 1:nrow(Gset)) { #traverse all genorypes
      
      #Assuming 1 known contribtion
      nAG = matrix(0,nrow=length(Aset),ncol=1,dimnames = list(Aset,"unknown")) #create allele counting matrix for unknown contributors
      nAG[,1] = table(factor(Gset[gind,],level=Aset)) #insert contribution
      
      mui <-  c(nAG[,1,drop=FALSE]*shape0) #expected contribution
      mui  <- mui*exp( (bpv-125)/100*log(beta) )  #in case of degradation
      
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
      if(length(psiDO)>0) vali = vali + sum(pgamma(ATv[loc],shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
      if(length(psiDI)>0) vali = vali + sum( dexp( Y[psiDI]-ATv[loc], rate=lamv[loc], log=TRUE) + log(pCv[loc]*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
      if(length(psiDI)==0) vali = vali + log(1-pCv[loc]) #in case of no dropin
      pEvid = pEvid + exp(vali)*Glist$Gprob[gind] #ind
    }    
    expect_equal(as.numeric(logLikv[loc]),as.numeric(log(pEvid)))
  }      
  
  #Excpected joint values:
  expect_equal(as.numeric(thhat),c( 921.88777210  ,    0.09885714   ,   0.77116993   ,   0.10356872  ,    0.05099952 )) #compare with manual derived
  expect_equal(mle$fit$loglik,-155.344653303007) #check
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-6,valid$ProbObs,c(0.71696899513367,0.975445161946618,0.663105946875243,0.374698748018734,0.578837036809496,0.611369826504089,0.0881900001996293,0.297801543974021,0.0846863235469782,0.371865205623908,0.433999182737134,0.0454415214924983,0.067135790031381,0.203763185593546,0.269355487664269,0.090056301694929,0.423165591958658,0.952950994193261,0.133731672879591,0.781069589536292,0.214296447013678,0.0294941771328225,0.643624576816175,0.988744077125576,0.262139516610987))

  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_equal(as.numeric(DC$table2[,2]),c(1,1,1,1,1,1,1))
})



test_that("check maximum likelihood Hd (sibling):", {
  ibd0 = c(1/4,1/2,1/4) #assuming the unknown is a sibling of poi
  mle = contLikMLE(nC=1,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=0,knownRef = 1, xi=NULL,prC=pCv,nDone=2,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=NULL,seed=seed0, knownRel = 1,ibd=ibd0 )
  thhat=mle$fit$thetahat #obtain maximum likelihood estimates
  logLikv = logLiki(mle) #obtain per marker resutls
  
  #Excpected joint values:
  expect_equal(sum(logLikv),mle$fit$loglik)
  expect_equal(as.numeric(thhat),c(921.81771878   ,   0.09883688   ,   0.77120574   ,   0.10356310  ,    0.05100558)) 
  expect_equal(mle$fit$loglik,-141.239022451048)
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
      
  #exctract params:  
  shape0 = 1/(thhat[2]^2) #get shape param
  scale0 = thhat[1]/shape0 #get scale param
  beta =  thhat[3] #degrad slope param
  xiB = thhat[4] #backwards stutter prop
  xiF = thhat[5] #forward stutter prop
  
  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
    #loc = locsCheck[2]
    
    freq = dat$popFreq[[loc]] 
    Aset = names(freq) #get allele outcome 
    Aset0 = head(Aset,-1) #keep non-dropouts
    Astutt = c(as.numeric(Aset0)-1,as.numeric(Aset0)+1)
    Aset = c(Aset,setdiff(Astutt,Aset0)) #insert potential stutters last
    
    #Prepare base pair info for kits:
    bpv = rep(NA,length(Aset)) #obtain base pairs for kit
    subkit <- subset(kitinfo,toupper(kitinfo$Marker)==loc) #control on dyer
    for(aa in Aset) { #for each allele
      ind <- which(subkit$Allele==aa)
      if(length(ind)==0) ind <- which.min(abs(as.numeric(subkit$Allele) - as.numeric(aa))) #nearest neighbour (in allele name)-> maximum of the alleles
      bpv[aa==Aset] <- subkit$Size[ind] #insert base pair
    }
    
    #Prepare PH data (same order as Aset):
    Y = rep(0,length(Aset))  #preparing PHs after extending allele vector
    Y[match(dat$samples[[1]][[loc]]$adata,Aset)] = dat$samples[[1]][[loc]]$hdata  #insert PHs
    
    #Get outcome of unknowns
    poi =  dat$refData[[loc]][[1]] #get alleles of POI (known under Hd)
    Glist = calcGjoint(freq=freq,nU=1,fst=fstv[loc],refK=poi,refR=poi,ibd=ibd0) #condition on genotype of poi
    Gset = Glist$G #get allele out come of unknowns
    nG = nrow(Gset) #get number of genotypes
    pEvid = 0 #obtain probabiliity of evidence (sum over all genotypes)
    for(gind in 1:nG) { #traverse all genorypes
      nAG = matrix(0,nrow=length(Aset),ncol=1,dimnames = list(Aset,"unknown")) #create allele counting matrix for unknown contributors
      nAG[,1] = table(factor(Gset[gind,],level=Aset)) #insert contribution
      
      mui <-  c(nAG[,1,drop=FALSE]*shape0) #expected contribution
      mui  <- mui*exp( (bpv-125)/100*log(beta) )  #in case of degradation
      
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
      if(length(psiDO)>0) vali = vali + sum(pgamma(ATv[loc],shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
      if(length(psiDI)>0) vali = vali + sum( dexp( Y[psiDI]-ATv[loc], rate=lamv[loc], log=TRUE) + log(pCv[loc]*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
      if(length(psiDI)==0) vali = vali + log(1-pCv[loc]) #in case of no dropin
      pEvid = pEvid + exp(vali + log(Glist$Gprob[gind]) ) #ind
    }    
    #sum(exp(rowSums(log(pMAT))))
    expect_equal(as.numeric(logLikv[loc]),as.numeric(log(pEvid)))
  } #end for each locus 
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-6,valid$ProbObs,c(0.717152792427557,0.97550656567467,0.663095378357748,0.37478619534339,0.578920711857432,0.611413161806438,0.0882103425865904,0.297730914703659,0.0846080057224975,0.37201124186449,0.434294475085395,0.0454639699742578,0.0671626403060103,0.203837623157615,0.269374333686294,0.0900733709011589,0.423294444239972,0.953011641147383,0.133778868047268,0.781169833342648,0.214369322047219,0.0295099408202252,0.643741447592566,0.988778229867342,0.262187688500659))
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_equal(as.numeric(DC$table2[,2]),c(1,1,1,1,1,1,1))
  
})

