#Testing that the numerical calculation of loglik is correct
#rm(list=ls());library(euroformix);library(testthat)
#library(euroformix);library(testthat)
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
  expect_equal(as.numeric(thhat),c(921.89699446,0.09882445,0.77113515,0.10358234,0.05100923)) #compare with manual derived
  expect_equal(mle$fit$loglik,-134.0967910) #check
    
  #Check cumulative probabilities (based on C++ code):
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-6,valid$ProbObs,c(0.716848735470755,0.975510422635757,0.663014597615598,0.374688874838355,0.579043517859287,0.611595010879323,0.0881941099266706,0.297987487625436,0.0847394063702697,0.371779616625961,0.43397321673334,0.0453934890944117,0.0671227331388472,0.203778649860947,0.269150314294627,0.0900432223429808,0.423216519265083,0.953051294066193,0.133749975065192,0.781300925973849,0.214125408308847,0.0294640903832653,0.643421265143422,0.988771424287719,0.262167100193662))
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
  expect_equal(as.numeric(thhat),c( 921.99621854,0.09880220,0.77111265,0.10355932,0.05100468 )) #compare with manual derived
  expect_equal(mle$fit$loglik,-155.34466387623) #check
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-6,valid$ProbObs,c(0.717060759455611,0.975444267208888,0.66305049890436,0.37482403573331,0.57880766299698,0.611316950226889,0.0882332700214559,0.297723105080007,0.0846093965943636,0.371903546642005,0.433520687461352,0.0452410884635613,0.0671514094443232,0.203416619167846,0.269223410605214,0.0898246808237023,0.423358097186722,0.952992425527182,0.133806028975264,0.78108990619852,0.214199187182661,0.0293524798504319,0.643540928996117,0.988747149952457,0.26220677774933))

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
  expect_equal(as.numeric(thhat),c(921.891318559959,0.0988551533869736,0.771167338641696,0.103570748609616,0.0510011348459797)) 
  expect_equal(mle$fit$loglik,-141.2390191251)
  
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
  expect_approx(1e-6,valid$ProbObs,c(0.716940777332816,0.975448777271414,0.663078541711507,0.374687899068804,0.578845171885752,0.611380032937213,0.0881878324813739,0.29781455217932,0.0846907068295644,0.371846607291081,0.433989465231522,0.0454371015194454,0.0671322662063131,0.203763765949529,0.269323370303125,0.0900555076339959,0.423159813731808,0.952956864151953,0.133729256050488,0.781083621979113,0.214273543997062,0.02949219087587,0.643586446611134,0.98874549852057,0.262134086866158))
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_equal(as.numeric(DC$table2[,2]),c(1,1,1,1,1,1,1))
  
})

