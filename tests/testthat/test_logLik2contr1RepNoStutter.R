#Testing that the numerical calculation of loglik is correct
#rm(list=ls());library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results
kit0 = "testkit" #name of selected kit

expect_approx = function(tol,x,y) { #helpfunction with specified tolerance
  expect_equal(as.numeric(x),as.numeric(y),tolerance = tol)
}

examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0("test_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("test_evids.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("test_refs.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))[1]#Use only 1st replicate
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


test_that("check maximum likelihood Hp: Ref1 + Ref2", {
  NOC = 2
  cond = c(1,2) #known contributors
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,xi=0,prC=pCv,nDone=2,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0,seed=seed0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #exctract params:  
  mx=thhat[1:NOC]
  shape0 = 1/(thhat[NOC+2]^2) #get shape param
  scale0 = thhat[NOC+1]/shape0 #get scale param
  beta =  thhat[NOC+3] #degrad slope param
  
  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
    freq = dat$popFreq[[loc]]
    Aset = names(freq) #get allele outcome 
    
    #Prepare base pair info for kits:
    bpv = rep(NA,length(Aset)) #obtain base pairs for kit
    subkit <- subset(kitinfo,toupper(kitinfo$Marker)==loc) #control on dyer
    for(aa in Aset) { #for each allele
      ind <- which(subkit$Allele==aa)
      if(length(ind)==0) ind <- which.min(abs(as.numeric(subkit$Allele) - as.numeric(aa))) #nearest neighbour (in allele name)-> maximum of the alleles
      bpv[aa==Aset] <- subkit$Size[ind] #insert base pair
    }
    
    #Assuming 2 known contribtions
    condRef =  dat$refData[[loc]][cond] #get alleles of conds
    nAG = matrix(0,nrow=length(Aset),ncol=NOC,dimnames = list(Aset)) #create allele counting matrix for unknown contributors
    for(ind in 1:length(condRef)) nAG[,ind] = table(factor(condRef[[ind]],level=Aset)) #insert contribution
    
    #Prepare PH data (same order as Aset):
    Ylist = list() #list per samples
    for(sample in names(dat$samples)) {
      Ylist[[sample]] = rep(0,length(Aset))  #preparing PHs after extending allele vector
      Ylist[[sample]][match(dat$samples[[sample]][[loc]]$adata,Aset)] = dat$samples[[sample]][[loc]]$hdata  #insert PHs
    }
    
    mui <-  c(nAG[,drop=FALSE]%*%mx)*shape0 #expected contribution
    mui  <- mui*exp( (bpv-125)/100*log(beta) )  #in case of degradation
    
    vali = 0
    for(sample in names(Ylist)) {
      #Divide set into dropin, contr and dropout
      psiDI <- which( Ylist[[sample]]>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
      psiYmu <- which( Ylist[[sample]]>0 & mui>0 )  #contributing to model and observed PH
      psiDO <- which( Ylist[[sample]]==0 & mui>0 )  #dropout elem
      
      if(length(psiYmu)>0) vali =  vali + sum(dgamma(Ylist[[sample]][psiYmu],shape=mui[psiYmu],scale=scale0,log=TRUE))   # CONTRIBUTION TO PH
      if(length(psiDO)>0) vali = vali + sum(pgamma(ATv[loc],shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
      if(length(psiDI)>0) vali = vali + sum( dexp( Ylist[[sample]][psiDI]-ATv[loc], rate=lamv[loc], log=TRUE) + log(pCv[loc]*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
      if(length(psiDI)==0) vali = vali + log(1-pCv[loc]) #in case of no dropin
    }
    expect_equal(as.numeric(logLikv[loc]),as.numeric(vali))
  }      
  
  #Excpected joint values:
  expect_equal(as.numeric(thhat),c( 0.723637935771717,0.276362064228283,688.465246180605,0.203496638155194,0.774111003592664)) #compare with manual derived
  expect_equal(mle$fit$loglik,-203.634478762721) #check
  
  #Check cumulative probabilities (based on C++ code):
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-6,valid$ProbObs,c(0.939649540948368,0.716111036213725,0.91421261622325,0.181269246922018,0.372020388994919,0.0949395497127479,0.45118836390597,0.966915299720932,0.424166394519259,0.154727468100078,0.0487705754992863,0.687059590644753,0.122710003692877,0.213047121429088,0.760399230864881,0.510318451414263,0.438384533825002,0.120346813651313,0.797026118838089,0.0691243829634655,0.315246578396999,0.548859774143454,0.854697828042385,0.279295354452448,0.257939040415481,0.838782355870223,0.961437757739786,0.776869839851568,0.456934409654664))
})


test_that("check maximum likelihood Hd: 2 unknown (unrelated):", {
  NOC = 2
  knownNonContrRefs = c(1,2) #known non-contributors
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,knownRef = knownNonContrRefs,xi=0,prC=pCv,nDone=2,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0,seed=seed0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #exctract params:  
  mx=thhat[1:NOC]
  shape0 = 1/(thhat[NOC+2]^2) #get shape param
  scale0 = thhat[NOC+1]/shape0 #get scale param
  beta =  thhat[NOC+3] #degrad slope param

  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
# loc=locsCheck[1]
    freq = dat$popFreq[[loc]]
    Aset = names(freq) #get allele outcome 
  
    #Prepare base pair info for kits:
    bpv = rep(NA,length(Aset)) #obtain base pairs for kit
    subkit <- subset(kitinfo,toupper(kitinfo$Marker)==loc) #control on dyer
    for(aa in Aset) { #for each allele
      ind <- which(subkit$Allele==aa)
      if(length(ind)==0) ind <- which.min(abs(as.numeric(subkit$Allele) - as.numeric(aa))) #nearest neighbour (in allele name)-> maximum of the alleles
      bpv[aa==Aset] <- subkit$Size[ind] #insert base pair
    }
    
    #Assuming NO known contributors
    nAG = matrix(0,nrow=length(Aset),ncol=NOC,dimnames = list(Aset)) #create allele counting matrix for unknown contributors

    #Prepare PH data (same order as Aset):
    Ylist = list() #list per samples
    for(sample in names(dat$samples)) {
      Ylist[[sample]] = rep(0,length(Aset))  #preparing PHs after extending allele vector
      Ylist[[sample]][match(dat$samples[[sample]][[loc]]$adata,Aset)] = dat$samples[[sample]][[loc]]$hdata  #insert PHs
    }

    nAG0 = nAG #temporary store contribution outcome
    Glist = calcGjoint(freq=freq,nU=2,fst=fstv[loc],sortComb=FALSE,refK=unlist( dat$refData[[loc]][knownNonContrRefs]) )
    Gset = Glist$G #get allele out come of unknowns
    
    #prod(dim(Glist$Gprob)) # number of combinations to traverse
    pEvid = 0 #obtain probabiliity of evidence (sum over all genotypes)
    for(c1 in 1:nrow(Gset)) { #contributor 1 genotypes
      for(c2 in 1:nrow(Gset)) { #contributor 2 genotypes
        nAG = nAG0
        nAG[,1] = table(factor(Gset[c1,],level=Aset)) #insert contribution
        nAG[,2] = table(factor(Gset[c2,],level=Aset)) #insert contribution
          
        mui <-  c(nAG[,drop=FALSE]%*%mx)*shape0 #expected contribution
        mui  <- mui*exp( (bpv-125)/100*log(beta) )  #in case of degradation
              
        vali = 0 
        for(sample in names(Ylist)) { #traversing all samples
          #Divide set into dropin, contr and dropout
          psiDI <- which( Ylist[[sample]]>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
          psiYmu <- which( Ylist[[sample]]>0 & mui>0 )  #contributing to model and observed PH
          psiDO <- which( Ylist[[sample]]==0 & mui>0 )  #dropout elem
          
          if(length(psiYmu)>0) vali =  vali + sum(dgamma(Ylist[[sample]][psiYmu],shape=mui[psiYmu],scale=scale0,log=TRUE))   # CONTRIBUTION TO PH
          if(length(psiDO)>0) vali = vali + sum(pgamma(ATv[loc],shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
          if(length(psiDI)>0) vali = vali + sum( dexp( Ylist[[sample]][psiDI]-ATv[loc], rate=lamv[loc], log=TRUE) + log(pCv[loc]*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
          if(length(psiDI)==0) vali = vali + log(1-pCv[loc]) #in case of no dropin
        }
        pEvid = pEvid + exp( vali + log(Glist$Gprob[c1,c2])) #sum up
      } #end for each genotype of C2
    } #end for each genotype of C1
    expect_equal(as.numeric(logLikv[loc]),as.numeric(log(pEvid)))
  }      
  
  #Excpected joint values:
  expect_approx(1e-3,thhat,c(0.718560679161608,0.281439320838392,699.018591591508,0.231186789294128,0.770818686051447))
  expect_equal(mle$fit$loglik,-234.095699644915) #check
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-3,valid$ProbObs,c(0.748907977796133,0.958605205699638,0.729244910702136,0.0182122947504831,0.674437437167819,0.260839819881155,0.259901440316275,0.951215393909042,0.46753958399596,0.236782114784974,0.0281382281720188,0.694747921504197,0.194628786595511,0.30312890660941,0.633833974741386,0.200111878822325,0.490582172976337,0.268457678901422,0.33307175170677,0.223354718274522,0.144210073837861,0.329811302193477,0.641182960658797,0.299482734815411,0.898110537317147,0.459356205884671,0.969224431217337,0.470075398724615,0.736372060891528))

  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_approx(1e-3,as.numeric(DC$table2[,5]),c(0.6269,0.5201,0.3149,0.752,0.7833,0.5257,0.5703))
})

