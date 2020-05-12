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

test_that("check maximum likelihood Hp:", {
  NOC = 2
  #CALC loglik under Hp:
  cond = c(1,2)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,xi=NULL,prC=pCv,nDone=2,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=NULL,seed=seed0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates

  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #exctract params:  
  mx=thhat[1:NOC]
  shape0 = 1/(thhat[NOC+2]^2) #get shape param
  scale0 = thhat[NOC+1]/shape0 #get scale param
  beta =  thhat[NOC+3] #degrad slope param
  xiB = thhat[NOC+4] #backwards stutter prop
  xiF = thhat[NOC+5] #forward stutter prop
    
  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
#loc=locsCheck[1]
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
    
    vali = 0
    for(sample in names(dat$samples)) {
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
  expect_approx(1e-3,thhat,c( 0.79607926  ,    0.20392074  ,  794.83757183  ,    0.15096867 ,     0.74678609, 0.11871766   ,   0.05499105) ) #compare with manual derived
  expect_equal(mle$fit$loglik,-330.464144513856) #check
    
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-3,valid$ProbObs,c(0.756510873302842,0.216297025971664,0.621738469465481,0.260427169952819,0.935774421587891,0.83278899396675,0.168054546649229,0,0.254901765847725,0.112572234185492,0.0643522204872648,0.193018972424357,0.591854354795273,0.351192344778938,0.989931931763211,0.934540441340693,0.723807069666164,0.971264061517507,0.344191901226563,0.0652547406514865,0.748126125582802,0.466995116390264,0.38257646689127,0.220289225840389,0.583033359761642,0.0958103756111627,0.824613902053727,0.950317390046415,0.357402727789437,0.241811058503077,0.629865594643927,0.31682373048975,0.870418900887443,0.030795914444826,0.0824290074425967,0.618726238543472,0.712491826359744,0.0161017641964499,0.958644955785059,0.627904697612822,0.507775305497841,0.696402183505275,0.876225821224122,0.653971272221824,0.164864485234383,0.24201450463526,0.740212319974726,0.739361388336834,0.112338174373275,0.765132932740444,0.622227546290014,0.974434467731492,0.907151620219564,0.554083278596624,0.581479913899373,0.343066048410637,0.568977176250329))
  
})

test_that("check maximum likelihood Hd (unrelated):", {
  NOC = 2
  #CALC loglik under Hp:
  cond = c(1,0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = 2,xi=NULL,prC=pCv,nDone=2,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=NULL,seed=seed0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #exctract params:  
  mx=thhat[1:NOC]
  shape0 = 1/(thhat[NOC+2]^2) #get shape param
  scale0 = thhat[NOC+1]/shape0 #get scale param
  beta =  thhat[NOC+3] #degrad slope param
  xiB = thhat[NOC+4] #backwards stutter prop
  xiF = thhat[NOC+5] #forward stutter prop
  
  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
    #loc=locsCheck[2]
    freq = dat$popFreq[[loc]]
    Aset = names(freq) #get allele outcome 
    Aset0 = head(Aset,-1) #keep non-dropouts
    Astutt = c(as.numeric(Aset0)-1,as.numeric(Aset0)+1)
    Aset = c(Aset,setdiff(Astutt,Aset0)) #insert potential stutters last
    FWstutt = match(Aset,as.character(as.numeric(Aset)+1)) #alleleinds to stutter from
    BWstutt = match(Aset,as.character(as.numeric(Aset)-1)) #alleleinds to stutter from
    
    #Prepare base pair info for kits:
    bpv = rep(NA,length(Aset)) #obtain base pairs for kit
    subkit <- subset(kitinfo,toupper(kitinfo$Marker)==loc) #control on dyer
    for(aa in Aset) { #for each allele
      ind <- which(subkit$Allele==aa)
      if(length(ind)==0) ind <- which.min(abs(as.numeric(subkit$Allele) - as.numeric(aa))) #nearest neighbour (in allele name)-> maximum of the alleles
      bpv[aa==Aset] <- subkit$Size[ind] #insert base pair
    }
    
    #Assuming 1 known contribtion
    condRef =  dat$refData[[loc]][cond] #get alleles of conds
    nAG = matrix(0,nrow=length(Aset),ncol=NOC,dimnames = list(Aset)) #create allele counting matrix for unknown contributors
    for(ind in 1:length(condRef)) nAG[,ind] = table(factor(condRef[[ind]],level=Aset)) #insert contribution
    
    #Prepare PH data (same order as Aset):
    Ylist = list() #list per samples
    for(sample in names(dat$samples)) {
      Ylist[[sample]] = rep(0,length(Aset))  #preparing PHs after extending allele vector
      Ylist[[sample]][match(dat$samples[[sample]][[loc]]$adata,Aset)] = dat$samples[[sample]][[loc]]$hdata  #insert PHs
    }

    nAG0 = nAG #temporary store
    Glist = calcGjoint(freq=freq,nU=1,fst=fstv[loc],refK=unlist( dat$refData[[loc]]))
    Gset = Glist$G #get allele out come of unknowns

    pEvid = 0 #obtain probabiliity of evidence (sum over all genotypes)
    for(gind in 1:nrow(Gset)) { #traverse all genorypes
      nAG = nAG0
      nAG[,2] = table(factor(Gset[gind,],level=Aset)) #insert contribution
        
      mui <-  c(nAG[,drop=FALSE]%*%mx)*shape0 #expected contribution
      mui  <- mui*exp( (bpv-125)/100*log(beta) )  #in case of degradation
      
      #Delegate stutterprod
      stuttB <- mui[BWstutt]*xiB #backward-stutter parts
      stuttF <- mui[FWstutt]*xiF #forward-stutter parts
      indBW = !is.na(stuttB)
      indFW = !is.na(stuttF)
      indLooseStutt = indBW & indFW #index of alleles which are assumed to not loose stutter product
      mui[indLooseStutt] = mui[indLooseStutt]*(1- (xiB+xiF)) #loose stutter products
      mui[indBW] = mui[indBW] + stuttB[indBW]
      mui[indFW] = mui[indFW] + stuttF[indFW]
      
      vali = 0 
      for(sample in names(dat$samples)) { #traversing all samples
        #Divide set into dropin, contr and dropout
        psiDI <- which( Ylist[[sample]]>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
        psiYmu <- which( Ylist[[sample]]>0 & mui>0 )  #contributing to model and observed PH
        psiDO <- which( Ylist[[sample]]==0 & mui>0 )  #dropout elem
        
        if(length(psiYmu)>0) vali =  vali + sum(dgamma(Ylist[[sample]][psiYmu],shape=mui[psiYmu],scale=scale0,log=TRUE))   # CONTRIBUTION TO PH
        if(length(psiDO)>0) vali = vali + sum(pgamma(ATv[loc],shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
        if(length(psiDI)>0) vali = vali + sum( dexp( Ylist[[sample]][psiDI]-ATv[loc], rate=lamv[loc], log=TRUE) + log(pCv[loc]*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
        if(length(psiDI)==0) vali = vali + log(1-pCv[loc]) #in case of no dropin
      }
      pEvid = pEvid + exp( vali + log(Glist$Gprob[gind])) #sum up
    } #end for each genotype
    expect_equal(as.numeric(logLikv[loc]),as.numeric(log(pEvid)))
  }      
  
  #Excpected joint values:
  expect_approx(1e-3,thhat,c(0.763477081556841,0.236522918443159,796.610875241322,0.148066384014808,0.739731359928258,0.119369627520945,0.0504508269641699))
  expect_equal(mle$fit$loglik,-349.572712290044) #check
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-3,valid$ProbObs,c(0.875682169517986,0.212936448616972,0.693555629194084,0.3099441082868,0.908184891213207,0.769788793824497,0.144147813536516,0,0.279310877237671,0.11710753920421,0.043269076720267,0.168835353332635,0.571970955600711,0.329849054374258,0.989513697686217,0.919137984396742,0.645794580976664,0.956080732454599,0.334039818720927,0.0708036638837956,0.813135778419479,0.551709734017975,0.401285283704567,0.102916610322478,0.39739303732468,0.128300822236649,0.871277964925281,0.898331232684839,0.211578543664533,0.276517614702183,0.677597505383126,0.384678008960847,0.908564060320084,0.0418766233402921,0.106942748677193,0.540708613996736,0.64429064008271,0.0228161729625381,0.974285237651567,0.505598869736413,0.444479010892354,0.678497041748205,0.886611062936997,0.637311848247842,0.215328761121263,0.305925672230426,0.886063921122954,0.991455579878122,0.701253334190676,0.690610355374979,0.525134919495781,0.779984906153653,0.528225418079361,0.533086097371343,0.561223842696954,0.411886625791345,0.642138471043274))

  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_approx(1e-3,as.numeric(DC$table2[,5]),c(0.9250 ,0.8745, 0.5431, 1.0000 ,0.9988 ,0.8069, 0.9833))
  
})


test_that("check maximum likelihood Hd (sibling):", {
  ibd0 = c(1/4,1/2,1/4) #assuming the unknown is a sibling of poi
  RelInd = 2 #ref index of related individual
  NOC = 2
  #CALC loglik under Hp:
  cond = c(1,0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = 2,xi=NULL,prC=pCv,nDone=2,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=NULL,seed=seed0,knownRel = RelInd,ibd=ibd0)
  thhat=mle$fit$thetahat2 #obtain maximum likelihood estimates
  
  #COMPARE PER MARKER RESULTS WITH MANUAL DERIVED:
  logLikv = logLiki(mle) #obtain per marker resutls
  expect_equal(sum(logLikv),mle$fit$loglik)
  
  #exctract params:  
  mx=thhat[1:NOC]
  shape0 = 1/(thhat[NOC+2]^2) #get shape param
  scale0 = thhat[NOC+1]/shape0 #get scale param
  beta =  thhat[NOC+3] #degrad slope param
  xiB = thhat[NOC+4] #backwards stutter prop
  xiF = thhat[NOC+5] #forward stutter prop
  
  locsCheck = names(logLikv) #loci to check values for
  for(loc in locsCheck) {
    #loc=locsCheck[3]
    freq = dat$popFreq[[loc]]
    Aset = names(freq) #get allele outcome 
    Aset0 = head(Aset,-1) #keep non-dropouts
    Astutt = c(as.numeric(Aset0)-1,as.numeric(Aset0)+1)
    Aset = c(Aset,setdiff(Astutt,Aset0)) #insert potential stutters last
    FWstutt = match(Aset,as.character(as.numeric(Aset)+1)) #alleleinds to stutter from
    BWstutt = match(Aset,as.character(as.numeric(Aset)-1)) #alleleinds to stutter from
    
    #Prepare base pair info for kits:
    bpv = rep(NA,length(Aset)) #obtain base pairs for kit
    subkit <- subset(kitinfo,toupper(kitinfo$Marker)==loc) #control on dyer
    for(aa in Aset) { #for each allele
      ind <- which(subkit$Allele==aa)
      if(length(ind)==0) ind <- which.min(abs(as.numeric(subkit$Allele) - as.numeric(aa))) #nearest neighbour (in allele name)-> maximum of the alleles
      bpv[aa==Aset] <- subkit$Size[ind] #insert base pair
    }
    
    #Assuming 1 known contribtion
    condRef =  dat$refData[[loc]][cond] #get alleles of conds
    nAG = matrix(0,nrow=length(Aset),ncol=NOC,dimnames = list(Aset)) #create allele counting matrix for unknown contributors
    for(ind in 1:length(condRef)) nAG[,ind] = table(factor(condRef[[ind]],level=Aset)) #insert contribution
    
    #Prepare PH data (same order as Aset):
    Ylist = list() #list per samples
    for(sample in names(dat$samples)) {
      Ylist[[sample]] = rep(0,length(Aset))  #preparing PHs after extending allele vector
      Ylist[[sample]][match(dat$samples[[sample]][[loc]]$adata,Aset)] = dat$samples[[sample]][[loc]]$hdata  #insert PHs
    }
    
    nAG0 = nAG #temporary store
    Glist = calcGjoint(freq=freq,nU=1,fst=fstv[loc],refK=unlist( dat$refData[[loc]]), refR=dat$refData[[loc]][[RelInd]],ibd=ibd0)
    Gset = Glist$G #get allele out come of unknowns
    
    pEvid = 0 #obtain probabiliity of evidence (sum over all genotypes)
    for(gind in 1:nrow(Gset)) { #traverse all genorypes
      nAG = nAG0
      nAG[,2] = table(factor(Gset[gind,],level=Aset)) #insert contribution
      
      mui <-  c(nAG[,drop=FALSE]%*%mx)*shape0 #expected contribution
      mui  <- mui*exp( (bpv-125)/100*log(beta) )  #in case of degradation
      
      #Delegate stutterprod
      stuttB <- mui[BWstutt]*xiB #backward-stutter parts
      stuttF <- mui[FWstutt]*xiF #forward-stutter parts
      indBW = !is.na(stuttB)
      indFW = !is.na(stuttF)
      indLooseStutt = indBW & indFW #index of alleles which are assumed to not loose stutter product
      mui[indLooseStutt] = mui[indLooseStutt]*(1- (xiB+xiF)) #loose stutter products
      mui[indBW] = mui[indBW] + stuttB[indBW]
      mui[indFW] = mui[indFW] + stuttF[indFW]
      
      vali = 0 
      for(sample in names(dat$samples)) { #traversing all samples
        #Divide set into dropin, contr and dropout
        psiDI <- which( Ylist[[sample]]>0 & mui==0 ) #drop-in (not a stutter). ANY DROPIN?
        psiYmu <- which( Ylist[[sample]]>0 & mui>0 )  #contributing to model and observed PH
        psiDO <- which( Ylist[[sample]]==0 & mui>0 )  #dropout elem
        
        if(length(psiYmu)>0) vali =  vali + sum(dgamma(Ylist[[sample]][psiYmu],shape=mui[psiYmu],scale=scale0,log=TRUE))   # CONTRIBUTION TO PH
        if(length(psiDO)>0) vali = vali + sum(pgamma(ATv[loc],shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
        if(length(psiDI)>0) vali = vali + sum( dexp( Ylist[[sample]][psiDI]-ATv[loc], rate=lamv[loc], log=TRUE) + log(pCv[loc]*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
        if(length(psiDI)==0) vali = vali + log(1-pCv[loc]) #in case of no dropin
      }
      pEvid = pEvid + exp( vali + log(Glist$Gprob[gind])) #sum up
    } #end for each genotype
    expect_equal(as.numeric(logLikv[loc]),as.numeric(log(pEvid)))
  }      
  
  #Excpected joint values:
  expect_approx(1e-3,as.numeric(thhat),c(0.781811348820458,0.218188651179542,794.667708593774,0.151698944113907,0.745636804143434,0.118944752922951,0.0513713668220271))
  expect_equal(mle$fit$loglik,-337.188990585992) #check
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-3,valid$ProbObs,c(0.812359252930424,0.21251152200153,0.651825785583529,0.286932405533184,0.92532840933404,0.809636267821101,0.16063160365626,0,0.25294916221782,0.109362903222284,0.0623229106296523,0.194615568023339,0.591946730269557,0.351061065824279,0.990935831840615,0.939387872442509,0.684830443116803,0.962724850150029,0.313015471411921,0.0660023898747749,0.768313829574371,0.495003434916729,0.387484887415543,0.165126923747297,0.502404111958594,0.111287267972911,0.842130595825601,0.929251279446608,0.29325049129632,0.261813303071063,0.654822326607357,0.345793092924021,0.884335365050262,0.0362989256091051,0.093693728591745,0.588737125672951,0.685486740893559,0.0199012550489418,0.9638549125583,0.572060836401672,0.464802513814793,0.660664819080248,0.88600502202015,0.670697945598854,0.182585557274025,0.263562459732762,0.749688183136658,0.985013979088397,0.366305827549567,0.724355979524796,0.562052195881761,0.940452017558372,0.791269767626697,0.546030505192734,0.573670129594725,0.371119723792328,0.597336170101742))
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_equal(as.numeric(DC$table2[,5]),c(0.9757, 0.9547, 0.9670, 1.0000, 0.9995, 0.9783, 0.5534))
  
})
