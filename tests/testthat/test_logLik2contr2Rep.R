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
  expect_approx(1e-3,thhat,c(0.79607720896432,0.20392279103568,794.836439003443,0.150972128575762,0.746779030469506,0.11871725859296,0.0549843239305395)) #compare with manual derived
  expect_equal(mle$fit$loglik,-330.46414456887) #check
    
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-3,valid$ProbObs,c(0.756508215977985,0.21630427174326,0.621741075982208,0.260436351115485,0.935784465176423,0.832812874423049,0.168056637039667,0,0.254918848676835,0.11258488524009,0.0643661284147846,0.193044840718972,0.591851210467705,0.35119099790422,0.989931907856857,0.934542823969781,0.723800567801081,0.971261236517766,0.344189691038075,0.0652556405675384,0.748141283661917,0.467021405613272,0.382572894333458,0.220280196982473,0.583014628091662,0.0958112491892617,0.82460283776591,0.950311575204053,0.357395946971504,0.24183481252569,0.62989257748154,0.316820818994879,0.870409861152667,0.0307989698986432,0.0824340656466703,0.618745912542384,0.712507252064016,0.0161043826436387,0.958642921620643,0.627890986607015,0.507776533735793,0.696399759710498,0.8762276905707,0.653981200378729,0.164878450739975,0.242030072496172,0.740222944281938,0.739336297412009,0.112330236255872,0.765143347103222,0.622243086412992,0.974430796943293,0.907143628498147,0.55412136567563,0.581517050055367,0.343063665044419,0.568969708300312))
  
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
  expect_approx(1e-3,thhat,c(0.763486321081856,0.236513678918144,796.592283453393,0.148089673966996,0.739757394864906,0.119377184621431,0.0504547157244887))
  expect_equal(mle$fit$loglik,-349.57271366966) #check
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-3,valid$ProbObs,c(0.875711193923581,0.213052285800921,0.693514761949751,0.309944348670644,0.908156684742772,0.769755012676495,0.144102077116148,0,0.279345058808352,0.117142211417019,0.0432772452468894,0.168832130012311,0.57191057077521,0.329809045050305,0.989508026594311,0.919139251447746,0.64576541521204,0.956058481070563,0.333982477052052,0.0707810249467666,0.813057706740672,0.551636304158464,0.401219622290983,0.102990444672623,0.397500182348655,0.128369197823146,0.871274970693453,0.898338216299065,0.21167295095782,0.276474784711449,0.677516418319843,0.384747939708633,0.908554193329498,0.041904299242138,0.106986583130699,0.540723589387731,0.64429093181231,0.0228343491733412,0.974270628178829,0.50563891027821,0.444517493231431,0.67850836791512,0.886543896861955,0.637214243074982,0.215346970632388,0.305934273432303,0.88604337270328,0.991452328491644,0.701276881632775,0.690547922565004,0.525077182881645,0.780060775024309,0.528353860372558,0.533028098737082,0.561163187828175,0.411946600474555,0.642164079823821))

  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_approx(1e-3,as.numeric(DC$table2[,5]),c(0.9249,0.8744,0.5435,1,0.9988,0.8067,0.9833))
  
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
  expect_approx(1e-3,as.numeric(thhat),c(0.781819963994145,0.218180036005855,794.658693174281,0.151696493827377,0.745639242593977,0.118952975359393,0.0513751749723349))
  expect_equal(mle$fit$loglik,-337.188990826475) #check
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-3,valid$ProbObs,c(0.81239096315799,0.212513981286232,0.651856382842464,0.286952908626159,0.925338743530147,0.809651184915082,0.160621005361242,0,0.252966160992566,0.109370913064564,0.0623282816075049,0.194631685436729,0.591937022902966,0.351051743575428,0.990937906992317,0.939396212044697,0.684869907566975,0.962735823090934,0.313045736293998,0.0659997266203354,0.768326996733106,0.495014547368226,0.387465277843515,0.165162885263044,0.502468178053952,0.11129737554708,0.842153823601571,0.929276531611937,0.293303793857651,0.261777949156753,0.65478232318429,0.345817721065986,0.884354428053382,0.0363007302222458,0.0936990874403453,0.588746085379983,0.685496872809319,0.019902096335696,0.963862381312286,0.572116786322504,0.464821958654795,0.660687061192617,0.886018473869222,0.670717383461362,0.182591831244872,0.263571999378962,0.749722433210035,0.985016030499462,0.366289568303248,0.724343558550917,0.562039052784012,0.940495575066176,0.791384258934212,0.545988327620577,0.573628827629681,0.371145509169973,0.59736659628644))
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_equal(as.numeric(DC$table2[,5]),c(0.9757,0.9547,0.967,1,0.9995,0.9783,0.5532))
  
})
