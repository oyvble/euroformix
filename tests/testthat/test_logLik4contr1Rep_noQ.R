#Testing that the numerical calculation of loglik is correct
#rm(list=ls());library(euroformix);library(testthat)
seed0 = 1 #important to get reproducible results
kit0 = "GlobalFiler" #name of selected kit

expect_approx = function(tol,x,y) { #helpfunction with specified tolerance
  expect_equal(as.numeric(x),as.numeric(y),tolerance = tol)
}

examples = paste(path.package("euroformix"),"examples",sep=.Platform$file.sep)
popfn = paste(examples,paste0("GF_freq.csv"),sep=.Platform$file.sep)
evidfn = paste(examples,paste0("GF_evid.csv"),sep=.Platform$file.sep)
reffn = paste(examples,paste0("GF_refs.csv"),sep=.Platform$file.sep)

#Obtain data (taken from runexample):
popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
samples = sample_tableToList(tableReader(evidfn))#,threshT)
refs = sample_tableToList(tableReader(reffn))
kitinfo = getKit(kit0) #get kitinfo

#Set specific settings
AT = 75
pC = 0.05
lam = 0.01
fst = 0.01

#plotEPG2(samples,kit0,refs, AT=AT) #plot data

dat = prepareData(samples,refs,popFreq,threshT=AT,minF=NULL, normalize = TRUE) #obtain data to use for analysis
NOC = 4

#cond = 1:NOC
#prepC = prepareC(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,kit=kit0,fst=fst,incS=TRUE,incFS=TRUE)
#nC=NOC;samples=dat$samples;popFreq=dat$popFreq;refData=dat$refData;condOrder=cond;xi=NULL;prC=pC;nDone=2;threshT=AT;fst=fst;lambda=lam;kit=kit0;xiFW=NULL;knownRel=NULL;ibd=NULL;incS=TRUE;incFS=TRUE;knownRef=NULL

test_that("check maximum likelihood Hp:", {
  cond = 1:NOC #condition on all contributors
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,xi=NULL,prC=pC,nDone=2,threshT=AT,fst=fst,lambda=lam,kit=kit0,xiFW=NULL,seed=seed0)
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
# loc=locsCheck[5]
    freq = dat$popFreq[[loc]]
    Aset <- Aset0 <- names(freq) #get allele outcome 
    if( "99"%in%Aset ) Aset0 = head(Aset,-1) #keep only non-dropouts
    Astutt = c(as.numeric(Aset0)-1,as.numeric(Aset0)+1) #get stutter alleles
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
      if(length(psiDO)>0) vali = vali + sum(pgamma(AT,shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
      if(length(psiDI)>0) vali = vali + sum( dexp( Ylist[[sample]][psiDI]-AT, rate=lam, log=TRUE) + log(pC*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
      if(length(psiDI)==0) vali = vali + log(1-pC) #in case of no dropin
    }
    expect_equal(as.numeric(logLikv[loc]),as.numeric(vali))
  }      

  #Excpected joint values:
  expect_approx(1e-3,thhat,c(0.0849197106883292,0.261582140403597,0.309627608620431,0.343870540287642,17579.6310701948,0.15792270257044,0.992617873375767,0.0763488258988656,0.00669115722332238))
  expect_equal(mle$fit$loglik,-1205.66237498629) #check
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-3,valid$ProbObs,c(0.157092363482124,0.126063627217922,0.193287936798544,0.55215594757013,0.574136791076054,0.754306670291347,0.594391616394651,0.685650974431689,0.967838321755042,0.763958285503071,0.976177899037668,0.525118384969115,0.445851622744666,0.599656499609508,0.533996465999201,0.715060110221713,0.483945523272939,0.265866821933372,0.459377363675086,0.389537016797711,0.911833050322112,0.958695671207372,0.947410451992007,0,0.476839304632617,0.114365360639507,0.340428318084126,0.0291393016403834,0.262472845215376,0.0203466898917488,0.546970080664599,0.848184300447472,0.854599420076798,0.805236226617257,0.654889386006718,0.850350906206946,0.807980886950636,0.227444050180073,0.800684231255438,0.931437746262491,0.384695860087853,0.242979659429404,0.107455848345705,0.0256023629323732,0.041305383421994,0.0823988590588677,0.401256380233557,0.840220002286824,0.895370024920971,0.810886549256312,0.880308521235186,0.992066019499911,0.381152062062264,0.612651397121306,0.819225806351106,0.667221977200717,0.694945719087542,0.917142578182176,0.864754308106087,0.591163076840966,0.11299306935693,0.357562340671511,0.259208100442622,0.558816275801595,0.554498931592885,0.57312906578412,0.639449903672061,0.790428650745787,0.253916736591829,0.667141239060134,0.435772219110316,0.496628967582104,0.442042978741613,0.552995130749183,0.166324616343266,0.90061407474756,0.474615494596373,0.36422759522691,0.115639985756667,0.0144093684135088,0.0168515025310058,0.270251163204855,0.373378387024017,0.467118255376492,0.541073970791788,0.731743100075752,0.720268890001429,0.600303514328857,0.738537758767662,0.0858013129428207,0.444158470943382,0.777975869516141,0.426646474678559,0.940052767691144,0.56668931324846,0.260832220491447,0.972949265729017,0.656997515424709,0.769564130249784,0.545077318063061,0.971968702055493,0.5271970988612,0.922623371393376,0.232588012807899,0.42248834186471,0.965048607434935,0.402224081273666,0.665581111649183,0.146021125525541,0.0344698733132603,0.365155007123083,0.0306583744332375,0.104572591821129,0.192635887618537,0.813588915441737,0.68930129259635,0.738147251969069,0.388538760592555,0.519891162732848,0.524627142871568,0.751834210017747,0.605177298391863,0.749776058906259,0.71917171601619,0.156615537890267,0.422256814118398,0.636859255430608,0.0323348723921149,0.010693925629425,0.0513636765025808,0.0203311331055038,0.00607951480164596,0.0570880183007099,0.360582183208324,0.123723004479805,0.114296374113723,0.0122880388614752,0.31241978845628,0.4670015609873,0.65421235150988,0.355107343111459,0.118727832078229,0.142782547734817))
})

unknownIndHd = 4 #index of unknown

test_that("check maximum likelihood Hd (unrelated):", {
  
  #CALC loglik under Hd:
  cond = 1:NOC
  cond[unknownIndHd] = 0
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = unknownIndHd,xi=NULL,prC=pC,nDone=2,threshT=AT,fst=fst,lambda=lam,kit=kit0,xiFW=NULL,seed=seed0)
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
#  loc=locsCheck[5]
    freq = dat$popFreq[[loc]]
    Aset <- Aset0 <- names(freq) #get allele outcome 
    if( "99"%in%Aset ) Aset0 = head(Aset,-1) #keep only non-dropouts
    Astutt = c(as.numeric(Aset0)-1,as.numeric(Aset0)+1) #get stutter alleles
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
    Glist = calcGjoint(freq=freq,nU=1,fst=fst,refK=unlist( dat$refData[[loc]] ))
    Gset = Glist$G #get allele out come of unknowns
    
    pEvid = 0 #obtain probabiliity of evidence (sum over all genotypes)
    for(gind in 1:nrow(Gset)) { #traverse all genorypes
      nAG = nAG0
      nAG[,unknownIndHd] = table(factor(Gset[gind,],level=Aset)) #insert contribution
      
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
        if(length(psiDO)>0) vali = vali + sum(pgamma(AT,shape=mui[psiDO],scale=scale0,log=TRUE)) #CONTR TO DROPOUT
        if(length(psiDI)>0) vali = vali + sum( dexp( Ylist[[sample]][psiDI]-AT, rate=lam, log=TRUE) + log(pC*freq[psiDI]) ) #CONTR TO DROPIN: Scaled for each dropin. THIS IS OK
        if(length(psiDI)==0) vali = vali + log(1-pC) #in case of no dropin
      }
      pEvid = pEvid + exp( vali + log(Glist$Gprob[gind])) #sum up
    } #end for each genotype
    expect_equal(as.numeric(logLikv[loc]),as.numeric(log(pEvid)))
  }      

  #Excpected joint values:
  expect_approx(1e-3,thhat,c(0.0854118055923308,0.26140389133384,0.309537682873213,0.343646620200616,17513.0413473543,0.159256212760716,0.996708231117273,0.0768791996353518,0.00671117127638308))
  expect_equal(mle$fit$loglik,-1272.43904020949) #check
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-3,valid$ProbObs,c(0.155382068289545,0.124939715624677,0.217528905547355,0.250067266653976,0.255282468905418,0.614083409950597,0.59311793245331,0.690780615817632,0.99735835072076,0.760216412103784,0.976266314504108,0.530698185130165,0.441504759072782,0.602919145173533,0.503199008080787,0.709007610623681,0.476366771348657,0.263282084445915,0.453376339181195,0.388796745647391,0.969896416308529,0.980236683214538,0.944568084481545,0,0.455361811008114,0.112538404567331,0.350168563183662,0.0228393222527843,0.243522725422181,0.0193259923711185,0.543077549652613,0.84399126905821,0.93196295409851,0.832803683785416,0.651144224387977,0.846544467572072,0.803161545077326,0.226222887221802,0.797126498122717,0.929018007056172,0.383695947650136,0.244351620355276,0.100761281374397,0.0221001788300563,0.0427938926740758,0.0689305611979076,0.399957899239541,0.84077612409975,0.895183459199456,0.807203172396284,0.878605377622459,0.996944082014627,0.378267789535483,0.610259320966297,0.819206171508748,0.666956198855077,0.691626607202529,0.927747231494925,0.864068246934173,0.591965307300707,0.113238089562777,0.366505739769041,0.257799857815624,0.55794333098306,0.563048597460421,0.62244306055421,0.665585381243964,0.788683570657415,0.26121847128989,0.664763324026329,0.43094244058313,0.490906727339698,0.438626329197937,0.547619554037332,0.164316981381354,0.90411272882921,0.555083875512695,0.35864787712788,0.124208839377118,0.0208983143562177,0.000495626800867346,0.184593758941606,0.0321731114238114,0.466135739766282,0.545381571387036,0.721113497667015,0.714194997884136,0.599673553109791,0.822918761934741,0.0848279088818242,0.442158505233935,0.780766042879573,0.72631433310876,0.973539163868424,0.564821131743443,0.254872119865882,0.970778326459433,0.652023832437339,0.761089084014592,0.538384182935732,0.994506095808145,0.52352909860031,0.921014579509413,0.231628036058138,0.791655721167759,0.953280413960719,0.400356051005551,0.966580678903848,0.138259017077563,0.0265906904176967,0.23530358520313,0.0278087789417335,0.0546423706580726,0.187004155720449,0.808372388503044,0.683124954260183,0.736804984655174,0.384566775722738,0.512070994115135,0.520151193237597,0.744834216288526,0.600157757522063,0.741333334301132,0.711855659398977,0.155222759269925,0.413983417230158,0.632526879486163,0.0310149008297766,0.00124138303267101,0.00192754567980665,0.00387350396105958,0.00384045977420651,0.0537052265386087,0.0176123920191309,0.0165259650216947,0.0662938975573624,0.0087332522468615,0.310709956610394,0.466613105588043,0.623128561432522,0.361341768420317,0.116029310393548,0.161693407606026))
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_approx(1e-3,as.numeric(DC$table2[,ncol(DC$table2)-1]),c(0.9731,0.9948,0.9993,0.9865,0.9991,1,0.9997,0.9992,1,1,1,0.9377,0.9967,0.9937,0.9916,0.9595,0.998,1,0.9973,0.9985,0.9997))

})









