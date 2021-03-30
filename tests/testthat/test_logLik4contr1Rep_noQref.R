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

#Force the reference to have a an allele not observed AND other observed spans the freq data
refs[["K48"]][["D16S539"]]$adata = c("7","13") #insert another allele to reference

#Set specific settings
AT = 75
pC = 0.05
lam = 0.01
fst = 0.01

#plotEPG2(samples,kit0,refs, AT=AT) #plot data
dat = prepareData(samples,refs,popFreq,threshT=AT,minF=NULL, normalize = TRUE) #obtain data to use for analysis

#dat$refData[["D16S539"]][["K48"]]
#The new
#dat$popFreq$D16S539  
#all(sapply(dat$popFreq,function(x) "99"%in%names(x)))
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
    indBW = !is.na(stuttB)
    indFW = !is.na(stuttF)
    
    indLooseStutt = indBW & indFW #index of alleles which are assumed to not loose stutter product
    mui[indLooseStutt] = mui[indLooseStutt]*(1- (xiB+xiF)) #loose stutter products
    
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
  expect_approx(1e-6,thhat,c(0.0709877463436237,0.269099895982247,0.312072443715005,0.347839913959124,17581.8975096388,0.166616668703413,0.992389760520995,0.0803986277103819,0.00715714253517648))
  expect_equal(mle$fit$loglik,-1220.27569260485) #check
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-6,valid$ProbObs,c(0.150521828930323,0.129017144846988,0.228467861760927,0.54845453232446,0.551025448572646,0.743603922306363,0.56670567065844,0.647078420516996,0.96552203563591,0.72944119321914,0.966315706325333,0.574966328826211,0.419825369471192,0.604351740670531,0.530779080317959,0.689028177437823,0.473366970133889,0.25956138281631,0.589555210069244,0.377700797128753,0.9081415102761,0.952197003861505,0.930257294134962,0,0.964001426693287,0.10920533625107,0.33530997762538,0.0340612120318416,0.293474380495426,0.0194269101084747,0.533843565425133,0.851229776832529,0.834769887603082,0.789004577230323,0.6299599818825,0.835661611225969,0.818805988014892,0.217879776149216,0.777318049250994,0.920933506949827,0.368747314805291,0.243856848168899,0.107623664401371,0.0270486923189509,0.0621494510618991,0.115999691326595,0.384658179743412,0.823051147301873,0.873469444032969,0.858437615214214,0.882426171909721,0.987195530947229,0.364253291721414,0.587609259766233,0.800939298219282,0.690162826630711,0.679819339876287,0.917923850515752,0.838681344126468,0.585603769042774,0.108831448018924,0.3648478423452,0.258051541471364,0.634713164566059,0.554567167034126,0.556454544465003,0.62695892631951,0.765704240827498,0.265292231302275,0.645735683574142,0.416713882987989,0.505963808143047,0.42398868669,0.544942215770793,0.174319311469953,0.883018168456399,0.474148264898436,0.499660375506832,0.11571638742421,0.0211783278118127,0.0225630818400991,0.259217784150249,0.379944042068196,0.443031127951348,0.519420319805187,0.79443176654201,0.699268721887667,0.569357145719495,0.749275397363761,0.082021699603366,0.4107619320649,0.751107799593617,0.459682316792766,0.923220066300423,0.679405947435602,0.254343088765757,0.969164537942175,0.630664883519954,0.841085502465694,0.507879816115286,0.963134227383912,0.496217691619821,0.903896664425648,0.221246254565253,0.414759800923582,0.953596067689215,0.538813384329159,0.652104598968089,0.147282660371762,0.0521566257374797,0.34071781457721,0.0316618285278525,0.105681210171069,0.319753130048458,0.801585159486496,0.653972173142554,0.709814806829492,0.384471344006776,0.645506320173346,0.503781894883697,0.741373188469528,0.580531554529097,0.729841342076075,0.689034487163202,0.156289657382338,0.55845082201917,0.613282505914422,0.0309109861862066,0.0126674286563066,0.0763320622334614,0.0237267155894622,0.0100038910658077,0.0687080988554974,0.49835456825794,0.204258532192271,0.112262050095937,0.0135503391049485,0.299573511037875,0.472656904735171,0.714274601558259,0.371393336584955,0.118231176986192,0.151743743639126))
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
  
#  paste(as.numeric(DC$table2[,ncol(DC$table2)-1]),collapse=",")
  
  expect_approx(1e-6,thhat,c(0.0718795847924063,0.268765260016843,0.312470218687185,0.346884936503565,17582.7610742919,0.169407047060386,0.992376285603184,0.0810860103697447,0.00727219542589593))
  expect_equal(mle$fit$loglik,-1286.86282716154) #check
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-6,valid$ProbObs,c(0.150667777756406,0.135684406532004,0.273328175271617,0.275128911225095,0.269357688714114,0.629771888638613,0.562950593361656,0.646903637048789,0.995641418915122,0.720276785559275,0.963396384688748,0.570693692703241,0.417503787749528,0.611023722744812,0.490369278575144,0.681269739465804,0.470026636952463,0.260818561509788,0.579860432694256,0.379421675788413,0.964937546134779,0.975737661381746,0.924591075212005,0,0.0470546328256264,0.107372403829423,0.761669204861663,0.0270116499363608,0.270988547126152,0.0181115700495052,0.529918958846438,0.84701017228698,0.933746005669016,0.828898038823402,0.62517515920215,0.831437392671024,0.814918976570205,0.21499017651114,0.771979493114926,0.917771562092293,0.366454572283092,0.2441848917711,0.104711848735701,0.0258220132418516,0.0658587242169975,0.105692969515725,0.382238243376732,0.818745536957146,0.86999597888226,0.850003636647242,0.87536555335164,0.995427254007237,0.356874590216717,0.584639318825447,0.800312654402129,0.687389501805001,0.675392461995144,0.926318898061658,0.835418921786238,0.584741833827869,0.108757644774108,0.36776238275041,0.254587398626465,0.626559487667932,0.554834734202588,0.618089335590123,0.658867856094704,0.760305336240191,0.269374563748496,0.639347260807439,0.415248723943062,0.507720467668478,0.421192049724362,0.543987750952979,0.178283514407582,0.893920264587571,0.562840518696976,0.494581618265675,0.128921799940487,0.0324550682470925,0.000898330075145238,0.19037502532434,0.0446240002111826,0.440418903165214,0.519391272096148,0.773097282183519,0.68110041433906,0.559445710508233,0.833384928643557,0.0801409622539663,0.407868292841542,0.74621704623658,0.721294747107602,0.962653670821861,0.671998843119438,0.248528832648162,0.965156678027908,0.625483538855598,0.825908097222906,0.50269903992736,0.992787395075306,0.490368105918689,0.896255353198215,0.220212101394028,0.767010022007997,0.932118868523459,0.524460515130642,0.951486838210987,0.140010630683135,0.0391848763624229,0.210711083526561,0.03062694772947,0.0559665064824078,0.310141524089663,0.798135705259929,0.650474705002356,0.722612059847482,0.377023919901509,0.638199293713971,0.500124571426681,0.738343977577655,0.577551929228537,0.730470347467843,0.688953613483702,0.154030712049479,0.552371478811467,0.60864480896565,0.0307448533415152,0.0022712096352937,0.00462547894871486,0.00705669468866013,0.00642694994362807,0.0690266078222601,0.0284984522746116,0.0288132598762784,0.0743722665546151,0.010597059465195,0.297206148097956,0.468268426624323,0.669371095917231,0.380523537259943,0.115833148577734,0.185739324081408))
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_approx(1e-6,as.numeric(DC$table2[,ncol(DC$table2)-1]),c(0.9649,0.9934,0.9984,0.9793,0.8579,1,0.9996,0.9974,1,1,1,0.9137,0.9948,0.9928,0.9831,0.9462,0.9954,0.9998,0.994,0.9963,0.9992))
})









