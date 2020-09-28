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
samples = sample_tableToList(tableReader(evidfn))[1]
refs = sample_tableToList(tableReader(reffn))
locs = names(popFreq)#[10:20]

#Set common settings for all markers:
ATv = 500 #put high to remove alleles in markers
pCv = 0.05
lamv = 0.01
fstv = 0.02


test_that("Test that wrong data format causes error:", {

#Filter data based on detection threshold
samples2 =samples[[1]] #copy data
for(loc in names(samples2)) {
  keep =  samples2[[loc]]$hdata>=ATv
  samples2[[loc]]$hdata = samples2[[loc]]$hdata[keep]
  samples2[[loc]]$adata = samples2[[loc]]$adata[keep]
}

errorOccured = FALSE 
tryCatch({
  mle = contLikMLE(nC=2,samples=list(Evid=samples2),popFreq=popFreq,xi=0,prC=pCv,nDone=2,threshT=ATv,fst=fstv,lambda=lamv,xiFW=0,seed=seed0)
},error = function(e) errorOccured <<- TRUE) #print(e))
expect_equal(errorOccured,TRUE) #CHECK that function call caused error
})

#plotEPG2(samples,kit0,refs, AT=ATv) #plot data
dat = prepareData(samples,refs,popFreq[locs],threshT=ATv,minF=NULL, normalize = TRUE) #obtain data to use for analysis
#plotMPS2(dat$samples,AT = ATv,refData = dat$refData,locYmax = FALSE)
#sapply(dat$refData,function(x) lapply(x, function(y) length(y)))

test_that("check maximum likelihood Hp (P1 + 1U):", {
  
  NOC = 2
  #CALC loglik under Hp:
  cond = c(1,0)
  mle = contLikMLE(nC=NOC,samples=dat$samples,popFreq=dat$popFreq,refData=dat$refData,condOrder=cond,knownRef = NULL,xi=0,prC=pCv,nDone=2,threshT=ATv,fst=fstv,lambda=lamv,kit=kit0,xiFW=0,seed=seed0)
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
# loc=locsCheck[1]
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
    Glist = calcGjoint(freq=freq,nU=1,fst=fstv,refK=unlist(condRef))
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
  expect_approx(1e-6,thhat,c(0.800719189324181,0.199280810675819,730.581491903075,0.631339725096235))
  expect_equal(mle$fit$loglik,-1409.45126218041) 
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE) #check 20 first only
  expect_approx(1e-6,valid$ProbObs[1:20],c(0.590179643627026,0.612445505503976,0.249519866401614,0.909564130104912,0.683045750865029,0.839829285178395,0.147909480600492,0.158977961867369,0.0286730633727604,0.665960598018581,0.721428156047013,0.134485356394851,0.435450094431388,0.454699113403926,0.0505784615059424,0.525840861483247,0.722767113367438,0.745758258760812,0,0.105644701649654))

  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_approx(1e-6,DC$table2[1:20,5],c(0.5145,0.4994,0.5758,0.6303,0.4802,0.5059,0.6409,0.4862,0.5197,0.4855,1,0.5061,0.5784,0.4644,0.487,0.5109,0.4976,0.5175,1,1))

})

test_that("check maximum likelihood Hd: 2 unknown (unrelated):", {
  
  NOC = 2
  knownNonContrRefs = 1 #known non-contributors
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
    
    #Assuming NO known contributors
    nAG = matrix(0,nrow=length(Aset),ncol=NOC,dimnames = list(Aset)) #create allele counting matrix for unknown contributors
    
    #Prepare PH data (same order as Aset):
    Ylist = list() #list per samples
    for(sample in names(dat$samples)) {
      Ylist[[sample]] = rep(0,length(Aset))  #preparing PHs after extending allele vector
      Ylist[[sample]][match(dat$samples[[sample]][[loc]]$adata,Aset)] = dat$samples[[sample]][[loc]]$hdata  #insert PHs
    }
    
    nAG0 = nAG #temporary store contribution outcome
    Glist = calcGjoint(freq=freq,nU=2,fst=fstv,sortComb=FALSE,refK=unlist( dat$refData[[loc]][knownNonContrRefs]) )
    Gset = Glist$G #get allele out come of unknowns
    
    #prod(dim(Glist$Gprob)) # number of combinations to traverse
    pEvid = 0 #obtain probabiliity of evidence (sum over all genotypes)
    for(c1 in 1:nrow(Gset)) { #contributor 1 genotypes
      for(c2 in 1:nrow(Gset)) { #contributor 2 genotypes
        nAG = nAG0
        nAG[,1] = table(factor(Gset[c1,],level=Aset)) #insert contribution
        nAG[,2] = table(factor(Gset[c2,],level=Aset)) #insert contribution
        
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
        pEvid = pEvid + exp( vali + log(Glist$Gprob[c1,c2])) #sum up
      } #end for each genotype of C2
    } #end for each genotype of C1
    expect_equal(as.numeric(logLikv[loc]),as.numeric(log(pEvid)))
  }      
  
  #Excpected joint values:
  expect_approx(1e-6,thhat,c(0.50000401203235,0.49999598796765,743.962888261923,0.668780008839223))
  expect_equal(mle$fit$loglik,-1470.39533062578) 
  
  valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,alpha=0.01,verbose = FALSE)
  expect_approx(1e-6,valid$ProbObs[1:20],c(0.572632045881159,0.574282162800037,0.364361130476897,0.753157582617469,0.833177268714064,0.856400073731605,0.141516149955756,0.128769132130379,0.0603426142223828,0.703634471949587,0.628205626639378,0.14985072753912,0.41144814210303,0.420137599252618,0.10911109858334,0.655060519061695,0.762131114749059,0.803781634408062,0,0.0725784849555745))
  
  #Calculate deconvolution (DC):
  DC = deconvolve(mle,verbose=FALSE)
  expect_approx(1e-6,DC$table2[1:20,5],c(0.5988,0.5135,0.5792,0.8076,0.587,0.5139,0.7924,0.5522,0.5988,0.5426,1,0.5012,0.7621,0.6611,0.4877,0.5657,0.5838,0.5835,1,1))
  
 })

