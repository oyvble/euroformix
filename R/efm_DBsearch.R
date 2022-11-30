#' @title efm_DBsearch
#' @author Oyvind Bleka
#' @description A EFM GUI helpfunction to run the Database search LR comparisons
#' @param mmTK The EFM GUI environment variable (CONTAINS DATA AND MODEL FOR CALCULATIONS)
#' @param ITYPE Type of LR interpretation (QUAN=EFM or QUAL=LRmix)
#' @export

efm_DBsearch = function(mmTK, ITYPE="QUAN") {
  set <- get("setDB",envir=mmTK) #get setup for DB
  popFreq <- set$popFreq #get original saved popfreq 
  popFreqQ <- set$popFreqQ #get popFreq saved by setup
  locs_hd <- names(popFreqQ) #get locus analysed under hd
  mixSel <- names(set$samples) #get name of selected mixtures
  refData <- set$refDataQ #take out selected references         
  fst = set$fst[1] #use only first fst element for qual model
  pC <- get("optDB",envir=mmTK)$QUALpC #get drop-in parameter for qual model from option
  
  if(ITYPE=="QUAN") {
    mleobj <- set$mlefit_hd #get object under hd
    opt <- get("optMLE",envir=mmTK)
    logLi_hd <- logLiki(mleobj) #already calculated, maxThreads=opt$maxThreads) #get log-probabilities for each locus (under Hd)
  }
  
  #LRmix settings:
  nsample <- get("optLRMIX",envir=mmTK)$nsample
  totA <-  sapply(  set$samples, function(x) sum(sapply(x,function(y) length(y$adata)) ) ) #number of alleles for each samples
  refHd <- NULL
  if( any(set$condOrder_hd>0) ) refHd <- lapply(set$refData ,function(x) x[set$condOrder_hd]) #must have the original refData!
  pDhat <- rep(0.1,length(totA)) #defualt value
  if(ITYPE=="QUAL") {
    #print("Obtaining drop-out estimate with maximum likelihood estimate")
    print("Estimating allele dropout probability based on MonteCarlo method...")
    for(ss in 1:length(totA)) { #for each sample (do dropout estimation) under Hd
      Msamp <- max(2000,25*totA[ss]) #number of samples for each vectorization
      DOdist <- simDOdistr(totA=totA[ss],nC=set$nC_hd,popFreq=popFreq,refData=refHd,minS=nsample,prC=pC,M=Msamp)
      pDhat[ss] <- quantile(DOdist ,0.5) #get median
    }
    print(paste0("Median(s) is given as pDhat=",pDhat))
    pDhat[is.na(pDhat)] <- 0.1 #impute 0.1
  }
  pDhat <- round(mean(pDhat),2) #used pDhat in database search
  #     print(paste0("Mean of the medians over each sample was given as pDhat=",pDhat))
  print(paste0("Dropout probability used for DB-search (Qual):",pDhat))
  
  #Prepare variables
  nU_hp <- set$nC_hp - sum(set$condOrder_hp>0) #number of unknowns under Hp                    
  nU_hd <- set$nC_hd - sum(set$condOrder_hd>0) #number of unknowns under Hd                    
  DBtab <- numeric()  #used to store table when searched
  locevid <- unique(unlist(unlist(lapply( set$samples, function(x) names(x) )))) #get unique locus names
  prim = .getPrim()  #get prime numbers. max length of alleles for reference-database storage is 244
  dbData <- get("dbData",envir=mmTK) #load all DB data
  for(dsel in set$dbSel) { #for each selected databases
    subD <- dbData[[dsel]] #get selected database
    dblocs <- toupper(colnames(subD)) #get database locs
    indD <- rownames(subD) #get individual names in database
    macD <- rep(0,length(indD)) #matching allele counter for each reference
    nlocs <- rep(0,length(indD)) #Number of loci which are used for calculating score - Note: Require that reference in DB has a genotype
    LR1 <- rep(1,nrow(subD)) #LRmix vec
    dblocs <- locevid[locevid%in%dblocs] #consider only loci within sample
    
    
    ########################################################################
    #step 1) convert allele-names of elements in database to one in popFreq
    for(loc in dblocs) { #for each locus in db     
      if(is.null(popFreq[[loc]])) next #skip to next locus
      Ainfo <- names(unlist(popFreq[[loc]])) #extract allele-info of frequncies
      #translate database to original genotypes
      Pinfo <- prim[1:length(Ainfo)]
      G = t(as.matrix(expand.grid(rep(list(Ainfo,Ainfo )))))
      GP = t(as.matrix(expand.grid(rep(list(Pinfo,Pinfo )))))
      keep = GP[2,]>=GP[1,] #unique genotypes
      G <- G[,keep]  #store genotypes
      GP <- GP[1,keep]*GP[2,keep] #get prime product
      G[!G%in%names(popFreqQ[[loc]])] <- "99" #Rename missing alleles in popFreqQ to "99":
      G0 <- paste0(G[1,],paste0("/",G[2,])) #make db-ref in one vector only
      
      newRow <- rep(NA,length(indD))  
      for(j in 1:length(GP)) { #for each genotype in population: Check in database
        #Always: Find corresponding genotype name by looking up stored primenumbers (same order as in popFreq!)
        rowind <- which(subD[,which(loc==colnames(subD))]==GP[j]) #samples with this genotype
        if(length(rowind)==0) next
        newRow[rowind] <- G0[j]
        #Get MAC of references:
        tmpmac <- 0
        #Count matching alleles over all mixtures:
        for(msel in mixSel) { #for each selected mixture
          evid0 <- set$samples[[msel]][[loc]]$adata
          tmpmac <- tmpmac + sum(G[,j]%in%evid0)
        } 
        macD[rowind] = macD[rowind] + tmpmac #add match counter  
        nlocs[rowind] <- nlocs[rowind] + 1 #counted only once!
      } #end for each genotype
      subD[,which(loc==colnames(subD))] <- newRow #force insertion of genotype-names
    }
    
    #########################################
    #STEP 2) CALCULATE Qualitative LR vals with dropout prob pDhat
    #step 2) qualitative LRs (always done alongside)
    LR1 <- rep(1,nrow(subD)) #LRmix vec
    for(loc in dblocs ) { #for each locus in db 
      if(is.null(popFreq[[loc]])) next #skip to next locus
      evidlist <- lapply( set$samples, function(x) x[[loc]]$adata ) #take out sample data:
      condRefHp <- unlist(refData[[loc]][set$condOrder_hp] ) #take out known refs under Hp 
      dbR <- subD[,which(loc==colnames(subD))] #take out DB-refs
      isNA <- is.na(dbR) #indicate references with missing loci
      #dbR[is.na(dbR)] <- 0 #insert zero
      #if(all(isNA)) next #skipt locus if none to calculate
      dbR2 <- matrix(NA,nrow=nrow(subD),ncol=2) #create a matrix with NA
      dbR2[!isNA,] <- t(matrix(unlist(strsplit(dbR[!isNA] , "/")) ,nrow=2)) #store into new matrix
      move = as.numeric(dbR2[,2]) < as.numeric(dbR2[,1])
      dbR2[move[!isNA],] <- dbR2[move[!isNA],c(2,1),drop=FALSE]   #sort they are same genotype           
      #dbR2 <- dbR2[!isNA,] #not removed
      #if(sum(!isNA)==1) dbR2 <- rbind(dbR2) #require conversion if one possible combination
      undbR <- unique(dbR2) #get unique genotypes
      for(j in 1:nrow(undbR)) { #for each unique reference profile
        DBrefGeno <- undbR[j,] #take out alleles of reference
        nU_hp2 <- nU_hp + sum(is.na(DBrefGeno[1])) #add an extra unknown if locus missing (alleles is NA)
        dbind <-  which(dbR2[,1]==DBrefGeno[1] & dbR2[,2]==DBrefGeno[2]) #get index of matching genotypes
        if(length(dbind)==0) {
          dbind <- which(is.na(dbR2[,1])) #index of miss markers
          condRefHp2 <- c(NULL,condRefHp ) #conditional references under hp
          DBrefGeno = NULL
        } else {
          condRefHp2 <- c(DBrefGeno,condRefHp ) #conditional references  under hp
        }
        Evid <- NULL
        for(ss in 1:length(evidlist)) { #for each evidence
          Ei <- evidlist[[ss]]	
          if(length(Ei)==0) Ei = "0" #take 
          if(ss>1 && ss<length(evidlist)) Ei <- c(Ei,"0") #insert zero if more replicates and more to go
          Evid <- c(Evid,Ei)
        } #end for each evidence
        hp0 <- calcQual( Evid,condRefHp2, NULL,nU_hp2, fst, pDhat, pC, popFreqQ[[loc]])
        if(j==1 || fst>0) hd0 <- calcQual( Evid, condRefHp, DBrefGeno ,nU_hd, fst, pDhat, pC, popFreqQ[[loc]])
        LR1[dbind] <- LR1[dbind]*hp0/hd0   #if more alleles than unknown
      }#end for each genotypes
    } #end for each locus
    LR1 = round(log10(LR1),2) #use log10 scale (and round)
    
    #STEP 2) CALCULATE LR vals (Depending on model)
    LR2 = rep(NA,length(LR1))
    if(ITYPE=="QUAN") {   #CONT LR calculation for each reference in table: FOR each database: calculate LR for each samples 
      print(paste0("Calculating quantitative LR for ",nrow(subD)," individual(s) in database ",dsel,"..."))
      #unsubD <- unique( subD ) #get unique values. Not in use
      for(rind in 1:length(indD)) { #for each individual in database
        Dind <- subD[rind,,drop=FALSE] #take out individual
        
        dblocs2 <- dblocs #take out loci which the reference in database have: #Update from v1.9: !NA is removed
        locevalind <- locs_hd%in%dblocs2
        loceval <- locs_hd[locevalind] #locus to evaluate 
        
        #setup for hp:
        #insert Dind to refData         
        if(is.null(refData)) { #
          refData2 <- list()
          condOrder_hp <- 1 #put conditional-index to model  
          condOrder_hd <- 0 #put conditional-index to model 
        } else { 
          refData2 <- refData[loceval] #take out only relevant loci to analyse
          condOrder_hp <- c(set$condOrder_hp,max(set$condOrder_hp)+1) #put conditional-index to model 
          condOrder_hd <- c(set$condOrder_hd,0) #put conditional-index to model 
        }
        nRefs <- length(condOrder_hp) #number of references in refData2
        for(loc in loceval) {
          refData2[[loc]]$ijoisdjskwa <- unlist(strsplit(Dind[ which(loc==colnames(Dind)) ], "/"))  #SOME BUG OBSERVED HERE: insert data into a new ref: name it with a random text to avoid similar with others
          if(is.na(refData2[[loc]]$ijoisdjskwa[1])) refData2[[loc]]$ijoisdjskwa <- numeric() #Update from v1.9: added line
        }
        samples <- lapply( set$samples, function(x) x[loceval] ) #take only relevant mixture data:
        
        logLi_hdeval <- logLi_hd[locevalind] #take out relevant values
        mlefit_hp   <- calcMLE(set$nC_hp+1,samples,popFreqQ[loceval],refData2,condOrder_hp,set$knownref_hp,set$kit,set$DEG,set$BWS,set$FWS, set$threshT, set$prC, set$lambda, set$fst, NULL,NULL,           set$minFreq,set$normalize,  opt$steptol, opt$nDone, opt$delta, opt$difftol, opt$seed, FALSE, set$priorBWS,set$priorFWS, opt$maxThreads, set$adjQbp)
        #mlefit_hp contLikMLE(set$nC_hp+1,samples,popFreqQ[loceval],refData2,condOrder_hp,set$knownref_hp,set$xi,set$prC,opt$nDone,set$threshT,set$fst,set$lambda,delta=opt$delta,pXi=set$pXi,kit=set$kit,verbose=FALSE,difftol=opt$difftol ,                                  xiFW=set$xiFW,pXiFW=set$pXiFW, maxThreads=opt$maxThreads,seed=opt$seed,steptol=opt$steptol)
        if(any(set$fst>0)) { #must calculate Hd once again (assume Rj is known)
          knownref_hd = c(set$knownref_hp,nRefs) #include DB-ref
          mlefit_hdj  <- calcMLE(set$nC_hd,samples,popFreqQ[loceval],refData2,condOrder_hd,    knownref_hd,set$kit,set$DEG,set$BWS,set$FWS, set$threshT, set$prC, set$lambda, set$fst, set$knownRel,set$ibd,set$minFreq,set$normalize,  opt$steptol, opt$nDone, opt$delta, opt$difftol, opt$seed, FALSE, set$priorBWS,set$priorFWS, opt$maxThreads, set$adjQbp)
          #mlefit_hdj<- contLikMLE(set$nC_hd,  samples,popFreqQ[loceval],refData2,condOrder_hd,          nRefs,set$xi,set$prC,opt$nDone,set$threshT,set$fst,set$lambda,delta=opt$delta,pXi=set$pXi,kit=set$kit,verbose=FALSE,difftol=opt$difftol,knownRel=set$knownRel,ibd=set$ibd ,xiFW=set$xiFW,pXiFW=set$pXiFW, maxThreads=opt$maxThreads,seed=opt$seed,steptol=opt$steptol)
          LR2[rind] <- mlefit_hp$fit$loglik - mlefit_hdj$fit$loglik #insert calculated LR adjusted by fst-correction
        } else {
          LR2[rind] <- mlefit_hp$fit$loglik - sum(logLi_hdeval) #insert calculated LR:
        }  
        if(rind%%10==0) print(paste0(round(rind/length(indD)*100),"% finished")) #Update each 10th
      } #end for each individual
      LR2 = round(LR2/log(10),2) #use log10 scale (and round)
    } #end type !QUAL
    LR2[is.na(LR2)] = "-" #updated in v1.9: NA causes error, replaced with "-"
    print(paste0(100,"% finished for database ",dsel))
    DBtab <- rbind(DBtab , cbind(indD,LR2,LR1,macD,nlocs) ) #add to DBtab
  } #end for each databases
  colnames(DBtab) <- c("Referencename","quanLR","qualLR","MAC","nMarkers")
  return(DBtab)
}