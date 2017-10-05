#set directory:
#setwd("/Users/KeithMacBookPro/Dropbox (Personal)/CSUEB/Grants/NIJ 2013 DNA grant/CSUEB Laboratory/PGS Evaluations/likeLTD/Test files")
require(euroformix); sessionInfo()

################
#help functions#
################
readFreq <- function(file) { #import popfrequencies:
 table <- read.table(file,header=TRUE,sep=",")
 locs <- toupper(colnames(table[-1]))
 popFreq <- list()
 for(i in 1:length(locs)) {
  freqs <- table[,i+1]
  popFreq[[i]] <-  table[!is.na(freqs),i+1]
  names(popFreq[[i]]) <-  table[!is.na(freqs),1]
 }
 names(popFreq) <- locs
 return(popFreq)
}

tableReader=function(filename) {
  tab <- read.table(filename,header=TRUE,sep="\t",stringsAsFactors=FALSE)
  tryCatch( {  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=",",stringsAsFactors=FALSE) } ,error=function(e) e) 
  tryCatch( {  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE) } ,error=function(e) e) 
  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE)
  return(tab) #need dataframe to keep allele-names correct!!
 }

sample_tableToList = function(X) {
  cn = colnames(X) #colnames 
  lind = grep("marker",tolower(cn),fixed=TRUE) #locus col-ind
  if(length(lind)==0) lind = grep("loc",tolower(cn),fixed=TRUE) #try another name
  sind = grep("sample",tolower(cn),fixed=TRUE) #sample col-ind
  if(length(sind)>1)  sind = sind[grep("name",tolower(cn[sind]),fixed=TRUE)] #use only sample name
  A_ind = grep("allele",tolower(cn),fixed=TRUE) #allele col-ind
  H_ind = grep("height",tolower(cn),fixed=TRUE) #height col-ind
  ln = unique(toupper(X[,lind])) #locus names: Convert to upper case
  sn = unique(as.character(X[,sind])) #sample names
  I = length(ln)
  Y = list() #insert non-empty characters:
  for(k in 1:length(sn)) { #for each sample in matrix
   Y[[sn[k]]] = list() #one list for each sample
   for(i in 1:I) { #for each locus
     xind = X[,sind]==sn[k] & toupper(X[,lind])==ln[i] #get index in X for given sample and locus
     if(sum(xind)==0) next
     keep <- which(!is.na(X[xind,A_ind]) & X[xind,A_ind]!="")
     if(length(A_ind)>0) Y[[sn[k]]][[ln[i]]]$adata = as.character(X[xind,A_ind][keep])
     if(length(H_ind)>0) Y[[sn[k]]][[ln[i]]]$hdata = as.numeric(as.character(X[xind,H_ind][keep]))
   }   
  }
  names(Y) <- sn
  return(Y)
}

getData <- function(mixData2,refData2,popFreq,stutter) { #Helpfunction to get data to analyse
 locs <- names(popFreq)
 mixData <- lapply(mixData2,function(x) return(x[locs])) #return selected loci
 refData <- list()
 for(loc in locs)  refData[[loc]] <- lapply(refData2,function(x) return(x[[loc]]$adata)) #return selected loci
 Qret <- Qassignate(samples=mixData, popFreq,  refData,incS=stutter,incR=FALSE)
 return(list(samples=mixData,refData=Qret$refData,popFreq=Qret$popFreq))
}
###################################################################

#SCRIPT STARTS HERE:

#Outfile to store results
outf <- "results.txt" #"/Users/KeithMacBookPro/Dropbox (Personal)/CSUEB/Grants/NIJ 2013 DNA grant/CSUEB Laboratory/PGS Evaluations/likeLTD/Test files/combinedResults.csv"


#get popfreq file:
databaseFile = "ESX17_Norway.csv"
popFreq <- readFreq(databaseFile) #import population freqs.
#names(popFreq) #loci to consider

#Get evidences (files)
evidfold <- "spor" #"/Users/KeithMacBookPro/Dropbox (Personal)/CSUEB/Grants/NIJ 2013 DNA grant/CSUEB Laboratory/PGS Evaluations/likeLTD/Test files/Multisample files"
files = list.files(evidfold)

#get references:
refFile <- "refs/refs.txt" #"/Users/KeithMacBookPro/Dropbox (Personal)/CSUEB/Grants/NIJ 2013 DNA grant/CSUEB Laboratory/PGS Evaluations/likeLTD/likeLTD reference files//NIJ B reference datafile.csv"
refData=sample_tableToList(tableReader(refFile)) #load references
rN <- names(refData) #names of references


#Model setup:
kit = "ESX17" #"Identifiler"
threshT = 150 #25 #detection threshold (rfu)
fst = 0.01
nC = 2 #assumed number of contributors

stutter=T##F #consider stutter model?
xi = NULL #stutter parameter: Default is that it is estimated (NULL)
if(!stutter) xi = 0  #otherwise set as 0 => no stutter

dropin=F #consider drop-in model?
pC<-lambda<-0
if(dropin) {
 pC<-0.05 #dropin frequency.
 lambda<-0.01 #dropin peak height param 
}

#optimization and MCMC iterations
nDone=4 #Number of required converged estimates 
niter=2000 #Number of iteration for LR sensitivity


# out object
out = NULL

# Loop over cases
begin=Sys.time() #start timer
for(i in 1:length(files)) {
	
  evidfile = paste0(evidfold,"/",files[i]) #evidence file
  mixData = sample_tableToList( tableReader( evidfile) ) #get sample to analyse
  refData2 <- refData[i] #consider only ref "i" 
  hpcond <- c(1) #Hp condition: ref i is contributor 1
  hdcond <- c(0) #Hd condition: ref i is not-contributor
  knownRefHd <- 1 #condition under Hd that ref i is a known non-contributors

  #plotEPG(Data=mixData,kitname=kit,threshT=threshT,refcond=refData2,showPH=TRUE) #plotting evidence with ref
  dat <- getData(mixData,refData2,popFreq,stutter=stutter) #process data for euroformix calculations

  #Perform calculatations
  #calculate LR (mle based) using continuous model:
  set.seed(1)
  hpfit <- contLikMLE(nC,dat$samples,popFreq=dat$popFreq,threshT=threshT,nDone=nDone,xi=xi,refData=dat$refData,prC=pC,lambda=lambda,fst=fst,condOrder=hpcond)
  #print(hpfit$fit$thetahat2) #estimated parameters
  #validMLEmodel(hpfit,kit=kit,plottitle="Hp") #valid model assumption under hp

  set.seed(1)
  hdfit <- contLikMLE(nC,dat$samples,popFreq=dat$popFreq,threshT=threshT,nDone=nDone,xi=xi,refData=dat$refData,prC=pC,lambda=lambda,fst=fst,condOrder=hdcond,knownRef=knownRefHd)
  #print(hdfit$fit$thetahat2) #estimated parameters
  #validMLEmodel(hdfit,kit=kit,plottitle="Hd") #valid model assumption under hd

  #WoE:
  LRmle <- exp(hpfit$fit$loglik - hdfit$fit$loglik)
  log10LRmle <- (hpfit$fit$loglik - hdfit$fit$loglik)/log(10)#

  #MCMC:
log10LRmleCons <- NA
if(0) { #if you want you can allow this section to calculate your "conservative" LR estimate
  qq0 <- 0.05
  delta=10  #variance of the MCMC proposal distribution
  set.seed(1)
  mcmchp <- contLikMCMC(hpfit,niter=niter,delta = delta)
  #validMCMC(mcmchp)
  set.seed(1)
  mcmchd <- contLikMCMC(hdfit,niter=niter,delta = delta)
  #validMCMC(mcmchd)
  log10LRmleCons <- quantile((mcmchp$postlogL-mcmchd$postlogL)/log(10),qq0)
}
	
  # update out object
  out = rbind(out, c(files[i],log10LRmle,log10LRmleCons))	 #notice on log10 scale = Bins
} #end loop for each evidence
end=Sys.time() #end timer
runtime=difftime(end,begin) #Calculate the total running time:
print(paste0("Routine spent ",round(as.numeric(runtime))," seconds"))

# Export overall results
colnames(out) = c("Filename","log10LR (ML)","log10LR (Cons)")
write.csv(out,outf,row.names=F)
