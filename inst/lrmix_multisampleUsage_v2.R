require(euroformix); sessionInfo()
setwd("C:\\Users\\oyvbl\\Dropbox\\Forensic\\LRmixStudio")
#rm(list=ls())
#source("lrmix_multisampleUsage.R")
library(forensim);library(euroformix)

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

tableReader2=function(filename) {
  tab <- read.table(filename,header=TRUE,sep="\t",stringsAsFactors=FALSE,row.names=NULL)
  tryCatch( {  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=",",stringsAsFactors=FALSE,row.names=NULL) } ,error=function(e) e) 
  tryCatch( {  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE,row.names=NULL) } ,error=function(e) e) 
  if(ncol(tab)==1) tab <- read.table(filename,header=TRUE,sep=";",stringsAsFactors=FALSE,row.names=NULL)
  return(tab) #need dataframe to keep allele-names correct!!
 }

sample_tableToList = function(X,threshT=NULL) {
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
     if(length(H_ind)>0) { #If peak heights are considered
      PH <- as.numeric(as.character(X[xind,H_ind][keep])) #get the peak heights
      if(!is.null(threshT)) keep = which(PH>=threshT) #keep only alleles above thrshold (if given)
      Y[[sn[k]]][[ln[i]]]$hdata = PH[keep]
     }
     if(length(A_ind)>0) {
      Y[[sn[k]]][[ln[i]]]$adata = as.character(X[xind,A_ind][keep])
     }
   }   
  }
  names(Y) <- sn
  return(Y)
}

getData <- function(mixData2,refData2,popFreq) { #Helpfunction to get data to analyse
 locs <- names(popFreq)
 mixData <- lapply(mixData2,function(x) return(x[locs])) #return selected loci
 refData <- list()
 for(loc in locs)  refData[[loc]] <- lapply(refData2,function(x) return(x[[loc]]$adata)) #return selected loci
 Qret <- Qassignate(samples=mixData, popFreq,  refData,incS=FALSE,incR=FALSE)  #NB: NOTICE THE CHANGE HERE OF inclS=FALSE even for stutter model (this has been updated in v2(
 return(list(samples=mixData,refData=Qret$refData,popFreq=Qret$popFreq))
}

calcLR <- function(pD) {
          LR<-1
          pDvec = rep(pD,nC)
          for(loc in names(dat$popFreq)) { #for each locus
           Ei <- NULL #get evidence
           for(ss in 1:length(dat$samples)) { #fix samples
            if(ss>1) Ei <- c(Ei,0) #seperate with 0  
            adata <-dat$samples[[ss]][[loc]]$adata
            if(length(adata)==0) adata=0 #is empty
            Ei <- c(Ei,adata)
           } 
           rdata <- dat$refData[[loc]] #reference data
           hpval <- likEvid( Ei,T=unlist(rdata),V=NULL,x=nC-1,theta=fst, prDHet=pDvec, prDHom=pDvec^2, prC=pC, freq=dat$popFreq[[loc]])
           hdval <- likEvid( Ei,T=NULL,V=unlist(rdata),x=nC,theta=fst, prDHet=pDvec, prDHom=pDvec^2, prC=pC, freq=dat$popFreq[[loc]])
           LR <- LR*hpval/hdval
         } #end for each markers
         return(LR)
}

###################################################################

#SCRIPT STARTS HERE:

#get popfreq file:
databaseFile = "popFreq.csv" # opt$database
#The allele frequency file
popFreq <- readFreq(databaseFile) #import population freqs
#names(popFreq) #loci to consider

#Get evidences (files)
evidfold <- "Evids" #opt$samples #The folder-name with files including evidence profiles
files = list.files(evidfold)

#get references:
refFile <- "refs.txt"# opt$ref #the file including references
refData=sample_tableToList(tableReader2(refFile)) #load references
rN <- names(refData) #names of references


#Model setup:
kit = "NGM"
threshT = 200 #opt$threshold #25 #detection threshold (rfu)
fst = 0.01
nC = 3# opt$unknowns #assumed number of contributors

dropin=TRUE #opt$doDropin #consider drop-in model?
pC=0
if(dropin) {
 pC=0.05 #dropin frequency.
}

#Outfile to store results
setup <- paste0("T",threshT,"_fst",fst,"_pC",pC,"_C",nC)
outf <- paste0("QualResults_",setup,".csv")

cn=c("EvidFile","POI","LR","Dropout")
out = matrix(nrow=0,ncol=length(cn))
colnames(out) = cn

# Loop over cases
begin=Sys.time() #start timer
for(i in 1:length(files)) { #for each evidence files (may contain replicates)
	
  evidfile = paste0(evidfold,"/",files[i]) #evidence files assumed to be looking in the evidfolder
  mixData = sample_tableToList( X=tableReader2(evidfile),threshT=threshT ) #get sample to analyse. NOTICE THAT THE PEAK HEIGHT THRESHOLD IS GIVEN AS ARGUMENT

  for(j in 1:length(rN)) { #for each reference
   refData2 <- refData[j] #consider only ref "j" as POI 
   hpcond <- c(1) #Hp condition: ref i is contributor 1. This example only consider 1 reference profile. With x reference profiles this must be a x long vector.
   hdcond <- c(0) #Hd condition: ref i is not-contributor .This example only consider 1 reference profile. With x reference profiles this must be a x long vector.
   knownRefHd <- 1 #condition under Hd that ref i is a known non-contributors. This is a vector specifying which of the i-th references that are known non-contributors under hd.

   #plotEPG(Data=mixData,kitname=kit,threshT=threshT,refcond=refData2,showPH=TRUE) #plotting evidence with ref
   dat <- getData(mixData,refData2,popFreq) #process data for euroformix calculations (NOTICE THE CHANGE HERE OF NOT INCLUDING STUTTERS)
   nS = length(dat$samples) #number of samples

   #Perform calculatations
   set.seed(1)
   totAv <- sapply(dat$samples, function(x) sum(sapply(x,function(y) length(y$adata)))) #get number of alleles
   refData3 <- list()
   for(loc in names(dat$popFreq)) refData3[[loc]] <- lapply(refData2,function(x) x[[loc]]$adata) #get format for simDOdistr

   dropqq <- c(0.05,0.95) #quantiles to estimate
   totA=floor(mean(totAv)) #round down (this is what LRmix Studio)
   niter = 1e4 #required number of samples

   #PERFORM MC simulations:
   dihp <- simDOdistr(totA=totA,nC=nC,popFreq,refData=refData3,minS=niter, prC=pC,M=2000) #consider only model under Hd
	tmpHp = quantile(dihp,dropqq)
     print("Hp quantiles"); print(tmpHp);
   dihd <- simDOdistr(totA=totA,nC=nC,popFreq,refData=NULL,minS=niter, prC=pC,M=2000) #consider only model under Hd
     tmpHd = quantile(dihd,dropqq)
     print("Hd quantiles"); print(tmpHd)
 
   div <- c(tmpHp ,tmpHd )
   LRmc <- Vectorize(calcLR)(div)
   dropind = which.min(LRmc)
   dhat = div[dropind]
   LR <- LRmc[dropind] #get conservative LR in LRmix

   # update out object
   out = rbind(out, c(files[i],rN[j],LR,signif(dhat,3)))	
 
   # Export overall results
   write.table(out,file=outf,row.names=FALSE)
  }
} #end loop for each evidence
end=Sys.time() #end timer
runtime=difftime(end,begin) #Calculate the total running time:
paste("Time taken: ", sprintf("%.2fmin", runtime))
