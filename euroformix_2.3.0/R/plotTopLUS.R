#' @title plotTopLUS
#' @author Oyvind Bleka
#' @description A dataplotter for MPS data with LUS format
#' @details  Plots the expected coverage reads of the top genotypes. The coverage reads for corresponding alleles (one sample) are superimposed.
#' @param MLEobj An object returned from contLikMLE
#' @param DCobj An object returned from devonvolve: Must be run with same object as MLEobj
#' @param threshT The detection threshold can be shown in gray in the plot.
#' @param LUSsymbol The separator for each allele. Assumed to be "_" separating sequence variant and LUS.
#' @export
plotTopLUS <- function(MLEobj,DCobj=NULL,threshT=0,LUSsymbol="_") {
  if(is.null(DCobj)) DCobj <- deconvolve(MLEobj,maxlist=1) #get top candidate profiles
  #extract info from DC (deconvolution) object
  topG <- sapply(DCobj$rankGi,function(x) x[1,-ncol(x),drop=F])
  if(is.null(nrow(topG))) topG <- t(topG) #consider as matrix
  pG <- as.numeric(sapply(DCobj$rankGi,function(x) x[1,ncol(x)])) #probabilities

  #estimates:
  thhat <- MLEobj$fit$thetahat2 #get estimates
  nContr <- MLEobj$model$nC #number of contributors
  mx <- thhat[1:nContr]   
  mu <- thhat[nContr+1]   
  beta <- 1
  if(!is.null(MLEobj$model$kit)) beta <- thhat[nContr+3]   
  Ccols <- c("cyan","green","coral","gold","hotpink","darkorange","lightgoldenrod") #c("black","gray","brown","darkorange") #contributor cols. MAX IS 7

  #fix order prior:
  kit <- NULL
  if(!is.null(MLEobj$model$kit)) kit <- getKit(MLEobj$model$kit)

  Data <- MLEobj$model$samples #get sample profiles
  refcond <- MLEobj$model$refData #get ref profiles (always 1.contributor)
  sn <- names(Data) #sample names
  refcond <- names(refcond[[1]])[MLEobj$model$condOrder>0] #ref names considered
  refData =  lapply(MLEobj$model$refData,function(x) x[refcond]) #get those to condition on

#Data is list with allele and height data. 
nS = length(Data)
for(ss in names(Data)) { #Create several plots if multiple samples
#  ss = names(Data)[1]
   mixData = Data[[ss]]
   locs <- toupper(names(mixData))
   nL <- length(locs) 

  tab <- cbind(locs,"","")
 for(l in 1:nL) {
  ind <- which(toupper(names(mixData))==locs[l])
  tmp <- mixData[[ind]]$adata
  keep <- !is.na(tmp) | length(tmp)>0 | tmp!=""
  if(!any(keep)) next
 
  av <- mixData[[ind]]$adata[keep] #alleles 
  hv <- as.numeric(mixData[[ind]]$hdata[keep]) #p.h. 
  av2 = unlist(strsplit(topG[,colnames(topG)==toupper(locs[l])],"/"))
  if(!is.null(refData)) {
   av2 <-  unique(av2,unlist(refData[[locs[l]]])) #get only unique
  }
  newa <- av2[!av2%in%av] #new variants
  if(length(newa)>0) {
    av <- c(av,newa)
    hv <- c(hv,rep(0,length(newa)))
  }  
  tab[l,2] <- paste0(av,collapse="/")
  tab[l,3] <- paste0(hv,collapse="/")
 } #end for each locus

 nR <- 5 #number of rows (is this maximum?)
 nC <- ceiling((nL+1)/nR) #number of columns
 mat <- t(matrix(1:(nR*nC),nrow=nC,ncol=nR)) #Get layout of markers
 mat[mat>(nL+1)] <- 0 #don't plot

 #CREATING PLOT HERE
 dev.new(width=25, height=10) 
 layout(mat, widths=rep(1, ncol(mat)), heights=rep(1, ncol(mat)))
 par(mar=c(2, 1.5, 1.5, 0), mgp=c(0, 0.5, 0), las=0)
 #for(l in 1:nL) plot(0,0,xlim=0:1,ylim=0:1)
 for(l in 1:nL) {
 #print(tab[l,1])
  loc = toupper(tab[l,1]) #locs[l] #get locus
  av <- unlist(strsplit(tab[l,2],"/"))
  hv <- as.numeric(unlist(strsplit(tab[l,3],"/")))
  if(length(av)==0) next
  isdupl = duplicated(av) #consider only unique alleles
  av = av[!isdupl]
  hv = hv[!isdupl]

 #PLOT ALL MARKERS
   tmp <-  strsplit(av,LUSsymbol)
   av1 <-  sapply(tmp,function(x) x[1]) #get seq variant
   av2 <-  sapply(tmp,function(x) x[2]) #get LUS variant
   avtab <- cbind(av1,av2) 
   av1 = av1[!is.na(av1)]
   av2 = av2[!is.na(av2)]
   unSEQ <- sort(as.numeric(unique(av1))) #get unique seq vairants
   unLUS <- sort(as.numeric(unique(av2))) #get unique lus vairants

   bpVec = rep(125,length(unSEQ)) #default values
   if(!is.null(kit)) {
    subtab = kit[toupper(kit$Marker)==loc,,drop=FALSE]
    for(a in unSEQ) {
      bpind = which.min( abs(as.numeric(subtab$Allele)-a))
      bpVec[which(unSEQ==a)] = subtab$Size[bpind]
    }
   }

   #INCLUDE MODEL INFO
   G = topG[,colnames(topG)==loc]
   EYlist = list()
   for (allel in av) { # Loop over all peaks.
     aind = which(allel==av) #get index
     contr <- sapply(strsplit(G,"/"),function(x) sum(x%in%allel)) #get contribution on the allele
     EYlist[[aind]] <- c(0,cumsum(contr*mx*mu*beta^((bpVec[av1[aind]==unSEQ]-125)/100))) #expected peak heights
   } 
   ymarg = 1.5 #increasing to have space for locus name
   ymax <- max(hv,unlist(EYlist))*ymarg
   xr <- 1:length(unSEQ) #number of seq alleles
   zr <- 1:length(unLUS) #number of LUS alleles
   nX <- length(xr) #number of seq alleles
   nZ <- length(zr) #number of LUS alleles

   d1 <- 0.5 #width of a seq-rectangle
   d2 <- d1/nZ #width of LUS rectangle
   shift <- 0.5 #=(xr[2]-xr[1])/2 shift should be middle of grid size

   seg1init <- (1:nX) - d1/2 #init index of seq variants
   seg2init <- cumsum(c(0,rep(d2,nZ-1))) #init index of LUS variants

   plot(0,0,xlim=c(-1,nX),ylim=c(0,ymax),ty="n",axes=FALSE,xlab="",ylab="")
   lines(c(-1,nX),c(0,0),lty=1)
   if(threshT>0) lines(c(0,nX),rep(threshT,2),lty=1,col="gray") 
   axis(1,at=xr-shift,labels=unSEQ, tck=0,lty=0,padj=-0.5,cex.lab=1.2)
   axis(2)#,padj=1.5)

   graycol = gray.colors(nZ,start=0.05,end=0.95) 
   for(i in 1:nrow(avtab)) { #for each alleles
    a1 <- avtab[i,1]
    a2 <- avtab[i,2]
    seg1 <- seg1init[a1==unSEQ] #what seq segment
    seg2 <- seg2init[a2==unLUS] #what LUS segment
    if(is.na(a2)) seg2 = mean(seg2init) #use the middle
    colind <- graycol[which(a2==unLUS)] #color of lus
    xpos <- seg1+seg2 #what x start to draw
    if(!is.na(a2)) rect(xpos-shift,0,xpos-shift+d2,hv[i],col=colind)

    #INCLUDE MODEL INFO
    EY = EYlist[[i]] 
    for(j in 1:nContr) { #for each contributor:
      rect(xpos-shift+d2/4,EY[j],xpos-shift+d2*3/4,EY[j+1],col=adjustcolor(Ccols[j],alpha.f=0.8) ) 
    }
   } #end for each alleles

   legend("topleft",legend=paste0("LUS=",unLUS),col=graycol,pch=15, bty="n",cex=0.9)

   if(l==1) {
     if(!is.null(refData)) legend("topright",legend=paste0("'",1:length(refcond),"\'=",refcond),col=1,bty="n",cex=1)
     legend("bottomleft",legend=paste0("Contr.",1:nContr," = ",signif(mx,2)),bty="n",pch=15,col=Ccols[1:nContr])
   }
  prob <- pG[match(loc,colnames(topG))] #get posterior probabilities
  colcode = "red"
  if(prob>0.90) colcode <- "orange"
  if(prob>0.95) colcode <- "forestgreen" #codes
  mtext(loc,3,line=-1,col=colcode )
 
  if(is.null(refData) || length(refData)==0 || length(refData[[loc]])==0) next #skip if no references

  #Add ref info:
  av <-  sapply(refData[[loc]],function(x) {
    if(length(x)>0) {
     return(x)
    } else {
     return(rep("",2))
    } 
  })
  cs <-  colnames(av) #get contributors
  csid <- which(refcond%in%cs) #get Contr id
 
  emt <- av=="" #check if empty
  csid2 <- col(av)[!emt] #belonging contributor
  av <- av[!emt] 
  tmp <-  strsplit(av,LUSsymbol)
  av1 <-  sapply(tmp,function(x) x[1]) #get seq variant
  av2 <-  sapply(tmp,function(x) x[2]) #get LUS variant
  avtab <- unique(cbind(csid2,av1,av2) ) #consider only uniques
  colnames(avtab) <- c("Contr","SEQ","LUS")
  avtab2 <- unique(avtab[,-1,drop=F]) #look at unique seq/lus independent of contr

  for(i in 1:nrow(avtab2)) { #for each alleles of references
   a1 <- avtab2[i,1]
   a2 <- avtab2[i,2]
   seg1 <- seg1init[a1==unSEQ] #what seq segment
   seg2 <- seg2init[a2==unLUS] #what LUS segment
   if(is.na(seg2)) seg2=0 #allele 99
   if(length(seg1)==0 | length(seg2)==0) next
   xpos <- seg1+seg2 #what x start to draw
   ctrs <- paste0(unique(avtab[avtab[,2]==a1 & avtab[,3]==a2,1]),collapse="/") #get contributors
   if(ctrs=="NA") ctrs <- paste0(unique(avtab[avtab[,2]==a1,1]),collapse="/") #get contributors

   #axis(1,at=xpos-shift+d2/2,labels=ctrs,tck=0,lty=0,padj=0.8,cex.lab=0.8)
   mtext(ctrs, side=1, line=1, at=xpos-shift+d2/2,cex=0.7,col=1)
  }
 } #end for each markers

 #Label the contributor colors:
 mtext(ss,3,outer=TRUE,line=-1.5)
 }#end for each samples
} #end function

