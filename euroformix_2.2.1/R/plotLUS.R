#' @title plotLUS
#' @author Oyvind Bleka
#' @description A dataplotter for MPS data
#' @details Plots peak height with corresponding allele variant for one sample
#' @param mixData List of adata- and hdata-elements. adata must be separated with LUSsymbol
#' @param sn The Sample name, will be plotted in title.
#' @param threshT The detection threshold can be shown in gray in the plot.
#' @param refData condition on a list$refname$locname$adata of reference alleles which are labeled in EPG
#' @param LUSsymbol The separator for each allele. Assumed to be "_" separating sequence variant and LUS.
#' @export
plotLUS <- function(mixData,sn="",refData=NULL,threshT=0,LUSsymbol="_") {
 #Data is list with allele and height data. Only one sample!
 locs <- toupper(names(mixData))
 locs <- locs[!grepl("AM",locs)]
 nL <- length(locs) 

 contr = NULL
 if(!is.null(refData)) contr <- names(refData) #get unique contributors

 tab <- cbind(locs,"","")
 for(l in 1:nL) {
  ind <- which(toupper(names(mixData))==locs[l])
  tmp <- mixData[[ind]]$adata
  keep <- !is.na(tmp) | length(tmp)>0 | tmp!=""
  if(!any(keep)) next
 
  av <- mixData[[ind]]$adata[keep] #alleles 
  hv <- as.numeric(mixData[[ind]]$hdata[keep]) #p.h. 
  if(!is.null(refData)) {
   av2 = lapply(refData,function(x) {
      ind = which(toupper(names(x))==locs[l]) #loci is not casesensetive for refs
      if(length(ind)>0) return(x[[ind]]$adata)
   })
   av2 <-  unique( unlist(av2)) #get only unique
   newa <- av2[!av2%in%av] #new variants
   if(length(newa)>0) {
     av <- c(av,newa)
     hv <- c(hv,rep(0,length(newa)))
   }
  }
  tab[l,2] <- paste0(av,collapse="/")
  tab[l,3] <- paste0(hv,collapse="/")
 }

 nR <- 5 #number of rows
 nC <- ceiling((nL+1)/nR) #number of columns
 mat <- t(matrix(1:(nR*nC),nrow=nC,ncol=nR))
 mat[mat>(nL+1)] <- 0 #don't plot

 layout(mat, widths=rep(1, ncol(mat)), heights=rep(1, ncol(mat)))
 par(mar=c(2, 1.5, 1.5, 0), mgp=c(0, 0.5, 0), las=0)
 #for(l in 1:nL) plot(0,0,xlim=0:1,ylim=0:1)
 for(l in 1:nL) {
 #print(tab[l,1])
  av <- unlist(strsplit(tab[l,2],"/"))
  hv <- as.numeric(unlist(strsplit(tab[l,3],"/")))
 
  if(length(av)==0) { #if empty plot nothing  
   plot(0,0,xlim=0:1,ylim=0:1,ty="n",axes=FALSE,xlab="",ylab="")
   text(0.5,0.5,"emtpy")
  } else {
   tmp <-  strsplit(av,LUSsymbol)
   av1 <-  sapply(tmp,function(x) x[1]) #get seq variant
   av2 <-  sapply(tmp,function(x) x[2]) #get LUS variant
   avtab <- cbind(av1,av2) 

   unSEQ <- sort(as.numeric(unique(av1))) #get unique seq vairants
   unLUS <- sort(as.numeric(unique(av2))) #get unique lus vairants

   ymarg = 1.5 #increasing to have space for locus name
   ymax <- max(hv)*ymarg
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

   for(i in 1:nrow(avtab)) { #for each alleles
    a1 <- avtab[i,1]
    a2 <- avtab[i,2]
    seg1 <- seg1init[a1==unSEQ] #what seq segment
    seg2 <- seg2init[a2==unLUS] #what LUS segment
    colind <- which(a2==unLUS) #color of lus
    xpos <- seg1+seg2 #what x start to draw
    rect(xpos-shift,0,xpos-shift+d2,hv[i],col=colind)
   }
   legend("topleft",legend=paste0("LUS=",unLUS),col=1:nZ,pch=15, bty="n",cex=0.9)
  }
  mtext(tab[l,1],3,line=-1,col="black")
 
  if(length(av)==0) next #skip if no data
  if(is.null(refData) || length(refData)==0) next #skip if no references

  #Add ref info:
  av <-  sapply(refData,function(x) {
    ind = which(toupper(names(x))==tab[l,1]) #loci is not casesensetive for refs
    if(length(ind)>0) {
     return(x[[ind]]$adata)
    } else {
     return(rep("",2))
    } 
  })
  cs <-  colnames(av) #get contributors
  csid <- which(contr%in%cs) #get Contr id
 
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
   if(length(seg1)==0 | length(seg2)==0) next
   colind <- which(a2==unLUS) #color of lus
   xpos <- seg1+seg2 #what x start to draw
   ctrs <- paste0(unique(avtab[avtab[,2]==a1 & avtab[,3]==a2,1]),collapse="/") #get contributors
   #axis(1,at=xpos-shift+d2/2,labels=ctrs,tck=0,lty=0,padj=0.8,cex.lab=0.8)
   mtext(ctrs, side=1, line=1, at=xpos-shift+d2/2,cex=0.7,col=colind)
  }
 } #end for each markers
 plot(0,0,ty="n",axes=FALSE,xlab="",ylab="")
 if(!is.null(refData)) legend("topleft",legend=paste0("'",1:length(contr),"\'=",contr),col=1,bty="n",cex=1)
 mtext(sn,3,outer=TRUE,line=-1.5)
} #end function

