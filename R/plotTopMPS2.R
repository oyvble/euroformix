#' @title plotTopMPS2
#' @author Oyvind Bleka
#' @description MPS data visualizer (interactive)
#' @details Plots the expected peak heights of the top genotypes. The peak heights for corresponding alleles (one sample) are superimposed.
#' @param MLEobj An object returned from contLikMLE
#' @param DCobj An object returned from devonvolve: Must be run with same object as MLEobj
#' @param AT A detection threshold can be shown in a dashed line in the plot (constant). Possibly a AT[[loc]] list.
#' @param ST A stochastic threshold can be shown in a dashed line in the plot (constant). Possibly a ST[[loc]] list.
#' @param grpsymbol A separator for each allele giving plot grouping. Useful for separating conventional repeat units (RU) and sequence variant.
#' @param locYmax A boolean of whether Y-axis should be same for all markers (FALSE) or not (TRUE this is default)
#' @param options A list of possible plot configurations. See comments below
#' @return sub A plotly widget
#' @export

plotTopMPS2 = function(MLEobj,DCobj=NULL,AT=NULL,ST=NULL,grpsymbol="_",locYmax=TRUE,options=NULL) {
 require(plotly) #required package
 Qallele = "99"
 if(is.null(options$h0)) { h0 = 300 } else { h0 = options$h0 } # 5500/nrows #standard height for each dye (depends on number of rows? No)
 if(is.null(options$w0)) { w0 = 1800 } else { w0 = options$w0 } # standard witdh when printing plot
 if(is.null(options$marg0)) { marg0 = 0.015 } else { marg0 = options$marg0 } #Margin between subplots
 if(is.null(options$txtsize0)) { txtsize0 = 12 } else { txtsize0 = options$txtsize0 } #text size for alleles
 if(is.null(options$locsize0)) { locsize0 = 20 } else { locsize0 = options$locsize0 } #text size for loci
 if(is.null(options$minY)) { minY = 100 } else { minY = options$minY } #default minimum Y-axis length
 if(is.null(options$ymaxscale)) { ymaxscale = 1.06 } else { ymaxscale = options$ymaxscale } #y-axis scaling to the locus name positions
 if(is.null(options$grptype)) { grptype="group" } else { grptype = options$grptype }#,"stack" "group" is default 

 mixData=MLEobj$model$samples #copy
 refData=MLEobj$model$refData #copy

 sn = names(mixData) #get samples names
 nS = length(sn) #number of replicates
 locs = names(mixData[[1]]) #get locus names
 nL = length(locs)

 nrefs = 0 
 if(!is.null(refData)) {
  refn = names(refData[[1]])
  refn = refn[MLEobj$model$condOrder>0] #ref names conditioned on
  nrefs = length(refn)
 }
 if(is.null(DCobj)) DCobj <- deconvolve(MLEobj,maxlist=1) #get top candidate profiles
 #extract info from DC (deconvolution) object
 topG <- sapply(DCobj$rankGi,function(x) x[1,-ncol(x),drop=F])
 if(is.null(nrow(topG))) topG <- t(topG) #consider as matrix
 pG <- as.numeric(sapply(DCobj$rankGi,function(x) x[1,ncol(x)])) #probabilities
 names(pG) = toupper(colnames(topG)) #assign loci names

 #estimates:
 thhat <- MLEobj$fit$thetahat2 #get estimates
 nC <- MLEobj$model$nC #number of contributors
 mx <- thhat[1:nC]   
 mu <- thhat[nC+1]   
 sigma <- thhat[nC+2]   
 beta <- 1
 kitinfo = NULL
 if(!is.null(MLEobj$model$kit)) {
  beta <- thhat[nC+3]   
  kitinfo = euroformix::getKit(MLEobj$model$kit) #names(kitinfo)
 }
#  Ccols <- c("cyan","green","coral","gold","hotpink","darkorange","lightgoldenrod") #c("black","gray","brown","darkorange") #contributor cols
 Ccols <- c("blue","green","coral","gold","hotpink","darkorange","lightgoldenrod") #c("black","gray","brown","darkorange") #contributor cols


 df = numeric() #store data: (sample,marker,allele,height)
 for(ss in sn) { #create a seperate EPG plot for each samples
  #locs = names(mixData[[ss]])
  for(loc in locs) {
   #loc=locs[1]

    edat = mixData[[ss]][[loc ]] #get evid data   
    if(nrefs==0) {
      rdat = NULL
    } else {
      rdat = refData[[loc]][refn] #get ref data (list) for loci given first
    }

    if(is.null(edat) && is.null(rdat)  ) next #skip if no data (evid or ref)

    av = edat$adata
    hv = edat$hdata
    av2 = unique(unlist(rdat))

    G = topG[,colnames(topG)==toupper(loc)] #get genotypes
    av2 = unique( c(av2,unlist(strsplit(G,"/"))) ) #get all alleles 
    adda = av2[!av2%in%av] 

    #add missing:
    if(length(adda)>0) {
     av = c(av,adda)
     hv  =  c(hv, rep(0,length(adda)) )
    }
 
    if(length(av)==0) { #add dummy variables if no alleles
      av <- ""
      hv <- 0 
      av1 <- av2 <- rep("",length(av))
    } else { #otherwise if observed:
      tmp = strsplit(av,grpsymbol) 
      av1 = sapply(tmp,function(x) x[1]) #extract RU allele
      av2 = sapply(tmp,function(x) paste0(x[-1],collapse=grpsymbol)) #collapse other levels if several

     #sort alleles increasingly (handle strings)
      suppressWarnings({ av1n = as.numeric(av1)})
      if(any(is.na(av1n))) av1n = av1
      ord = order(av1n) 
      av = av[ord]
      hv = hv[ord]
      av1 = av1[ord]
      av2 = av2[ord]
    }

    #ref text under each allele (follow original av)
    reftxt <- rep("",length(av))
    if(nrefs>0) {
     for(rr in 1:nrefs) { #for each ref
      indadd = which(av%in%unlist(rdat[[rr]])) #index of alleles to add to text
      hasprevval = indadd[nchar(reftxt[indadd])>0] #indice to add backslash (sharing alleles)
      reftxt[ hasprevval ] = paste0(reftxt[ hasprevval ],"/")      
      reftxt[indadd] = paste0( reftxt[indadd], rr)
     }
    }

    if(!is.null(kitinfo)) {
     tmp = kitinfo[toupper(kitinfo$Marker)==loc,]
     ind = match(av1,tmp$Allele) #get index to extract
     bv = tmp$Size[ind]  #get sizes directly from lookup
     isna = which(is.na(ind)) #which alleles are missing?
	if(length(isna)>0) avuse = as.numeric(tmp$Allele) #alleles available (called only once)

     for(missind in isna) {#impute missing bp:
      if( av[missind]==Qallele) {
       bv[missind] = max(tmp$Size)
      } else {
       newa = as.numeric(av1[missind ])
       impuse = which.min(abs(newa - avuse )) #index of closest allele 
       newa2 =  avuse[impuse] #closest allele 
       diff = newa - newa2
       bpadd1 = floor(diff)*tmp$Repeat[impuse]  #integer add
       bpadd2 = (diff-floor(diff))*10  #decimal add
       bv[missind] = tmp$Size[impuse] + bpadd1 + bpadd2  #estimate bp to insert
      }
     }
    }

    #GET EXPECTation (and std) of PH
    EYv = rep("",length(av))
    SDv = rep("",length(av))

    for (aa in av) { # Loop over all alleles in locus
       ind = which(av==aa)
       contr <- sapply(strsplit(G,"/"),function(x) sum(x%in%aa)) #get contribution
       if(!is.null(kitinfo)) {
        EY <- cumsum(contr*mx*mu*beta^((bv[ind]-125)/100)) #expected peak heights for each contributors
       } else {
        EY <- cumsum(contr*mx*mu) #expected peak heights for each contributors
       }
       EYv[ind] <- paste0(round(EY),collapse="/") #notice rounding to whole number!!
       SDv[ind] <- paste0(round(EY*sigma),collapse="/") #notice rounding to whole number!!
    }
    df = rbind(df, cbind(ss,loc,av,hv,reftxt,av1,av2,EYv,SDv) )
  } #end for each loci
 } #end for each samples
 df = data.frame(Sample=df[,1],Marker=df[,2],Allele=df[,3],Height=as.numeric(df[,4]),reftxt=df[,5],Allele1=df[,6], Allele2=df[,7],EXP=df[,8],SD=df[,9],stringsAsFactors=FALSE)
 #df[,-c(3,7)]
 #colnames(df)
 EYmax <- max(as.numeric(unlist( strsplit(df$EXP,"/") )))
 ymax1 <- ymaxscale*max(minY,EYmax,df$Height) #global max y

#GRAPHICAL SETUP BASED ON SELECTED KIT:
ncols = 5 #number of locs per row (depend on amax)?
maxA = max(aggregate(df$Height,by=list(df$Sample,df$Marker),length)$x)
if(maxA<=2 && nL>40) ncols=10 #Number of cols should depend on allele outcome (SNP vs STR)
if(!is.null(options$ncols)) ncols=options$ncols #use number of column specified in options 
h1 = h0*ncols #standard height for each graph (DEPEND ON NUMBER COLS)	
nrows0=ceiling((nL+1)/ncols) #number of rows to use (use 10 per column)
if(nrefs>0) nrows0=ceiling((nL+1)/ncols) #number of rows to use (use 10 per column)

hline <- function(y = 0, color = "black",xr=0:1) {
  list(
    type = "line", 
    x0 = xr[1], 
    x1 = xr[2], 
    y0 = y, 
    y1 = y, 
    line = list(color = color,dash = 'dot',width=2)
  )
}
  
getTcol <- function(color, deg = 50) {
  rgb.val <- col2rgb(color)
  return( rgb(rgb.val[1], rgb.val[2], rgb.val[3],max = 255, alpha = (1-deg)*255))
}
transdeg = .8

#SEPARATE PLOTS

for(ss in sn) {
 #ss = sn[1]
 #locs = locs[1:30]
 plist = list() #create plot object for each marker
 for(loc in locs) {
#loc  = locs[21]
   AT1 <- AT #temporary on analytical threshold
   ST1 <- ST #temporary on stochastic threshold
   if(!is.null(AT) && is.list(AT) ) AT1 = AT[[dye]] #ignores dye if not found
   if(!is.null(ST) && is.list(ST) ) ST1 = ST[[dye]] #ignores dye if not found

   dfs = df[df$Sample==ss & df$Marker%in%loc,] #extract subset 
   dfs$Allele = as.character(dfs$Allele)
   dfs$Allele1 = as.character(dfs$Allele1)
   dfs$Allele2 = as.character(dfs$Allele2)
  
   nA = length(dfs$Allele)
   av1 = unique(dfs$Allele1) #get unique alleles
   nA1 = length(av1) #number of unique
   xpos = 0:(nA-1) #position of alleles (standard)
   reptab = table(dfs$Allele1)
   nR = max(reptab) #number of layers
   repcols = rep("black",nR) #gray.colors(nR, start = 0.3, end = 0.9) #color level will be adjust regarding #layers

   repcol = rep(NA,nA)#get layer index
   for(aa in av1) {
    ind = which(dfs$Allele1==aa)
    repcol[ind] = 1:length(ind)
   }
   xpos1 = rep(NA,nA)
   for(i in 1:nA) xpos1[i] = which(dfs$Allele1[i]==av1)-1
   xpos0 = 0:(nA1-1)
   atxtL = nchar(av1) #get allele length

   p = plotly::plot_ly(dfs,height=h1,showlegend = TRUE,colors=repcols[1:nR] )
   p = plotly::add_trace(p,type = "bar", x = xpos1,y=~Height,name=repcol,hoverinfo="y+text",hoverlabel=list(font=list(size=12),namelength=1000),text =~Allele,color=as.factor(repcol))

   #ADD EXPECTATIONS
   shifts = seq(-0.25,0.25,l=nR)
   if(nR==1) shifts = 0

   Ev = strsplit(dfs$EXP,"/") #extract
   shapeList = list()  #add opacity shapes
   cc = 1 #counter
   for(i in 1:length(Ev)) {
    x0 = xpos1[i] #get allele position (centered)
    x1 = x0+shifts[repcol[i]] #layer decides layer pos
    Ev2 = c(0,Ev[[i]])
    for(j in 1:length(Ev[[i]])) { #for each contributor 
     if(Ev2[j+1]==Ev2[j]) next #skip if equal
      shapeList[[cc]] = list(type = "rect",fillcolor = Ccols[j], line = list(color = Ccols[j],width=0.1), opacity = transdeg,x0 =x1-1/(2*(nR+1)), x1 = x1+1/(2*(nR+1)),y0 = Ev2[j], y1 = Ev2[j+1])#,xref="x", yref = "y")
      cc = cc + 1
    }
   }

   #Add threshold lines to shapes
   if(!is.null(AT))  {
      shapeList[[cc]] = hline(AT1,xr=c(-0.5,nA1-0.5))
      cc = cc + 1
   }
   if(!is.null(ST)) shapeList[[cc]] = hline(ST1,xr=c(-0.5,nA1-0.5))
   if(locYmax)  ymax1 = ymaxscale*max(minY,AT1,ST1,as.numeric(unlist(Ev)),dfs$Height)  #get max 

   dye2 = "black" 
    prob = pG[names(pG)==toupper(loc)]
    if( length(prob)>0 ) {
     if( prob>=.95 ) {
      dye2 = "forestgreen"
     } else if(prob>=.9) {
      dye2 = "orange"
     } else {
      dye2 = "red"
     }
     p = plotly::add_annotations(p, x=(nA1-1)/2 ,y=ymax1,text=loc,showarrow=FALSE,font = list(color = dye2,family = 'Gravitas One',size = locsize0))  #ADD LOCI NAME
    }
   if(max(atxtL)<=5) p = plotly::add_annotations(p, x=xpos0 ,y=rep(0,nA1),text=av1 ,showarrow=FALSE,font = list(color = 1,family = 'sans serif',size = txtsize0),yshift=-10)  #ADD ALLELE NAMES
   #Rotate alleles if many alleles? Using tickangle: layout(xaxis=list(tickangle=-45))
   #p = plotly::add_annotations(p, x=(nA1-1)/2 ,y=ymax1,text=loc,showarrow=FALSE,font = list(color = 1,family = 'Gravitas One',size = locsize0))  #ADD LOCI NAMES 

   if(nrefs>0) { #need to take care of different layers (shift on x-axis) + missing PHs
     refuse = which(dfs$reftxt!="")  #!duplicated(dfs$Allele1) #only put reference under relevant alleles
     for(rr in refuse) { #for each refs (some are missing PH)
      x0 = xpos1[rr] #get allele position (centered)
      x1 = x0+shifts[repcol[rr]] #layer decides layer pos
      p = plotly::add_annotations(p, x=x1,y=0,text=dfs$reftxt[rr],showarrow=FALSE,font = list(color = repcols[repcol[rr]],family = 'sans serif',size = txtsize0),yshift=-25)  #ADD ALLELE NAMES
      if(dfs$Height[rr]==0)  p = plotly::add_trace(p,type = "scatter",mode="markers", x = x1,y=0,name=repcol[rr],hoverinfo="y+text",hoverlabel=list(font=list(size=12),namelength=1000),text=dfs$Allele[rr],color=as.factor(repcol[rr]))  #add a point to missing alleles (to get hovering)
     } #end for each refuse
   } #end if references
   p = plotly::layout(p,xaxis = list(showticklabels = FALSE,title = ""),yaxis=list(range=c(0,ymax1), showline = TRUE,title = ""),shapes=shapeList)
   plist[[loc]] <- p
 } #end for each loc
# subplot(plist, nrows = nrows0, shareX = FALSE, shareY = FALSE,margin=marg0,titleY= TRUE)%>%plotly::layout(title=ss,barmode = grptype)%>%hide_legend()

 #In last plot we will show the Mx and the contributors 
 p = plotly::plot_ly(x = xpos,y=rep(0,length(xpos)),height=h1,showlegend = FALSE,type="scatter",mode="markers")
 p = plotly::add_annotations(p, x=tail(xpos,1),y= c(ymax1-ymax1/6*(1:length(mx))),text= paste0("Contr.C",1:length(mx),"(",Ccols[1:length(mx)],")=",signif(mx,3)),showarrow=FALSE,font = list(family = 'sans serif',size = 15),xshift=0,xanchor = 'right') 

 if(nrefs>0) p = plotly::add_annotations(p, x=0,y=c(ymax1-ymax1/6*(1:nrefs)-ymax1/12),text= paste0("Label ",1:nrefs,": ",refn),showarrow=FALSE,font = list(family = 'sans serif',size = 15),xshift=0,xanchor = 'left')  #ADD ALLELE NAMES
 p = plotly::layout(p,xaxis = list(showline=FALSE, showticklabels = FALSE,title = ""),yaxis=list(range=c(0,ymax1), showline=FALSE,showticklabels = FALSE,title = ""))#,colorway =dye2) 
 plist[[nL+1]] = p

 sub = plotly::subplot(plist, nrows = nrows0, shareX = FALSE, shareY = FALSE,margin=marg0,titleY= TRUE)%>%plotly::layout(title=ss,barmode = grptype)%>%plotly::hide_legend()%>%plotly::config(scrollZoom=TRUE, displaylogo=FALSE,modeBarButtonsToRemove=c("lasso2d","select2d","hoverClosestCartesian","hoverCompareCartesian","toggleSpikelines"),toImageButtonOptions=list(width=w0))

 print(sub)
 } #end for each samples
 return(sub) #return last created

} #end function
