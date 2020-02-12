#' @title plotMPS2
#' @author Oyvind Bleka
#' @description MPS data visualizer (interactive)
#' @details Plots intensities with corresponding allele variant for one sample. Does not yet handle replicates. Can handle RU grouping with separator grpsymbol.
#' @param mixData List of mixData[[ss]][[loc]] =list(adata,hdata), with samplenames ss, loci names loc, allele vector adata (can be strings or numeric), intensity vector hdata (must be numeric) 
#' @param refData List of refData[[rr]][[loc]] or refData[[loc]][[rr]] to label references (flexible). Visualizer will show dropout alleles. 
#' @param AT A detection threshold can be shown in a dashed line in the plot (constant). Possibly a AT[[loc]] list.
#' @param ST A stochastic threshold can be shown in a dashed line in the plot (constant). Possibly a ST[[loc]] list.
#' @param grpsymbol A separator for each allele giving plot grouping. Useful for separating conventional repeat units (RU) and sequence variant.
#' @param locYmax A boolean of whether Y-axis should be same for all markers (FALSE) or not (TRUE this is default)
#' @param options A list of possible plot configurations. See comments below
#' @return sub A plotly widget
#' @export

plotMPS2 = function(mixData,refData=NULL,AT=NULL,ST=NULL,grpsymbol="_",locYmax=TRUE,options=NULL) {
 require(plotly) #required package
 if(is.null(options$h0)) { h0 = 300 } else { h0 = options$h0 } # 5500/nrows #standard height for each dye (depends on number of rows? No)
 if(is.null(options$w0)) { w0 = 1800 } else { w0 = options$w0 } # standard witdh when printing plot
 if(is.null(options$marg0)) { marg0 = 0.015 } else { marg0 = options$marg0 } #Margin between subplots
 if(is.null(options$txtsize0)) { txtsize0 = 12 } else { txtsize0 = options$txtsize0 } #text size for alleles
 if(is.null(options$locsize0)) { locsize0 = 20 } else { locsize0 = options$locsize0 } #text size for loci
 if(is.null(options$minY)) { minY = 100 } else { minY = options$minY } #default minimum Y-axis length
 if(is.null(options$ymaxscale)) { ymaxscale = 1.06 } else { ymaxscale = options$ymaxscale } #y-axis scaling to the locus name positions
 if(is.null(options$grptype)) { grptype="group" } else { grptype = options$grptype }#,"stack" "group" is default 

 sn = names(mixData) #get samples names
 nS = length(sn) #number of replicates
 locs = names(mixData[[1]]) #get locus names
 nL = length(locs)

 locFirst = FALSE #boolean of refData[[loc]][[rr]] (or refData[[rr]][[loc]])
 nrefs = 0 
 if(!is.null(refData)) {
  refn =  names(refData) #default structure (same as old) 
  if(any(refn%in%locs)) { #convert data structure of reference
    refn = names(refData[[1]])
    locFirst = TRUE
  }
  nrefs = length(refn)
 }

 df = numeric() #store data: (sample,marker,allele,height)
 for(ss in sn) { #create a seperate EPG plot for each samples
  #locs = names(mixData[[ss]])
  for(loc in locs) {
   #loc=locs[1]

    edat = mixData[[ss]][[loc ]] #get evid data   
    if(is.null(refData)) {
      rdat = NULL
    } else {
      if(locFirst) rdat = refData[[loc]] #get ref data (list) #get ref data (list)      
      if(!locFirst) rdat = lapply(refData,function(x) x[[loc]]) #get ref data (list)      
    }

    if(is.null(edat) && is.null(rdat)  ) next #skip if no data (evid or ref)

    av = edat$adata
    hv = edat$hdata
    av2 = unique(unlist(rdat))
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
      av1 = sapply(tmp,function(x) x[1])
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
  df = rbind(df, cbind(ss,loc,av,hv,reftxt,av1,av2) )
  } #end for each loci
 } #end for each samples
df = data.frame(Sample=df[,1],Marker=df[,2],Allele=df[,3],Height=as.numeric(df[,4]),reftxt=df[,5],Allele1=df[,6], Allele2=df[,7],stringsAsFactors=FALSE)
#df[,-c(3,7)]
#colnames(df)
ymax1 <- ymaxscale*max(minY,df$Height) #global max y

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
repcols = c("black","red","blue","green","orange","purple") 

#SEPARATE PLOTS FOR EACH SINGLE SAMPLES

for(ss in sn) {
 #ss = sn[1]
 #locs = locs[1:30]
 plist = list() #create plot object for each marker
 for(loc in locs) {
#loc  = locs[13]
   AT1 <- AT #temporary on analytical threshold
   ST1 <- ST #temporary on stochastic threshold
   if(!is.null(AT) && is.list(AT) ) AT1 = AT[[dye]] #ignores dye if not found
   if(!is.null(ST) && is.list(ST) ) ST1 = ST[[dye]] #ignores dye if not found

   dfs = df[df$Sample==ss & df$Marker%in%loc,] #extract subset 
   dfs$Allele = as.character(dfs$Allele)
   dfs$Allele1 = as.character(dfs$Allele1)
   dfs$Allele2 = as.character(dfs$Allele2)
   if(locYmax)  ymax1 = ymaxscale*max(minY,AT1,ST1,dfs$Height)  #get max 

   nA = length(dfs$Allele)
   av1 = unique(dfs$Allele1) #get unique alleles
   nA1 = length(av1) #number of unique alleles
   xpos = 0:(nA-1) #position of alleles (standard)
   reptab = table(dfs$Allele1)
   nR = max(reptab) #number of layers

   repcol = rep(NA,nA)#get layer index
   for(aa in av1) {
    ind = which(dfs$Allele1==aa)
    repcol[ind] = 1:length(ind)
   }
   xpos1 = rep(NA,nA)
   for(i in 1:nA) xpos1[i] = which(dfs$Allele1[i]==av1)-1
   xpos0 = 0:(nA1-1)
   atxtL = nchar(av1) #get allele length

   p = plotly::plot_ly(dfs,height=h1,showlegend = TRUE,colors=repcols[1:nR] )%>%plotly::layout(xaxis = list(showticklabels = FALSE,title = ""),yaxis=list(range=c(0,ymax1), showline = TRUE,title = ""))
   if(!is.null(AT))   p = plotly::layout(p,shapes = list(hline(AT, xr=c(-0.5,nA1-0.5) ) ))
   if(!is.null(ST))   p = plotly::layout(p,shapes = list(hline(ST, xr=c(-0.5,nA1-0.5) ) ))
   p = plotly::add_trace(p,type = "bar", x = xpos1,y=~Height,name=repcol,hoverinfo="y+text",hoverlabel=list(font=list(size=12),namelength=1000),text =~Allele,color=as.factor(repcol))

   if(max(atxtL)<=5) p = plotly::add_annotations(p, x=xpos0 ,y=rep(0,nA1),text=av1 ,showarrow=FALSE,font = list(color = 1,family = 'sans serif',size = txtsize0),yshift=-10)  #ADD ALLELE NAMES
   #Rotate alleles if many alleles? Using tickangle: layout(xaxis=list(tickangle=-45))
   p = plotly::add_annotations(p, x=(nA1-1)/2 ,y=ymax1,text=loc,showarrow=FALSE,font = list(color = 1,family = 'Gravitas One',size = locsize0))  #ADD LOCI NAMES 

   if(nrefs>0) { #need to take care of different layers (shift on x-axis) + missing PHs
     shifts = seq(-0.25,0.25,l=nR)
     if(nR==1) shifts = 0
     refuse = which(dfs$reftxt!="")  #!duplicated(dfs$Allele1) #only put reference under relevant alleles
     for(rr in refuse) { #for each refs (some are missing PH)
      x0 = xpos1[rr] #get allele position (centered)
      x1 = x0+shifts[repcol[rr]] #layer decides layer pos
      p = plotly::add_annotations(p, x=x1,y=0,text=dfs$reftxt[rr],showarrow=FALSE,font = list(color = repcols[repcol[rr]],family = 'sans serif',size = txtsize0),yshift=-25)  #ADD ALLELE NAMES
      if(dfs$Height[rr]==0)  p = plotly::add_trace(p,type = "scatter",mode="markers", x = x1,y=0,name=repcol[rr],hoverinfo="y+text",hoverlabel=list(font=list(size=12),namelength=1000),text=dfs$Allele[rr],color=as.factor(repcol[rr]))  #add a point to missing alleles (to get hovering)
     } #end for each refuse
   } #end if references
   plist[[loc]] <- p
 } #end for each loc

 if(nrefs>0) {  #In last plot we will show the contributors 
  p = plotly::plot_ly(x = xpos,y=rep(0,length(xpos)),height=h1,showlegend = FALSE,mode="scatter",type="scatter")
  p = plotly::add_annotations(p, x=0,y=c(ymax1-ymax1/5*(1:nrefs)),text= paste0("Label ",1:nrefs,": ",refn),showarrow=FALSE,font = list(colors = "black",family = 'sans serif',size = 15),xshift=0,xanchor = 'left')  #ADD ALLELE NAMES
  p = plotly::layout(p,xaxis = list(showline=FALSE, showticklabels = FALSE,title = ""),yaxis=list(range=c(0,ymax1), showline=FALSE,showticklabels = FALSE,title = ""))#,colorway =dye2) 
  plist[[nL+1]] = p
 }
 #sub = subplot(plist, nrows = nrows0, shareX = FALSE, shareY = FALSE,margin=marg0,titleY= TRUE)%>%hide_legend()
 sub = plotly::subplot(plist, nrows = nrows0, shareX = FALSE, shareY = FALSE,margin=marg0,titleY= TRUE)%>%plotly::layout(title=ss,barmode = grptype)%>%plotly::hide_legend()%>%plotly::config(scrollZoom=TRUE, displaylogo=FALSE,modeBarButtonsToRemove=c("lasso2d","select2d","hoverClosestCartesian","hoverCompareCartesian","toggleSpikelines"),toImageButtonOptions=list(width=w0))
 print(sub)
 } #end for each samples
 return(sub) #return last created
} #end function
