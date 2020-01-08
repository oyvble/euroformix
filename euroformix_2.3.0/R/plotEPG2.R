#' @title plotEPG2
#' @author Oyvind Bleka
#' @description EPG data visualizer (interactive)
#' @details Plots peak height with corresponding allele for sample(s) for a given kit.
#' @param mixData List of mixData[[ss]][[loc]] =list(adata,hdata), with samplenames ss, loci names loc, allele vector adata (can be strings or numeric), intensity vector hdata (must be numeric) 
#' @param refData List of refData[[rr]][[loc]] or refData[[loc]][[rr]] to label references (flexible). Visualizer will show dropout alleles. 
#' @param kitname Short name of kit: See supported kits with getKit()
#' @param AT A detection threshold can be shown in a dashed line in the plot (constant). Possibly a AT[[dye]] list.
#' @param ST A stochastic threshold can be shown in a dashed line in the plot (constant). Possibly a ST[[dye]] list.
#' @param dyeYmax A boolean of whether Y-axis should be same for all markers (FALSE) or not (TRUE this is default)
#' @param plotRepsOnly A boolean of whether only replicate-plot is shown in case of multiple samples (TRUE is default) 
#' @param options A list of possible plot configurations. See comments below
#' @return sub A plotly widget
#' @export

plotEPG2 = function(mixData,kit,refData=NULL,AT=NULL,ST=NULL,dyeYmax=TRUE,plotRepsOnly=TRUE,options=NULL) {
 #AT (analyitcal threshold),ST (stochastic threshold). Can be given marker/dye specific
 require(plotly) #required package

 sn = names(mixData) #get samples names
 nS = length(sn) #number of replicates
 locs = names(mixData[[1]]) #get locus names

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

 #GRAPHICAL SETUP BASED ON SELECTED KIT:
 kitinfo = euroformix::getKit(kit) #names(kitinfo)
 if( is.na(kitinfo)[1] ) {
  print("The kit name was not recognized by getKit!")
  return()
 }
 dyes <- dyes2 <- unique(kitinfo$Color) #get dyes
 dyes2[dyes=="yellow"] = "orange" #exchange col because of illcondtioned
 dyes2[dyes=="green"] = "forestgreen" #exchange col because of illcondtioned

 nrows = length(dyes) #number of dyes/rows
 if(is.null(options$h0)) { h0 = 1200  } else { h0 = options$h0 }  # 5500/nrows #standard height for each dye (depends on number of rows? No)
 if(is.null(options$w0)) { w0 = 1800 } else { w0 = options$w0 }  # standard witdh when printing plot
 if(is.null(options$marg0)) { marg0 = 0.02  } else { marg0 = options$marg0 } 
 if(is.null(options$txtsize0)) { txtsize0 = 15 } else { txtsize0 = options$txtsize0 }  
 if(is.null(options$locsize0)) { locsize0 = 20 } else { locsize0 = options$locsize0 } 
 if(is.null(options$minY)) { minY = 100 } else { minY = options$minY }  #default minimum Y-axis length
 if(is.null(options$ymaxscale)) { ymaxscale = 1.05 } else { ymaxscale = options$ymaxscale }  #default minimum Y-axis length

 #Create list with dye,marker,bp (for observed data)
 bprng = range(kitinfo$Size) #get range (same range for all plots)
# bprng[1] = bprng[1]/2 #widen out on left?
 POS = aggregate(kitinfo$Size,by=list(kitinfo$Color,kitinfo$Marker),FUN=mean) #get marker positions (bp)

#Create dataset (per dye info with bp)
df = numeric() #store data: (sample,marker,allele,height,bp)
for(dye in dyes) {
#   dye=dyes[1]
  loctab = POS[POS[,1]==dye,-1,drop=FALSE]
  locs = toupper(as.character(loctab[,1])) #get locs (upper case variant)
  for(ss in sn) { #create a seperate EPG plot for each samples
#ss=sn[1]
   for(loc in locs) {
   #loc=locs[2]
    edat = mixData[[ss]][[loc]] #get evid data   
    if(is.null(refData)) {
      rdat = NULL
    } else {
      if(locFirst) rdat = refData[[loc]] #get ref data (list) #get ref data (list)      
      if(!locFirst) rdat = lapply(refData,function(x) x[[loc]]) #get ref data (list)      
    }
    av = edat$adata
    if( is.null(edat) && is.null(rdat)  ) next #skip if no data (evid or ref)
    hv = edat$hdata
    av2 = unique(unlist(rdat))
    adda = av2[!av2%in%av] 
    adda <- adda[!is.na(adda)] #remove NAs

    #add missing:
    if(length(adda)>0) {
     av = c(av,adda)
     hv  =  c(hv, rep(0,length(adda)) )
    }
    if(length(av)==0) next #skip if still no data

    #ref text under each allele
    reftxt = rep("",length(av))
    if(nrefs>0) {
     for(rr in 1:nrefs) { #for each ref
      indadd = which(av%in%unlist(rdat[[rr]])) #index of alleles to add to text
      hasprevval = indadd[nchar(reftxt[indadd])>0] #indice to add backslash (sharing alleles)
      reftxt[ hasprevval ] = paste0(reftxt[ hasprevval ],"/")      
      reftxt[indadd] = paste0( reftxt[indadd], rr)
     }
    }

    tmp = kitinfo[toupper(kitinfo$Marker)==loc,]
    ind = match(av,tmp$Allele) #get index to extract
    bv = tmp$Size[ind]  #get sizes directly from lookup
    isna = which(is.na(ind)) #which alleles are missing?
    if(length(isna)>0) avuse = as.numeric(tmp$Allele) #alleles available (called only once)
    for(missind in isna) {#impute missing bp:
      newa = as.numeric(av[missind ])
      impuse = which.min(abs(newa - avuse )) #index of closest allele 
      newa2 =  avuse[impuse] #closest allele 
      diff = newa - newa2
      bpadd1 = floor(diff)*tmp$Repeat[impuse]  #integer add
      bpadd2 = (diff-floor(diff))*10  #decimal add
      bv[missind] = tmp$Size[impuse] + bpadd1 + bpadd2  #estimate bp to insert
    }
    df = rbind(df, cbind(ss,loc,av,hv,bv,reftxt) )
   } #end for each samples
  } #end for each loci
} #end for each dye
df = data.frame(Sample=df[,1],Marker=df[,2], Allele=df[,3],Height=as.numeric(df[,4]),bp=as.numeric(df[,5]),reftxt=df[,6],stringsAsFactors=FALSE)
ymax1 =  ymaxscale*max(minY,df$Height) #global max y

#SEPARATE PLOTS
if(nS==1 || !plotRepsOnly) { #plot separate plot only in this case

for(ss in sn) { #create a seperate EPG plot for each samples
# ss =sn[1]
 plist = list() #create plot object for each color

 for(dye in dyes) {
#   dye=dyes[1]
  AT1 <- AT #temporary on analytical threshold
  ST1 <- ST #temporary on stochastic threshold
  if(!is.null(AT) && is.list(AT) ) AT1 = AT[[dye]] #ignores dye if not found
  if(!is.null(ST) && is.list(ST) ) ST1 = ST[[dye]] #ignores dye if not found

  dyeind = which(dyes==dye)
  dye2 = dyes2[dyeind]
  loctab = POS[POS[,1]==dye,-1,drop=FALSE]
  locs = toupper(as.character(loctab[,1])) #get locs 
  poslocs = loctab[,2] #get corresponding positions

  dfs = df[df$Sample==ss & df$Marker%in%locs,] #extract subset 
  if(dyeYmax) ymax1 = 1.05*max(minY,AT1,ST1,dfs$Height)  #get max 

  p = plotly::plot_ly(colors=dye2,mode="lines",height=h0) #df,x = ~bp,y=~Height,type="scatter",mode="markers",colors=dye2,name=~Allele)
  if(!is.null(AT1)) p <- plotly::add_lines(p,x = bprng, y = rep(AT1,2),color=factor(1),line=list(dash = 'dot',width=2),showlegend = FALSE)
  if(!is.null(ST1)) p <- plotly::add_lines(p,x = bprng, y = rep(ST1,2),color=factor(1),line=list(dash = 'dash',width=2),showlegend = FALSE)
  for(j in 1:nrow(dfs))  p = plotly::add_trace(p,x =dfs$bp[j] + 1*c( -1/4,0,1/4),y =c(0,dfs$Height[j],0 ),name=as.character(dfs$Allele[j]),type = "scatter" , mode = "lines", fill = "tozeroy",fillcolor=dye2,showlegend = FALSE,color=factor(1))
  p = plotly::add_annotations(p, x=poslocs ,y=ymax1,text=locs,showarrow=FALSE,font = list(color = dye2,family = 'Gravitas One',size = locsize0))  #ADD LOCI NAMES
  p = plotly::add_annotations(p, x=dfs$bp,y=rep(0,nrow(dfs)),text=dfs$Allele,showarrow=FALSE,font = list(color = 1,family = 'sans serif',size = txtsize0),yshift=-10)  #ADD ALLELE NAMES
  if(nrefs>0) {
     p = plotly::add_annotations(p, x=dfs$bp,y=rep(0,nrow(dfs)),text=dfs$reftxt,showarrow=FALSE,font = list(color = 1,family = 'sans serif',size = txtsize0),yshift=-25)  #ADD ALLELE NAMES
     if(dyeind ==1) p = plotly::add_annotations(p, x=rep(bprng[2],2),y=c(ymax1-ymax1/10*(1:nrefs)),text= paste0("Label ",1:nrefs,": ",refn),showarrow=FALSE,font = list(colors = 1,family = 'sans serif',size = 15),xshift=0,xanchor = 'right')  #ADD ALLELE NAMES
  }
  p = plotly::layout(p,xaxis = list(range = bprng,showticklabels = FALSE,title = ""),yaxis=list(range=c(0,ymax1), showline = TRUE,title = "Heights (RFU)"))#,colorway =dye2) 
  plist[[dye]] = p
 }
 sub = plotly::subplot(plist, nrows = nrows, shareX = TRUE, shareY = FALSE,margin=marg0,titleY= TRUE)
 sub = plotly::layout(sub ,title=ss,barmode = 'group',xaxis = list(title = ""))%>%plotly::config(scrollZoom=TRUE, displaylogo=FALSE,modeBarButtonsToRemove=c("hoverClosestCartesian","hoverCompareCartesian","toggleSpikelines"),toImageButtonOptions=list(width=w0)) 
 print(sub) 

 if(nS==1) return(sub) #return function if no replicates
}

} #end if
repcols = c("black","red","blue","forestgreen","orange","purple") [1:nS]

#REPS IN SAME PLOT
 plist = list() #create plot object for each color
 for(dye in dyes) {
#   dye=dyes[1]
  AT1 <- AT #temporary on analytical threshold
  ST1 <- ST #temporary on stochastic threshold
  if(!is.null(AT) && is.list(AT) ) AT1 = AT[[dye]] #ignores dye if not found
  if(!is.null(ST) && is.list(ST) ) ST1 = ST[[dye]] #ignores dye if not found

  dyeind = which(dyes==dye)
  dye2 = dyes2[dyeind]
  loctab = POS[POS[,1]==dye,-1,drop=FALSE]
  locs = toupper(as.character(loctab[,1])) #get locs 
  poslocs = loctab[,2] #get corresponding positions
  dfs = df[df$Marker%in%locs,] #extract subset 
  dfs1 = unique( subset(dfs,select=c("Marker","Allele","bp","reftxt") ) ) #extract unique info (to label loci/alleles etc)

  if(dyeYmax) ymax1 = ymaxscale*max(minY,AT1,ST1,dfs$Height)  #get max 
  p = plotly::plot_ly(dfs,type = "bar",height=h0, colors=repcols,showlegend = FALSE)
  if(!is.null(AT1)) p = plotly::add_segments(p,x = bprng[1], xend = bprng[2], y = AT1, yend = AT1,color=I(dye2),line=list(dash = 'dot',width=2))#,inherit=FALSE)
  if(!is.null(ST1)) p <- plotly::add_lines(p,x = bprng, y = rep(ST1,2),color=factor(1),line=list(dash = 'dash',width=2),showlegend = FALSE)
  p = plotly::add_trace(p,x=~bp,y=~Height,name=~Sample,showlegend = FALSE,color=~Sample,hoverlabel=list(font=list(size=12),namelength=1000),text =~Allele) #dfs,x = ~bp,y=~Height,type="scatter",mode="markers",colors=dye2,name=~Allele)
  p = plotly::add_annotations(p, x=poslocs ,y=ymax1,text=locs,showarrow=FALSE,font = list(color = dye2,family = 'Gravitas One',size = locsize0))  #ADD LOCI NAMES
  p = plotly::add_annotations(p, x=dfs1$bp,y=rep(0,nrow(dfs1)),text=dfs1$Allele,showarrow=FALSE,font = list(color = 1,family = 'sans serif',size = txtsize0),yshift=-10)  #ADD ALLELE NAMES
  if(nrefs>0) {
     p = plotly::add_annotations(p, x=dfs1$bp,y=rep(0,nrow(dfs1)),text=dfs1$reftxt,showarrow=FALSE,font = list(color = 1,family = 'sans serif',size = txtsize0),yshift=-25)  #ADD ALLELE NAMES
     if(dyeind ==1) p = plotly::add_annotations(p, x=rep(bprng[2],2),y=c(ymax1-ymax1/10*(1:nrefs)),text= paste0("Label ",1:nrefs,": ",refn),showarrow=FALSE,font = list(colors = 1,family = 'sans serif',size = 15),xshift=0,xanchor = 'right')  #ADD ALLELE NAMES
  }
  p = plotly::layout(p,xaxis = list(range = bprng,showticklabels = FALSE,title = ""),yaxis=list(range=c(0,ymax1), showline = TRUE,title = "Heights (RFU)"))#,autosize=FALSE,width=10)#,colorway =dye2) 
  plist[[dye]] = p
 }
 sub = plotly::subplot(plist, nrows = nrows, shareX = FALSE, shareY = FALSE,margin=marg0,titleY= TRUE) 
 sub = plotly::layout(sub ,title=paste0(sn,collapse="/"),barmode = 'group')%>%plotly::config(scrollZoom=TRUE, displaylogo=FALSE,modeBarButtonsToRemove=c("lasso2d","select2d","hoverClosestCartesian","hoverCompareCartesian","toggleSpikelines"),toImageButtonOptions=list(width=w0)) 
 print(sub)

 return(sub)
} #end function





  