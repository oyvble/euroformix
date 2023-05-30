#' @title getSampleType
#' @description A function for determining data format
#' @details returns "EPG" (CE),"LUS" (MPS with LUS format) or "MPS" (else it is string)
#' @param samples A list of samples (evidence) with structure [[samplename]][[locus]] = list(adata,...)
#' @param kit Kitname returned from getKit()
#' @param LUSsymbol Symbol to be recognized as a LUS allele

#helpfunction to determine type (
getSampleType = function(samples,kit=NULL,LUSsymbol="_") { 
  isMPS = TRUE #can always be plotted
  isLUS = all(unlist(sapply(samples,function(x)  sapply(x,function(y) all(grepl(LUSsymbol,y$adata),na.rm=TRUE)) )))  #ADDED: check if alleles are given on LUS format 
  isEPG = FALSE 
  if(!is.null(kit) && !is.na(kit) && kit!="" ) { #kit must be proper name
    kitinfo = getKit(kit) #check kit
    if(length(kitinfo)>1) isEPG = TRUE
    #print("Please specify kit in order to show EPG. Otherwise default MPS format is shown!")
  }  
  
  type = "MPS" #default
  if(isLUS) {
    type = "LUS"
  } else if(isEPG) {
    type = "EPG"
  } 
  if(type=="EPG") {  #Special case if kit is given
    colorinfo = euroformix::getKit(kit,"COLOR") #Extract color info of kit
    if( length( unique(colorinfo[,2]) )==1 ) type = "MPS" #SET AS MPS IF ONLY 1 COLOR
  }
  return(type)
}
