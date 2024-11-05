#' @title sample_tableToList
#' @author Oyvind Bleka
#' @description Converting table to list format
#' @param table A table with data (Evid or refs): A list on form [[name]][[locus]]
#' @return outL A data list
#' @export

sample_tableToList = function(table) {
  colname = colnames(table) #colnames in file
  lind = grep("marker",tolower(colname),fixed=TRUE) #locus col-ind
  if(length(lind)==0) lind = grep("loc",tolower(colname),fixed=TRUE) #try another name
  sind = grep("sample",tolower(colname),fixed=TRUE) #sample col-ind
  if(length(sind)>1)  sind = sind[grep("name",tolower(colname[sind]),fixed=TRUE)] #use only sample name
  A_ind = grep("allele",tolower(colname),fixed=TRUE) #allele col-ind
  H_ind = grep("height",tolower(colname),fixed=TRUE) #height col-ind
  locs = unique(toupper(table[,lind])) #locus names: Use uniques and Convert to upper case
  samplenames = unique(as.character(table[,sind])) #sample names
  outL = list() #Init outList (insert non-empty characters):
  for(samplename in samplenames) { #for each sample in matrix
    outL[[samplename]] = list() #one list for each sample
    for(loc in locs) { #for each locus
      rowinds = which(table[,sind]==samplename & toupper(table[,lind])==loc) #get row index in table for given sample and locus
      if(sum(rowinds)==0) next #no data found (skip marker)
      
      alleles <- heights <- numeric() #init vector with data
      for(rowind in rowinds) { #for each row indices
		    av <- table[rowind,A_ind] #obtain raw alleles
        keep <- !is.na(av) & !(av%in%c(""," ")) #get boolean vector of alleles to keep
        alleles <- c(alleles, table[rowind,A_ind[keep]] ) #extract allele(s), be be levels
        if(length(H_ind)>0) {
          heights <- c(heights, table[rowind,H_ind[keep]]) #extract peak heights (intensities)
        }
      } #end for each rows
      
      insAlleles <- character() #default if no data
      insHeights <- numeric() #default  if no data
      if(length(alleles)>0) { #if any data to store:
        if(length(A_ind)>0) insAlleles = as.character(alleles)
        if(length(H_ind)>0) insHeights = as.numeric(as.character(heights))
      }
      if(length(A_ind)>0) outL[[samplename]][[loc]]$adata = insAlleles
      if(length(H_ind)>0) outL[[samplename]][[loc]]$hdata = insHeights
    } #end for each loci
  } #end for each samples
  return(outL)
}