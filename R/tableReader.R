#' @title tableReader
#' @author Oyvind Bleka
#' @description  A robust function for reading tables
#' @details This function reads text based tables and recognizes which of the three common separates that are used.
#' @param filename Name of file to be read.
#' @param header Whether header should be included or not
#' @param delims delimiter options to use when reading table
#' @export

 tableReader=function(filename, header=TRUE, delims = c("\t",",",";")) {
   
  reader = function(delim) read.table(filename,header=header,sep=delim,as.is=TRUE,row.names=NULL,check.names = FALSE,comment.char = "") 
  #NOTES: check.names=FALSE important for avoiding . for spaces.  comment.char is remove to handle '#'
  for(delim in delims) {
    tab = reader(delim)
    if(ncol(tab)>1) break #stop if proper table
  }
  return(tab) 
 }