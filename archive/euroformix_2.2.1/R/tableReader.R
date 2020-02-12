#' @title tableReader
#' @author Oyvind Bleka
#' @description  A robust function for reading tables
#' @details This function reads text based tables and recognizes which of the three common separates that are used.
#' @export
 tableReader=function(filename,header=TRUE) {
  readF <- function(...) { #read function
   #require(data.table) #used to read very large files fast
   #fread(...)
   read.table(...)
  }
  tab <- readF(filename,header=header,sep="\t",stringsAsFactors=FALSE, check.names=FALSE)
  tryCatch( {  if(ncol(tab)==1) tab <- readF(filename,header=header,sep=",",stringsAsFactors=FALSE, check.names=FALSE) } ,error=function(e) e) 
  tryCatch( {  if(ncol(tab)==1) tab <- readF(filename,header=header,sep=";",stringsAsFactors=FALSE, check.names=FALSE) } ,error=function(e) e) 
  return(tab) 
 }