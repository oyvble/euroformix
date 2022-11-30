#' @title tableSaver
#' @author Oyvind Bleka
#' @description  A robust function for saving tables (interacts with gfile)
#' @details This function prompts the user to obtain a filename and saves text based tables and ensures that a proper extension is used (txt or csv)
#' @param tab Table to be saved
#' @param sep Separator to use for saving table (decided by the program)
#' @export


tableSaver = function(tab,sep="txt") {
  tabfile  = efm_gfile(text="Save table",type="save") #csv is correct format!
  if(length(tabfile)==0) return()
  
  #Ensure that a proper extension is given to the text file
  ext <- strsplit(basename(tabfile), split="\\.")[[1]] #get extension
  addExt = FALSE #whether to add extension to full file name
  if(length(ext)==1) { #IF extension is needed (YES if length was 1)
    addExt = TRUE
  } else if(length(ext)>1) { #obtain extension candidates: Note that . may be in filename
    if(!ext[length(ext)]%in% c("txt","csv")) addExt = TRUE #check if last extension is valid. Add extension IF NOT
  }
  if(addExt) tabfile = paste0(tabfile,".",sep) #then simply add suggested extension
  
  if(sep=="txt" | sep=="tab") write.table(tab,file=tabfile,quote=FALSE,sep="\t",row.names=FALSE) 
  if(sep=="csv") write.table(tab,file=tabfile,quote=FALSE,sep=";",row.names=FALSE) 
  #if(sep=="csv2") write.table(tab,file=tabfile,quote=FALSE,sep=";",row.names=FALSE) 
  print(paste("Table saved in file: ",tabfile,sep=""))
} #end file

