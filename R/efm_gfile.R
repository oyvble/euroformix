#' @title efm_gfile
#' @author Oyvind Bleka
#' @description A helpfunction to obtain file names from browser
#' @details Ensures that correct coding is applied to the file name returned from gWidgets2::gfile
#' @param text Input argument to gfile
#' @param type Input argument to gfile
#' @param filter Input argument to gfile
#' @param initf Input argument to gfile
#' @export

#Helpfunction to interact with file I/O interface
#This function is written since the encoding in  gWidgets2::gfile is fixed to UTF-8 which doesn't handle special letters
efm_gfile <- function(text,type,filter=list(),initf=NULL) { #Possible bug: text ignored when type="selectdir"
  file <- gWidgets2::gfile(text=text,type=type,filter=filter,initial.filename=initf)
  Encoding(file) <- options()$encoding #Set to local encoder: Handle special cases.
  return(file)
}

