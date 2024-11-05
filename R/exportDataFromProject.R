#' @title exportDataFromProject
#' @author Oyvind Bleka
#' @description Exports data from a EFM-project file into multiple files in EFM format
#' @param projectFile File name (full path) of the project file
#' @param evidfn File name for the exported evidence data
#' @param reffn File name for the exported reference data
#' @param freqfn File name for the exported frequency data
#' @export 

#library(euroformix)
#exportDataFromProject(envirfile)
exportDataFromProject = function(projectFile,evidfn = "evids.csv", reffn = "refs.csv",  freqfn = "freqs.csv") {
  mmTK = NULL #declared
  load(projectFile) #loading environment
  if(is.null(mmTK)) stop("Couldnt read project")
  
  #export evidence data:
  evidList = mmTK$mixData
  if(is.null(evidList)) {
    print("Couldnt find evidence data..")
  } else {
    evidTable = sample_listToTable(evidList)
    write.table(evidTable,file=evidfn,quote=FALSE,sep=";",row.names = FALSE)
  }
  
  #export refernce data:
  refList = mmTK$refData
  if(is.null(evidList)) {
    print("Couldnt find reference data..")
  } else {
    refTable = sample_listToTable(refList)
    write.table(refTable,file=reffn,quote=FALSE,sep=";",row.names = FALSE)
  }
  
  #export frequency data:
  popFreq = mmTK$popFreq
  if(is.null(evidList)) {
    print("Couldnt find frequency data..")
  } else {
    freqTable = freqs_listToTable(popFreq)
    write.table(freqTable,file=freqfn,quote=FALSE,sep=";",row.names = FALSE)
  }
}
  