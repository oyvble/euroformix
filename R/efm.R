#' @title efm
#' @author Oyvind Bleka
#' @description efm (EuroForMix) is a GUI wrapper for euroformix
#' @details The function is a graphical layer for the functions in the package euroformix. See ?euroformix for more information.
#' @param envirfile A Rdata file including a saved environment of a project
#' @export

#library(euroformix);envirfile=NULL#;efm()
efm = function(envirfile=NULL) {
 LUSsymbol = "_" #Added in version 1.11.0 defined as constant. Used for showing MPS based STRs in RU_LUS format (LUS must be numeric)
 MPSsymbol = ":" #Added in version 2.2.0 defined as constant. Used for showing MPS based SNPs/STRs in RU:something format (both can be strings)

 #size of main window
 mwH <- 500
 mwW <- 1000

 #type of gwidgets-kit
 #library(gWidgetstcltk)
 options(guiToolkit="tcltk")

 #version:
 version = utils::packageVersion("euroformix") #follows same version as package number

 #software name:
 softname <- paste0("EuroForMix v",version)

 #GUI-Restriction on maximum number of unknowns
 maxUsetup <- 4 

 #Spacing between widgets
 spc <- 10
 
 #Strider link (can be customized in toolbar)
 striderlink = "https://strider.online/frequencies_xml/download" #"https://strider.online/frequencies/download"

 #Marker specific setting names (objects)
 objnameMarkerSettings = c("threshv","pCv","lamv","fstv") #get object names for marker settings (optMarkerSetup)
 
 #####################
 #create environment #
 #####################
  pgkPath <- path.package("euroformix", quiet = FALSE) # Get package path.
  .sep <- .Platform$file.sep # Platform dependent path separator. 
  deffreq <- paste(pgkPath,"FreqDatabases",sep=.sep) #default path to freq-files

  #the file used to store system settings (opt-settings)
  optSetupFile <- paste(pgkPath,"configSetup",sep=.sep)  
  if(!file.exists(optSetupFile)) {  #use default values if not existing
    optSetup = list(easyMode=FALSE,adjQbp=FALSE,thresh0=50,fst0=0,pC0=0.05,lam0=0.01,pXi="dbeta(x,1,1)",pXiFW="dbeta(x,1,1)")
    #adjQbp: adjusted Q-allele fragment length (bp)
    #pXi: prior density function of stutter proportion parameter
    #thresh0: default value of detection threshold value
  } else {
    optF <- scan(file=optSetupFile,what=character(),quiet=TRUE)
    optSetup = list(easyMode=optF[1]=="TRUE",adjQbp=optF[2]=="TRUE",thresh0=as.numeric(optF[3]),fst0=as.numeric(optF[4]),pC0=as.numeric(optF[5]),lam0=as.numeric(optF[6]),pXi=optF[7],pXiFW=optF[8])
  }

  #the file used to store per-marker settings (opt-settings)
  optMarkerSetupFile <- paste(pgkPath,"configMarkerSetup",sep=.sep)  
  if(!file.exists(optMarkerSetupFile)) {  #use default values if not existing
    optMarkerSetup = NULL #object is empty if not specified 
  } else {
    optF <- readLines(optMarkerSetupFile) #read saved data
    nLocs = as.integer(optF[1]) #first element is number of markers
    mat = matrix(optF[-1] ,nrow=nLocs) #obtain matrix with locus information
    
    optL = list() #storing objects from matrix
    for(c in 2:ncol(mat)) {
      optL[[c-1]] = as.numeric(mat[,c] ) #convert to numbers
      names(optL[[c-1]]) = mat[,1]  #insert locus names
    }
    names(optL) = objnameMarkerSettings
    optMarkerSetup = optL
  }
  
  
 if(is.null(envirfile)) {
  mmTK = new.env( parent = globalenv() ) #create new envornment object (must be empty)

  #Toolbar options: can be changed any time by using toolbar
  assign("optSetup",optSetup,envir=mmTK) 
  assign("optMarkerSetup",optMarkerSetup,envir=mmTK) #default is common marker settings
  
  assign("optFreq",list(freqsize=0,wildsize=5,minF=NULL,normalize=1),envir=mmTK) #option when new frequencies are found (size of imported database,minFreq), and missmatch options
  assign("optMLE",list(nDone=3,delta=1,dec=4,difftol=0.01,maxThreads=0,seed=NULL,steptol=1e-3,alpha=0.01),envir=mmTK) #options when optimizing,validation (nDone,delta)
  assign("optMCMC",list(delta=2,niter=2000,quantile=0.05,seed=1),envir=mmTK) #options when running MCMC-simulations (delta, niter,seed=1)
  assign("optINT",list(reltol=0.1,maxeval=20000,dev=3),envir=mmTK) #options when integrating (reltol and settings for boundaries)
  assign("optDC",list(alphaprob=0.99,maxlist=20),envir=mmTK) #options when doing deconvolution (alphaprob, maxlist)
  assign("optDB",list(maxDB=10000,QUALpC=0.05,ntippets=10),envir=mmTK)  #options when doing database search (maxDB)
  assign("optLRMIX",list(range=0.6,nticks=31,nsample=2000,alpha=0.05),envir=mmTK) #options when doing LRmix
  assign("STRidER",list(url=striderlink),envir=mmTK) #path for importing STRidER data
  
  #initializing environment variables
  assign("workdir",NULL,envir=mmTK) #assign working directory to mmTK-environment
  assign("popList",NULL,envir=mmTK) #Added in v1.10.0, assign imported popfreq list to mmTK-environment
  assign("selPopKitName",NULL,envir=mmTK) #selected kit and population for popFreq

  #imported data:
  assign("popFreq",NULL,envir=mmTK) #assign to mmTK-environment
  assign("mixData",NULL,envir=mmTK) #assign to mmTK-environment
  assign("refData",NULL,envir=mmTK) #assign to mmTK-environment
  assign("dbData",NULL,envir=mmTK) #assign dbData: matrix referenceDatabase to mmTK-environment (encoded)

  #models: stored setups for model specification
  assign("setEVID",NULL,envir=mmTK) #assign model (evidence weighting)
  assign("setEVIDbag",NULL,envir=mmTK) #create bag of calculated WOE
  assign("setDB",NULL,envir=mmTK) #assign model (database search)
  assign("setDC",NULL,envir=mmTK) #assign model (deconvolution)
  assign("setGEN",NULL,envir=mmTK) #assign model (generation)

  #results: stored results after calculations
  assign("resDB",NULL,envir=mmTK) #assign database search results (i.e. ranked table of results)
  assign("resEVID",NULL,envir=mmTK) #assign evidence weighting results (i.e. calculated LR with MLE estimates)
  assign("resDC",NULL,envir=mmTK) #assign deconvolved results (i.e. ranked tables of results)
  assign("resEVIDLRMIX",NULL,envir=mmTK) #assign evidence weighting results - Based on LRmix
 } else { #restore from file
    load(envirfile) #loading environment
     
    #MAKING EFM-PROJECTS BACKWARD COMPATIBLE
     
    #v0.6.0 
    if( is.null( mmTK$optSetup )) assign("optSetup",optSetup,envir=mmTK)  
     
    #1.9.4
    if( is.null( mmTK$optMLE$difftol ))  mmTK$optMLE$difftol  <- 0.01 #set default value
    
    #v1.10.0
    if( is.null( mmTK$popList)) assign("popList",NULL,envir=mmTK) #Added in v1.10.0, 
     
    #v3.0.0
    if( is.null( mmTK$optFreq$normalize )) { #whether allele frequencies should be normalized after including new alleles
      mmTK$optFreq$normalize <- 1 #set default value (YES)
      assign("optFreq",mmTK$optFreq,envir=mmTK)   #store again
    }
    if( is.null( mmTK$optSetup$pXiFW )) mmTK$optSetup$pXiFW="dbeta(x,1,1)" #set default prior distribution of beta
    if( is.null( mmTK$optMarkerSetup )) assign("optMarkerSetup",NULL,envir=mmTK) #default is common marker settings
    if( is.null( mmTK$resEVID$LRupper )) mmTK$resEVID$LRupper = NA
    if( is.null( mmTK$resEVID$adjLRmle )) mmTK$resEVID$adjLRmle = NA

    #v3.2.0
    if( is.null( mmTK$STRidER )) {
      mmTK$STRidER$url  <- striderlink #set default value
      assign("STRidER",mmTK$STRidER,envir=mmTK)   #store again
    }
    
    #v4.0.0
    if( !is.null(mmTK$optSetup$maxloc)) { #v4.0
      mmTK$optSetup$maxloc = FALSE #adjQbp
      names(mmTK$optSetup)[names(mmTK$optSetup)%in%"maxloc"] = "adjQbp" #rename variable
    }
    if( is.null( mmTK$optINT$dev )) {
      mmTK$optINT$dev <- 3 #set default value
      assign("optINT",mmTK$optINT,envir=mmTK)   #store again
    }
    assign("optSetup",mmTK$optSetup,envir=mmTK)   #store agaimn
    
    
    #from v3.0.0 onwards
    if( mmTK$optMLE$delta==10 ) mmTK$optMLE$delta=1 #be sure that randomizer is reduced (if loaded old projects)
    if( is.null( mmTK$optMLE$maxThreads )) mmTK$optMLE$maxThreads  <- 0 #set default value
    if( is.null( mmTK$optMLE$steptol ))  mmTK$optMLE$steptol  <- 1e-3 #set default value
    if( is.null( mmTK$optMLE$alpha ))  mmTK$optMLE$alpha  <- 0.01 #set default value
    if( is.null( mmTK$optMCMC$quantile )) mmTK$optMCMC$quantile=0.05 #set default quantile for mcmc
    if( is.null( mmTK$optMCMC$seed )) mmTK$optMCMC$seed=1 #set default seed for mcmc
    assign("optMLE",mmTK$optMLE,envir=mmTK)   #store again
    assign("optMCMC",mmTK$optMCMC,envir=mmTK)   #store again
    assign("resEVID",mmTK$resEVID,envir=mmTK) #assign evidence weighting results (i.e. calculated LR with MLE estimates)
    
  }

 ####################################
 #auxiliary functions and variables:#
 ####################################
  
 NAtoSign <- function(x) {
  x[is.na(x)] <- "-" #NB: New version of gtable does not accept NA values
  return(x)
 }
# helptext = function(obj,txt) { gWidgets2::addHandlerRightclick(obj,handler=function(h,...) { gWidgets2::gmessage(txt,title="Detailed information") }) }
  helptext = function(obj,txt) { gWidgets2::tooltip(obj) = txt } 
 
 #helpfunction to return minimum frequency (used for new alleles)
 getminFreq = function() {
  popFreq <- get("popFreq",envir=mmTK) #get selected popFreq
  freqsize <- get("optFreq",envir=mmTK)$freqsize #get selected size of frequence-database
  minfreq <- get("optFreq",envir=mmTK)$minF #get selected size of frequence-database
  if(freqsize>0) {
   return(5/(2*freqsize))
  } else if(is.null(minfreq)) {
   return(min(unlist(popFreq))) #minumum observed frequency was used 
  } else {
   return(minfreq) #specified observed frequency was used 
  }
 }

 #Function to get data from environment
 #sel used to select a specific datasubset
 getData = function(type,sel=NULL) {
   Data <- NULL
   if(type=="mix") Data <- get("mixData",envir=mmTK) #assign kit to mmTK-environment
   if(type=="ref") Data <- get("refData",envir=mmTK) #assign kit to mmTK-environment 
   if(type=="db") Data <- get("dbData",envir=mmTK) #assign kit to mmTK-environment 
   if(!is.null(sel)) return(Data[sel]) #returns only selected datasubset
   return(Data)
 }

 #function for inserting sample/ref/db-names into existing gWidgets2::gcheckboxgroup
 getDataNames_type = function(type) {
  subD <- getData(type)
  if(!is.null(subD)) { return( names(subD))
  } else { return("") }
 }

 #Function which takes rownames and adds to first column
 addRownameTable = function(tab) {
  tmp <- colnames(tab)
  tab <- cbind(rownames(tab),tab)
  colnames(tab) <- c("X",tmp)
  return(tab)
 }



###################################################################
#############GUI HELPFUNCTIONS#####################################
###################################################################

 
 #Helpfunction to get focus 
 getFocus = function() {
   gWidgets2::visible(mainwin) <- TRUE
   gWidgets2::focus(mainwin) <- TRUE
 }
 
 #Menu bar file-lists:
 f_setwd = function(h,...) {
  dirfile = efm_gfile(text="Select folder",type="selectdir")
  if(length(dirfile)==0) return()
  setwd(dirfile)
  assign("workdir",dirfile,envir=mmTK) #assign working directory
 }
  
 f_openproj = function(h,...) {
  projfile = efm_gfile(text="Open project",type="open", filter=list("Project"=list(patterns=list("*.Rdata")) , "All files"=list(patterns=list("*"))))
  if(length(projfile)==0) return()
  gWidgets2::dispose(mainwin)
  efm(projfile) #send environment into program
 }
 
 f_saveproj = function(h,...) {
  projfile = efm_gfile(text="Save project",type="save")
  if(length(projfile)==0) return()
  
   if(length(unlist(strsplit(projfile,"\\.")))==1) projfile = paste0(projfile,".Rdata") #add extension if missing
   tmp = sort(sapply(mmTK,object.size)/1e6,decreasing=TRUE)
   print(paste0("Size of stored objects (in MB): ",round(sum(tmp),2))) #prints size of each stored object
   print(tmp) #prints size of each stored object
   #save(mmTK,file=projfile,compress="xz") #save environment, #Dont compress!
   save(mmTK,file=projfile,compress="xz",eval.promises=FALSE,precheck=FALSE,compression_level=2)
   print(paste("Project saved in ",projfile,sep=""))
  
 }

 f_settings = function(h,...) { #creates a table, with different settings
   opt <- get("optSetup",envir=mmTK) 
   setwin <- gWidgets2::gwindow(paste0("Settings"),visible=FALSE)
   tabval = gWidgets2::glayout(spacing=0,container=(setwin)) 
   w0 <- 15 #width of text box
   #Model parameters: 
   tabval[1,1] <- gWidgets2::glabel(text="Easy mode:",container=tabval)
   tabval[1,2] <- gWidgets2::gradio(items=c("NO","YES"), selected = as.numeric(opt$easyMode==TRUE)+1, horizontal = TRUE,container=tabval)
   tabval[2,1] <- gWidgets2::glabel(text="Analytical threshold (AT)",container=tabval)
   tabval[2,2] <- gWidgets2::gedit(text=opt$thresh0,width=w0,container=tabval)
   tabval[3,1] <- gWidgets2::glabel(text="Fst-correction (theta)",container=tabval)
   tabval[3,2] <- gWidgets2::gedit(text=opt$fst0,width=w0,container=tabval)
   tabval[4,1] <- gWidgets2::glabel(text="Probability of drop-in (PrC)",container=tabval)
   tabval[4,2] <- gWidgets2::gedit(text=opt$pC0,width=w0,container=tabval)
   tabval[5,1] <- gWidgets2::glabel(text="Drop-in hyperparam (lambda)",container=tabval)
   tabval[5,2] <- gWidgets2::gedit(text=opt$lam0,width=w0,container=tabval)
   tabval[6,1] <- gWidgets2::glabel(text="Prior: BW stutter-prop. \n function(x)=",container=tabval)
   tabval[6,2] <- gWidgets2::gedit(text=opt$pXi,width=w0,container=tabval)
   tabval[7,1] <- gWidgets2::glabel(text="Prior: FW Stutter-prop. \n function(x)=",container=tabval)
   tabval[7,2] <- gWidgets2::gedit(text=opt$pXiFW,width=w0,container=tabval)
   tabval[8,1] <- gWidgets2::glabel(text="Adjust fragmentlength \nof Q-allele: ",container=tabval)
   tabval[8,2] <- gWidgets2::gradio(items=c("NO","YES"), selected = as.numeric(opt$adjQbp==TRUE)+1, horizontal = TRUE,container=tabval)
   tabval[9,1] <- gWidgets2::gbutton("Save", container=tabval,handler = function(h, ...) { 

    #Delete previous marker settings if clicking save 
    if(file.exists(optMarkerSetupFile)) file.remove(optMarkerSetupFile) #delete file if exists
    assign("optMarkerSetup",NULL,envir=mmTK) #default is common marker settings
    opt$easyMode <- gWidgets2::svalue(tabval[1,2])=="YES" #easy Mode?
    opt$thresh0 <- as.numeric(gWidgets2::svalue(tabval[2,2]))
    opt$fst0 <- as.numeric(gWidgets2::svalue(tabval[3,2])) 
    opt$pC0 <- as.numeric(gWidgets2::svalue(tabval[4,2])) 
    opt$lam0 <- as.numeric(gWidgets2::svalue(tabval[5,2]))
    opt$pXi <- gWidgets2::svalue(tabval[6,2]) #get pXi function
    opt$pXiFW <- gWidgets2::svalue(tabval[7,2]) #get pXiFW function
    opt$adjQbp <- gWidgets2::svalue(tabval[8,2])=="YES" #adjusted fragmenthlength of Q-allele?
    if( any( sapply(opt,is.na) ) ) { 
       gWidgets2::gmessage("Invalid input in settings, please set another value!",title="Wrong input",icon="error")
       return()
    }
   assign("optSetup",opt,envir=mmTK)  #assign user-value to opt-list
   opt2 =   c(opt$easyMode,opt$adjQbp,opt$thresh0,opt$fst0,opt$pC0,opt$lam0,opt$pXi,opt$pXiFW) #ensure correct order
   write(opt2,file=optSetupFile)    #save to file in installation folder (must be correct order)
   gWidgets2::dispose(setwin)
  } )
  gWidgets2::visible(setwin) <- TRUE
 }

 #creates a table with different settings for each markers (based on popfreq)
 f_markerSettings = function(h,...) { 
   optL = get("optMarkerSetup",envir=mmTK)  #may contain stored marker based info
   opt0 <- get("optSetup",envir=mmTK) #contains the default settings (common marker settings)
   popFreq = get("popFreq",envir=mmTK)
   if(is.null(popFreq)) {
     gWidgets2::gmessage("Please select population frequencies to get marker specific settings!\n\nSelecting kit will also give you dye (color) information. ",title="Missing data",icon="error")
     return() 
   }
   locs <- toupper(names(popFreq)) #obtain marker names (same as in popFreq). Force to capital letters
     
   kitname <- get("selPopKitName",envir=mmTK)[1]  #get name of selected kit
   dyes = NULL #default is no colors (dyes)
   if(!is.null(kitname)) {
     kitdyes = getKit(kitname,"COLOR") #get kitinfo from selected kit
     if(!is.na(kitdyes) && length(kitdyes)>1) { #if kitinformatation found
       locDyeTab = kitdyes[toupper(kitdyes$Marker)%in%locs,,drop=FALSE] #obtain list (overlap in loci)
       locs = toupper(locDyeTab[,1]) #update loci order (force to be capital letter)
       dyes = locDyeTab[,2] #get corrsponding colors
     }
   } 
   nLocs = length(locs) #get number of loci
   if( nLocs > 50  ) {
     gWidgets2::gmessage("The number of markers was large. This may take some time to obtain!",title="Large data",icon="error")
   }
   nObj = length(objnameMarkerSettings) #number of objects
   
   getvec = function(x) { #helpfunction to assign locus name on vector
     names(x) = locs
     return(x)
   }

   setToDefault = function(init=FALSE,usecommon=FALSE) { #helpfunctions to set values back to default (init is variable declaration)
     opt0default = c(opt0$thresh0,opt0$pC0,opt0$lam0,opt0$fst) #default values (if element not found in optMarkerSetup )
     for(m in 1:nLocs) {
       if(init) tabval[m+1,1] <- gWidgets2::glabel(text=locs[m],container=tabval) #init markers
       for(c in 1:nObj) { #for each objects
         insval = opt0default[c] #default value to insert
         if(!usecommon && !is.null(optL)) { #if not using common values (from settings) restoring data 
           tmp= optL[[c]]
           tmp = tmp[names(tmp)==locs[m]] #extract value of relevant marker
           if(length(tmp)==1) insval = tmp
         } 
         if(init) {
           tabval[m+1,c+1] <- gWidgets2::gedit(text=insval,width=w0,container=tabval) #insert value
         } else {
           gWidgets2::svalue(tabval[m+1,c+1]) <- insval #insert default value
         }
       } #end for each objects       
       if(!is.null(dyes) && init)  tabval[m+1,6] <- gWidgets2::glabel(text=dyes[m],container=tabval) #init dyes
     } #end for each loci
   }
   
   #CREATE GUI WINDOW:
   setwin <- gWidgets2::gwindow(paste0("Marker specific settings"),visible=FALSE, width=750, height=500)
   grouplay <- gWidgets2::ggroup(spacing=2,container=(setwin),horizontal = FALSE, use.scrollwindow = TRUE)  #set group layout
   
   tabsel = gWidgets2::glayout(spacing=10,container=(grouplay),horizontal = TRUE)  #set group (will contain buttons)
   tabsel[1,1] = gWidgets2::gbutton(text="Restore",container=tabsel, handler = function(h,...) {
     gWidgets2::visible(setwin) = FALSE #hide window
     setToDefault(FALSE) #dont init 
     gWidgets2::visible(setwin) = TRUE #show window
     gWidgets2::focus(setwin) #set top
     
   }) 
   tabsel[1,2] = gWidgets2::gbutton(text="Set to default",container=tabsel, handler = function(h,...) {
     gWidgets2::visible(setwin) = FALSE #hide window
     setToDefault(FALSE,TRUE) #dont init 
     gWidgets2::visible(setwin) = TRUE #show window
     gWidgets2::focus(setwin) #set top
     
   }) 
   tabsel[1,3] = gWidgets2::gbutton(text="Fill out dye info",container=tabsel, handler = function(h,...) {
     gWidgets2::visible(setwin) = FALSE #hide window) 
     #filling in values over all dyes (using those specified in  top dye )
     unDyes = unique(dyes) #obtain unique colors
     for(c in 1:nObj) { #for each objects
       for(dye in unDyes) {  #for each dye
         rowind = which(dyes==dye) #obtain row indices of selected dyes
         vals = sapply(rowind,function(x) gWidgets2::svalue(tabval[x+1,c+1])) #obtain values
         useval = vals[vals!=""] #obtain value in cells
         if(length(useval)>0) {
           for(ii in rowind) gWidgets2::svalue(tabval[ii+1,c+1]) <- useval[1] #fill with first value only
         }
       }
     }
     gWidgets2::visible(setwin) = TRUE #show window
     gWidgets2::focus(setwin) #set top
   })
   if(is.null(dyes)) gWidgets2::enabled(tabsel[1,3]) = FALSE #deactivate if no dye information found

   tabsel[1,4] = gWidgets2::gbutton(text="Empty all",container=tabsel, handler = function(h,...) {
     gWidgets2::visible(setwin) = FALSE #hide window
     for(m in 1:nLocs) {
       for(c in 1:nObj) {
         gWidgets2::svalue(tabval[m+1,c+1]) = "" #insert emtpy cell
       }
     }
     gWidgets2::visible(setwin) = TRUE #show window
     gWidgets2::focus(setwin) #set top
   } ) 
   
   #extracting values from table
   tabsel[1,5] = gWidgets2::gbutton(text="Save settings",container=tabsel, handler = function(h,...) {
     optL2 = list()
     for(c in 1:nObj) {
       vecExtract = rep(NA,nLocs) #vector to extract
       for(m in 1:nLocs) {
         tmp = as.numeric( gWidgets2::svalue(tabval[m+1,c+1]) )  #convert from text to number
         checkPositive(tmp,"The cell value") #check if pos value
         vecExtract[m] = tmp
       }
       optL2[[c]] = getvec(vecExtract)
     }  
     names(optL2) = objnameMarkerSettings #naming objects
     assign("optMarkerSetup",optL2,envir=mmTK)  #store to internal object
     
     #Store to file:     
     write(c(nLocs,locs,unlist(optL2)),file=optMarkerSetupFile)    #save to file in installation folder
     gWidgets2::dispose(setwin) #close window when successful
   }) 
   
   tabval = gWidgets2::glayout(spacing=0,container=(grouplay))  #create grid layout
   w0 <- 15 #width of textbox

   tabval[1,1] <- gWidgets2::glabel(text="Marker",container=tabval) 
   tabval[1,2] <- gWidgets2::glabel(text="Analyt. thresh (AT)",container=tabval) 
   tabval[1,3] <- gWidgets2::glabel(text="Dropin prob. (pC)",container=tabval) 
   tabval[1,4] <- gWidgets2::glabel(text="Hyperparam (lambda)",container=tabval) 
   tabval[1,5] <- gWidgets2::glabel(text="Fst-correction (theta)",container=tabval) 
   if(!is.null(dyes)) tabval[1,6] <- gWidgets2::glabel(text="Dye (color)",container=tabval) 
   
   setToDefault(TRUE) #set to default values (initiate)
   gWidgets2::visible(setwin) <- TRUE
 } #end gui-function # f_markerSettings()
 
 f_quitproj = function(h,...) {
  ubool <- gWidgets2::gconfirm("Do you want to save project?",title="Quit Program",icon="info")
  if(ubool) {
    f_saveproj(h)
  } else { 
   print("Program terminated without saving")
  }
  gWidgets2::dispose(mainwin) #remove window!
 }

 #helpfunction to get value in from user and store
 setValueUser <- function(what1,what2,txt,allowNULL=FALSE,allowText=FALSE) {
   listopt <- get(what1,envir=mmTK) #get object what 1.
   val <- listopt[[what2]]
   if(is.null(val)) val ="" #gwidgets2 does not handle NULL, must use empty string
   sw <- gWidgets2::gwindow(title="User input",visible=FALSE, width=300,height=50)
   grid <- gWidgets2::glayout(spacing=0,container=sw )
   grid[1,1] <- gWidgets2::glabel(txt, container=grid)
   grid[1,2] <- gWidgets2::gedit(text=val,container=grid,width=30)
   grid[2,1] <- gWidgets2::gbutton("OK", container=grid,handler = function(h, ...) { 
    GUIval = gWidgets2::svalue(grid[1,2]) #obtain GUI value
    if(allowNULL && GUIval=="") { #if accepting empty string
      tmp = NULL #Insert NULL
    } else {
      tmp <- GUIval
      if(!allowText) {
        tmp <- as.numeric(GUIval) #insert new value
        if(is.na(tmp)) {
          NAerror(what2)
          return()
        }
      }
    }
    listopt[[what2]] <- tmp
    assign(what1,listopt,envir=mmTK) #assign user-value to opt-list
    gWidgets2::dispose(sw)
   } )
   grid[2,2] <- gWidgets2::gbutton("Cancel", container=grid,handler = function(h, ...) { gWidgets2::dispose(sw) } )
   gWidgets2::visible(sw) <- TRUE
 }

 #helpfunction to get value from user
 getValueUser <- function(txt="",val=0) {
  val2 <- gWidgets2::ginput(txt, text=val, title="User input",icon="question")
  return(val2)   
 }

 NAerror <- function(what) {
  gWidgets2::gmessage(paste0(what," must be specified as a valid value"),title="Wrong input",icon="error")
  #stop("Wrong user-input")
 }

 #helpfunction which checks that at value is in interval of [0,1]
 checkProb = function(x,what) {
  if(is.na(x)) NAerror(what)
  if(x < 0 || x>1) {
   gWidgets2::gmessage(paste0(what," must be specified in interval [0,1] "),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
 }
 checkPositive = function(x,what,strict=FALSE) {
  if(is.na(x)) NAerror(what)
  if(x < 0 ) {
   gWidgets2::gmessage(paste0(what," cannot be a negative number"),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
  if(strict && x==0) {
   gWidgets2::gmessage(paste0(what," cannot be zero"),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
 }
 checkPosInteger = function(x,what) {
  if(is.na(x)) NAerror(what)
  if(x < 1 || round(x)!=x) {
   gWidgets2::gmessage(paste0(what," must be a positive integer"),title="Wrong input",icon="error")
   stop("Wrong user-input")
  }
 }

 #helpfunction for printing evidence sample to terminal
 printEvid = function(subD) {
  locs <- names(subD) #get unique loci
  mixtab <- matrix(ncol=2,nrow=length(locs))
  for(loc in  locs) { #for each locus
        mixA <- subD[[loc]]$adata
        mixH <- subD[[loc]]$hdata
        if(!is.null(mixA)) mixtab[which(loc==locs),1] <- paste0(mixA ,collapse="/")
        if(!is.null(mixH)) mixtab[which(loc==locs),2] <- paste0(mixH ,collapse="/")
  }
  rownames(mixtab) <- locs
  colnames(mixtab) <- c("Allele","Height")
  .printTable(mixtab)
 }  

 #helpfunction for printing reference sample to terminal
 printRefs = function(refData,refSel=NULL) {
   popFreq <- get("popFreq",envir=mmTK) #get frequencies
   fst <- get("optSetup",envir=mmTK)$fst0 #get Fst-correction
   nRefs <- length(refSel) #number of selected references
   locs <- unique(unlist(lapply(refData,names))) #get unique loci
   reftab <- matrix(ncol=nRefs,nrow=length(locs)) #last row is RMP
   RMPs <- matrix(1,ncol=nRefs,nrow=2) #RMP and 1/RMP
   colnames(RMPs) <- refSel
   rownames(RMPs) <- c("RMP","log10(1/RMP)")
   for(rsel in refSel) {
    for(loc in  locs) { #for each locus
      refA <-refData[[rsel]][[loc]]$adata
      if(!is.null(refA)) {
       reftab[which(loc==locs),which(rsel==refSel)] <- paste0(refA ,collapse="/")
       if(!is.null(popFreq) && loc%in%names(popFreq) ) {
        tmp <- popFreq[[loc]][ names(popFreq[[loc]])%in%refA ]
        prob <- 0 #default: Some of alleles not found in freqtable
        if(length(refA)!=2) next #skip if not 2 alleles
        if(refA[1]==refA[2]) { #if homozygous
         if(length(tmp)==1) prob <- (2*fst+(1-fst)*tmp)/(1+fst)*(3*fst+(1-fst)*tmp)/(1+2*fst) #tmp^2 
        } else {#if heterozygous
         if(length(tmp)==2) prob <- 2*(fst+(1-fst)*tmp[1])/(1+fst)*(fst+(1-fst)*tmp[2])/(1+2*fst) #2*prod(tmp) 
        }
        RMPs[1,which(refSel==rsel)] <- RMPs[1,which(refSel==rsel)]*prob 
       }
      }
    }
   }
   RMPs[2,] <- round(log10(1/RMPs[1,]),2) #log10 of invert RMP (LR)
   RMPs[1,] <- signif(RMPs[1,],2)
   
   rownames(reftab) <- locs
   colnames(reftab) <- refSel 
   .printTable(reftab)
   if(!is.null(popFreq)) {
     print(paste0("Calculation of random match probability and its inverse for fst=",fst))
     .printTable(RMPs)
     misslocs <- setdiff(locs,names(popFreq)) #missing loci compared to refs loci
     if(length(misslocs)>0) print( paste0("RMP not using loci ", paste0(misslocs,collapse="/") ))
   }
 }

 #Helpfunction to load data in a file with values and fit drop-in model
 f_fitdropin = function(h,...) { 
   impfile = efm_gfile(text="Find drop-in data",type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab")),"all"=list(patterns=list("*"))))
   if(length(impfile)==0) return()
   data = tableReader(impfile,header=FALSE) #load profile
   x = as.numeric(unlist(data)) #convert to a numerical vector
   x = x[!is.na(x)] #dont consider NAs

   threshT <- get("optSetup",envir=mmTK)$thresh0 #get specified threshold from settings
   x = x[x>=threshT] #UPDATED v2.1.1: consider dropinPH above threhsold
   if(length(x)==0) { #if no drop-in observed
    gWidgets2::gmessage("No PHs above detection threshold found. \nNo estimation could be done!",title="Drop-in model estimation", icon="info")
    return()
   }
   
   print(paste0(x,collapse=","))
   lambda = signif(length(x)/sum(x-threshT),2) #estimated hyperparam with 2 significant numbers
   print(paste0("Estimated hyperparameter lambda=",lambda))

   optSetup <- get("optSetup",envir=mmTK) #get settings
   optSetup$lam0=lambda #update parameter in settings
   assign("optSetup",optSetup,envir=mmTK)  #store settings again
    
   #PLOT DROPIN MODEL:
   par(mfrow=c(1,1))
   h = hist(x,breaks=25,probability=TRUE,xlab="RFU",main=paste0("Estimated lambda=",lambda))
   lines(h$mid,dexp(h$mid-threshT ,rate=lambda ),lty=1,lwd=2)
   legend("topright",legend=paste0("Drop-in P.H. distribution\nw/lambda=",lambda),lty=1)
 }

 
##################################### 
###########GUI WINDOW STARTS#########
##################################### 
suppressWarnings({
 
 ##########
 #Menu bar#
 ##########
 mblst = list( #project saving and so on
  File=list(  
    gWidgets2::gaction('Set directory',handler=f_setwd),
    gWidgets2::gaction('Open project',handler=f_openproj),
    gWidgets2::gaction('Save project',handler=f_saveproj),
    gWidgets2::gaction('Settings',handler=f_settings),
    gWidgets2::gaction('Marker settings',handler=f_markerSettings), #updated v3.0.0
    gWidgets2::gaction('Quit',handler=f_quitproj,icon="close")
  ),
  Frequencies=list(
    gWidgets2::gaction('Set size of frequency database',handler=function(h,...) {  
      setValueUser(what1="optFreq",what2="freqsize",txt="Set size of imported freq database \n(Min observed used if not spesified):") 
    }),
    gWidgets2::gaction('Set minimum frequency',handler=function(h,...) {  
      setValueUser(what1="optFreq",what2="minF",txt="Set minimum freq for new alleles\n(Min observed used if not spesified):") 
    }),
    gWidgets2::gaction('Set whether to normalize frequencies',handler=function(h,...) {  
      setValueUser(what1="optFreq",what2="normalize",txt="Should frequencies always add up to one\nafter including rare alleles? (1=YES,0=NO)") 
    }),
#   'Select stutter model',handler=function(h,...) {  
#      setValueUser(what1="optFreq",what2="stutterMod",txt="Select prefered stutter model \n(1=Before v2,2=From v2, <v1.12 used 1)") 
#    }),
    gWidgets2::gaction('Set number of wildcards in false positives match',handler=function(h,...) {  
      setValueUser(what1="optFreq",what2="wildsize",txt="Set number of missmatches (wildcards) in false positive match:") 
    }),
    gWidgets2::gaction('Set URL for STRidER import',handler=function(h,...) {  
      setValueUser(what1="STRidER",what2="url",txt="Set URL path:",allowText=TRUE) 
    })
  ),
  Optimization=list(
    gWidgets2::gaction('Set number of successful optimizations',handler=function(h,...) {  
      setValueUser(what1="optMLE",what2="nDone",txt="Set required number of (identical) optimizations:") 
    }),
    gWidgets2::gaction('Set variation of randomizer',handler=function(h,...) {  
      setValueUser(what1="optMLE",what2="delta",txt="Set variance of start point of MLE randomizer:") 
    }),
    gWidgets2::gaction('Set difference tolerance',handler=function(h,...) {  
      setValueUser(what1="optMLE",what2="difftol",txt="Set tolerance of successive MLEs (logLik):") 
    }),
    gWidgets2::gaction('Set seed of randomizer',handler=function(h,...) { 
      setValueUser(what1="optMLE",what2="seed",txt="Set seed of start value randomizer:",allowNULL=TRUE) 
    }),
    gWidgets2::gaction('Set accuracy of optimization',handler=function(h,...) { 
     setValueUser(what1="optMLE",what2="steptol",txt="Set accuracy of MLE optimization (steptol, see ?nlm):") 
    }),
    gWidgets2::gaction('Set significance level of validation',handler=function(h,...) { 
      setValueUser(what1="optMLE",what2="alpha",txt="Set significance level for model validation:") 
    }),
    gWidgets2::gaction('Set maximum threads for computation',handler=function(h,...) {  
      setValueUser(what1="optMLE",what2="maxThreads",txt="Set max number of threads used for parallelisation \n(0 means using all threads):") 
    })
  ),
  MCMC=list(
    gWidgets2::gaction('Set number of samples',handler=function(h,...) {  
      setValueUser(what1="optMCMC",what2="niter",txt="Set number of sample iterations:") 
    }),
    gWidgets2::gaction('Set variance of randomizer',handler=function(h,...) {  
      setValueUser(what1="optMCMC",what2="delta",txt="Set variation of randomizer:") 
    }),
    gWidgets2::gaction('Set quantile',handler=function(h,...) {  
      setValueUser(what1="optMCMC",what2="quantile",txt="Set quantile of conservative LR:") 
    }),
    gWidgets2::gaction('Set seed of randomizer',handler=function(h,...) {
      setValueUser(what1="optMCMC",what2="seed",txt="Set seed of MCMC randomizer:",allowNULL=TRUE) 
    })
  ),
  Integration=list(
    gWidgets2::gaction('Set relative error requirement',handler=function(h,...) {  
      setValueUser(what1="optINT",what2="reltol",txt="Set relative error:") 
    }),
    gWidgets2::gaction('Set maximum number of evaluations',handler=function(h,...) {  
      setValueUser(what1="optINT",what2="maxeval",txt="Set maximum number of evaluations for calculating integral:") 
    }),
    gWidgets2::gaction('Set deviation scale',handler=function(h,...) {  
      setValueUser(what1="optINT",what2="dev",txt="Set deviation scale for obtaining parameter boundaries (see getParamLimits)") 
    })
  ),
  Deconvolution=list(
    gWidgets2::gaction('Set required summed probability',handler=function(h,...) {  
      setValueUser(what1="optDC",what2="alphaprob",txt="Set required summed posterior genotype-probability of list:") 
    }),
    gWidgets2::gaction('Set max listsize',handler=function(h,...) {  
      setValueUser(what1="optDC",what2="maxlist",txt="Set size of maximum elements in deconvoluted list:") 
    })
  ),
  'Database search'=list(
    gWidgets2::gaction('Set maximum view-elements',handler=function(h,...) {  
      setValueUser(what1="optDB",what2="maxDB",txt="Set max size of view elements from database:") 
    }),
    gWidgets2::gaction('Set drop-in probability for qualitative model',handler=function(h,...) {  
      setValueUser(what1="optDB",what2="QUALpC",txt="Set allele drop-in probability for qualitative model:") 
    }),
    gWidgets2::gaction('Set number of non-contributors',handler=function(h,...) {  
      setValueUser(what1="optDB",what2="ntippets",txt="Set number of non-contributor samples in non-contributor plot:") 
    })
  ),
  'Qual LR'=list(
    gWidgets2::gaction('Set upper range for sensitivity',handler=function(h,...) {  
      setValueUser(what1="optLRMIX",what2="range",txt="Set upper limit of dropout in sensitivity analysis:") 
    }),
    gWidgets2::gaction('Set nticks for sensitivity',handler=function(h,...) {  
      setValueUser(what1="optLRMIX",what2="nticks",txt="Set number of ticks in sensitivity analysis:") 
    }),
    gWidgets2::gaction('Set required samples in dropout distr.',handler=function(h,...) {  
      setValueUser(what1="optLRMIX",what2="nsample",txt="Set required number number of samples in dropout distribution:") 
    }),
    gWidgets2::gaction('Set significance level in dropout distr.',handler=function(h,...) {  
      setValueUser(what1="optLRMIX",what2="alpha",txt="Set significance level (quantiles) in dropout distribution:") 
    })
  )
 )

##################################################################################################
########### Program starts #######################################################################
##################################################################################################

 #change working directory to the one stored in mmTK-environment
 wd=get("workdir",envir=mmTK) #assign working directory to mmTK-environment
 if(!is.null(wd) && file.exists(wd) ) setwd(wd)
 
 #Main window:
 mainwin <- gWidgets2::gwindow(softname, visible=FALSE, width=mwW,height=mwH)
 gWidgets2::gmenu(mblst,container=mainwin)
 nb = gWidgets2::gnotebook(container=mainwin)
 tabGEN = gWidgets2::glayout(spacing=spc,container=nb,label="Generate data") #tab1: Generates data(with peak heights) for a given model (plot EPG in addition)
 tabimport = gWidgets2::ggroup(horizontal=FALSE,spacing=10,container=nb,label="Import data") #tab2: (imports all files)
 tabmodel = gWidgets2::glayout(spacing=spc,container=nb,label="Model specification") #tab3: specify model used in weight-of-evidence (INT/MLE) or in a Database search 
 tabMLE = gWidgets2::glayout(spacing=spc,container=nb,label="MLE fit")#,expand=TRUE,fill=TRUE) #results from MLE
 tabDC = gWidgets2::ggroup(horizontal=FALSE,spacing=spc,container=nb,label="Deconvolution") #results from a deconvolution
 tabDB= gWidgets2::ggroup(horizontal=FALSE, spacing=spc,container=nb,label="Database search") #results from a database search
 tabLRMIX <- gWidgets2::glayout(spacing=spc,container=nb,label="Qual. LR") #LRmix


######################################################
###############Tab 1: Generate Data:##################
######################################################

 refreshTabGEN = function(thlist=list(mu=1000,sigma=0.15,xi=0.1,beta=1,mx=NULL,xiFW=0) ) { #can be repeated
  gWidgets2::visible(mainwin) <- FALSE 
  tabGENtmp <- gWidgets2::glayout(spacing=0,container=(tabGEN[1,1,expand=TRUE] <- gWidgets2::ggroup(container=tabGEN))) 
  
  #load/save helpfunctions for generated data
  f_openprof = function(h,...) {
    proffile = efm_gfile(text="Open profile",type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab")),"all"=list(patterns=list("*"))))
    if(length(proffile)==0) return() 
     Data = tableReader(proffile) #load profile
     .printTable(Data)
     setDataToGUI(sample_tableToList(Data)[[1]] ,h$action) #convert from table to list and load into GUI
    
  }
  f_saveprof = function(h,...) {
    Data <- list(getDataFromGUI(h$action)) #get data from GUI
    if(h$action==0) names(Data) <- paste0("stain",sample(100,1))
    if(h$action>0) names(Data) <- paste0("ref",h$action)
    Data = sample_listToTable(Data) #convert from table to list
    tableSaver(Data,"csv")
  }

  #helpfunction to set Data to GUI
  setDataToGUI <- function(Data,type) {
    dloc <- names(Data) #locus in Data
    for(i in 1:length(locs)) {
     dind <- grep(toupper(locs[i]),toupper(dloc)) #get dataindex
     if(type==0) { #if evidence
       gWidgets2::svalue(tabGENb2[i,1]) <- gWidgets2::svalue(tabGENb2[i,2]) <- ""
       if(length(dind)>0) { #if locus found
         if(!is.null(Data[[dind]]$adata)) gWidgets2::svalue(tabGENb2[i,1]) <- paste0(Data[[dind]]$adata,collapse=",") #insert alleles
         if(!is.null(Data[[dind]]$hdata)) gWidgets2::svalue(tabGENb2[i,2]) <- paste0(Data[[dind]]$hdata,collapse=",") #insert alleles
       }
     } else { #if reference
       gWidgets2::svalue(tabGENb3[i,type]) <- ""
       if( length(dind)>0 && !is.null(Data[[dind]]$adata) ) gWidgets2::svalue(tabGENb3[i,type]) <- paste0(Data[[dind]]$adata,collapse=",") #insert
     }
   } #end for each locus
  } #end function

  #helpfunction to get data from GUI (type=0 for evidence, otherwise its reference)
  getDataFromGUI <- function(type) {
    outloc <- locs #store locs
    Data <- list()
    for(i in 1:length(locs)) {
     outloc[i] <- gWidgets2::svalue(tabGENb1[i,1]) #get new loc-names
     if(type==0) { #store evidence
      Data[[outloc[i]]] <- list( adata=unlist(strsplit(gWidgets2::svalue(tabGENb2[i,1]),",")) , hdata=as.numeric(unlist(strsplit(gWidgets2::svalue(tabGENb2[i,2]),","))) )
     } else { #store reference
      Data[[outloc[i]]] <- list( adata=unlist(strsplit(gWidgets2::svalue(tabGENb3[i,type]),",")) )
     }    
    } #end for each locus
    return(Data)
  }

  #layout: 
  tabGENtop = gWidgets2::glayout(spacing=3,container=(tabGENtmp[2,1] <-gWidgets2::gframe(spacing=3,container=tabGENtmp))) 
  tabGENa = gWidgets2::glayout(spacing=0,container=(tabGENtop[1,1] <-gWidgets2::gframe("Parameters",container=tabGENtop))) 
  tabGENd = gWidgets2::glayout(spacing=3,container=(tabGENtop[2,1] <-gWidgets2::gframe("Further action",container=tabGENtop))) 
  #tabGENc = gWidgets2::glayout(spacing=3,container=(tabGENtop[3,1] <-gWidgets2::gframe("Import/Export profile",container=tabGENtop)))  

  #layout of data  
  tabGENc = gWidgets2::glayout(spacing=3,container=(tabGENtmp[1,2] <-gWidgets2::gframe("Import/Export profile",container=tabGENtmp)))  
  tabGENb = gWidgets2::glayout(spacing=0,container=(tabGENtmp[2,2] <-gWidgets2::gframe("Edit",container=tabGENtmp))) 
  
  #number of contributors
  set <- get("setGEN",envir=mmTK) #get stored setup
  nC <- set$nC_hd #number of contributors
 
  #default values
  if(is.null(thlist$mx)) thlist$mx <- round((nC:1)/sum(nC:1),3) #default value

  width=8 #fixed width of edit
  #user input:
  tabGENa[1,1] <- gWidgets2::glabel("P.H.expectation",container=tabGENa)
  tabGENa[1,2] <- gWidgets2::gedit(thlist$mu,container=tabGENa)
  tabGENa[2,1] <- gWidgets2::glabel("P.H.variability",container=tabGENa)
  tabGENa[2,2] <- gWidgets2::gedit(thlist$sigma,container=tabGENa)
  tabGENa[3,1] <- gWidgets2::glabel("Degrad.slope",container=tabGENa)
  tabGENa[3,2] <- gWidgets2::gedit(thlist$beta,container=tabGENa)
  tabGENa[4,1] <- gWidgets2::glabel("BW stutter-prop.",container=tabGENa)
  tabGENa[4,2] <- gWidgets2::gedit(thlist$xi,container=tabGENa)
  tabGENa[5,1] <- gWidgets2::glabel("FW stutter-prop.",container=tabGENa)
  tabGENa[5,2] <- gWidgets2::gedit(thlist$xiFW,container=tabGENa)
  gWidgets2::size(tabGENa[1,2]) <- gWidgets2::size(tabGENa[2,2]) <- gWidgets2::size(tabGENa[3,2]) <- gWidgets2::size(tabGENa[4,2]) <- gWidgets2::size(tabGENa[5,2]) <- width
  
  for(k in 1:nC) { #for each contributors
   tabGENa[k+5,1] <- gWidgets2::glabel( paste0("Mix-prop.",k," (mx",k,")") ,container=tabGENa)
   tabGENa[k+5,2] <- gWidgets2::gedit(thlist$mx[k],container=tabGENa)
   gWidgets2::size(tabGENa[k+5,2]) = width
  }

  #Generate:
  kit <- get("selPopKitName",envir=mmTK)[1]  #get selected kit
  if(is.na(kit)) kit = NULL
  simdata <- genDataset(nC, popFreq=set$popFreq, mu=thlist$mu, sigma=thlist$sigma,beta=thlist$beta,sorted=FALSE,threshT=set$threshT, refData=set$refData, mx=thlist$mx/sum(thlist$mx),nrep=1, stutt = thlist$xi, prC=set$prC,lambda=set$lambda,kit=kit, stuttFW=thlist$xiFW)

  #insert data in GUI
  mixData <- simdata$samples[[1]]
  refData <- simdata$refData
  locs <- names(mixData) #get locus names

  #show Loci
  tabGENb1 <-  gWidgets2::glayout(spacing=0,container=(tabGENb[1,1] <-gWidgets2::gframe("Loci",container=tabGENb))) 
  for(i in 1:length(locs))  {
    tabGENb1[i,1] = gWidgets2::gedit(locs[i],container=tabGENb1)
    gWidgets2::size(tabGENb1[i,1]) <- nchar(locs[i])+3
  }

  #show allele,heights
  tabGENb2 <-  gWidgets2::glayout(spacing=0,container=(tabGENb[1,2] <-gWidgets2::gframe("Evidence (allele,heights)",container=tabGENb))) 
  for(i in 1:length(locs)) {
   adata <- mixData[[locs[i]]]$adata
   hdata <- round(mixData[[locs[i]]]$hdata)
   tabGENb2[i,1] = gWidgets2::gedit(paste0(adata,collapse=","),container=tabGENb2)
   tabGENb2[i,2] = gWidgets2::gedit(paste0(hdata,collapse=","),container=tabGENb2,width=sum(nchar(hdata)) + length(hdata))
   gWidgets2::size(tabGENb2[i,1]) <- sum(nchar(adata)) + length(adata)
   gWidgets2::size(tabGENb2[i,2]) <- sum(nchar(hdata)) + length(hdata)
  }

  #show references:
  tabGENb3 <-  gWidgets2::glayout(spacing=0,container=(tabGENb[1,3] <-gWidgets2::gframe("Reference(s)",container=tabGENb))) 
  for(k in 1:nC) {
   for(i in 1:length(locs)) {
    adata <- refData[[locs[i]]][[k]]
    tabGENb3[i,k] = gWidgets2::gedit(paste0(adata,collapse=","),container=tabGENb3)
    gWidgets2::size(tabGENb3[i,k]) <- sum(nchar(adata)) + length(adata)
   }
  }

  #storage buttons:
  tabGENc[1,1] <- gWidgets2::gbutton(text="Store evidence",container=tabGENc,handler=f_saveprof,action=0)
  for(k in 1:nC) tabGENc[1,1+k] <- gWidgets2::gbutton(text=paste0("Store ref",k),container=tabGENc,handler=f_saveprof,action=k)
  tabGENc[2,1] <- gWidgets2::gbutton(text="Load evidence",container=tabGENc,handler=f_openprof,action=0)
  for(k in 1:nC) tabGENc[2,1+k] <- gWidgets2::gbutton(text=paste0("Load ref",k),container=tabGENc,handler=f_openprof,action=k)

  #further action
  tabGENd[1,1] <- gWidgets2::gbutton(text="Generate again",container=tabGENd,handler=function(h,...) { 
   mx <- rep(0,nC) #mixture proportions
   for(k in 1:nC)  mx[k] <- as.numeric(gWidgets2::svalue(tabGENa[5+k,2])) 
   mu0 = as.numeric(gWidgets2::svalue(tabGENa[1,2])) #exp. PH
   sigma0 = as.numeric(gWidgets2::svalue(tabGENa[2,2])) #var. PH
   beta0 = as.numeric(gWidgets2::svalue(tabGENa[3,2])) #degrad param
   xiBW = as.numeric(gWidgets2::svalue(tabGENa[4,2])) #BW stutter prop.
   xiFW = as.numeric(gWidgets2::svalue(tabGENa[5,2])) #FW stutter prop.
   refreshTabGEN(list(mx=mx/sum(mx),mu=mu0,sigma=sigma0,beta=beta0,xi=xiBW,xiFW=xiFW) ) #refresh GUI
  })

  tabGENd[2,1] <- gWidgets2::gbutton(text="Plot EPG",container=tabGENd,handler=function(h,...) {
   kit <- get("selPopKitName",envir=mmTK)[1] #get
   if(is.null(kit) || is.na(kit)) {
     gWidgets2::gmessage("Select a valid kit to show data")
     return()
   }
   Data <- getDataFromGUI(0) #get data from GUI
   if(requireNamespace("plotly")) { #visualise epg with reference
     refData = list()
     if(nC>0) {
       for(k in 1:nC) refData[[ paste0("ref",k) ]] <- getDataFromGUI(k) #get reference data from GUI
     } else {
       refData = NULL  
     }
     plotEPG2(list(stain=Data),kit,refData,AT=set$threshT)  #Plot epg if kit was recognized
   } else { #shown if plotly not found
     plotEPG(list(stain=Data),kitname=kit,threshT=0) #plot epg's
   }
   gWidgets2::focus(mainwin) <- TRUE
  })

  gWidgets2::visible(mainwin) <- TRUE #INCREASE SPEED
  gWidgets2::focus(mainwin) <- TRUE
 } #end refreshTabGE

####################################################
###############Tab 2: Import Data:##################
####################################################

  #helpfunction to delete selected profiles 
  f_deleteprof = function(h,...) {
    type=h$action #get type of profile
    col = 1 + 1*sum(type=="ref") + 2*sum(type=="db")  #column in table to use
    sel <- gWidgets2::svalue(tabimportB[3,col])  #get selected profiles (can be several)
    Data2 <- getData(type) #get data from mmTK-environment
    if(length(Data2)==0 || length(sel)==0)  return() #return if no data to drop
    keep = setdiff(names(Data2),sel)
    Data2 = Data2[keep]
    tab = cbind(Profiles=keep) #table to insert
    if(type=="db") tab = cbind(Databases=keep) #table to insert
    tabimportB[3,col][] <- tab #update
    assign(paste0(type,"Data"),Data2,envir=mmTK) #assign data to mmTK-environment again    
  }
  
 
 #b) load/save profiles/database: Supports any filetype
 f_importprof = function(h,...) {
  type=h$action #get type of profile
#  proffile = efmgfile(text=paste("Open ",type,"-file",sep=""),type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab"))))
  proffile = efm_gfile(text=paste("Open ",type,"-file",sep=""),type="open",filter=list("text"=list(patterns=list("*.txt","*.csv","*.tab")),"all"=list(patterns=list("*"))))
  if(length(proffile)==0) return() 
  
   Data = tableReader(proffile) #load profile
   #print("Raw fil import:")
   #.printTable(Data[1:min(nrow(Data),100),]) #print raw-input data
   
  ###################
  ##DATABASE IMPORT##
  ###################
   if(type=="db") { #importing database file
    prim = .getPrim()  #get prime numbers. max length of alleles for reference-database storage is 244
     
    popFreq <- get("popFreq",envir=mmTK) 
    if(is.null(popFreq)) {
     gWidgets2::gmessage("Population frequencies needs to be imported for database search",title="Error")
    } else {
     minFreq <- getminFreq()
     #saving MEMORY by convert database here!
     #Data <- dbToGenoTable(Data) #convert database table to genotype matrix
     
     #Assumption: No reference has same samplename. All samplenames are rowed sequentially
     cn = colnames(Data) #colnames 
     lind = grep("marker",tolower(cn),fixed=TRUE) #locus col-ind
     sind = grep("sample",tolower(cn),fixed=TRUE) #sample col-ind
     aind = grep("allele",tolower(cn),fixed=TRUE) #allele col-ind
     aind <- aind[1:2] #assume only 2 possible alles in reference profiles
     locsDB <- unique(Data[,lind]) #unique markers in DB
     locsPop <- toupper(names(popFreq)) #markers in population
     sn <- unique(Data[,sind]) #unique samples in DB
     newData <- numeric() #encoded reference table
     dbDatalocs <- numeric() #locus in newData
     for(loc in locsDB) { #for each marker in DB:
      locind <- grep(toupper(loc),toupper(names(popFreq)),fixed=TRUE) #get position in popFreq
      if(length(locind)==0) next #if locus not existing in popFreq the locus in not imported
      newCol <- rep(NA,length(sn)) #new column in newData
      subA <- Data[Data[,lind]==loc,aind] #extract allele data with matching marker
      subS <- Data[Data[,lind]==loc,sind] #extract sample names with matching marker
      isHom <- which(subA[,2]=="" | is.na(subA[,2])) #if homozygote assignation:
      if(length(isHom)>0) subA[isHom,2] <- subA[isHom,1] #copy first allele
      okSind <- which(sn%in%subS) #samples to consider for particular marker
      if(length(okSind)==0) next #if no samples exists
      print(paste0("Import data from locus: ",loc)) 
      newCol[okSind] <- 1 #init existing as 1. NA for missing allele-info
      doneEncode <- matrix(FALSE,ncol=2,nrow=length(okSind)) #matrix checking which we are finished with encoding a genotype
      Afreqs <- names(popFreq[[locind]]) #get allele-names. Update for each column
      for(k in 1:length(aind)) { #for each allele-column
       for(j in 1:length(Afreqs)) { #for each unique allele in popFreq:
        okAind <- which(subA[,k]==Afreqs[j]) #find matching alleles in subA[okSind]
        if(length(okAind)==0) next
        doneEncode[okAind,k] = TRUE #check that is finished
        newCol[okSind][okAind] <- newCol[okSind][okAind]*prim[j] #multiply with primenumber
       } #end for each allele j

       #THREAT NEW ALLELES,MISSTYPOS ETC:
       if(any(!doneEncode[,k])) { #if any not encoded
        newA <- unique(subA[!doneEncode[,k],k]) #get new alleles
        newA <- newA[!is.na(newA)] #important to remove NA's
        if(length(newA)==0) next
        tmp <- popFreq[[locind]]
        popFreq[[locind]] <- c(tmp, rep(minFreq,length(newA)))
        names(popFreq[[locind]]) <- c(names(tmp),newA) #add unique
        print(paste0("In locus ",loc,": Allele(s) ",newA," was inserted with min. frequency ",prettyNum(minFreq)))
        for(j in 1:length(newA)) { #for each unique allele in popFreq:
         okAind <- which(subA[,k]==newA[j]) #find matching alleles in subA[okSind]
         if(length(okAind)==0) next
         newCol[okSind][okAind] <- newCol[okSind][okAind] * prim[ which(names(popFreq[[locind]])==newA[j]) ] #multiply with EXTENDED primenumber
        } #end for each allele j
       } #end if not encoded 
       Afreqs <- names(popFreq[[locind]]) #get allele-names again!
      } #end for each column k
      dbDatalocs <- c(dbDatalocs,toupper(names(popFreq)[locind])) #all locus
      newData <- cbind(newData,newCol) #add column
     } #end for each locus loc
   
     #RESCALE popFreq?
     for(i in 1:length(popFreq)) {
      popFreq[[i]] <- popFreq[[i]]/sum(popFreq[[i]])
     }
     print("Frequencies was normalized")
     colnames(newData) <- dbDatalocs
     rownames(newData) <- sn #insert sample names
     print(paste0("Database successfully imported with ",nrow(newData)," samples."))

     tmp <- unlist(strsplit(proffile,"/",fixed=TRUE)) #just label the file
     fname <- tmp[length(tmp)] #get filename
     fname <- unlist(strsplit(fname,"\\."))[1]

     dbData <- get("dbData",envir=mmTK) #get already existing databases
     if(is.null(dbData)) dbData <- list() #create list-object
     dbData[[fname]] <- newData #add database to list

     assign("dbData",dbData,envir=mmTK) #store matrix in environment for later use
     assign("popFreq",popFreq,envir=mmTK) #assign updated popFreq
     
     tabimportB[3,3][] <- cbind(Databases=names(dbData))
    } #end if popFreq exist
    
   } else { #if not DB
     #ORDINARY IMPORT (EVID OR REF)
    Data = sample_tableToList(Data) #convert from table to list 

    #Block updated in v4.0.0: Check for data
    tooMany <- FALSE #indicator whether OL allele was observed
    hasOL <- FALSE #indicator of whether reference has more than 2 alleles
    for(kn in names(Data)) { #for each profile
      for(loc in names(Data[[kn]])) { #for each profile
        av = Data[[kn]][[loc]]$adata #obtain alleles (som may be missing)
        if(length(av)>0) { #need observations
          indOL = av%in%"OL"
          if(any(indOL)) {
            hasOL = TRUE #flag
            av = av[!indOL] #update alleles
            if(type=="mix") {
               Data[[kn]][[loc]]$adata = av  #update alleles
               Data[[kn]][[loc]]$hdata = Data[[kn]][[loc]]$hdata[!indOL] #update peak heights
            }
          }
        }
        if(type=="ref") { #special handling refs
          if( length(av)==1) Data[[kn]][[loc]]$adata = rep(av,2) #duplicated for homozygous
          if( length(av)>2 ) tooMany = TRUE #flag
          Data[[kn]][[loc]]$hdata = NULL #ensure that hdata is not given
        }
      }
    }
    if(tooMany) {
      gWidgets2::gmessage("Too many alleles where given for a genotype in a reference profile. Must be a maximum of 2.","Error in data",icon="error")
      return(NULL)
    }
    if(hasOL) { 
      bool = gWidgets2::gconfirm("OL allele observed in any of the markers. Please check data and import again. Do you still want to import the data where these are automatically removed?","Error in data")
      if(!bool) return(NULL)
    }

    #get already stored data:
    Data2 <- getData(type) #get data from mmTK-environment

    if(is.null(Data2)) { #if no previous already there
     Data2 <- Data
    } else {
     for(kn in names(Data)) { #for each profile
      Data2[[kn]] <- Data[[kn]] #insert dataframe
     }
    }
    assign(paste0(type,"Data"),Data2,envir=mmTK) #assign data to mmTK-environment
    newTab = cbind(Profiles=names(Data2)) #table to insert
    tabimportB[3,1 + sum(type=="ref")][] <- newTab
   }
 }

 #Helpfunction to assign population when clicked in list
 f_selectedPop = function(h,...) {
   if(!is.null(get("dbData",envir=mmTK))) {
     print("You can not change selected population after loading a database!") 
     return()
   }
   
   dbList <- get("popList",envir=mmTK) #get population frequency info from mmTK-environment
   if(is.null(dbList)) return()
   
   popsel <- gWidgets2::svalue(tabimportA[3,2]) #get selected population
   popList <- dbList[[popsel]] #get selected frequencies
   assign("popFreq",popList,envir=mmTK) #assign popFreq get("popFreq",envir=mmTK)
   
   popkitname <- get("selPopKitName",envir=mmTK) 
   if(is.null(popkitname)) popkitname <- rep(NA,2) #initiate kit-pop selection
   popkitname[2] <- popsel #insert selected population
   assign("selPopKitName",popkitname,envir=mmTK)  #store selected popname
   
   if(!is.null(popList)) { #NEW in v1.9: Gives a warning when allele frequencies are 'odd'
     locs <- names(popList)
     check1 <- locs[sapply(popList,length)==0] #loci which was empty
     check2 <- locs[sapply(popList,function(x) any(x<0))] #any loci containing negative frequencies
     check3 <- locs[sapply(popList,function(x) round(sum(x),2)!=1) ] #any loci containing negative frequencies
     txt <- paste0("Selected population frequency file: ",paste0(popkitname,collapse="-"))
     
     giveWarning = function(loci,what) {
       loctext <- ifelse(length(loci)>1,"Loci: ","Locus: ")
       txt2 <- paste0(txt,"\n\n",loctext, paste0(loci,collapse="/"),"\n",what,".\n\nPlease check the allele frequency file!")
       gWidgets2::gmessage(txt2,title="Incorrect allele frequencies observed!",icon="warning")
     }
     if(length(check1)>0)  giveWarning(check1,what="did not have any allele frequencies")
     if(length(check2)>0)  giveWarning(check2,what="had atleast one negative allele frequency")
     if(length(check3)>0)  giveWarning(check3,what="had summed allele frequencies different from one (with 2 decimal roundoff)")
   }
 }
 
 #Helpfunction to export frequency data (only frequency file currently implemented)
 f_exportfreq = function(h,...) {
   popFreq <- get("popFreq",envir=mmTK) #get frequencies
   if(is.null(popFreq)) {
     gWidgets2::gmessage("Please import and select population frequencies!",icon="info")
     return()
   } else {
     nL <- length(popFreq)
     unAchr <- unique(unlist(lapply( popFreq,names) )) #get character alleles
     ord <- order(as.numeric(unAchr))  #all alleles must be able to convert to numeric
     unAchr <- unAchr[ord]  #make increasing order

     outtab = matrix("",ncol=nL,nrow=length(unAchr)) #table to save to file
     for(i in 1:nL) { #go through all markers
       freqs <- popFreq[[i]] #get frequencies
       outtab[ match(names(freqs),unAchr) ,i ] = freqs #insert frequencies
     } 
     outtab = cbind(unAchr,outtab)
     colnames(outtab) = c("Allele",names(popFreq)) #insert marker names
     tableSaver(outtab,"csv") #save table (with csv)
   }
 }

 #prints evidence, references, EPG, databases and population frequencies
 f_viewdata = function(h,...) {
  help_gmessage = function(what) gWidgets2::gmessage(paste0("Please import and select ",what),icon="info")
  
  #types: freq,mix,ref,db
  evidD <- getData("mix") #get selected data
  mixSel  <- numeric()
  if(!is.null(evidD)) mixSel <- gWidgets2::svalue(tabimportB[3,1])  #get selected mixtures
  threshT <- get("optSetup",envir=mmTK)$thresh0 #get specified threshold from settings

  if(h$action=="freq") { #prints frequencies
   wildsize <- get("optFreq",envir=mmTK)$wildsize
   popFreq <- get("popFreq",envir=mmTK) #get frequencies
   if(is.null(popFreq)) {
     help_gmessage("population frequencies.")
   } else {
    locs <- names(popFreq)
    nL <- length(locs)
    unAchr <- unique(unlist(lapply( popFreq,names) )) #get character alleles
    ord <- order(as.numeric(unAchr)) 
    unAchr <- unAchr[ord]  #make increasing order
    outD <- unAchr

    matsi <- matrix(NA,nrow=nL,ncol=length(mixSel))  #vector with summed alleles for each marker (used for random match prob calcs)
    rownames(matsi) <- locs
    colnames(matsi) <- mixSel
    for(i in 1:nL) {
     freqs <- popFreq[[i]] #get frequencies
     newRow <- rep(NA,length(unAchr))
     for(j in 1:length(freqs)) {
      rowind <- which(unAchr==names(freqs[j] ))
      newRow[rowind] <- freqs[j]
     }
     outD <- cbind(outD,newRow)
     for(msel in mixSel) {
      adata <- evidD[[msel]][[locs[i]]]$adata #selected samples
      if(length(adata)>0) matsi[i,which(msel==mixSel)] <- sum(freqs[names(freqs)%in%adata])
     } #end for each samples
    } #end for each locus

    #plot table with frequncies
    colnames(outD) = c("Allele",locs) 
    dbwin <- gWidgets2::gwindow("Population frequencies", visible=FALSE)#,height=mwH)
    tab <- gWidgets2::gdf(NAtoSign(outD) ,container=dbwin) #create table
    gWidgets2::visible(dbwin) <- TRUE

    #Calculate random match probability of matching for each samples
    for(msel in mixSel) { 
     cind <- which(msel==mixSel)
     print(paste0("Calculating false positive MAC probability for sample ",msel,"... "))
     si <- matsi[!is.na(matsi[,cind]),cind]
     macdist <- exactMACdistr(si)
     macdist <- tail(macdist,wildsize+1) #allow wildcardsize missmatches
     cumprob <- rev(cumsum(rev(macdist)))
     ymax <- max(cumprob)
     delta <- 0.03
     if(msel!=mixSel[1]) dev.new() #create new plot for next EPG
     barplot(cumprob,main="Random match probability having number of allele matches>=k",xlab="k number of allele matches",space=0,ylim=c(0,ymax+ymax*2*delta))
     text(1:length(cumprob)-0.5,cumprob+ymax*delta,labels=format(cumprob,digits=2))
     mtext(paste0("Sample: ",msel))
    }
   } #end if popFreq
  } #end freq

  if(h$action=="mix") { #prints EPG
    nS <- length(mixSel) #number of selected samples
    if(nS==0) {
     help_gmessage("evidence profile(s).")
     return() #return if no samples selected
    }
    refData <- getData("ref") #get selected references
    refSel <- numeric()
    if(!is.null(refData))  refSel <- gWidgets2::svalue(tabimportB[3,2])  #get selected references
    kitname <- get("selPopKitName",envir=mmTK)[1]  #gWidgets2::svalue(tabimportA[3,1]) #get kit name. Must be same name as in generateEPG
    #first: Print evidences:
    
    evidDsel <- list()
    for(msel in mixSel) { #PRINT DATA
      subD <- evidDsel[[msel]] <- evidD[[msel]] #selected samples
      print("------------------------------------")
      print(paste("Samplename: ",msel,sep=""))
      printEvid(subD)
    }
    plotSumPH(evidDsel,kitname) #Fitting SumPH and show

    #SHOW EPG   
    #determine whether it is EPG/MPS(LUS)/(strings): Check if "_" is used. Otherwise check with all alleles are strings
    sampletype = getSampleType(evidDsel,kit=kitname,LUSsymbol=LUSsymbol) #get sample type (EPG/LUS/MPS)
    hasPlotly = requireNamespace("plotly") #check if having plotly
    
    if(hasPlotly) { #if using browswer to plot profile
      switch(sampletype,
             "EPG" = plotEPG2(evidDsel,kitname,refData[refSel],AT=threshT),
             "LUS" = plotMPS2(evidDsel,refData[refSel],AT=threshT,grpsymbol=LUSsymbol),
             "MPS" = plotMPS2(evidDsel,refData[refSel],AT=threshT,grpsymbol=MPSsymbol))
    } else {
      switch(sampletype, 
             "EPG" = plotEPG(Data=evidDsel,kitname,refcond=refData[refSel],threshT=threshT),
             print("Data format not supported without plotly installed!"))
    }
  } #end show mix

  if(h$action=="ref") { #print tables only
     refData <- getData("ref")
     refSel <- numeric()
     if(!is.null(refData))  refSel <- gWidgets2::svalue(tabimportB[3,2])  #get selected references
     nRefs <- length(refSel)
     if(nRefs==0) {
       help_gmessage("reference profile(s).")
       return()
     }
     printRefs(refData,refSel)
     #second: Print #matches to selected samples: Condition on samples inside evids
     for(msel in mixSel) { #for each selected evidence 
      print(paste0("Number of matching alleles with samplename ",msel,":"))
      subD <- evidD[[msel]] #selected samples
      locs <- names(subD)
      tab <- matrix(NA,ncol=nRefs,nrow=length(locs))
      rownames(tab) <- locs
      colnames(tab) <- refSel 
      for(loc in  locs) { #for each locus
        if( length(grep("AM",loc))>0) next
        if( length(subD[[loc]]$adata)==0) next
        for(rsel in refSel) {
          refA <- refData[[rsel]][[loc]]$adata
          if(is.null(refA) || length(refA)==0) next #skip if no data
          tab[which(loc==locs),which(rsel==refSel)] <- sum(refA%in%subD[[loc]]$adata)
        }
      } 
      MAC <- colSums(tab,na.rm=TRUE) #remove NA's
      nLocs <- colSums(!is.na(tab))
      nMiss = 2*nLocs - MAC
      tab2 <- rbind(tab,MAC,nLocs,MatchRate=round(MAC/(2*nLocs),2),nMissing=nMiss)
      .printTable(tab2)

     } #end for each mix    
     print("------------------------------------")
   } #end if references

   if(h$action=="db") {  #view imported databases (reference)
    popFreq <- get("popFreq",envir=mmTK)
    dbSel <- gWidgets2::svalue(tabimportB[3,3])  #get selected database(s)
    if(length(dbSel)==0 || nchar(dbSel[1])==0) { #if none selected
      help_gmessage("reference database(s).")
      return()
    }
    prim = .getPrim()  #get prime numbers. max length of alleles for reference-database storage is 244
    dbData <- get("dbData",envir=mmTK) #load all DB data
    for(dsel in dbSel) { #for each selected databases
     subD <- dbData[[dsel]] #get selected data
     dblocs <- toupper(colnames(subD)) #get database locs
     outD <- rownames(subD) #will be character
     macD <- matrix(0,nrow=length(outD),ncol=length(mixSel)) #matching allele counter for each reference
     nlocs <- rep(0,length(outD))
     for(i in 1:length(dblocs)) { #for each locus in db     
      Ainfo <- names(unlist(popFreq[[dblocs[i]]])) #extract allele-info
      #translate to genotypes
      Pinfo <- prim[1:length(Ainfo)]
      G = t(as.matrix(expand.grid(rep(list(Ainfo,Ainfo )))))
      GP = t(as.matrix(expand.grid(rep(list(Pinfo,Pinfo )))))
      keep = GP[2,]>=GP[1,] #unique genotypes
      G0 <- G[,keep]  #store genotypes

      G <- paste0(G0[1,],paste0("/",G0[2,]))
      GP <- GP[,keep]  #store genotypes
      GP <- GP[1,]*GP[2,]
      newRow <- rep(NA,nrow(subD))  
      tmpmac <- matrix(0,nrow=length(GP),ncol=length(mixSel)) #genotype match for each samples
      for(j in 1:length(GP)) { #for each genotype
       #Always: Find corresponding genotype name by looking up stored primenumbers (same order as in popFreq!)
       rowind <- which(subD[,i]==GP[j]) #samples with this genotype
       if(length(rowind)==0) next
       newRow[rowind] <- G[j]
       #Optional (if mixtures selected): Count matching alleles:
       hasloc <- FALSE #total number locus checked over evidence
       for(msel in mixSel) { #for each selected mixture
        evid0 <- evidD[[msel]][[dblocs[i]]]$adata
        if(!is.null(evid0)) hasloc <- TRUE
        tmpmac[j,which(msel==mixSel)] <- sum(G0[,j]%in%evid0)
       } 
       if(ncol(macD)>0) { #if compared with evidence
        macD[rowind,] = macD[rowind,] + tmpmac[j,] #add match counter 
        nlocs[rowind] <- nlocs[rowind] + as.integer(hasloc)
       } #end if evidence had locus
      } #end for each genotype
      outD <- cbind(outD,newRow) #extend
     } #end for each locus
     colnames(outD) <- c("Reference",dblocs)
     nmax <- min(nrow(outD), get("optDB",envir=mmTK)$maxDB) #max limit of number in showing dropdown
     if(nrow(outD)>nmax ) print(paste0("Database contained more than ",nmax," samples. Showing first ",nmax," samples only!"))


     #if selected Mixtures: show ranked matches in database-reference:
     if(ncol(macD)>0) { #
      ord <- order(rowSums(macD),decreasing=TRUE)
      outD2 <- cbind(outD[ord,1],macD[ord,],nlocs[ord])
      colnames(outD2) <- c("Reference",mixSel,"nMarkers") 
      dbwin2 <- gWidgets2::gwindow(paste0("Number of sample matching alleles in references in database ",dsel), visible=FALSE)#,height=mwH)
      if(nmax<=1) gWidgets2::gtable(NAtoSign(outD2),container=dbwin2) #create table
      if(nmax>1) gWidgets2::gtable(NAtoSign(outD2[1:nmax,]),container=dbwin2) #create table
      gWidgets2::visible(dbwin2) <- TRUE
     }
     dbwin <- gWidgets2::gwindow(paste0("References in imported database ",dsel), visible=FALSE)#,height=mwH)
     if(nmax<=1) gWidgets2::gdf(NAtoSign(outD),container=dbwin) #create table
     if(nmax>1) gWidgets2::gdf(NAtoSign(outD[1:nmax,]),container=dbwin) #create table
     gWidgets2::visible(dbwin) <- TRUE        

    } #end for each databases
   } #end if db -case
 }  #end viewdata

 
###################
#start GUI layout:#
###################
 tabimportA = gWidgets2::glayout(spacing=5,container=gWidgets2::gframe("Step 1) Import and select Population frequencies",container=tabimport)) #kit and population selecter
 #tabimportA2 = gWidgets2::glabel("",container=tabimport) #evidence,ref dataframe
 tabimportB = gWidgets2::glayout(spacing=5,container=gWidgets2::gframe("Step 2) Import and select Evidence, Reference, Database",container=tabimport,expand=TRUE,fill=TRUE),expand=TRUE,fill=TRUE) #evidence,ref dataframe
 #tabimportB2 = gWidgets2::glabel("",container=tabimport) #evidence,ref dataframe
 tabimportC = gWidgets2::glayout(spacing=10,container=gWidgets2::gframe("Step 3) Select Interpretation",container=tabimport)) #Tasks button

 setPopsToList = function(dbList) {
   assign("popList",dbList,envir=mmTK) #assign population frequency info to mmTK-environment
   popNames = names(dbList)
   tabimportA[3,2][] <- popNames #update list
   gWidgets2::svalue(tabimportA[3,2]) <-  popNames[1] #select first population (this will be stored by using f_selectedPop)
 }

 #Choose box and import button
 tabimportA[1,1] = gWidgets2::gbutton(text="Import from file",container=tabimportA,handler=
  function(h,...) {
   f = efm_gfile(text="Select file",type="open",filter = list(`All files` = list(patterns = c("*"))))
   if(length(f)==0) return()
    setPopsToList(freqImport(f,url=FALSE,xml=FALSE))
    f_selectedPop(h)    
  }
 )
 helptext(tabimportA[1,1],paste0("Choose a frequency file.\nExample files are found in directory\n",deffreq))
 
 tabimportA[1,2] = gWidgets2::gbutton(text="Import from Inst.",container=tabimportA,handler=
  function(h,...) {
    setPopsToList(freqImport(list.files(deffreq,full.names=TRUE),url=FALSE,xml=FALSE))
    f_selectedPop(h)
  }
 )
 helptext(tabimportA[1,2],paste0("Loading frequency files stored in ",deffreq))
 
 tabimportA[1,3] = gWidgets2::gbutton(text="Import from STRidER",container=tabimportA,handler=
  function(h,...) {
   striderlink1 = get("STRidER",envir=mmTK)$url #obtain strider link URL
   setPopsToList(freqImport(striderlink1,url=TRUE,xml=TRUE))
   f_selectedPop(h)
  }
 )
 helptext(tabimportA[1,3],paste0("Click to import directly from the STRidER database.\nLink: ",striderlink))

 tabimportA[2,1] <- gWidgets2::glabel(text="Select STR kit:",container=tabimportA)
 tabimportA[2,2] <- gWidgets2::glabel(text="Select population:",container=tabimportA)


 #previous stored kit/pop-names:
 initPops <- names(get("popList",envir=mmTK)) #get population list already saved from mmTK-environment
 initKits <- c("NONE",getKit()) #get kit names (shortname)
 selKit <- 0 #default selection
 selPop <- 0 #default selection
 if(is.null(initPops)) initPops <- ""
 popkitname <- get("selPopKitName",envir=mmTK) 
 if(!is.null(popkitname)) { 
  if(!is.na(popkitname[1])) {
   selKit <- which(popkitname[1]==initKits)
   if(length(selKit)==0) selKit <- 0
  }
  if(!is.na(popkitname[2])) {
   selPop <- which(popkitname[2]==initPops)
   if(length(selPop)==0) { #if not found, just show poptext (the popFreq object is already stored)
    selPop <- 1
    initPops <- popkitname[2]
   }
  }
 }

 tabimportA[3,1] <- gWidgets2::gcombobox(items=initKits, width=100, selected =selKit, editable = FALSE, container = tabimportA, handler=
    function(h,...) {
     if(!is.null(get("dbData",envir=mmTK))) {
      print("You can not change selected kit after loading a database!") 
     } else {
      kitsel <- gWidgets2::svalue(tabimportA[3,1]) #get selected kit
      popkitname <- get("selPopKitName",envir=mmTK) 
      if(is.null(popkitname)) popkitname <- rep(NA,2) #initiate kit-pop selection
      
      if(kitsel!="NONE") {
        popkitname[1] <- kitsel #insert selected kit (if not none)
      } else {
        popkitname[1] = NA #set to none
      }
      assign("selPopKitName",popkitname,envir=mmTK) #store to environment
     } #end else
    })
 
 #population-selection
 tabimportA[3,2] <- gWidgets2::gcombobox(items=initPops, width=100, selected = selPop, editable = FALSE , container = tabimportA, handler=f_selectedPop)

 tabimportA[2,3] <-  gWidgets2::gbutton(text="Export frequencies",container=tabimportA,handler=f_exportfreq)  #view popFreq-data
 helptext(tabimportA[2,3],"Exports a frequency file (in EuroForMix/LRmixStudio format) for selected population.")
 
 tabimportA[3,3] <-  gWidgets2::gbutton(text="View frequencies",container=tabimportA,handler=f_viewdata,action="freq")  #view popFreq-data
 helptext(tabimportA[3,3],"Shows the selected population frequencies from the drop-down menu. \n\nIf evidence(s) is selected, the probability for a random profile to have more than k number of allele matches to the sample is given in a plot.")

 #Set import buttons
 tabimportB[1,1] = gWidgets2::gbutton(text="Import evidence",container=tabimportB,handler=f_importprof,action="mix")
 tabimportB[1,2] = gWidgets2::gbutton(text="Import reference",container=tabimportB,handler=f_importprof,action="ref")
 tabimportB[1,3] = gWidgets2::gbutton(text="Import database",container=tabimportB,handler=f_importprof,action="db")

 helptext(tabimportB[1,1],"Imports sample profile(s) from a selected file into the software. \n\nThe column names must contain 'sample..', 'marker', 'allele..'. Optional: 'height..'")
 helptext(tabimportB[1,2],"Imports reference profile(s) from a selected file into the software. \n\nThe column names must contain 'sample..', 'marker', 'allele..'.")
 helptext(tabimportB[1,3],"Imports reference profile(s) from a selected file into the software. \n\nThe column names must contain 'sample..', 'marker', 'allele..'.")
 
 #Set View buttons
 tabimportB[2,1] = gWidgets2::gbutton(text="View evidence",container=tabimportB,handler=f_viewdata,action="mix")
 tabimportB[2,2] = gWidgets2::gbutton(text="View references",container=tabimportB,handler=f_viewdata,action="ref")
 tabimportB[2,3] = gWidgets2::gbutton(text="View database",container=tabimportB,handler=f_viewdata,action="db")
 gWidgets2::size(tabimportB[2,1]) <- gWidgets2::size(tabimportB[2,2]) <- gWidgets2::size(tabimportB[2,3]) <- 30

 helptext(tabimportB[2,1],"Shows the data in the selected evidence(s) both in terminal and in an epg-like plot (shown in the browser if plotly is installed). \n\nIf selected kit is recognized by the software and peak heights are imported, the sum of the peak heights per marker are fitted against fragment length assuming a gamma-model. \n\nThe p-value for an extreme marker is based on the fitted model when leaving out the marker.")
 helptext(tabimportB[2,2],"Shows the allele data for each selected reference(s) in terminal only. \n\nShows number of alleles for selected references matching against selected evidence(s).")
 helptext(tabimportB[2,3],"Shows the allele data for each reference(s) in selected database(s). \n\nShows number of alleles for each references in selected database(s) matching against selected evidence(s).")

 #Create tables for data selection
 tabimportB[3,1] = gWidgets2::gcheckboxgroup(items= getDataNames_type("mix"), container = tabimportB, use.table=TRUE)
 tabimportB[3,2] = gWidgets2::gcheckboxgroup(items= getDataNames_type("ref"), container = tabimportB, use.table=TRUE)
 tabimportB[3,3] = gWidgets2::gcheckboxgroup(items= getDataNames_type("db"), container = tabimportB, use.table=TRUE)
 gWidgets2::size(tabimportB[3,1]) <- gWidgets2::size(tabimportB[3,2]) <- gWidgets2::size(tabimportB[3,3]) <- c(20,150)
 #gWidgets2::add(tabimportB2a,tabimportB_mix)
 #gWidgets2::add(tabimportB2,tabimportB_ref,expand=TRUE,fill=TRUE)
 #gWidgets2::add(tabimportB2,tabimportB_db,expand=TRUE,fill=TRUE)

 #Set Delete buttons
 tabimportB[4,1] = gWidgets2::gbutton(text="Delete evidence",container=tabimportB,handler=f_deleteprof,action="mix")
 tabimportB[4,2] = gWidgets2::gbutton(text="Delete reference",container=tabimportB,handler=f_deleteprof,action="ref")
 tabimportB[4,3] = gWidgets2::gbutton(text="Delete database",container=tabimportB,handler=f_deleteprof,action="db")
 
 #helpfunction used to extract selected importdata-elements to further model-setup
 selectDataToModel <- function(h,....) {
   #All: popFreq must be imported!
   #GEN: No other requirement
   #EVID: must have both mixture and reference profiles
   #DB: Database search require both mixture and database, reference profiles is optional
   #DC: Deconvolution requires only mixtures. Reference profiles is optional
   popFreq <- get("popFreq",envir=mmTK)
   mixSel <- refSel <- dbSel <- numeric()
   if(length(tabimportB[3,1][])>0) mixSel <- gWidgets2::svalue(tabimportB[3,1])  #get selected mixtures
   if(length(tabimportB[3,2][])>0) refSel <- gWidgets2::svalue(tabimportB[3,2])  #get selected references
   if(length(tabimportB[3,3][])>0) dbSel <- gWidgets2::svalue(tabimportB[3,3])  #get selected databases

   if(is.null(popFreq)) {
    gWidgets2::gmessage("No frequencies was specified!\n Please import table with population frequencies.")
   } else if(h$action!="GEN" && length(mixSel)==0) {
    gWidgets2::gmessage("Please import and select mixture-profile!")
   } else if(h$action=="EVID" && length(refSel)==0) {
    gWidgets2::gmessage("Please import and select reference-profiles for weight of evidence!")
   } else if(h$action=="DB" && length(dbSel)==0) {
    gWidgets2::gmessage("Please import and select reference-database for database search!")
   } else {
    refreshTabModel(mixSel,refSel,dbSel,h$action) #refresh table with selected data
    gWidgets2::svalue(nb) <- 3 #change tab of notebook
   }
 }

 #Button-choices further:
 tabimportC[1,1] = gWidgets2::gbutton(text="Weight-of-Evidence",container=tabimportC,handler=selectDataToModel,action="EVID")
 helptext(tabimportC[1,1],"A module for calculating the Likelihood Ratio for selected sample(s) (treated as replicates). \n\nSelected reference(s) can later be conditioned on in the hypotheses.")

 tabimportC[1,2] = gWidgets2::gbutton(text="Deconvolution",container=tabimportC,handler=selectDataToModel,action="DC")
 helptext(tabimportC[1,2],"A module for ranking the most likely profiles of the unknown contributors for selected sample(s) (treated as replicates).\n\nThe user will first need to fit a gamma-model based on maximum likelihood estimation.")

 tabimportC[1,3] = gWidgets2::gbutton(text="Database search",container=tabimportC,handler=selectDataToModel,action="DB")
 helptext(tabimportC[1,3],"A module for calculating the Likelihood Ratio for selected sample(s) (treated as replicates) for each profile(s) in the selected database(s). \n\nSelected reference(s) can later be conditioned on in the hypotheses.")

 tabimportC[2,1] = gWidgets2::gbutton(text="Fit dropin data",container=tabimportC,handler=f_fitdropin)
 helptext(tabimportC[2,1],"A module for fitting the drop-in P.H. model based on data from an imported text file.")
 
 tabimportC[2,2] = gWidgets2::gbutton(text="Generate sample",container=tabimportC,handler=selectDataToModel,action="GEN")
 helptext(tabimportC[2,2],"A module for generating alleles based on selected population frequencies and corresponding peak heights based on the gamma-model.")

 tabimportC[2,3] = gWidgets2::gbutton(text="RESTART",container=tabimportC,handler=function(h,...) {
   gWidgets2::dispose(mainwin) #remove window!
   efm() #and open EFM again
 })
 helptext(tabimportC[2,3],"A button to restart EuroForMix")

#ADD EASY MODE
 if(get("optSetup",envir=mmTK)$easyMode) {
    gWidgets2::enabled(tabimportC[1,2]) <- FALSE  #deactivate deconvolution
    gWidgets2::enabled(tabimportC[1,3]) <- FALSE  #deactivate database search
    gWidgets2::enabled(tabimportC[2,1]) <- FALSE  #deactivate fit drop-in data
    gWidgets2::enabled(tabimportC[2,2]) <- FALSE  #deactivate generate sample
    gWidgets2::enabled(tabimportB[1,3]) <- FALSE #deactivete import database
    gWidgets2::enabled(tabimportB[2,3]) <- FALSE #deactivete database view
    gWidgets2::enabled(tabimportB[4,3]) <- FALSE #deactivete database delete
    
    gWidgets2::enabled(tabimportA[2,3]) <- FALSE #deactivete export frequencies
 }

####################################################################################################################
#######################################Tab 3: Model setup:##########################################################
#####################################################################################################################
 
 #Helpfunction to do model search (opens new windows)
 doModelSearch = function() {
   #prepare settings:
   opt <- get("optMLE",envir=mmTK) #options when optimizing (nDone,delta)
   checkPositive(opt$delta,"Variance parameter of randomizer")
   checkPosInteger(opt$nDone,"Number of random startpoints")
   
   set = get("setEVID",envir=mmTK) #get model setup 
   maxNOC = set$nC_hd #get number of contributors under hd  (this is max)   

   knownRefPOI = setdiff(set$knownref_hd, set$knownref_hp) #get index of POI (in Hd but not in Hp)
   if(length(knownRefPOI)!=1) {
     gWidgets2::gmessage("Only one POI can be selected. Please respecify hypotheses to proceed!")
     return(NULL) #return from function
   }
   POInames = names(set$refData[[1]])[knownRefPOI] #get name of POI to use
   condnames = names(set$refData[[1]])[which(set$condOrder_hd>0)] #get names 
   minNOC = length(condnames) + 1 #minimum number of contrs
   
   booltxt = c("YES","NO")
   defaultNOC = minNOC:maxNOC #set default range of NOC
   defaultBWS <- defaultFWS <- defaultDEG <- c(FALSE,TRUE) #booltxt[2] #default is no model 
   if(set$DEG) defaultDEG[1] = TRUE  #add possible to traverse Degradation model
   if(set$BWS) defaultBWS[1] =  TRUE  #add possible to traverse Backward stutter model
   if(set$FWS) defaultFWS[1] = TRUE #add possible to traverse Forward stutter model
   
   #user can change settings to traverse:
   NOCtxt = defaultNOC
   if(length(NOCtxt)>1) NOCtxt = paste0(min(NOCtxt),"-",max(NOCtxt))
   subwin = gWidgets2::gwindow("Select models to compare",visible=FALSE)  
   tabLay <- gWidgets2::glayout(spacing=spc,container=subwin)
   tabLay[1,1] =  gWidgets2::glabel( text="Select outcome for number of contributors (NOC):",container=tabLay)
   tabLay[1,2] =  gWidgets2::gedit( paste0(NOCtxt,collapse=",") ,container=tabLay)
   tabLay[2,1] =  gWidgets2::glabel( text="Select outcome for degradation:",container=tabLay)
   tabLay[2,2] =  gWidgets2::gcheckboxgroup(items=booltxt,checked =defaultDEG , horizontal=TRUE,container=tabLay )
   #tabLay[2,2] =  gWidgets2::gradio(items=booltxt,selected =defaultDEG , horizontal=TRUE,container=tabLay )
   #if(is.null(set$kit)) gWidgets2::enabled(tabLay[2,2]) = FALSE #can't change if kit (degrad) is not selected
   tabLay[3,1] =  gWidgets2::glabel( text="Select outcome for backward stutter:",container=tabLay)
   tabLay[3,2] =  gWidgets2::gcheckboxgroup(items=booltxt,checked = defaultBWS , horizontal=TRUE,container=tabLay )
   tabLay[4,1] =  gWidgets2::glabel( text="Select outcome for forward stutter:",container=tabLay)
   tabLay[4,2] =  gWidgets2::gcheckboxgroup(items=booltxt,checked = defaultFWS , horizontal=TRUE,container=tabLay )
   tabLay[5,2] =  gWidgets2::gbutton("Quit",container=tabLay, handler=function(h,...) gWidgets2::dispose(subwin) )
   tabLay[5,1] =  gWidgets2::gbutton("Evaluate",container=tabLay, handler=function(h,...) {
     #Set range of number of contributors (NOC):
     NOC = gWidgets2::svalue(tabLay[1,2]) #get number of outcome
     if(grepl("-",NOC)) { #if separator is included
       tmp = as.integer(unlist(strsplit(NOC,"-")))
       NOC = as.integer(tmp[1]:tmp[2])
     } else if(grepl(",",NOC)) { #if otehr separator is included
       NOC = as.integer(unlist(strsplit(NOC,",")))
     }
     if(length(NOC)==1) NOC = as.integer(NOC)
     if( length(NOC)==0 || !is.integer(NOC) || any(NOC<minNOC) || any(NOC>5) ) {
       gWidgets2::gmessage("The number of contributors was not correctly specified (or exceeded 5). Please respecify!")
       return(NULL) #return from function
     }
     #Set outcome of optional models:
     bool = c(TRUE,FALSE) #connects to booltxt = (YES/NO)
     modList = list() #boolean list for each opt models:
     for(i in 2:4) { #traverse models
       tmp = gWidgets2::svalue(tabLay[i,2]) 
       if(length(tmp)==0) {
         gWidgets2::gmessage("Missing model option. Please respecify!")
         return(NULL) #user must specify value!
       } 
       modList[[i-1]] = rev(bool[booltxt%in%tmp]) #find corresponding boolean (reverse to start with FALSE=simpler model)
     }
     modelDEG <- modList[[1]]
     modelBWS <- modList[[2]]
     modelFWS <- modList[[3]]
     
     txt1 = paste0("POI=",POInames)
     if(minNOC>1) txt1 = paste0(txt1,"\nConditionals=",paste0(condnames,collapse="/"))
     txt1 = paste0(txt1,"\nNumber of contributors={",paste(NOC,collapse=","),"}")
     txt1 = paste0(txt1,"\n\nModel combinations:") #
     txt1 = paste0(txt1,"\nDegrad: ",paste(booltxt[bool%in%modelDEG],collapse="/") ) 
     txt1 = paste0(txt1,"\nBW stutter: ",paste(booltxt[bool%in%modelBWS],collapse="/") ) 
     txt1 = paste0(txt1,"\nFW stutter: ",paste(booltxt[bool%in%modelFWS],collapse="/") ) 
     
     gWidgets2::dispose(subwin) #dispose window after extracting values (checking, and continue)
     evalBool = gWidgets2::gconfirm(paste0("Following setup to be evaluated:\n",txt1,"\n\nDo you want to continue?"),title="Searching for optimal model",icon="info")
     if(!evalBool) return(NULL) #return if not evaluating    
     
     #FIND OPTIMAL MODEL (use settings under Hd):
     searchList <- contLikSearch(NOC,modelDEG,modelBWS,modelFWS,set$samples,set$popFreq,set$refData,set$condOrder_hd, knownRefPOI,set$prC, opt$nDone,set$threshT,set$fst,set$lambda,opt$delta, set$kit,TRUE,opt$difftol,set$knownRel,set$ibd, set$minFreq, set$normalize, set$priorBWS, set$priorFWS, opt$seed, opt$maxThreads, opt$alpha, opt$steptol, set$adjQbp)
     AIC = searchList$outtable[,2] #obtain crietion
     optimInd = which.max(AIC)[1] #get index of optimal model. Use simpler model if "Tie"
     
     #CREATE AND STORE RESULT TABLE
     searchTable = cbind(searchList$modoutcome, searchList$outtable) #get search table
     
     #Show as table in GUI
     tabSearch_win = gWidgets2::gwindow("Model comparison results",visible = FALSE,width=550)
     tabSearch_GUI = gWidgets2::gtable(searchTable ,container=tabSearch_win) #show table to user
     gWidgets2::size(tabSearch_GUI) = list(column.widths=c(35,rep(55,ncol(searchTable)-1)))
     gWidgets2::svalue(tabSearch_GUI,index=1) = optimInd #optimInd #highlight best model
     gWidgets2::visible(tabSearch_win) = TRUE
     
     #store optimal model to SET
     set$mlefit_hp = searchList$hpfitList[[optimInd]]
     set$mlefit_hd = searchList$hdfitList[[optimInd]]
     set$SEARCH = searchTable #Store search table here
     
     #NEED TO UPDATE the SET object (NOC for each hyp, and 
     storeLRvalues(set) #store LR values (needed before refreshing model etc)
     assign("setEVID",set,envir=mmTK) #store setup for EVID
     
     #Important to update results:
     refreshTabMLE("START") #Show final model in MLE tab
     gWidgets2::svalue(nb) <- 4 #go to mle-fit window when finished
   }) #end evaluate button
   gWidgets2::visible(subwin) = TRUE
 } #end function
 
  #Helpfunction to do data selection (opens new windows)
  selectDataFromUser = function(SELDATA) {
    #SELDATA = get("SELDATA",envir=mmTK) 
    mixData = SELDATA$mixData
    refData = SELDATA$refData
    locs = SELDATA$locs
    mixSel = names(mixData)
    refSel = names(refData)
    nSamples = length(mixSel)
    sellocs = NULL #Loci to be selected
    
    #winsize = c(500,1000)
    #CREATING A NEW WINDOW
    selwin = gWidgets2::gwindow("Data selection",visible=FALSE)  
    panel  = gWidgets2::ggroup(horizontal = FALSE,container = selwin, use.scrollwindow = TRUE,expand=TRUE,fill=TRUE)
    tab  = gWidgets2::glayout(spacing=3,container = panel)#,expand=TRUE,fill=TRUE)
    
    #Show loci names
    for(loc in locs) { #insert locus names from popFreq
       tab[1+which(loc==locs),1] <- loc  #insert loc-name
    }
  
    #SELECTION OF EVIDS
    for(msel in mixSel) { #for each selected mixture
      tab[1,1 + which(msel==mixSel)] <- gWidgets2::glabel(text=msel,container=tab) #get selected mixturenames
      for(loc in locs) { #for each locus
       exist <- !is.null(mixData[[msel]][[loc]]$adata) ##&& !is.null(mixData[[msel]][[loc]]$hdata) #check if exist alleles!
       tab[1+which(loc==locs),1 + which(msel==mixSel)]  <- gWidgets2::gcheckbox(text="",container=tab,checked=exist)
       if(!exist) gWidgets2::enabled(tab[1+which(loc==locs),1 + which(msel==mixSel)]) <- FALSE #deactivate non-existing locus
      }
    }  
    
    #SELECTION OF REFS
    for(rsel in refSel) { #for each selected reference
      tab[1,1 + nSamples + which(rsel==refSel)] <- gWidgets2::glabel(text=rsel,container=tab) #name of reference
      for(loc in locs) { #for each locus
       adata = unlist(refData[[rsel]][[loc]])#$adata
       exist <- !is.null(adata) && length(adata)==2 #check if allele exists (assume duploid!)
       tab[1+which(loc==locs),1 + nSamples + which(rsel==refSel)]  <- gWidgets2::gcheckbox(text="",container=tab,checked=exist)
       if(!exist) gWidgets2::enabled(tab[1+which(loc==locs),1 + nSamples + which(rsel==refSel)]) <- FALSE #deactivate non-existing locus
      }
    }
    
    #Include button which saves user-changes:
    tab[1,1] <- gWidgets2::gbutton("Confirm",container=tab, handler=function(h,...) {
      
      bool = gWidgets2::gconfirm("Do you want to apply the data selection changes?",icon="question")
      if(!bool) return(NULL) 
      
      for(loc in locs) { #for each locus
        isOK <- TRUE
        for(msel in mixSel) { #Locus is not evaluated if any of the evid profiles are turned off
          isOK <- isOK && gWidgets2::svalue(tab[1+which(loc==locs),1 + which(msel==mixSel)])  #check if locus checked for stains
        }
        for(rsel in refSel) {
          bool <- gWidgets2::svalue(tab[1+which(loc==locs),1 + nSamples + which(rsel==refSel)]) #check if locus checked for references
          if(!bool || is.null(refData[[rsel]][[loc]])) refData[[rsel]][[loc]]$adata = numeric() #INSERT numeric() if missing or not checked
        }
        if(isOK) sellocs <- c(sellocs,loc) #locus can be evaluated
      } #end for each locus
    
      if(length(sellocs)==0) { #don't do anything if no loci will be evaluated
        gWidgets2::gmessage("No loci are evaluated! Be sure that all selected data have valid data in their loci.",title="No loci found!",icon="error")
      }
      assign("SELDATA",list(mixData=mixData,refData=refData,locs=sellocs),envir=mmTK)   #STORE SELECTED DATA
      gWidgets2::dispose(selwin) #close window
    })
    gWidgets2::visible(selwin) <- TRUE
  } #end data selection function
 

#  mixSel=names( getData("mix"));refSel=names(getData("ref") );dbSel=numeric();type="EVID"
  #THIS IS A MAIN FUNCTION TO SHOW Model SPECIFICATION
  refreshTabModel = function(mixSel,refSel,dbSel,type) { 
    #type={"GEN","EVID","DB","DC"}
    gWidgets2::visible(mainwin) <- FALSE # dispose(tabmodel) 
    
    #a nested main "helpfunction" which takes GUI settings and stores them in "set'type'"
    storeSettings = function(lrtype="PLOT") {
      #lrtype="CONT","QUAL","PLOT"
      optL = get("optSetup",envir=mmTK) #get setup list
      SELDATA = get("SELDATA",envir=mmTK) #obtain selected data (potentially modified with selectData function) 
      mixData = SELDATA$mixData
      refData = SELDATA$refData
      locs = SELDATA$locs #some loci may have been de-selected (if any evid profiles de-selected) 
      mixSel = names(mixData)
      refSel = names(refData)
      
      #HERE: GET USER SPECIFIED DATA SELECTION     
      popFreq <- popFreq[locs] #consider only relevant loci in popFreq
      print(c("Locs to be evaluated:",paste0(locs,collapse=",")))
      
      #Check if samples have peak heights if cont. model is used:
      if(lrtype=="CONT") {
        hdatas <- sapply( mixData, function(x) sum(sapply(x,function(y) length(y$hdata)) ) ) #number of alleles for each samples
        if(any(hdatas==0)) {
          gWidgets2::gmessage(paste0("The sample(s) ", paste0(mixSel[hdatas==0],collapse=",")," did not have peak heights! \nEvaluate with qualitative LR model"),title="Wrong data input!",icon="error")
          stop("No evaluation done.")
        }
      }
      
      #Check for duplicated alleles in loci
      for(loc in locs) {
        tmp <- popFreq[[loc]] #get allele-names in popFreq
        newA <- numeric()
        for(ss in mixSel) { #for each stain sample
          adata <- mixData[[ss]][[loc]]$adata
          if(length(adata)!=length(unique(adata))) {
            gWidgets2::gmessage(paste0("The sample ", ss," has duplicated alleles in locus ",loc,"\nPlease remove the duplicate from the text-file!"),title="Wrong data input!",icon="error")
            stop("No evaluation done.")
          }
        }
      } #end for each loci
      
      #READ FROM GUI TO DECIDE MODELS:
      DEG <- BWS <- FWS <- FALSE
      if(type!="GEN") { #if not generating data
        DEG <- gWidgets2::svalue(tabmodelMODopt[1,2])=="YES" #get boolean of degradation
        BWS <- gWidgets2::svalue(tabmodelMODopt[2,2])=="YES" #get boolean of degradation
        FWS <- gWidgets2::svalue(tabmodelMODopt[3,2])=="YES" #get boolean of degradation
        if(!BWS && FWS) {
          gWidgets2::gmessage("Can't turn on Forward Stutter model without turning on Backward stutter model.",title="Wrong model choice",icon="error")
          stop("No evaluation done.")
        }
      }

      #READ FROM SETTINGS (default values):
      threshT <- optL$thresh0 #get specified threshold
      prC <-   optL$pC0 #prC is probability of dropin
      lambda <-   optL$lam0 #lambda is hyperparameter to dropin-peak height model
      fst <- optL$fst0 #theta-correction
      priorBWS = eval( parse(text= paste("function(x)",optL$pXi)) , envir=mmTK ) #get prior from settings. 
      priorFWS = eval( parse(text= paste("function(x)",optL$pXiFW)) , envir=mmTK ) #get prior from settings. 
     
      #Receive marker specific settings (Replace the constant settings:
      optLmarker = get("optMarkerSetup",envir=mmTK) #this is marker specific 
      if(!is.null(optLmarker)) { #if specified
        for(c in 1:length(optLmarker)) { #for each type USE IBUILT FUNCTION setVecRightOrder ??
          tmp = setVecRightOrder(optLmarker[[c]],locs) #get right order of threshv,pCv,lamv,fstv
          if(c==1) threshT = tmp #get AT threshold per-marker values
          if(c==2) prC = tmp  #get dropin prob per-marker values
          if(c==3) lambda = tmp #get lambda per-marker values
          if(c==4) fst = tmp #get fst per-marker values
        }
      }
      
      checkPrior <- function(txt="") { #throw message 
        gWidgets2::gmessage(txt,title="Wrong data input!",icon="error")
        stop(txt)
      }
      tryCatch({priorBWS(0.1)}, error = function(e) checkPrior("BW stutter prior was wrongly specified.") ) 
      tryCatch({priorFWS(0.1)}, error = function(e) checkPrior("FW stutter prior was wrongly specified.") ) 
      if(!is.numeric(priorBWS(0.1))) checkPrior()
      
      #assign model of degradation and stutter models:
      kit <- get("selPopKitName",envir=mmTK)[1] #get selected kit
      if(is.na(kit)) kit <- NULL
      
      #CHECK VARIABLES:
      sapply(threshT,function(x) checkPositive(x,"Threshold"))
      sapply(prC,function(x) checkProb(x,"Allele drop-in probability"))
      if(any(prC>0) && lrtype=="CONT") sapply(lambda, function(x) checkPositive(x,"Drop-in peak height hyperparameter",strict=TRUE) )
      sapply(fst, function(x) checkProb(x,"fst-correction"))
      
      #prepare data for function in euroformix! Data-object must have correct format!
      #go through selected data from GUI:
      samples <- list()
      refData2 <- NULL
      if(nRefs>0) refData2 <- list()
      for(loc in locs) { #for each selected locus in popFreq
        for(msel in mixSel) { #for each mixture samples
          subD <- mixData[[msel]][[loc]] #get selected data
          if(!is.null(subD$hdata)) { #if PH data given
            threshT0 = getMarkerVal(threshT,loc) #Find correct threshold to use
            keep <- subD$hdata>=threshT0 #allele to keep (above threshold)
            subD$hdata <- subD$hdata[keep]
            subD$adata <- subD$adata[keep]
          }
          samples[[msel]][[loc]] <- subD #insert samples
        }
        if(nRefs>0) {
          refData2[[loc]] <- list()
          for(rsel in refSel) {
            av = unlist(refData[[rsel]][[loc]])
            if(is.null(av)) av = character() #important to avoid NULL
            refData2[[loc]][[rsel]] <-  av #$adata #insert references: Note the changing format!!
          }
        }
      } #end for each locus
      
      #get specified hypotheses  (contributors)
      condOrder_hp <- condOrder_hd <- rep(0,nRefs)
      if(type=="DC") condOrder_hp <- NULL
      for(rsel in refSel) { #for each reference under hp and hd
        if(!type%in%c("DC","GEN")) {
          valhp <- as.integer(gWidgets2::svalue(tabmodelMODhp[which(rsel==refSel),1])) 
          condOrder_hp[which(rsel==refSel)] <- valhp +  valhp*max(condOrder_hp)
        }
        valhd <- as.integer(gWidgets2::svalue(tabmodelMODhd[which(rsel==refSel),1])) 
        condOrder_hd[which(rsel==refSel)] <- valhd + valhd*max(condOrder_hd)
      }
      
      #get specified hypotheses (typed non-contributors but): 
      #UPDATE IN v3.4.0: INCLUDE REFS NOT PUT FORWARD IN ANY HYPOTHESIS
      knownref_hp <- which(condOrder_hp==0) #those references not condition on under hp
      knownref_hd <- which(condOrder_hd==0) #those references not condition on under hd
      if(length(knownref_hp)==0) knownref_hp = NULL
      if(length(knownref_hd)==0) knownref_hd = NULL
      
      #number of contributors in model:
      nC_hp <- NULL
      if(!type%in%c("DC","GEN")) {
        nC_hp <-  as.integer(gWidgets2::svalue(tabmodelMODhp[nRefs+1,2])) + sum(condOrder_hp>0)
        checkPosInteger(nC_hp + sum(type=="DB"),"Number of contributors under Hp")
      }
      nC_hd <-  as.integer(gWidgets2::svalue(tabmodelMODhd[nRefs+1,2])) + sum(condOrder_hd>0)
      checkPosInteger(nC_hd,"Number of contributors under Hd")
      
      #get model parameters: This is also done for SNPs!
      Qdes <- TRUE #always true in EFM calculation (drastic speedup)
      if(lrtype=="GEN") Qdes <- FALSE #Q-designation turned off when generating data
      optfreq <- get("optFreq",envir=mmTK) #get frequency option
      minFreq <- getminFreq()
      normalize = as.logical(optfreq$normalize) #get bool of whether to normalize
      ret <- Qassignate(samples, popFreq, refData2,doQ=Qdes,incS=FALSE,incR=FALSE,minF=minFreq,normalize=normalize) #call function in euroformix
      popFreqQ <- ret$popFreq
      refDataQ <- ret$refData
      
      ############################################
      #################RELATEDNESS################
      ############################################
      #get relationship type of last unknown under Hd
      rel_type <- gWidgets2::svalue( tabmodelMODhd[nRefs+3,1]) #get selected relationship (name)
      rel_ibd <- relatednessIBD[relname==rel_type,] #get selected ibd
      names(rel_ibd) <- rel_type
      rel_refName <- gWidgets2::svalue( tabmodelMODhd[nRefs+5,1]) #get name of related individual
      if(rel_type==relname[1] || rel_refName==refSelItems[1] || length(refSel)==0) { #refSel has lenght zero if not selected
        rel_ref = NULL
      } else {
        rel_ref <-  which(refSel==rel_refName) #get index
        names(rel_ref) <- rel_refName
      }
      #knownref_hd
      #get input to list: note: "fit_hp" and "fit_hd" are list-object from fitted model
      #FINALLY STORE IN SETTINGS LIST:
      set <-list(samples=samples,refData=refData2,popFreq=popFreq,refDataQ=refDataQ,popFreqQ=popFreqQ, minFreq=minFreq,normalize=normalize,    #DATA
          nC_hp=nC_hp,nC_hd=nC_hd,condOrder_hp=condOrder_hp,condOrder_hd=condOrder_hd,knownref_hp=knownref_hp,knownref_hd=knownref_hd,ibd=rel_ibd,knownRel=rel_ref, #propositions
          prC=prC,threshT=threshT,fst=fst,lambda=lambda,priorBWS=priorBWS,priorFWS=priorFWS,kit=kit,DEG=DEG,BWS=BWS,FWS=FWS, adjQbp= optL$adjQbp) #Model
      if(type=="DB") set$dbSel <- dbSel #store name of selected databases
      assign(paste0("set",type),set,envir=mmTK) #store setup for relevant type
    } #end store settings from GUI to environment
    ####################################################
    
    #PREPARING DATA (finds relevant locs to use: etc)
    popFreq <- get("popFreq",envir=mmTK)
    locs <- names(popFreq)  #get names of loci for imported population frequencies. used as default in list
    kitname <- get("selPopKitName",envir=mmTK)  #get selected kit to use
    hasKit = FALSE
    if(!is.null(kitname) && !is.na(kitname[1])) {
      kitinfo = getKit(kitname[1],"Marker")
      if(!is.na(kitinfo[1])) {
        hasKit = TRUE #has a recognized kit
        locs <- intersect(locs,toupper(kitinfo)) #use only overlapping allele
      }
    }
    nRefs = length(refSel) #number of ref-profiles
    mixData=getData("mix",mixSel)
    #create temporary variable storage for data selection (own function which syncs with data object):
    SELDATA = list(mixData=mixData,refData=getData("ref",refSel),locs=locs) #obtain selection data object
    assign("SELDATA",SELDATA,envir=mmTK)  #also store here first time
    
    #CHECK WHAT TYPE IT IS: INDICATE WHETHER THE DATA FORMAT IS SOMETHING ELSE THAN CE/EPG: MPS OR LUS. And also maybe SNP:
    isMPS = FALSE
    isLUS = all(unlist(sapply(mixData,function(x)  sapply(x,function(y) all(grepl(LUSsymbol,y$adata),na.rm=TRUE)) )))  #ADDED: check if alleles are given on LUS format 
    if(!isLUS) suppressWarnings( { isMPS = all(unlist(sapply(mixData,function(x)  sapply(x,function(y) all(is.na( as.numeric(y$adata) )))))) }) #ADDED: check if alleles are given on string-format (i.e MPS)grepl(y$adata),na.rm=TRUE)) )))    
    isSNP <- all(sapply(popFreq,length)==2) #check if SNP data: I.e. 2 allele outcome for all markers - stutter/deg models not possible
    
    #############################################################################
    #GUI PROCEEDS: Orgainizing model options and data for the user to select:
    tabmodeltmp <- gWidgets2::glayout(spacing=spc,container= tabmodel[1,1] <- gWidgets2::ggroup(container=tabmodel) ) 
    
    #Left and right table layout
    tabmodelMOD = gWidgets2::glayout(spacing=5,container=(tabmodeltmp[1,1] <-gWidgets2::gframe("Model specification",container=tabmodeltmp))) 
    tabmodelSELCALC = gWidgets2::glayout(spacing=10,container=(tabmodeltmp[1,2] <-gWidgets2::gframe(spacing=10,container=tabmodeltmp)))  
    
    #Upper and lower table layout (for right table layout)
    tabmodelSEL = gWidgets2::glayout(spacing=0,container=(tabmodelSELCALC[1,1] <-gWidgets2::gframe("Data",container=tabmodelSELCALC)))  
    tabmodelCALC = gWidgets2::glayout(spacing=5,container=(tabmodelSELCALC[2,1] <-gWidgets2::gframe("Calculations",container=tabmodelSELCALC)))  
    
    #Hypothesis selection: subframe of A
    txt <- "Contributor(s) under Hp:"
    if(type=="DB") txt <- paste0(txt, "\n(DB-reference already included)")
    if(type%in%c("DB","EVID")) tabmodelMODhp = gWidgets2::glayout(spacing=0,container=(tabmodelMOD[2,1] <-gWidgets2::gframe(txt,container=tabmodelMOD))) 
    tabmodelMODhd = gWidgets2::glayout(spacing=0,container=(tabmodelMOD[3,1] <-gWidgets2::gframe("Contributor(s) under Hd:",container=tabmodelMOD)))
    if(type!="GEN")  tabmodelMODopt = gWidgets2::glayout(spacing=0,container=(tabmodelMOD[4,1] <-gWidgets2::gframe("Model options",container=tabmodelMOD))) 
    
    #specify references under hypotheses
    for(rsel in refSel) {
    if(type%in%c("DB","EVID")) tabmodelMODhp[which(rsel==refSel),1]  <- gWidgets2::gcheckbox(text=rsel,container=tabmodelMODhp,checked=TRUE) #Check as default under Hp
    tabmodelMODhd[which(rsel==refSel),1]  <- gWidgets2::gcheckbox(text=rsel,container=tabmodelMODhd,checked=!(type=="EVID")) #references under Hd (not checked if evidnece)
    }
    
    Krange <- 0:4 #default Contr range
    #specify number of unknowns
    if(!type%in%c("DC","GEN")) {
    tabmodelMODhp[nRefs+1,1] <- gWidgets2::glabel(text="#unknowns (Hp): ",container=tabmodelMODhp)
    tabmodelMODhp[nRefs+1,2] <- gWidgets2::gcombobox(items=Krange,selected=2,editable=TRUE,container=tabmodelMODhp)
    gWidgets2::size(tabmodelMODhp[nRefs+1,2]) = 4 #set with
    #    tabmodelMODhp[nRefs+1,2] <- gWidgets2::gedit(text="1",container=tabmodelMODhp,width=4)
    }
    tabmodelMODhd[nRefs+1,1] <- gWidgets2::glabel(text="#unknowns (Hd): ",container=tabmodelMODhd)
    #tabmodelMODhd[nRefs+1,2] <- gWidgets2::gedit(text="2",container=tabmodelMODhd,width=4)
    tabmodelMODhd[nRefs+1,2] <- gWidgets2::gcombobox(items=Krange,selected=3,editable=TRUE,container=tabmodelMODhd)
    gWidgets2::size(tabmodelMODhd[nRefs+1,2]) = 4 #set with
    
    #BLOCK added in version 2.0.0: Relatedness
    tabmodelMODhd[nRefs+2,1] <-  gWidgets2::glabel(text="\nLast unknown is",container=tabmodelMODhd)
    relatednessIBD = rbind( c(1,0,0), t(replicate(2,c(0,1,0))) , c(1/4,1/2,1/4),  t(replicate(5,c(1/2,1/2,0))) ,c(3/4,1/4,0), c(0,0,1) )  #Defined at https://strbase.nist.gov/pub_pres/Lewis-Towson-kinship-Apr2010.pdf
    relname <- rownames(relatednessIBD) <- c("Unrelated" , "Parent","Child", "Sibling" , "Uncle","Nephew","Grandparent","Grandchild","Half-sibling" , "Cousin", "Twin (ident.)")
    tabmodelMODhd[nRefs+3,1] <-  gWidgets2::gcombobox(items=rownames(relatednessIBD),container=tabmodelMODhd,editable=FALSE) 
    tabmodelMODhd[nRefs+4,1] <-  gWidgets2::glabel(text="to",container=tabmodelMODhd)
    refSelItems <- c("",refSel) #selection list of references
    tabmodelMODhd[nRefs+5,1] <-  gWidgets2::gcombobox(items=refSelItems,container=tabmodelMODhd,editable=FALSE)
    gWidgets2::size(tabmodelMODhd[nRefs+3,1]) <- gWidgets2::size(tabmodelMODhd[nRefs+5,1]) <- 10
    #gWidgets2::enabled(  tabmodelMODhd[nRefs+3,1] ) <- FALSE #not implemented
    
    #Model parameters: 
    if(type!="GEN") {
      tabmodelMODopt[1,1] <- gWidgets2::glabel(text="Degradation:",container=tabmodelMODopt)
      tabmodelMODopt[1,2] <- gWidgets2::gradio(items=c("YES","NO"), selected = ifelse(hasKit,1,2), horizontal = TRUE,container=tabmodelMODopt)
      tabmodelMODopt[2,1] <- gWidgets2::glabel(text="BW Stutter:",container=tabmodelMODopt)
      tabmodelMODopt[2,2] <- gWidgets2::gradio(items=c("YES","NO"), selected = 2, horizontal = TRUE,container=tabmodelMODopt)
      tabmodelMODopt[3,1] <- gWidgets2::glabel(text="FW Stutter:",container=tabmodelMODopt) #added version 3.0.0
      tabmodelMODopt[3,2] <- gWidgets2::gradio(items=c("YES","NO"), selected = 2, horizontal = TRUE,container=tabmodelMODopt)
      if(!hasKit) gWidgets2::enabled( tabmodelMODopt[1,2] ) <- FALSE  #is kit defined? If not then don't consider degradation model
      
      #added v1.9.3:
      if(isSNP) { #deactivate model options if SNPs
       gWidgets2::enabled( tabmodelMODopt[1,2] ) <- FALSE
       gWidgets2::enabled( tabmodelMODopt[2,2] ) <- FALSE
       gWidgets2::enabled( tabmodelMODopt[3,2] ) <- FALSE
      } else {
       #Stutters accepted for numerical vairants (STR-RU/LUS)   
       if(isMPS) {
         gWidgets2::enabled( tabmodelMODopt[2,2] ) <- FALSE  #Deactivate stutter if alleles are strings (BW)
         gWidgets2::enabled( tabmodelMODopt[3,2] ) <- FALSE  #Deactivate stutter if alleles are strings (FW)
       }
      } #end if not SNP 
    } #end if not GEN
    
    #View evaluating evidence/databases 
    #Plot EPG of selected data (calls storeSettings first and then plotEPG by importing selected data)
    if(type!="GEN") {
      #Data selection: CREATES A SEPARATE WINDOW
      tabmodelSELevid = gWidgets2::glayout(spacing=0,container=(tabmodelSEL[1,1] <-gWidgets2::gframe("Evidence(s)",container=tabmodelSEL))) 
      for(msel in mixSel) {
        tabmodelSELevid[which(msel==mixSel),1] <- gWidgets2::gcheckbox(text=msel,container=tabmodelSELevid,checked=TRUE)
        gWidgets2::enabled(tabmodelSELevid[which(msel==mixSel),1]) <- FALSE #cant select here
      }
      tabmodelSEL[2,1] = gWidgets2::gbutton(text="Select data",container=tabmodelSEL,handler= function(h,...) { selectDataFromUser(SELDATA) })
      tabmodelSEL[3,1] = gWidgets2::gbutton(text="Show selected",container=tabmodelSEL,handler= function(h,...) { 
        storeSettings("PLOT") #store settings
        #loads each selections and plot epg:
        set = get(paste0("set",type),set,envir=mmTK) #store setup for relevant type
        print("Assumed population frequencies:")
        print(unlist(set$popFreqQ))
        print("Considered references:")
        .printTable(t(sapply(set$refDataQ,function(x) {  sapply(x,function(y) { paste0(y,collapse="/") } ) })))
        print("Considered Evidence samples:")
        for(sn  in names(set$samples)) {
         print("------------------------------------")
         print(paste("Samplename: ",sn,sep=""))
         printEvid(set$samples[[sn]])
        }
      })
      if(type=="DB") { #add database-names if included:
      tabmodelSELdb = gWidgets2::glayout(spacing=0,container=(tabmodelSEL[3,1] <-gWidgets2::gframe("Database(s) to search",container=tabmodelSEL))) 
      for(dsel in dbSel) tabmodelSELdb[which(dsel==dbSel),1] =  gWidgets2::glabel(text=dsel,container=tabmodelSELdb)
      }
    }

    #Calculation button:  
    if(type=="GEN") {
      tabmodelCALC[1,1] = gWidgets2::gbutton(text="Generate sample",container=tabmodelCALC,handler= function(h,...) { 
       storeSettings("CONT") #store settings
       refreshTabGEN() #generate a dataset based on stored settings
       gWidgets2::svalue(nb) <- 1 #go to data generation-window
      })
    } else {
      tabmodelCALC[1,1] = gWidgets2::gbutton(text="Quantitative LR\n(Maximum Likelihood based)",container=tabmodelCALC,handler=
      	function(h,...) {
            storeSettings("CONT") #store settings
            refreshTabMLE(type) #refresh MLE fit tab (i.e. it fits the specified model)
            gWidgets2::svalue(nb) <- 4 #go to mle-fit window (for all cases) when finished
        }) #end cont. calculation button
        
      if(type=="EVID") { # Insert own button for model searching for LR
        tabmodelCALC[2,1] = gWidgets2::gbutton(text="Optimal quantitative LR\n(automatic model search)",container=tabmodelCALC,handler=function(h,...) {
          storeSettings("CONT") #store settings FIRST
          doModelSearch() #excecure model search         
        })
      } #end if EVID

      if(type%in%c("EVID","DB")) {
        txt = "Qualitative LR\n(semi-continuous)"
        tabmodelCALC[3,1] = gWidgets2::gbutton(text="Qualitative LR\n(semi-continuous)",container=tabmodelCALC,handler=function(h,...) {
          storeSettings("QUAL") #store model-settings (use other input)
          if(type=="DB") {
            doDBsearch("QUAL")
          } else {
            refreshTabLRMIX() #refresh LRmix calculation tab (i.e. it fits the specified model)
            gWidgets2::svalue(nb) <- 7 #go to LRmix tab
          }
        }) #end cont. calculation button
  
        if(get("optSetup",envir=mmTK)$easyMode) { #CHECK IF EASY MODE:
          gWidgets2::enabled(tabmodelCALC[3,1]) <- FALSE  #deactivate qualitative model
        }
      } #end if evid or db
   } #end if not gen
    gWidgets2::visible(mainwin) <- TRUE
    gWidgets2::focus(mainwin) <- TRUE
  } #end refresh setup tab-frame


##########################################################################################
############################Tab 4: MLE estimation:########################################
##########################################################################################
#WE DO MLE-FITTING HERE, and also DECONVOLUTION-function AND DATABASE SEARCHING is implemented here (saves memory usage?)
  
  #Helpfunction to carry out database searching
  doDBsearch <- function(ITYPE="QUAN") {  #take a given inference-type {"MLE","INT","QUAL"} #qual means that only qualitative LR is done
    DBtab = efm_DBsearch(mmTK,ITYPE)
    assign("resDB",DBtab ,envir=mmTK) #assign deconvolved result to environment
    if(ITYPE=="QUAN") refreshTabDB(1) #cont LR is order to sort with
    if(ITYPE=="QUAL") refreshTabDB(2) #qual LR is order to sort with
    gWidgets2::svalue(nb) <- 6 #go to database search results window when finished     
    gWidgets2::focus(mainwin) <- TRUE
  } #end doDBsearch
  
# helpfunction to perform Bayesian integration (Bayes Factor)
  doINT <- function(type="EVID") { #Used either with EVID (or DB
    dig=4 #number of digits in output
    set <- get(paste0("set",type),envir=mmTK)
    if(is.null(set)) return()
    
    bool = gWidgets2::gconfirm("Are you sure you want to proceed with the numerical integration of the Full Bayesian LR (Bayes Factor)?\nIt may take a while!",title="Quantitative LR (Bayesian based)")#,icon="question")
    if(!bool) return()   
    
    optint <- get("optINT",envir=mmTK) #options when integrating (reltol and boundaries)
    if(is.null(optint$dev)) optint$dev=3 #use default
    print(paste0("Calculating integrals with relative error ",optint$reltol))
    print("This may take a while...")
    
    #Use wrapper function to calculate:
    LRint = calcLRint(set$mlefit_hp,set$mlefit_hd, reltol = optint$reltol, maxEval=optint$maxeval, dev=optint$dev,verbose = TRUE)
    
    #Print a GUI message:
    txt <- paste0("The Bayes Factor (Bayesian based LR)\nwas calculated as \nlog10LR=",format(LRint$log10LR,digits=dig)," [",format(LRint$log10LRerror[1],digits=dig)," , ",format(LRint$log10LRerror[2],digits=dig),"]")
    #print(txt)
    
    bool = gWidgets2::gconfirm(paste0(txt,"\nDo you want to store the results?"),title="Quantitative LR (Bayesian based)")#,icon="question")
    if(bool) {
      resEVID = get("resEVID",envir=mmTK) #assign integrated LR to session, stored in report
      resEVID$INT = list(LR=LRint$LR, LRerror = LRint$LRerror, reltol=optint$reltol, maxEval=optint$maxeval, dev=optint$dev) #insert Bayesian calculated LR to evid results 
      assign("resEVID",resEVID,envir=mmTK) #assign integrated LR to session, stored in report
    }
  } #end Integration

  #helpfunction ran when call deconvolution
  doDC <- function(mleobj) {
     dcopt <- get("optDC",envir=mmTK) #options when Deconvolution
     dcobj <- deconvolve(mlefit=mleobj,alpha=dcopt$alphaprob,maxlist=dcopt$maxlist) 
     DCtable1<-addRownameTable(dcobj$table2)
     colnames(DCtable1)[1] <- "Locus"
     DCtable2<-dcobj$table1
     DCtable3<-dcobj$table3
     DCtable4<-dcobj$table4
     assign("resDC",list(DCtable1,DCtable2,DCtable3,DCtable4),envir=mmTK) #assign deconvolved result to environment
     refreshTabDC() #update table with deconvolved results
     gWidgets2::svalue(nb) <- 5 #go to deconvolution results window (for all cases) when finished     
     gWidgets2::focus(mainwin) <- TRUE
   }

  #helpfunction ran when call MCMC
  doMCMC <- function(mleobj,showValid=TRUE,seed=1,delta=NULL,niter=NULL) { 
     optlist <- get("optMCMC",envir=mmTK)  #options for MCMC 
     #optint <- get("optINT",envir=mmTK) #get boundaries
     if(any(is.na(mleobj$fit$thetaSigma))) return();
     if(is.null(delta)) delta = optlist$delta #obtain delta from toolbar
     if(is.null(niter)) niter = optlist$niter #obtain number of iterations from toolbar
     
     print(paste0("Sampling ",niter," samples with variation ",delta,":"))
     print("Note: You can change default number of iterations in toolbar menu.")
#     mcmcfit <- contLikMCMC(mleobj,uppermu=optint$maxmu,uppersigma=optint$maxsigma,upperxi=optint$maxxi,optlist$niter,optlist$delta)
#     mcmcfit <- contLikMCMC(mleobj,optlist$niter,optlist$delta)
     mcmcfit <- contLikMCMC(mleobj,niter,delta,seed=seed) 
     print(paste0("Sampling acceptance rate=",round(mcmcfit$accrat,2),". This should be around 0.25 [0.15-0.35]"))
     print(paste0("Estimation of the marginalized likelihood (log)=",mcmcfit$logmargL))
     if(showValid) validMCMC(mcmcfit,acf=FALSE) #don't plot acf
     return(mcmcfit)
  }

  #Simulating LR over the parameter space
  #mlehp=set$mlefit_hp;mlehd=set$mlefit_hd
  simLR = function(mlehp,mlehd) {
     if(any(is.na(mlehp$fit$phiSigma)) || any(is.na(mlehd$fit$phiSigma)) ) return();
     opt <- get("optMCMC",envir=mmTK)  #options for MCMC 
     rec_niters = .getMinESS( max(length(mlehp$fit$phihat),length(mlehd$fit$phihat))) #obtain number of recommended samples
     bool = gWidgets2::gconfirm(paste0("Are you sure you want to proceed with\nthe LR sensitivity analysis?\n\nThe number of iterations is set to ",opt$niter,"\nRecommended minimum number of samples is ",rec_niters,"\n\nYou can change this number under the MCMC settings."))#,icon="question")
     if(!bool) return()
     
     #Use new function directly
     res <- get("resEVID",envir=mmTK) #get all setup-object 
     delta = opt$delta
     mcmcObjList = NULL
     if(!is.null(res$MCMC)) { #if already run
       delta = res$MCMC$deltaTuned #obtain tuned delta from ealier run
       mcmcObjList = res$MCMC$mcmcObjList #obtain results from earlier simulations
     }
     
     qq0 = opt$quantile #quantile to calculate conservative LR for
     obj = calcLRmcmc(mlehp, mlehd, opt$niter, opt$delta,quantile=qq0,seed=opt$seed,mcmcObjList=mcmcObjList)#, accRateGolden=0.25,diffTol=0.1, niterTune=200, diffSeed=999, verbose=TRUE
     dev.new()      #ENSURE TO GET TRACEPLOT IN ADDITION

     #Obtain output:
     #LRcons = obj$LRcons #obtain estimated conservative LR
     log10LRconsCI = obj$log10LRconsCI #obtain 95% CI
     log10LRcons = obj$log10LRcons
     log10LRbayes = obj$log10LRbayes
     log10LR <- obj$log10LRdistr ##Obtain Random values of log10LR:
     
     print("----Final results from MCMC sampling----")
     print("Estimation of the Bayes Factor (approximate):")
     print(paste0("log10LR=",log10LRbayes))
     #print(paste0("LR=",10^log10LR_bayesian))
     print("Different quantiles of the log10LR distribution:")
     qqs <- c(0.01,0.025,0.05,0.10,0.25,0.5) #quantiles to calculate
     print(quantile(log10LR,qqs))

     #DRAW FIGURE:
     d <- density(log10LR)
     plot(d,xlab="log10 LR",ylab="log10LR distr",main="Sensitivity of LR")
     mtext(paste0("Number of samples: ",obj$niter))
     abline(v=log10(obj$LRmle),lty=2) #MLE based
     #lines(d$x, dnorm(d$x,mean=mean(log10LR),sd=sd(log10LR)),lty=2,col="gray")
     abline(v=log10LRcons,col=4,lty=1)
     abline(v=log10LRconsCI,col=4,lty=2) #also plot 95% CI
     abline(v=log10LRbayes,col=2,lty=2) #show bayesian estimate
     
     txt1 = paste0("MLE=",round(log10(obj$LRmle),2))
     txt2 = paste0("BayesFactor=",round(log10LRbayes,2))
     txt3 = paste0(qq0*100,"% percentile = ",round(log10LRcons,2))
     txt4 = paste0("95% CI=[",paste0(round(log10LRconsCI,2),collapse=","),"]")
     legtxt = c(txt1,txt2,txt3,txt4)
     legend("topright",legend=legtxt, col=c(1,2,4,4),lty=c(2,2,1,2))

     #Update MCMC results and store in setEVID
     res$MCMC <- list(nSamples=obj$niter,log10LRcons=log10LRcons,log10LRconsCI=log10LRconsCI,quantile=qq0,log10LRbayes=log10LRbayes, seed=opt$seed, delta=obj$delta,deltaTuned=obj$deltaTuned, mcmcObjList=obj$mcmcObjList) #Storing niter and seed as specified. storing LR on log10 scale
     assign("resEVID",res,envir=mmTK) #store again
  }

  #Tippet-analysis frame: Draw random (possibly related) non-contributors (specified under Hd) 
  doTippet <- function(tipind,type) { #tipref is index in refData to exchange with random man from population
    #type="MLE"
    niter <- get("optDB",envir=mmTK)$ntippets
    bool = gWidgets2::gconfirm(paste0("Are you sure you want to proceed with\nthe non-contributor analysis?\nThe number of samples is set to ",niter,"\n\nYou can change this number under the Database settings!"))#,icon="question")
    if(!bool) return()
    
    #Obtain calculation options:
    set = get("setEVID",envir=mmTK)
    res = get("resEVID",envir=mmTK)
    MLEopt <- get("optMLE",envir=mmTK) 
    INTopt <- get("optINT",envir=mmTK) 
    LRobs = NULL    
    if(type=="INT") LRobs = res$intLR$log10LR #obtain LR based from Bayes Factor

    #obtain results and store tippet statistics
    tipobj <- calcTippet(tipind,set$mlefit_hp,set$mlefit_hd,niter, type, LRobs, MLEopt, INTopt,seed = NULL)
    res$TIPPET = list(stat = tipobj$stat, type=tipobj$type) #store only stats
    assign("resEVID",res,envir=mmTK) 
  } #end Tippet function

  #helpfunction to translate fitted MLE objects, get and store LR values
  storeLRvalues = function(set) { 
    #obtain calculated LR values with helpfunction (different variants based on MLE)
    obj = calcLRmle(set$mlefit_hp,set$mlefit_hd) 
    resEVID <- list(LRmle=obj$LR,LRlap=obj$LRlap,LRi=obj$LRmarker,LRupper=obj$LRupper,adjLRmle=obj$LRadj) 
    assign("resEVID",resEVID,envir=mmTK) #store EVID calculations for showing later (also report)
    return(resEVID)
  }
  
  #helpfunction to print msg to screen
  #modelfitmsg =function() gWidgets2::gmessage("The one-sample Kolmogorov-Smirnov test\nrejected the peak height model assumption\n(with significance level 0.05)",title="Rejection of model assumption",icon="info")
  doValidMLE = function(mlefit,hyp) { #function to get significance level in Validation plot
    txt = paste0("PP-plot under ",hyp)
    alpha <- get("optMLE",envir=mmTK)$alpha  #as.numeric(getValueUser("Set significance level \nin model validation:",0.01))
    #checkPositive(alpha,"The significance level",strict=TRUE)
    valid = validMLEmodel(mlefit,NULL,txt,alpha=alpha)
    
    #Store validation result here (For further in report)
    nFailed = sum(valid$Significant) #store number of failed points
    objname = paste0("nFailed",hyp)
    resEVID <- get("resEVID",envir=mmTK) 
    resEVID[[objname]] = nFailed
    assign("resEVID",resEVID,envir=mmTK) 
  }   
  
  
  #helpfunction to show topEPG/topMPS fits
  plotTop = function(mlefit){
    kit <- get("selPopKitName",envir=mmTK)[1] #get selected kit: Used in modelvalidation
    sampletype = getSampleType(mlefit$model$samples,kit=kit,LUSsymbol=LUSsymbol) #get sample type (EPG/LUS/MPS) 
    hasPlotly = requireNamespace("plotly") #check if having plotly
    
    if(hasPlotly) { #if using browswer to plot profile
      switch(sampletype,
             "EPG" = plotTopEPG2(mlefit,kit=kit),
             "LUS" = plotTopMPS2(mlefit,grpsymbol=LUSsymbol),
             "MPS" = plotTopMPS2(mlefit,grpsymbol=MPSsymbol))
    } else {
      switch(sampletype, 
             "EPG" = plotTopEPG(mlefit,kitname=kit),
             "LUS" = plotTopLUS(mlefit,LUSsymbol=LUSsymbol),
             "MPS" = gWidgets2::gmessage("Install the package plotly to show MPS plot!",title="Package not found",icon="info"))
    }
  } #end if plotting top
  
###########################################################
#FUNCTION WHICH REFRESHES THE MLEfit panel (After new fit)#
###########################################################
  
  #type="EVID"
  refreshTabMLE = function(type) { 
    #type={"EVID","DB","DC","START"}
    gWidgets2::visible(mainwin) <- FALSE 
    
#################################################################
#helpfunction used to show MLE fit for a given object (Hp or Hd)#
    #mlefit=get("setEVID",envir=mmTK)$mlefit_hp #get setup for EVID
    tableMLE = function(mlefit,tabmleX,sig0=2) {
      SUM_PH  =  unlist(lapply(mlefit$model$samples,function(x)  lapply(x,function(y) sum(y$hdata)))) #get SUM PH for all markers across all replicates
      AVG_PH = mean(SUM_PH) #obtain average PH 
      
      mle <- cbind(mlefit$fit$thetahat2,sqrt(diag(mlefit$fit$thetaSigma2)))
      pnames2 <- rownames(mle) #parameter names
      tab <- cbind(pnames2,format(mle,digits=sig0))
      colNames <- colnames(tab) <-  c("Param.","MLE","Std.Err.")
    
      isCondInd = which(mlefit$model$condOrder>0) #get index of conditional
      nCond = length(isCondInd) #number of conditional contr
      
      NOC = mlefit$model$nC #number of contr
      nRefsMarkers = sapply(mlefit$model$refData,length)
      refNames = names(mlefit$model$refData[[which.max(nRefsMarkers)]]) #obtain reference names

      #show results in table:  
      tabmleX1 = gWidgets2::glayout(spacing=1,container=(tabmleX[1,1] <-gWidgets2::gframe("Parameter estimates:",container=tabmleX,expand=TRUE,fill=TRUE)),expand=TRUE,fill=TRUE) 
      #gWidgets2::gdf(tab,container=tabmleX1,expand=TRUE,fill=TRUE)#,noRowsVisible=TRUE) #add to frame
      #tabmleX1 = gWidgets2::ggroup(spacing=0,container=tabmleX,expand=TRUE,fill=TRUE)
      
      tabmleX1[1,1] = gWidgets2::glabel(colNames[1],container=tabmleX1)
      tabmleX1[1,2] = gWidgets2::glabel(colNames[2],container=tabmleX1)
      tabmleX1[1,3] = gWidgets2::glabel(colNames[3],container=tabmleX1)
      gWidgets2::font(tabmleX1[1,1]) <- gWidgets2::font(tabmleX1[1,2]) <- gWidgets2::font(tabmleX1[1,3]) <- list(weight="bold",size=11)
      
      for(j in 1:nrow(tab)) {
        for(i in 1:3) {
          tabmleX1[j+1,i] = gWidgets2::glabel(tab[j,i],container=tabmleX1)
          
          if(i%in%1:2 && j<=NOC) {
            Mx = as.numeric(tab[j,2]) #obtain Mx
            Mxtxt = paste0("Mix-prop=",format( Mx*100,digits=2),"%")
            Mxtxt2 = paste0("RFU*Mx=",round( AVG_PH*Mx,2))
            if(j<=nCond) {
              helptext(tabmleX1[j+1,i],paste0("Contr: ",refNames[isCondInd[j]],"\n",Mxtxt,"\n",Mxtxt2) )
            } else {
              helptext(tabmleX1[j+1,i],paste0("Contr: Unknown ",j-nCond,"\n",Mxtxt,"\n",Mxtxt2) )
            }
          } #end if first 
        }
        #gWidgets2::font(tabmleX1[j+1,1]) <- list(style="italic")
      }

      tabmleX2 = gWidgets2::glayout(spacing=0,container=(tabmleX[2,1] <-gWidgets2::gframe("Maximum Likelihood value",container=tabmleX))) 
      tabmleX2[1,1] =  gWidgets2::glabel(text="logLik=",container=tabmleX2)
      tabmleX2[1,2] =  gWidgets2::glabel(text=round(mlefit$fit$loglik,sig0),container=tabmleX2)
      tabmleX2[2,1] =  gWidgets2::glabel(text="adj.loglik=",container=tabmleX2) #show adj.loglik=-AIC/2, where AIC= -2*logLik + 2*nparam -AIC/2
      tabmleX2[2,2] =  gWidgets2::glabel(text=round((mlefit$fit$loglik - length(mlefit$fit$thetahat)),sig0),container=tabmleX2)
      #tabmleX2[3,1] =  gWidgets2::glabel(text="Lik=",container=tabmleX2)
      #tabmleX2[3,2] =  gWidgets2::glabel(text= getSmallNumber(mlefit$fit$loglik,sig0),container=tabmleX2)
    }#end function

    #CALCULATIONS BEGIN
    #optimizing options
    opt <- get("optMLE",envir=mmTK) #options when optimizing (nDone,delta)
    dec <- opt$dec #number of significant numbers to have in MLE print
    checkPositive(opt$delta,"Variance parameter of randomizer")
    checkPosInteger(opt$nDone,"Number of random startpoints")
    
    if(type=="START") { #loads already calculated results if program starts
      set <- get("setEVID",envir=mmTK) #get setup for EVID
      if(is.null(set)) set <- get("setDC",envir=mmTK) #get setup for DC if EVID not found
      mlefit_hd <- set$mlefit_hd
      mlefit_hp <- set$mlefit_hp
      if(is.null(mlefit_hd)) return(); #EVID/DC was not prev calculated, return out of function!
    } else { #otherwise, function was called to make new calculations
      #type="EVID"
      #take out relevant parameters from stored list and prepare calculations:
      set <- get( paste0("set",type) ,envir=mmTK) #get setup for EVID/DB/DC
      
      #fit under hp: (only for evidence)
      mlefit_hp <- NULL #not used otherwise
      if(type=="EVID") {
        print("Calculating under Hp...")
        time <- system.time({  mlefit_hp <- calcMLE(set$nC_hp,set$samples,set$popFreq,set$refData,set$condOrder_hp,set$knownref_hp,set$kit,set$DEG,set$BWS,set$FWS, set$threshT, set$prC, set$lambda, set$fst, NULL,NULL,              set$minFreq, set$normalize,  opt$steptol, opt$nDone, opt$delta, opt$difftol, opt$seed, TRUE, set$priorBWS,set$priorFWS, opt$maxThreads, set$adjQbp) })[3]      
        print(paste0("Optimizing under Hp took ",format(time,digits=5),"s"))
        if(!is.null(set$mlefit_hp) && set$mlefit_hp$fit$loglik>mlefit_hp$fit$loglik )  mlefit_hp <- set$mlefit_hp #the old model was better
      }
      
      #fit under hd: (does it for all type of calculations)
      nUhp <- set$nC_hp-sum(set$condOrder_hp>0) #number of unknowns	 
      print("Calculating under Hd...")
      time <- system.time({    mlefit_hd <- calcMLE(set$nC_hd,set$samples,set$popFreq,set$refData,set$condOrder_hd,set$knownref_hd,set$kit,set$DEG,set$BWS,set$FWS, set$threshT, set$prC, set$lambda, set$fst, set$knownRel, set$ibd,  set$minFreq, set$normalize,  opt$steptol, opt$nDone, opt$delta, opt$difftol, opt$seed, TRUE, set$priorBWS,set$priorFWS, opt$maxThreads, set$adjQbp) })[3]
      print(paste0("Optimizing under Hd took ",format(time,digits=5),"s"))
      if(!is.null(set$mlefit_hd) && set$mlefit_hd$fit$loglik>mlefit_hd$fit$loglik )  mlefit_hd <- set$mlefit_hd #the old model was better
      
      #store MLE results:
      set$mlefit_hp=mlefit_hp #store fitted mle-fit
      set$mlefit_hd=mlefit_hd #store fitted mle-fit
      if(type=="EVID") {
        storeLRvalues(set) #store LR valeus based on fitted mle-fit
        assign("setEVID",set,envir=mmTK) #store setup for EVID
      }
      if(type=="DB") assign("setDB",set,envir=mmTK) #store setup for DB
      if(type=="DC") assign("setDC",set,envir=mmTK) #store setup for DC
    } #END OBAINING mlefit
	
#######################
#GUI (common under Hd)#
#######################
    #PROVIDES TWO TABLES:
    tabMLEmeta <- gWidgets2::glayout(spacing=1,container=(tabMLE[1,1] <-gWidgets2::gframe("Evaluation",container=tabMLE))) #main panel
    #tabMLEmeta[1,1] = gWidgets2::gcombobox(1:3,container = tabMLESEL)
    
    #Provide meta info about samples and hypotheses (NOC:
    txtSamples = paste0("Sample(s): ",paste0(names(set$samples),collapse="/"))
    refNames = names(set$refData[[1]])
    condHp = refNames[which(set$condOrder_hp>0)]
    condHd = refNames[which(set$condOrder_hd>0)]
    
    tabMLEmeta[1,1] = gWidgets2::glabel(txtSamples,container = tabMLEmeta)
    if(type!="DC") {
      NOC = set$mlefit_hp$model$nC 
      cond = condHp
      if(type=="DB") {
        NOC = NOC + 1       #need to add one extra if datbase search.
        cond = c("DBref",cond) #add POI first 
      }
      txtHp = paste0("Hp: NumContr=",NOC,". Conditional ref(s): ",ifelse(length(cond),paste0(cond,collapse="/"),"none"))
      tabMLEmeta[2,1] = gWidgets2::glabel(txtHp,container = tabMLEmeta)
    }
    txtHd = paste0("Hd: NumContr=",set$mlefit_hd$model$nC,". Conditional ref(s): ",ifelse(length(condHd),paste0(condHd,collapse="/"),"none"))
    if(!is.null(set$knownRel)) txtHd = paste0(txtHd,". Last unknown related to ",refNames[set$knownRel]) #include last unknown info
    tabMLEmeta[3,1] = gWidgets2::glabel(txtHd,container = tabMLEmeta)
    
    tabMLEtmp <- gWidgets2::glayout(spacing=1,container=(tabMLE[2,1,expand=TRUE] <- gWidgets2::ggroup(container=tabMLE))) #main panel
    tabmleA = gWidgets2::glayout(spacing=0,container=(tabMLEtmp[1,1] <- gWidgets2::gframe("Estimates under Hd",container=tabMLEtmp,expand=TRUE,fill=TRUE)),expand=TRUE,fill=TRUE) 
    tableMLE(mlefit_hd,tabmleA)
    tabmleA3 = gWidgets2::glayout(spacing=0,container=(tabmleA[3,1] <-gWidgets2::gframe("Further Action",container=tabmleA))) 
    tabmleA3[1,1] <- gWidgets2::gbutton(text="MCMC simulation",container=tabmleA3,handler=function(h,...) { doMCMC(mlefit_hd) } )
    tabmleA3[2,1] <- gWidgets2::gbutton(text="Deconvolution",container=tabmleA3,handler=function(h,...) { doDC(mlefit_hd) }  )
    tabmleA3[3,1] <- gWidgets2::gbutton(text="Model validation",container=tabmleA3,handler=function(h,...) { doValidMLE(mlefit_hd,"Hd") } )
    tabmleA3[4,1] <- gWidgets2::gbutton(text="Model fitted P.H.",container=tabmleA3,handler=function(h,...) { 
      plotTop(mlefit_hd) #UPDATED v2.2.0
    } )

#ADD EASY MODE
    if(get("optSetup",envir=mmTK)$easyMode) {
     gWidgets2::enabled(tabmleA3[1,1]) <- FALSE #deactivate MCMC
     #gWidgets2::enabled(tabmleA3[4,1]) <- FALSE #deactivate Model fitted PH
    }
    if( !is.null(mlefit_hp) ) { #used only for weight-of-evidence
     tabmleB = gWidgets2::glayout(spacing=0,container=(tabMLEtmp[1,2] <-gWidgets2::gframe("Estimates under Hp",container=tabMLEtmp,expand=TRUE,fill=TRUE)),expand=TRUE,fill=TRUE) 
     tableMLE(mlefit_hp,tabmleB)
     tabmleB3 = gWidgets2::glayout(spacing=0,container=(tabmleB[3,1] <-gWidgets2::gframe("Further Action",container=tabmleB))) 
     tabmleB3[1,1] <- gWidgets2::gbutton(text="MCMC simulation",container=tabmleB3,handler=function(h,...) { doMCMC(mlefit_hp) } )
     tabmleB3[2,1] <- gWidgets2::gbutton(text="Deconvolution",container=tabmleB3,handler=function(h,...) {  doDC(mlefit_hp) }  )
     tabmleB3[3,1] <- gWidgets2::gbutton(text="Model validation",container=tabmleB3,handler=function(h,...) { doValidMLE(mlefit_hp,"Hp") } )
     tabmleB3[4,1] <- gWidgets2::gbutton(text="Model fitted P.H.",container=tabmleB3,handler=function(h,...) { 
       plotTop(mlefit_hp) #UPDATED v2.2.0
     } )

    #ADD EASY MODE
     if(get("optSetup",envir=mmTK)$easyMode) {
      gWidgets2::enabled(tabmleB3[1,1]) <- FALSE #deactivate MCMC
      #gWidgets2::enabled(tabmleB3[4,1]) <- FALSE #deactivate Model fitted PH
     }
    }


    fixmsg <- "The specified model could not explain the data.\nPlease re-specify the model."
    if(is.infinite(mlefit_hd$fit$loglik)) gWidgets2::gmessage(fixmsg,title="Wrong model specification",icon="error")

    if(type=="EVID")  if(!is.infinite(mlefit_hd$fit$loglik) && is.infinite(mlefit_hp$fit$loglik)) gWidgets2::gmessage(fixmsg,title="Wrong model specification",icon="error")

    if( !is.null(mlefit_hp) ) { #used only for weight-of-evidence
     tabmleC = gWidgets2::glayout(spacing=5,container=(tabMLEtmp[1,3] <-gWidgets2::gframe("",container=tabMLEtmp))) 
     resLR <- get("resEVID",envir=mmTK) #get EVID calculations 
     
      #CREATING NEW LAYOUT:
     lrobs = format(log10(resLR$LRmle),digits=dec) #obtain LR-value to show
     lrupper = format(log10(resLR$LRupper),digits=dec)
     tabmleC1 = gWidgets2::glayout(spacing=0,container=(tabmleC[1,1] <-gWidgets2::gframe("Joint LR",container=tabmleC))) 
     tabmleC1[1,1] =  gWidgets2::glabel(text=paste0("log10LR=",lrobs),container=tabmleC1)
  
     lrobs_subsourcetxt = format(log10(resLR$adjLRmle),digits=dec) #obtain sub-source LR
     #helptext(tabmleC1[1,1],paste0("Sub-source log10LR=",lrobs_subsourcetxt) )      
     tabmleC1[2,1] =  gWidgets2::glabel(text=paste0("Upper boundary=",lrupper),container=tabmleC1)

     #LR-per locus layout (create button): SHOW IN separate window
     tabmleC1[3,1] =  gWidgets2::gbutton(text="Show LR per-marker",container=tabmleC1, handler=function(h,...) {
       sig0 = 2
       LRmarker = resLR$LRi #obtain Locus specific LR:
       outD = cbind(Marker=names(LRmarker), LR=format(LRmarker,digits=sig0),log10LR=format(log10(LRmarker),digits=sig0))
       dbwin <- gWidgets2::gwindow("LR per-marker", visible=FALSE)#,height=mwH)
       tab <- gWidgets2::gdf(NAtoSign(outD) ,container=dbwin) #create table
       gWidgets2::visible(dbwin) <- TRUE
     })
     
     #We show weight-of-evidence
     tabmleD = gWidgets2::glayout(spacing=5,container=(tabmleC[3,1] <-gWidgets2::gframe("Further",container=tabmleC))) 
     #tabmleD[1,1] <- gWidgets2::gbutton(text="Optimize more",container=tabmleD,handler=function(h,...) { refreshTabMLE(type)  } )
     tabmleD[1,1] <- gWidgets2::gbutton(text="LR sensitivity",container=tabmleD,handler=function(h,...) { simLR(mlefit_hp,mlefit_hd) } ) 
     tabmleD[2,1] <- gWidgets2::gbutton(text="Bayes Factor",container=tabmleD,handler=function(h,...) { doINT() } )  #trigger integer alculation
     if(get("optSetup",envir=mmTK)$easyMode) gWidgets2::enabled(tabmleD[2,1]) <- FALSE #deactivate 
     helptext(tabmleD[1,1],"Calculates the conservative LR based on MCMC simulations. Also giving an estimate of Bayes Factor")
     helptext(tabmleD[2,1],"Calculates the Full Bayesian LR (Bayes Factor) using the Integrated approach")
     
     #postanalysis
     tabmleF = gWidgets2::glayout(spacing=3,container=(tabmleC[2,1] <-gWidgets2::gframe("Non-contributor analysis",container=tabmleC))) 
     tippets <- setdiff(set$knownref_hd,set$knownref_hp)  #known non-contributors under Hd but not Hp (index)
     if(!is.null(tippets)) {
       tN <- names(set$refData[[1]][tippets]) #tippet names
       tabmleF[1,1] <- gWidgets2::glabel( "Select reference to\nreplace with non-contributor:",container=tabmleF)
       tabmleF[2,1] <- gWidgets2::gcombobox( items=tN ,container=tabmleF)
       tabmleF[3,1] <- gWidgets2::gbutton(text="Sample MLE based",container=tabmleF,handler=function(x) {
         doTippet(tipind=tippets[which(tN==gWidgets2::svalue(tabmleF[2,1]))],type="MLE")  #get tip-index in refData
       })
       tabmleF[4,1] <- gWidgets2::gbutton(text="Sample Bayesian based",container=tabmleF,handler=function(x) { 
         doTippet(tipind=tippets[which(tN==gWidgets2::svalue(tabmleF[2,1]))],type="INT")  #get tip-index in refData
       })
       helptext(tabmleF[3,1],"Replaces the POI with a random man and calculates the LR using the MLE approach")
       helptext(tabmleF[4,1],"Replaces the POI with a random man and calculates the Bayes Factor using the Integrated approach")
       if(get("optSetup",envir=mmTK)$easyMode) gWidgets2::enabled(tabmleF[4,1]) <- FALSE #deactivate tippets for bayesian approach
     } #end if tippets
     
     tabmleD[3,1] <- gWidgets2::gbutton(text="Create report",container=tabmleD,handler=function(h,...) { efm_createReport(mmTK,type)})
    }
    if(type=="DB") tabmleA3[2,2] <- gWidgets2::gbutton(text="Search Database",container=tabmleA3,handler=function(h,...) { doDBsearch("QUAN")} )

    gWidgets2::visible(mainwin) <- TRUE
    gWidgets2::focus(mainwin) <- TRUE #focus window after calculations are done
  } #end refresh tab-frame of MLE-fit

  refreshTabMLE(type="START") #Show already calculted evidence-results when program starts


##############################################################
###############Tab 5: Deconvolution results:##################
##############################################################

 f_savetableDC = function(h,...) {
   if(is.null(DCtable)) {
    gWidgets2::gmessage("There is no deconvolution results available!")
   } else {
    tableSaver(DCtable[], "txt") #save deconvolution results
   }
 }
  
  f_saveDCasRef= function(h,...) {
    DCtables <- get("resDC",envir=mmTK) #get deconvolved results
    topRanked = DCtables[[1]] #get top ranked table
    ncol = ncol(topRanked)
    numRefs = (ncol-1)/3 #number of references to export
    outtab = NULL #get outtabe
    for(r in seq_len(numRefs)) {
      genos = topRanked[,3*(r-1) + 2] #obtain genotypes
      genos =  t(matrix(unlist(strsplit(genos,"/")),nrow=2)) #obtain alleles per genotyes
      new = cbind(paste0("C",r),topRanked[,1],genos)
      outtab = rbind(outtab,new)
    }
    colnames(outtab) = c("SampleName","Marker","Allele 1","Allele 2")
    tableSaver(outtab,"csv")
  }
  
  refreshTabDC = function(dctype=1) { #1=table1 (top marginal results),2=table2 (joint results), 3=table3 (all marginal results per genotype), 4=table4 (all marginal results per alleles)
   DCtables <- get("resDC",envir=mmTK) #get deconvolved results
   if(!is.null(DCtables)) DCtable[] = NAtoSign(DCtables[[dctype]])  #update Table
 }

 #CREATE DECONV-GUI 
 tabDCa = gWidgets2::glayout(spacing=1,container=tabDC) #table layout
 tabDCb = gWidgets2::ggroup(spacing=1,container=tabDC,expand=TRUE,fill=TRUE)
 itemvecDC = c("Top Marginal","All Joint","All Marginal (G)","All Marginal (A)")
 tabDCa[1,1] <- gWidgets2::glabel("Select layout:",container=tabDCa)
 tabDCa[1,2] <-  gWidgets2::gradio(items=itemvecDC,selected=1,horizontal=TRUE,container=tabDCa,handler=function(x) {
   refreshTabDC( which(itemvecDC==gWidgets2::svalue(tabDCa[1,2])) )
 })
 tabDCa[2,1] <- gWidgets2::gbutton(text="Save table",container=tabDCa,handler=f_savetableDC)  
 tabDCa[2,2] <- gWidgets2::gbutton(text="Save Top Ranked Genotype as Reference",container=tabDCa,handler=f_saveDCasRef)  
 helptext(tabDCa[2,2],"Store top ranked genotypes to a reference file (EFM-format)")
 
 #ADD DC TABLE
 DCtable = gWidgets2::gtable(items="",multiple = TRUE,container = tabDCb,expand=TRUE,fill=TRUE)
 gWidgets2::add(tabDCb,DCtable,expand=TRUE,fill=TRUE)
 refreshTabDC() #open results when program starts



##############################################################
###############Tab 6: Database search:########################
##############################################################

 f_savetableDB = function(h,...) {
   if(is.null(DBtable)) {
     gWidgets2::gmessage("There is no deconvolution results available!")
   } else {
     tableSaver(DBtable[], "txt") #save deconvolution results
   }
 }

 refreshTabDB = function(ranktype=1) {
   DBtable2 <- get("resDB",envir=mmTK) #get database result 
   if(!is.null(DBtable2)) {
     if(nrow(DBtable2)<=1 )  {
       DBtable[] <- DBtable2
     } else { #several candidates
       colRank = ranktype+1 #column used to rank
       if(DBtable2[1,colRank]=="-") return() #don't rank
       
       vals = as.numeric(DBtable2[,colRank]) #obtain values
       ord <- order(vals,decreasing=TRUE) #need to convert to numeric!
       sizeDB = min(length(ord), get("optDB",envir=mmTK)$maxDB) #max list
       ord = ord[seq_len(sizeDB)] #trim ordered idx
       DBtable2 = cbind(Rank=seq_len(sizeDB),DBtable2[ord,,drop=FALSE]) #add rank to table
       DBtable[] <- DBtable2  
     }
   }
 }

 #Create table:
 tabDBa = gWidgets2::glayout(spacing=1,container=tabDB) #table layout
 tabDBb = gWidgets2::ggroup(spacing=1,container=tabDB,expand=TRUE,fill=TRUE) #table layout
 tabDBa[1,1] <- gWidgets2::glabel("Sort table:",container=tabDBa)
 itemvecDB <- c("quanLR","qualLR","MAC","nMarkers")
 tabDBa[1,2] <- gWidgets2::gradio(items=itemvecDB,selected=1,horizontal=TRUE,container=tabDBa,handler=function(x) {
   refreshTabDB( which(itemvecDB==gWidgets2::svalue(tabDBa[1,2])) )
 })
 tabDBa[2,1] <- gWidgets2::gbutton(text="Save table",container=tabDBa,handler=f_savetableDB)  
 
 #ADD DB TABLE
 DBtable = gWidgets2::gtable(items="",multiple = TRUE,container = tabDBb,expand=TRUE,fill=TRUE)
 gWidgets2::add(tabDBb,DBtable,expand=TRUE,fill=TRUE) #add table
 refreshTabDB(2) #when program starts: Consider qual-rank

 
 
###############################################################
###############Tab 7: LRmix module:############################
###############################################################
 #uses only qualitative information

  f_savetableEVIDLRMIX = function(h,...) { #function for storing LR
    LRi <- get("resEVIDLRMIX",envir=mmTK) #get EVID calculations when GUI starts
    if(is.null(LRi)) {
      gWidgets2::gmessage("There was no Weight-of-Evidence results available!")
    } else {
      tab <- c(LRi,prod(LRi))
      tab <- cbind(c(names(LRi),"Joint"),tab,log10(tab))
      colnames(tab) <- c("Locus","LR","log10LR")
      tableSaver(tab, "txt") 
    }
  }
 
  #Refreshing LR-GUI  
  refreshTabLRMIX = function() {
    tabLRMIXtmp <- gWidgets2::glayout(spacing=spc,container=(tabLRMIX[1,1,expand=TRUE] <- gWidgets2::ggroup(container=tabLRMIX))) 
    gWidgets2::visible(mainwin) <- FALSE
    sig0 = 4
    
    #OBTAIN DATA/Hypo/Settings FOR ANALYSIS
    set <- get("setEVID",envir=mmTK) #get setup for EVID
    mod <- set$model #obtain model options
    par <- set$param  #obtain model parameters
    
    #helpfunction to update LR calculed value
    updateQualLR = function(jointLR, sig=4) { 
      log10LR = sum(log10(jointLR))
      tabLRmixB[1,2] = gWidgets2::glabel(text=format(10^log10LR,digits=sig),container=tabLRmixB)
      tabLRmixB[2,2] = gWidgets2::glabel(text=round(log10LR,2),container=tabLRmixB)
    }
  
    #helpfunction for calculating LR for each given dropout pD (takes a numeric): REQUIRE DATA TO BE PREPARED!
    doLR = function(pD) {
      hpvec <- hdvec <- rep(1,length(locs))
      for(loc in locs) {
        fst0 = getMarkerVal(set$fst,loc) #extract per marker based fst
        pC0 = getMarkerVal(set$prC,loc) #extract per marker based drop-in prob 
        hpvec[which(loc==locs)] <- calcQual( evidList[[loc]], refList_hp[[loc]]$Ri, refList_hp[[loc]]$Ki, set$nC_hp-refList_hp[[loc]]$nRefs, fst0,  pD, pC0, popFreqList[[loc]])
        hdvec[which(loc==locs)] <- calcQual( evidList[[loc]], refList_hd[[loc]]$Ri, refList_hd[[loc]]$Ki, set$nC_hd-refList_hd[[loc]]$nRefs, fst0,  pD, pC0, popFreqList[[loc]])
      }
      LRi <- hpvec/hdvec
      names(LRi) <- locs
      return(LRi)
    }
  
    #Function to calculate the MLE based LR (qual):
    f_calcQualLRmle = function(h,...) { #added function v1.11
      mlehp = calcQualMLE(set$nC_hp, set$samples, set$popFreqQ, set$refDataQ, set$condOrder_hp, set$knownref_hp, set$prC,set$fst)
      mlehd = calcQualMLE(set$nC_hd, set$samples, set$popFreqQ, set$refDataQ, set$condOrder_hd, set$knownref_hd, set$prC,set$fst)
      pDhat <- c(mlehp$pDhat,mlehd$pDhat)#, 1/(1+exp(-c(foohp$est,foohd$est))) #get estimated dropouts
      logLikHp = mlehp$loglik
      logLikHd = mlehd$loglik
      log10LRmle =   (logLikHp - logLikHd)/log(10) #get LR on log10 scale
      
      print(paste0("logLik_hp=",logLikHp)) #Print maximum likelihood under Hp 
      print(paste0("logLik_hd=",logLikHd)) #Print maximum likelihood under Hp 
      v1 <- signif(log10LRmle,digits=sig0)
      v2 <- signif(10^log10LRmle,digits=sig0)
      v3 <- signif(pDhat[1],digits=sig0)
      v4 <- signif(pDhat[2],digits=sig0)
      vectxt =  c("Results of MLE based LR:",paste0("LR=",v2),paste0("log10LR=",v1),paste0("Estimated pD under Hp=",v3),paste0("Estimated pD under Hd=",v4))
      txt = paste0(paste0(vectxt,collapse="\n"), "\n\nDo you want to export the results to a text file?")
      bool <- gWidgets2::gconfirm(txt,title="Maximum Likelihood based qualitative LR")
      if(bool) tableSaver(vectxt,sep="txt")
    } #end qualLRmleBased
    
    f_calcSensitivity = function(h,...) {
      optLRMIX <- get("optLRMIX",envir=mmTK) #options when integrating (reltol and boundaries)
      range <- optLRMIX$range
      nTicks <- optLRMIX$nticks
      checkProb(range,"The range")
      checkPosInteger(nTicks, "The number of ticks")  
      pDvec <- seq(1e-6,range,l=nTicks) #set very small dropout probability
      LRivec <- Vectorize(doLR) (pDvec) #return LR(pD,loc)
      logLR <- colSums(log10(LRivec)) #joint LR
      ylims <- range(logLR[!(is.nan(logLR) | is.infinite(logLR))]) 
      if(ylims[2]>-1) ylims[1] <- -1  #lower limit
      plot(pDvec ,logLR ,xlab="Pr(Dropout)",ylab="log_10(LR)",ylim=ylims + c(0,0.5),main="Sensitivity plot")
      lines(spline(pDvec ,logLR,n=nTicks )) #use spline to interpolate points
      abline(h=0,col="gray")
    }
    
    f_calcCons = function(h,...) {
      noSamples = function(hyp,M) { #helpfunction for tell user that wrong model assumption was used.
        gWidgets2::gmessage(paste0("No samples was accepted out of the first ",M," samples.\nPlease retry sampling or change hypothesis ",hyp),title="Wrong model specification",icon="error")
        return()
      }  
      optLRMIX <- get("optLRMIX",envir=mmTK) 
      nsample <- optLRMIX$nsample
      alpha <- optLRMIX$alpha
      qq <- c(alpha,0.5,1-alpha) #Dropout quantiles to consider 
      totA <-  sapply(  set$samples, function(x) sum(sapply(x,function(y) length(y$adata)) ) ) #number of alleles for each samples
      print("Total number of observed alleles for sample(s):")
      .printTable(totA)
      refHp <- refHd <- NULL
      if( any(set$condOrder_hp>0) ) refHp <- lapply(set$refData ,function(x) x[set$condOrder_hp]) #must have the original refData!
      if( any(set$condOrder_hd>0) ) refHd <- lapply(set$refData ,function(x) x[set$condOrder_hd]) #must have the original refData!
      
      for(ss in 1:nS) { #for each sample (do dropout calculation)
        print(paste0("For evidence ",names(set$samples)[[ss]],":"))
        print("Estimating quantiles from allele dropout distribution under Hp...")
        Msamp <- max(2000,25*totA[ss]) #number of samples for each vectorization
        DOdist <- simDOdistr(totA=totA[ss],nC=set$nC_hp,popFreq=set$popFreq,refData=refHp,minS=nsample,prC=set$prC,M=Msamp)
        if(length(DOdist)==0) noSamples("Hp",Msamp)
        qqhp <- quantile(DOdist ,qq) #get estimated quantiles
        .printTable(qqhp)
        print("Estimating quantiles from allele dropout distribution under Hd...")
        DOdist <- simDOdistr(totA=totA[ss],nC=set$nC_hd,popFreq=set$popFreq,refData=refHd,minS=nsample,prC=set$prC,M=Msamp)
        if(length(DOdist)==0) noSamples("Hd",Msamp)
        qqhd <- quantile(DOdist ,qq) #get estimated quantiles
        .printTable(qqhd)
        if(ss==1) consDO <- rbind(qqhp[-2],qqhd[-2])
        if(ss>1) {
          if(qqhp[1]<consDO[1,1]) consDO[1,1] <- qqhp[1]
          if(qqhp[3]>consDO[1,2]) consDO[1,2] <- qqhp[3]
          if(qqhd[1]<consDO[2,1]) consDO[2,1] <- qqhd[1]
          if(qqhd[3]>consDO[2,2]) consDO[2,2] <- qqhd[3]
        }
      }
      consDO <- signif(consDO,2) #2 significant numbers to use
      print(consDO) #print out drop-out values
      
      #calculate corresponding LR's
      consLRi <- Vectorize(doLR) (c(consDO))
      conslog10LR <- colSums(log10(consLRi)) #get joint log10LR
      selInd <- which.min(conslog10LR) #select reporting LR (which gives min LR)
      selDO <- consDO[selInd] #take out selected dropout
      LRi <- consLRi[,selInd] #take out selected LR-values
      gWidgets2::svalue(tabLRmixA2[1,2]) <- selDO  #update selected dropout
      assign("resEVIDLRMIX",LRi,envir=mmTK) #assign evidence weighting results - Based on LRmix
      updateQualLR(LRi) #update table with calculated LR
    }
    
    #Helpfunction for calculating non-contributor LR:
    calcNonContr = function(tipsel,pD) {
      #First step: calculating LR for all genotypes in original popFreq.
      Glist <- getGlist(popFreqList) #get random man-Glist (use encoded values)
      print("Precalculating for non-contributor plot...")
  
      #calculate LRs directly here: 
      tipind <- set$knownref_hd[tipsel] #get tip-ind in refData
      modtipind <- set$condOrder_hp[tipind] #get position in system of tippet. Necessary for QUAL model
      for(loc in locs) { #Calcualte for each locus:
        fst0 = getMarkerVal(set$fst,loc) #extract per marker based fst
        pC0 = getMarkerVal(set$prC,loc) #extract per marker based drop-in prob 
        nG <- length(Glist[[loc]]$Gprob) #number of genotypes
        Glist[[loc]]$LR <- rep(NA,nG) #init space for LR
        refhptmp <- refList_hp[[loc]]$Ri  #take out contributing replicates under Hp
        nrefhdtmp <- refList_hd[[loc]]$Ki  #take out non-contributing replicates under Hd
        for(j in 1:nG) { #for each genotypes
          refhptmp[ 2*modtipind -c(1,0) ] <- Glist[[loc]]$G[j,] #insert genotype to reference
          nrefhdtmp[ 2*tipsel-c(1,0) ] <- Glist[[loc]]$G[j,] #insert genotype to reference (noncontributor)
          nUhp = set$nC_hp-length(refhptmp)/2
          nUhd = set$nC_hd-refList_hd[[loc]]$nRefs
          hp0 <- calcQual( evidList[[loc]], refhptmp,            refList_hp[[loc]]$Ki, nUhp ,fst0,pD,pC0,popFreqList[[loc]])
          hd0 <- calcQual( evidList[[loc]], refList_hd[[loc]]$Ri,nrefhdtmp,            nUhd ,fst0,pD,pC0,popFreqList[[loc]])
          Glist[[loc]]$LR[j] <- hp0/hd0 #store LR
        }
      } #end for each locus
      
      #Obtain LR for non-contributors and POI
      nSamples <- get("optDB",envir=mmTK)$ntippets #get number of tippets to simulate
      print(paste0("Simulating ",nSamples," non-contributors..."))
      lr0 <- NULL #log10 LR
      LRi <- get("resEVIDLRMIX",envir=mmTK) #get already calcualted LR values
      if(!is.null(LRi))  lr0 <- sum(log10(LRi))
      RMLR <- rep(0,nSamples) #vector of tippets
      for(loc in locs) {
        if(all(round(Glist[[loc]]$LR,8)==1)) next #skip if only LR=1
        RMLR <- RMLR + sample(log10(Glist[[loc]]$LR),nSamples,replace=TRUE,prob=Glist[[loc]]$Gprob)
      }
      stat = plotTippet(RMLR,lr0,mtxt=paste0("Qualitative-based with pD=",format(pD,digits=3)))
      print(stat) #print tippet statistics
    } #end calc TIPPETs
    
    #helpfunction to get conditional refs under a hypothesis
    getConds <- function(condOrder,knownref) {
      cond <- which(condOrder>0) #ind of conditional refs (they are increasingly sorted)
      Ri <- Ki <- NULL
      for(rr in cond ) Ri <- c(Ri,set$refDataQ[[loc]][[rr]])
      for(rr in knownref) Ki <- c(Ki,set$refDataQ[[loc]][[rr]])
      return(list(Ri=Ri,Ki=Ki,nRefs=length(Ri)/2 ))  #Updated in v1.9: consider number of alleles divided by 2
    }
  
    #Prepare Evidence and refs under each hypothesis: MUST BE GLOBALLY AVAILABLE
    #UPDATED FROM v1.11: STRINGS ARE ENCODED AS integers 1:n (with n as number of alleles), since LRmix require numerics. 
    locs <- names(set$popFreqQ) #get analysing loci
    nS <- length(set$samples) #number of samples
    popFreqList = set$popFreqQ #make copy of frequencies (IMPORTANT!)
    evidList <- list()
    refList_hp <- list()
    refList_hd <- list()
    for(loc in locs) {
      Ei <- NULL #get evidence
      avOld <- names(popFreqList[[loc]]) #old names
      names(popFreqList[[loc]]) <- seq_along(popFreqList[[loc]])  #encoded alleles (used for order of alleles) (popFreq2)
      for(ss in 1:nS) {
       if(ss>1) Ei <- c(Ei,0) #seperate with 0  
       adata <- set$samples[[ss]][[loc]]$adata
       adata <- match(adata,avOld) #update names to be index of popFreq2
       if(length(adata)==0) adata=0 #is empty
       Ei <- c(Ei,adata)
      } 
      evidList[[loc]] <- Ei
      refList_hp[[loc]] <- getConds(condOrder=set$condOrder_hp,knownref=set$knownref_hp) #under hp
      refList_hd[[loc]] <- getConds(set$condOrder_hd,set$knownref_hd) #under hd
    
      if(!is.null(refList_hp[[loc]]$Ri)) refList_hp[[loc]]$Ri <-  match(refList_hp[[loc]]$Ri,avOld) #update names to be index of popFreq2
      if(!is.null(refList_hd[[loc]]$Ri)) refList_hd[[loc]]$Ri <-  match(refList_hd[[loc]]$Ri,avOld) #update names to be index of popFreq2
      if(!is.null(refList_hp[[loc]]$Ki)) refList_hp[[loc]]$Ki <-  match(refList_hp[[loc]]$Ki,avOld) #update names to be index of popFreq2
      if(!is.null(refList_hd[[loc]]$Ki)) refList_hd[[loc]]$Ki <-  match(refList_hd[[loc]]$Ki,avOld) #update names to be index of popFreq2
      if(length(refList_hp[[loc]]$Ri)==0) refList_hp[[loc]]$Ri <- NULL #Updated in v1.9: Missing markers are supported for LRmix
      if(length(refList_hd[[loc]]$Ri)==0) refList_hd[[loc]]$Ri <- NULL #Updated in v1.9: Missing markers are supported for LRmix
    }
  
    #GUI:
    tabLRmixA = gWidgets2::glayout(spacing=0,container=(tabLRMIXtmp[1,1] <-gWidgets2::gframe("Analysis of qualitative LR",container=tabLRMIXtmp))) 
    tabLRmixA1 = gWidgets2::glayout(spacing=0,container=(tabLRmixA[1,1] <-gWidgets2::gframe("Preanalysis",container=tabLRmixA)))  
    tabLRmixA2 = gWidgets2::glayout(spacing=0,container=(tabLRmixA[2,1] <-gWidgets2::gframe("Calculation",container=tabLRmixA))) 
    tabLRmixA3 = gWidgets2::glayout(spacing=0,container=(tabLRmixA[3,1] <-gWidgets2::gframe("Non-contributor analysis",container=tabLRmixA))) 
    tabLRmixA4 = gWidgets2::glayout(spacing=0,container=(tabLRmixA[4,1] <-gWidgets2::gframe("MLE based",container=tabLRmixA)))  #Frame added in v1.11

    #analysis buttons
    tabLRmixA1[1,1] <- gWidgets2::gbutton(text="Sensitivity",container=tabLRmixA1,handler= f_calcSensitivity) #get range from options under toolbar:
    tabLRmixA1[2,1] <- gWidgets2::gbutton(text="Conservative LR",container=tabLRmixA1,handler= f_calcCons)
    tabLRmixA2[1,1] <-  gWidgets2::glabel(text="Dropout prob:",container=tabLRmixA2) #direct LR-Calculation
    tabLRmixA2[1,2] <-  gWidgets2::gedit(text="0.05",container=tabLRmixA2) #this is updated after dropout distr is ran
    gWidgets2::size(tabLRmixA2[1,2]) <- 8
  
    tabLRmixA2[2,1] <-  gWidgets2::gbutton(text="Calculate LR",container=tabLRmixA2,handler=function(x) {
      pD <- as.numeric(gWidgets2::svalue(tabLRmixA2[1,2]))
      checkProb(pD,"The allele dropout probability")
      LRi <- doLR(pD) 
      assign("resEVIDLRMIX",LRi,envir=mmTK) #assign evidence weighting results - Based on LRmix
      updateQualLR(LRi) #update table with calculated LR
    }) 
    tabLRmixA2[2,2] <-  gWidgets2::gbutton(text="Save table",container=tabLRmixA2,handler=f_savetableEVIDLRMIX)
  
    #Tippet-analysis frame:
    tippets <- set$knownref_hd #known non-contributors under Hd
    if(!is.null(tippets)) {
      tippetNames <- names(set$refData[[1]][tippets]) #tippet names
      tabLRmixA3[1,1] <- gWidgets2::glabel( "Select reference to\nreplace with non-contributor:",container=tabLRmixA3)
      tabLRmixA3[2,1] <- gWidgets2::gcombobox( items=tippetNames ,container=tabLRmixA3)
      gWidgets2::size(tabLRmixA3[2,1]) <- 10
      tabLRmixA3[3,1] <- gWidgets2::gbutton(text="Sample non-contributors",container=tabLRmixA3,handler= function(h,...) {
        pD <- as.numeric(gWidgets2::svalue(tabLRmixA2[1,2])) #take dropout-value as given in GUI
        tipref <- gWidgets2::svalue(tabLRmixA3[2,1]) #get name of reference to tippet
        tipsel <- which(tippetNames==tipref) #index of tippet to select     
        calcNonContr(tipsel,pD) #calculate non-contributor
      })
    } #end if not tippet possible
    tabLRmixA4[1,1] <-  gWidgets2::gbutton("Calculate LR",container=tabLRmixA4,handler=f_calcQualLRmle) 

    #Show calculated LR:    
    tabLRmixB = gWidgets2::glayout(spacing=0,container=(tabLRMIXtmp[1,2] <-gWidgets2::gframe("Joint LR",container=tabLRMIXtmp))) 
    tabLRmixB[1,1] =  gWidgets2::glabel(text="LR=",container=tabLRmixB) #direct LR-Calculation
    tabLRmixB[2,1] =  gWidgets2::glabel(text="log10LR=",container=tabLRmixB) #direct LR-Calculation

    tabLRmixB[3,1] =  gWidgets2::gbutton("Show\nper-marker",container=tabLRmixB, handler=function(h,...) {
      LRmarker <- get("resEVIDLRMIX",envir=mmTK) #get already calcualted LR values
      if(is.null(LRmarker)) return() #not obtained  
      outD = cbind(Marker=names(LRmarker), LR=round(LRmarker,2),log10LR=round(log10(LRmarker),2))
      dbwin <- gWidgets2::gwindow("LR per-marker", visible=FALSE)#,height=mwH)
      tab <- gWidgets2::gdf(NAtoSign(outD) ,container=dbwin) #create table
      gWidgets2::visible(dbwin) <- TRUE
    })
    
    getFocus()
  } #end refresh funtion
   
  #END OF PROGRAM
  for(nbvisit in 6:2) gWidgets2::svalue(nb) <- nbvisit #visit tables 
  getFocus()

 })
} #end funcions
