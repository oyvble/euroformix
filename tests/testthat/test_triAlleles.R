#Testing that the numerical calculation of loglik is correct
#rm(list=ls());library(euroformix);library(testthat)
nDone0=3
steptol0=1e-6
s0 = 3 #signif of checking

kit0 = "testkit" #name of selected kit
examples = paste(path.package("euroformix"),"examples_triAllele",sep=.Platform$file.sep)
testfolds = list.dirs(examples,recursive = FALSE)

#Store expected results:
expectedResults_logLik = c(-190.869,-197.709,-205.752,-205.185,-141.419,-346.799,-351.533)
expectedResults_valid = list()
expectedResults_valid[[1]] = c(0.031,0.586,0.91,0.209,0.956,0.157,0.389,0.073,0.662,0.978,0.475,0.569,0.964,0.217,0.071,0.594,0.729,0.682,0.316,0.294,0.304,0.684,0.646,0.241,0.328)
expectedResults_valid[[2]] = c(0.028,0.588,0.913,0.536,0.204,0.959,0.153,0.383,0.068,0.666,0.979,0.474,0.567,0.967,0.213,0.066,0.596,0.73,0.682,0.313,0.291,0.298,0.683,0.649,0.237,0.326)
expectedResults_valid[[3]] = c(0.029,0.635,0.909,0.533,0.229,0.963,0.099,0.18,0.382,0.07,0.696,0.977,0.511,0.563,0.971,0.244,0.067,0.629,0.725,0.677,0.349,0.326,0.298,0.677,0.679,0.27,0.362)
expectedResults_valid[[4]] = c(0.027,0.601,0.918,0.542,0.207,0.964,0.361,0.153,0.386,0.067,0.675,0.982,0.482,0.573,0.97,0.215,0.064,0.605,0.737,0.689,0.317,0.294,0.299,0.693,0.658,0.239,0.33)
expectedResults_valid[[5]] = c(0.614,0.813,0.145,0.881,0.182,0.076,0.205,0.765,0.2,0.822,0.819,0.401,0.903,0.366,0.809,0.509,0.934,0.781,0.512,0.714,0.781,0.106,0.073,0.373)
expectedResults_valid[[6]] = c(0.782,0.263,0.619,0.271,0.94,0.846,0.281,0.133,0.079,0.215,0.195,0.023,0.401,0.448,0.601,0.372,0.988,0.928,0.721,0.967,0.355,0.097,0.745,0.476,0.262,0.606,0.103,0.812,0.949,0.402,0.399,0.717,0.32,0.859,0.604,0.037,0.093,0.654,0.74,0.02,0.953,0.678,0.548,0.721,0.869,0.652,0.179,0.257,0.776,0.741,0.129,0.785,0.653,0.969,0.896,0.598,0.623,0.346,0.563)
expectedResults_valid[[7]] = c(0.681,0.208,0.597,0.283,0.984,0.955,0.264,0.134,0.059,0.16,0.161,0.018,0.525,0.567,0.531,0.319,0.978,0.906,0.668,0.944,0.322,0.088,0.727,0.484,0.232,0.536,0.116,0.769,0.933,0.421,0.554,0.819,0.318,0.818,0.522,0.048,0.106,0.725,0.793,0.036,0.937,0.623,0.482,0.647,0.853,0.654,0.196,0.268,0.838,0.66,0.116,0.796,0.678,0.947,0.859,0.726,0.746,0.342,0.536)
expectedResults_DC = list()
expectedResults_DC[[1]] = list(D3S1358=c(1,"16/18/19"))
expectedResults_DC[[2]] = list(D3S1358=c(1,"16/18/19"),VWA=c(2,"16/18/20"))
expectedResults_DC[[3]] = list(D3S1358=c(1,"16/18/19"),VWA=c(2,"16/18/20"),D16S539=c(1,"10/11/12"))
expectedResults_DC[[4]] = list(D3S1358=c(1,"16/18/19"),VWA=c(2,"16/18/20"),D16S539=c(1,"10/11/12"),D16S539=c(2,"10/11/12"))
expectedResults_DC[[5]] = list(VWA=c(2,"16/18/99"))
expectedResults_DC[[6]] <- expectedResults_DC[[7]] <- list(D16S539=c(1,"10/11/12"))


for(testIdx in seq_along(testfolds)) {
  test_that(paste0("TriAllele example nr ",testIdx,":"), {

#    testIdx = 6
    testfold = testfolds[testIdx]
    isRepStuttFold = basename(testfold)=="test_triAlleleRepsStutters"
    #read data:
    evidfn = "test_evid.csv"
    if(isRepStuttFold)    evidfn = "test_evids.csv"
    popfn = paste(testfold,paste0("test_freq.csv"),sep=.Platform$file.sep)
    evidfn = paste(testfold,evidfn,sep=.Platform$file.sep)
    reffn = paste(testfold,paste0("test_refs.csv"),sep=.Platform$file.sep)
    
    #Obtain data (taken from runexample):
    popFreq = freqImport(popfn)[[1]] #obtain list with population frequencies
    samples = sample_tableToList(tableReader(evidfn))#,threshT)
    refData = sample_tableToList(tableReader(reffn))
  
    #Set settings (default):
    AT = 50
    pC = 0.05
    lambda = 0.01
    fst = 0
    NOC = 2
    condOrder = 1:2
    
    makeTest = function(DEG=FALSE,BWS=FALSE,FWS=FALSE) {
      
#  BWS=FALSE;FWS=FALSE;DEG=FALSE
#  BWS=TRUE;FWS=TRUE;DEG=TRUE
#  BWS=TRUE;FWS=FALSE;DEG=TRUE
      mle = calcMLE(nC=NOC,samples=samples,popFreq=popFreq,refData=refData,condOrder=condOrder,pC=pC,AT=AT,fst=fst,lambda=lambda,kit=kit0,steptol=steptol0,nDone=nDone0,DEG=DEG,BWS=BWS,FWS=FWS)

      #mle$prepareC$triAlleles
      useResIdx = testIdx
      if(BWS && !FWS) useResIdx = testIdx+ 1 #shift with one
      expect_equal(round(mle$fit$loglik,s0),expectedResults_logLik[useResIdx]   ) 
      
      #CHECK deconvolution (DC):
      DC = deconvolve(mle)
      expectedResults_DC_idx = expectedResults_DC[[useResIdx]]
      for(eidx in seq_along(expectedResults_DC_idx)) {
#        eidx = 1
        DCelem =expectedResults_DC_idx[eidx]
        marker = names(DCelem)
        contrIdx = as.integer(DCelem[[1]][1])
        genoNameExp = DCelem[[1]][2] #expected 
        genoNameObs = DC$toprankGi[[marker]][1,contrIdx] #observedf
        expect_equal(genoNameExp,genoNameObs)
      }
      
      #CHECK validation
      valid = validMLEmodel(mle,kit=kit0,createplot=FALSE,verbose = FALSE)
      expect_equal(round(valid$ProbObs,s0),expectedResults_valid[[useResIdx]])    
      #print(paste0(round(valid$ProbObs,s0),collapse = ","))

    }
    
    if(isRepStuttFold) {
      makeTest(TRUE,TRUE,TRUE)
      makeTest(TRUE,TRUE,FALSE)
    } else {
      makeTest() #assume no DEG, FWS, nor BWS 
    }
  })

} #end for each folder
    
