
# EuroForMix <img src="man/figures/efmlogo.png" align="right" height=140/>

Software for the interpretation of complex DNA profiles based on the
information in peak heights. The model can take degradation and stutters
into account.

## Installation

Installation for Windows: Open R (v4.2 or newer) and write

``` r
install.packages('https://github.com/oyvble/euroformix/releases/download/v4.2.4/euroformix_4.2.4.zip',repos=NULL,type='win.binary')
install.packages(c('numDeriv','gWidgets2tcltk','cubature','XML','curl','plotly'))

#Run EuroForMix (GUI)
library(euroformix);efm() 
```

Also possible to download zip file from
<https://github.com/oyvble/euroformix/releases> (compiled for Windows
only)

Alternative installation directly from source through GitHub (requires
R-tools):

``` r
install.packages("remotes")
remotes::install_github("oyvble/euroformix")
```

Info about suggested R-packages (all are optional to install):  
- gWidgets2tcltk: Used for R GUI  
- plotly: Used to vizualiaze DNA-profiles  
- cubature: Used for numerical calculation of Bayes Factor  
- XML/curl: Used to read allele frequencies from Strider population
databases

## Open the GUI with data from manual (saved project)

This will open the EuroForMix GUI with a tutorial project:

``` r
library(euroformix) #load package
pkg = path.package("euroformix") #get package install folder
proj = paste0(pkg,"/tutorialdata/EFM4proj_imported.Rdata") #obtain project file
library(euroformix);efm(proj)
```

## Command line tutorial

### Part 1: Calculating MLE based LR

#### Step 1: Import and visualize profiles

``` r
library(euroformix) #load package
pkg = path.package("euroformix") #get package install folder

kit = "ESX17" #defining kit to use (must be defined in getKit())
AT = 50 #analytical threshold used (global for all markers)

#Importing allele frequencies
freqFile = paste0(pkg,"/FreqDatabases/",kit,"_Norway.csv") #frequency file to use
popFreq =  freqImport(freqFile)[[1]] #need to select 1st population


#Importing evidence and reference profiles:
evidfn = paste0(pkg,"/examples/",kit,"_3p.csv")
reffn = paste0(pkg,"/examples/",kit,"_refs.csv")
evidData = sample_tableToList(tableReader(evidfn))
refData = sample_tableToList(tableReader(reffn))
```

``` r
plotEPG2(evidData,kit,refData) #Show in graphical interface
```

#### Step 2: Specify hypotheses for interpretation

Hypothesis sets: Ref3 as person of interest (POI)  
Hp: Ref1 + ref3 + 1 unknown  
Hd: Ref2 + 2 unknowns (all unrelated)

``` r
#Set up hypothesis (contributors)
POIidx = 3 #index of POI (in refData)
#Must construct a 'contribution vector' for each hypothesis:
condHp = c(1,0,2) #C1=Ref1, C2=Ref3
condHd = c(1,0,0) #C1=Ref1
knownRefhp = NULL #No known non-contributor reference under Hp
knownRefhd = POIidx #known non-contributor reference under Hd
NOC = 3 #assumed number of contributors
```

#### Step 3: Model fit of Hp and Hd

``` r
#We keep degradation and back-stutter models on (default), but turns off forward stutter model:
mleHp = calcMLE(NOC,evidData,popFreq,refData, condHp, knownRefhp, kit, FWS=FALSE) 
mleHd = calcMLE(NOC,evidData,popFreq,refData, condHd, knownRefhd, kit, FWS=FALSE) 
```

#### Step 4: Obtain calculated LR based on MLE

``` r
MLEresult = calcLRmle(mleHp,mleHd)
LRmle = MLEresult$log10LR #get LR on log10 scale
LRmleMarkers = MLEresult$log10LRmarker #get LR per markers
upperLR = MLEresult$log10LRupper #get theoretical upper LR
```

#### Step 5: Perform model validation

``` r
validhp = validMLEmodel(mleHp,"Hp")
validhd = validMLEmodel(mleHd,"Hd")
nSignifHp = sum(validhp$Significant) #numbers outside envelope
nSignifHd = sum(validhp$Significant) #numbers outside envelope
```

#### Step 6: Provide deconvolution

``` r
DCtableHp = deconvolve(mleHp)$table2 #top ranked genotypes (Hp)
DCtableHd = deconvolve(mleHd)$table2 #top ranked genotypes (Hd)
```

#### Step 7: Show model fit (expectation vs observations)

``` r
plotTopEPG2(mleHp)
```

#### Step 8: Provide non-contributor simulations:

``` r
nTippets = 10 #this typically takes a while (depending on the number)
tippets = calcTippet(POIidx,mleHp,mleHd,nTippets,seed = 1234) 
```

### Part 2: Calculating Bayesian based LR

#### Step 9: Calculate conservative LR (and estimate Bayes Factor)

``` r
#obtain 10% quantile as 'conservative LR'
mcmc = calcLRmcmc(mleHp,mleHd, 5000,quantile = 0.10,seed=1234)
LRcons_mcmc = mcmc$log10LRcons #estimated conservative LR
LRbayes_mcmc = mcmc$log10LRbayes #estimated Bayes factor (LR) using MCMC

#Possible to increase the number of samples (extending the MCMC chain):
mcmcExtended = calcLRmcmc(mleHp,mleHd, 5000,quantile = 0.10,seed=1234,,mcmcObjList=mcmc$mcmcObj)
```

#### Step 10: Calculate Bayes Factor with numerical integration

``` r
#The following calculation would take some time
int = calcLRint(mleHp,mleHd)

#Obtaining calculated LR with relative errors:
LRbayes_int = int$log10LR 
LRbayes_intError = int$log10LRerror
```
