SUGGESTIONS:
 - Also perform model search under Hp (alternative choice).
 - Adding new functionality: EFM remembers all used settings when closed (as for CaseSolver). 
   Following could be stored: Options in toolbar,  Selected Kit, Selected population frequencies, more??
 - For PHexp plots (plotTopEPG2):
	- Add allele frequency to hover labels (easy to hightlight)
	- Add LR value of each marker as text label (under marker name, or after)

KNOWN ISSUES (TROUBLESHOOT):
 - Upgrade to latest Rcpp package (1.0.9) to fix issue "function 'Rcpp_precious_remove' not provided by package 'Rcpp'". Notice that this package must be installed and loaded non-virtually (in order to work with Citrix)

Future version:
 - Add the possibility to define a prior for the PHvar param.
 - Make a scroller for generateData panel to enable the possibility to generate data with many markers.

- Faster Non-contributor tests by optimizing the calcMLE function (update genoWeights, skip preSearch, modelValid and DCcalc)


EuroForMix v4.2.3 (Release date: 2025-02-07)
=============================================
 - Improved how LR is presented when the model under Hp or/and Hd is not fitting the data.
   - Always showing two decimal points for log10LR.
   - In calcMLE-L321: Set logLik=-Inf if logLik reached .Machine$double.xmax.
   - Modifying .getSmallNumber to include special situations (LR=1, LR=0 and LR=NaN).
   
 - Fixing crash in calcMLE if all start values gave logLik=-Inf values (happens if model does not fit the data):
   - calcMLE-L163: Adding "if(length(maxL)==0) maxL = -Inf"
   - calcMLE-L172: Adding !is.na(valdiff) 
   - calcMLE-L182: Use seq_len(nrow(bestPreSearchParams))
    
EuroForMix v4.2.2 (Release date: 2025-02-06)
=============================================
 - Fixing crash in calcMLE when stutter-proportion was very close to zero and required a impute
   - Caused by wrong maxPhi impute when it was infinite (need to separate on Mx and stutter prop for the impute).

 - Fixing bugs in contLikSearch function:
   - Component of POI was not correct when extracting the Mx (could cause wrong Mx values).
   - At L56 where knownNonContrHp was set as NULL even if not supposed to. 
	 - This caused differences against ordinary calculations when including references that are not conditioned on.
   - Allowing the model search to go up to 6 contributors.
   - resttol is now also included as argument to function.

 - Fixing an issue in efm when having triAlleles of a reference:
	- The reference was deactivated in the data selection window (this is now fixed by allowing at least 2 alleles: efm:L1738).

EuroForMix v4.2.1 (Release date: 2025-01-24)
=============================================
 - Fixing an issue when the model is not converging and throws an error for the model validation .
	- Fixed in validMLEmodel (L41-L43) by allowing NaN values (forced as 0).
	- Also causing the deconvolution to fail (fixed in deconvolve:L132)
	- The issue was triggered when doing Automatic model search when setting a too low NOC and having replicates.

EuroForMix v4.2.0 (Release date: 2025-01-10)
=============================================
 - Adding the possibility to condition on references with tri-alleles:
	- Storing alleles which are exceeding allele index 2 is stored in a triAllele vector (prepareC:L321-L328)
		- Structuring the information of triAlleles in a vector with three elements: markerIdx,alleleIdx,contributorIdx.
	- Modification in EFMfastengine.cpp and getContributionIndices.cpp
		- Structuring triAlles to all markers (L790-L812),
		- Adding allele-contribution for triAlleles for known references (L172-L179). Utilizes a new function: getContributionIndicesAllele.		
	- Other modifications to support triAlleles:
		- In prepareC-L319: Removing restriction that "References can't have more than two alleles!"
		- In efm:L1108: Removing restriction that imported references cannot have more than two alleles (comment out). 
		- In deconvolve:L96-L114: Adding triAlleles to top genotype of known reference 
		- In efm:L2579-L2603, updating helpfunction f_saveDCasRef.
	- Note: 
		- The triAllele will not provide any stutter product when the forward stutter model is assumed (technical limiation).
		- The genotype including triAllele(s) cannot have more than 2 equal alleles (not possible for references with "homozygous tri-alleles")
		- It is not possible condition on references with triAlleles when calculating Qual.LR or generating new samples.
	- Added tests: Including 1-4 triAlleles (two contributors). Also shared.
 - In the Joint LR result of EFM GUI: Also including the non-logged version of the LR 
 - Fixing the issue that very small LR results causes log10LR=-Inf:
	- Updated calcLRmle function by utilizing the .getSmallNumber helpfunction
	- In efm:L2288: Adding log10LRmle from calcLRmle function	
 - Fixed color of markers in plotTopEpg2 (L90+L158): Did not obtain marker specific genotypeProbs
 - Changing maximum number of contributors to specify in a hypothesis to 5 instead of 4 (efm-L36)

EuroForMix v4.1.1 (Release date: 2024-11-18)
=============================================
 - In calcMLE function (L230-L244):
	- More robust handling when Mx-solution is on zero-boundary (phiMax modified). Avoids crash.

EuroForMix v4.1.0 (Release date: 2024-11-05)
=============================================
 - A new presearch algorithm is carried out before MLE optimization (and INT and MCMC calculations).
	- The presearch spans a set of mixture proportion parameter vectors and calculates the likelihood for each outcome.
	- The top "nDone" (req. optimizations) sets from the presearch is used as starting points for the optimization. 
		- This causes it to avoid random startpoints and hence "random results". Seed will not have any effect.
	- Notes: 
		- Less than "nDone" optimizations can happen when no other likely startpoints are found.
		- The optimization is now more robust (more likely to find the global solution).
	
 - Faster likelihood calculations: 
	- A new implementation restricts the genotype outcome based on information from the presearch 
	- The restriction criterion is based on a threshold (default=1e-6, selected through development validation).
		- The threshold can be changed under Optimization->"Set threshold for genotype restriction".
		- Threshold value set by user added under "Optimalisation setting" in the Report text file.
		- A warning is given in the R console if it is recommended to lower the threshold value.
		- Threshold value=0 means no restriction (full genotype outcome used).

 - Faster and more effiecient deconvolution and model validation calculations:
	- Carried out in EFMfastengine C++ code utilizing the restriction and the v4 algorithm.
	- The algorithm utilizes full parallelization thanks to "declaring the reduction" into vectors.
	- Note: Joint Deconvolved profiles no longer accessible. List elements "Table 1" and "rankGi" from deconvolve is now empty.

 - Modified MCMC method:
	- Restricts the mixture proportions for the unknown (unrelated) contributors.
	- A parameter will be assumed as known (fixated) if its corresponding StdErr<0.001, and will not be part of the MCMC. 
		- These are typically parameters estimated close to zero (boundary estimates).
 
 - Modified integral method (Bayes Factor) where the limits of the mixture proportion integrals are correctly defined.
	- The defined limit for the integration was set as a deviation of 2 instead of 3 (too wide earlier).
	- Scale returned from getParamLimits is now based on the maximum likelihood value (and not the absolute value of this)
	- Integration of a stutter proportion (and mixture proportion) variable will not be carried out if the width of the integral is less than 0.01.
	- NOTE: Integral method is not always calculated correctly when defining relationship in the hypotheses (depends on size of component).

 - Added functions:
    - exportDataFromProject: Exporting evidence and reference profiles and frequency data in project to EFM compatible files.
	- freqs_listToTable: Convert frequencies in list format to a freqency table.
	- prepareCobj: Wrapper function to create and prepare the C++ object (fill data and prepare structure)
	- Helpfunctions (hidden):
		- .preSearch: Performs the presearch of Mixture proportion (Mx) set.
		- .getMxOutcome: Obtaining the outcome of Mx proportion sets to traverse in the presearch.
		- .getMxValid: Make sure that the Mx vector is correct order (restricting the unknown unrelated).
		- .checkMxValid: Check that the Mx vector is correct order (restricting the unknown unrelated).
		- .getFittedParams: Prepares output of fitted parameters in MLE (after optimization).
		- .getDataToPlotProfile: Obtaining data for showing model contribution in "plotTop" plots.

 - Replaced functions due to being outdated (now requires plotly installed):
	- plotTopEPG now calls plotTopEPG2
	- plotTopLUS now calls plotTopMPS

 - Minor notes:
	- URL link to STRidER was changed to https://strider.online/frequencies/xml.
	- Degradation and stutter parameters are now log-transformed instead of logit transformed in the optimization domain.
		- The degradation parameter can now have estimated values greater than one.
		- The stutter proportion parameters was restricted to not exceeding one.
	- The test environment now calls the C++ function "calcGenoWeightsMax" instead of "loglik" because loglik requires restriction function to be called first.	
	- Adding a test for the integral expression for mixture proportions (test_integrals.R).
	- In calcMLE function: Arguments Seed and delta will no longer have any effect since a pre-search used instead of random start points.
	- Removing a few options in toolbar of the EFM GUI which are no longer in use.
	- In efm function: Adding version to a project. Also printing the version of a project when loading it (if found).
	- In create report (efm_createReport.R):
		- Time of created report now given in whole seconds instead of decimals of seconds.
		- Get only version number of R when printing out the R-version.
	
 - Fixed bug: Only for v4.0 (not earlier).
	- AMEL was wrongly calculated when using a stutter model: Fixed in (EFMfastengine.cpp) by applying NumStutterModelsMAX to m_possibleContributionsPerContributor.

EuroForMix v4.0.9 (Release date: 2024-05-24)
=============================================
 - Changed STRIDER-path to "https://strider.online/frequencies/xml"
 - For relationship calculations:
    - Require that all markers are typed for a related reference (throws a warning if not)
	- Introduction of new relationship model: If ibd is given but no related reference index is given, the last two unknowns will be assumed specified relationship.

EuroForMix v4.0.8 (Release date: 2023-08-21)
=============================================
 - Fixed bug with non-contributor test, occuring when having rare alleles of known non-contributors under Hd.
	- Added lines in calcTippet:L75-L82.
 - Improved README file: Easier information to get started. Inspired by work of https://github.com/magnusdv/forrel

EuroForMix v4.0.7 (Release date: 2023-05-30)
=============================================
 A major issue caused by new R-version (v4.3.0) was fixed (caused program to crash):
 - Fixed bug when trying to specify marker specific settings. At efm:L325; "!is.na" was removed.
 - Also modifying getSampleType:L15.

 Minor changes:
 - Removing printout when using calcLRmcmc with verbose=FALSE argument.
 - Fixed small bug in plotLUS:L26. 
 - Tidying up NEWS.


EuroForMix v4.0.6 (Release date: 2023-04-28)
=============================================
 - Fixed bug when not having full reference profile of related individual (prepareC:L377). Caused program to crash.
 - Fixed bug when calculating non-contributor tests for related individuals with partial profile (tidying up calcGjoint function).
 - An error message is given to explain a crash when following is happening: 
	- Evidence profile covers all alleles in frequency database and a conditional reference has a rare allele (not in database).


EuroForMix v4.0.5 (Release date: 2023-03-27)
=============================================
 - Fixed bug when having full marker dropout across all replicates (when evaluating at least two replicates):
	- In prepareC:L216. Drop-in model vector not correctly assigned for Q-allele of replicates. 
	- Consequence: Lead to exclusionary LR and model validation could trigger an error.
 - Fixed "min" typo on fitgammamodel:L54. Also modify "min" to "minimum" at gammamodel:L90 to avoid issues.
 - Removed "DEG" as argument to function prepareC, since not in use. Function calls are adjusted.


EuroForMix v4.0.4 (Release date: 2023-01-19)
=============================================
 - An issue was found when evaluating the relationship of the major of a clear major/minor profile (Thanks to Damir Tesanovic for highlighting this). 
	- The issue causes the optimizer to struggle in reaching global maximum (and causing wrong MCMC results). 
	- Modifying code in helpfunctions.R->.paramrandomizer: Including boolean whether last unknown is related (Mx for last unknown must vary freely).
	- The issue causes the BayesFactor integration to be wrong.
	- Modifying integration limits in calcINT.R: The related now counts as a known contributor (Mx for last unknown must vary freely).
	- An extra variable from prepareC was included (hasKinship) which eases the recognizion.
	- Including additional test in test_logLik2contr2Rep: Unknown is sibling of major contributor.
  
 Minor changes:
 - Avoids error in GUI when forgetting to select related individual: Added check at efm:L1924.
 - sortComb argument in calcGjoint was removed (not used anymore).
 - Error estimate for numerical integration does no longer scale with combination factor (calcINT:L192 is ignored).


EuroForMix v4.0.3 (Release date: 2023-01-05)
=============================================
 - The modification of the EFMfastengine code in previous version (v4.0.2) led to worse speed: The implementation is reversed.


EuroForMix v4.0.2 (Release date: 2023-01-04)
=============================================
 - Slight speed improvement by modifying the EFMfastengine code:
	- Precalculation of weights are now parallelized on both allele and outcome traversion.
	- Removed unused arguments in getContributionIndices C++ function (doesn't affect speed).   
 - Fixed bug when calculating non-contributor LRs for qualitative model (efm:L2796).
 - Fixed bug when calculating model validation for R v3.6.x and earlier (data.frame in validMLEmodel.R caused wrong format)
 - Fixed crash issue with earlier R-versions (v3.5.x) which used earlier versions of RcppArmadillo (uvec vectors not initiated as zeros)


EuroForMix v4.0.1 (Release date: 2022-12-02)
=============================================
 - Fixed bug causing non-contributor panel in GUI to vanish.
 - Fixed issue causing R versions earlier than 4.0 to compile:
	- Removing shared variables for OpenMP operation (EFMfastengine-L220): 


EuroForMix v4.0.0 (Release date: 2022-11-28)
=============================================
 Special thanks to: Victor Saragoni (tester), Damir Tesanovic (tester) and Jerry Hoogenboom (algorithm).

 Major changes:
 - Speedup of qualitative likelihood calculations: calcQual is an alternative implementation of forensim::likEvid (can be ~10x faster)
 - Speedup of quantitative likelihood calculations: MLE, INT, MCMC based calculations utilizes this.
 - Modified underlying C++ code for deconvolution and validMLE (utilizing EFMrep development).
 - All calculation functions are back-compatible: Same arguments are kept.
 - The old functions calls the new ones indirectly: contLikINT->calcINT, contLikMLE->calcMLE

 - Modified fragment length designations:
	- Non-observed alleles now linearly interpolated instead of using closest.
	- New model option for assigning the fragment length of Q-allele: Weighted average of non-observed allele frequencies.
	- Can be turned on in Settings or given as argument in following functions: calcMLE,calcINT,contLikSearch.
  
 - Updated feature for Conservative LR calculation (LR sensitivity): 
	- The calculation can be done accumulatively. In the GUI: Pressing the button again will extended the MCMC chain.
	- A 95% CI is included to take into account the MCMC simulation uncertainty.
	- A traceplot is included to give the user an idea of accuracy as a function of number of iterations.
	- The Bayes Factor estimate based on MCMC is also superimposed in both plots.
    
 - Added functions:
	- Main functions: tableSaver, efm_gfile, efm_DBsearch (database search module is extracted outside of efm function)
	- Helpfunctions (hidden) used for calculations:  .convBack, .paramrandomizer, .calcJacobian, .secondToTimeformat, .printTable, .plotTippet, .getFragLength
	- plotSumPH placed outside efm function (used to show degradation trends)
  
 - New functions to simplified calculations where input is fitted models under Hp and Hd.
	- calcLRmle: Returns MLE based LRs.
	- calcLRmcmc: Based on MCMC. Provides both conservative and estimated full bayesian.
	- calcLRint: First utilizes getParamLimits based on fitted model under Hd to obtain param limits.
	- calcTippet: Non-contributor analysis
  
 - Updated Upper boundary LR: Now takes into account following (results may differ due to this):
	- All typed individuals (under Hp) are included for calculations. 
	- Missing markers of POI should not affect RMP and scaling.

 - In GUI:
	- Created button for showing LR-per markers in results
	- Imported Profiles/Databases can now be removed from the selection list (Data panel).
	- In "MLEfit" panel: Also show which sample(s) and hypotheses that are evaluated.
	- A new button for storing top ranked genotypes as references (after Deconvolution), stores as EFM-format.
	- Integration based LR (Bayes Factor) is no longer possible for database search.
	- Changed options for Integration, MCMC
	- All calculated results are stored in resEVID: MLE,INT,MCMC (and validHp/validHd). SEARCH is stored in setEVID.
	- The created reported is now done in a separate function: createReport (taking mmTK environment as input).
	- Additional comparison values of Ref->Evid added when "View Reference".
	- When importing profiles: Gives user option whether to still import data if OL alleles are observed (automatically removed)
	- The data selection in Model panel is moved to a separate window (activated with button). Also accepting any number of markers
	- In settings: "Max number of loci" is removed. This is no longer a limitation (data selection and LR-per marker modified)
	- Peak height plots: EPG/Bar-plots are shown in Browser only if plotly is installed.
	- The sumPeakHeight visualiziation against fragment length is modified: 
		- A smoothing curve is shown for each replicate (colored). Also no longer showing p-values.
  
 - The calcLikMLE argument "maxIter" is replaced with "difftol": tolerance for being exact in log-likelihood value (relevant when nDone>1).
	- Function noncontrMLE replaced with a new function calcTippet (supports both MLE/INT).
	- Additional information added to the report:
		- Detailed information about samples and references are included into the report.
		- The number of evaluations and time usage for the MLE method is included to the report ("Number of evaluations")

 - Update in plotTopEPG2/plotTopMPS2 plots: Stutter products are now also shown

 - Minor changes: 
	- Tidy up RMP printout for references.
	- Fixed issue of showing MCMC when assuming one contributor.
	- Degradation option can now be included for the model search.
	- Plot functions plotEPG2/plotMPS2
	- plotTopEPG2/plotTopMPS2 now uses information from fitted model (prepareC)
	- Fixed crash when replicates without kit specification is visualized (L1257: Use lapply instead of sapply).
	- Database search LR results now shown on log10 scale
	- Better error handling regarding the markerSpecific settings. 
		- Still error if user forget to use "Save settings".
		- Added guidance text to "setVecRightOrder" function.
 
 
EuroForMix v3.4.1 (Release date: 2022-03-22)
=============================================
- efm function updates:
	- Fixed issue with possible wrong reference name under Hd when hovering Mix-prop parameters.
	- Additional statistics given when hovering Mix-prop parameters: "Average rfu (across markers) multiplied with mixture proportion"
	- Avoid crash when directory of stored project was not found.
	- Restoring MLE results from a stored project when deconvolution module was used.
- plotLUS updated: accepts no adata attribute for refData.
- Added license file: GNU Lesser General Public License v3.0


EuroForMix v3.4.0 (Release date: 2022-02-05)
=============================================
 The STRidER import issue is fixed:
 - Connection is using curl instead of Rcurl to parse text from file directly (package dependency changed).

 Important decision changes:
 - In efm:L2140: References which are put forward to the Hypothesis specification GUI, but not included as contributors, are included as typed non-contributors. 
 - Model validation calculations slightly modified to integrate from 'AT-1' instead of 'AT'; calcloglik.cpp:L(782,785,793,795).
	- The reason is that peak heights are discrete and not continuous.
	- This avoids zero pvalue when observed peak=AT .
 - The unknown related individual under Hd is now last contributor.

 Minor updates:
 - Avoiding hardcoded flags in Makevars file enables better system compatibility for Linux/MAC (thanks to Simon Urbanek for this suggestion).
 - The numerical test code has been tidied up (testthat).
 - Fixed bug in plotMPS2: Marker specific AT/ST argument was not used if provided.
 - prepareData now does no longer require that alleles of references are stored in element name "adata".
 - The fitgammamodel function was updated: Make pre-fitting of optimization start points even better and more robust (based on linear regression):
	- The MLE approach calculations may give slightly different values, depending on the steptol parameter.
 - Fixed a formula error for kinship calculations when Fst>0 and at least two unknowns (thanks to Jerry Hoogenboom for highlighting this).
	- In calcloglik.cpp: The relatedness genotype probability is now calculated last (OK when having one relationship).
	- calcGjoint function was updated.
	- Several kinship hypotheses was added to test_logLik2contr1RepNoStutter.

EuroForMix v3.3.3 (Release date: 2021-10-14)
=============================================
 Faster computations (~10% increase).
 - Now using Rmath::pgamma instead of boost::gamma_p function in the C-code (the BH R-package is no longer required).


EuroForMix v3.3.2 (Release date: 2021-10-06)
=============================================
 Fixing Import from STRidER issue giving "SSL certificate problem.." crash:
 - The getURL call in freqImport function was updated.

 Fixing a whitespace issue which occured for some textfiles:
 - The sample_tableToList now also ignores " " alleles.

 Minor update:
 - Included MiniFiler as a kit (Thanks to Oskar Hansson to providing this).


EuroForMix v3.3.1 (Release date: 2021-05-26)
=============================================
 Important update:
 - The URL to the STRidER connection can be customized in toolbar (Frequencies->Set URL for STRidER import).
	- The default URL path has been changed (earlier path did not work). Thanks to Marcin Wozniak for making aware of this issue.

 Minor update:
 - The Lik value in report (Create report) now avoids underflow. A helpfunction "getSmallNumber" was created.
 - Suppressing gWidgets2 warnings in the efm function (suppressWarnings wrapped around the GUI code). 

 Bug fixes:
 - In efm function (L1209,L1396): gWidgets2::gmessage gave error because of wrong arguments (messages are put first).


EuroForMix v3.3.0 (Release date: 2021-03-30)
=============================================
 Important updates:
 - The mix-proportion parameter of the last contributor (C_K) is now shown in the MCMC-results (function validMCMC updated).
 - Wrong calculations were obtained when the alleles in the evidence profile covers the ones defined in the frequency data, for at least one of the markers, and at the same time a stutter-model is considered (BW or/and FW).  
	- This bug was introduced in v3.0.0 and remains until v3.2.0. See "Important Update" section on the website and "bug fixes" below for more details. Thanks to Kevin Cheng, John Buckleton and Jo-Anne Bright for providing the example leading to this discovery.
 
 Minor updates:
 - The Optimal model search now also uses the seed from Optimization setting. efm:L2303.
 - Function genDataset (L50) now throws a warning instead of printing when the frequencies is "not a valid simplex".
 - In the example code for some of the external functions: now uses the right path to the population frequency file in the installation folder.
 - The insertion of rare alleles are now printed to the console.
 - The function contLikSearch now also contains "steptol" as argument (obtaining identical results as ordinary optimization). Thanks to Damir Tesanovic for noticing this.

 Bug fixes:
 - A bug was found when not having Q-allele (evidence alleles fully overlaps with frequency alleles) and considering stutter model:
	- prepareC:L244-L245; The dummy index "-1" should not be included. 
	- Adding 2 blocks to src/calcloglik:(L350-357 and L777-784); Last allele is also traversed if not a Q-allele. 	
 - If the reference had a dropout allele not defined in population frequency table and the observed evidence fully overlapped with the frequency table, the program earlier crashed.
	- The program now force the insertion of a Q-allele with minimum frequency to the frequency table (Qassignate:L77-82).
	- Care must be taken if this happens only under Hp (then there will be allele frequency missmatch under Hp and Hd).
 - Not using encoding when working with gfile caused crash for non-UTF8 characters (resolved by applying a mygfile helpfunction)
 - Setting empty Seed after putting number in "Set seed of randomizer" is now possible.

 Added tests:
 - test_logLik1contr_noQ: Example where population frequency table is same outcome as evidence profile for one marker.
 - test_logLik4contr1Rep_noQ: Same as above but a larger example considering a 4-person mixture (frequency of D16 has same outcome as evidence profile).
 - test_logLik4contr1Rep_noQref: Same example as above but we force the reference to have a drop-out allele for the D16 marker (causing inclusion of Q-allele).


EuroForMix v3.2.0 (Release date: 2021-02-17)
=============================================
 Important updates:
 - GUI R-package changed to gWidgets2 instead of gWidgets (gWidgets was removed from CRAN). Remember to install gWidgets2 and gWidgets2tcltk if not already done.
	- Helptext may be available by hovering the mouse over button or text.
	- The user can now hover the mix-proportion parameter estimates to highlight the name of the corresponding contributor.
	- The list of sample profiles are now stored in a table (can list a large number).
 - Importing population frequency data from selected folder will only use a folder from the euroformix installation path (~euroformix/FreqDatabases)
 - The number of failed PH-validation points are now included to the exported report (this is calculated when exported).
 - The MCMC simulations are automatically calibrated by iteratively running 100 samples under Hp (fastest). 
 - Function contLikMCMC now return logMargL instead of margL to avoid underflow (log-scale).
 - Modified MCMC results in report: 
	- The report now also contains the estimated marginalized LR (Bayesian based on MCMC)
	- The variation of the randomizer parameter is now also mentioned (needed for full reproducibility).
 - The order of contributors under Hp is now same for both automatic model search and sole LR calculations (earlier the POI was always put last).
 - The calculation of the upper boundary LR has now been modified when there are at least 2 unknown contributors under Hd and fst>0: The random match probability (RMP) is scaled (with value greater or equal than 1) with a formula depending only on selected fst value. Notice: The LR of a POI fitting a clear major unknown components could exceed 1/RMP(POI) when fst>0 when assuming at least 2 unknowns under Hd.
 - It is no longer required to install external R-packages in order to load euroformix. This would be useful for those who only wants to use CaseSolver.

 Minor updates:
 - Changed name in "File->Setup" removing "Default" word.
 - Transparant level on barcolors in plotTopEPG2 has been increased (better visualization).
 - When importing data from text files: A previous printout to R console is no longer performed.
 - The degradation model is no longer a choice for automatic model selection (due to new GUI implementation).
 - A message is shown when trying to select "View evidence/reference/database" without any selection.
 - In genDataset function: Frequencies are rescaled to sum to 1 instead of throwing the "not a simplex" error.
 - There are no longer any restrictions of number of markers to visualize in "Marker settings".

 Bug fixes:
 - Marker settings are reset when saving "global settings" (File->Setup). Avoids earlier crash when using only global settings after Marker settings are saved.
 - Fixed problem in plotTopMPS2/plotTopEPS2 which caused marker names to not be shown when assuming 1 contributor.
 - Handle that references may have zero alleles in a marker (fixed bug in "View References" and made inactive in "Data selection")


EuroForMix v3.1.1 (Release date: 2020-09-08)
=============================================
 Bug fix:
 - Problem loading projects and seting working directory when path involves non UTF8-characters.


EuroForMix v3.1.0 (Release date: 2020-09-08)
=============================================
 Important updates:
 - In GUI: It is now possible to export a frequency file from the 'Import data' panel (based on "Selected population").
 - In GUI and contLikMLE function: The user can now change the steptol parameter used in the optimization. This makes it possible to increase accuracy of the results (steptol=0.001 is default from v3, whereas it was 1e-6 in earlier versions).
 - MLE Optimizer strategy: Default number of optimizations set to 3 (instead of 2). This will make the MLE more robust. The user can change this at any time.
 - Default drop-in probability is now "0.05" instead of "0".
 - A new R-function called "qualLikMLE" has been included to the package, enabling more flexibility for inferring the qualitative model: Different drop-out probability for each contributor, in combination with fixed values, is possible. NB: This is not yet implemented in the GUI.
 - Enabling degradation model for MPS-STR sequences with "RU:sequence" format, when kit is selected (only ForenSeq included): Example is "10:[ATCG]10", where "10" is extracted for fragment length lookup.
 - Adding an error message if data is not properly structured before running contLikMLE: The function prepareC now throws the error "Please use the prepareData function to structure the data correctly." if last element.

 Bug fixes:
 - In function deconvolution: The startMarker indices was not correctly updated when no unknowns was given. This earlier caused crash.
 - In function contLikMLE-L349: "maxPhi" replaced with "phi" (only syntax change).
 - In function contLikSearch-L86: knownRel replaced with NULL. This bug caused the Optimal model search to fail when conditioning on known references under Hp/Hd and having relatedness under Hd. Thanks to Kent Harman for disovering this.

 Added tests:
 - test_logLik2contr1RepSNP -> Numerical testing with 1 and 2 unknowns for 1 replicate (ensures that data with SNP frequency data becomes correct).
 - test_logLik2contr1RepNoStutter -> Numerical testing with 2 unknowns for 1 replicate (ensures that C++ function calcLogLikGammaMarkerUnknown is tested).
 - test_QualLogLik2contr2Rep -> Numerical testing for drop-out probability inference with qualitative model, with different variants.
 - test_logLikcontr2RepSNP renamed as test_logLik2contr2RepSNP.

 Minor updates:
 - Added error message in the prepareC function (L185-L189) if there are alleles in the frequency outcome which was not observed (except "99" allele): "In marker ....: There are alleles in the frequencies data are not found in the observed data. This is not allowed. Please use the prepareData function beforehand! Program stops.". Thanks to Guro Dorum for providing an example with this.
 - The "Optimal quantitative LR (automatic model search)" module can only be run if hypothesis Hp and Hd only differ if a POI under Hp is replaced with an unknown under Hd.
 - Better error print out when allele frequencies does not sum to one (prepareC:L65-L70 and calcGjoint:L31-L32).
 - Added error message in function getData if any marker names between the data types are inconsistent.
 - Small text changes in report and GUI regarding optimization settings.


EuroForMix v3.0.3 (Release date: 2020-05-06)
=============================================
 Important bug fix in the fitgammamodel function (bug introduced in v3.0.0):
 - The returned degradation parameter value was not correct (always overestimated). This influences the degradation slope plot in the EFM GUI. Also, the fix may cause faster convergence of MLE. 

 Other fixed issue: 
 - GUI crashed when loading settings after saving settings with old projects (created with a version before v3.0.0). Thanks to Peter Gill for pointing out this issue. 


EuroForMix v3.0.2 (Release date: 2020-05-06)
=============================================
 Issue fix:
 - The MCMC simulations failed if first random startpoint gave a likelihood of zero. Thanks to Eugenio Alladio and Maria Martin Agudo for reporting this issue.


EuroForMix v3.0.1 (Release date: 2020-04-24)
=============================================
 Bug fix:
 - Markers containing space caused problem when loading Marker settings after saving (reading from the configMarkerSetup file failed).


EuroForMix v3.0.0 (Release date: 2020-04-01)
=============================================
 Summary of version 3 updates:
 - Faster calculations due to paralellization on deeper level; causes speed to scale with the number of logical CPU cores (threads).
 - Progress bar for all types of calculations. Time estimates given for the MLE method (if time exceeds 10s).
 - Automatic model selection for the MLE method (number of contributors, Degradation (ON/OFF), BW-stutter (ON/OFF), FW-stutter (ON/OFF) decided based on the AkaikeInformationCriterion)
 - Model extensions: Forward stutter model and Marker/Dye specific settings for analytical thresholds, dropin model and fst.
 - Validation plots created (almost) instantly
 - Seed settings: The MCMC and conservative LR values are now reproducible. Also optional seed can be given for the MLE method (reproducible start values, useful for time comparisons).
 - Upper boundary of Likelihood Ratio provided in GUI (1/'Random match probability').

 Details of the update:
 - GUI changes (efm):
	- In 'Generate data' panel: New layout. Includes Forward stutter. 
	- In 'Model specification' panel: Possibility of forward stutter model added.
	- In 'MLE fit' panel: New layout for better per-marker space. 
		- 'Only LR results' removed. 
		- adj.loglik added to loglik: adj.loglik = loglik - nparam
	- In create report this was added: 1) Frequency for rare alleles and 2) information of whether the frequencies were normalized after imputing rare alleles.
	
 - Bug fixes:
	- In function deconvolve: The genotype of known contributors were not taken into account when doing deconvolution with theta-correction (only known non-contributors).
	- In function plotTopMPS2: Inserting base pair of closest allele (for missing bp)

 - New implementation of the C++ engine (Armadillo/Rcpp not used any more):
	 - The C++ engine performs parallelization on genotype marginalisation over one for-loop which makes the speed scale with the number of CPU cores.

 - Model extensions:
	 - Forward stutter model.
	 - Marker dependent AT/pC/lambda/Fst parameters (also implemented for qualitative model).

 - New features:
 	- Monitoring progress (a progressbar is provided in console)
	- Automatic model selection (based on the AIC): Number of contributors, Degradation (ON/OFF), BW-stutter (ON/OFF), FW-stutter (ON/OFF)
 	- Seeds added to following functions: contLikMCMC,contLikMLE 
	- Optional 'normalization' handling of rare alleles (missing alleles in population frequency data): 
		- Can be changed under "Frequencies -> Set whether to normalize frequencies" (1=Yes, 0=No) 
		- Allele frequencies are by default normalized after including new alleles. In EFM v2 the default was not to normalize - set option equal 0=No for concordant results.

 - Following functions has been modified and updated because of the new C++ implementation:
	- prepareC
	- contLikMLE
	- contLikINT
	- contLikMCMC (depends on contLikMLE output)
	- validMLEmodel (depends on contLikMLE output)
	- deconvolve (depends on contLikMLE output)
	- logLiki (depends on contLikMLE output)

 - Following functions are not used anymore (but still there): contLikMCMCpara (calls contLikMCMC), contLikMLEpara (calls contLikMLE), iszerolik

 - Added functions: 
	- runexample: Illustration of using the most important methods in euroformix.
	- calcRMP: Calculates the random match probability of a specific profile.
	- conLikSearch: Searches through all model outcomes and shows the optimal model in the MLE panel. Table is visualised and also stored in report file.
	- sample_listToTable: Converts list structure to table (used by exporting)
	- sample_tableToList: Converts tables to list (used by importing)
	- prepareData: preparing data with Qassignation and filter data wrt analytical thresholds (possibly per marker)

 - Modified functions:
	- genDataset: Also includes forward stutter simulation and marker based analytical threshold
	- getKit: The user can specify other kit files than the one by default (more flexibility)
	- simDOdistr: The function now accepts dropin probability per marker (if provided)
 
 - Includes automated testing of the numerical calculations by C++ using 'testthat' (compares with manual R implementations)
	- See separete presentation "EFM3numtestPres" for a description
 
 - Other notes:
	- Degradation slope param can never exceed 1. This was possible before 
	- New optimization strategy for MLE (contLikMLE): 
		- New strategy: Require nDone number of obtained local maximum values (set to 2 as default).
		- Random start points are better chosen: Mixture proportions (flat, but sorted for unknowns), PHexp/PHvar/stuttprop/degradslope have lower variation.
		- Faster convergence in nlm, but about same accuracy (minimum step tolerance (steptol) set to 1e-3 instead of 1e-6).
	- New sampling strategy for MCMC: No block-wise sampling provided and only one chain is run.
	- The parallel versions contLikMLEpara and contLikMCMCpara are no longer used (calls the original functions)
	- The likelihood of underflow numbers are now shown in MLE results.
	- CaseSolver v1.7 or newer should be used with this version.

