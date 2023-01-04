//WRITTEN BY OYVIND BLEKA (2022). ALGORITHM PROVIDED BY JERRY HOOGENBOOM
//IMPLEMENTATION OF THE IMPORTANT getContributionIndices function

#include <RcppArmadillo.h> //necessary for functions find, intercept
using namespace arma;

uvec getContributionIndices(int genotypeCombinationIndex, int numberOfContributorIterations,  Mat<uword> outG1allele, int alleleIndexOfDropouts, 
	int numStutterModels, int *BWto, int *FWto, int numberOfGenotypes, uvec contributionIndices, Col<int> contributorPower ) {
	
	//Rcpp::Rcout << outG1allele << "-----------\n";  
	int genoIdx; //used to indicate genotype index
    uword alleleIndex1; //used as help variable
    bool hasNumStutterModels1 = numStutterModels > 0; //save as variable since used twice
    bool hasNumStutterModels2 = numStutterModels > 1; //save as variable since used twice
    uword alleleIndexOfDropouts0 = alleleIndexOfDropouts; //convert to uword
  
    int genotypePower = 1;
	for(int contributorIndex=0; contributorIndex < numberOfContributorIterations; contributorIndex++) { //loop over all contributors		
		genoIdx =  (genotypeCombinationIndex/genotypePower) % numberOfGenotypes; // //get genotype for partifular contributor
		urowvec alleleIndices = outG1allele.row(genoIdx);  //alleleIndices are index for alleles in a genotype
		//Rcpp::Rcout << alleleIndices << "\n";
		
		//THIS IS INNER FUNCTION (updateContributionCombinationIndices)
		for(int aind=0; aind<2; aind++) {			
		  alleleIndex1 = alleleIndices(aind);
		  contributionIndices[alleleIndex1] += contributorPower[contributorIndex];
		  if(hasNumStutterModels1 && (alleleIndex1 != alleleIndexOfDropouts0) ) {
			if(hasNumStutterModels2) { // Performance trick: compiler combines >=1 and >1 into a single comparison
			  contributionIndices[ FWto[alleleIndex1] ] += 6 * contributorPower[contributorIndex];
			}
			contributionIndices[ BWto[alleleIndex1] ] +=  3 * contributorPower[contributorIndex];
		  }
		}
		if (hasNumStutterModels1 && alleleIndices(0)==alleleIndices(1)) { // Need to make a small correction for homozygotes...
		  if (alleleIndex1 != alleleIndexOfDropouts0) {
  			contributionIndices[ BWto[alleleIndex1] ] -= contributorPower[contributorIndex];
			if (hasNumStutterModels2) {
				contributionIndices[ FWto[alleleIndex1] ] -= 4 * contributorPower[contributorIndex];
			}
		  }
		} //END INNER FUNCTION
		genotypePower *= numberOfGenotypes;
	 }
  
  return contributionIndices;	
}
