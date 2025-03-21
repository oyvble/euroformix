#include <RcppArmadillo.h> //necessary for functions find, intercept
using namespace arma;

extern uvec getContributionIndices(size_t genotypeCombinationIndex, int numberOfContributorIterations,  Mat<uword> outG1allele, int alleleIndexOfDropouts, 
	int numStutterModels, int *BWto, int *FWto, size_t numberOfGenotypes, uvec contributionIndices, Col<int> contributorPower );
	
extern uvec getContributionIndicesAllele(uword alleleIndex1, int alleleIndexOfDropouts, int numStutterModels, int *BWto, int *FWto, uvec contributionIndices, Col<int> contributorPower );
