#include <RcppArmadillo.h> //necessary for functions find, intercept
using namespace arma;

extern uvec getContributionIndices(int genotypeCombinationIndex, int numberOfContributorIterations,  Mat<uword> outG1allele, int alleleIndexOfDropouts, 
	int numStutterModels,int nPotentialStutter, int *BWto, int *FWto, int numberOfGenotypes, int possibleContributionsPerContributor, 
	uvec contributionIndices, Col<int> contributorPower );
	