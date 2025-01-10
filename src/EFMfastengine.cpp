//WRITTEN BY OYVIND BLEKA (2022). ALGORITHM INSPIRED BY JERRY HOOGENBOOM
//This script contains the fast implementation of EFM (new algorithm).
//The script was updated in June 2023 to enable restricted calculations.
//Used for following calculations: MLE (optimization), MCMC (simulations), INT (numerical integral)

#include <vector> //vector storage
#include <thread> //used to obtain number of logical processes
#include <RcppArmadillo.h> //necessary for functions find, intercept
#include <Rmath.h> //includes pgamma,qgamma (note that this has namespace 'R')
#ifdef _OPENMP
#include <omp.h> //parallelization
#endif
#include "getContributionIndices.h" //necessary for function getContributionIndices
#include "genoProbFuns.h" //helpfunction for genotype probs
using namespace arma;
using namespace std;

//ORGANISATION:  ExposedClass -> EFMclass -> EFMmarker
/* RCPP_MODULE is used to access internal C++ functions through R (done through the "ExposedClass"). Makes it possible to access stored objects in C++ memory.
//Step 1: filldata - used to initiate data objects (not depending on the NOC=number of contributors)

//Step 2: prepare - used to prepare variables for the likelihood function call that doesn't change for new parameters: Drop-in weights and genotype probs.

//Step 3 (optional): calcGenoWeightsMax - Update genotype weight vector for given params. Consists of sub-steps A and B:
//	Substep A: calculation of the Likelihood matrix EvidWeights (carried out across all data).
//	Substep B: calculation of the total sum by traversing all genotypes (carried out for each marker).

//Step 4 (optional): restrictGenos - Construct vector of genotypes to traverse (may give a restricted set, depending on the input threshold)

//Step 5: loglik - used to calculate the likelihood function for a given parameter set. Consists of sub-steps A and B (traversing only restricted set).

//Step 6: logliki - Returning likelihood value per marker
*/

class EFMmarker { //Each marker is treated separately
	private: 
	int m_NumGenos1p; //number of genotypes for 1 contributor 
	umat m_outG1allele; //allele index per genotype outcome	

	//DATA VARIABLES
	int m_NumPotStutters; //number of unobserved potential stutters (BW+FW)
	int m_NumAlleles; //total number of alleles (both alleles and potential stutters)	
	int m_NumAllelesTot; //total number of alleles (both alleles and potential stutters)
	int m_NumStutterModels; //model type (0=No stutt, 1=BW, 2=BWFW)

	mat m_peaks; //Intensity vector (copied number)
	mat m_dropinweight; //Dropin weights vector (pC*freq*dexp(Y-AT,lambda)	logged
	vector<double> m_freqs; //Frequency vector
	vector<double> m_NumTypedAllele; //counting number of typed alleles of type a
	double m_NumTypedTot; //counting number of typed alleles (total)
	double m_AT; //analytical threshold
	double m_fst; //theta correction 
	double m_pC; //drop-in prob
	double m_lambda; //drop-in lambda (necessary for model validation)
	
	int m_QalleleIndex; //index with Q-allele (-1 means none is Qallele)
	vector<int> m_BWto; //Indicating which alleles backward/forward stutter to: Needed for getContributionIndices function
	vector<int> m_FWto; //Indicating which alleles backward/forward stutter to: Needed for getContributionIndices function

	//Konwn contributor variables: Threated as the first contributors (Mx1,..,MxNOK etc)
	int m_NOK; //number of known contributors
	int m_NOU; //number of unknown contributors
	int m_NumReps; //number of replicates	
	
	Col<int> m_knownContributorGenotypes; //The genotype index for the known contributors (m_NOK length)
	Col<int> m_unknownContributorPower; //indicates which contributorPower value to begin iterate with (used in getContributionIndices)	
	int m_possibleContributionsPerContributor; //this is the (3/6/10)^NOC outcome of compact matrix
	
	//STORAGE OF LARGE OBJECTS:	
	fvec m_pGvec; //Genotype probabilities (nJointCombFULL); //float gives ~1e-6 in precision, double gives 1e-16 in precision

	//Objects for restricting genotype information
	uvec m_jointCombUse; //index storage of which genotypes to use (m_nJointCombREST)	
	fvec m_genoWeightsHighest; //highest calculated genotype probabilities seen so far (nJointCombFULL); dynamical updated (fetch after DC)
		
	//Storage of relationship params
	int m_relGind;
	vector<double> m_ibd;
	
	//Include special contribution for conditional references: 
	umat m_triAlleles; //A matrix specifying a 3rd allele. col1=Alleleposition and col2=contributor number. Implemented from v4.1.2
	
	//The constructor prepares the variables with necessary data (arguments in constructor)
	public:
	size_t m_nJointCombREST; //restricted genotype outcome (<=nGeno1p^m_NOU)
	size_t m_nJointCombFULL; //full genotype outcome (nGeno1p^m_NOU)
	uvec m_contributionIndicesConditionedKnowns; //(numberOfModeledAlleles+nPotentialStutter); //this is the lookup indices to return

	EFMmarker(int markerIndex, int NumReps, int NumAlleles, int startIndMarker_nAlleles, int startIndMarker_nAllelesReps, vector<double> peaksLong, vector<double> freqsLong, vector<double> dropLong, 
		double nTypedLong, vector<double> maTypedLong, vector<int> BWtoLong, vector<int> FWtoLong, 
		int NumPotStutters, int startIndMarkerTot, int QalleleIndex, double pC, double fst, double AT, double lambda,
		int NOK, vector<int> GindKnownLong, int relGindLong, vector<double> ibdLong, umat triAlleles) {

		//Gprob = *Gprob1U; //set genotype probability of 1st unknown
		m_NumReps = NumReps; //set number of replicates
		m_NOK = NOK; //set number of knowns
		m_NumPotStutters = NumPotStutters; //set number of potential stutters
		m_NumAlleles = NumAlleles; //number of alleles (not with potential stutters)
		m_NumAllelesTot = NumAlleles + NumPotStutters; //total number of alleles (includes both allele set and non-observed stutters
		m_NumTypedTot = nTypedLong; //set number of total prev. typed
		m_pC = pC;
		m_lambda = lambda;
		m_fst = fst;
		m_AT = AT;
		m_QalleleIndex = QalleleIndex;
		m_ibd = ibdLong;
		m_relGind = relGindLong;
		m_triAlleles = triAlleles; 
		
		//Following are NumAlleles*NumReps long
		m_peaks.zeros(NumReps,NumAlleles); //init vector
		m_dropinweight.zeros(NumReps,NumAlleles); //init vector
		
		//init vectors:
		for(int i=0; i<NumAlleles; i++) { //for each alleles (observed + dropout)
			m_freqs.push_back( freqsLong[startIndMarker_nAlleles + i] ); //copy value
			m_NumTypedAllele.push_back( maTypedLong[startIndMarker_nAlleles + i] ); //copy value
			m_BWto.push_back( BWtoLong[startIndMarker_nAlleles + i] ); //copy value
			m_FWto.push_back( FWtoLong[startIndMarker_nAlleles + i] ); //copy value
			
			for(int r=0; r<NumReps; r++) { //for each alleles (observed + dropout)
				int idx = r + NumReps*i; //obtain index position adjusted for replicate
				m_peaks(r,i) = peaksLong[startIndMarker_nAllelesReps + idx]; //copy value
				m_dropinweight(r,i) = dropLong[startIndMarker_nAllelesReps + idx]; //copy value			
			}
		}
		//Rcpp::Rcout << m_dropinweight << "\n";

		//Create contribution matrix (1 genotype):
		m_NumGenos1p = int(NumAlleles*(NumAlleles + 1) / 2); //get number of Genotype outcome		
		m_outG1allele.set_size(m_NumGenos1p,2); //init nG1x2 matrix (allele names as indices 0,...,nA-1)
		int cc = 0; //counter oveer all genotypes
		for (int i = 0; i < NumAlleles; i++) {
			for (int j = i; j < NumAlleles; j++) { 
				m_outG1allele(cc,0) = i; //include index
				m_outG1allele(cc,1) = j; //include index
				cc++; //iterate to next genotype outcome
			}
		}		
		//Rcpp::Rcout << outG1contr << "-----------\n" << m_outG1allele;	
		//Note: Can condition on references with missing markers: Make sure that these are sorted: The references with most markers comes first (assume subsets)!
		m_knownContributorGenotypes.zeros(m_NOK);
		for (int k = 0; k < m_NOK; k++) {
			m_knownContributorGenotypes[k] = GindKnownLong[markerIndex*m_NOK + k]; //copy over value (-1 if missing)
		}		
	} //end constructor

	//Function which prepares variables for the repeated likelihood calculations
	//Goes through all combinations ones to store static (unchanged) variables, such as genotype combinations and drop-in model
	void prepareMarker(int NOC, int NumStutterModels, int NumStutterModelsMAX) {		
		m_NumStutterModels = NumStutterModels;//copy: this variable is marker specific
		
		//Prepare index vector for known contributors first:		
		m_possibleContributionsPerContributor = 3; //assuming no stutters
		for(int t = 0; t < NumStutterModelsMAX; t++) {
			m_possibleContributionsPerContributor += 3 + t; //formula for number of outcome per contributor
		}		
		
		//Handling known references here (marker may be missing)
		m_contributionIndicesConditionedKnowns.zeros(m_NumAllelesTot); //declare		
		m_NOU = NOC-m_NOK; //obtain number of unknowns (can also be higher if missing markers)
		m_unknownContributorPower.zeros(NOC); //init Potentially all are unknown
		Col<int> contributorPower = {1}; //Updated for each new conditional (always)		
		int knownsCounter = 0; //count number of non-missing refs 
		int unknownsCounter = 0; //count number of unknowns
		for (int k = 0; k < NOC; k++) {  //traverse all contributors 		
			//Informing the known contributors (update allele lookup indices for these)		
			if( k < m_NOK && m_knownContributorGenotypes[k] > -1) { //obtain value ("-1" means missing marker)
				m_contributionIndicesConditionedKnowns = getContributionIndices(m_knownContributorGenotypes[k], 1,m_outG1allele,m_QalleleIndex,m_NumStutterModels, &(m_BWto[0]),&(m_FWto[0]), 
															m_NumGenos1p, m_contributionIndicesConditionedKnowns, contributorPower);
							
				//ADDING (OPTIONAL) TRI-ALLELE CONTRIBUTION HERE:
				for(uword row=0; row<m_triAlleles.n_rows; row++) {
					int contrIdx = m_triAlleles(row,1); //obtain contribution index
					if(contrIdx==k) { //when same contribution:
						m_contributionIndicesConditionedKnowns = getContributionIndicesAllele(m_triAlleles(row,0),m_QalleleIndex,m_NumStutterModels, &(m_BWto[0]),&(m_FWto[0]), m_contributionIndicesConditionedKnowns, contributorPower);
						//Rcpp::Rcout << m_contributionIndicesConditionedKnowns << "\n";
					}							
				}
				knownsCounter++; //count if non-missing
			} else {				
				m_unknownContributorPower[unknownsCounter] = pow(m_possibleContributionsPerContributor,k); //aligning position of unknown here
				unknownsCounter++; //update number of unknowns
			}			
			contributorPower[0] *= m_possibleContributionsPerContributor; //always update
		}
		m_NOU = unknownsCounter; //update the number of unknown contributors
		m_NOK = knownsCounter; //update the number of known contributors (non-missing marker)				
		m_nJointCombFULL = (size_t) pow(m_NumGenos1p,m_NOU); //obtain full genotype outcome size		
		m_pGvec.set_size(m_nJointCombFULL); //init genotyping vector (+ dropin weights) that is scaled with evidence matrix
		
		//New in version 4.1.x: Allocate extra vector for keeping values
		m_nJointCombREST = m_nJointCombFULL; //init as same (can avoid restriction steps)			
		m_genoWeightsHighest.set_size(m_nJointCombFULL); //init genotyping vector (+ dropin weights) that is scaled with evidence matrix
		m_genoWeightsHighest.fill(-1000000); //insert some very small number (genotype weights must be higher than this)
		//m_genoWeightsHighest.ones( m_nJointCombFULL ); //init dynamic genotyping vector (used in peeling off genotypes). Zero is lowest (non-logged)

		//Store index of where Y is positive (may be irregular)
		vector<uvec> hasPosY_list; //init vector with uvec objects. 
		for(int r=0; r < m_NumReps; r++) { 
			uvec hasPosY = find( m_peaks.row(r) >= m_AT); //obtain indices where Y is positive (per replicate)
			hasPosY_list.push_back(hasPosY); //add to vector (list)
		}
		
		#pragma omp parallel for  //variables are all shared (access)
		for(size_t gind=0; gind < m_nJointCombFULL; gind++) { //traverse all genoypes
			rowvec maTypedvec = m_NumTypedAllele; //Creating copy of counter (per allele)
			double nTyped = m_NumTypedTot; //creating copy of total counter
			
			int genotypePower = 1; //don't think about known contributors here
			double genoProd = 1; //genotype probability
			int UgindPrev=0;
			int aindR=0;
			int bindR=0;			
			for(int k=0; k < m_NOU; k++) { //loop over all unknown contributors				
				int Ugind = (gind/genotypePower) % m_NumGenos1p; // //get genotype for partifular contributor	
				genotypePower *= m_NumGenos1p; //update genotype comb ready for next contributors
				int aindU = m_outG1allele(Ugind,0);
				int bindU = m_outG1allele(Ugind,1);
				
				//Handle related individual (placed last)
				int Rgind = -1; //genotype of related unknown (default)
	
				//Rcpp::Rcout << m_pGvec << "\n";	
				if(m_ibd[0]<0.999 && k==(m_NOU-1)) { //If relationship given and IF last unknown
					if(m_relGind > -1) { //check if related individual was given (do as normal)
						Rgind = m_relGind; //copy					
						aindR = m_outG1allele(Rgind,0); 
						bindR = m_outG1allele(Rgind,1); 
					} else if(m_NOU>=2) { //otherwise we evaluate another type of relationship (requires at least 2 unknowns)
						Rgind = UgindPrev; //set related genotype same as previous unknown genotype
						aindR = m_outG1allele(Rgind,0); 
						bindR = m_outG1allele(Rgind,1); 
					}
				}
				genoProd *= prob_relUnknown(aindU, bindU, Ugind, &(m_freqs[0]), m_fst, &(maTypedvec[0]), nTyped, aindR, bindR, Rgind, &(m_ibd[0]));				
				maTypedvec[aindU] += 1; //update allele count (allele 1)
				maTypedvec[bindU] += 1; //update allele count (allele 2)
				nTyped += 2; //update total count
				
				UgindPrev = Ugind; //store previous unknown allele1
			}
			m_pGvec[gind] = log(genoProd); //store genotype-prob outcome
			
			//get contribution indices for each allele: stored as vector: Known contributors are taken into account through "m_contributionIndicesConditionedKnowns"
			uvec indiceAlleles = getContributionIndices(gind, m_NOU, m_outG1allele,m_QalleleIndex,m_NumStutterModels,  &(m_BWto[0]),&(m_FWto[0]), 
									m_NumGenos1p, m_contributionIndicesConditionedKnowns, m_unknownContributorPower);

			//Obtain info about drop-in alleles: indiceAlleles with value zero means no contribution			
			uvec noContribution = find(indiceAlleles==0); //positions with no contribution
			for(int r=0; r < m_NumReps; r++) { 
				rowvec dropinweight = m_dropinweight.row(r); //obtain row
				uvec dropinInd = intersect(hasPosY_list[r],noContribution); //obtain drop-in allele positions (assume nRep=1)
				int nDropIn = dropinInd.n_elem; //get number of dropin elems
				if( nDropIn>0 ) {
					m_pGvec[gind] +=  sum(dropinweight(dropinInd));
				} else {
					m_pGvec[gind] +=  log(1-m_pC);
				}
			} //end for each replicates
		} //end for all genotypes
		//Rcpp::Rcout << m_pGvec << "\n";	
	} //end indexing function
	
	//Perform restriction based on calculations
	void restrictGenosMarker(double restGenoThresh) {
		//restGenoThresh: Must be in [0,1]. No restriction=0, Full restriction = 1.
		float log_pGvecMax = m_genoWeightsHighest.max(); //obtain max element (already logged)
		m_jointCombUse = find( m_genoWeightsHighest - log_pGvecMax >= log(restGenoThresh)); //return indices which satisfy criterion
		m_nJointCombREST = m_jointCombUse.size(); //obtain size of element		
	}


	//Calculate the genotype weights and check if higher than previous: Always traversing all genotypes
	double calcGenoWeightsMarkerStoreMax(mat EvidWeightsSub) {
		double likEvidMarker = 0.0;  //total sum over all genotypes	(likelihood of evidence for a marker)
		#pragma omp parallel for reduction(+:likEvidMarker) //default(none)
		for(size_t gind=0; gind < m_nJointCombFULL; gind++) { //traverse all genoypes	
			double log_genoWeight = m_pGvec[gind]; //includes prG and drop-in prob							
			//obtaining outcome indices directly for each allele		
			uvec indiceAlleles = getContributionIndices(gind, m_NOU, m_outG1allele,m_QalleleIndex,m_NumStutterModels, &(m_BWto[0]),&(m_FWto[0]),
							m_NumGenos1p, m_contributionIndicesConditionedKnowns, m_unknownContributorPower);			
			for(int aa=0; aa < m_NumAllelesTot; aa++) { //traverse each allele (allele outcome)
				if(indiceAlleles(aa)>0) log_genoWeight += EvidWeightsSub(indiceAlleles[aa],aa); //look up evidence weight based on index
			} //end for each allele	
			if(log_genoWeight > m_genoWeightsHighest[gind]) m_genoWeightsHighest[gind] = log_genoWeight; //insert if higher		
			likEvidMarker += exp(log_genoWeight); //add to global sum			
		}
		return log(likEvidMarker); //return logged value
	}
	
	//Engine to calculate the log-likelihood function: Traverse all genotypes and use calculated lookup matrix:
	double calcLogLikMarker(mat EvidWeightsSub) {
		double likEvidMarker = 0.0;  //total sum over all genotypes	(likelihood of evidence for a marker)		
		//ALWAYS APPLY ASSUMING RESTRICTION (MUST BE EXECUTED AFTER RUNNING restrictGenosMarker otherwise it gives crash...)
		#pragma omp parallel for reduction(+:likEvidMarker) //default(none)
		for(size_t rgind=0; rgind < m_nJointCombREST; rgind++) { //traverse all genoypes			
			uword gind = m_jointCombUse[rgind]; //obtain index of genotype to use (decided beforehand)
			double log_genoWeight = m_pGvec[gind]; //includes prG and drop-in prob							
			uvec indiceAlleles = getContributionIndices(gind, m_NOU, m_outG1allele,m_QalleleIndex,m_NumStutterModels, &(m_BWto[0]),&(m_FWto[0]),
							m_NumGenos1p, m_contributionIndicesConditionedKnowns, m_unknownContributorPower);	
			//look up evidence values for each allele (replicates already accounted for per allele)
			for(int aa=0; aa < m_NumAllelesTot; aa++) { //traverse each allele (allele outcome)
				if(indiceAlleles[aa]>0) log_genoWeight += EvidWeightsSub(indiceAlleles[aa],aa); //look up evidence weight based on index				
			} //end for each allele	
			likEvidMarker += exp(log_genoWeight); //add to global sum
		}
		//Rprintf("%f\n",likEvidMarker);
		return log(likEvidMarker); //return logged value
	}
	
	//THIS IS A NEW FUNCTION FOR MAKING DECONV FAST (not supporting top joint which is never used)
	mat calcMargDCmarker(mat EvidWeightsSub) {		
		#pragma omp declare reduction(vec_double_plus : vector<double> : transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
		vector<double> sumVector(m_NumGenos1p*m_NOU,0);
		#pragma omp parallel for reduction(vec_double_plus : sumVector)		
		for(size_t rgind=0; rgind < m_nJointCombREST; rgind++) { //traverse all genoypes			
			uword gind = m_jointCombUse[rgind]; //obtain index of genotype to use (decided beforehand)
			double log_genoWeight = m_pGvec[gind]; //includes prG and drop-in prob			
			uvec indiceAlleles = getContributionIndices(gind, m_NOU, m_outG1allele,m_QalleleIndex,m_NumStutterModels, &(m_BWto[0]),&(m_FWto[0]),
										m_NumGenos1p, m_contributionIndicesConditionedKnowns, m_unknownContributorPower);						
			for(int aa=0; aa < m_NumAllelesTot; aa++) { //traverse each allele (allele outcome)
				if(indiceAlleles(aa)>0) log_genoWeight += EvidWeightsSub(indiceAlleles[aa],aa); //look up evidence weight based on index
			}							
			double LikVal = exp(log_genoWeight); //obtain likelihood value
			//NOTE: Here it is also possible to remove for genotype probabilities (NB: adjust for drop-in model first)
			
			//lAST PART IS TO INSERT LIKELIHOOD VALUE TO THE MARGINAL MATRIX FOR THE UNKNOWNS
			int genotypePower = 1;
			for(int contributorIndex=0; contributorIndex < m_NOU; contributorIndex++) { //loop over all contributors		
				int genoIdx = (gind/genotypePower) % m_NumGenos1p; // //get genotype for partifular contributor	
				//marginalGenoWeightsUnknowns(contributorIndex,genoIdx) += LikVal; //adding likelihood value to correct genotype category
				sumVector[genoIdx*m_NOU + contributorIndex] += LikVal; //adding likelihood value to correct genotype category
				genotypePower *= m_NumGenos1p;
			}
		}
		
		//Copy values over to the arma structure afterwards:
		mat marginalGenoWeightsUnknowns(m_NOU,m_NumGenos1p,fill::zeros); 
		for(int genoIdx=0; genoIdx < m_NumGenos1p; genoIdx++) {
			for(int contributorIndex=0; contributorIndex < m_NOU; contributorIndex++) { //loop over all contributors		
				marginalGenoWeightsUnknowns(contributorIndex,genoIdx) = sumVector[genoIdx*m_NOU + contributorIndex];
			}
		}
		return(marginalGenoWeightsUnknowns);
	}				

	//THIS IS A NEW FUNCTION FOR MAKING VALIDATION FAST 
	mat calcValidMLEmarker(mat EvidWeightsSub, mat contrMatPerAlleleSub, umat lookUpMat, vector<double> modelValidParams, int *BWfrom, int *FWfrom) {				
		//Unwrap the parameters
		double shape0 = modelValidParams[0];
		double scale0 = modelValidParams[1];
		double xiB = modelValidParams[2];
		double xiF = modelValidParams[3];
		double peakMaxModel = R::qgamma(0.99999, 2*shape0,scale0, 1,0); //max observation in theory (used only for model validation)
		//Rcpp::Rcout << peakMaxModel << "\n";
	
		//Calculate the cumulative expressions for each allele (And each replicate, these cannot be combined)		
		#pragma omp declare reduction(vec_double_plus : vector<double> : transform(omp_out.begin(), omp_out.end(), omp_in.begin(), omp_out.begin(), plus<double>())) initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))
		vector<double> sumVector(2*m_NumAlleles*m_NumReps,0);
		#pragma omp parallel for reduction(vec_double_plus : sumVector)			
		for(size_t rgind=0; rgind < m_nJointCombREST; rgind++) { //traverse all (restricted) genoypes			
			uword gind = m_jointCombUse[rgind]; //obtain index of genotype to use (decided beforehand)
			double log_genoWeight = m_pGvec[gind]; //includes prG and drop-in prob			
			uvec indiceAlleles = getContributionIndices(gind, m_NOU, m_outG1allele,m_QalleleIndex,m_NumStutterModels, &(m_BWto[0]),&(m_FWto[0]),
								m_NumGenos1p, m_contributionIndicesConditionedKnowns, m_unknownContributorPower);		

			for(int alleleIdx=0; alleleIdx < m_NumAllelesTot; alleleIdx++) { //traverse each allele (allele outcome)
				//If having contribution and not beeing allele 'outerAllele' (calculated
				if(indiceAlleles(alleleIdx)>0) log_genoWeight += EvidWeightsSub(indiceAlleles[alleleIdx],alleleIdx); //look up evidence weight based on index
			}							

			//traverse each observed alleles (replicates are handled under each of the two situations)
			double logLikVal_OBS,logLikVal_MAX; //obtain the two necessary statistics				
			for(int outerAllele=0; outerAllele < m_NumAlleles; outerAllele++) { 
				uword contrIdx = indiceAlleles[outerAllele]; //obtain outcome index for the traversed (outer allele)
				
				//Two situations will happen: Peak is dropin OR a non-drop 
				if(contrIdx==0) {  //if the traversed allele 'outerAllele' is dropin (i.e. no contribution)				
					for(int outerRep=0; outerRep < m_NumReps; outerRep++) { //Traverse each replicate	
				
						//Part 1: Remove already integrated drop-in part
						double log_genoWeightCopy = log_genoWeight; //make a copy
						//log_genoWeightCopy -= m_dropinweight(outerRep,outerAllele); //remove particular drop-in weight
						double peak = m_peaks(outerRep,outerAllele); //Obtain peak height			
						if(peak < m_AT) continue; //in case of dropout
					
						//Part 2: Calculate the cumulative drop-in parts
						if( fabs(peak-m_AT) < 1 ) peak += 1; //ensure that peak is greater than AT 
						double ldexpObs = log(m_lambda) - m_lambda*(peak - m_AT); //obtain log of dexp(y-T)
						double pexp1 = 1-exp(-m_lambda*( peak - m_AT)); //cdf for observed peaks
						double pexp2 = 1; //cdf for max peaks (always one)
						logLikVal_OBS = log_genoWeightCopy + log(pexp1) - ldexpObs;
						logLikVal_MAX = log_genoWeightCopy + log(pexp2) - ldexpObs;
						int insIdx = outerAllele*m_NumReps+outerRep; 
						sumVector[insIdx] += exp(logLikVal_OBS); //Insert obs part
						sumVector[insIdx + m_NumAlleles*m_NumReps] += exp(logLikVal_MAX); //Insert max part						
					}
					
				} else { //Otherwise it is an observed allele that must be evaluated with cumulative prob
					//Note: Adjusting the drop-in weight part is not necessary
					
					//Part 1: Obtain shape value for the particular allele:
					//Step 1: non-stutter part
					double shape1 = contrMatPerAlleleSub(lookUpMat(contrIdx,0),outerAllele); //obtain shape argument from contrMat
					shape1 *= (1-xiB-xiF); //loose stutter products (always since this is an observed allele)
					
					//Step 2: BW stutter part:
					if( BWfrom[outerAllele] > -1 ) {
						shape1 += xiB*contrMatPerAlleleSub(lookUpMat(contrIdx,1),BWfrom[outerAllele]);
					}
					//Step 3: FW stutter part:
					if( FWfrom[outerAllele] > -1 ) {
						shape1 += xiF*contrMatPerAlleleSub(lookUpMat(contrIdx,2),FWfrom[outerAllele]);
					}


					//Part 2: Calculate pdf weights for each replicates	(necessary)
					vec weightReps(m_NumReps, fill::zeros);
					for(int innerRep=0; innerRep < m_NumReps; innerRep++) { //Traverse each replicate	
						double peak = m_peaks(innerRep,outerAllele); //Obtain peak height for other replicate
						if(peak < m_AT) { //in case of dropout
							weightReps[innerRep] = R::pgamma(m_AT, shape1,scale0, 1, 1); 
						} else {
							weightReps[innerRep] = R::dgamma(peak, shape1,scale0, 1); 									
						}
					}
						
					//Part 3: Calculate the cumulative pdf weights
					double lowerIntegral = R::pgamma(m_AT-1, shape1,scale0, 1, 0); //common for both terms					
					double maxIntegral = R::pgamma(peakMaxModel, shape1,scale0, 1, 0); //common for any replicates
					
					for(int outerRep=0; outerRep < m_NumReps; outerRep++) { //Traverse each replicate						
						double log_genoWeightCopy = log_genoWeight; //make a copy					
						double peak = m_peaks(outerRep,outerAllele); //Obtain peak height	
						if(peak < m_AT) continue; //dont consider in case of dropout					
						double obsIntegral = R::pgamma(peak, shape1,scale0, 1, 0);				
						
						log_genoWeightCopy -= weightReps[outerRep]; //subtract weight of particular replicate
						logLikVal_OBS = log_genoWeightCopy + log(obsIntegral - lowerIntegral);
						logLikVal_MAX = log_genoWeightCopy + log(maxIntegral - lowerIntegral);						
						int insIdx = outerAllele*m_NumReps+outerRep; 
						sumVector[insIdx] += exp(logLikVal_OBS); //Insert obs part
						sumVector[insIdx + m_NumAlleles*m_NumReps] += exp(logLikVal_MAX); //Insert max part
					}					
				} //end situation (drop-in or non-drop)
			} //end each allele
		} //end each genotype outcome
		
		//Copy values over to the arma structure afterwards:
		mat UaCum(2,m_NumAlleles*m_NumReps,fill::zeros); //int storage	
		for(int outerAllele=0; outerAllele < m_NumAlleles; outerAllele++) { 
			for(int outerRep=0; outerRep < m_NumReps; outerRep++) { //Traverse each replicate	
				int insIdx = outerAllele*m_NumReps+outerRep;
				UaCum(0,insIdx) = sumVector[insIdx];					
				UaCum(1,insIdx) = sumVector[insIdx + m_NumAlleles*m_NumReps];
			}
		}		
		return(UaCum);
	}			
};

//Class to calculate logLik over all markers (also used for optimization options)
class EFMclass {
	private: 
	int m_NumMarkers; //number of markers to traverse
	int m_NOC; //number of contributors
	vector<EFMmarker> EFMmarkerv; //note the type <EFMmarker> and not <*EFMmarker>
	
	//store full vectors of data also here (used for rapid evidWeigh calculations)
	vector<double>  m_PeaksLong;	
	vector<int> 	m_BWfromLong;
	vector<int> 	m_FWfromLong;				
	vector<int> 	m_BWtoLong; //used to recognize if stutters are given away
	vector<int>		m_startIndMarker_nAlleles; //only ordinary alleles (not Potential stutters)
	vector<int>		m_NumPotStutters;
	vector<int>		m_startIndMarker_nAllelesTot;
	vector<int>		m_startIndMarker_nAllelesReps;
	vector<int>		m_NumAlleles;
	vector<int>		m_NumRepsMarkers;
	vector<double>	m_AT;
	Row<double> 	m_basepairLong;

	//helpvariables:
	int m_TotalNumAllelesPS; //total number of alleles (including potential stutters)
	int m_TotalNumAllelesWithReps; //total number of alleles (including replicates for each allele)
	vector<int>		m_lookupAlleleIdx; //used to obtain allele index for given marker and allele-iteration (includes Q-allele)
	vector<int>		m_lookupMarkerIdx; //used to obtain marker index for given allele-iteration (includes Q-allele)
	vector<int>  	m_lookupCumulativeAlleleIdx; //helpfunction look up correct non-potential alleles when traversing all alleles (including potential stutters)

	//ELEMENTS CONCERNING INDEX LOOKUP
	int m_NumCombOut0; //number of elements in contrVec0_key (outcome '0'). "0" meaning that no stutter model is applied
	int m_NumCombOut1; //number of rows in contrMat1_key (outcome '1'). "1" meaning that a stutter model may have been applied
	umat m_contrMat0; //contribution matrix (per contributor for outcome0). 
	umat m_lookUpMat; //lookup index BETWEEN STUTTER OUTCOME AND non-stutter outcome
	int m_NumStutterModelsMAX; //model type (0=Only DEG, 1=BW, 2=BWFW)	
	vector<int> m_NumStutterModels; //vector to take into account possible none-stutter loci (AMEL)

	//Dynamically updated loglikelihood matrices of each outcome (used for calculating logLik and modelValid)
	mat m_EvidWeights; //Evidence weight is stored globally here (also includes potential stutters)
	mat m_contrMatPerAllele;  //contribution matrix per allele (stored here for use in modelValidation as well)
	vector<double> m_modelValidParams; //Copy of different kinds of necessary parameters used in modelValidation)
		
	public: 
	vec m_loglikMarkers; //stored for later access (directly to R)	
	//vec m_loglikMarkersGenos; //stored for later access (directly to R).

	//Main function to make restriction on genotype outcome (based on a threshold)
	double restrictGenos(double restGenoThresh) { 		
		#pragma omp parallel for //Parallel restriction per marker
		for(int m=0; m<m_NumMarkers; m++) {
			(EFMmarkerv[m]).restrictGenosMarker(restGenoThresh); //perform restriction per marker
		}

		//Last: counting up number genotypes to travers (RESTRICTED VS FULL)
		double nJointCombREST = 0; //NB: Note the convertion to double to overcome potential overflow
		double nJointCombFULL = 0;
		for(int m=0; m<m_NumMarkers; m++) {
			nJointCombREST += (EFMmarkerv[m]).m_nJointCombREST;  //add up
			nJointCombFULL += (EFMmarkerv[m]).m_nJointCombFULL;  //add up
		}
		return(nJointCombREST/nJointCombFULL); //return average proportion
	}

	//This is the main function for calculating the EvidenceWeight matrix (done before joint genotype iterations)
	void calcEvidWeights(vector<double> theta) {
		//Fetch and prepare params
		vec Mx(m_NOC);
		for(int k=0; k<m_NOC; k++) Mx[k] = theta[k]; //init
		double mu = theta[m_NOC]; //obtain PHexp param
		double omega = theta[m_NOC+1]; //obtain PHvar param
		double loggedBeta = log(theta[m_NOC+2]); //obtain slope param (take log to easy further calculations)				
		double shape0 = 1/(omega*omega); //get shape param
		double scale0 = mu/shape0;// get scale param
		double const1 = 1/(scale0); //constant 1
	 	double const2 = log(scale0); //constant 2
		double xiB = theta[m_NOC+3]; //obtain stutter prop (BW) param
		double xiF = theta[m_NOC+4]; //obtain stutter prop (FW) param
		m_modelValidParams = {shape0,scale0,xiB, xiF}; //make additional outer copy of these params used for model validation
		
		//PREPARE EVIDENCE MATRIX
		//Step 1: Prepare contribution matrix without stutters (Mix-prop and degradation)		
		vec contrWithMx = m_contrMat0*Mx; //get Mx contribution (no degrad)		
		
		
		Row<double> degScale = shape0 * exp(m_basepairLong * loggedBeta); //vectorized ARMA scale with shape0
		m_contrMatPerAllele = contrWithMx * degScale; //degScale; //an outer matrix product (gives matrix). Contribution per allele (column-wise)		
		//Note: contrMatPerAllele is a (out0 x m_NumAllelesTot) matrix (inclues Q-allele but not potential stutters)
		//Rcpp::Rcout << m_basepairLong << "\n";
		//Rcpp::Rcout << contrWithMx << "\n";
				
		//Step 2: Modify contribution matrix wrt stutters (produce contrMat1 using m_lookUpMat): In same step: Calculate weight-of-evidence matrix
		//Parallelize over all alleles in total (utilizing many threads)		
		#pragma omp parallel for //default(none) shared(EvidWeights,m_TotalNumAllelesPS, m_contrMat0,m_lookUpMat, xiB, xiF,m_NumStutterModels, BWfrom,FWfrom,m_NumCombOut1,m_AT,scale0,const1,const2 ) 
		for(int aindLong=0; aindLong < m_TotalNumAllelesPS; aindLong++) { //Traversing all alllees (also potential stutters)
			//Rcpp::Rcout << "allele=" << aindLong << "\n";
			vec colVEC; //help variable	
			vec contrVEC(m_NumCombOut1,fill::zeros);//help variable, init as zero (contribution vector for all outcome)
			
			//look up alleles (not including potentialstutters) and marker index
			int cumAlleleIdx = m_lookupCumulativeAlleleIdx[aindLong]; //obtain allele to consider: -1 if part of "potential stutte	
			int alleleIdx = m_lookupAlleleIdx[aindLong]; //obtain allele to consider: -1 if part of "potential stutte	
			int markerIdx = m_lookupMarkerIdx[aindLong]; //obtain marker index
			int NumReps = m_NumRepsMarkers[markerIdx]; //number of replicates for given marker			
			int startIndMarker_nAllelesReps = m_startIndMarker_nAllelesReps[markerIdx]; //obtain start index of marker
			double AT = m_AT[markerIdx]; //obtain assumed analytical threshold
			if(cumAlleleIdx > -1) { //Traverse the set of alleles within possible genotypes (also Q-allele)
				colVEC = m_contrMatPerAllele.col( cumAlleleIdx ); 
				contrVEC += colVEC( m_lookUpMat.col(0) ); //contributor alleles
				if( m_BWtoLong[cumAlleleIdx] > -1 ) { //If giving away stutter (better check than using Q-allele index)
					contrVEC *= (1-xiB-xiF); //then loose stutter products (scaling)
				}			
			}
			//Continue and check where to add stutter products
			if( m_BWfromLong[aindLong] > -1 ) {
				colVEC = m_contrMatPerAllele.col( m_startIndMarker_nAlleles[markerIdx] + m_BWfromLong[aindLong]); //obtain allele which contrAllele get BW stutter from. Remember SHIFT
				contrVEC += xiB*colVEC( m_lookUpMat.col(1) ); //modify based on stutter products
			}
			if( m_FWfromLong[aindLong] > -1) {
				colVEC = m_contrMatPerAllele.col( m_startIndMarker_nAlleles[markerIdx] + m_FWfromLong[aindLong]); //obtain allele which contrAllele get BW stutter from. Remember SHIFT
				contrVEC += xiF*colVEC( m_lookUpMat.col(2) ); //modify based on stutter products							
			}
			
			//PREPARE WEIGHT EVIDENCE: 
			double shape1, weight, shapeConst;
			vec colVECcalc(m_NumCombOut1); //help variable (for calculation)
			for(int combOutIdx=0; combOutIdx < m_NumCombOut1; combOutIdx++) { //traverse each unique shape param
				shape1 = contrVEC[combOutIdx]; //get unique shape param
				weight = 0; //reset weight again				
				if( cumAlleleIdx < 0) { //Assuming dropout for all replicates
					weight+= NumReps *(R::pgamma(AT, shape1,scale0, 1, 1)); //calculate for all replicates (identical)
				} else { //Peak heights are expected (but perhaps not for all replicates
					shapeConst = -lgamma(shape1) - shape1*const2; //calculate shape constant only one time //Note: lgamma is from std
					for(int repIdx=0; repIdx < NumReps; repIdx++) {	//Traverse all replicates (for given allele)
						double peak = m_PeaksLong[startIndMarker_nAllelesReps + alleleIdx*NumReps + repIdx]; //obtain peak height wrt long-vector (takes into account replicates)
						if(peak < AT) {
							weight += R::pgamma(AT, shape1,scale0, 1, 1); //calculate for certain replicates
						} else {
							weight += shapeConst + (shape1-1)*log(peak) - peak * const1; 					
						}
					} //end for each replicate					
				} //end if expecting peak heights
				colVECcalc[combOutIdx] = weight;
			}
			m_EvidWeights.col(aindLong) = colVECcalc; //insert calculcated vector			
		} //END CALC		
	}
	
	//This is the main function for calculating the logLik for each marker 
	double calcLogLik(vector<double> theta, bool calcGenoWeights = false) {
		calcEvidWeights(theta);
				
		//Finally calculate per marker
		double jointloglik = 0.0;
		m_loglikMarkers.set_size(m_NumMarkers);
		for(int m=0; m<m_NumMarkers; m++) {
			int fromCol = m_startIndMarker_nAllelesTot[m];
			int toCol = m_startIndMarker_nAllelesTot[m+1] - 1;
			//Rcpp::Rcout << "from=" << fromCol << " to=" << toCol << "\n";			
			if(calcGenoWeights) {
				m_loglikMarkers[m] = (EFMmarkerv[m]).calcGenoWeightsMarkerStoreMax(m_EvidWeights.cols(fromCol,toCol)); //insert marker information
			} else {
				m_loglikMarkers[m] = (EFMmarkerv[m]).calcLogLikMarker(m_EvidWeights.cols(fromCol,toCol) ) ; //insert marker information
			}
			jointloglik += m_loglikMarkers[m];
		}
		return jointloglik; //return joint result
	}
	
	//Obtain list with marginal DC per unknowns
	Rcpp::List calcMargDC() {	
		Rcpp::List retList(m_NumMarkers);
	//	#pragma omp parallel for //Note the parallelization here since this cannot be done within the function
		for(int m=0; m<m_NumMarkers; m++) {
			int fromCol = m_startIndMarker_nAllelesTot[m];
			int toCol = m_startIndMarker_nAllelesTot[m+1] - 1;
			mat margDCmarker = (EFMmarkerv[m]).calcMargDCmarker(m_EvidWeights.cols(fromCol,toCol) ); //calculate for marker
			retList[m] = margDCmarker; //insert matrix
		}	
		return retList;
	}	
	
	//Obtain list with validMLE per marker
	Rcpp::List calcValidMLE() {	
		Rcpp::List retList(m_NumMarkers);
		//#pragma omp parallel for //Note the parallelization here since this cannot be done within the function
		for(int m=0; m<m_NumMarkers; m++) {
			int fromCol = m_startIndMarker_nAllelesTot[m]; //note that this also includes potential stutters (tot)
			int toCol = m_startIndMarker_nAllelesTot[m+1] - 1;
			mat contrMatPerAlleleSub = m_contrMatPerAllele.cols(m_startIndMarker_nAlleles[m],m_startIndMarker_nAlleles[m+1] - 1); //make a copy of evidWeight
			mat validMLEmarker = (EFMmarkerv[m]).calcValidMLEmarker(m_EvidWeights.cols(fromCol,toCol), contrMatPerAlleleSub, m_lookUpMat, m_modelValidParams, &(m_BWfromLong[fromCol]), &(m_FWfromLong[fromCol]) ); 
			retList[m] = validMLEmarker; //insert matrix
		}
		return retList;
	}	
		
	//The main prepare function (must be done before likelihood calculations!!)
	//Here the allele-outcome matrices are prepared. Only depending on the number of contributors.
	//A contribution lookup-matrix is created to link unstuttered (contrMat0) and stuttered relation (i.e., contrMat1 which is not needed since "getContributionIndices" function already gives index)
	void prepare(int NOC) {
		m_NOC = NOC;
		uvec out0 = { 0, 1, 2 }; //defined outcomes for key0		
		int nOut0 = 3; //size of outcome (no stutter)
		int nOut1; //size of outcome (with possible stutters)
		int nCols; //number of columns used for additional stutters (1 + BW + FW)
		umat out1; //nCols x nOutcome (transposed)
		if(m_NumStutterModelsMAX==0) {
			out1 = {{0,1,2}}; //defined outcomes for key1 (given per column)	
		} else if(m_NumStutterModelsMAX==1) {
			out1 = { {0,1,2,0,1,0},		   {0,0,0,1,1,2}}; //defined outcomes for key1 (given per row)			
		} else { //(numStutterModels0==2) {
			out1 = { {0,1,2,0,1,0,0,1,0,0},{0,0,0,1,1,2,0,0,0,1},{0,0,0,0,0,0,1,1,2,1} }; //defined outcomes for key1 (given per column)			
		}
		nOut1 = out1.n_cols; //obtain number of outcome
		nCols = out1.n_rows; //obtain number of stutter outcome (notice that out1 is transposed)
		//Note that: nCols==numStutterModels0
		//Rcpp::Rcout << out1 << "\n";
		//Rcpp::Rcout << nCols << "\n";
		//Rcpp::Rcout << nOut1 << "\n";
		
		//Expaning outcome wrt number of contributors
		m_NumCombOut0 = (int) pow(nOut0,NOC); //number of outcomes (no stutter effect)
		m_NumCombOut1 = (int) pow(nOut1,NOC); //number of outcomes (with BW+FW stutter effects)			
		m_contrMat0.zeros(m_NumCombOut0,NOC);  //contribution matrix
		m_lookUpMat.zeros(m_NumCombOut1,nCols); //lookup matrix (used to connect non-stutter outcome and stuttered outcome)		
		uvec contrMat0_key(m_NumCombOut0,fill::zeros);  //Key vector (no stutter)		
		
		//INIT BIG EvidWeight matrix (globally for all markers)
		m_EvidWeights.zeros(m_NumCombOut1 , m_TotalNumAllelesPS); //Init Weight of evidence matrix (important with init different from zero)

		//CONSTRUCTION OF lookUpMat (same for all markers).
		//CONSTRUCT Key0: No STUTTER
		int contrPower; //used to loop combinations
		int combind; //combination index outcome
		int base10;//used to indicate which 10-base letter (used to create key)
		for(int ind=0; ind < m_NumCombOut0; ind++) { //traverse all combinations	
			base10 = 1; 
			contrPower = 1;
			for(int k=0; k<NOC; k++) {
				combind =  (ind/contrPower) % nOut0; // get outcome index
				uword contr = out0[combind]; //get contribution				
				m_contrMat0(ind,k) = contr; //notice change of order (and index -1 to start from 0)
				contrMat0_key[ind] += contr*base10; //notice the base here (add to number)
				base10 *= 10; 
				contrPower *= nOut0; //update for combination traversion
			}	
		} //end for all combinations
		//Rcpp::Rcout << m_contrMat0 << "\n";
		//Rcpp::Rcout << contrMat0_key << "\n";	

		//DEFINED lookUpMat BY CONSTRUCTING Key1(With STUTTER) on the go
		//Must build up a key across all stutter types (nCols)
		for(int ind=0; ind < m_NumCombOut1; ind++) { //traverse all combinations	
			base10 = 1; 
			contrPower = 1;
			uvec keyVec(nCols,fill::zeros); //init key-vector as zero (numStutter models + 1)
			//IDEA IS TO FILL keyVec UP (all 3 values for all outcomes), STEP 1
			//AND MAP THE KEY TO contrMat_key index, STEP 2
			for(int k=0; k<NOC; k++) {
				combind =  (ind/contrPower) % nOut1; // get outcome index			
				uvec contrs = out1.col(combind); //get specific column (contribution)
				keyVec += contrs*base10; //notice the base here (add to number)
				base10 *= 10; //update base afterwards
				contrPower *= nOut1; //update for combination traversion
			}			
			//Rcpp::Rcout << keyVec.t() << "\n";
					
			//LAST: LOOK UP CORRESPONDING INDEX OF contrMat0_key
			for(int j=0; j<nCols; j++) {
				uvec indVec = find( contrMat0_key==keyVec[j]); //store index (starts from 0)
				m_lookUpMat(ind,j) = indVec[0]; //store index 				
			}		
		} //end for all combinations
		//Rcpp::Rcout << m_lookUpMat << "\n";
		
		//EXECUTES PRE-CALCULATIONS FOR EACH MARKER (this is data specific)
		for(int m=0; m<m_NumMarkers; m++) { 
			//send marker specific 'number of stutter-model'
			(EFMmarkerv[m]).prepareMarker(NOC, m_NumStutterModels[m],m_NumStutterModelsMAX); //Note: Also include max 
		}
		//Rcpp::Rcout << m_EvidWeights << "\n";
	}

	//Constructor:
	EFMclass(vector<int>  NumStutterModels, int NumMarkers, vector<int> NumRepsMarkers, vector<int> NumAlleles, vector<int> startIndMarker_nAlleles, vector<int> startIndMarker_nAllelesReps, 
		vector<double> peaksLong, vector<double> freqsLong, vector<double> dropLong, 
		vector<double> nTypedLong, vector<double> maTypedLong, vector<double> basepairLong, 
		vector<int> BWfromLong, vector<int> FWfromLong, vector<int> BWtoLong, vector<int> FWtoLong,
		vector<int> NumPotStutters, vector<int> startIndMarker_nAllelesTot, vector<int> QalleleIndex, vector<double> pC, vector<double> fst, vector<double> AT, vector<double> lambda,
		int NOK, vector<int> GindKnownLong,  vector<int> relGindLong, vector<double> ibdLong, vector<int> triAlleles) {
				
		m_NumStutterModels = NumStutterModels;
		m_NumMarkers = NumMarkers; 
		m_PeaksLong = peaksLong;
		m_BWfromLong = BWfromLong;
		m_FWfromLong = FWfromLong;
	 	m_BWtoLong = BWtoLong;
		m_startIndMarker_nAlleles = startIndMarker_nAlleles;
		m_startIndMarker_nAllelesTot = startIndMarker_nAllelesTot;
		m_startIndMarker_nAllelesReps = startIndMarker_nAllelesReps;
		m_AT = AT;
		m_NumRepsMarkers = NumRepsMarkers;
		
		//PREPARE allele lookup
		m_NumStutterModelsMAX = 0; //obtain upper stutter model (to take into account that some marker may not use stutter model)
		m_TotalNumAllelesPS = 0; //Total number of alleles (both allele set and non-observed potential stutters) 
		m_TotalNumAllelesWithReps = 0; //Total number of alleles (only allele set + replicates)
		int TotalNumAlleles = 0; //Total number of alleles (only allele set). Use only locally here		
		int counter=0; //allele counter (not counting potential stutter alleles)
		for(int m=0; m<m_NumMarkers; m++) {	 //traverse each marker
			for(int a=0; a < NumAlleles[m]; a++) {
				m_lookupCumulativeAlleleIdx.push_back(counter++); //include counter for allele outcome (include Q-allele but not PotentialStutters)
				m_lookupAlleleIdx.push_back(a); //simply include allele index
			}
			for(int a=0; a < NumPotStutters[m]; a++) {
				m_lookupCumulativeAlleleIdx.push_back(-1); //insert -1 for PotentialStutters				
				m_lookupAlleleIdx.push_back(-1); //insert -1 for PotentialStutters				
			}
			for(int a=0; a < (NumAlleles[m]+NumPotStutters[m]); a++)  m_lookupMarkerIdx.push_back( m ); //append marker index						
			TotalNumAlleles += NumAlleles[m]; 
			m_TotalNumAllelesPS += (NumAlleles[m]+NumPotStutters[m]);
			m_TotalNumAllelesWithReps += (NumAlleles[m]*NumRepsMarkers[m]);
			if( NumStutterModels[m] > m_NumStutterModelsMAX) m_NumStutterModelsMAX = NumStutterModels[m]; //set as highest observed
		}
		
		//m_basepairLong(basepairLong);
		m_basepairLong.set_size(TotalNumAlleles); //init
		for(int a=0; a < TotalNumAlleles; a++)  m_basepairLong[a] = basepairLong[a]; //copy to vector (need other vector type)
		
		//structure information about tri-alleles: 
		vector<umat> triAllelesPerMarker; //create vector structure to store triAlleles matrix per marker
		umat triAllelesEmpty(0,2); //init empty object (in case of none)
		umat triAlleles1row(1,2); //init empty object (in case of none)		
		int nTriAlleles = triAlleles.size()/3; //obtain size of vector and adjust for number of elements	
		int elemIdxCounter = 0; //Used to iterate over the elements
		for(int m=0; m< NumMarkers; m++) {
			umat triAlleleTmp = triAllelesEmpty; //set as empty by default if no data
			for(int idx=elemIdxCounter; idx<nTriAlleles; idx++) { //traverse each row in new matrix				
				int marker = triAlleles[3*idx]; //Obtain marker index
				if(marker==m) { //Store object 
					umat newRow = triAlleles1row; 
					newRow(0,0) = triAlleles[3*idx+1];
					newRow(0,1) = triAlleles[3*idx+2];					
					triAlleleTmp = join_cols(triAlleleTmp,newRow); //append row
					elemIdxCounter++; //remember to update this as well
				} else if(marker>m) {
					break; //stop loop if marker exceed
				}
			}
			//Rcpp::Rcout << triAlleleTmp << "\n";
			triAllelesPerMarker.push_back(triAlleleTmp); //insert to vector
		}
			
		//Prepare EFMmarker object
		for(int m=0; m< NumMarkers; m++) {
			EFMmarker *marker = new EFMmarker(m, NumRepsMarkers[m], NumAlleles[m], startIndMarker_nAlleles[m],startIndMarker_nAllelesReps[m], peaksLong, freqsLong, dropLong, nTypedLong[m], maTypedLong, BWtoLong, FWtoLong, NumPotStutters[m], startIndMarker_nAllelesTot[m], QalleleIndex[m], pC[m], fst[m], AT[m], lambda[m], NOK, GindKnownLong,  relGindLong[m], ibdLong, triAllelesPerMarker[m]); //prepare new marker
			EFMmarkerv.push_back(*marker); //insert marker information to vector (push back gives correct marker order)
		}
		
		
	}
};


//REQUIRED FOR EXPOSED IMPLEMENTATION:
//This class is used as interface to EFMclass with functions
class ExposedClass { 
	//private:
		//vector<EFMclass> efmStorage; //must store data in vector object
	EFMclass *efmStorage;
	int m_maxThreads = 0;
	
	public:      
		//exposed function to fill data structure
		void filldata(vector<int> NumStutterModels, int NumMarkers, vector<int> NumRepsMarkers, vector<int> NumAlleles, 
		    vector<int> startIndMarker_nAlleles, vector<int> startIndMarker_nAllelesReps, 
			vector<double> peaksLong, vector<double> freqsLong, vector<double> dropLong, 
			vector<double> nTypedLong, vector<double> maTypedLong, vector<double> basepairLong, 			
			vector<int> BWfromLong, vector<int> FWfromLong, vector<int> BWtoLong, vector<int> FWtoLong,
			vector<int> NumPotStutters, vector<int> startIndMarker_nAllelesTot, vector<int> QalleleIndex, 
			vector<double> pC, vector<double> fst, vector<double> AT, vector<double> lambda, int NOK, vector<int> GindKnownLong, 
			vector<int> relGindLong, vector<double> ibdLong, vector<int> triAlleles, int maxThreads) {
				
			m_maxThreads = maxThreads; //copy
			efmStorage = new EFMclass(NumStutterModels, NumMarkers, NumRepsMarkers, NumAlleles, startIndMarker_nAlleles, startIndMarker_nAllelesReps, peaksLong, freqsLong, dropLong, nTypedLong, maTypedLong, basepairLong, 
			BWfromLong, FWfromLong, BWtoLong, FWtoLong, NumPotStutters, startIndMarker_nAllelesTot, QalleleIndex, pC, fst, AT, lambda, NOK, GindKnownLong, relGindLong, ibdLong, triAlleles);			
			//efmStorage.push_back(*efm); //store EFMclass object			
		}
		
		//exposed function to create indexing: Drop-in weights already calculated (but still need pC)
		void prepare(int NOC) {		
			#ifdef _OPENMP
				int numThreads = thread::hardware_concurrency();
				if(m_maxThreads>0) numThreads = min(numThreads,m_maxThreads); 
				omp_set_num_threads(numThreads); //preparing CPU parallelization	
			#endif
			if(efmStorage) efmStorage->prepare(NOC); 
		}
		
		//exposed function to restrict genotype outcome, done by comparing weights with input threshold
		double restrictGenos(double restGenoThresh) { 
			#ifdef _OPENMP
				int numThreads = thread::hardware_concurrency();
				if(m_maxThreads>0) numThreads = min(numThreads,m_maxThreads); 
				omp_set_num_threads(numThreads); //preparing CPU parallelization	
			#endif

			if(efmStorage) {
				return  efmStorage->restrictGenos(restGenoThresh); //perform restriction			
			} else {
				return NAN;
			}
		}

		//exposed function to calculate genotype weights (keep max)
		double calcGenoWeightsMax(vector<double> theta) { 
			#ifdef _OPENMP
				int numThreads = thread::hardware_concurrency();
				if(m_maxThreads>0) numThreads = min(numThreads,m_maxThreads); 
				omp_set_num_threads(numThreads); //preparing CPU parallelization	
			#endif
			if(efmStorage) {
				return efmStorage->calcLogLik(theta, true); //calculate genotype weights and also return likelihood function			
			} else {
				return NAN;
			}
		}

		//exposed function to caclulate the likelihood func
		double loglik(vector<double> theta) { 
			#ifdef _OPENMP
				int numThreads = thread::hardware_concurrency();
				if(m_maxThreads>0) numThreads = min(numThreads,m_maxThreads); 
				omp_set_num_threads(numThreads); //preparing CPU parallelization	
			#endif
			if(efmStorage) {
				return efmStorage->calcLogLik(theta); //calculate and return likelihood function			
			} else {
				return NAN;
			}
		}

		//exposed function to caclulate the likelihood func: Need AT for cdf calculations
		vec loglikMarkers() { 
			if(efmStorage) {				
				return( efmStorage->m_loglikMarkers );
			} else {
				return NULL;
			}
		}
		
		//exposed function to calculate marginal genotype weights of unknowns (run loglik to set evidWeights first)
		Rcpp::List calcMargDC() { //resturn
			#ifdef _OPENMP
				int numThreads = thread::hardware_concurrency();
				if(m_maxThreads>0) numThreads = min(numThreads,m_maxThreads); 
				omp_set_num_threads(numThreads); //preparing CPU parallelization	
			#endif
			if(efmStorage) {
				return efmStorage->calcMargDC(); //calculate and return	
			} else {
				return NAN;
			}
		}

		//exposed function to calculate model validation statistics (run loglik to set evidWeights first)
		Rcpp::List calcValidMLE() { //resturn
			#ifdef _OPENMP
				int numThreads = thread::hardware_concurrency();
				if(m_maxThreads>0) numThreads = min(numThreads,m_maxThreads); 
				omp_set_num_threads(numThreads); //preparing CPU parallelization	
			#endif
			if(efmStorage) {
				return efmStorage->calcValidMLE(); //calculate and return		
			} else {
				return NAN;
			}
		}
		
		//exposed function to remove storage
		void close() {
			if(efmStorage) {
				delete efmStorage; //free memory		
				efmStorage = NULL;
			}
		}
		
};


RCPP_MODULE(mod){ //mod is name of module (needed for loading exposed class)
    using namespace Rcpp ;

    class_<ExposedClass>( "ExposedClass" )
	.constructor()

    // Get access to functions: Must be called sequentially
	.method("filldata", &ExposedClass::filldata , "Fill object with data") //Step 1
	.method("prepare", &ExposedClass::prepare , "Prepare data for calculations") //Step 2
	.method("calcGenoWeightsMax", &ExposedClass::calcGenoWeightsMax , "Calculate genotype weights (stored in vector, only maximum)") //Step 3
	.method("restrictGenos", &ExposedClass::restrictGenos , "Construct restricted genotype vector") //Step 4
	.method("loglik", &ExposedClass::loglik , "Calculate and obtain likelihood value") //Step 5
	.method("logliki", &ExposedClass::loglikMarkers , "Obtain marker specific likelihood values") //Step 6	
	.method("calcMargDC", &ExposedClass::calcMargDC, "Obtain marginal likelihood for all outcome") //Step 7 (DC)
	.method("calcValidMLE", &ExposedClass::calcValidMLE, "Obtain model validation results") //Step 8 (Valid)
	.method("close", &ExposedClass::close , "Free memory")
	;
}
