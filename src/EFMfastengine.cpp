//WRITTEN BY OYVIND BLEKA (2022). ALGORITHM INSPIRED BY JERRY HOOGENBOOM
//This script contains the fast implementation of EFM (new algorithm)
//Used for following calculations: MLE (optimization), MCMC (simulations), INT (numerical integral)

#include <vector> //vector storage
#include <thread> //used to obtain number of logical processes
#include <RcppArmadillo.h> //necessary for functions find, intercept
#include <Rmath.h> //includes pgamma (note that this has namespace 'R')
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
//Step 3: loglik - used to calculate the likelihood function for a given parameter set. Consists of 2 steps:
//	Step 3a: calculation of the Likelihood matrix EvidWeights (carried out across all data).
//	Step 3b: calculation of the total sum by traversing all genotypes (carried out for each marker).
*/

class EFMmarker { //Each marker is treated separately
	private: 
	int m_NumGenos1p; //number of genotypes for 1 contributor 
	Mat<uword> m_outG1allele; //allele index per genotype outcome	

	//DATA VARIABLES
	int m_NumPotStutters; //number of unobserved potential stutters (BW+FW)
	int m_NumAllelesTot; //total number of alleles (both alleles and potential stutters)

	mat m_peaks; //Intensity vector (copied number)
	mat m_dropinweight; //Dropin weights vector (pC*freq*dexp(Y-AT,lambda)	logged
	vector<double> m_freqs; //Frequency vector
	vector<double> m_NumTypedAllele; //counting number of typed alleles of type a
	double m_NumTypedTot; //counting number of typed alleles (total)
	double m_AT; //analytical threshold
	double m_fst; //theta correction 
	double m_pC; //drop-in prob
	int m_QalleleIndex; //index with Q-allele (-1 means none is Qallele)
	vector<int> m_BWto; //Indicating which alleles backward/forward stutter to: Needed for getContributionIndices function
	vector<int> m_FWto; //Indicating which alleles backward/forward stutter to: Needed for getContributionIndices function

	//Konwn contributor variables: Threated as the first contributors (Mx1,..,MxNOK etc)
	int m_NOK; //number of known contributors
	int m_NOU; //number of unknown contributors
	int m_NumReps; //number of replicates	
	
	Col<int> m_knownContributorGenotypes; //The genotype index for the known contributors (m_NOK length)
	Col<int> m_unknownContributorPower; //indicates which contributorPower value to begin iterate with (used in getContributionIndices)
	uvec m_contributionIndicesConditionedKnowns; //(numberOfModeledAlleles+nPotentialStutter); //this is the lookup indices to return
	int m_possibleContributionsPerContributor; //this is the (3/6/10)^NOC outcome of compact matrix
	
	//STORAGE OF LARGE OBJECTS:
	int m_nJointCombFULL; //full genotype oucome (nGeno1p^m_NOU)
	Col<float> m_pGvec; //Genotype probabilities (nJointCombFULL); //float gives ~1e-6 in precision, double gives 1e-16 in precision
	int m_NumStutterModels; //model type (0=Only DEG, 1=BW, 2=BWFW)

	//Storage of relationship params
	int m_relGind;
	vector<double> m_ibd;
	//The constructor prepares the variables with necessary data (arguments in constructor)
	public:
	EFMmarker(int markerIndex, int NumReps, int NumAlleles, int startIndMarker_nAlleles, int startIndMarker_nAllelesReps, vector<double> peaksLong, vector<double> freqsLong, vector<double> dropLong, 
		double nTypedLong, vector<double> maTypedLong, vector<int> BWtoLong, vector<int> FWtoLong, 
		int NumPotStutters, int startIndMarkerTot, int QalleleIndex, double pC, double fst, double AT,
		int NOK, vector<int> GindKnownLong, int relGindLong, vector<double> ibdLong) {

		//Gprob = *Gprob1U; //set genotype probability of 1st unknown
		m_NumReps = NumReps; //set number of replicates
		m_NOK = NOK; //set number of knowns
		m_NumPotStutters = NumPotStutters; //set number of potential stutters
		m_NumAllelesTot = NumAlleles + NumPotStutters; //total number of alleles (includes both allele set and non-observed stutters
		m_NumTypedTot = nTypedLong; //set number of total prev. typed
		m_pC = pC;
		m_fst = fst;
		m_AT = AT;
		m_QalleleIndex = QalleleIndex;
		m_ibd = ibdLong;
		m_relGind = relGindLong;
		
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

	//Function which prepares variables for the repeated genotype calculations 
	//Goes through all combinations ones to store static (unchanged) variables, such as genotype combinations and drop-in model
	void prepareMarker(int NOC, int NumStutterModels) {		
		m_NumStutterModels = NumStutterModels;//copy: this variable is marker specific
		
		//Prepare index vector for known contributors first:		
		m_possibleContributionsPerContributor = 3; //assuming no stutters
		for(int t = 0; t < m_NumStutterModels; t++) {
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
				knownsCounter++; //count if non-missing
			} else {				
				m_unknownContributorPower[unknownsCounter] = pow(m_possibleContributionsPerContributor,k); //aligning position of unknown here
				unknownsCounter++; //update number of unknowns
			}			
			contributorPower[0] *= m_possibleContributionsPerContributor; //always update
		}
		m_NOU = unknownsCounter; //update the number of unknown contributors
		m_NOK = knownsCounter; //update the number of known contributors (non-missing marker)				
		m_nJointCombFULL = (int) pow(m_NumGenos1p,m_NOU); //obtain full genotype outcome size		
		m_pGvec.set_size(m_nJointCombFULL); //init genotyping vector (+ dropin weights) that is scaled with evidence matrix

		//Store index of where Y is positive (may be irregular)
		vector<uvec> hasPosY_list; //init vector with uvec objects. 
		for(int r=0; r < m_NumReps; r++) { 
			uvec hasPosY = find( m_peaks.row(r) >= m_AT); //obtain indices where Y is positive (per replicate)
			hasPosY_list.push_back(hasPosY); //add to vector (list)
		}
		
		#pragma omp parallel for  //variables are all shared (access)
		for(int gind=0; gind < m_nJointCombFULL; gind++) { //traverse all genoypes
			rowvec maTypedvec = m_NumTypedAllele; //Creating copy of counter (per allele)
			double nTyped = m_NumTypedTot; //creating copy of total counter
			
			int genotypePower = 1; //don't think about known contributors here
			double genoProd = 1; //genotype probability
			for(int k=0; k < m_NOU; k++) { //loop over all unknown contributors				
				int Ugind = (gind/genotypePower) % m_NumGenos1p; // //get genotype for partifular contributor	
				genotypePower *= m_NumGenos1p; //update genotype comb ready for next contributors
				int aindU = m_outG1allele(Ugind,0);
				int bindU = m_outG1allele(Ugind,1);
				
				//Handle related individual (placed last)
				int Rgind = -1; //genotype of related unknown
				int aindR = 0;
				int bindR = 0;
				if(k==(m_NOU-1) && m_relGind > -1) { //check if last visited
					Rgind = m_relGind; //copy
					aindR = m_outG1allele(Rgind,0); 
					bindR = m_outG1allele(Rgind,1); 
				}
				genoProd *= prob_relUnknown(aindU, bindU, Ugind, &(m_freqs[0]), m_fst, &(maTypedvec[0]), nTyped, aindR, bindR, Rgind, &(m_ibd[0]));				
				maTypedvec[aindU] += 1; //update allele count (allele 1)
				maTypedvec[bindU] += 1; //update allele count (allele 2)
				nTyped += 2; //update total count				
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

	
	//Engine to calculate the log-likelihood function
	double calcLogLikMarker(mat EvidWeightsSub) {
		//Traverse all genotypes and use calculated lookup matrix:
		double bigsum = 0.0;  //total sum over all genotypes		
		//return log(bigsum); //return logged value
		#pragma omp parallel for reduction(+:bigsum) //default(none) shared(m_NumAllelesTot,EvidWeightsSub,m_pGvec,m_NOU,m_outG1allele,m_QalleleIndex,m_NumStutterModels,m_NumPotStutters,m_BWto,m_FWto,m_NumGenos1p,m_contributionIndicesConditionedKnowns,m_unknownContributorPower) 
		for(int gind=0; gind < m_nJointCombFULL; gind++) { //traverse all genoypes	
			double logevidProb = m_pGvec[gind]; //includes prG and drop-in prob							
			//obtaining outcome indices directly for each allele		
			uvec indiceAlleles = getContributionIndices(gind, m_NOU, m_outG1allele,m_QalleleIndex,m_NumStutterModels, &(m_BWto[0]),&(m_FWto[0]),
							m_NumGenos1p, m_contributionIndicesConditionedKnowns, m_unknownContributorPower);
			
			//look up evidence values for each allele (replicates already accounted for per allele)
			for(int aa=0; aa < m_NumAllelesTot; aa++) { //traverse each allele (allele outcome)
				if(indiceAlleles[aa]==0) continue; //note that this could be avoided to use? Using skip is a bit faster than not using it.			
				logevidProb += EvidWeightsSub(indiceAlleles[aa],aa); //look up evidence weight based on index
			} //end for each allele	
			//Rprintf("gind=%i, loglik=%f\n",gind,logevidProb);			
			//Rprintf( "Thread %d works with idx %d\n", omp_get_thread_num(), gind);
			bigsum += exp(logevidProb); //add to global sum
		}
		//Rprintf("%f\n",bigsum);
		return log(bigsum); //return logged value
	}
	
	/*Calculating the log-likelihood function per genotype (return large vector). NOT USED
	Col<float> calcLogLikGenos(mat EvidWeightsSub) {
		//Traverse all genotypes and use calculated lookup matrix:		
		Col<float> logLikGeno(m_nJointCombFULL); //init full vector
		for(int gind=0; gind < m_nJointCombFULL; gind++) { //traverse all genoypes	
			double logevidProb = m_pGvec[gind]; //includes prG and drop-in prob											
			logLikGeno[gind] = logevidProb; //add to global sum
		}
		return logLikGeno; //return logged value
	}
	*/
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
	vector<int>		m_startIndMarker_nAlleles;
	vector<int>		m_NumPotStutters;
	vector<int>		m_startIndMarker_nAllelesTot;
	vector<int>		m_startIndMarker_nAllelesReps;
	vector<int>		m_NumAlleles;
	vector<int>		m_NumRepsMarkers;
	vector<double>	m_AT;
	Row <double> m_basepairLong;

	//helpvariables:
	int m_TotalNumAllelesPS; //total number of alleles (including potential stutters)
	vector<int>		m_lookupAlleleIdx; //used to obtain allele index for given marker and allele-iteration (includes Q-allele)
	vector<int>		m_lookupMarkerIdx; //used to obtain marker index for given allele-iteration (includes Q-allele)
	vector<int>  	m_lookupCumulativeAlleleIdx; //helpfunction look up correct non-potential alleles when traversing all alleles (including potential stutters)

	//ELEMENTS CONCERNING INDEX LOOKUP
	int m_NumCombOut0; //number of elements in contrVec0_key (outcome '0'). "0" meaning that no stutter model is applied
	int m_NumCombOut1; //number of rows in contrMat1_key (outcome '1'). "1" meaning that a stutter model may have been applied
	Mat<uword> m_contrMat0; //contribution matrix (per contributor for outcome0). 
	Mat<uword> m_lookUpMat; //lookup index BETWEEN STUTTER OUTCOME AND non-stutter outcome
	int m_NumStutterModelsMAX; //model type (0=Only DEG, 1=BW, 2=BWFW)	
	vector<int> m_NumStutterModels; //vector to take into account possible none-stutter loci (AMEL)
	mat m_EvidWeights; //Evidence weight is stored globally here: Also used as a cache to indicate which outcomes to calc (values>0)
	
	public: 
	Col<double> m_loglikMarkers; //stored for later access (directly to R)	
	//Col<double> m_loglikMarkersGenos; //stored for later access (directly to R).
	
	//calculating (joint) loglik val for given param
	//In this class the preparation of the EvidWeight matrix is already made (looping across all alleles_tot).
	double calcLogLik(vector<double> theta) {
		//Fetch and prepare params
		Col<double> Mx(m_NOC);
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

		//PREPARE EVIDENCE MATRIX
		//Step 1: Prepare contribution matrix without stutters (Mix-prop and degradation)		
		Col<double> contrWithMx = m_contrMat0*Mx; //get Mx contribution (no degrad)		
		Row <double> degScale = shape0 * exp(m_basepairLong * loggedBeta); //vectorized ARMA scale with shape0
		mat contrMatPerAllele = contrWithMx * degScale; //degScale; //an outer matrix product (gives matrix). Contribution per allele (column-wise)
		//Note: contrMatPerAllele is a (out0 x m_NumAllelesTot) matrix (inclues Q-allele but not potential stutters)
		//Rcpp::Rcout << m_basepairLong << "\n";
		//Rcpp::Rcout << contrWithMx << "\n";
		
		//Step 2: Modify contribution matrix wrt stutters (produce contrMat1 using m_lookUpMat): In same step: Calculate weight-of-evidence matrix
		//Parallelize over all alleles in total (utilizing many threads)		
		#pragma omp parallel for //default(none) shared(EvidWeights,m_TotalNumAllelesPS, m_contrMat0,m_lookUpMat, xiB, xiF,m_NumStutterModels, BWfrom,FWfrom,m_NumCombOut1,m_AT,scale0,const1,const2 ) 
		for(int aindLong=0; aindLong < m_TotalNumAllelesPS; aindLong++) { //Traversing all alllees (also potential stutters)
			//Rcpp::Rcout << "allele=" << aindLong << "\n";
			Col <double> colVEC; //help variable	
			Col <double> contrVEC(m_NumCombOut1,fill::zeros);//help variable, init as zero (contribution vector for all outcome)
			
			//look up alleles (not including potentialstutters) and marker index
			int cumAlleleIdx = m_lookupCumulativeAlleleIdx[aindLong]; //obtain allele to consider: -1 if part of "potential stutte	
			int alleleIdx = m_lookupAlleleIdx[aindLong]; //obtain allele to consider: -1 if part of "potential stutte	
			int markerIdx = m_lookupMarkerIdx[aindLong]; //obtain marker index
			int NumReps = m_NumRepsMarkers[markerIdx]; //number of replicates for given marker			
			int startIndMarker_nAllelesReps = m_startIndMarker_nAllelesReps[markerIdx]; //obtain start index of marker
			double AT = m_AT[markerIdx]; //obtain assumed analytical threshold
			if(cumAlleleIdx > -1) { //Traverse the set of alleles within possible genotypes (also Q-allele)
				colVEC = contrMatPerAllele.col( cumAlleleIdx ); 
				contrVEC += colVEC( m_lookUpMat.col(0) ); //contributor alleles
				if( m_BWtoLong[cumAlleleIdx] > -1 ) { //If giving away stutter (better check than using Q-allele index)
					contrVEC *= (1-xiB-xiF); //then loose stutter products (scaling)
				}			
			}
			//Continue and check where to add stutter products
			if( m_BWfromLong[aindLong] > -1 ) {
				colVEC = contrMatPerAllele.col( m_startIndMarker_nAlleles[markerIdx] + m_BWfromLong[aindLong]); //obtain allele which contrAllele get BW stutter from. Remember SHIFT
				contrVEC += xiB*colVEC( m_lookUpMat.col(1) ); //modify based on stutter products
			}
			if( m_FWfromLong[aindLong] > -1) {
				colVEC = contrMatPerAllele.col( m_startIndMarker_nAlleles[markerIdx] + m_FWfromLong[aindLong]); //obtain allele which contrAllele get BW stutter from. Remember SHIFT
				contrVEC += xiF*colVEC( m_lookUpMat.col(2) ); //modify based on stutter products							
			}
			
			//PREPARE WEIGHT EVIDENCE: 
			//Note:Obtaining the unique shape1 values is SLOW (hence removed)
			double weight; //Calculated weight for a given allele in allele set
			double peak; //temporary peak heigt 
			double shapeConst; //temporary variable for shape constnat
			int repIdx; //replicate variable
			Col <double> colVECcalc(m_NumCombOut1); //help variable (for calculation)
			for(int combOutIdx=0; combOutIdx < m_NumCombOut1; combOutIdx++) { //traverse each unique shape param
				double shape1 = contrVEC[combOutIdx]; //get unique shape param
				weight = 0; //reset weight again				
				if( cumAlleleIdx < 0) { //Assuming dropout for all replicates
					weight+= NumReps *(R::pgamma(AT, shape1,scale0, 1, 1)); //calculate for all replicates (identical)
				} else { //Peak heights are expected (but perhaps not for all replicates
					shapeConst = -lgamma(shape1) - shape1*const2; //calculate shape constant only one time
					for(repIdx=0; repIdx < NumReps; repIdx++) {	//Traverse all replicates (for given allele)
						peak = m_PeaksLong[startIndMarker_nAllelesReps + alleleIdx*NumReps + repIdx]; //obtain peak height wrt long-vector (takes into account replicates)
						if(peak < AT) {
							weight += R::pgamma(AT, shape1,scale0, 1, 1); //calculate for certain replicates
						} else {
							weight += shapeConst + (shape1-1)*log(peak) - peak * const1; //Note: lgamma is from std							
						}
					} //end for each replicate					
				} //end if expecting peak heights
				colVECcalc[combOutIdx] = weight;
			}
			m_EvidWeights.col(aindLong) = colVECcalc; //insert calculcated vector			
		} //END CALC		

		//Finally calculate per marker
		m_loglikMarkers.set_size(m_NumMarkers); 		
		double jointloglik = 0.0;				
		for(int m=0; m<m_NumMarkers; m++) {
			int fromCol = m_startIndMarker_nAllelesTot[m];
			int toCol = m_startIndMarker_nAllelesTot[m+1] - 1;
			//Rcpp::Rcout << "from=" << fromCol << " to=" << toCol << "\n";			
			m_loglikMarkers[m] = (EFMmarkerv[m]).calcLogLikMarker(m_EvidWeights.cols(fromCol,toCol) ) ; //insert marker information to vector			
			jointloglik += m_loglikMarkers[m];
		}
		return jointloglik; //return joint result
	}
	
	//exposed function to prepare indexing for given NOC
	//Here the allele-outcome matrices are prepared. Only depending on the number of contributors.
	//A contribution lookup-matrix is created to link unstuttered (contrMat0) and stuttered relation (i.e., contrMat1 which is not needed since "getContributionIndices" function already gives index)
	void prepare(int NOC) {
		m_NOC = NOC;
		uvec out0 = { 0, 1, 2 }; //defined outcomes for key0		
		int nOut0 = 3; //size of outcome (no stutter)
		int nOut1; //size of outcome (with possible stutters)
		int nCols; //number of columns used for additional stutters (1 + BW + FW)
		Mat<uword> out1; //nCols x nOutcome (transposed)
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
		Col<uword> contrMat0_key(m_NumCombOut0,fill::zeros);  //Key vector (no stutter)		
		
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
			Col<uword> keyVec(nCols,fill::zeros); //init key-vector as zero (numStutter models + 1)
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
			(EFMmarkerv[m]).prepareMarker(NOC, m_NumStutterModels[m]); 
		}
		//Rcpp::Rcout << m_EvidWeights << "\n";
	}

	//Constructor:
	EFMclass(vector<int>  NumStutterModels, int NumMarkers, vector<int> NumRepsMarkers, vector<int> NumAlleles, vector<int> startIndMarker_nAlleles, vector<int> startIndMarker_nAllelesReps, 
		vector<double> peaksLong, vector<double> freqsLong, vector<double> dropLong, 
		vector<double> nTypedLong, vector<double> maTypedLong, vector<double> basepairLong, 
		vector<int> BWfromLong, vector<int> FWfromLong, vector<int> BWtoLong, vector<int> FWtoLong,
		vector<int> NumPotStutters, vector<int> startIndMarker_nAllelesTot, vector<int> QalleleIndex, vector<double> pC, vector<double> fst, vector<double> AT,
		int NOK, vector<int> GindKnownLong,  vector<int> relGindLong, vector<double> ibdLong) {
				
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
			if( NumStutterModels[m] > m_NumStutterModelsMAX) m_NumStutterModelsMAX = NumStutterModels[m]; //set as highest observed
		}
		
		//m_basepairLong(basepairLong);
		m_basepairLong.set_size(TotalNumAlleles); //init
		for(int a=0; a < TotalNumAlleles; a++)  m_basepairLong[a] = basepairLong[a]; //copy to vector (need other vector type)
		
		//Prepare EFMmarker object
		for(int m=0; m< NumMarkers; m++) {
			EFMmarker *marker = new EFMmarker(m, NumRepsMarkers[m], NumAlleles[m], startIndMarker_nAlleles[m],startIndMarker_nAllelesReps[m], peaksLong, freqsLong, dropLong, nTypedLong[m], maTypedLong, BWtoLong, FWtoLong, NumPotStutters[m], startIndMarker_nAllelesTot[m], QalleleIndex[m], pC[m], fst[m], AT[m], NOK, GindKnownLong,  relGindLong[m], ibdLong ); //prepare new marker
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
			vector<double> pC, vector<double> fst, vector<double> AT, int NOK, vector<int> GindKnownLong, 
			vector<int> relGindLong, vector<double> ibdLong, int maxThreads) {
				
			m_maxThreads = maxThreads; //copy
			efmStorage = new EFMclass(NumStutterModels, NumMarkers, NumRepsMarkers, NumAlleles, startIndMarker_nAlleles, startIndMarker_nAllelesReps, peaksLong, freqsLong, dropLong, nTypedLong, maTypedLong, basepairLong, 
			BWfromLong, FWfromLong, BWtoLong, FWtoLong, NumPotStutters, startIndMarker_nAllelesTot, QalleleIndex, pC, fst, AT, NOK, GindKnownLong, relGindLong, ibdLong);			
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
		Col<double> loglikMarkers() { 
			if(efmStorage) {				
				return( efmStorage->m_loglikMarkers );
			} else {
				return NULL;
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

    // Get access to functions:
	.method("filldata", &ExposedClass::filldata , "Fill object with data")
	.method("prepare", &ExposedClass::prepare , "Prepare data for calculations")
	.method("loglik", &ExposedClass::loglik , "Calculate and obtain likelihood value")
	.method("logliki", &ExposedClass::loglikMarkers , "Obtain marker specific likelihood values")
	.method("close", &ExposedClass::close , "Free memory")
	;
}
