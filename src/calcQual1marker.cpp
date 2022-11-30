//SCRIPT FOR CALCULATING QUALITATIV MODEL (ONE MARKER)

#include <vector> //vector storage
#include <cmath> //includes lgamma
#include <RcppArmadillo.h> //necessary for functions find, intercept
#ifdef _OPENMP
#include <omp.h> //parallelization
#endif

using namespace arma;
using namespace std;

class recursiveClass {
	public: 

	double m_pD;
	double m_pC;
	double m_fst;
	int m_nUnknowns;
	int m_nAlleles;
	int m_nReps;
	vector<double> m_freq;		
	vector<int> m_nTyped;
	Row<double> m_contrAlleles;
	double m_permFactorConstant; //used as constant (numerator for permutation scale (logged))
	int m_nTotTyped; //count number of total typed
	Row<double> m_likvalVEC; //likelihood value
	Mat<int> m_hasEvidMat;
						
	//prepare for genotype combinations
	int m_nGenos;
	Mat<double> m_contrMat;
	Mat<int> m_alleleMat;   

	//helpfunction to get permutuation factor
	double getPermutationFactor(vector<int> genoJointVec) {
		//sort(genoJointVec); //Note: genoJointVec is already sorted
		double permFactor = m_permFactorConstant; //calculate permutation factor				
		int counter = 1; //init counter
		for(int k = 1; k < m_nUnknowns; k++) { //start traverse from 2nd element
			if( genoJointVec[k]==genoJointVec[k-1]) { //CHECK IF PREVIOUS GENOTYPE WAS EQUAL OR NOT
				counter++; //increment if identical as previous value									
				if(k >= (m_nUnknowns-1)) permFactor /= tgamma(counter+1); //calculate factor if last element
			} else {
				permFactor /= tgamma(counter+1); 
				counter = 1; //reset counter
			}	
		}
		return( permFactor );
	}

	//helpfunction to calculate genotype prob
	double getGenoProb(int genoIdx, int *ptr_nTyped, int *ptr_nTotTyped) {
		double genoProb = 1;
		for (int i = 0; i < 2; i++) { //traverse both alleles in genotype
			int alleleIdx = m_alleleMat(genoIdx,i); //obtain allele idx for genotype
			genoProb *=  (ptr_nTyped[alleleIdx]*m_fst + (1-m_fst)*m_freq[alleleIdx]) / (1+(ptr_nTotTyped[0]-1)*m_fst);
			ptr_nTyped[alleleIdx]++; //add counter
			ptr_nTotTyped[0]++;
		}
		if(m_alleleMat(genoIdx,0)!=m_alleleMat(genoIdx,1)) genoProb *= 2; //If heterygous
		return(genoProb);
	}
	
	//Obtain data likelihood (can be multiple replicates)
	double getEvidProb(Row<double> contrAlleles) {
		double likEvid = 1;
		for(int r=0; r< m_nReps; r++) {
			bool hasDropin = false;
			for(int a=0; a < m_nAlleles; a++) { //traverse each alllee
				if( contrAlleles(a) > 0.5 ) { //IF CONTRIBUTOR:
					double dropProb = pow(m_pD, contrAlleles(a)); //dropout contribution
					if(m_hasEvidMat(r,a) > 0.5) {
						likEvid *= (1 - dropProb);
					} else {
						likEvid *= dropProb;
					}
				} else if(m_hasEvidMat(r,a)> 0.5) { //IF NO CONTRIBUTOR
					hasDropin = true;
					likEvid *= m_pC*m_freq[a];
				}
			}
			if(!hasDropin) likEvid *= (1-m_pC); //no dropin
		}
		return(likEvid);
	}

	//DEFINING THE RECURSIVE FUNCTION
	void recFun(int genoIdx1, int contrIdx, vector<int> genoJointIdx_rec, double genoProb_rec, Row<double> contrAlleles_rec, vector<int> nTyped_rec, int nTotTyped_rec) {
		//genoIdx1: the index of 1st unknown genotype contributor (due to parallel process)
		for(int genoIdx=0; genoIdx<m_nGenos; genoIdx++ ) { //Traverse each genotype comb		
			if( (contrIdx>0) && (genoIdx > genoJointIdx_rec[contrIdx-1])) break;
			
			//Re-copy variables for each new genotype (to be updated
			double genoProb_upd = genoProb_rec;
			vector<int> genoJointIdx_upd = genoJointIdx_rec;
			Row<double> contrAlleles_upd = contrAlleles_rec;
			vector<int> nTyped_upd = nTyped_rec;
			int nTotTyped_upd = nTotTyped_rec;
				
			//Update variables
			genoJointIdx_upd[contrIdx] = genoIdx; //insert genotype contribution	
			contrAlleles_upd += m_contrMat.row(genoIdx); //update contribution
			//Rcpp::Rcout << contrIdx << "\n";
			
			//STEP 1 (ALWAYS): Calculate genotype prob
			genoProb_upd *= getGenoProb(genoIdx, &(nTyped_upd[0]), &nTotTyped_upd);
			
			//Post calculations: Traverse deeper or calculate the likelihood of evidince
			if( contrIdx < (m_nUnknowns-1) ) { //IF BEFORE LAST CONTRIBUTOR WE RECURSIVE DEEPER
				recFun(genoIdx1, contrIdx+1 , genoJointIdx_upd, genoProb_upd, contrAlleles_upd, nTyped_upd, nTotTyped_upd); 
			} else {					
				double likEvid = getEvidProb(contrAlleles_upd); //Obtain data likelihood						
				int permFactor = getPermutationFactor(genoJointIdx_upd); //Obtain permutation factor						
				m_likvalVEC[genoIdx1] += likEvid*genoProb_upd*permFactor; //Calculate final likelilhood
			}
		} //end for each genotype
	} //End recursive function
	
	//Main function for calculations:
	double calcLik() {
		#pragma omp parallel for
		for(int genoIdx1=0; genoIdx1 < m_nGenos; genoIdx1++ ) { //Traverse each genotype comb					
			double genoProb_upd = 1.0; //init new
			vector<int> genoJointIdx_upd(m_nUnknowns, 0); //init new
			Row<double> contrAlleles_upd = m_contrAlleles;
			vector<int> nTyped_upd = m_nTyped;
			int nTotTyped_upd = m_nTotTyped;
				
			//Update variables
			genoJointIdx_upd[0] = genoIdx1; //insert genotype contribution	
			contrAlleles_upd += m_contrMat.row(genoIdx1); //update contribution
			
			//STEP 1 (ALWAYS): Calculate genotype prob
			genoProb_upd *= getGenoProb(genoIdx1, &(nTyped_upd[0]), &nTotTyped_upd);
	
			if(m_nUnknowns > 1) {
				recFun(genoIdx1, 1, genoJointIdx_upd, genoProb_upd, contrAlleles_upd, nTyped_upd, nTotTyped_upd);
			} else {
				m_likvalVEC[genoIdx1] = genoProb_upd*getEvidProb(contrAlleles_upd); //Obtain data likelihood	
			}
		}
		double likval = sum(m_likvalVEC); //calculate sum		
		return likval;
	}

	recursiveClass(double *pD, double *pC, double *fst, int *nUnknowns,  int *nAlleles, int *nReps, int *contrAlleles,
					int *nTyped, double *freq,  int *hasEvidVec ) {						
		m_pD = *pD;
		m_pC = *pC;
		m_fst = *fst;
		m_nUnknowns = *nUnknowns;
		m_nAlleles = *nAlleles;
		m_nReps = *nReps;
		m_permFactorConstant = tgamma(m_nUnknowns+1); 	
		
		//Init genotype structure:
		m_nGenos = int(m_nAlleles*(m_nAlleles + 1) / 2); //get Genotype outcome		
		m_likvalVEC.zeros(m_nGenos); //init for parall computing		
		m_contrMat.zeros(m_nGenos,m_nAlleles);
		m_alleleMat.zeros(m_nGenos,2);
		m_contrAlleles.zeros(m_nAlleles);
		m_nTotTyped = 0;
		int cc = 0; //counter oveer all		
		for (int i = 0; i < m_nAlleles; i++) {
			for (int j = i; j < m_nAlleles; j++) { 
				m_alleleMat(cc,0) = i; //include index
				m_alleleMat(cc,1) = j; //include index
				m_contrMat(cc,i) += 1; //insert contr at allele i
				m_contrMat(cc,j) += 1; //insert contr at allele j
				cc++; //iterate to next genotype outcome
			}			
			if(contrAlleles[i]>0) m_contrAlleles(i) = contrAlleles[i]; //copy
			m_nTyped.push_back(nTyped[i]); //copy
			m_freq.push_back(freq[i]); //copy
			m_nTotTyped += nTyped[i]; 
		}
		
		//Need to store indices of alleles in replicates that are observed
		m_hasEvidMat.set_size(m_nReps,m_nAlleles);
		for(int r=0; r < m_nReps; r++) {		
			for(int a=0; a < m_nAlleles; a++) {
				m_hasEvidMat(r,a) = hasEvidVec[m_nAlleles*r + a]; //fill in
			}
		}
		//Rcpp::Rcout << m_hasEvidMat << "\n";
	} //Done constructor
	
};


extern "C" {
	void calcQual1marker(double *lik, double *pD, double *pC, double *fst, int *nUnknowns,  
		int *nAlleles, int *nReps, int *contrAlleles, int *nTyped, double *freq,  int *hasEvidVec) {

		recursiveClass recObj( pD, pC, fst, nUnknowns, nAlleles, nReps, contrAlleles, nTyped, freq, hasEvidVec);
		*lik = recObj.calcLik(); //calculcate
	}
} //end external

