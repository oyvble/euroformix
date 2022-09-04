//This script contains functions for calculating the likelihoods function.
//AUTHOR: Oyvind Bleka, April 2020
/*ABOUT:

Contains two main functions which are called from R: 
1) loglikgammaC: Calculating the likelihood value for model (given parameter values) 
2) cumvalgammaC: Calculating the cumulative probability of PH (instead of numerical integration)

Both functions first structure data in the EFMmarker class before calculating the 'large-sum' (utilizing parallelization) for each markers separately:
1) loglikgammaC calls one of three functions (depending on situation):
	a) calcLogLikGammaMarkerUnknown: 1 replicate, only unknown contributors, no relatedness, last allele is dropout, no stutter
	b) calcLogLikGammaMarkerUnknownStutter: 1 replicate, only unknown contributors, no relatedness, last allele is dropout, only stutter
	c) calcLogLikGammaMarker: All other situations than a or b (known contributors, relatedness, diallelic markers, replicates).	

2) cumvalgammaC calls 'calccumval' (all model types): Calculates the cumulative probabilities of AT (analytical threshold) and maxY (max obtainable PH) for each alleles separately (conditioned on the others)
*/

#include <vector> //vector storage
#include <cmath> //includes lgamma
#include <thread> //used to obtain number of logical processes
#include <Rmath.h> //includes pgamma
#ifdef _OPENMP
#include <omp.h> //parallelization
#endif

using namespace std;

class EFMmarker {
	private: 
  	const double m_dSmalltol = 1.0e-30; //a tiny number > 0 (avoiding zero roundoff errors)
  	int cc,i,j; //iterators
  	int m_nNumGenos1p; //number of genotypes for 1 contributor 
  	vector< vector<double> > m_vOutG1contr;
  	vector< vector<int> > m_vOutG1mat;

  	//DATA VARIABLES
  	int m_nNumRep; //number of replicates (for marker)
  	int m_nNumAlleles; //number of alleles in marker
  	vector<double> m_vIntensity; //Intensity vector (copied number)
  	vector<double> m_vFreq; //Frequency vector
  	vector<double> m_vBasePairs; //base pair information from the kitinfo (for each alleles). Used for degradation model
  	vector<int> m_vBackstutters,m_vForwardStutters; //Indices for Backward/forward stutter vector (indicating which alleles backward/forward stutter to whom)
  	vector<double> m_vAlleleCounts; //counting number of typed alleles of type a
  	double m_nNumTyped; //counting number of typed alleles (total)
  	int m_nNumPotentUnobsStutters; //number of unobserved potential stutters (BW+FW)
  	vector<int> m_vPotentBackStutters,m_vPotentForwardStutters; //Indices for Backward/Forward potential unseen stutter vector (indicating which alleles backward/forward stutter to whom)

  	//Known contributor variables	
  	vector<int> m_vKnownGenoIdx; //genotype index of known contributors (follows the m_nNOK)
  	vector<int> m_vKnownContribIdx; //which contributor index the knowns have
  	vector<int> m_vUnknownContribIdx; //which contributor index the unknown have
  	int m_nNOC; //number of contributuros
  	int m_nNOK; //number of known contributors
  	int m_nNOU; //number of unknown contributors
      //The constructor prepares the variables with necessary data (arguments in constructor)
  
  	//added for related individuals
  	//vector<int> kindRel; //which contributor index the related individual has
  	int m_nUnknownContribRelIdx = -1; //(unknown) contributor index of related (k=1,..,m_nNOU)
  	int m_nRelGenoIdx = -1; //genotype index of related
	
	public: 
	  vector< vector<double> > cumprobvals; //m_nNumAlleles x 2 matrix with cumulative probabilities (used in function	calccumval)
	
	EFMmarker(int *NOC2, 
           int *NOK2, 
           int *knownGind, 
           int *nR, 
           int *nA, 
           double *peaksLong, 
           double *freqsLong, 
           double *m_nNumTypedLong, 
           double *maTypedLong, 
           double *basepairLong, 
           int *BWstuttindLong, 
           int *FWstuttindLong, 
           int *nPS, 
           int *BWPstuttindLong, 
           int *FWPstuttindLong, 
           int *relGind) { 
		//NOC2 = number of contrs (known + unknown)
		//knownGind = m_nNOC long vector with genotype index of known contributors
		//nR = number of replicates (needed to traverse peak height data)
		//nA = number of alleles. NB: Assuming same length for all replicates!!!!
		//peaksLong = longvector of Intensity values
		//freqsLong = longvector of frequency values
		//m_nNumTypedvec = longvector number of typed alleles (total)
		//maTypedLong = longvector of number of typed alleles of type a
		//basepairLong = vectorized vector of the already adjusted base pairs (x-150)/100 
		//BWstuttindLong = vectorized vector of backward stutter indices (gives index of what allele it receive stutter from  (Allele 1=index0))
		//FWstuttindLong = vectorized vector of forward stutter indices (gives index of what allele it receive stutter from  (Allele 1=index0))
		//nPS = number of potential stutters per marker
		//BWPstuttindLong = vectorized vector of backward stutter indices (gives index of what allele it receive stutter from  (Allele 1=index0))
		//FWPstuttindLong = vectorized vector of forward stutter indices (gives index of what allele it receive stutter from  (Allele 1=index0))
		//relGind = m_nNOC long vector with genotype index of related contributors
		
		m_nNOC = *NOC2; //copy number of contributors
		m_nNOK = *NOK2; //copy number of contributors
		m_nNOU = m_nNOC - m_nNOK; //number of unknowns
		m_nNumRep = *nR; //set number of replicates
		m_nNumAlleles = *nA; //set number of alleles
		m_nNumPotentUnobsStutters = *nPS; //set number of potential stutters (must be gloabal)
		m_nNumTyped = *m_nNumTypedLong; //set number of total prev. typed
		m_vIntensity.assign(m_nNumRep*m_nNumAlleles,0.0); //init PH-vector
		m_vFreq.assign(m_nNumAlleles,0.0); //init vector
		m_vBackstutters.assign(m_nNumAlleles,0); //init vector
		m_vForwardStutters.assign(m_nNumAlleles,0); //init vector
		m_vAlleleCounts.assign(m_nNumAlleles,0); //init vector
		m_vBasePairs.assign(m_nNumAlleles,0); //init vector
		for(i=0; i<m_nNumAlleles; i++) { //for each alleles (observed + dropout)
			m_vIntensity[i] = peaksLong[i]; //copy value
			m_vFreq[i] = freqsLong[i]; //copy value
			m_vBackstutters[i] = BWstuttindLong[i]; //copy value
			m_vForwardStutters[i] = FWstuttindLong[i]; //copy value
			m_vAlleleCounts[i] = maTypedLong[i]; //copy value
			m_vBasePairs[i] = basepairLong[i]; //copy value
		}
		for(i=m_nNumAlleles; i<(m_nNumRep*m_nNumAlleles); i++) { //include remaining PHs
			m_vIntensity[i] = peaksLong[i]; //copy value
		}
		
		m_vPotentBackStutters.assign(m_nNumPotentUnobsStutters,0); //init vector
		m_vPotentForwardStutters.assign(m_nNumPotentUnobsStutters,0); //init vector
		for(i=0; i<m_nNumPotentUnobsStutters; i++) { //for each unobserved potential stutters (BW+FW)
			m_vPotentBackStutters[i] = BWPstuttindLong[i]; //copy value
			m_vPotentForwardStutters[i] = FWPstuttindLong[i]; //copy value
		}		
		
		//Create contribution matrix (1 genotype): Could be as part of constructor to save time
		m_nNumGenos1p = int(m_nNumAlleles*(m_nNumAlleles + 1) / 2); //get Genotype outcome
		m_vOutG1contr.assign(m_nNumGenos1p, vector<double>(m_nNumAlleles, 0.0)); //init nG1xnA matrix (contribution matrix). Indicating what alleles that are contributoed
		m_vOutG1mat.assign(m_nNumGenos1p, vector<int>(2, 0)); //init nG1x2 matrix (allele names as indices 0,...,nA-1)
		cc = 0; //counter oveer all
		for (i = 0; i < m_nNumAlleles; i++) {
			for (j = i; j < m_nNumAlleles; j++) { 
				m_vOutG1mat[cc][0] = i; //include index
				m_vOutG1mat[cc][1] = j; //include index
				m_vOutG1contr[cc][i] += 1.0; //insert contr at allele i
				m_vOutG1contr[cc][j] += 1.0; //insert contr at allele j
				cc++; //iterate to next genotype outcome
			}
		}
		
		//Prepare vector for known contributors (and also unknown): Need to know positions!
		m_vKnownGenoIdx.assign(m_nNOK,0); //genotype index in vector 
		m_vKnownContribIdx.assign(m_nNOK,0); //contributor index in vector
		m_vUnknownContribIdx.assign(m_nNOU,-1); //contributor index in vector
		//kindRel.assign(m_nNOC,-1); //genotype index in vector, one index for each contributors  (-1 means no related)
		
		cc = 0; //counters for known
		j = 0; //counter for unknowns
		for(i=0; i<m_nNOC; i++) { //for each contributors:
			if(knownGind[i]>=0) { //If contributor is known (genotype given)
				m_vKnownGenoIdx[cc] = knownGind[i]; //copy genotype index
				m_vKnownContribIdx[cc] = i; //insert contributor index
				cc++; //update counter for knowns
			} else { //if contributor is unknown (genotype not given)
				if(relGind[i] != -1 ) { //in case of related
					m_nUnknownContribRelIdx = j; //insert contributor index of unknowns
					m_nRelGenoIdx = relGind[i]; //insert genotype index
				}					
				m_vUnknownContribIdx[j] = i; //insert contributor index of unknown. 
				j++; //update counter for unknowns
				
			}
		}		
	} //end constructor
	
		
	//////////////////////////////////////////////////////////
	//HELPFUNCTION FOR CALCULATING RELATEDNESS PROBABILITIES//
	//////////////////////////////////////////////////////////
	
	//helpfunction to get allele probability:
	double prob_a(double Pa, double mm, double nn, double *fst2) {
		return( (mm*(*fst2) + (1-(*fst2))*Pa)/(1+(nn-1)*(*fst2)) ); 
	}

	double prob_relUnknown(int Ugind, int Rgind, double *ibd2, double *fst2, double *m_vAlleleCounts2, double *m_nNumTyped2) {
		//Ugind = genotype index of unknown indiviudal 
		//Rgind = genotype index of related indiviudal 		
		//ibd2 = ibd vector (m_nNOC long)
		double genoSum; //used to sum the genotype probability tfor unknowns		
		int aind = m_vOutG1mat[Ugind][0]; //get allele index of genotype g_1
		int bind = m_vOutG1mat[Ugind][1]; //get allele index of genotype g_2
		bool Uhom = aind==bind; //boolean of whether unknown genotype is homozygote
		
		//First step: Calculate random match probability of unrelated situation:
		genoSum = prob_a(m_vFreq[aind],m_vAlleleCounts2[aind],*m_nNumTyped2,fst2); //init with prob 1st allele  
		if(Uhom) { //if unknown is homozygote					
			genoSum *= prob_a(m_vFreq[aind],m_vAlleleCounts2[aind]+1,*m_nNumTyped2+1,fst2); //calculate random match prob (always used) 					
		} else { //if unknown is heterozygote variant
			genoSum *= 2*prob_a(m_vFreq[bind],m_vAlleleCounts2[bind],*m_nNumTyped2+1,fst2); //calculate prob 2st allele  (and scale with 2)
		}
			
		//Extension with kappa-coefficient: SEE FORMULAS IN TABLE A.3 in Book "A forensic practicioners guide...."
		if( Rgind != -1 ) { //if related is specified (not -1 index)
			genoSum *= ibd2[0]; //multiply with kappa0
			if( Ugind==Rgind ) { //if Unknown and Related are same genotype
				genoSum+=ibd2[2]; //sum with kappa2
				
				if(Uhom) { //if unknown is homozygote
					genoSum+= prob_a(m_vFreq[aind],m_vAlleleCounts2[aind],*m_nNumTyped2,fst2)*ibd2[1] ; //multiply with kappa1 
				} else { //if unknown is heterozygote variant
					genoSum+= (prob_a(m_vFreq[aind],m_vAlleleCounts2[aind],*m_nNumTyped2,fst2)+prob_a(m_vFreq[bind],m_vAlleleCounts2[bind],*m_nNumTyped2,fst2))*ibd2[1]/2 ; //multiply with kappa1 						
				}	
				
			} else { //if not the same genotype we need to check overlap (a,b)~(c,d)
				bool A1eq = aind==m_vOutG1mat[Rgind][0]; //check if a=c 
				bool A2eq = aind==m_vOutG1mat[Rgind][1]; //check if a=d 
				bool B1eq = bind==m_vOutG1mat[Rgind][0]; //check if b=c 
				bool B2eq = bind==m_vOutG1mat[Rgind][1]; //check if b=d 

				if( A1eq || A2eq || B1eq || B2eq) { //if any overlap (one shared allele): THERE ARE 3 OUTCOME!!
					if(Uhom) { //if Unknown is homozygote we know what allele to use
						genoSum += prob_a(m_vFreq[aind],m_vAlleleCounts2[aind],*m_nNumTyped2,fst2)*ibd2[1]/2; //multiply with kappa1 	
					} else { //if Unknown is heterozygote we need to found the non-overlapping allele 
						bool Rhom = m_vOutG1mat[Rgind][0]==m_vOutG1mat[Rgind][1]; //check if related genotype is homozygote
						int cind; //temporary variable to select non-overlapping allele
						if(A1eq || A2eq) { //Allele bind is non-overlapping
							cind = bind;
						} else {
							cind = aind; //Allele aind is non-overlapping
						}
						if(Rhom) { //should not divide by 2 if related geno is homozygote
							genoSum += prob_a(m_vFreq[cind],m_vAlleleCounts2[cind],*m_nNumTyped2,fst2)*ibd2[1]; //multiply with kappa1 																			
						} else { //should divide by 2 if related geno is heterozygote
							genoSum += prob_a(m_vFreq[cind],m_vAlleleCounts2[cind],*m_nNumTyped2,fst2)*ibd2[1]/2; //multiply with kappa1 																												
						}
					}
				}	//otherwise none shared allele (and we are done (kappa0 already included)												
			}
		}
		return( genoSum ); //return genotype prob
	}
	
	
	//Engine to calculate the log-likelihood function
	//LIKELIHOOD FUNCTION (POSSIBLE KNOWN CONTRIBUTORS AND RELATED INDIVIDUALS)
	double calcLogLikGammaMarker(double *fst,double *AT, double *lambda, double *prC, double *mixprop, double *mu, double *sigma, double *beta, double *xiB0, double *xiF0, double *ibd) {

		//MODEL PARAMETER VARIABLES:
		//fst = sub-population correction
		//AT = Analytical threshold (assumed constant per markers)
		//mixprop = mixture proportions (length must be *m_nNOU)
		//shape = 1/sig^2, where sigma=coefficient of PH amount
		//scale = mu*sig^2, where mu=expected PH amount
		//beta is degradation slope parameters
		//xiB0 is baseline of expected stutter proportion (backward stutter)
		//xiF0 is baseline of expected stutter proportion (forward stutter)
	    //ibd = kappa-coefficient for each contributors (vectorised)
		int	nGK = pow(m_nNumGenos1p, m_nNOU); //get total number of combined genotypes (for unknown contrs)
		double invshape = (*sigma)*(*sigma);//get inverse of shape parameter
		double shape = 1/invshape; //get shape parameter
		double scale = (*mu)*invshape; //get scale parameter
		const double const1 = 1/(scale); //constant 1
	 	const double const2 = log(scale); //constant 2
			
		//Prepare variables:
		int cind,a,r;
		vector<double> shapevK(m_nNumAlleles, 0.0); //Precalculation for known contributors:
		vector<double> shapev0(m_nNumAlleles, 0.0); //degrad scaling of shape
		vector<double> dropin0(m_nNumAlleles*m_nNumRep, 0.0); //drop-in contribution per peak height
		vector<double> constY0(m_nNumAlleles*m_nNumRep, 0.0); //const1 times m_vIntensity	
		for (a = 0; a < m_nNumAlleles; a++) { //traverse each observed alleles (indices), also those PH=0
			shapev0[a] = exp( log(shape) + m_vBasePairs[a]*log(*beta) ); //scaling with degradation model (assumed already scaled)
			for(r = 0; r < m_nNumRep; r++) { //traverse each replicates (observed alleles indicated by PH)
				cind = a*m_nNumRep + r; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
				dropin0[cind] = log(*lambda) - (*lambda)*(m_vIntensity[cind]- (*AT) ) + log(*prC) + log(m_vFreq[a]); //add to dropin prob to sum	
				constY0[cind] = m_vIntensity[cind]*const1; 
			}
		}

		//Sum up contribution for each alleles (Taking into account mix proportions): ONLY CALCULATED FOR KNOWN CONTRIBUTORS INITIALLY
		int aa, kk;
		for (kk = 0; kk < m_nNOK; kk++) { //for each known contributors 
			for (aa = 0; aa < m_nNumAlleles; aa++) { //traverse each alleles (indices)
				shapevK[aa] += m_vOutG1contr[m_vKnownGenoIdx[kk]][aa] * mixprop[m_vKnownContribIdx[kk]]; //contr from contr k to allele a
			}
		}
		//CALCULATING LARGE SUM:
		double bigsum = 0.0;  //total sum over all genotypes
		
		#pragma omp parallel for reduction(+:bigsum) //shared(shapev0, dropin0,constY0) //perform parallel on summing up bigsum
		for (int largeIter = 0; largeIter < nGK; largeIter++) { //for each combined/joint genotype outcome

			//Following variables must be declared for each iteration (non-shared):
			int a,k,r; //used to traverse alleles(a), contributors(k) and replicates(r)
			int aind; //used as index for alleles in genotypes
			vector<int> jointGind(m_nNOU, 0); //This will be the permuation contribution index for the unknown inds (directly corresponds to indices of Gmarg)
			vector<double> shapev = shapevK; //make a copy of existing shapevector
			vector<double> m_vAlleleCounts2 = m_vAlleleCounts; //Creating copy of counter (per allele)
			double m_nNumTyped2 = m_nNumTyped; //creating copy of total counter

			//jointGind = digits(largeIter, base = m_nNumGenos1p, pad = m_nNOU); #Equivalent operaion
			double genoProd = 1.0; //calculating the genotype probability of the unknowns
			int modrest = largeIter; //used to keep remained after modulo (init as iter number)
			bool inserted = false; //boolean of whether all digits are inserted (Rest is zero padded)
			for (k = 0; k < m_nNOU; k++) { //for each unknown contributors (summing up wrt both contr (m_vOutG1contr) and mx (mixprop)): Need each contr to derive shapev
				if (!inserted) { //if not all digits inserted
					if ( k>0 ) {
						modrest = int((modrest - jointGind[k - 1]) / m_nNumGenos1p); //extract remaining, divide to get to next digit (necessary with int converion?)
					}
					jointGind[k] = modrest % m_nNumGenos1p; //INSERT NUMBER: convert number to "m_nNumGenos1p" basis 	
					if (modrest < m_nNumGenos1p) { //check if rest is smaller than base (only run if not inserted)
						inserted = true;  //then all digits are inserted
					}	
				} //else { jointGind[k] = 0; //zero pad	}
		
				//Sum up contribution for each alleles (Taking into account mix proportions)
				for (a = 0; a < m_nNumAlleles; a++) { //traverse each alleles (indices)
					shapev[a] += m_vOutG1contr[jointGind[k]][a] * mixprop[m_vUnknownContribIdx[k]]; //contr from contr k to allele a. NOTICE THE k+*m_nNOK shift!
				}
								
				//CALCULATE GENOTYPE PROBS OF UNKNOWNS (MAY BE RELATED != -1)
				if( m_nUnknownContribRelIdx != k ) { //if unknown is not related
					for(a = 0;a<2; a++) {
						aind = m_vOutG1mat[jointGind[k]][a]; //get allele index of genotype g_a
						genoProd *= ((*fst)*m_vAlleleCounts2[aind] + (1-*fst)*m_vFreq[aind]) / (1 + (m_nNumTyped2-1)*(*fst)); 
						m_vAlleleCounts2[aind] += 1; //update allele count for particular genotype
						m_nNumTyped2 += 1; //update total count
					}
					if(m_vOutG1mat[jointGind[k]][0]!=m_vOutG1mat[jointGind[k]][1]) { // heterozygous variant
						genoProd *=2; //multiply by 2  if het. variant
					} 
				} 
			} //end for each unknown unrelated contributors

			//AN UNKNOWN RELATED IS CALCULATED LAST (AFTER ALLELE COUNT UPDATES)
			if( m_nUnknownContribRelIdx != -1 ) { //if an unknown was related (then k=m_nUnknownContribRelIdx)
				genoProd *= prob_relUnknown( jointGind[m_nUnknownContribRelIdx], m_nRelGenoIdx, (ibd + 3*m_vUnknownContribIdx[m_nUnknownContribRelIdx]), fst, &(m_vAlleleCounts2[0]), &m_nNumTyped2);	 //Note: scale with 3 because ibd is a '3-long vector' per contributor						
			}

			//////////////////////////////////////////////
			//Calculating the inner sum-part begins here//
			//////////////////////////////////////////////

			//Scaling shape parameter with degrad model
			for (a = 0; a < m_nNumAlleles; a++) { //traverse each observed alleles (indices), having PH>0
				shapev[a] *= shapev0[a]; //scaling with degradation model (assumed already scaled)
			}
			//Scaling shape parameter with stutter model:
			vector<double> shapev2 = shapev;  //make a copy of existing shapevector
			double shape1; //modified shape param (caused by stutters)
			
			//Stutter model calculations (for given allele)
			if( (*xiB0>m_dSmalltol) || (*xiF0>m_dSmalltol) ) { //only run if stutterparam is positive
				for (a = 0; a < (m_nNumAlleles-1); a++) { //traverse each alleles. Potential non-Qallele is traversed below
					shape1 = 0.0; //copy original contribution
										
					//CONSIDERING STUTTERS: Always running (with param equal zero)
					if(m_vForwardStutters[a] != -1) { //if a-1 stutter of a exists (found in FWstutter info)
						shape1 = (*xiF0)*shapev[m_vForwardStutters[a]]; //proportion obtained from allele a-1
					} 
					if(m_vBackstutters[a] != -1) { //if a+1 stutter of a exists (found in BWstutter info)
						shape1 += (*xiB0)*shapev[m_vBackstutters[a]]; //proportion obtained from allele a-1
					}
					shapev2[a] = (1- *xiB0 - *xiF0)*shapev[a] + shape1; //updated shape parameter (summing fraction lost and attained), proportion xiB/xiF always lost (possible to unseen potential stutters see later)
				} 
				
				//HANDLING SPECIAL SCENARIO WHEN LAST ALLELE IS NOT Q-allele
				a = (m_nNumAlleles-1); //set allele index as last allele
				if(m_vForwardStutters[a] != -1 || m_vBackstutters[a] != -1) {  //IF THERE WAS STUTTER INDEX DEFINED
					shape1 = 0.0; //copy original contribution									
					if(m_vForwardStutters[a] != -1) { shape1 = (*xiF0)*shapev[m_vForwardStutters[a]]; } //proportion obtained from allele a-1
					if(m_vBackstutters[a] != -1) { shape1 += (*xiB0)*shapev[m_vBackstutters[a]]; } //proportion obtained from allele a+1
					shapev2[a] = (1- *xiB0 - *xiF0)*shapev[a] + shape1; //updated shape parameter (summing fraction lost and attained), proportion xiB/xiF always lost (possible to unseen potential stutters see later)				
				} //end special scenario
			}
			
			//Summing up contribution of each alleles:
			vector<bool> anyDropin(m_nNumRep,false); //use vector to indicate drop-in
			double logevidProb = 0.0; //evidence probability (weight)
			int cind; //cumulative allele indexing for PHs (for traversing over all replicates)
			for (a = 0; a < m_nNumAlleles; a++) { //traverse each observed alleles (also consider the last allele since Q-allele may have PH last allee)
				for(r = 0; r < m_nNumRep; r++) { //traverse each replicates (observed alleles indicated by PH)
					cind = a*m_nNumRep + r; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
					if( m_vIntensity[cind] > m_dSmalltol ){ //if PH>0 (this is same as PH>=AT since threshold has already been applied)
						if(shapev2[a] > m_dSmalltol) { //If contribution and PH>0  					//CONTRIBUTION SET (A)
							logevidProb += -lgamma(shapev2[a]) - shapev2[a]*const2 + (shapev2[a]-1)*log(m_vIntensity[cind]) - m_vIntensity[cind]*const1; //gamma/lgamma is found in cmath							
						} else { //If contribution and PH>0 	//DROPIN SET (C)
							anyDropin[r] = true; //there was at least one dropin
							//evidProb *= lambda*exp(lambda*(m_vIntensity[a]-AT))*prC*m_vFreq[a]; //add to dropin prob to sum			
							logevidProb += dropin0[cind]; //log(*lambda) - (*lambda)*(m_vIntensity[cind] - (*AT) ) + log(*prC) + log(m_vFreq[a]); //add to dropin prob to sum										
						}							
					} else if(shapev2[a] > m_dSmalltol) { //IF PH missing  and contribution(it't a dropout //CONTRIBUTION SET (B)						
						logevidProb += pgamma(*AT, shapev2[a],scale, 1, 1); //log( gamma_p(shapev2[a],(*AT)*const1) );  //Add log(dropout-probability)							
					}
				} //end for each replicates (r)
			} //end for each observed alleles
			
			//Need to add this block to calc unobserved potential stutters (done ones for all replicates, since this can be scaled)
			if( (*xiB0>m_dSmalltol) || (*xiF0>m_dSmalltol) ) { //only run if param is positive
				for (a = 0; a < m_nNumPotentUnobsStutters; a++) { //traverse for each unobserved potential stutters
					shape1 = 0.0; //set to zero by default
					if(m_vPotentForwardStutters[a] != -1) { //if a-1 stutter of a exists (found in FWPstutter info)
						shape1 = (*xiF0)*shapev[m_vPotentForwardStutters[a]]; //proportion obtained from allele a-1
					} 
					if(m_vPotentBackStutters[a] != -1) { //if a+1 stutter of a exists (found in BWPstutter info)
						shape1 += (*xiB0)*shapev[m_vPotentBackStutters[a]]; //proportion obtained from allele a+1
					}
					if( shape1>m_dSmalltol) { //if there was a contribution, we need to calculate the prob for dropout						
							logevidProb += m_nNumRep*pgamma(*AT, shape1,scale, 1, 1); //log( gamma_p(shape1,(*AT)*const1) ); //scaling with number replicates
					}
				}
			}
			
			//Liklihood NO DROP-IN: Traverse each replicates (observed alleles indicated by PH)
			for(r = 0; r<m_nNumRep; r++) {
				if(!anyDropin[r]) { //if no drop-in occured for a replicate
					logevidProb += log(1-*prC); //
				}
			}
			bigsum += exp(logevidProb)*genoProd; //add temp res to big sum
		} //end for each combinations (bigsum)
	
		return log(bigsum); //return function
	}
	
	
	double calcLogLikGammaMarkerUnknown(double *fst,double *AT, double *lambda, double *prC, double *mixprop, double *mu, double *sigma, double *beta) {
		int	nGK = pow(m_nNumGenos1p, m_nNOU); //get total number of combined genotypes (for unknown contrs)
		double invshape = (*sigma)*(*sigma);//get inverse of shape parameter
		double shape = 1/invshape; //get shape parameter
		double scale = (*mu)*invshape; //get scale parameter
		const double const1 = 1/(scale); //constant 1
	 	const double const2 = log(scale); //constant 2
			
		//Precalculation for known contributors:
		vector<double> shapevK(m_nNumAlleles, 0.0); //contribution vector of each alleles after taking into account mx

		//Prepare variables not needed in parallel:
		vector<double> shapev0(m_nNumAlleles, 0.0); //contribution vector of each alleles after taking into account mx
		vector<double> dropin0(m_nNumAlleles, 0.0); //drop-in contribution per peak height
		vector<double> constY0(m_nNumAlleles, 0.0); //PH vector multiplied with constant
		for (int a = 0; a < m_nNumAlleles; a++) { //traverse each observed alleles (indices), having PH>0
			shapev0[a] = exp( log(shape) + m_vBasePairs[a]*log(*beta) ); //scaling with degradation model (assumed already scaled)
			dropin0[a] = log(*lambda) - (*lambda)*(m_vIntensity[a]- (*AT) ) + log(*prC) + log(m_vFreq[a]); //add to dropin prob to sum	
			constY0[a] = m_vIntensity[a]*const1; 
		}

		//CALCULATING LARGE SUM:
		double bigsum = 0.0;  //total sum over all genotypes
		
		#pragma omp parallel for reduction(+:bigsum) //shared(shapev0, dropin0,constY0) //perform parallel on summing up bigsum
		for (int largeIter = 0; largeIter < nGK; largeIter++) { //for each combined/joint genotype outcome

			//Following variables must be declared for each iteration (non-shared):
			int a; //used to traverse alleles
			int k; //used to traverse contributors
			int aind; //used as index for alleles in genotypes
			vector<int> jointGind(m_nNOU, 0); //This will be the permuation part (directly corresponds to indices of Gmarg)
			vector<double> shapev(m_nNumAlleles, 0.0); //contribution vector of each alleles after taking into account mx
			vector<double> m_vAlleleCounts2 = m_vAlleleCounts; //Creating copy of counter (per allele)
			double m_nNumTyped2 = m_nNumTyped; //creating copy of total counter

			//jointGind = digits(largeIter, base = m_nNumGenos1p, pad = m_nNOC); #Equivalent operaion
			double genoProd = 1.0; //calculating the genotype probability of the unknowns
			int modrest = largeIter; //used to keep remained after modulo (init as iter number)
			bool inserted = false; //boolean of whether all digits are inserted (Rest is zero padded)
			for (k = 0; k < m_nNOU; k++) { //for each contributor (summing up wrt both contr (m_vOutG1contr) and mx (mixprop)
				if (!inserted) { //if not all digits inserted
					if ( k>0 ) {
						modrest = int((modrest - jointGind[k - 1]) / m_nNumGenos1p); //extract remaining, divide to get to next digit (necessary with int converion?)
					}
					jointGind[k] = modrest % m_nNumGenos1p; //INSERT NUMBER: convert number to "m_nNumGenos1p" basis 						
					if (modrest < m_nNumGenos1p) { //check if rest is smaller than base (only run if not inserted)
						inserted = true;  //then all digits are inserted
					}	
				} //else { jointGind[k] = 0; //zero pad	}

				//Sum up contribution for each alleles (Taking into account mix proportions)
				for (a = 0; a < m_nNumAlleles; a++) { //traverse each alleles (indices)
					shapev[a] += m_vOutG1contr[jointGind[k]][a] * mixprop[k]; //contr from contr k to allele a
				}
				
				//Calculing the genotype probababilities where fst is taken into account (Balding-Nicholsen correction): (fst*ma  + (1-fst)*pa)/(1+ (n-1)*fst)
				for(a = 0;a<2; a++) {
					aind = m_vOutG1mat[jointGind[k]][a]; //get allele index of genotype g_a
					genoProd *= ((*fst)*m_vAlleleCounts2[aind] + (1-*fst)*m_vFreq[aind]) / (1 + (m_nNumTyped2-1)*(*fst)); 
					m_vAlleleCounts2[aind] += 1; //update allele count
					m_nNumTyped2 += 1; //update total count
				}

				if(m_vOutG1mat[jointGind[k]][0]!=m_vOutG1mat[jointGind[k]][1]) { // heterozygous variant
					genoProd *=2; //multiply by 2  if het. variant
				}
			} //end for each contributor

			//Calculating the inner sum-part begins here:
			//Scaling shape parameter with degrad model:
			for (a = 0; a < m_nNumAlleles; a++) { //traverse each observed alleles (indices), having PH>0
				shapev[a] *= shapev0[a]; //scaling with degradation model (assumed already scaled)
			}
			
			//Summing up contribution of each alleles:
			bool anyDropin = false;
			double logevidProb = 0.0; //evidence probability (weight)
			double shape1; //modified shape param (caused by stutters)
			for (a = 0; a < (m_nNumAlleles-1); a++) { //traverse each observed alleles (indices), having PH>0 (no one are dropouts)				
				shape1 = shapev[a]; //copy original contribution	
				if( shape1>m_dSmalltol ) { //If contribution 
					logevidProb += -lgamma(shape1) - shape1*const2 + (shape1-1)*log(m_vIntensity[a]) - constY0[a]; //gamma/lgamma is found in cmath							
				} else { //no contribution and PH present
					anyDropin = true; //there was at least one dropin
					logevidProb += dropin0[a]; //add to dropin prob to sum			
				} 
			} //end for each observed alleles
						
			//Calculate for the dropout allele (also including potential stutters:
			if(shapev[m_nNumAlleles-1]>m_dSmalltol) { //only if any contributors
					logevidProb += pgamma(*AT, shapev[m_nNumAlleles-1],scale, 1, 1); //log( gamma_p(shapev[m_nNumAlleles-1],(*AT)*const1) ); 
			}
			
			if(!anyDropin) { //if no drop-in occured
				logevidProb += log(1-*prC); //
			}

			//bigsum += evidProb*genoProd; //add temp res to big sum
			bigsum += exp(logevidProb)*genoProd; //add temp res to big sum
		} //end for each combinations (bigsum)	
		return log(bigsum); //return function
	}
	
	double calcLogLikGammaMarkerUnknownStutter(double *fst,double *AT, double *lambda, double *prC, double *mixprop, double *mu, double *sigma, double *beta, double *xiB0, double *xiF0) {
		int	nGK = pow(m_nNumGenos1p, m_nNOU); //get total number of combined genotypes (for unknown contrs)
		double invshape = (*sigma)*(*sigma);//get inverse of shape parameter
		double shape = 1/invshape; //get shape parameter
		double scale = (*mu)*invshape; //get scale parameter
		const double const1 = 1/(scale); //constant 1
	 	const double const2 = log(scale); //constant 2
			
		//Precalculation for known contributors:
		vector<double> shapevK(m_nNumAlleles, 0.0); //contribution vector of each alleles after taking into account mx

		//Prepare variables not needed in parallel:
		vector<double> shapev0(m_nNumAlleles, 0.0); //contribution vector of each alleles after taking into account mx
		vector<double> dropin0(m_nNumAlleles, 0.0); //drop-in contribution per peak height
		vector<double> constY0(m_nNumAlleles, 0.0); //PH vector multiplied with constant
		for (int a = 0; a < m_nNumAlleles; a++) { //traverse each observed alleles (indices), having PH>0
			shapev0[a] = exp( log(shape) + m_vBasePairs[a]*log(*beta) ); //scaling with degradation model (assumed already scaled)
			dropin0[a] = log(*lambda) - (*lambda)*(m_vIntensity[a]- (*AT) ) + log(*prC) + log(m_vFreq[a]); //add to dropin prob to sum	
			constY0[a] = m_vIntensity[a]*const1; 
		}

		//CALCULATING LARGE SUM:
		double bigsum = 0.0;  //total sum over all genotypes
		
		#pragma omp parallel for reduction(+:bigsum) //shared(shapev0, dropin0,constY0) //perform parallel on summing up bigsum
		for (int largeIter = 0; largeIter < nGK; largeIter++) { //for each combined/joint genotype outcome

			//Following variables must be declared for each iteration (non-shared):
			int a; //used to traverse alleles
			int k; //used to traverse contributors
			int aind; //used as index for alleles in genotypes
			vector<int> jointGind(m_nNOU, 0); //This will be the permuation part (directly corresponds to indices of Gmarg)
			vector<double> shapev(m_nNumAlleles, 0.0); //contribution vector of each alleles after taking into account mx
			vector<double> m_vAlleleCounts2 = m_vAlleleCounts; //Creating copy of counter (per allele)
			double m_nNumTyped2 = m_nNumTyped; //creating copy of total counter

			//jointGind = digits(largeIter, base = m_nNumGenos1p, pad = m_nNOC); #Equivalent operaion
			double genoProd = 1.0; //calculating the genotype probability of the unknowns
			int modrest = largeIter; //used to keep remained after modulo (init as iter number)
			bool inserted = false; //boolean of whether all digits are inserted (Rest is zero padded)
			for (k = 0; k < m_nNOU; k++) { //for each contributor (summing up wrt both contr (m_vOutG1contr) and mx (mixprop)
				if (!inserted) { //if not all digits inserted
					if ( k>0 ) {
						modrest = int((modrest - jointGind[k - 1]) / m_nNumGenos1p); //extract remaining, divide to get to next digit (necessary with int converion?)
					}
					jointGind[k] = modrest % m_nNumGenos1p; //INSERT NUMBER: convert number to "m_nNumGenos1p" basis 						
					if (modrest < m_nNumGenos1p) { //check if rest is smaller than base (only run if not inserted)
						inserted = true;  //then all digits are inserted
					}	
				} //else { jointGind[k] = 0; //zero pad	}

				//Sum up contribution for each alleles (Taking into account mix proportions)
				for (a = 0; a < m_nNumAlleles; a++) { //traverse each alleles (indices)
					shapev[a] += m_vOutG1contr[jointGind[k]][a] * mixprop[k]; //contr from contr k to allele a
				}
				
				//Calculing the genotype probababilities where fst is taken into account (Balding-Nicholsen correction): (fst*ma  + (1-fst)*pa)/(1+ (n-1)*fst)
				for(a = 0;a<2; a++) {
					aind = m_vOutG1mat[jointGind[k]][a]; //get allele index of genotype g_a
					genoProd *= ((*fst)*m_vAlleleCounts2[aind] + (1-*fst)*m_vFreq[aind]) / (1 + (m_nNumTyped2-1)*(*fst)); 
					m_vAlleleCounts2[aind] += 1; //update allele count
					m_nNumTyped2 += 1; //update total count
				}

				if(m_vOutG1mat[jointGind[k]][0]!=m_vOutG1mat[jointGind[k]][1]) { // heterozygous variant
					genoProd *=2; //multiply by 2  if het. variant
				}
			} //end for each contributor

			//Calculating the inner sum-part begins here:
			//Scaling shape parameter with degrad model:
			for (a = 0; a < m_nNumAlleles; a++) { //traverse each observed alleles (indices), having PH>0
				shapev[a] *= shapev0[a]; //scaling with degradation model (assumed already scaled)
			}
			
			//Summing up contribution of each alleles:
			bool anyDropin = false;
			double logevidProb = 0.0; //evidence probability (weight)
			double shape1; //modified shape param (caused by stutters)
			for (a = 0; a < (m_nNumAlleles-1); a++) { //traverse each observed alleles (indices), having PH>0 (no one are dropouts). LAST INDEX IS A Q-allele (PH=0)
				
				//Stutter model calculations (for given allele)
				shape1 = 0.0; //copy original contribution
				//CONSIDERING STUTTERS: Always running (with param equal zero)?? 
				if(m_vForwardStutters[a]!=-1) { //if a-1 stutter of a exists (found in FWstutter info)
					shape1 = (*xiF0)*shapev[m_vForwardStutters[a]]; //proportion obtained from allele a-1
				} 
				if(m_vBackstutters[a]!=-1) { //if a+1 stutter of a exists (found in BWstutter info)
					shape1 += (*xiB0)*shapev[m_vBackstutters[a]]; //proportion obtained from allele a-1
				}
				shape1 += (1- *xiB0 - *xiF0)*shapev[a]; //updated shape parameter (summing fraction lost and attained), proportion xiB/xiF always lost (possible to unseen potential stutters see later)
	
				if( shape1>m_dSmalltol ) { //If contribution 
					logevidProb += -lgamma(shape1) - shape1*const2 + (shape1-1)*log(m_vIntensity[a]) - constY0[a]; //gamma/lgamma is found in cmath							
				} else { //no contribution and PH present
					anyDropin = true; //there was at least one dropin
					logevidProb += dropin0[a]; //add to dropin prob to sum			
				} 
			} //end for each observed alleles
						
			//Calculate for the dropout allele (also including potential stutters:
			if(shapev[m_nNumAlleles-1]>m_dSmalltol) { //only if any contributors
					logevidProb += pgamma(*AT, shapev[m_nNumAlleles-1],scale, 1, 1); //log( gamma_p(shapev[m_nNumAlleles-1],(*AT)*const1) ); 
			}
			//Need to add this block to calc unobserved potential stutters
			for (a = 0; a < m_nNumPotentUnobsStutters; a++) { //traverse for each unobserved potential stutters
				shape1 = 0.0; //set to zero by default
				if(m_vPotentForwardStutters[a] != -1) { //if a-1 stutter of a exists (found in FWPstutter info)
					shape1 += (*xiF0)*shapev[m_vPotentForwardStutters[a]]; //proportion obtained from allele a-1
				} 
				if(m_vPotentBackStutters[a] != -1) { //if a+1 stutter of a exists (found in BWPstutter info)
					shape1 += (*xiB0)*shapev[m_vPotentBackStutters[a]]; //proportion obtained from allele a+1
				}
				if( shape1 > m_dSmalltol ) { //if there was a contribution, we need to calculate the prob for dropout
						logevidProb += pgamma(*AT, shape1,scale, 1, 1); //log( gamma_p(shape1,(*AT)*const1) ); //
				}
			}
			
			if(!anyDropin) { //if no drop-in occured
				logevidProb += log(1-*prC); //
			}

			//bigsum += evidProb*genoProd; //add temp res to big sum
			bigsum += exp(logevidProb)*genoProd; //add temp res to big sum
		} //end for each combinations (bigsum)	
		return log(bigsum); //return function
	}


	//Calculate cumulative probability for all observations in marker 
	void calccumval(double *maxY, double *fst,double *AT, double *lambda, double *prC, double *mixprop, double *mu, double *sigma, double *beta, double *xiB0, double *xiF0, double *ibd) {
		//maxY is max value of PH outcome
		int	nGK = pow(m_nNumGenos1p, m_nNOU); //get total number of combined genotypes (for unknown contrs)
		double invshape = (*sigma)*(*sigma);//get inverse of shape parameter
		double shape = 1/invshape; //get shape parameter
		double scale = (*mu)*invshape; //get scale parameter
		const double const1 = 1/(scale); //constant 1
	 	const double const2 = log(scale); //constant 2
			
		//Precalculation for known contributors:
		vector<double> shapevK(m_nNumAlleles, 0.0); //contribution vector of each alleles after taking into account mx

		//Sum up contribution for each alleles (Taking into account mix proportions): ONLY CALCULATED FOR KNOWN CONTRIBUTORS INITIALLY
		int aa, kk;
		for (kk = 0; kk < m_nNOK; kk++) { //for each known contributors 
			for (aa = 0; aa < m_nNumAlleles; aa++) { //traverse each alleles (indices)
				shapevK[aa] += m_vOutG1contr[m_vKnownGenoIdx[kk]][aa] * mixprop[m_vKnownContribIdx[kk]]; //contr from contr k to allele a
			}
		}

		//For each observations we traverse the modified "large sum" 
		cumprobvals.assign(m_nNumRep*m_nNumAlleles, vector<double>(2, -1.0)); //init (m_nNumRep*m_nNumAlleles) x 2 matrix (col1=yA,col2=T)
		
		//double i; //indicating which allele (for a replicate) that is used 
		for(i=0; i<(m_nNumRep*m_nNumAlleles); i++) { //Traverse all PHs (all reps for each alleles)
			if(m_vIntensity[i] < *AT) {
				continue; //skip if not observed
			}
			for(j=0; j<2; j++) { //traverse 0-2 (0=peak height considered, 1=AT considered,2=MaxY)
			
				//CALCULATING LARGE SUM:
				double bigsum = 0.0;  //total sum over all genotypes
				#pragma omp parallel for reduction(+:bigsum) //perform parallel on summing up bigsum
				for (int largeIter = 0; largeIter < nGK; largeIter++) { //for each combined/joint genotype outcome
					//Following variables must be declared for each iteration (non-shared):
					int a,k,r; //used to traverse alleles(a), contributors(k) and replicates(r)
					int aind; //used as index for alleles in genotypes
					vector<int> jointGind(m_nNOU, 0); //This will be the permuation contribution index for the unknown inds (directly corresponds to indices of Gmarg)
					vector<double> shapev = shapevK; //make a copy of existing shapevector
					vector<double> m_vAlleleCounts2 = m_vAlleleCounts; //Creating copy of counter (per allele)
					double m_nNumTyped2 = m_nNumTyped; //creating copy of total counter

					//jointGind = digits(largeIter, base = m_nNumGenos1p, pad = m_nNOU); #Equivalent operaion
					double genoProd = 1.0; //calculating the genotype probability of the unknowns
					int modrest = largeIter; //used to keep remained after modulo (init as iter number)
					bool inserted = false; //boolean of whether all digits are inserted (Rest is zero padded)
					for (k = 0; k < m_nNOU; k++) { //for each unknown contributors (summing up wrt both contr (m_vOutG1contr) and mx (mixprop)): Need each contr to derive shapev
						if (!inserted) { //if not all digits inserted
							if ( k>0 ) {
								modrest = int((modrest - jointGind[k - 1]) / m_nNumGenos1p); //extract remaining, divide to get to next digit (necessary with int converion?)
							}
							jointGind[k] = modrest % m_nNumGenos1p; //INSERT NUMBER: convert number to "m_nNumGenos1p" basis 	
							if (modrest < m_nNumGenos1p) { //check if rest is smaller than base (only run if not inserted)
								inserted = true;  //then all digits are inserted
							}	
						} //else { jointGind[k] = 0; //zero pad	}
		
						//Sum up contribution for each alleles (Taking into account mix proportions)
						for (a = 0; a < m_nNumAlleles; a++) { //traverse each alleles (indices)
							shapev[a] += m_vOutG1contr[jointGind[k]][a] * mixprop[m_vUnknownContribIdx[k]]; //contr from contr k to allele a. NOTICE THE k+*m_nNOK shift!
						}
						
						//CALCULATE GENOTYPE PROBS OF UNKNOWNS (MAY BE RELATED != -1)
						if( m_nUnknownContribRelIdx != k ) {
							for(a = 0;a<2; a++) {
								aind = m_vOutG1mat[jointGind[k]][a]; //get allele index of genotype g_a
								genoProd *= ((*fst)*m_vAlleleCounts2[aind] + (1-*fst)*m_vFreq[aind]) / (1 + (m_nNumTyped2-1)*(*fst)); 
								m_vAlleleCounts2[aind] += 1; //update allele count for particular genotype
								m_nNumTyped2 += 1; //update total count
							}
							if(m_vOutG1mat[jointGind[k]][0]!=m_vOutG1mat[jointGind[k]][1]) { // heterozygous variant
								genoProd *=2; //multiply by 2  if het. variant
							} 
						} 
						
						//AN UNKNOWN RELATED IS CALCULATED LAST (AFTER ALLELE COUNT UPDATES)
						if( m_nUnknownContribRelIdx != -1 ) { //if an unknown was related (k=m_nUnknownContribRelIdx)
							genoProd *= prob_relUnknown( jointGind[m_nUnknownContribRelIdx], m_nRelGenoIdx, (ibd + 3*m_vUnknownContribIdx[m_nUnknownContribRelIdx]), fst, &(m_vAlleleCounts2[0]), &m_nNumTyped2);	 //Note: scale with 3 because ibd is a '3-long vector' per contributor						
						}						
					} //end for each unknown contributor

					//////////////////////////////////////////////
					//Calculating the inner sum-part begins here//
					//////////////////////////////////////////////


					//Scaling shape parameter with degrad model:
					for (a = 0; a < m_nNumAlleles; a++) { //traverse each observed alleles (indices), having PH>0
						shapev[a] *= exp( log(shape) + m_vBasePairs[a]*log(*beta) ); //scaling with degradation model (assumed already scaled)
					}
					//Scaling shape parameter with stutter model:
					vector<double> shapev2 = shapev;  //make a copy of existing shapevector
					double shape1; //modified shape param (caused by stutters)
					
					//Stutter model calculations (for given allele)
					if( (*xiB0>m_dSmalltol) || (*xiF0>m_dSmalltol) ) { //only run if stutterparam is positive
						for (a = 0; a < (m_nNumAlleles-1); a++) { //traverse each alleles. Potential non-Qallele is traversed below
							shape1 = 0.0; //copy original contribution
							
							//CONSIDERING STUTTERS: Always running (with param equal zero)
							if(m_vForwardStutters[a] != -1) { //if a-1 stutter of a exists (found in FWstutter info)
								shape1 = (*xiF0)*shapev[m_vForwardStutters[a]]; //proportion obtained from allele a-1
							} 
							if(m_vBackstutters[a] != -1) { //if a+1 stutter of a exists (found in BWstutter info)
								shape1 += (*xiB0)*shapev[m_vBackstutters[a]]; //proportion obtained from allele a-1
							}
							shapev2[a] = (1- *xiB0 - *xiF0)*shapev[a] + shape1; //updated shape parameter (summing fraction lost and attained), proportion xiB/xiF always lost (possible to unseen potential stutters see later)
						} 
						
						//HANDLING SPECIAL SCENARIO WHEN LAST ALLELE IS NOT Q-allele
						a = (m_nNumAlleles-1); //set allele index as last allele
						if(m_vForwardStutters[a] != -1 || m_vBackstutters[a] != -1) {  //IF THERE WAS STUTTER INDEX DEFINED
							shape1 = 0.0; //copy original contribution									
							if(m_vForwardStutters[a] != -1) { shape1 = (*xiF0)*shapev[m_vForwardStutters[a]]; } //proportion obtained from allele a-1
							if(m_vBackstutters[a] != -1) { shape1 += (*xiB0)*shapev[m_vBackstutters[a]]; } //proportion obtained from allele a+1
							shapev2[a] = (1- *xiB0 - *xiF0)*shapev[a] + shape1; //updated shape parameter (summing fraction lost and attained), proportion xiB/xiF always lost (possible to unseen potential stutters see later)				
						} //end special scenario
					}
					
					//Summing up contribution of each alleles:
					vector<bool> anyDropin(m_nNumRep,false); //use vector to indicate drop-in
					double logevidProb = 0.0; //evidence probability (weight)
					int cind; //cumulative allele indexing for PHs (for traversing over all replicates)
					for (a = 0; a < m_nNumAlleles; a++) { //traverse each observed alleles (also consider the last allele since Q-allele may have PH last allee)
						for(r = 0; r < m_nNumRep; r++) { //traverse each replicates (observed alleles indicated by PH)
							cind = a*m_nNumRep + r; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
						
							if( cind!=i ) { //proceed as before if different index
								if( m_vIntensity[cind] > m_dSmalltol ){ //if evaluating index was not ii and PH>0 
									if( shapev2[a] > m_dSmalltol) { //If contribution and PH>0  					//CONTRIBUTION SET (A)
										logevidProb += -lgamma(shapev2[a]) - shapev2[a]*const2 + (shapev2[a]-1)*log(m_vIntensity[cind]) - m_vIntensity[cind]*const1; //gamma/lgamma is found in cmath							
									} else { //If contribution and PH>0 	//DROPIN SET (C)
										anyDropin[r] = true; //there was at least one dropin
										logevidProb += log(*lambda) - (*lambda)*(m_vIntensity[cind] - (*AT) ) + log(*prC) + log(m_vFreq[a]); //add to dropin prob to sum										
									}							
								}  else if(shapev2[a]>m_dSmalltol) { //IF PH missing  and contribution(it't a dropout //CONTRIBUTION SET (B)
									logevidProb += pgamma(*AT, shapev2[a],scale, 1, 1); //log( gamma_p(shapev2[a],(*AT)*const1) );  //Add log(dropout-probability)
								}
							} else { //modify if allele index is same (assures that m_vIntensity[cind]>=(*AT) (see early in loop)
								if(shapev2[a]>m_dSmalltol) { //If contribution and PH>0  					//CONTRIBUTION SET (A)
									if(j==0) { //CASE OF PH
										logevidProb += log( pgamma(m_vIntensity[cind], shapev2[a],scale, 1, 0) - pgamma(*AT-1, shapev2[a],scale, 1, 0));
										//logevidProb += log( gamma_p(shapev2[a],m_vIntensity[cind]*const1) -  gamma_p(shapev2[a],(*AT)*const1) );  //evaluate for upper limit - lower limit
									} else if(j==1) { //in case of j=1:  //CASE OF maxY
										logevidProb += log( pgamma(*maxY, shapev2[a],scale, 1, 0) - pgamma(*AT-1, shapev2[a],scale, 1, 0));
										//logevidProb += log( gamma_p(shapev2[a],(*maxY)*const1) - gamma_p(shapev2[a],(*AT)*const1) );  //evaluate for upper limit - lower limit										
									} 
									
								} else { //If contribution and PH>0 //DROPIN SET (C)
									anyDropin[r] = true; //there was at least one dropin
									logevidProb += log(*prC) + log(m_vFreq[a]); //add to cumulative dropin prob to sum	
									if(j==0) { //CASE OF PH		
										logevidProb += log( 1 - exp(- (*lambda)*( m_vIntensity[cind] - (*AT-1) )) ); //calc logarithm of cumulative expression
									} else if(j==1) { //CASE OF maxY
										logevidProb += log( 1 - exp(- (*lambda)*(*maxY-(*AT-1)))); //calc logarithm of cumulative expression
									} 							
								}							
							}							
						} //end for each replicates (r)
					} //end for each observed alleles
						
					//Need to add this block to calc unobserved potential stutters (done ones for all replicates, since this can be scaled)
					if( (*xiB0>m_dSmalltol) || (*xiF0>m_dSmalltol) ) { //only run if param is positive
						for (a = 0; a < m_nNumPotentUnobsStutters; a++) { //traverse for each unobserved potential stutters
							shape1 = 0.0; //set to zero by default
							if(m_vPotentForwardStutters[a] != -1) { //if a-1 stutter of a exists (found in FWPstutter info)
								shape1 = (*xiF0)*shapev[m_vPotentForwardStutters[a]]; //proportion obtained from allele a-1
							} 
							if(m_vPotentBackStutters[a] != -1) { //if a+1 stutter of a exists (found in BWPstutter info)
								shape1 += (*xiB0)*shapev[m_vPotentBackStutters[a]]; //proportion obtained from allele a+1
							}
							if( shape1 > m_dSmalltol ) { //if there was a contribution, we need to calculate the prob for dropout
								logevidProb += m_nNumRep*pgamma(*AT, shape1,scale, 1, 1); //log( gamma_p(shape1,(*AT)*const1) ); //scaling with number replicates
							}
						}
					}
					
					//Liklihood NO DROP-IN: Traverse each replicates (observed alleles indicated by PH)
					for(r = 0; r<m_nNumRep; r++) {
						if(!anyDropin[r]) { //if no drop-in occured for a replicate
							logevidProb += log(1-*prC); //
						}
					}
					
					//bigsum += evidProb*genoProd; //add temp res to big sum
					bigsum += exp(logevidProb)*genoProd; //add temp res to big sum
				} //end for each combinations (bigsum)
				cumprobvals[i][j] = bigsum; //storing non-logged probability value 
			} //end for each iteration j = 0-2
		} //end for each PH iteration i
	} //end function
	
};


/*MAIN IMPORTS DATA AN STRUCTURES THE DATA FOR MARKER BASED FUNCTIONS*/
extern "C" {
		

//The P(E|theta)= sum_g P(E|g,theta)*P(g) expression
void loglikgammaC(double *logLik, int *NOC, int *NOK, int *knownGind, double *mixprop, double *mu, double *sigma, double *beta, double *xiB,double *xiF, double *AT, double *pC, double *lambda, double *fst,  int *m_nNumRep, int *nM, int *nA, double *peaksLong, double *freqsLong, double *m_nNumTypedLong, double *maTypedLong, double *basepairLong, int *BWstuttindLong, int *FWstuttindLong, int *nPS, int *BWPstuttindLong, int *FWPstuttindLong, int *maxThreads, int *isPhi, int *anyRel, int *relGind, double *ibd) {
	//return items:
	//logLik = the joint log likelihood over all markers
	
	//MOdel params:
	//NOC = number of known contributors (constant)
	//NOK = number of unknown contributors (per marker)	
	//knownGind = genotype index of the known contributors
	//mixprop = mixture proportions (NOC-1 long vector)
	//mu = expectation of peak heights (numMarkers long vector)
	//sigma = coefficient-of variation of peak heights (numMarkers long vector)
	//beta = degradation slope (numMarkers long vector)
	//xiB = Backward stutter proportion (expectation)
	//xiF = Forward stutter proportion (expectation)
	//ibd = kappa-coefficients to specify the relatedness between unknowns and related NOC*3 vector long
	
	//Pre-specified params:
	//AT = Analytical threshold (LOD)  (numMarkers long vector)
	//pC = Dropin prob. parameter  (numMarkers long vector)
	//lambda = Dropin hyperparameter for exponential distr (numMarkers long vector)
	//fst = Sub-population structure parameter (numMarkers long vector)
	
	//Data variables:
	//m_nNumRep = number of replicates (a vector per marker, different markers may have different number of replicates)
	//nM = number of markers (per sample): M_s_m M_1,1
	//nA = number of alleles per markers (population alleles)
	//peaksLong = vectorized vector of intensities (PeakHeights/Coverage), vectorized over all replicates (missing are indicated as zero)
	//freqsLong = vectorized vector of frequencies
	//m_nNumTypedLong = number of total number of previous typed alleles (used for fst). Vectorized
	//maTypedLong = number of total number of previous typed alleles (used for fst). Vectorized
	//basepairLong = vectorized vector of the already adjusted base pairs (x-150)/100 
	//BWstuttindLong = vectorized vector of backward stutter indices (gives index of what allele it receive stutter from  (Allele 1=index0))
	//FWstuttindLong = vectorized vector of forward stutter indices (gives index of what allele it receive stutter from  (Allele 1=index0))
	//nPS = number of potential stutters per marker
	//startIndMarkersPS = start index per marker (potential stutters)
	//BWPstuttindLong = vectorized vector of backward stutter indices (gives index of what allele it receive stutter from  (Allele 1=index0))
	//FWPstuttindLong = vectorized vector of forward stutter indices (gives index of what allele it receive stutter from  (Allele 1=index0))
	//relGind = vectorized indices of genotypes of the related individuals (nM*NOC long vector)

	//Other vars:
	//maxThreads = max number of threads used for paralellisation
	//isPhi = boolean of whether Real domain of mixture proportion variable should be considered (must transform back)
	//anyRel = boolean of whether any of the unknowns are related 
	#ifdef _OPENMP
	int numThreads = thread::hardware_concurrency();
	int useThreads = min(numThreads,*maxThreads);
	omp_set_num_threads(useThreads);  //set number of threads to use 
	#endif
	
	//Baseline of Backward (BW) expected stutter param per marker
	//Slope of Backward (BW) expected stutter param per marker
	
	//TRANSFORMING MIXTURE PROPORTIONS:
	vector<double> mixprop2(*NOC,1.0); //Vector with mixture proportions
	if(*NOC>1)  {
		double cs = 0.0; //init cumulative sum
		for(int i=0; i< (*NOC-1); i++) { //for each mix props
			if(*isPhi==1) { //if transformed variable
				mixprop2[i] = (1-cs)/(1+exp(-mixprop[i])); //transform further (see formula for explanation)
			} else {
				mixprop2[i] = mixprop[i]; //no transform
			}
			cs += mixprop2[i]; //add mixture proportion to cumulative sum
		} //end for each 
		mixprop2[*NOC-1] = 1.0 - cs; //Insert last contributor
    }

	//Calculating the loglik (over all markers), which is returned
	double jointloglik = 0.0;
	int startIndPS = 0; //start marker index for potential stutters (own vectors)
	int startIndMarker1=0;//start marker index (1 rep)
	int startIndMarker2=0;//start marker index (m_nNumRep[m] reps)
	const double m_dSmalltol = 1.0e-30; //a tiny number > 0 (avoiding zero roundoff errors)
	for(int m=0; m< *nM; m++) {		
		EFMmarker *marker = new EFMmarker( NOC,(NOK+m),(knownGind + (*NOC)*m), (m_nNumRep+m), (nA+m), (peaksLong+startIndMarker2),(freqsLong+startIndMarker1), (m_nNumTypedLong+m), (maTypedLong+startIndMarker1),(basepairLong+startIndMarker1),(BWstuttindLong+startIndMarker1),(FWstuttindLong+startIndMarker1), (nPS+m), (BWPstuttindLong + startIndPS),(FWPstuttindLong + startIndPS), (relGind+(*NOC)*m) ); //prepare new marker

		//if only unknown unrelated with 1 replicate observed, last allele is dropout: RUN FASTEST CODE VARIANTS (LESS TO RUN)
		if(NOK[m]==0 && m_nNumRep[m]==1 && *anyRel==0 && peaksLong[startIndMarker1 + nA[m]-1 ] < m_dSmalltol ) {  
			if( xiB[m]>m_dSmalltol || xiF[m]>m_dSmalltol) { //IF stutter model version (assumes BW/FW stutter prop is positive
				jointloglik += marker->calcLogLikGammaMarkerUnknownStutter((fst+m), (AT+m), (lambda+m) ,(pC+m), &(mixprop2[0]), (mu + m),(sigma + m),(beta+m),(xiB+m),(xiF+m)) ; //insert marker information to vector			
			} else { //run model without stutter			
				jointloglik += marker->calcLogLikGammaMarkerUnknown((fst+m), (AT+m), (lambda+m) ,(pC+m), &(mixprop2[0]), (mu + m),(sigma + m),(beta+m)) ; //insert marker information to vector			
			}
		} else { //OTHERWISE A SLOWER CODE IS RUN
			jointloglik += marker->calcLogLikGammaMarker((fst+m), (AT+m), (lambda+m) ,(pC+m), &(mixprop2[0]), (mu + m),(sigma + m),(beta+m),(xiB+m),(xiF+m),ibd) ; //insert marker information to vector			
		}
		if( isinf(jointloglik) ) {
			break; //stop calc if -INF
		}
		//*(logLikv+m) = jointloglik; //store likelihood value for marker
		startIndMarker1 += nA[m]; //get start position of marker m+1 (1 rep)
		startIndMarker2 += nA[m]*m_nNumRep[m]; //get start position of marker m+1 (m_nNumRep), used only for Peaks only
		startIndPS += nPS[m]; //add number of potential stutters
		delete marker; //delete after each declaration
	}
	*logLik = jointloglik;
} //end main function



//Obtaining cumulative probabilities for all observations (used for PP plots)
void cumvalgammaC(double *pvalPH, double *pvalMAX, double *maxY, int *NOC, int *NOK, int *knownGind, double *mixprop, double *mu, double *sigma, double *beta, double *xiB,double *xiF, double *AT, double *pC, double *lambda, double *fst,  int *m_nNumRep, int *nM, int *nA, double *peaksLong, double *freqsLong, double *m_nNumTypedLong, double *maTypedLong, double *basepairLong, int *BWstuttindLong, int *FWstuttindLong, int *nPS, int *BWPstuttindLong, int *FWPstuttindLong, int *maxThreads, int *relGind, double *ibd) {
	//pval PH/AT/MAX  cumulative value vectors with size m_nNumReps*m_nNumAlleles' 
	//maxY	is maximum value to use

	#ifdef _OPENMP
	int numThreads = thread::hardware_concurrency();
	int useThreads = min(numThreads,*maxThreads);
	omp_set_num_threads(useThreads);  //set number of threads to use 
	#endif
	
	//Baseline of Backward (BW) expected stutter param per marker
	//Slope of Backward (BW) expected stutter param per marker
	
	//TRANSFORMING MIXTURE PROPORTIONS:
	vector<double> mixprop2(*NOC,1.0); //Vector with mixture proportions
	if(*NOC>1)  {
		double cs = 0.0; //init cumulative sum
		for(int i=0; i< (*NOC-1); i++) { //for each mix props
			mixprop2[i] = mixprop[i]; //no transform
			cs += mixprop2[i]; //add mixture proportion to cumulative sum
		} //end for each 
		mixprop2[*NOC-1] = 1.0 - cs; //Insert last contributor
    }

	//Calculating the pvalues for each observations (provided peaksLong)
	int startIndPS = 0; //start marker index for potential stutters (own vectors)
	int startIndMarker1=0;//start marker index (1 rep)
	int startIndMarker2=0;//start marker index (m_nNumRep[m] reps)
	int i,j; //used to traverse alleles and replicates
	for(int m=0; m< *nM; m++) {		
		EFMmarker *marker = new EFMmarker( NOC,(NOK+m),(knownGind + (*NOC)*m), (m_nNumRep+m), (nA+m), (peaksLong+startIndMarker2),(freqsLong+startIndMarker1), (m_nNumTypedLong+m), (maTypedLong+startIndMarker1),(basepairLong+startIndMarker1),(BWstuttindLong+startIndMarker1),(FWstuttindLong+startIndMarker1), (nPS+m), (BWPstuttindLong + startIndPS),(FWPstuttindLong + startIndPS), (relGind+(*NOC)*m) ); //prepare new marker
		marker->calccumval(maxY, (fst+m), (AT+m), (lambda+m) ,(pC+m), &(mixprop2[0]), (mu + m),(sigma + m),(beta+m),(xiB+m),(xiF+m),ibd) ; //calculate: ONLY ONE CODE VARIANT FOR THIS
		for(i=0; i < nA[m]*m_nNumRep[m]; i++) { //For each alleles in marker
			j = startIndMarker2 + i; //get index over set
			*(pvalPH+j) = marker->cumprobvals[i][0]; //copy value for cumulative on PH
			*(pvalMAX+j) = marker->cumprobvals[i][1]; //copy value for cumulative on MAX
		}		
		startIndMarker1 += nA[m]; //get start position of marker m+1 (1 rep)
		startIndMarker2 += nA[m]*m_nNumRep[m]; //get start position of marker m+1 (m_nNumRep), used only for Peaks only
		startIndPS += nPS[m]; //add number of potential stutters
		delete marker; //delete after each declaration
	}
} //end main function

} //end external
