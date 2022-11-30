//This script was originally written to calculate the likelihood value for a replicate model (EFMrep).
//AUTHOR: Oyvind Bleka, July 2021, modified September 2022
/*ABOUT:
- Contains ordinary implementation of EFM:
	Used for Deconvolution. Genotypes are traversed only once.	
- No structuring of data is needed (already done in prepC)
- All markers are looped outside the inner large-sum loop.
*/


#include <vector> //vector storage
#include <cmath> //includes lgamma
#include <thread> //used to obtain number of logical processes
#include <Rmath.h> //includes pgamma
#include "genoProbFuns.h" //helpfunction for genotype probs
#ifdef _OPENMP
#include <omp.h> //parallelization
#endif

using namespace std;

/*MAIN IMPORTS DATA AN STRUCTURES THE DATA FOR MARKER BASED FUNCTIONS*/
extern "C" {
		

//The P(E|theta)= P(E|g,theta)*P(g) expression is returned
void loglikGamma_allcomb2( double *loglikVEC, int *startIndMarker_nJointGenos, int *nJointGenos, int *NOC, int *NOK, int *nKnowns, 
	double *mixProp, double *PHexp, double *PHvar, double *DEG, double *stutt,  //stutt is stutter proportion vector (1 element for each kind)
	double *AT, double *fst, double *dropinProb, double *dropinWeight, 
    int *nMarkers, int *nRepMarkers, int *nAllelesVEC, int *NumPotStutters, int *startIndMarker_nAlleles, int *startIndMarker_nAllelesReps,	
	double *peaksLong, double *freqsLong, double *nTypedLong, double *maTypedLong, double *basepairLong, int *nStutters, int *stuttFromIndVEC, 
	int *stuttToIndVEC, int *stuttParamIndVEC, int *startIndMarker_nStutters, int *knownGind, int *relGind, double *ibd, int *maxThreads) {	
	
	#ifdef _OPENMP
		int numThreads = thread::hardware_concurrency();
		if(*maxThreads>0) numThreads = min(numThreads,*maxThreads); 
		omp_set_num_threads(numThreads); //preparing CPU parallelization	
	#endif
	
	int aa, kk, rr; //indices for alleles, contributors and replicates
	int nReps = 1; //fix to 1
	//Prepare transformation of parameters (per replicate) to save computation:
	vector<double> shape0(nReps,0.0);
	vector<double> scale0(nReps,0.0);
	vector<double> const1(nReps,0.0);
	vector<double> const2(nReps,0.0);
	for(rr=0; rr< nReps; rr++) { //Traverse each replicate vars (param potential different for each)
		double PHvarSq = PHvar[rr]*PHvar[rr]; //Square param
		shape0[rr] = 1/PHvarSq; // //theta_Am[locind]/theta_omegasq; //shape for 'full het allele'    
		scale0[rr] = PHexp[rr]*PHvarSq; //obtain scale param 
		const1[rr] = 1/scale0[rr]; //constant 1 (used in pgamma,dgamma)
		const2[rr] = log(scale0[rr]); //constant 2 (used in dgamma)		
	}
	
	//Calculating over all markers
	const double smalltol = 1.0e-30; //a tiny number > 0 (avoiding zero roundoff errors)
	for(int locind=0; locind< *nMarkers; locind++) {	//for each marker:
	
		//OBTAIN CONSTANTS FOR SPECIFIC MARKER:
	
		//default settings
		double fst0 = fst[locind]; //theta-correction param	(marker specific)
		double dropinProb0 = dropinProb[locind]; //obtain dropin prob (marker specific, but same for all replicates)
		double AT0 = AT[locind]; //obtain AT (marker specific, but same for all replicates)
			
		//Prepare dimensions and data vectors
		int NOK0 = nKnowns[locind]; //obtain number of contributors (may be different for different markers)
		int NOU = *NOC - NOK0; //number of unknowns (may be different for markers)
		int nReps0 = nRepMarkers[locind]; //number of replicates for specific marker (may be different for different markers)
		
		//Prepare dimen
		int nAlleles = nAllelesVEC[locind]; //number of alleles
		int nPS = NumPotStutters[locind]; //number of potential stutters
		int nAlleles2 = nPS + nAlleles; //number of alleles (including potential stutters)

		//Create contribution matrix (1 genotype):
		int numGenos1p = int(nAlleles*(nAlleles + 1) / 2); //get Genotype outcome
		vector<double> outG1contr(nAlleles*numGenos1p, 0); //init nG1xnA matrix (contribution matrix). Indicating what alleles that are contributoed
		vector<int> outG1allele(2*numGenos1p,0);  //init nG1x2 matrix (allele names as indices 0,...,nA-1)
		int cc = 0; //counter oveer all genotypes
		for(int i = 0; i < nAlleles; i++) {
			for(int j = i; j < nAlleles; j++) { 
				outG1allele[2*cc + 0] = i; //include index
				outG1allele[2*cc + 1] = j; //include index
				outG1contr[nAlleles*cc + i] += 1.0; //insert contr at allele i
				outG1contr[nAlleles*cc + j] += 1.0; //insert contr at allele j
				cc++; //iterate to next genotype outcome
			}
		}
		
		//Obtain start indices (because of traversing all markers)
		int SI_nAlleles0 = startIndMarker_nAlleles[locind]; //index start number of alleles
		int SI_nAllelesReps0 = startIndMarker_nAllelesReps[locind]; //index start number of alleles (taking into account Reps)
		int SI_nJointGenos0 = startIndMarker_nJointGenos[locind];
		//int SI_nRepMarkers0 = startIndMarker_nRepMarkers[locind]; //index start number of alleles (taking into account Reps)
		//int SI_outG1allele0 = startIndMarker_outG1allele[locind];
		//int SI_outG1contr0 = startIndMarker_outG1contr[locind];
		int SI_nStutters0 = startIndMarker_nStutters[locind];
		
		//Prepare vector for known contributors (and also unknown): Need to know positions!
		vector<int> GindKnown(NOK0,0); //genotype index in vector 
		vector<int> kindKnown(NOK0,0); //contributor index in vector
		vector<int> kindUnknown(NOU,0); //contributor index in vector		
		int kindRel = relGind[locind]; //Put last in unknown traversion (-1 means no related)
	
		cc = 0; //counters for known
		int jj = 0; //counter for unknowns
		int knownGind0; //obtain genotype index of known 
		for(kk=0; kk< *NOC; kk++) { //for each contributors:
			//kindRel[kk] = relGind[ (*NOC)*locind + kk ]; //copy genotype index
			knownGind0 = -1; //init
			if(kk<*NOK) knownGind0 = knownGind[ (*NOK)*locind + kk ]; //copy genotype index  (KNOWN)			
			if(knownGind0>=0) { //If contributor is known (genotype given)
				GindKnown[cc] = knownGind0; //copy genotype index
				kindKnown[cc] = kk; //insert contributor index
				cc++; //update counter for knowns
			} else { //if contributor is unknown (genotype not given)
				kindUnknown[jj] = kk; //insert contributor index
				jj++; //update counter for unknowns
			}
		}
				
		//Prepare variables (before sum-iterations): Save computations
		double nTyped = nTypedLong[locind]; //copy total number of typed alleles
		vector<double> maTypedvec(nAlleles,0.0); //copy number of typed alleles (each type)
		vector<double> shapevK(nAlleles2, 0.0); //Precalculation for known contributors (vectorized over all replicates of particular marker)
		vector<double> shapev0(nAlleles2, 0.0); //degrad scaling of shape (init vector) (vectorized over all replicates of particular marker)	
		//int aaind; //vectorization index (see info below)
		int repIDmarker=0; //Used to indicate correct replicate (some markers may have fewer replicates!!!)
		//int mixPropContrInd; //used to indicate correct mixture proportion vector (of particular replicate)
		for (aa = 0; aa < nAlleles; aa++) { //traverse each observed alleles (indices), also the Q-allele. Potential not necessary!
			maTypedvec[aa] = maTypedLong[ SI_nAlleles0 + aa ]; //copy previously typed alleles
			//aaind = aa; //get index vectorized vector which includes (allele,rep) indices as [(1,1), (1,2), (1,3), (2,1), (2,2) etc]
			shapev0[aa] = exp( log( shape0[repIDmarker] ) + basepairLong[SI_nAlleles0 + aa ]*log( DEG[repIDmarker] ) ); //scaling with degradation model (assumed already scaled)
			
			//Sum up contribution for each alleles (Taking into account mix proportions): ONLY CALCULATED FOR KNOWN CONTRIBUTORS INITIALLY
			for (kk = 0; kk < NOK0; kk++) { //for each known contributors 	
				//mixPropContrInd = kindKnown[kk] + repIDmarker*(*NOC); //obtain correct index of mixture proportion (vectorized across replicates)
				shapevK[aa] += outG1contr[ nAlleles*GindKnown[kk] + aa] * mixProp[ kindKnown[kk] ]; //contr from contr k to allele aa
			} //end for each known contributor
		}
		
		#pragma omp parallel for //shared(shapev0, dropin0,constY0) //perform parallel on summing up bigsum
		for (int iter = 0; iter < nJointGenos[locind]; iter++) { //for each combined/joint genotype outcome
		
			//int iter = combUseVEC[ SI_nJointGenos0 + iter2]; //obtain correct iteration index to consider

			//Following variables PHexpst be declared for each iteration (non-shared):
			int a,k,r; //used to traverse alleles(a), contributors(k) and replicates(r)
			int aind; //used as index for alleles in genotypes
			int repID=0; //Used to indicate correct replicate (some markers may have fewer replicates!!!)
			//int mixPropInd; //used to indicate correct mixture proportion vector (of particular replicate)
			//vector<int> jointGind(NOU, 0); //This will be the contribution index for the unknown inds (directly corresponds to indices of Gmarg)
			vector<double> shapev = shapevK; //make a copy of existing shapevector
			vector<double> maTypedvec2 = maTypedvec; //Creating copy of counter (per allele)
			double nTyped2 = nTyped; //creating copy of total counter

			//jointGind = digits(iter, base = numGenos1p, pad = NOU); #Equivalent operaion
			//int startInd_alleles; //init for allele index
			double genoProd = 1.0; //calculating the genotype probability of the unknowns
			int genotypePower = 1; 
			int genoidx; //genotype index		
			int allele_contr; //allele contribution (decided by genotype combinations)
			for (k = 0; k < NOU; k++) { //for each unknown contributors (summing up wrt both contr (outG1contr) and mx (mixProp)): Need each contr to derive shapev
				genoidx = (iter/genotypePower) % numGenos1p; // //get genotype for partifular contributor	
				genotypePower *= numGenos1p; //update genotype comb	 ready for next contr

				//Sum up contribution for each allele:rep (Taking into account mix proportions for each replicate)
				for (a = 0; a < nAlleles; a++) { //traverse each alleles (indices)				
					allele_contr = outG1contr[ nAlleles*genoidx + a]; //+ SI_outG1contr0 obtain allele contribution
					if( allele_contr==0 ) continue; //skip allele if no contribution from genotypes (SPEEDUP???)
					//aind = a; //get index vectorized vector which includes (allele,rep) indices as [(1,1), (1,2), (1,3), (2,1), (2,2) etc]
					//repID = repIDvec[ SI_nRepMarkers0 + r ]; //obtain correct repID for marker		
					//mixPropInd = kindUnknown[k] + repID*(*NOC); //obtain correct index of mixture proportion (vectorized across replicates)
					shapev[a] +=  allele_contr*mixProp[  kindUnknown[k] ]; //contr from contr k to allele a. NOTICE THE k+*NOK0 shift!
				}
		
				//CALCULATE GENOTYPE PROBS OF UNKNOWNS (MAY BE RELATED:  kindRel!= -1)
				//Sending pointer where marker indices are starting for outG1 and freqLong vector
				//SI_outG1allele0 + 2*genoidx; //obtain genotype index of unknown
				int aindU = outG1allele[2*genoidx];
				int bindU = outG1allele[2*genoidx+1];
					
				//Handle related individual (placed last)
				int Rgind = -1; //genotype of related unknown
				int aindR = 0;
				int bindR = 0;
				if(k==(NOU-1) && kindRel > -1) { //check if last visited
					Rgind = kindRel; //copy
					aindR = outG1allele[2*Rgind ]; 
					bindR = outG1allele[2*Rgind+1]; 
				}
				genoProd *= prob_relUnknown(aindU,bindU, genoidx, freqsLong + SI_nAlleles0, fst0,  &(maTypedvec2[0]), nTyped2,  aindR, bindR, Rgind, ibd );	 //Note: scale with 3 because ibdLong is a '3-long vector' per contributor						
						
				//LAST: UPDATE COUNTERS FOR ALLELES (a and b)
				maTypedvec2[ outG1allele[2*genoidx  ] ] +=1; //update allele count for particular allele (1)
				maTypedvec2[ outG1allele[2*genoidx+1] ] +=1; //update allele count for particular allele (2)					
				nTyped2 += 2; //update total count
			} //end for each unknown contributor

			//////////////////////////////////////////////
			//Calculating the inner sum-part begins here//
			//////////////////////////////////////////////

			//Scaling shape parameter with degrad model (could be done elsewhere?  done last in prev operatoin)
			for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (indices), having PH>0
				//aind = a;//*nReps0 + r; //get index vectorized vector which includes (allele,rep) indices as [(1,1), (1,2), (1,3), (2,1), (2,2) etc]
				shapev[a] *= shapev0[a]; //scaling with degradation model (assumed already scaled)					
			}
			
			//Scaling shape parameter with stutter model:			
			vector<double> shapev2 = shapev;  //make a copy of existing shapevector
			int stuttind,stuttind2,stuttFromInd,stuttToInd,stuttParamInd; //init vars
			//repID = repIDvec[ SI_nRepMarkers0 + r ]; //obtain correct repID for marker
			for(stuttind=0; stuttind < nStutters[locind]; stuttind++) { 
				stuttind2 = SI_nStutters0 + stuttind; //obtain index of stutter
				stuttFromInd = stuttFromIndVEC[stuttind2]; //*nReps0 + r;  //NB: CAREFUL WITH INDEX
				if( shapev[stuttFromInd]>smalltol) { //ONLY NECESSARY TO PROVIDE MODIFICATION IF shapeval>0
					stuttToInd = stuttToIndVEC[stuttind2]; //*nReps0 + r;  //NB: CAREFUL WITH INDEX
					stuttParamInd = stuttParamIndVEC[stuttind2]; // + 2*repID;  //Obtain parameter index (NOTE: ALWAYS using 2 elements per replicate). NB: CAREFUL WITH INDEX					
					shapev2[stuttToInd] += stutt[stuttParamInd] * shapev[stuttFromInd]; // #OBTAINED stutters
					shapev2[stuttFromInd] -= stutt[stuttParamInd] * shapev[stuttFromInd]; // #SUBTRACTED stutters
				}
			}

			//Summing up likelihood contribution of each alleles:reps:
			vector<int> nDropin(nReps0,0); //count number of drop-in for each replicate
			double logevidProb = 0.0; //evidence probability (weight)
			
			double peak; //obtaining observeed peak height 
			for(r = 0; r < nReps0; r++) { //traverse each replicates (observed alleles indicated by PH)
				//repID = repIDvec[ SI_nRepMarkers0 + r ]; //obtain correct repID for marker
				//AT0 = AT[SI_nRepMarkers0 + r]; //obtain AT to use (follows same as nReps per marker)
				for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (also consider the last allele since Q-allele may have PH last allee)				
					aind = a*nReps0 + r; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
					peak = peaksLong[SI_nAllelesReps0 + aind]; //Notice the cumulative locus shift (Long vector)
					
					//IF NOT DROPOUT
					if( peak >= AT0 ){ //if PH>0 (this is same as PH>=AT since threshold has already been applied)
						if(shapev2[a] > smalltol) { //If contribution and PH>0  					//CONTRIBUTION SET (A)
							logevidProb += -lgamma(shapev2[a]) - shapev2[a]*const2[repID] + (shapev2[a]-1)*log(peak) - peak*const1[repID]; //gamma/lgamma is found in cmath							
						} else { //If contribution and PH>0 	//DROPIN SET (C)							
							logevidProb += dropinWeight[ SI_nAllelesReps0 + aind ]; //likelihood for observed dropin (both PH and lambda are included)							
							nDropin[r] += 1; //count dropin for particular replicate															   
						}							
						
					//OTHERWISE IT IS DROPOUT
					} else if(shapev2[a] > smalltol) { //IF PH missing  and contribution(it't a dropout //CONTRIBUTION SET (B)
						logevidProb += pgamma(AT0, shapev2[a],scale0[repID], 1, 1); //Add log(dropout-probability)
						//log( gamma_p(shapev2[aind],AT0*const1[repID]) );  
					}
					
				} //end for each replicates (r)
			} //end for each observed alleles
			
			//Weight potential drop-outs (stutters)
			 if(nPS>0) {
				//repID = repIDvec[ SI_nRepMarkers0 + r ]; //obtain correct repID for marker
				//AT0 = AT[SI_nRepMarkers0 + r]; //obtain AT to use (follows same as nReps per marker)
				for(a=nAlleles; a<nAlleles2; a++) { //traverse remaining alleles
				  //aind = a;//*nReps0 + r; //get index vectorized vector which includes (allele,rep) indices as [(1,1), (1,2), (1,3), (2,1), (2,2) etc]
				  if(shapev2[a]>smalltol) { //IF DROPOUT						
				    //Note the scaling here
					logevidProb += nReps0*pgamma(AT0, shapev2[a],scale0[repID], 1, 1); //Add log(dropout-probability).
				  }
				}
			 }
						
			//Liklihood OF NUMBER OF NOISE/DROP-IN: Traverse each replicates (observed alleles indicated by PH)			
			for(r = 0; r<nReps0; r++) {
				//double dropinProb0 = dropinProb[SI_nRepMarkers0 + r]; //obtain dropin prob (follows same as nReps per marker)
				if(nDropin[r]==0) { //in case of no droping
					logevidProb += log(1-dropinProb0);
				} 
			}
			
			//FINAL INSERTION OF VALUES:
  		   loglikVEC[SI_nJointGenos0 + iter] = logevidProb + log(genoProd); //calculate P(E|gj)P(gj)			   		   
		} //end for each combination iterations (bigsum)
		//logLik[0] += log(bigsum); //calculate logLik by adding log(innerSUM)	
	} //end for each marker
} //end main function

} //end external