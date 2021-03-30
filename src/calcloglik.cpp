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


#include <iostream> //print to terminal
#include <fstream> ////READ FROM DATA
#include <vector> //vector storage
#include <cmath> //includes lgamma
#include <ctime> //calculate time
#include <thread> //used to obtain number of logical processes
#include <string> //
#include <stdio.h> //indludes printf-function
#include <omp.h> //parallelization
#include <boost/math/special_functions/gamma.hpp> 

using namespace std;
using namespace boost::math;
//using namespace dlib; Conflict with boost::math it seems...

class EFMmarker {
	private: 
	const double smalltol = 1.0e-30; //a tiny number > 0 (avoiding zero roundoff errors)
	int cc,i,j; //iterators
	int numGenos1p; //number of genotypes for 1 contributor 
	vector< vector<double> > outG1contr;
	vector< vector<int> > outG1mat;

	//DATA VARIABLES
	int nRep; //number of replicates (for marker)
	int nAlleles; //number of alleles in marker
	vector<double> Yvec; //Intensity vector (copied number)
	vector<double> Fvec; //Frequency vector
	vector<double> bpvec; //base pair information from the kitinfo (for each alleles). Used for degradation model
	vector<int> BWvec,FWvec; //Indices for Backward/forward stutter vector (indicating which alleles backward/forward stutter to whom)
	vector<double> maTypedvec; //counting number of typed alleles of type a
	double nTyped; //counting number of typed alleles (total)
	int nPstutters; //number of unobserved potential stutters (BW+FW)
	vector<int> BWPvec,FWPvec; //Indices for Backward/Forward potential unseen stutter vector (indicating which alleles backward/forward stutter to whom)

	//Known contributor variables	
	vector<int> GindKnown; //genotype index of known contributors (follows the NOK)
	vector<int> kindKnown; //which contributor index the knowns have
	vector<int> kindUnknown; //which contributor index the unknown have
	int NOC; //number of contributuros
	int NOK; //number of known contributors
	int NOU; //number of unknown contributors
    //The constructor prepares the variables with necessary data (arguments in constructor)

	//added for related individuals
	vector<int> kindRel; //which contributor index the related individual has

	public: 
	vector< vector<double> > cumprobvals; //nAlleles x 2 matrix with cumulative probabilities (used in function	calccumval)
	
	EFMmarker(int *NOC2, int *NOK2, int *knownGind, int *nR, int *nA, double *peaksLong, double *freqsLong, double *nTypedLong, double *maTypedLong, double *basepairLong, int *BWstuttindLong, int *FWstuttindLong, int *nPS, int *BWPstuttindLong, int *FWPstuttindLong, int *relGind) { 
		//NOC2 = number of contrs (known + unknown)
		//knownGind = NOC long vector with genotype index of known contributors
		//nR = number of replicates (needed to traverse peak height data)
		//nA = number of alleles. NB: Assuming same length for all replicates!!!!
		//peaksLong = longvector of Intensity values
		//freqsLong = longvector of frequency values
		//nTypedvec = longvector number of typed alleles (total)
		//maTypedLong = longvector of number of typed alleles of type a
		//basepairLong = vectorized vector of the already adjusted base pairs (x-150)/100 
		//BWstuttindLong = vectorized vector of backward stutter indices (gives index of what allele it receive stutter from  (Allele 1=index0))
		//FWstuttindLong = vectorized vector of forward stutter indices (gives index of what allele it receive stutter from  (Allele 1=index0))
		//nPS = number of potential stutters per marker
		//BWPstuttindLong = vectorized vector of backward stutter indices (gives index of what allele it receive stutter from  (Allele 1=index0))
		//FWPstuttindLong = vectorized vector of forward stutter indices (gives index of what allele it receive stutter from  (Allele 1=index0))
		//relGind = NOC long vector with genotype index of related contributors
		
		NOC = *NOC2; //copy number of contributors
		NOK = *NOK2; //copy number of contributors
		NOU = NOC-NOK; //number of unknowns
		nRep = *nR; //set number of replicates
		nAlleles = *nA; //set number of alleles
		nPstutters = *nPS; //set number of potential stutters (must be gloabal)
		nTyped = *nTypedLong; //set number of total prev. typed
		Yvec.assign(nRep*nAlleles,0.0); //init PH-vector
		Fvec.assign(nAlleles,0.0); //init vector
		BWvec.assign(nAlleles,0); //init vector
		FWvec.assign(nAlleles,0); //init vector
		maTypedvec.assign(nAlleles,0); //init vector
		bpvec.assign(nAlleles,0); //init vector
		for(i=0; i<nAlleles; i++) { //for each alleles (observed + dropout)
			Yvec[i] = peaksLong[i]; //copy value
			Fvec[i] = freqsLong[i]; //copy value
			BWvec[i] = BWstuttindLong[i]; //copy value
			FWvec[i] = FWstuttindLong[i]; //copy value
			maTypedvec[i] = maTypedLong[i]; //copy value
			bpvec[i] = basepairLong[i]; //copy value
		}
		for(i=nAlleles; i<(nRep*nAlleles); i++) { //include remaining PHs
			Yvec[i] = peaksLong[i]; //copy value
		}
		
		BWPvec.assign(nPstutters,0); //init vector
		FWPvec.assign(nPstutters,0); //init vector
		for(i=0; i<nPstutters; i++) { //for each unobserved potential stutters (BW+FW)
			BWPvec[i] = BWPstuttindLong[i]; //copy value
			FWPvec[i] = FWPstuttindLong[i]; //copy value
		}		
		
		//Create contribution matrix (1 genotype): Could be as part of constructor to save time
		numGenos1p = int(nAlleles*(nAlleles + 1) / 2); //get Genotype outcome
		outG1contr.assign(numGenos1p, vector<double>(nAlleles, 0.0)); //init nG1xnA matrix (contribution matrix). Indicating what alleles that are contributoed
		outG1mat.assign(numGenos1p, vector<int>(2, 0)); //init nG1x2 matrix (allele names as indices 0,...,nA-1)
		cc = 0; //counter oveer all
		for (i = 0; i < nAlleles; i++) {
			for (j = i; j < nAlleles; j++) { 
				outG1mat[cc][0] = i; //include index
				outG1mat[cc][1] = j; //include index
				outG1contr[cc][i] += 1.0; //insert contr at allele i
				outG1contr[cc][j] += 1.0; //insert contr at allele j
				cc++; //iterate to next genotype outcome
			}
		}
		
		//Prepare vector for known contributors (and also unknown): Need to know positions!
		GindKnown.assign(NOK,0); //genotype index in vector 
		kindKnown.assign(NOK,0); //contributor index in vector
		kindUnknown.assign(NOU,0); //contributor index in vector
		kindRel.assign(NOC,-1); //genotype index in vector, one index for each contributors  (-1 means no related)
		
		cc = 0; //counters for known
		j = 0; //counter for unknowns
		for(i=0; i<NOC; i++) { //for each contributors:
			kindRel[i] = relGind[i]; //copy genotype index 		
			if(knownGind[i]>=0) { //If contributor is known (genotype given)
				GindKnown[cc] = knownGind[i]; //copy genotype index
				kindKnown[cc] = i; //insert contributor index
				cc++; //update counter for knowns
			} else { //if contributor is unknown (genotype not given)
				kindUnknown[j] = i; //insert contributor index
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

	double prob_relUnknown(int Ugind, int Rgind, double *ibd2, double *fst2, double *maTypedvec2, double *nTyped2) {
		//Ugind = genotype index of unknown indiviudal 
		//Rgind = genotype index of related indiviudal 		
		//ibd2 = ibd vector (NOC long)
		double genoSum; //used to sum the genotype probability tfor unknowns		
		int aind = outG1mat[Ugind][0]; //get allele index of genotype g_1
		int bind = outG1mat[Ugind][1]; //get allele index of genotype g_2
		bool Uhom = aind==bind; //boolean of whether unknown genotype is homozygote
		
		//First step: Calculate random match probability of unrelated situation:
		genoSum = prob_a(Fvec[aind],maTypedvec2[aind],*nTyped2,fst2); //init with prob 1st allele  
		if(Uhom) { //if unknown is homozygote					
			genoSum *= prob_a(Fvec[aind],maTypedvec2[aind]+1,*nTyped2+1,fst2); //calculate random match prob (always used) 					
		} else { //if unknown is heterozygote variant
			genoSum *= 2*prob_a(Fvec[bind],maTypedvec2[bind],*nTyped2+1,fst2); //calculate prob 2st allele  (and scale with 2)
		}
			
		//Extension with kappa-coefficient: SEE FORMULAS IN TABLE A.3 in Book "A forensic practicioners guide...."
		if( Rgind != -1 ) { //if related is specified (not -1 index)
			genoSum *= ibd2[0]; //multiply with kappa0
			if( Ugind==Rgind ) { //if Unknown and Related are same genotype
				genoSum+=ibd2[2]; //sum with kappa2
				
				if(Uhom) { //if unknown is homozygote
					genoSum+= prob_a(Fvec[aind],maTypedvec2[aind],*nTyped2,fst2)*ibd2[1] ; //multiply with kappa1 
				} else { //if unknown is heterozygote variant
					genoSum+= (prob_a(Fvec[aind],maTypedvec2[aind],*nTyped2,fst2)+prob_a(Fvec[bind],maTypedvec2[bind],*nTyped2,fst2))*ibd2[1]/2 ; //multiply with kappa1 						
				}	
				
			} else { //if not the same genotype we need to check overlap (a,b)~(c,d)
				bool A1eq = aind==outG1mat[Rgind][0]; //check if a=c 
				bool A2eq = aind==outG1mat[Rgind][1]; //check if a=d 
				bool B1eq = bind==outG1mat[Rgind][0]; //check if b=c 
				bool B2eq = bind==outG1mat[Rgind][1]; //check if b=d 

				if( A1eq || A2eq || B1eq || B2eq) { //if any overlap (one shared allele): THERE ARE 3 OUTCOME!!
					if(Uhom) { //if Unknown is homozygote we know what allele to use
						genoSum += prob_a(Fvec[aind],maTypedvec2[aind],*nTyped2,fst2)*ibd2[1]/2; //multiply with kappa1 	
					} else { //if Unknown is heterozygote we need to found the non-overlapping allele 
						bool Rhom = outG1mat[Rgind][0]==outG1mat[Rgind][1]; //check if related genotype is homozygote
						int cind; //temporary variable to select non-overlapping allele
						if(A1eq || A2eq) { //Allele bind is non-overlapping
							cind = bind;
						} else {
							cind = aind; //Allele aind is non-overlapping
						}
						if(Rhom) { //should not divide by 2 if related geno is homozygote
							genoSum += prob_a(Fvec[cind],maTypedvec2[cind],*nTyped2,fst2)*ibd2[1]; //multiply with kappa1 																			
						} else { //should divide by 2 if related geno is heterozygote
							genoSum += prob_a(Fvec[cind],maTypedvec2[cind],*nTyped2,fst2)*ibd2[1]/2; //multiply with kappa1 																												
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
		//mixprop = mixture proportions (length must be *NOU)
		//shape = 1/sig^2, where sigma=coefficient of PH amount
		//scale = mu*sig^2, where mu=expected PH amount
		//beta is degradation slope parameters
		//xiB0 is baseline of expected stutter proportion (backward stutter)
		//xiF0 is baseline of expected stutter proportion (forward stutter)
	    //ibd = kappa-coefficient for each contributors (vectorised)
		int	nGK = pow(numGenos1p, NOU); //get total number of combined genotypes (for unknown contrs)
		double invshape = (*sigma)*(*sigma);//get inverse of shape parameter
		double shape = 1/invshape; //get shape parameter
		double scale = (*mu)*invshape; //get scale parameter
		const double const1 = 1/(scale); //constant 1
	 	const double const2 = log(scale); //constant 2
			
		//Prepare variables:
		int cind,a,r;
		vector<double> shapevK(nAlleles, 0.0); //Precalculation for known contributors:
		vector<double> shapev0(nAlleles, 0.0); //degrad scaling of shape
		vector<double> dropin0(nAlleles*nRep, 0.0); //drop-in contribution per peak height
		vector<double> constY0(nAlleles*nRep, 0.0); //const1 times Yvec	
		for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (indices), also those PH=0
			shapev0[a] = exp( log(shape) + bpvec[a]*log(*beta) ); //scaling with degradation model (assumed already scaled)
			for(r = 0; r < nRep; r++) { //traverse each replicates (observed alleles indicated by PH)
				cind = a*nRep + r; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
				dropin0[cind] = log(*lambda) - (*lambda)*(Yvec[cind]- (*AT) ) + log(*prC) + log(Fvec[a]); //add to dropin prob to sum	
				constY0[cind] = Yvec[cind]*const1; 
			}
		}

		//Sum up contribution for each alleles (Taking into account mix proportions): ONLY CALCULATED FOR KNOWN CONTRIBUTORS INITIALLY
		int aa, kk;
		for (kk = 0; kk < NOK; kk++) { //for each known contributors 
			for (aa = 0; aa < nAlleles; aa++) { //traverse each alleles (indices)
				shapevK[aa] += outG1contr[GindKnown[kk]][aa] * mixprop[kindKnown[kk]]; //contr from contr k to allele a
			}
		}
		//CALCULATING LARGE SUM:
		double bigsum = 0.0;  //total sum over all genotypes
		
		#pragma omp parallel for reduction(+:bigsum) //shared(shapev0, dropin0,constY0) //perform parallel on summing up bigsum
		for (int largeIter = 0; largeIter < nGK; largeIter++) { //for each combined/joint genotype outcome

			//Following variables must be declared for each iteration (non-shared):
			int a,k,r; //used to traverse alleles(a), contributors(k) and replicates(r)
			int aind; //used as index for alleles in genotypes
			vector<int> jointGind(NOU, 0); //This will be the permuation contribution index for the unknown inds (directly corresponds to indices of Gmarg)
			vector<double> shapev = shapevK; //make a copy of existing shapevector
			vector<double> maTypedvec2 = maTypedvec; //Creating copy of counter (per allele)
			double nTyped2 = nTyped; //creating copy of total counter

			//jointGind = digits(largeIter, base = numGenos1p, pad = NOU); #Equivalent operaion
			double genoProd = 1.0; //calculating the genotype probability of the unknowns
			int modrest = largeIter; //used to keep remained after modulo (init as iter number)
			bool inserted = false; //boolean of whether all digits are inserted (Rest is zero padded)
			for (k = 0; k < NOU; k++) { //for each unknown contributors (summing up wrt both contr (outG1contr) and mx (mixprop)): Need each contr to derive shapev
				if (!inserted) { //if not all digits inserted
					if ( k>0 ) {
						modrest = int((modrest - jointGind[k - 1]) / numGenos1p); //extract remaining, divide to get to next digit (necessary with int converion?)
					}
					jointGind[k] = modrest % numGenos1p; //INSERT NUMBER: convert number to "numGenos1p" basis 	
					if (modrest < numGenos1p) { //check if rest is smaller than base (only run if not inserted)
						inserted = true;  //then all digits are inserted
					}	
				} //else { jointGind[k] = 0; //zero pad	}
		
				//Sum up contribution for each alleles (Taking into account mix proportions)
				for (a = 0; a < nAlleles; a++) { //traverse each alleles (indices)
					shapev[a] += outG1contr[jointGind[k]][a] * mixprop[kindUnknown[k]]; //contr from contr k to allele a. NOTICE THE k+*NOK shift!
				}
								
				//CALCULATE GENOTYPE PROBS OF UNKNOWNS (MAY BE RELATED != -1)
				if( kindRel[kindUnknown[k]] == -1 ) {
					for(a = 0;a<2; a++) {
						aind = outG1mat[jointGind[k]][a]; //get allele index of genotype g_a
						genoProd *= ((*fst)*maTypedvec2[aind] + (1-*fst)*Fvec[aind]) / (1 + (nTyped2-1)*(*fst)); 
						maTypedvec2[aind] += 1; //update allele count for particular genotype
						nTyped2 += 1; //update total count
					}
					if(outG1mat[jointGind[k]][0]!=outG1mat[jointGind[k]][1]) { // heterozygous variant
						genoProd *=2; //multiply by 2  if het. variant
					} 
				} else { //IF RELATED WE CALL ANOTHER FUNCTION
					// Multiply the genotype probability of unknown (potential) related individual (related index = kindUnknown[k])
					genoProd *= prob_relUnknown( jointGind[k], kindRel[kindUnknown[k]], (ibd + 3*kindUnknown[k]), fst, &(maTypedvec2[0]), &nTyped2);	 //Note: scale with 3 because ibd is a '3-long vector' per contributor						
					//LAST: UPDATE COUNTERS FOR ALLELES (a and b)
					maTypedvec2[outG1mat[jointGind[k]][0]] += 1; //update allele count for particular allele
					maTypedvec2[outG1mat[jointGind[k]][1]] += 1; //update allele count for particular allele
					nTyped2 += 2; //update total count
				}
			} //end for each unknown contributor

			//////////////////////////////////////////////
			//Calculating the inner sum-part begins here//
			//////////////////////////////////////////////

			//Scaling shape parameter with degrad model
			for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (indices), having PH>0
				shapev[a] *= shapev0[a]; //scaling with degradation model (assumed already scaled)
			}
			//Scaling shape parameter with stutter model:
			vector<double> shapev2 = shapev;  //make a copy of existing shapevector
			double shape1; //modified shape param (caused by stutters)
			
			//Stutter model calculations (for given allele)
			if( (*xiB0>smalltol) || (*xiF0>smalltol) ) { //only run if stutterparam is positive
				for (a = 0; a < (nAlleles-1); a++) { //traverse each alleles. Potential non-Qallele is traversed below
					shape1 = 0.0; //copy original contribution
										
					//CONSIDERING STUTTERS: Always running (with param equal zero)
					if(FWvec[a] != -1) { //if a-1 stutter of a exists (found in FWstutter info)
						shape1 = (*xiF0)*shapev[FWvec[a]]; //proportion obtained from allele a-1
					} 
					if(BWvec[a] != -1) { //if a+1 stutter of a exists (found in BWstutter info)
						shape1 += (*xiB0)*shapev[BWvec[a]]; //proportion obtained from allele a-1
					}
					shapev2[a] = (1- *xiB0 - *xiF0)*shapev[a] + shape1; //updated shape parameter (summing fraction lost and attained), proportion xiB/xiF always lost (possible to unseen potential stutters see later)
				} 
				
				//HANDLING SPECIAL SCENARIO WHEN LAST ALLELE IS NOT Q-allele
				a = (nAlleles-1); //set allele index as last allele
				if(FWvec[a] != -1 || BWvec[a] != -1) {  //IF THERE WAS STUTTER INDEX DEFINED
					shape1 = 0.0; //copy original contribution									
					if(FWvec[a] != -1) { shape1 = (*xiF0)*shapev[FWvec[a]]; } //proportion obtained from allele a-1
					if(BWvec[a] != -1) { shape1 += (*xiB0)*shapev[BWvec[a]]; } //proportion obtained from allele a+1
					shapev2[a] = (1- *xiB0 - *xiF0)*shapev[a] + shape1; //updated shape parameter (summing fraction lost and attained), proportion xiB/xiF always lost (possible to unseen potential stutters see later)				
				} //end special scenario
			}
			
			//Summing up contribution of each alleles:
			vector<bool> anyDropin(nRep,false); //use vector to indicate drop-in
			double logevidProb = 0.0; //evidence probability (weight)
			int cind; //cumulative allele indexing for PHs (for traversing over all replicates)
			for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (also consider the last allele since Q-allele may have PH last allee)
				for(r = 0; r < nRep; r++) { //traverse each replicates (observed alleles indicated by PH)
					cind = a*nRep + r; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
					if( Yvec[cind] > smalltol ){ //if PH>0 (this is same as PH>=AT since threshold has already been applied)
						if(shapev2[a] > smalltol) { //If contribution and PH>0  					//CONTRIBUTION SET (A)
							logevidProb += -lgamma(shapev2[a]) - shapev2[a]*const2 + (shapev2[a]-1)*log(Yvec[cind]) - Yvec[cind]*const1; //gamma/lgamma is found in cmath							
						} else { //If contribution and PH>0 	//DROPIN SET (C)
							anyDropin[r] = true; //there was at least one dropin
							//evidProb *= lambda*exp(lambda*(Yvec[a]-AT))*prC*Fvec[a]; //add to dropin prob to sum			
							logevidProb += dropin0[cind]; //log(*lambda) - (*lambda)*(Yvec[cind] - (*AT) ) + log(*prC) + log(Fvec[a]); //add to dropin prob to sum										
						}							
					} else if(shapev2[a] > smalltol) { //IF PH missing  and contribution(it't a dropout //CONTRIBUTION SET (B)
						try{ //high values of shape1 with low values of const1 causes error -> Need to handle the exception to avoid crash
							logevidProb += log( gamma_p(shapev2[a],(*AT)*const1) );  //Add log(dropout-probability)
						} catch(...) {
							logevidProb += log(0); //Add with -INF
						}
					}
				} //end for each replicates (r)
			} //end for each observed alleles
			
			//Need to add this block to calc unobserved potential stutters (done ones for all replicates, since this can be scaled)
			if( (*xiB0>smalltol) || (*xiF0>smalltol) ) { //only run if param is positive
				for (a = 0; a < nPstutters; a++) { //traverse for each unobserved potential stutters
					shape1 = 0.0; //set to zero by default
					if(FWPvec[a] != -1) { //if a-1 stutter of a exists (found in FWPstutter info)
						shape1 = (*xiF0)*shapev[FWPvec[a]]; //proportion obtained from allele a-1
					} 
					if(BWPvec[a] != -1) { //if a+1 stutter of a exists (found in BWPstutter info)
						shape1 += (*xiB0)*shapev[BWPvec[a]]; //proportion obtained from allele a+1
					}
					if( shape1>smalltol) { //if there was a contribution, we need to calculate the prob for dropout
						try{ //high values of shape1 with low values of const1 causes error -> Need to handle the exception to avoid crash
							logevidProb += nRep*log( gamma_p(shape1,(*AT)*const1) ); //scaling with number replicates
						} catch(...) {
							logevidProb += log(0); //Add with -INF
						}
					}
				}
			}
			
			//Liklihood NO DROP-IN: Traverse each replicates (observed alleles indicated by PH)
			for(r = 0; r<nRep; r++) {
				if(!anyDropin[r]) { //if no drop-in occured for a replicate
					logevidProb += log(1-*prC); //
				}
			}
			bigsum += exp(logevidProb)*genoProd; //add temp res to big sum
		} //end for each combinations (bigsum)
	
		return log(bigsum); //return function
	}
	
	
	double calcLogLikGammaMarkerUnknown(double *fst,double *AT, double *lambda, double *prC, double *mixprop, double *mu, double *sigma, double *beta) {
		int	nGK = pow(numGenos1p, NOU); //get total number of combined genotypes (for unknown contrs)
		double invshape = (*sigma)*(*sigma);//get inverse of shape parameter
		double shape = 1/invshape; //get shape parameter
		double scale = (*mu)*invshape; //get scale parameter
		const double const1 = 1/(scale); //constant 1
	 	const double const2 = log(scale); //constant 2
			
		//Precalculation for known contributors:
		vector<double> shapevK(nAlleles, 0.0); //contribution vector of each alleles after taking into account mx

		//Prepare variables not needed in parallel:
		vector<double> shapev0(nAlleles, 0.0); //contribution vector of each alleles after taking into account mx
		vector<double> dropin0(nAlleles, 0.0); //drop-in contribution per peak height
		vector<double> constY0(nAlleles, 0.0); //PH vector multiplied with constant
		for (int a = 0; a < nAlleles; a++) { //traverse each observed alleles (indices), having PH>0
			shapev0[a] = exp( log(shape) + bpvec[a]*log(*beta) ); //scaling with degradation model (assumed already scaled)
			dropin0[a] = log(*lambda) - (*lambda)*(Yvec[a]- (*AT) ) + log(*prC) + log(Fvec[a]); //add to dropin prob to sum	
			constY0[a] = Yvec[a]*const1; 
		}

		//CALCULATING LARGE SUM:
		double bigsum = 0.0;  //total sum over all genotypes
		
		#pragma omp parallel for reduction(+:bigsum) //shared(shapev0, dropin0,constY0) //perform parallel on summing up bigsum
		for (int largeIter = 0; largeIter < nGK; largeIter++) { //for each combined/joint genotype outcome

			//Following variables must be declared for each iteration (non-shared):
			int a; //used to traverse alleles
			int k; //used to traverse contributors
			int aind; //used as index for alleles in genotypes
			vector<int> jointGind(NOU, 0); //This will be the permuation part (directly corresponds to indices of Gmarg)
			vector<double> shapev(nAlleles, 0.0); //contribution vector of each alleles after taking into account mx
			vector<double> maTypedvec2 = maTypedvec; //Creating copy of counter (per allele)
			double nTyped2 = nTyped; //creating copy of total counter

			//jointGind = digits(largeIter, base = numGenos1p, pad = NOC); #Equivalent operaion
			double genoProd = 1.0; //calculating the genotype probability of the unknowns
			int modrest = largeIter; //used to keep remained after modulo (init as iter number)
			bool inserted = false; //boolean of whether all digits are inserted (Rest is zero padded)
			for (k = 0; k < NOU; k++) { //for each contributor (summing up wrt both contr (outG1contr) and mx (mixprop)
				if (!inserted) { //if not all digits inserted
					if ( k>0 ) {
						modrest = int((modrest - jointGind[k - 1]) / numGenos1p); //extract remaining, divide to get to next digit (necessary with int converion?)
					}
					jointGind[k] = modrest % numGenos1p; //INSERT NUMBER: convert number to "numGenos1p" basis 						
					if (modrest < numGenos1p) { //check if rest is smaller than base (only run if not inserted)
						inserted = true;  //then all digits are inserted
					}	
				} //else { jointGind[k] = 0; //zero pad	}

				//Sum up contribution for each alleles (Taking into account mix proportions)
				for (a = 0; a < nAlleles; a++) { //traverse each alleles (indices)
					shapev[a] += outG1contr[jointGind[k]][a] * mixprop[k]; //contr from contr k to allele a
				}
				
				//Calculing the genotype probababilities where fst is taken into account (Balding-Nicholsen correction): (fst*ma  + (1-fst)*pa)/(1+ (n-1)*fst)
				for(a = 0;a<2; a++) {
					aind = outG1mat[jointGind[k]][a]; //get allele index of genotype g_a
					genoProd *= ((*fst)*maTypedvec2[aind] + (1-*fst)*Fvec[aind]) / (1 + (nTyped2-1)*(*fst)); 
					maTypedvec2[aind] += 1; //update allele count
					nTyped2 += 1; //update total count
				}

				if(outG1mat[jointGind[k]][0]!=outG1mat[jointGind[k]][1]) { // heterozygous variant
					genoProd *=2; //multiply by 2  if het. variant
				}
			} //end for each contributor

			//Calculating the inner sum-part begins here:
			//Scaling shape parameter with degrad model:
			for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (indices), having PH>0
				shapev[a] *= shapev0[a]; //scaling with degradation model (assumed already scaled)
			}
			
			//Summing up contribution of each alleles:
			bool anyDropin = false;
			double logevidProb = 0.0; //evidence probability (weight)
			double shape1; //modified shape param (caused by stutters)
			for (a = 0; a < (nAlleles-1); a++) { //traverse each observed alleles (indices), having PH>0 (no one are dropouts)				
				shape1 = shapev[a]; //copy original contribution	
				if( shape1>smalltol ) { //If contribution 
					logevidProb += -lgamma(shape1) - shape1*const2 + (shape1-1)*log(Yvec[a]) - constY0[a]; //gamma/lgamma is found in cmath							
				} else { //no contribution and PH present
					anyDropin = true; //there was at least one dropin
					logevidProb += dropin0[a]; //add to dropin prob to sum			
				} 
			} //end for each observed alleles
						
			//Calculate for the dropout allele (also including potential stutters:
			if(shapev[nAlleles-1]>smalltol) { //only if any contributors
				try{ //high values of shape1 with low values of const1 causes error -> Need to handle the exception to avoid crash
					logevidProb += log( gamma_p(shapev[nAlleles-1],(*AT)*const1) ); 
				} catch(...) {
					logevidProb += log(0); //Add with -INF
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
	
	double calcLogLikGammaMarkerUnknownStutter(double *fst,double *AT, double *lambda, double *prC, double *mixprop, double *mu, double *sigma, double *beta, double *xiB0, double *xiF0) {
		int	nGK = pow(numGenos1p, NOU); //get total number of combined genotypes (for unknown contrs)
		double invshape = (*sigma)*(*sigma);//get inverse of shape parameter
		double shape = 1/invshape; //get shape parameter
		double scale = (*mu)*invshape; //get scale parameter
		const double const1 = 1/(scale); //constant 1
	 	const double const2 = log(scale); //constant 2
			
		//Precalculation for known contributors:
		vector<double> shapevK(nAlleles, 0.0); //contribution vector of each alleles after taking into account mx

		//Prepare variables not needed in parallel:
		vector<double> shapev0(nAlleles, 0.0); //contribution vector of each alleles after taking into account mx
		vector<double> dropin0(nAlleles, 0.0); //drop-in contribution per peak height
		vector<double> constY0(nAlleles, 0.0); //PH vector multiplied with constant
		for (int a = 0; a < nAlleles; a++) { //traverse each observed alleles (indices), having PH>0
			shapev0[a] = exp( log(shape) + bpvec[a]*log(*beta) ); //scaling with degradation model (assumed already scaled)
			dropin0[a] = log(*lambda) - (*lambda)*(Yvec[a]- (*AT) ) + log(*prC) + log(Fvec[a]); //add to dropin prob to sum	
			constY0[a] = Yvec[a]*const1; 
		}

		//CALCULATING LARGE SUM:
		double bigsum = 0.0;  //total sum over all genotypes
		
		#pragma omp parallel for reduction(+:bigsum) //shared(shapev0, dropin0,constY0) //perform parallel on summing up bigsum
		for (int largeIter = 0; largeIter < nGK; largeIter++) { //for each combined/joint genotype outcome

			//Following variables must be declared for each iteration (non-shared):
			int a; //used to traverse alleles
			int k; //used to traverse contributors
			int aind; //used as index for alleles in genotypes
			vector<int> jointGind(NOU, 0); //This will be the permuation part (directly corresponds to indices of Gmarg)
			vector<double> shapev(nAlleles, 0.0); //contribution vector of each alleles after taking into account mx
			vector<double> maTypedvec2 = maTypedvec; //Creating copy of counter (per allele)
			double nTyped2 = nTyped; //creating copy of total counter

			//jointGind = digits(largeIter, base = numGenos1p, pad = NOC); #Equivalent operaion
			double genoProd = 1.0; //calculating the genotype probability of the unknowns
			int modrest = largeIter; //used to keep remained after modulo (init as iter number)
			bool inserted = false; //boolean of whether all digits are inserted (Rest is zero padded)
			for (k = 0; k < NOU; k++) { //for each contributor (summing up wrt both contr (outG1contr) and mx (mixprop)
				if (!inserted) { //if not all digits inserted
					if ( k>0 ) {
						modrest = int((modrest - jointGind[k - 1]) / numGenos1p); //extract remaining, divide to get to next digit (necessary with int converion?)
					}
					jointGind[k] = modrest % numGenos1p; //INSERT NUMBER: convert number to "numGenos1p" basis 						
					if (modrest < numGenos1p) { //check if rest is smaller than base (only run if not inserted)
						inserted = true;  //then all digits are inserted
					}	
				} //else { jointGind[k] = 0; //zero pad	}

				//Sum up contribution for each alleles (Taking into account mix proportions)
				for (a = 0; a < nAlleles; a++) { //traverse each alleles (indices)
					shapev[a] += outG1contr[jointGind[k]][a] * mixprop[k]; //contr from contr k to allele a
				}
				
				//Calculing the genotype probababilities where fst is taken into account (Balding-Nicholsen correction): (fst*ma  + (1-fst)*pa)/(1+ (n-1)*fst)
				for(a = 0;a<2; a++) {
					aind = outG1mat[jointGind[k]][a]; //get allele index of genotype g_a
					genoProd *= ((*fst)*maTypedvec2[aind] + (1-*fst)*Fvec[aind]) / (1 + (nTyped2-1)*(*fst)); 
					maTypedvec2[aind] += 1; //update allele count
					nTyped2 += 1; //update total count
				}

				if(outG1mat[jointGind[k]][0]!=outG1mat[jointGind[k]][1]) { // heterozygous variant
					genoProd *=2; //multiply by 2  if het. variant
				}
			} //end for each contributor

			//Calculating the inner sum-part begins here:
			//Scaling shape parameter with degrad model:
			for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (indices), having PH>0
				shapev[a] *= shapev0[a]; //scaling with degradation model (assumed already scaled)
			}
			
			//Summing up contribution of each alleles:
			bool anyDropin = false;
			double logevidProb = 0.0; //evidence probability (weight)
			double shape1; //modified shape param (caused by stutters)
			for (a = 0; a < (nAlleles-1); a++) { //traverse each observed alleles (indices), having PH>0 (no one are dropouts). LAST INDEX IS A Q-allele (PH=0)
				
				//Stutter model calculations (for given allele)
				shape1 = 0.0; //copy original contribution
				//CONSIDERING STUTTERS: Always running (with param equal zero)?? 
				if(FWvec[a]!=-1) { //if a-1 stutter of a exists (found in FWstutter info)
					shape1 = (*xiF0)*shapev[FWvec[a]]; //proportion obtained from allele a-1
				} 
				if(BWvec[a]!=-1) { //if a+1 stutter of a exists (found in BWstutter info)
					shape1 += (*xiB0)*shapev[BWvec[a]]; //proportion obtained from allele a-1
				}
				shape1 += (1- *xiB0 - *xiF0)*shapev[a]; //updated shape parameter (summing fraction lost and attained), proportion xiB/xiF always lost (possible to unseen potential stutters see later)
	
				if( shape1>smalltol ) { //If contribution 
					logevidProb += -lgamma(shape1) - shape1*const2 + (shape1-1)*log(Yvec[a]) - constY0[a]; //gamma/lgamma is found in cmath							
				} else { //no contribution and PH present
					anyDropin = true; //there was at least one dropin
					logevidProb += dropin0[a]; //add to dropin prob to sum			
				} 
			} //end for each observed alleles
						
			//Calculate for the dropout allele (also including potential stutters:
			if(shapev[nAlleles-1]>smalltol) { //only if any contributors
				try{ //high values of shape1 with low values of const1 causes error -> Need to handle the exception to avoid crash
					logevidProb += log( gamma_p(shapev[nAlleles-1],(*AT)*const1) ); 
				} catch(...) {
					logevidProb += log(0); //Add with -INF
				}
			}
			//Need to add this block to calc unobserved potential stutters
			for (a = 0; a < nPstutters; a++) { //traverse for each unobserved potential stutters
				shape1 = 0.0; //set to zero by default
				if(FWPvec[a] != -1) { //if a-1 stutter of a exists (found in FWPstutter info)
					shape1 += (*xiF0)*shapev[FWPvec[a]]; //proportion obtained from allele a-1
				} 
				if(BWPvec[a] != -1) { //if a+1 stutter of a exists (found in BWPstutter info)
					shape1 += (*xiB0)*shapev[BWPvec[a]]; //proportion obtained from allele a+1
				}
				if( shape1 > smalltol ) { //if there was a contribution, we need to calculate the prob for dropout
					try{ //high values of shape1 with low values of const1 causes error -> Need to handle the exception to avoid crash
						logevidProb += log( gamma_p(shape1,(*AT)*const1) ); //
					} catch(...) {
						logevidProb += log(0); //Add with -INF
					}
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
		int	nGK = pow(numGenos1p, NOU); //get total number of combined genotypes (for unknown contrs)
		double invshape = (*sigma)*(*sigma);//get inverse of shape parameter
		double shape = 1/invshape; //get shape parameter
		double scale = (*mu)*invshape; //get scale parameter
		const double const1 = 1/(scale); //constant 1
	 	const double const2 = log(scale); //constant 2
			
		//Precalculation for known contributors:
		vector<double> shapevK(nAlleles, 0.0); //contribution vector of each alleles after taking into account mx

		//Sum up contribution for each alleles (Taking into account mix proportions): ONLY CALCULATED FOR KNOWN CONTRIBUTORS INITIALLY
		int aa, kk;
		for (kk = 0; kk < NOK; kk++) { //for each known contributors 
			for (aa = 0; aa < nAlleles; aa++) { //traverse each alleles (indices)
				shapevK[aa] += outG1contr[GindKnown[kk]][aa] * mixprop[kindKnown[kk]]; //contr from contr k to allele a
			}
		}

		//For each observations we traverse the modified "large sum" 
		cumprobvals.assign(nRep*nAlleles, vector<double>(2, -1.0)); //init (nRep*nAlleles) x 2 matrix (col1=yA,col2=T)
		
		//double i; //indicating which allele (for a replicate) that is used 
		for(i=0; i<(nRep*nAlleles); i++) { //Traverse all PHs (all reps for each alleles)
			if(Yvec[i] < *AT) {
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
					vector<int> jointGind(NOU, 0); //This will be the permuation contribution index for the unknown inds (directly corresponds to indices of Gmarg)
					vector<double> shapev = shapevK; //make a copy of existing shapevector
					vector<double> maTypedvec2 = maTypedvec; //Creating copy of counter (per allele)
					double nTyped2 = nTyped; //creating copy of total counter

					//jointGind = digits(largeIter, base = numGenos1p, pad = NOU); #Equivalent operaion
					double genoProd = 1.0; //calculating the genotype probability of the unknowns
					int modrest = largeIter; //used to keep remained after modulo (init as iter number)
					bool inserted = false; //boolean of whether all digits are inserted (Rest is zero padded)
					for (k = 0; k < NOU; k++) { //for each unknown contributors (summing up wrt both contr (outG1contr) and mx (mixprop)): Need each contr to derive shapev
						if (!inserted) { //if not all digits inserted
							if ( k>0 ) {
								modrest = int((modrest - jointGind[k - 1]) / numGenos1p); //extract remaining, divide to get to next digit (necessary with int converion?)
							}
							jointGind[k] = modrest % numGenos1p; //INSERT NUMBER: convert number to "numGenos1p" basis 	
							if (modrest < numGenos1p) { //check if rest is smaller than base (only run if not inserted)
								inserted = true;  //then all digits are inserted
							}	
						} //else { jointGind[k] = 0; //zero pad	}
		
						//Sum up contribution for each alleles (Taking into account mix proportions)
						for (a = 0; a < nAlleles; a++) { //traverse each alleles (indices)
							shapev[a] += outG1contr[jointGind[k]][a] * mixprop[kindUnknown[k]]; //contr from contr k to allele a. NOTICE THE k+*NOK shift!
						}
						
						//CALCULATE GENOTYPE PROBS OF UNKNOWNS (MAY BE RELATED != -1)
						if( kindRel[kindUnknown[k]] == -1 ) {
							for(a = 0;a<2; a++) {
								aind = outG1mat[jointGind[k]][a]; //get allele index of genotype g_a
								genoProd *= ((*fst)*maTypedvec2[aind] + (1-*fst)*Fvec[aind]) / (1 + (nTyped2-1)*(*fst)); 
								maTypedvec2[aind] += 1; //update allele count for particular genotype
								nTyped2 += 1; //update total count
							}
							if(outG1mat[jointGind[k]][0]!=outG1mat[jointGind[k]][1]) { // heterozygous variant
								genoProd *=2; //multiply by 2  if het. variant
							} 
						} else { //IF RELATED WE CALL ANOTHER FUNCTION
							// Multiply the genotype probability of unknown (potential) related individual (related index = kindUnknown[k])
							genoProd *= prob_relUnknown( jointGind[k], kindRel[kindUnknown[k]], (ibd + 3*kindUnknown[k]), fst, &(maTypedvec2[0]), &nTyped2); //Note: scale with 3 because ibd is a '3-long vector' per contributor	
							//LAST: UPDATE COUNTERS FOR ALLELES (a and b)
							maTypedvec2[outG1mat[jointGind[k]][0]] += 1; //update allele count for particular allele
							maTypedvec2[outG1mat[jointGind[k]][1]] += 1; //update allele count for particular allele
							nTyped2 += 2; //update total count
						}

						
					} //end for each unknown contributor

					//////////////////////////////////////////////
					//Calculating the inner sum-part begins here//
					//////////////////////////////////////////////


					//Scaling shape parameter with degrad model:
					for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (indices), having PH>0
						shapev[a] *= exp( log(shape) + bpvec[a]*log(*beta) ); //scaling with degradation model (assumed already scaled)
					}
					//Scaling shape parameter with stutter model:
					vector<double> shapev2 = shapev;  //make a copy of existing shapevector
					double shape1; //modified shape param (caused by stutters)
					
					//Stutter model calculations (for given allele)
					if( (*xiB0>smalltol) || (*xiF0>smalltol) ) { //only run if stutterparam is positive
						for (a = 0; a < (nAlleles-1); a++) { //traverse each alleles. Potential non-Qallele is traversed below
							shape1 = 0.0; //copy original contribution
							
							//CONSIDERING STUTTERS: Always running (with param equal zero)
							if(FWvec[a] != -1) { //if a-1 stutter of a exists (found in FWstutter info)
								shape1 = (*xiF0)*shapev[FWvec[a]]; //proportion obtained from allele a-1
							} 
							if(BWvec[a] != -1) { //if a+1 stutter of a exists (found in BWstutter info)
								shape1 += (*xiB0)*shapev[BWvec[a]]; //proportion obtained from allele a-1
							}
							shapev2[a] = (1- *xiB0 - *xiF0)*shapev[a] + shape1; //updated shape parameter (summing fraction lost and attained), proportion xiB/xiF always lost (possible to unseen potential stutters see later)
						} 
						
						//HANDLING SPECIAL SCENARIO WHEN LAST ALLELE IS NOT Q-allele
						a = (nAlleles-1); //set allele index as last allele
						if(FWvec[a] != -1 || BWvec[a] != -1) {  //IF THERE WAS STUTTER INDEX DEFINED
							shape1 = 0.0; //copy original contribution									
							if(FWvec[a] != -1) { shape1 = (*xiF0)*shapev[FWvec[a]]; } //proportion obtained from allele a-1
							if(BWvec[a] != -1) { shape1 += (*xiB0)*shapev[BWvec[a]]; } //proportion obtained from allele a+1
							shapev2[a] = (1- *xiB0 - *xiF0)*shapev[a] + shape1; //updated shape parameter (summing fraction lost and attained), proportion xiB/xiF always lost (possible to unseen potential stutters see later)				
						} //end special scenario
					}
					
					//Summing up contribution of each alleles:
					vector<bool> anyDropin(nRep,false); //use vector to indicate drop-in
					double logevidProb = 0.0; //evidence probability (weight)
					int cind; //cumulative allele indexing for PHs (for traversing over all replicates)
					for (a = 0; a < nAlleles; a++) { //traverse each observed alleles (also consider the last allele since Q-allele may have PH last allee)
						for(r = 0; r < nRep; r++) { //traverse each replicates (observed alleles indicated by PH)
							cind = a*nRep + r; //get index of PH (Y_rep,allele is vectorised as : Y11,Y21,Y31,Y12,Y22,Y32,... 
						
							if( cind!=i ) { //proceed as before if different index
								if( Yvec[cind] > smalltol ){ //if evaluating index was not ii and PH>0 
									if( shapev2[a] > smalltol) { //If contribution and PH>0  					//CONTRIBUTION SET (A)
										logevidProb += -lgamma(shapev2[a]) - shapev2[a]*const2 + (shapev2[a]-1)*log(Yvec[cind]) - Yvec[cind]*const1; //gamma/lgamma is found in cmath							
									} else { //If contribution and PH>0 	//DROPIN SET (C)
										anyDropin[r] = true; //there was at least one dropin
										logevidProb += log(*lambda) - (*lambda)*(Yvec[cind] - (*AT) ) + log(*prC) + log(Fvec[a]); //add to dropin prob to sum										
									}							
								}  else if(shapev2[a]>smalltol) { //IF PH missing  and contribution(it't a dropout //CONTRIBUTION SET (B)
									logevidProb += log( gamma_p(shapev2[a],(*AT)*const1) );  //Add log(dropout-probability)
								}
							} else { //modify if allele index is same (assures that Yvec[cind]>=(*AT) (see early in loop)
								if(shapev2[a]>smalltol) { //If contribution and PH>0  					//CONTRIBUTION SET (A)
									if(j==0) { //CASE OF PH
										logevidProb += log( gamma_p(shapev2[a],Yvec[cind]*const1) -  gamma_p(shapev2[a],(*AT)*const1) );  //evaluate for upper limit - lower limit
									} else if(j==1) { //in case of j=1:  //CASE OF maxY
										logevidProb += log( gamma_p(shapev2[a],(*maxY)*const1) - gamma_p(shapev2[a],(*AT)*const1) );  //evaluate for upper limit - lower limit										
									} 
									
								} else { //If contribution and PH>0 //DROPIN SET (C)
									anyDropin[r] = true; //there was at least one dropin
									logevidProb += log(*prC) + log(Fvec[a]); //add to cumulative dropin prob to sum	
									if(j==0) { //CASE OF PH		
										logevidProb += log( 1 - exp(- (*lambda)*(Yvec[cind]-*AT))); //calc logarithm of cumulative expression
									} else if(j==1) { //CASE OF maxY
										logevidProb += log( 1 - exp(- (*lambda)*(*maxY-*AT))); //calc logarithm of cumulative expression
									} 							
								}							
							}							
						} //end for each replicates (r)
					} //end for each observed alleles
						
					//Need to add this block to calc unobserved potential stutters (done ones for all replicates, since this can be scaled)
					if( (*xiB0>smalltol) || (*xiF0>smalltol) ) { //only run if param is positive
						for (a = 0; a < nPstutters; a++) { //traverse for each unobserved potential stutters
							shape1 = 0.0; //set to zero by default
							if(FWPvec[a] != -1) { //if a-1 stutter of a exists (found in FWPstutter info)
								shape1 = (*xiF0)*shapev[FWPvec[a]]; //proportion obtained from allele a-1
							} 
							if(BWPvec[a] != -1) { //if a+1 stutter of a exists (found in BWPstutter info)
								shape1 += (*xiB0)*shapev[BWPvec[a]]; //proportion obtained from allele a+1
							}
							if( shape1 > smalltol ) { //if there was a contribution, we need to calculate the prob for dropout
								logevidProb += nRep*log( gamma_p(shape1,(*AT)*const1) ); //scaling with number replicates
							}
						}
					}
					
					//Liklihood NO DROP-IN: Traverse each replicates (observed alleles indicated by PH)
					for(r = 0; r<nRep; r++) {
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
void loglikgammaC(double *logLik, int *NOC, int *NOK, int *knownGind, double *mixprop, double *mu, double *sigma, double *beta, double *xiB,double *xiF, double *AT, double *pC, double *lambda, double *fst,  int *nRep, int *nM, int *nA, double *peaksLong, double *freqsLong, double *nTypedLong, double *maTypedLong, double *basepairLong, int *BWstuttindLong, int *FWstuttindLong, int *nPS, int *BWPstuttindLong, int *FWPstuttindLong, int *maxThreads, int *isPhi, int *anyRel, int *relGind, double *ibd) {
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
	//nRep = number of replicates (a vector per marker, different markers may have different number of replicates)
	//nM = number of markers (per sample): M_s_m M_1,1
	//nA = number of alleles per markers (population alleles)
	//peaksLong = vectorized vector of intensities (PeakHeights/Coverage), vectorized over all replicates (missing are indicated as zero)
	//freqsLong = vectorized vector of frequencies
	//nTypedLong = number of total number of previous typed alleles (used for fst). Vectorized
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
	
	int numThreads = thread::hardware_concurrency();
	int useThreads = min(numThreads,*maxThreads);
	omp_set_num_threads(useThreads);  //set number of threads to use 
	
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
	int startIndMarker2=0;//start marker index (nRep[m] reps)
	const double smalltol = 1.0e-30; //a tiny number > 0 (avoiding zero roundoff errors)
	for(int m=0; m< *nM; m++) {		
		EFMmarker *marker = new EFMmarker( NOC,(NOK+m),(knownGind + (*NOC)*m), (nRep+m), (nA+m), (peaksLong+startIndMarker2),(freqsLong+startIndMarker1), (nTypedLong+m), (maTypedLong+startIndMarker1),(basepairLong+startIndMarker1),(BWstuttindLong+startIndMarker1),(FWstuttindLong+startIndMarker1), (nPS+m), (BWPstuttindLong + startIndPS),(FWPstuttindLong + startIndPS), (relGind+(*NOC)*m) ); //prepare new marker

		//if only unknown unrelated with 1 replicate observed, last allele is dropout: RUN FASTEST CODE VARIANTS (LESS TO RUN)
		if(NOK[m]==0 && nRep[m]==1 && *anyRel==0 && peaksLong[startIndMarker1 + nA[m]-1 ] < smalltol ) {  
			if( xiB[m]>smalltol || xiF[m]>smalltol) { //IF stutter model version (assumes BW/FW stutter prop is positive
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
		startIndMarker2 += nA[m]*nRep[m]; //get start position of marker m+1 (nRep), used only for Peaks only
		startIndPS += nPS[m]; //add number of potential stutters
		delete marker; //delete after each declaration
	}
	*logLik = jointloglik;
} //end main function



//Obtaining cumulative probabilities for all observations (used for PP plots)
void cumvalgammaC(double *pvalPH, double *pvalMAX, double *maxY, int *NOC, int *NOK, int *knownGind, double *mixprop, double *mu, double *sigma, double *beta, double *xiB,double *xiF, double *AT, double *pC, double *lambda, double *fst,  int *nRep, int *nM, int *nA, double *peaksLong, double *freqsLong, double *nTypedLong, double *maTypedLong, double *basepairLong, int *BWstuttindLong, int *FWstuttindLong, int *nPS, int *BWPstuttindLong, int *FWPstuttindLong, int *maxThreads, int *relGind, double *ibd) {
	//pval PH/AT/MAX  cumulative value vectors with size nReps*nAlleles' 
	//maxY	is maximum value to use

	int numThreads = thread::hardware_concurrency();
	int useThreads = min(numThreads,*maxThreads);
	omp_set_num_threads(useThreads);  //set number of threads to use 
	
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
	int startIndMarker2=0;//start marker index (nRep[m] reps)
	int i,j; //used to traverse alleles and replicates
	for(int m=0; m< *nM; m++) {		
		EFMmarker *marker = new EFMmarker( NOC,(NOK+m),(knownGind + (*NOC)*m), (nRep+m), (nA+m), (peaksLong+startIndMarker2),(freqsLong+startIndMarker1), (nTypedLong+m), (maTypedLong+startIndMarker1),(basepairLong+startIndMarker1),(BWstuttindLong+startIndMarker1),(FWstuttindLong+startIndMarker1), (nPS+m), (BWPstuttindLong + startIndPS),(FWPstuttindLong + startIndPS), (relGind+(*NOC)*m) ); //prepare new marker
		marker->calccumval(maxY, (fst+m), (AT+m), (lambda+m) ,(pC+m), &(mixprop2[0]), (mu + m),(sigma + m),(beta+m),(xiB+m),(xiF+m),ibd) ; //calculate: ONLY ONE CODE VARIANT FOR THIS
		for(i=0; i < nA[m]*nRep[m]; i++) { //For each alleles in marker
			j = startIndMarker2 + i; //get index over set
			*(pvalPH+j) = marker->cumprobvals[i][0]; //copy value for cumulative on PH
			*(pvalMAX+j) = marker->cumprobvals[i][1]; //copy value for cumulative on MAX
		}		
		startIndMarker1 += nA[m]; //get start position of marker m+1 (1 rep)
		startIndMarker2 += nA[m]*nRep[m]; //get start position of marker m+1 (nRep), used only for Peaks only
		startIndPS += nPS[m]; //add number of potential stutters
		delete marker; //delete after each declaration
	}
} //end main function

} //end external