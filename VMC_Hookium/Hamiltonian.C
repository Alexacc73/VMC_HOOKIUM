/*! \mainpage VMC HOOKIUM
 * \section intro_sec Introduction

*/

#include "myRandom.H"

#include <fstream>
#include <time.h>
#include <set>
#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <vector>
#include <math.h>
#include <iomanip>
#include <boost/math/special_functions/factorials.hpp>

/** Number of expansion coefficients to use*/
const int EXPAND = 5;

const int numWalkers = 100;
/**Thermalisation Phase */
const int numEquilSteps = 10000;

/** Accumulations of observables (Local energy) Phase */
const int numAccumSteps = 50000; 

/**Spring constant in the Hooke's potential in the Hamiltonian - is necessarily set to be 1/4 since this particular
value os the spring constant gives rise to the only analytically soluble Hamiltonian ground state. The exact analytical 
Energy is 2, for this k = 0.25.  */
const double kspring = 0.2500000000;

const double SQRT2 = 1.4142135623730950488 ;

/**Expansion coefficients for the eigenfunctions of the Hamiltonian (which uses a hermite polynomial basis set, which are
eigenfunctions of the quantum harmoinc oscillator). 5 terms are used to approximate the Hartree Fock single-particle Wavefunction. */
//const double WFCoeff [EXPAND] = {1.0};
const double WFCoeff [EXPAND] = {0.9940786692, 0.10824442, -0.00939263, 0.00156292, -0.000284251 } ;
// expand = 4 :::: {0.994077953, 0.10825043, -0.00940132, 0.00157574} ;

/**This coefficients multiplies the gaussian random number generator which is used to impose trial moves on the electrons.
Smaller values of DELTA result in smaller steps, and a higher acceptance probability of Monte-Carlo moves. DELTA is this tuned
in order to acheive an acceptance probability of approximately 50% . */
const double DELTA1 = 1.5;
const double DELTA2 = 1.0;



/**
* The nth Hermite polynomial:
* Simply defined via a recursion relation:
\f[
H_{n+1}(r) = 2xH_n (r) - 2nH_{n-1}(r) 
\f]
*/
double hermite_Nth(int N, double r){

	double H0 = 1.0;
	double H1 = 2*r;
	double H_n = 0;        // H_n
	double H_nm1 = 0;      // H_n_minus1
	double H_np1 = 0;      // H_n_plus1 
	if(N==0){
		return H0;
	}
	if(N==1){
		//std::cout << "INSIDE HERMITE ROUTINE, H1: r = " << r << std::endl;
		return H1;
	}
	H_nm1 = H0;
	H_n = H1;    // Currently only computed up to n = 2;
	double nfloat;
	for(int n = 1; n < N; n++ ){
		nfloat = n;
		H_np1 = 2*r*H_n - 2*nfloat*H_nm1 ;
		H_nm1 = H_n;
		H_n = H_np1;
	}
	return H_np1;
}


/**
* The first derivative of the nth Hermite polynomial, easily defined from the nth Hermite:
\f[
 H'_n (r) = 2nH_{n-1}(r)  
\f]
*/
double hermite_diff1_Nth(int N, double r){
	double result;
	double floatN = N;
	result = 2.0*floatN*hermite_Nth(N-1, r );
	return result;
}


/**
* The second derivative of the nth Hermite polynomial, again, easily defined from the nth Hermite:
\f[
  H''_n (r) = \frac{\partial}{\partial r} 2nH_{n-1}(r)  = 4n(n-1)H_{n-2} (r) 
\f]
*/
double hermite_diff2_Nth(int N, double r){
	double result;
	double floatN = N;
	result = 4*floatN*(floatN-1)*hermite_Nth(N-2, r);
	return result;
}


/**
* This is simply a constant coefficient for f(r), or "FofR_Kth":
\f[
 \beta _k = \frac{ \sqrt{2}}{ 2^k \sqrt{ (2k-1)! } (2\pi)^{ \frac{3}{4} } }
\f]
*/
double beta_Kth(int K){
	double beta;
	double denom;
	double Kfloat = K;
	double factorial = boost::math::factorial<double>(2.0*Kfloat - 1);
	beta = sqrt(2);
	denom = pow(2, Kfloat) * sqrt(factorial) * pow( (2.0*M_PI) , 0.75 );
	beta *= (1.0/denom);
	return beta;
}


/**
* This function defines \f$ f_k(r) \f$, which is itself a function of \f$ H_{2k-1}(r), \ \beta _{k} \f$ :
\f[
 f_k (r) = H_{2k-1}( \frac{r}{ \sqrt{2} } ) \frac {\beta _k}{r} 
\f]
*/
double FofR_Kth(int K, double r){
	//if(K == 1){
	//	std::cout << "----- 1 ----- r value = " << r << std::endl;
	//}
	int N = 2*K - 1;
	double Rinput = (r/SQRT2);
	double hermite_2km1 = hermite_Nth(N, Rinput );
	double beta = beta_Kth(K);
	double result = (hermite_2km1*beta) / r;
	return result;

}


/**
* This function returns the first derivative of \f$ f_k(r) \f$, which has the form:
\f[
 f'_k(r) = \beta _k \Big\{ \frac{1}{r} H'_{2k-1} (\frac{r}{\sqrt{2}}) - \frac{1}{r^2} H_{2k-1} (\frac{r}{\sqrt{2}}) ) \Big\} 
\f]
*/
double FofR_diff1_Kth(int K, double r){
	int N = 2*K - 1;
	double beta = beta_Kth(K);
	double Rinput = (r/SQRT2);
	double hermite_diff1 = hermite_diff1_Nth(N, Rinput );
	double hermite_2km1 = hermite_Nth(N, Rinput );
	double result = (beta/r)*( hermite_diff1 - (hermite_2km1/r) );
	return result;
}


/**
* This function returns the second derivative of \f$ f_k(r) \f$, which has the form:
\f[
f''_k(r) = \frac{\beta _k}{r^3} \Big\{  r^2  H''_{2k-1} (\frac{r}{\sqrt{2}}) - 2rH'_{2k-1}(\frac{r}{\sqrt{2}})
 + 2H_{2k-1} (\frac{r}{\sqrt{2}}) ) \Big\} 
\f]
*/
double FofR_diff2_Kth(int K, double r){
	int N = 2*K - 1;
	double beta = beta_Kth(K);
	double Rinput = (r/SQRT2);
	double hermite_2km1 = hermite_Nth(N, Rinput);
	double hermite_diff1 = hermite_diff1_Nth(N, Rinput);
	double hermite_diff2 = hermite_diff2_Nth(N, Rinput);
	double result = (beta/(r*r*r)) * ( r*r*hermite_diff2 - 2*r*hermite_diff1 + 2*hermite_2km1 );
	return result;

}


/**
* This function returns \f$ \phi _k (r) \f$, which represents the individual basisfunctions which sum to give the total 
* single-particle wave function. It has the form:
\f[
\phi _k (r) = f_k (r) exp\{ -\frac{r^2}{4} \} 
\f]
* Where the total wave-function isgiven by:
\f[
 \Psi (r) \approx \Psi _N (r) = \sum_{k = 1}^{N} c_k \phi _k (r) 
\f]
*/
double PHI_Kth(int K, double r){
	double FofR = FofR_Kth(K, r);
	double exponent = exp( -((r*r)/4.0) );
	double result = FofR*exponent;
	return result;
}


/**
* This function returns the second derivative of \f$ \phi _k (r) \f$, called the laplacian in this context of the Hamiltonian.
* It has the functional form:
\f[
\phi ''_k (r) = exp\{ -\frac{r^2}{4} \} \Big\{ f''_k(r) - 4rf'_k(r) + [4r^2 - 2]f(r) \Big\}
\f]
*/
double PHI_laplace_Kth(int K, double r){
	double exponent = exp(-(r*r)/4);
	double FofR = FofR_Kth(K, r);
	double FofR_diff1 = FofR_diff1_Kth(K, r);
	double FofR_diff2 = FofR_diff2_Kth(K, r);
	double result = 0.25*exponent*(4*FofR_diff2 - 4*r*FofR_diff1 + ( (r*r - 2.0)*FofR ) );
	return result;
}


/**
* This returns a WF, evaluated for a single electron (of the two) in the Hookium system, \f$ \Psi _1 (r_1) \f$.
* The total trail-WF is evaluated simply by the product of two of these functions :
\f[
\Psi _T (\{R_i\}) = \Psi _1 (r_1) \Psi _2 (r_2)
\f]
*/
double singleParticleWF(const int numTerms, double r ){
	double totalWF = 0;
	double partialPhiWF = 0;
	int K  = 0;
	for(int i = 0; i < numTerms; i++){
		K = i+1;
		partialPhiWF = PHI_Kth(K, r);
		//std::cout << "Partial Phi term " << K << " = " << partialPhiWF << std::endl;
		totalWF += partialPhiWF*WFCoeff[i];
	}
	//std::cout << "Total single-P WF = " << totalWF << '\n' << std::endl;
	return totalWF;
}


/**
* This function returns the probabilty weighting, upon moving an electron in space via a stochastic move.
* It returns:
\f[
\rho(R) = | \frac{\Psi(r')}{\Psi(r)} |^2
\f]
*/
double probabiltyWeight(double r1, double r2, double r1Trial, double r2Trial){
	double WFup = singleParticleWF(EXPAND, r1);
	double WFdown = singleParticleWF(EXPAND, r2);
	double WFtotal = WFup * WFdown;

	double WFupTrial = singleParticleWF(EXPAND, r1Trial);
	double WFdownTrial = singleParticleWF(EXPAND, r2Trial);
	double WFtotalTrial = WFupTrial * WFdownTrial;
	//std::cout << "CURRENT WF = " << WFtotal << std::endl;
	//std::cout << "TRIAL WF = " << WFtotalTrial << std::endl;
	double probRatio = (WFtotalTrial/WFtotal);
	probRatio *= probRatio;
	return probRatio;
}


/**
* Returns the instantaneous Hamiltonian eigenvalue for a specific local state of the wavefunction. Returns:
\f[
\hat{H} \Psi _T (R)
\f]
*/
double hamiltonianHookium(const int numTerms, double r1, double r2, double r_12){
	double WFup = singleParticleWF(numTerms, r1);
	//std::cout << "WFup = " << WFup << std::endl;
	double WFdown = singleParticleWF(numTerms, r2);
	//std::cout << "Wdown = " << WFdown << std::endl;
	double WFtotal = WFup*WFdown;

	//std::cout << "WFtotal = " << WFtotal << std::endl;
	double laplaceWFup = 0;
	double partialUp = 0;
	double laplaceWFdown = 0;
	double partialDown = 0;
	int K = 0;
	for(int i = 0; i < numTerms; i++){
		K = i+1;
		partialUp = partialDown = 0;
		partialUp = PHI_laplace_Kth(K, r1) * WFCoeff[i];
		partialDown = PHI_laplace_Kth(K, r2) * WFCoeff[i];
		laplaceWFup += partialUp;
		laplaceWFdown += partialDown;
	}

	double result;
	result = -0.5*(laplaceWFup*WFdown + laplaceWFdown*WFup) + WFtotal*(0.5*kspring*(r1*r1 + r2*r2) + 1.0/r_12) ;
	//std::cout << "Two laplacian components = " << -0.5*(laplaceWFup*WFdown + laplaceWFdown*WFup) << std::endl;
	//std::cout << "Electron interaction parts = " << 0.5*kspring*(r1*r1 + r2*r2) + 1/r_12 << std::endl;
	//std::cout << "SUMMATION of both - result = " << result << std::endl;
	return result;
}


/**
* Function used to determine the Local Energy, \f$ E_L \f$, which gives the stochastic energy estimator summed over mutiple
* Monte Carlo Cycles, it has the simple form:
\f[
E_L = \frac{ \hat{H} \Psi _T (R)}{\Psi _T (R)}
\f]
*/
double localEnergy(const int numTerms, double r1, double r2, double r_12){
	double hamil = 0; 
	double WFtotal = 0; 
	double result = 0;
    hamil = hamiltonianHookium(numTerms, r1, r2, r_12);
	WFtotal = singleParticleWF(numTerms, r1) * singleParticleWF(numTerms, r2) ;
	result = hamil/WFtotal;
	//std::cout << "Hamiltonian = " << hamil << '\n';
	//std::cout << "TOTAL WF = " << WFtotal << '\n';
	//std::cout << "RESULT of ratio = " << result <<std::endl;
	return result;
}





/**
------------------------------------------------
V M C  -  R O U T I N E S 
------------------------------------------------
*/

double getRvector(double x, double y, double z){
	double R;
	R = sqrt(x*x + y*y + z*z);
	return R;
}

double getR12vector(double x1, double y1, double z1, double x2, double y2, double z2){
	double r_12 = 0;
	r_12 = (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2);
	r_12 = sqrt(r_12);
	return r_12;
}

void initialiseWalker(std::vector< std::tuple<double, double, double> >& vector1,
	                  std::vector< std::tuple<double, double, double> >& vector2){
	double x1, y1, z1;
	double x2, y2, z2;
	x1 = (std_rand() - 0.5)*2 ;
	y1 = (std_rand() - 0.5)*2 ;
	z1 = (std_rand() - 0.5)*2 ;
	x2 = (std_rand() - 0.5)*2 ;
	y2 = (std_rand() - 0.5)*2 ;
	z2 = (std_rand() - 0.5)*2 ;
	vector1.push_back( std::make_tuple(x1, y1, z1) );
	vector2.push_back( std::make_tuple(x2, y2, z2) );
}


void metropStep(std::vector< std::tuple<double, double, double> >& vector1,
	            std::vector< std::tuple<double, double, double> >& vector2,
	            std::vector<double>& energyList,
	            int& successCounter,
	            int& walkerIDX,
	            bool& accum){
	double DELTA;
	if(accum == true){
		DELTA = DELTA2;
	}
	if(accum == false){
		DELTA == DELTA1;
	}
	double x1, y1, z1, x2, y2, z2;
	double x1trial, y1trial, z1trial, x2trial, y2trial, z2trial;
	double r1Trial, r2Trial;
	double r_12, r_12Trial;
	double r1Current, r2Current;
	double probTrial;
	double currentEnergy = 0;
	double newEnergy = 0;
	double weightedEnergy = 0;
	double metropolisRand = std_rand();
    x1 = std::get<0>(vector1[walkerIDX]);
    y1 = std::get<1>(vector1[walkerIDX]);
    z1 = std::get<2>(vector1[walkerIDX]);
    r1Current = getRvector(x1, y1, z1);
    x2 = std::get<0>(vector2[walkerIDX]);
    y2 = std::get<1>(vector2[walkerIDX]);
    z2 = std::get<2>(vector2[walkerIDX]);
    r2Current = getRvector(x2, y2, z2);
    r_12 = getR12vector(x1, y1, z1, x2, y2, z2);

    /*STEP 1: Propose a move from R --> R' */
	x1trial = x1 + DELTA*(std_rand() - 0.5 ) ;
	y1trial = y1 + DELTA*(std_rand() - 0.5 ) ;
	z1trial = z1 + DELTA*(std_rand() - 0.5 ) ;
	r1Trial = getRvector(x1trial, y1trial, z1trial);
	x2trial = x2 + DELTA*(std_rand() - 0.5 ) ;
	y2trial = y2 + DELTA*(std_rand() - 0.5 ) ;
	z2trial = z2 + DELTA*(std_rand() - 0.5 ) ;
	r2Trial = getRvector(x2trial, y2trial, z2trial);
	r_12Trial = getR12vector(x1trial, y1trial, z1trial, x2trial, y2trial, z2trial );

    /*STEP 2: Calculate the probability weighting:*/
	probTrial = probabiltyWeight(r1Current, r2Current, r1Trial, r2Trial);

    /*STEP 3 (Accumulation phase only): Accumulate the contribution to the local energy
    weighted by the metropolis acceptance criteria*/
	if(accum == true){
		currentEnergy = localEnergy(EXPAND, r1Current, r2Current, r_12);
		newEnergy = localEnergy(EXPAND, r1Trial, r2Trial, r_12Trial);
		weightedEnergy = (probTrial*newEnergy) + (1.0-probTrial)*currentEnergy;
		energyList.push_back(weightedEnergy);
	}
	//std::cout << "Prob Ratio test = " << probTrial << std::endl;
	//std::cout<< "ProbTrial Value = " << probTrial << std::endl;

	/*STEP 4: Accept or Reject the new move*/
	if( probTrial > metropolisRand ){ // RANDOM SINGLE ELECTRON MOVE ACCEPTED
		//std::cout << "RANDOM SINGLE ELECTRON MOVE ACCEPTED" << std::endl;
		successCounter += 1;
		// UPDATE POSITIONS FROM THE SUCCESSFUL TRIAL MOVES
		std::get<0>(vector1[walkerIDX]) = x1 = x1trial;
		std::get<1>(vector1[walkerIDX]) = y1 = y1trial;
		std::get<2>(vector1[walkerIDX]) = z1 = z1trial;
		std::get<0>(vector2[walkerIDX]) = x2 = x2trial;
		std::get<1>(vector2[walkerIDX]) = y2 = y2trial;
		std::get<2>(vector2[walkerIDX]) = z2 = z2trial;

		if(accum = false){ // Only necessary to update the following values in the equilibration phase
			r1Current = r1Trial;
    		r2Current = r2Trial;
    		r_12 = r_12Trial;
    	}
	}
	

	////std::cout << "r_12 = " << r_12 << std::endl;
	if(accum == false){
		newEnergy = localEnergy(EXPAND, r1Current, r2Current, r_12 );
		energyList.push_back(newEnergy);
	}

	if(newEnergy > 10000){
		std::cout << "!!! HIGH ENERGY !!!" << std::endl;
		std::cout << "Energy = " << newEnergy << std::endl;
		std::cout << "|r2-r1| distance = " << r_12 << std::endl;
		std::cout << "r1 = " << r1Current << ", and r2 = " << r2Current << '\n' << std::endl;
	}

}





/**
------------------------------------------------
M A I N  -  S T A R T S  -  H E R E 
------------------------------------------------
*/
int main(void){
	srand(492830);
	std::cout << std::setprecision(10);
	clock_t start;
	clock_t end;
	

	/*
	start = clock();
	end = clock();
	double TIME = end-start;
	std::cout << "LOOP TIME = " << TIME/CLOCKS_PER_SEC << std::endl;
	*/

    

    std::vector<double> localEn_Equilibrate;
    std::vector<double> meanEnergiesEquilibrate;

	/**Local Energy:  \f$ E_l = \frac{H \Psi}{\Psi} \f$*/
    std::vector<double> localEnergyAccumulator;
    /**Square od the local energy (used to calculate the variance, etc.) */
    std::vector<double> E_lAccumSqrd;

    std::vector<double> meanEnergiesAccum;

	std::vector< std::tuple<double, double, double> > r1 ; 
	std::vector< std::tuple<double, double, double> > r2 ; 


    /** -------------> INITIALISE N WALKERS <---------------- */
    for(int i = 0; i<numWalkers; i++){ 
		initialiseWalker(r1, r2);
	}
	std::cout << "[x y z] = [" << std::get<0>(r1[0]) <<" "<< std::get<1>(r1[0]) <<" "<< std::get<2>(r1[0]) << "]" << std::endl;


	/**
	------------------------------------------------
	E Q U I L I B R A T I O N   -   S T E P S
	------------------------------------------------
	*/
	bool ACCUM1 = false;
    int success = 0;
    double walkerMeanEnergy = 0;


    for(int idx = 0; idx < numWalkers; idx++){
    	for(int i = 0; i<numEquilSteps; i++){ // Perform TWO ELECTRONS at a time!!!
    		metropStep(r1, r2, localEn_Equilibrate, success, idx, ACCUM1);
    		walkerMeanEnergy += localEn_Equilibrate[i];
    	}
    	walkerMeanEnergy *= ( 1.0/double(localEn_Equilibrate.size()) );
    	meanEnergiesEquilibrate.push_back(walkerMeanEnergy);
    	walkerMeanEnergy = 0;
    	localEn_Equilibrate.clear();
    }


	std::cout << "[x y z] = [" << std::get<0>(r1[0]) <<" "<< std::get<1>(r1[0]) <<" "<< std::get<2>(r1[0]) << "]" << std::endl;


    std::ofstream energyEquil;
    energyEquil.open("ENERGIES_VMC");
	double instantE;
	double MEAN = 0;
	for(int i = 0; i < meanEnergiesEquilibrate.size(); i++){
		instantE = meanEnergiesEquilibrate[i];
		energyEquil << meanEnergiesEquilibrate[i] << std::endl;
		MEAN += meanEnergiesEquilibrate[i] ;
	}


	energyEquil.close();
	std::cout << " ----- EQUILIBRATION COMPLETE -> DATA WRITTEN TO FILE ----- " << '\n' << std::endl;
	double sizefloat = meanEnergiesEquilibrate.size() ;
	MEAN *= (1.0/sizefloat);
	std::cout << "Mean Stochastic Local Energy after Equilibration = " << MEAN << std::endl;
	std::cout << "Number of successful monte-carlo moves = " << success << std::endl;
	std::cout << "Proportion of moves which were successful : " << success / ((float)(numEquilSteps*numWalkers)) << '\n' << std::endl; 

	/**
	------------------------------------------------
	A C C U M U L A T I O N   -   S T E P S
	------------------------------------------------
	*/
	std::cout << " ----- BEGIN ACCUMULATION PHASE -----" << std::endl;

	
	bool ACCUM2 = true;
	int successAccum = 0;
	double accumMean = 0;
	for(int idx = 0; idx < numWalkers; idx++){
		for(int mv = 0; mv < numAccumSteps; mv++){
			metropStep(r1, r2, localEnergyAccumulator, successAccum, idx, ACCUM2);
			accumMean += localEnergyAccumulator[mv];
		}
		accumMean = (accumMean / double(localEnergyAccumulator.size()) ) ;
		meanEnergiesAccum.push_back(accumMean);
		accumMean = 0;
		localEnergyAccumulator.clear(); 
	}

	std::ofstream energyAccumfile;
    energyAccumfile.open("ENERGIES_ACCUM_VMC");
	double instantEAccum;
	double globalAccumMean = 0;
	for(int i = 0; i < meanEnergiesAccum.size(); i++){
		instantEAccum = meanEnergiesAccum[i];
		energyAccumfile << meanEnergiesAccum[i] << std::endl;
		globalAccumMean += meanEnergiesAccum[i] ;
	}
	energyAccumfile.close();
	double numAccums = meanEnergiesAccum.size();
	globalAccumMean *= (1.0/numAccums);
	std::cout << "Mean Stochastic Local Energy after ACCUMULATION = " << globalAccumMean << std::endl;
	std::cout << "Number of successful monte-carlo moves in ACCUMULATION phase = " << successAccum << std::endl;
	

	

	//std::cout << "c_1 = " << pow((1 - (WFCoeff[1]*WFCoeff[1] + WFCoeff[2]*WFCoeff[2] + WFCoeff[3]*WFCoeff[3] + WFCoeff[4]*WFCoeff[4] )), 0.5) << std::endl;
    
    return 1 ; 
}