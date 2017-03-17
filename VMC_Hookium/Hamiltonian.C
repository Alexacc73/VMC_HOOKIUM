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

const int EXPAND = 5;
const int numWalkers = 50;
const int numEquilSteps = 5000000;

const double kspring = 0.2500000000;
const double SQRT2 = 1.4142135623730950488 ;
const double WFCoeff [EXPAND] = {0.9940786692, 0.10824442, -0.00939263, 0.00156292, -0.000284251 } ;
// expand = 4 :::: {0.994077953, 0.10825043, -0.00940132, 0.00157574} ;
const double DELTA = 0.3;



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
	result = 2*N*hermite_Nth(N-1, r );
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
	result = 4*N*(N-1)*hermite_Nth(N-2, r);
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
	double factorial = boost::math::factorial<double>(2*Kfloat - 1);
	beta = sqrt(2);
	denom = pow(2, Kfloat) * sqrt(factorial) * pow(2*M_PI, 0.75 );
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
	double hermite_2km1 = hermite_Nth(N, (r/SQRT2) );
	double beta = beta_Kth(K);
	double result = (hermite_2km1*beta) / r;
	if(K == 1){
		//std::cout << "for phi_1 beta = " << beta << std::endl;
		//std::cout << "H_2k-1 for phi_1 = " << hermite_2km1 << std::endl;
		//std::cout << "----- 2 -----r value = " << r << std::endl;
	}
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
	double hermite_diff1 = hermite_diff1_Nth(N, (r/SQRT2));
	double hermite_2km1 = hermite_Nth(N, (r/SQRT2));
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
	double hermite_2km1 = hermite_Nth(N, (r/SQRT2));
	double hermite_diff1 = hermite_diff1_Nth(N, (r/SQRT2));
	double hermite_diff2 = hermite_diff2_Nth(N, (r/SQRT2));
	double result = (beta/(r*r*r))*( r*r*hermite_diff2 - 2*r*hermite_diff1 + 2*hermite_2km1 );
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
	double result = exponent*(FofR_diff2 - r*FofR_diff1 + (( (r*r/4) - 0.5)*FofR) );
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

	//std::cout << "_>_>_>_>_>__> FIRST UP DONE _>_>_>_>_>" << std::endl;


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
	double r_12;
	r_12 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);
	r_12 = sqrt(r_12);
	return r_12;
}

void initialiseWalker(std::vector< std::tuple<double, double, double> >& vector1,
	                  std::vector< std::tuple<double, double, double> >& vector2){
	double x1, y1, z1;
	double x2, y2, z2;
	x1 = (std_rand() - 0.5)*0.2 ;
	y1 = (std_rand() - 0.5)*0.2 ;
	z1 = (std_rand() - 0.5)*0.2 ;
	x2 = (std_rand() - 0.5)*0.2 ;
	y2 = (std_rand() - 0.5)*0.2 ;
	z2 = (std_rand() - 0.5)*0.2 ;
	vector1.push_back( std::make_tuple(x1, y1, z1) );
	vector2.push_back( std::make_tuple(x2, y2, z2) );
}


void metropStep(std::vector< std::tuple<double, double, double> >& vector1,
	            std::vector< std::tuple<double, double, double> >& vector2,
	            std::vector<double>& energyList,
	            int& successCounter,
	            int& walkerIDX){
	double x1, y1, z1, x2, y2, z2;
	double x1trial, y1trial, z1trial, x2trial, y2trial, z2trial;
	double r1Trial, r2Trial;
	double r_12;
	double r1Current, r2Current;
	double probTrial;
	double newEnergy = 0;
	double metropolisRand = std_rand();
    x1 = std::get<0>(vector1[walkerIDX]);
    y1 = std::get<1>(vector1[walkerIDX]);
    z1 = std::get<2>(vector1[walkerIDX]);
    r1Current = getRvector(x1, y1, z1);
    x2 = std::get<0>(vector2[walkerIDX]);
    y2 = std::get<1>(vector2[walkerIDX]);
    z2 = std::get<2>(vector2[walkerIDX]);
    r2Current = getRvector(x2, y2, z2);

	x1trial = x1 + DELTA*gasdev() ;
	y1trial = y1 + DELTA*gasdev() ;
	z1trial = z1 + DELTA*gasdev() ;
	r1Trial = getRvector(x1trial, y1trial, z1trial);
	x2trial = x2 + DELTA*gasdev() ;
	y2trial = y2 + DELTA*gasdev() ;
	z2trial = z2 + DELTA*gasdev() ;
	r2Trial = getRvector(x2trial, y2trial, z2trial);


	probTrial = probabiltyWeight(r1Current, r2Current, r1Trial, r2Trial);
	//std::cout << "Prob Ratio test = " << probTrial << std::endl;
	//std::cout<< "ProbTrial Value = " << probTrial << std::endl;
	if( probTrial > metropolisRand ){ // RANDOM SINGLE ELECTRON MOVE ACCEPTED
		//std::cout << "RANDOM SINGLE ELECTRON MOVE ACCEPTED" << std::endl;
		successCounter += 1;
		std::get<0>(vector1[walkerIDX]) = x1 = x1trial;
		std::get<1>(vector1[walkerIDX]) = y1 = y1trial;
		std::get<2>(vector1[walkerIDX]) = z1 = z1trial;
		std::get<0>(vector2[walkerIDX]) = x2 = x2trial;
		std::get<1>(vector2[walkerIDX]) = y2 = y2trial;
		std::get<2>(vector2[walkerIDX]) = z2 = z2trial;
        // UPDATE POSITIONS FROM THE SUCCESSFUL TRIAL MOVES
    	r1Current = getRvector(x1, y1, z1);
    	r2Current = getRvector(x2, y2, z2);
	}
	

	r_12 = getR12vector(x1, y1, z1, x2, y2, z2);
	////std::cout << "r_12 = " << r_12 << std::endl;
	newEnergy = localEnergy(EXPAND, r1Current, r2Current, r_12 );
	energyList.push_back(newEnergy);
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
	srand(492831);
	std::cout << std::setprecision(10);
	clock_t start;
	clock_t end;
	

	/*
	start = clock();
	end = clock();
	double TIME = end-start;
	std::cout << "LOOP TIME = " << TIME/CLOCKS_PER_SEC << std::endl;
	*/
    std::vector<double> localEnergyAccumulator;
    std::vector<double> meanEnergies;
	std::vector< std::tuple<double, double, double> > r1 ; 
	std::vector< std::tuple<double, double, double> > r2 ; 

    for(int i = 0; i<numWalkers; i++){ // Initialise N walkers
		initialiseWalker(r1, r2);
	}
	std::cout << "[x y z] = [" << std::get<0>(r1[0]) <<" "<< std::get<1>(r1[0]) <<" "<< std::get<2>(r1[0]) << "]" << std::endl;


	/**
	------------------------------------------------
	E Q U I L I B R A T I O N   -   S T E P S
	------------------------------------------------
	*/
    int success = 0;
    double walkerMeanEnergy = 0;
    for(int idx = 0; idx < numWalkers; idx++){
    	for(int i = 0; i<numEquilSteps; i++){ // Perform TWO ELECTRONS at a time!!!
    		//std::cout << "STEP " << i << " ------------------------->" << std::endl; 

    		metropStep(r1, r2, localEnergyAccumulator, success, idx);
    		//std::cout << " " << std::endl;
    	}
    	for(int k = 0; k<localEnergyAccumulator.size(); k++){
    		walkerMeanEnergy += localEnergyAccumulator[k];
    	}
    	walkerMeanEnergy *= ( 1/float(localEnergyAccumulator.size()) );
    	meanEnergies.push_back(walkerMeanEnergy);
    	walkerMeanEnergy = 0;
    	localEnergyAccumulator.clear();
    }
	std::cout << "[x y z] = [" << std::get<0>(r1[0]) <<" "<< std::get<1>(r1[0]) <<" "<< std::get<2>(r1[0]) << "]" << std::endl;


	/**
	------------------------------------------------
	A C C U M U L A T I O N   -   S T E P S
	------------------------------------------------
	*/
	

    std::ofstream energyAccum;
    energyAccum.open("ENERGIES_VMC");
	double MEAN = 0;
	double instantE;
	for(int i = 0; i < meanEnergies.size(); i++){
		instantE = meanEnergies[i];
		energyAccum << meanEnergies[i] << std::endl;
		MEAN += meanEnergies[i] ;
	}
	float sizefloat = meanEnergies.size() ;
	MEAN *= (1/sizefloat);
	std::cout << "Mean Stochastic Local Energy = " << MEAN << std::endl;
	std::cout << "Number of successful monte-carlo moves = " << success << std::endl;
	std::cout << "Proportion of moves which were successful : " << success / ((float)(numEquilSteps*numWalkers)) << std::endl; 
	

	//std::cout << "c_1 = " << pow((1 - (WFCoeff[1]*WFCoeff[1] + WFCoeff[2]*WFCoeff[2] + WFCoeff[3]*WFCoeff[3] + WFCoeff[4]*WFCoeff[4] )), 0.5) << std::endl;
    
    return 1 ; 
}