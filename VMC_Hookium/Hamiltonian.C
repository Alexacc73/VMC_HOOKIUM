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

const double kspring = 0.25;
const double SQRT2 = 1.4142135623730950488 ;
const double WFCoeff [4] = {1, 0.10825043, -0.00940132, 0.00157574} ;
const double DELTA = 0.05;
const int EXPAND = 4;

/**
* The nth Hermite polynomial:
* Simply defined via a recursion relation:
\f[
H_{n+1}(r) = 2xH_n (r) - 2nH_{n-1}(r) 
\f]
*/
double hermite_Nth(int N, double r){
	double H0 = 1.0;
	double H1 = 2.0*r;
	double H_n = 0;        // H_n
	double H_nm1 = 0;      // H_n_minus1
	double H_np1 = 0;      // H_n_plus1 
	if(N==0){
		return H0;
	}
	if(N==1){
		//std::cout << "H1 USED!" << std::endl;
		//std::cout << H1 << std::endl;
		return H1;
	}
	H_nm1 = H0;
	H_n = H1;    // Currently only computed up to n = 2;
	double nfloat;
	for(int n = 1; n < N; n++ ){
		nfloat = n;
		H_np1 = 2.0*r*H_n - 2.0*nfloat*H_nm1 ;
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
	//std::cout << "FACTORIAL USED = " << factorial << std::endl;
	beta = pow(2, 0.5);
	denom = pow(2, Kfloat) * pow( factorial ,0.5) * pow(2*M_PI, (3.0/4.0) );
	beta *= 1/denom;
	return beta;
}


/**
* This function defines \f$ f_k(r) \f$, which is itself a function of \f$ H_{2k-1}(r), \ \beta _{k} \f$ :
\f[
 f_k (r) = H_{2k-1}( \frac{r}{ \sqrt{2} } ) \frac {\beta _k}{r} 
\f]
*/
double FofR_Kth(int K, double r){
	int N = 2*K - 1;
	double hermite_2km1 = hermite_Nth(N, (r/SQRT2) );
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
	double exponent = exp(-(r*r)/4);
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
	double result = exponent*(FofR_diff2 - 4*r*FofR_diff1 + ((4*r*r - 2)*FofR) );
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
		totalWF += partialPhiWF*WFCoeff[i];
	}
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
	double probRatio = (WFtotalTrial/WFtotal)*(WFtotalTrial/WFtotal) ;
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
	int K;
	for(int i = 0; i < numTerms; i++){
		K = i+1;
		partialUp = PHI_laplace_Kth(K, r1) * WFCoeff[i];
		partialDown = PHI_laplace_Kth(K, r2) * WFCoeff[i];
		laplaceWFup += partialUp;
		laplaceWFdown += partialDown;
	}
	//std::cout << "LAPLACE up = " << laplaceWFup << std::endl;
	//std::cout << "LAPLACE down = " << laplaceWFdown << std::endl;
	double result;
	//std::cout << "LASt TERM = " << WFtotal/(fabs(r2-r1)) << std::endl;
	result = -0.5*(laplaceWFup*WFdown + laplaceWFdown*WFup) + 0.5*kspring*WFtotal*(r1*r1 + r2*r2) + WFtotal/r_12;
	//std::cout << "END HAMIL RESULT = " << result << std::endl;
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
	double hamil = hamiltonianHookium(numTerms, r1, r2, r_12);
	double WFtotal = singleParticleWF(numTerms, r1) * singleParticleWF(numTerms, r2) ;
	double result = hamil/WFtotal;
	return result;
}





/**
------------------------------------------------
V M C  -  R O U T I N E S 
------------------------------------------------
*/

double getRvector(double x, double y, double z){
	double R;
	R = pow((x*x + y*y + z*z), 0.5);
	return R;
}

double getR12vector(double x1, double y1, double z1, double x2, double y2, double z2){
	double r_12;
	r_12 = (x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1);
	r_12 = pow(r_12, 0.5);
	return r_12;
}

void initialiseWalker(std::vector< std::tuple<double, double, double> >& vector1,
	                  std::vector< std::tuple<double, double, double> >& vector2){
	double x1, y1, z1;
	double x2, y2, z2;
	x1 = std_rand() - 0.5;
	y1 = std_rand() - 0.5;
	z1 = std_rand() - 0.5;
	x2 = std_rand() - 0.5;
	y2 = std_rand() - 0.5;
	z2 = std_rand() - 0.5;
	vector1.push_back( std::make_tuple(x1, y1, z1) );
	vector2.push_back( std::make_tuple(x2, y2, z2) );
}


void metropStep(std::vector< std::tuple<double, double, double> >& vector1,
	            std::vector< std::tuple<double, double, double> >& vector2,
	            std::vector<double>& energyList){
	double x1, y1, z1, x2, y2, z2;
	double xtrial, ytrial, ztrial;
	double rTrial;
	double r_12;
	double r1Current, r2Current;
	double probTrial;
	double newEnergy = 0;
	double metropolisRand = std_rand();
    x1 = std::get<0>(vector1[0]);
    y1 = std::get<1>(vector1[0]);
    z1 = std::get<2>(vector1[0]);
    r1Current = getRvector(x1, y1, z1);
    x2 = std::get<0>(vector2[0]);
    y2 = std::get<1>(vector2[0]);
    z2 = std::get<2>(vector2[0]);
    r2Current = getRvector(x2, y2, z2);

	xtrial = x1 + DELTA*gasdev() ;
	ytrial = y1 + DELTA*gasdev() ;
	ztrial = z1 + DELTA*gasdev() ;
	rTrial = getRvector(xtrial, ytrial, ztrial);
	probTrial = probabiltyWeight(r1Current, r2Current, rTrial, r2Current);
	//std::cout<< "ProbTrial Value = " << probTrial << std::endl;
	if( probTrial > metropolisRand ){ // RANDOM SINGLE ELECTRON MOVE ACCEPTED
		//std::cout << "RANDOM SINGLE ELECTRON MOVE ACCEPTED" << std::endl;
		std::get<0>(vector1[0]) = xtrial;
		std::get<1>(vector1[0]) = ytrial;
		std::get<2>(vector1[0]) = ztrial;
	}
	x1 = std::get<0>(vector1[0]);
    y1 = std::get<1>(vector1[0]);
    z1 = std::get<2>(vector1[0]);
    r1Current = getRvector(x1, y1, z1);
	r_12 = getR12vector(x1, y1, z1, x2, y2, z2);
	////std::cout << "r_12 = " << r_12 << std::endl;
	newEnergy = localEnergy(EXPAND, r1Current, r2Current, r_12 );
	energyList.push_back(newEnergy);
	if(newEnergy > 100){
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
	srand(492832);
	std::cout << std::setprecision(10);
	clock_t start;
	clock_t end;

	std::cout << "Hamiltonian Element : " << hamiltonianHookium(EXPAND, 1, 0.5, 0.5) << std::endl; 
	/*
	start = clock();
	end = clock();
	double TIME = end-start;
	std::cout << "LOOP TIME = " << TIME/CLOCKS_PER_SEC << std::endl;
	*/
    std::vector<double> localEnergyAccumulator;

	std::vector< std::tuple<double, double, double> > r1 ; 
	std::vector< std::tuple<double, double, double> > r2 ; 
	initialiseWalker(r1, r2);
	std::cout << "[x y z] = [" << std::get<0>(r1[0]) <<" "<< std::get<1>(r1[0]) <<" "<< std::get<2>(r1[0]) << "]" << std::endl;


    std::ofstream energyAccum;
    energyAccum.open("ENERGIES_VMC");
    for(int i = 0; i<100000; i++){
    	metropStep(r1, r2, localEnergyAccumulator);
    	metropStep(r2, r1, localEnergyAccumulator);
    }
	std::cout << "[x y z] = [" << std::get<0>(r1[0]) <<" "<< std::get<1>(r1[0]) <<" "<< std::get<2>(r1[0]) << "]" << std::endl;
	 
	for(int i = 0; i < localEnergyAccumulator.size(); i++){
		energyAccum << localEnergyAccumulator[i] << std::endl;
	}
	return 1;
}