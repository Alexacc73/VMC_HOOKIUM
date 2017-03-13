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
#include <boost/math/special_functions/factorials.hpp>

const double SQRT2 = 1.4142135623730950488 ;

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
 \beta = \frac{ \sqrt{2}}{ 2^k \sqrt{ (2k-1)! } (2\pi)^{ \frac{3}{4} } }
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
 f_k (r) = H_{2k-1}( \frac{r}{ \sqrt{2} } ) \frac {\beta _{k}}{r} 
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








double hookiumWF(double&r ){

}


double laplacianHookiumWF(double& r){

}


double hamiltonianHookium(double& Ekinetic, double& Epot, double& exchange){

}







int main(void){
	std::cout << exp(4) << std::endl;
	double phi_N = PHI_Kth(1, 1);
	std::cout << "0th phi =  " << phi_N << std::endl;
	return 1;
}