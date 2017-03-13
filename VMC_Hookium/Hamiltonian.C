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


/**
* The nth Hermite polynomial:
* Simply defined via a recursion relation:
* \f$  H_{n+1}(r) = 2xH_n (r) - 2nH_{n-1}(r)  \f$
*/
double hermite_Nth(int N, double r){
	double H1 = 1;
	double H2 = 2*r;
	double H_n;        // H_n
	double H_nm1;      // H_n_minus1
	double H_np1;      // H_n_plus1 
	if(N==1){
		return H1;
	}
	if(N==2){
		return H2;
	}
	H_nm1 = H1;
	H_n = H2;    // Currently only computed up to n = 2;
	for(int n = 2; n < N; n++ ){
		H_np1 = 2*r*H_n - 2*n*H_nm1 ;
		H_nm1 = H_n;
		H_n = H_np1;
	}
	return H_np1;
}

/**
* The first derivative of the nth Hermite polynomial, easily defined from the nth Hermite:
* \f$ H'_n (r) = 2nH_{n-1}(r)  \f$
*/
double hermite_diff1_Nth(int N, double r){
	double result;
	result = 2*N*hermite_Nth(N-1, r );
	return result;
}


/**
* The second derivative of the nth Hermite polynomial, again, easily defined from the nth Hermite:
* \f$  H''_n (r) = \frac{\partial}{\partial r} 2nH_{n-1}(r)  = 4n(n-1)H_{n-2} (r)  \f$
*/
double hermite_diff2_Nth(int N, double r){
	double result;
	result = 4*N*(N-1)*hermite_Nth(N-2, r);
	return result;
}


/**
* This is simply a constant coefficient for f(r), or "FofR_Kth":
* \f$ \beta = \frac{ \sqrt{2}}{ 2^k \sqrt{ (2k-1)! } (2\pi)^{ \frac{3}{4} } } \f$
*/
double beta_Kth(int K){
	double beta;
	double denom;
	double Kfloat = K;
	double factorial = boost::math::factorial<double>(2*Kfloat - 1);

	beta = pow(2, 0.5);
	denom = pow(2, Kfloat) * pow( factorial ,0.5) * pow(2*M_PI, (3.0/4.0) );
	beta *= 1/denom;
	return beta;


}

double FofR_Kth(){

}

double FofR_diff1_Kth(){

}

double FofR_diff2_Kth(){

}

double THETA_Kth(){

}

double THETA_laplace_Kth(){

}








double hookiumWF(double&r ){

}


double laplacianHookiumWF(double& r){

}


double hamiltonianHookium(double& Ekinetic, double& Epot, double& exchange){

}







int main(void){
	std::cout << exp(4) << std::endl;

}