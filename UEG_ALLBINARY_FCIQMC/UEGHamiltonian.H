#include "GLOBAL_ORBITAL.H"

#include <iostream>
#include <stdlib.h>
#include <iostream>
#include <stdlib.h>
#include <algorithm> 
#include <boost/math/special_functions/binomial.hpp>
#include <bitset>  
#include <climits> 
#include <vector>



void Di_H_Dj(const double&, double (&)[ORB_SIZE][3], int&, int&, int&, bool&, double&  );
void Di_H_Di(const double&, const int&, long int&, long int&, double (&)[ORB_SIZE][3], double& );
void excitationAlpha_iaBeta_jb(int&, int&, int&, int&, const int&, const int&, long int&, long int&, double (&)[ORB_SIZE][3], int& );
void excitationSameSpinij_ab(int&, int&, int&, int&, const int&, const int&, long int&, double (&)[ORB_SIZE][3], int& );