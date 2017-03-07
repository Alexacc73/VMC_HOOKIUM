#include "UEGHamiltonian.H" 
#include "binaryManip.H"

#include <iostream>
#include <stdlib.h>
#include <algorithm> 

    /*printy stuff:
    std::cout << "s_i --> s_a = " << spin1_i << " --> " << spin1_a << std::endl;
    std::cout << "i --> a = [" << KEsortedKpoints[spin1_i][0] << KEsortedKpoints[spin1_i][1] << KEsortedKpoints[spin1_i][2] ;
    std::cout << "] --> [" << KEsortedKpoints[spin1_a][0] << KEsortedKpoints[spin1_a][1] << KEsortedKpoints[spin1_a][2] <<"]"<< std::endl;
    std::cout << "s_j --> s_b = " << spin1_j << " --> " << spin1_b << std::endl;
    std::cout << "j --> b = [" << KEsortedKpoints[spin1_j][0] << KEsortedKpoints[spin1_j][1] << KEsortedKpoints[spin1_j][2] ;
    std::cout << "] --> [" << KEsortedKpoints[spin1_b][0] << KEsortedKpoints[spin1_b][1] << KEsortedKpoints[spin1_b][2] <<"]"<< std::endl;
    */


inline  void INLdecimalToBinary(long int& decimal, std::string& binaryNum)
{
    //std::string binaryNum;
    binaryNum = std::bitset<ORB_SIZE>(decimal).to_string() ;
    //return binaryNum;
}

void Di_H_Dj(const double& cellLength, double (&KEsortedList)[ORB_SIZE][3], int& i, int& a, int& b, bool& ibSpinDifferent, double& RESULT  ){
	/* - - - RESULT: REUTURNS THE DOUBLE INTEGRAL < ij || ab >  =  <D_I | H | D_J> - - - */

  /* ->From the Determinant given by kpoints in detKlist: "i" is the index 
        from which one (alpha spin) electron is excited into orbital index "a".
     ->A second (beta spin) electron is chosen from a filled orbital, j, and is
        excited into the unfilled orbital index "b".
     ->It is assumed that D_J has already been chosen as a coupled DOUBLE excitation of D_I
  such that **CRYSTAL MOMENTUM and SPIN** are conserved.  */
  const double PI = 3.141592653589793;
  double coulomb = 0 ;
  double exchange = 0 ;
  double HElement_ij = 0 ;
  double qDotCoulomb = 0 ;
  double qDotExchange = 0 ;

  qDotCoulomb =  (KEsortedList[i][0] - KEsortedList[a][0])*(KEsortedList[i][0] - KEsortedList[a][0]);
  qDotCoulomb += (KEsortedList[i][1] - KEsortedList[a][1])*(KEsortedList[i][1] - KEsortedList[a][1]);
  qDotCoulomb += (KEsortedList[i][2] - KEsortedList[a][2])*(KEsortedList[i][2] - KEsortedList[a][2]);
  coulomb = 1/(PI*cellLength*qDotCoulomb);
  ////std::cout << "Coulomb part = " << coulomb << std::endl;

  if(ibSpinDifferent == false){
  	qDotExchange =  (KEsortedList[i][0] - KEsortedList[b][0])*(KEsortedList[i][0] - KEsortedList[b][0]);
  	qDotExchange += (KEsortedList[i][1] - KEsortedList[b][1])*(KEsortedList[i][1] - KEsortedList[b][1]);
  	qDotExchange += (KEsortedList[i][2] - KEsortedList[b][2])*(KEsortedList[i][2] - KEsortedList[b][2]);
  	exchange = 1/(PI*cellLength*qDotExchange);
  	////std::cout << "Exchange part = " << exchange << std::endl;
  }
  else{
  	exchange = 0.0;
  }

  HElement_ij = coulomb - exchange;
  RESULT = HElement_ij;
}



void Di_H_Di(const double& cellLength, const int& numElectrons, 
	long int& alphaBin, long int& betaBin, double (&KElist)[ORB_SIZE][3], double& RESULT){

	const double PI = 3.141592653589793;
	double kineticEnergy = 0;
	double exchangeEnergy = 0;
	double totalHamiltonian = 0;
	int alphaEl [numElectrons/2] ;
	int alphaAdd = 0;
	int betaEl [numElectrons/2] ;
	int betaAdd = 0;

    std::string ALPHA(ORB_SIZE, ' ');
    std::string BETA(ORB_SIZE, ' ');
    INLdecimalToBinary(alphaBin, ALPHA);
    INLdecimalToBinary(betaBin, BETA);


	for(int i = 0; i<ORB_SIZE; i++){
		if( ALPHA[i] == '1' ){
			alphaEl[alphaAdd] = ORB_SIZE-1-i ;
			alphaAdd += 1;
		}
		if( BETA[i] == '1' ){
			betaEl[betaAdd] = ORB_SIZE-1-i ;
			betaAdd += 1;
		}
	}

	for(int j = 0; j<(numElectrons/2); j++){
		kineticEnergy += KElist[alphaEl[j]][0]*KElist[alphaEl[j]][0] ;
		kineticEnergy += KElist[alphaEl[j]][1]*KElist[alphaEl[j]][1] ;
		kineticEnergy += KElist[alphaEl[j]][2]*KElist[alphaEl[j]][2] ;

		kineticEnergy += KElist[betaEl[j]][0]*KElist[betaEl[j]][0] ;
		kineticEnergy += KElist[betaEl[j]][1]*KElist[betaEl[j]][1] ;
		kineticEnergy += KElist[betaEl[j]][2]*KElist[betaEl[j]][2] ;
	}
	kineticEnergy *= 2*(PI/cellLength)*(PI/cellLength) ;

    double tempExchangeAlpha = 0;
    double tempExchangeBeta = 0;
	for(int k = 0; k<(numElectrons/2); k++){
		for(int m = k+1; m<(numElectrons/2); m++){
          
			tempExchangeAlpha = (KElist[alphaEl[k]][0] - KElist[alphaEl[m]][0])*(KElist[alphaEl[k]][0] - KElist[alphaEl[m]][0]) ;
			tempExchangeAlpha += (KElist[alphaEl[k]][1] - KElist[alphaEl[m]][1])*(KElist[alphaEl[k]][1] - KElist[alphaEl[m]][1]) ;
			tempExchangeAlpha += (KElist[alphaEl[k]][2] - KElist[alphaEl[m]][2])*(KElist[alphaEl[k]][2] - KElist[alphaEl[m]][2]) ;
			exchangeEnergy += (1/tempExchangeAlpha) ;

			tempExchangeBeta = (KElist[betaEl[k]][0] - KElist[betaEl[m]][0])*(KElist[betaEl[k]][0] - KElist[betaEl[m]][0]) ;
			tempExchangeBeta += (KElist[betaEl[k]][1] - KElist[betaEl[m]][1])*(KElist[betaEl[k]][1] - KElist[betaEl[m]][1]) ;
			tempExchangeBeta += (KElist[betaEl[k]][2] - KElist[betaEl[m]][2])*(KElist[betaEl[k]][2] - KElist[betaEl[m]][2]) ;
			exchangeEnergy += (1/tempExchangeBeta) ;
		}
	}
	exchangeEnergy *= (1/(PI*cellLength)) ;
	

	totalHamiltonian = kineticEnergy - exchangeEnergy ;
    RESULT = totalHamiltonian;
}


/* - - - - - - - - - - - - - EXCITATION OPERATOR - - - - - - - - - - - - - */
    /*
    1) Pick a filled alpha obrital *RANDOMLY*, i
    2) choose an UNfilled alpha determinant *RANDOMLY*, a
       ---> record change in K values from i-->a, delk_ia (for each n,m, and l)
    3) choose a filled beta orbital *RANDOMLY*, j
    4) Run through ALL unfilled beta orbitals, and record those which conserve K, given delk_ia (i.e, give -delk_ia ??)
       ---> Of those which conserve K, choose an unfilled one at *RANDOM* 
    5) Hopefully, the excitation ij ---> ab conserves K with 100% gauruntee 
    6) give the new determinant (J, with empty ij orbits nut filled ab orbits) unique binary labels for alpha and beta parts
    7) Calculate the Hamiltonian between determinants I and new J, < D_I | H | D_J >
       ---> given {i,j,a,b} as inputs, with acces to ordered KPOINTS, obviously. This is < D_I | H | D_J > is simply < ij || ab > 
    8) use < D_I | H | D_J > to calculate spawning prob! - algorithm continues!!! 
    */


void excitationAlpha_iaBeta_jb(int& alpha_i, int& alpha_a, int& beta_j, int& beta_b, const int& numElectrons, const int& numOrbitals,
long int& alphaDetBin, long int& betaDetBin, double (&KEsortedList)[ORB_SIZE][3], int& sign ){
	/* - - - EXCITATION OPERATOR - - - 
    1) Pick a filled alpha obrital *RANDOMLY*, i
    2) choose an UNfilled alpha determinant *RANDOMLY*, a
       ---> record change in K values from i-->a, delk_ia (for each n,m, and l)
    3) choose a filled beta orbital *RANDOMLY*, j
    4) Run through ALL unfilled beta orbitals, and record those which conserve K, given delk_ia (i.e, give -delk_ia ??)
       ---> Of those which conserve K, choose an unfilled one at *RANDOM* 
    5) Hopefully, the excitation ij ---> ab conserves K with 100% gauruntee */
    int filledAlphaOrbits [numElectrons/2] ;
    int UNfilledAlphaOrbits [numOrbitals - numElectrons/2] ;
    int add_fill = numElectrons/2 -1 ;
    int add_unfill = numOrbitals - numElectrons/2 -1 ;

    std::string alphaDetString(ORB_SIZE, ' ');
    std::string betaDetString(ORB_SIZE, ' ');

    INLdecimalToBinary(alphaDetBin, alphaDetString);
    INLdecimalToBinary(betaDetBin, betaDetString);
    
    for(int index = 0; index < numOrbitals; index++){ 
    	int korbit = numOrbitals-1-index;
    	if(alphaDetString[index] == '1'){
    		filledAlphaOrbits[add_fill] = korbit;
    		add_fill -= 1;
    	}
    	else{
    		UNfilledAlphaOrbits[add_unfill] = korbit;
    		add_unfill -= 1;
    	}
    }
    /* 1) Pick a filled alpha obrital *RANDOMLY*, i */
    int index_i = rand()%(numElectrons/2) ; // This pickes a radnom num between 0 and 6, (for 14 electrons, 7 alpha electrons)
    alpha_i = filledAlphaOrbits[index_i] ; // This is the obrital, i, out of which we will excite an electron
    /*2) choose an UNfilled alpha determinant *RANDOMLY*, a*/
    int index_a = rand()% (numOrbitals - numElectrons/2) ;
    alpha_a = UNfilledAlphaOrbits[index_a] ;
    int delta_K [3] = {0};
    /*   ---> record change in K values from i-->a (k_i - k_a), delk_ia (for each n,m, and l)*/
    delta_K[0] = KEsortedList[alpha_i][0] - KEsortedList[alpha_a][0] ;
    delta_K[1] = KEsortedList[alpha_i][1] - KEsortedList[alpha_a][1] ;
    delta_K[2] = KEsortedList[alpha_i][2] - KEsortedList[alpha_a][2] ;
    /*3) choose a filled beta orbital *RANDOMLY*, j */
    int filledBetaOrbits [numElectrons/2] ;
    int UNfilledBetaOrbits [numOrbitals - numElectrons/2] ;
    int beta_fill = numElectrons/2 -1;
    int beta_unfill = numOrbitals - numElectrons/2 -1;

    for(int index = 0; index < numOrbitals; index++){
    	int korbit = numOrbitals-1-index;
    	if(betaDetString[index] == '1'){
    		filledBetaOrbits[beta_fill] = korbit;
    		beta_fill -= 1;
    	}
    	else{
    		UNfilledBetaOrbits[beta_unfill] = korbit;
    		beta_unfill -= 1;
    	}
    }
    /*   ---> Find K of beta_j orbital */
    int index_j = rand()%(numElectrons/2) ;
    beta_j = filledBetaOrbits[index_j] ;
    int beta_j_k [3];
    beta_j_k[0] = KEsortedList[beta_j][0] ;
    beta_j_k[1] = KEsortedList[beta_j][1] ;
    beta_j_k[2] = KEsortedList[beta_j][2] ;

    int b = 0;
    bool foundBeta_b = false ;
    beta_b = -1;
    while( (b<(numOrbitals - numElectrons/2)) && (foundBeta_b == false)  ){
    	int unfilled_index = UNfilledBetaOrbits[b];
    	int Kj_n = beta_j_k[0] - KEsortedList[unfilled_index][0];
    	int Kj_m = beta_j_k[1] - KEsortedList[unfilled_index][1];
    	int Kj_l = beta_j_k[2] - KEsortedList[unfilled_index][2];
    	if( (Kj_n == -delta_K[0]) && (Kj_m == -delta_K[1]) && (Kj_l == -delta_K[2]) ){
    		// std::cout << "KConserving beta orbit, b = " << unfilled_index << std::endl;
    		beta_b = unfilled_index;
    		foundBeta_b = true;
    	}
    	b += 1 ;
    }
    if(foundBeta_b==false){
    	////std::cout << "E1: Did not find a K conserving Beta_b orbital" << std::endl;
    }

    if(foundBeta_b==true){
    	int alphai_excite;
    	int betaj_excite;
    	int temp_excite_alpha[numElectrons/2];
    	int temp_excite_beta[numElectrons/2];
    	for(int i = 0; i<numElectrons/2; i++){
    		temp_excite_alpha[i] = filledAlphaOrbits[i] ;
    		temp_excite_beta[i] = filledBetaOrbits[i] ;
    	}
    	temp_excite_alpha[index_i] = alpha_a;
    	temp_excite_beta[index_j] = beta_b;
    	std::sort(temp_excite_alpha, temp_excite_alpha + numElectrons/2) ;
    	std::sort(temp_excite_beta, temp_excite_beta + numElectrons/2) ;
    	for(int j = 0; j<numElectrons/2; j++){
    		if(temp_excite_alpha[j] == alpha_a){
    			alphai_excite = j;
    		}
    		if(temp_excite_beta[j] == beta_b){
    			betaj_excite = j;
    		}
    	}
    	int indexSum = 0;
    	indexSum = index_i + index_j + alphai_excite + betaj_excite;
    	//std::cout << "idx i = " << index_i << std::endl;
    	//std::cout << "idx_j = " << index_j<< std::endl;
    	//std::cout << "alphai_excite = " << alphai_excite << std::endl;
    	//std::cout << "betaj_excite = " << betaj_excite << std::endl;
    	//std::cout << "indexSum = " << indexSum << std::endl;
    	if(indexSum%2 == 0){
    		sign = 1;
    	}
    	else{
    		sign = -1;
    	}
    }
}





void excitationSameSpinij_ab(int& spin1_i, int& spin1_a, int& spin1_j, int& spin1_b, const int& numElectrons, const int& numOrbitals,
 long int& spin1DetBin, double (&KEsortedList)[ORB_SIZE][3], int& sign ){

 	int filledSpinOrbits [numElectrons/2]  ;
    int UNfilledSpinOrbits [numOrbitals - numElectrons/2] ;
    int add_fill = numElectrons/2 - 1 ;
    int add_unfill = numOrbitals - numElectrons/2 - 1 ;

    std::string spin1DetString(ORB_SIZE, ' ');
    INLdecimalToBinary(spin1DetBin, spin1DetString);
    for(int index = 0; index < numOrbitals; index++){
    	int korbit = numOrbitals-1-index;
    	if(spin1DetString[index] == '1'){
    		filledSpinOrbits[add_fill] = korbit;
    		add_fill -= 1;
    	}
    	else{
    		UNfilledSpinOrbits[add_unfill] = korbit;
    		add_unfill -= 1;
    	}
    }

    //std::cout << filledSpinOrbits[0] << " " << filledSpinOrbits[6] << std::endl;

    int index_i = rand()%(numElectrons/2) ; // This pickes a radnom num between 0 and 6, (for 14 electrons, 7 alpha electrons)
    /* Pick index_j to be DIFFERENT to index_i - find k value of spin1_j spinOrbital*/
    int index_j = index_i ;
    while( index_j == index_i ){ 
    	index_j = rand()%(numElectrons/2) ; 
    	//std::cout << "finished yet?" << std::endl;
    }
    if(index_i > index_j){ //Assert that i < j for calculating correct Hij sign purposes
    	int temp_j;
    	temp_j = index_j;
    	index_j = index_i;
    	index_i = temp_j;
    }

    spin1_i = filledSpinOrbits[index_i] ;
    spin1_j = filledSpinOrbits[index_j] ;
    int spin1_j_k [3];

    spin1_j_k[0] = KEsortedList[spin1_j][0] ;
    spin1_j_k[1] = KEsortedList[spin1_j][1] ;
    spin1_j_k[2] = KEsortedList[spin1_j][2] ;
    /*Choose index a at random*/
    int index_a = rand()% (numOrbitals - numElectrons/2 ) ;
    spin1_a = UNfilledSpinOrbits[index_a] ;

    int deltaK_ia [3] = {0};
    /*   ---> record change in K values from i-->a (k_i - k_a), delk_ia (for each n,m, and l)*/
    deltaK_ia[0] = KEsortedList[spin1_i][0] - KEsortedList[spin1_a][0] ;
    deltaK_ia[1] = KEsortedList[spin1_i][1] - KEsortedList[spin1_a][1] ;
    deltaK_ia[2] = KEsortedList[spin1_i][2] - KEsortedList[spin1_a][2] ;

    int b = 0 ;
    int unfilled_index = 0;
    int Kj_n = 0;
    int Kj_m = 0;
    int Kj_l = 0;
    bool foundspin1_b = false ;
    spin1_b = -1 ;
    while( (b<(numOrbitals - numElectrons/2)) && (foundspin1_b == false)  ){
    	unfilled_index = UNfilledSpinOrbits[b];
    	if(unfilled_index != spin1_a){ /*IMPORTANT!! CANNOT HAVE spin1_a == spin1_b*/
    		Kj_n = spin1_j_k[0] - KEsortedList[unfilled_index][0];
    		Kj_m = spin1_j_k[1] - KEsortedList[unfilled_index][1];
    		Kj_l = spin1_j_k[2] - KEsortedList[unfilled_index][2];
    		if( (Kj_n == -deltaK_ia[0]) && (Kj_m == -deltaK_ia[1]) && (Kj_l == -deltaK_ia[2]) ){
    			////std::cout << "KConserving beta orbit, b = " << unfilled_index << std::endl;
    			if(index_a > b){
    				int temp_b;
    				temp_b = b;
    				b = index_a;
    				index_a = temp_b;
    				spin1_a = UNfilledSpinOrbits[index_a] ;
    			}
    			spin1_b = UNfilledSpinOrbits[b] ;
    			foundspin1_b = true;

    		}
    	}
    	b += 1 ;
    }
    if(foundspin1_b==false){
    	////std::cout << "E2: Did not find a K conserving Beta_b orbital" << std::endl;
    }

    if(foundspin1_b==true){ /*This part is necessary to determine the sign of the Hij Matrix element!!*/
    	int temp_filled [numElectrons/2] ;
    	int temp_filled_excite [numElectrons/2] ;
    	int index_i_excite;
    	int index_j_excite;
    	for(int i = 0; i<numElectrons/2; i++){
    		temp_filled_excite[i] = filledSpinOrbits[i];
    	}
    	temp_filled_excite[index_i] = spin1_a;
    	temp_filled_excite[index_j] = spin1_b;
    	std::sort(temp_filled_excite, temp_filled_excite + numElectrons/2 );
    	for(int j = 0; j<numElectrons/2; j++){
    		if(temp_filled_excite[j] == spin1_a){
    			index_i_excite = j ;
    		}
    		if(temp_filled_excite[j] == spin1_b){
    			index_j_excite = j ;
    		}
    	}
    	int indexSum = index_i + index_j + index_i_excite + index_j_excite;
    	if(indexSum%2 == 0){// is ODD!
    		sign = 1;
    	}
    	else{
    		sign = -1;
    	}
    }
}
