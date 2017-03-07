/*! \mainpage FCIQMC HEAD PAGE
 * \section intro_sec Introduction
 *  ---> INCLUDED IN HEADER < GLOBAL_ORBITAL.H > is " const int ORB_SIZE = 125 " <--- 
    ---> NB - - NUMBER OF ORBITALS MUST BE KNOWN A PRIORI FOR "binaryManip.C* AND for KEsortedKpoints[size][3]

*/


#include "binaryManip.H"    /* Inludes  <bitset>  <climits>  <boost/math/special_functions/binomial.hpp>  <vector> */
#include "planeWaves.H"     /* Includes <iostream>  <set>  <math.h>  <stdlib.h>  <algorithm> */
#include "UEGHamiltonian.H" /*Includes initialisation: const int ORB_SIZE = 125;*/


#include <fstream>
#include <time.h>
#include <set>
#include <map>






/*
*Parameters for initial setup of UEG 
*/

/** A number, nobody knows what it means */
const double PI = 3.141592653589793;

/** rs controls density of the Electron Gas. 
* "rs" is the radius of the sphere whose volume is that of the cell divided by the number of electrons*/
const double rs = 1.0; 

/** Total number of electrons -> half allocated Alpha spin, the other half allocated Beta Spin.*/
const double numElectrons = 14; 
const int INTelectrons = numElectrons ;

/** Kc_CUTTOFF is the kinetic energy cutoff for the plane wave basis orbitals.
E.g, a cutoff of "2" will allow the orbital [4 0 0] but not [5 0 0]. Set cutoff = 2.4 for 57 Orbitals (114 Spin Orbitals) */
const double Kc_CUTTOFF = 2 ; 



/*
 *-----> Parameters for Main FCIQMC Algorithm <----- 
 */

/** delt is the Imaginary timestep for the propogation of the "walker" population */
const double delt = 0.003 ;

/** Zeta is a damping parameter which controls the agressiveness of the "shift" in the variable shift mode of the algorithm */
const double zeta = 0.04 ;

/** AShift controls how frequently the shift is changed in response to the population in the variable shift mode (AShift = 1 means every step) */
const int AShift = 2 ;

/** Number of steps after which to terminate the algorithm*/
const int numSteps = 500000;

/** After "walker critical" walkers have been spawned after a complete cycle (post annihilation) the variable shift mode is turned on */
const int walkerCritical = 500000;

/** initRefWalkers is the number of wlakers which are initially placed on the reference (i.e Hartree Fock) determinant to begin the spawning */
int initRefWalkers = 50;


long int pow2Array [ORB_SIZE];


inline constexpr std::uint64_t INLpow2 (std::uint64_t i)
{
    return std::uint64_t(1) << i;
}

inline const int INLgetPositionInList(std::pair<long int, long int>& uniqueDet, std::vector<long int>& alphaDetList, std::vector<long int>& betaDetList)
{
    int DETLIST_size = alphaDetList.size();
    int pos ;
    int j = 0;
    bool not_found = true ;

    long int alphaBIN = uniqueDet.first ;
    //binaryToDecimal(ALPHA, alphaBIN);
    long int betaBIN = uniqueDet.second ;
    //binaryToDecimal(BETA, betaBIN);
    while( (j<DETLIST_size) && (not_found == true)  ){
        if( (alphaBIN == alphaDetList[j]) && (betaBIN == betaDetList[j]) ){
            pos = j ;
            not_found = false ;
        }
        j++ ;
    }
    if(not_found == true){
        //pos = -1;
        std::cout << "OOPS POTISION FIND WENT WRONG!!!!!" << std::endl;
    }
    return pos ;
}




void SPAWN( const double cellLength, 
            std::vector<int>& trueWalkerList, 
            std::vector<int>& posWalkerList, 
            std::vector<int>& negWalkerList,
            double (&KEsortedList)[ORB_SIZE][3], 
            std::set< std::pair<long int, long int> >& uniqueDetSet, 
            std::pair< std::set< std::pair<long int, 
            long int> >::iterator , bool >& result, 
            std::vector<long int>& alphaDets, 
            std::vector<long int>& betaDets ){

    int walkerNum = 0;
    int intpSpawn = 0;
    int index_i = 0, index_j = 0;
    int index_a = 0, index_b = 0;
    int randChooseExcite = 0;
    int SIGN = 1;
    int position = 0;
    double pSpawn = 0.0; /*Probability of selecting certain determinant*/
    double remainderpSpawn = 0.0;
    double pGen = 0.0;
    double HijElement = 0.0;
    double metropolisRand = 0.0; /*Randomly generated number upon [0,1] to compare with pSelect (metropolis criteria)*/
    long int alphaDet;
    long int betaDet;
    //std::string alphaDetSTR;
    //std::string betaDetSTR;
    //std::string combinedAlphaBeta(ORB_SIZE*2 + 1, ' ');

    std::pair<long int, long int> iter ;

    bool ib_spinDifferent;

    int numDets = uniqueDetSet.size() ;
    for(int det = 0; det < numDets; det++){
        walkerNum = trueWalkerList[det] ;
        if(walkerNum != 0){

            for(int num = 0; num < walkerNum; num++){
                alphaDet = alphaDets[det] ;
                betaDet = betaDets[det] ;
                index_i = 0 ;
                index_j = 0 ;
                index_a = 0 ;
                index_b = 0 ;
                randChooseExcite = rand()%4 ;
                if(randChooseExcite == 0){ /* i-->a AND j-->b are both ALPHA spin*/
                    excitationSameSpinij_ab(index_i, index_a, index_j, index_b, INTelectrons, ORB_SIZE, alphaDet, KEsortedList, SIGN);
                    bool ib_spinDiffAlpha = false;
                    pGen = (1.0/( (numElectrons/2.0) *((numElectrons/2.0)-1.0)) ) * ( 2/(ORB_SIZE - (numElectrons/2)) ) ;
                    pGen = pGen/4.0 ;
                }
                if(randChooseExcite == 1){/* i-->a AND j-->b are both BETA spin*/
                    excitationSameSpinij_ab(index_i, index_a, index_j, index_b, INTelectrons, ORB_SIZE, betaDet, KEsortedList, SIGN);
                    bool ib_spinDiffAlpha = false;
                    pGen = (1.0/( (numElectrons/2.0) *((numElectrons/2.0)-1.0)) ) * ( 2/(ORB_SIZE - (numElectrons/2)) ) ;
                    pGen = pGen/4.0 ;
                }
                if( (randChooseExcite == 2)  || (randChooseExcite == 3) ){ /* i-->a is ALPHA,  j-->b is BETA*/
                    excitationAlpha_iaBeta_jb(index_i, index_a, index_j, index_b, INTelectrons, ORB_SIZE, alphaDet, betaDet, KEsortedList, SIGN);
                    bool ib_spinDifferent = true;
                    pGen = (1.0/((numElectrons/2.0)*(numElectrons/2.0)) ) * ( 1/(ORB_SIZE - (numElectrons/2)) ) ;
                    pGen = pGen/2.0 ;
                }
    
                if(index_b != -1){ /*Finding a k consevring orbital was SUCCESSFUL*/
                    /*1) Update newly excited determinants*/
                    // decimalToBinary(alphaDet, alphaDetSTR);
                    // decimalToBinary(betaDet, betaDetSTR);

                    if(randChooseExcite==0){
                        alphaDet -= pow2Array[  index_i ] ; //= '0' ;
                        alphaDet -= pow2Array[  index_j ] ;//= '0' ;
                        alphaDet += pow2Array[  index_a ] ;//= '1' ;
                        alphaDet += pow2Array[  index_b ] ;//= '1' ;

                        //std::cout << "New Alpha Det = " << alphaDet << std::endl;
                    }
                    if(randChooseExcite==1){
                        betaDet -= pow2Array[ index_i ] ;// = '0' ;
                        betaDet -= pow2Array[ index_j ] ;// = '0' ;
                        betaDet += pow2Array[ index_a ] ;// = '1' ;
                        betaDet += pow2Array[ index_b ] ;// = '1' ;
                    }
                    if( (randChooseExcite == 2)  || (randChooseExcite == 3) ){
                        alphaDet -= pow2Array[ index_i] ; // = '0' ;
                        alphaDet += pow2Array[ index_a] ; // = '1' ;
                        betaDet  -= pow2Array[ index_j] ; //= '0' ;
                        betaDet  += pow2Array[ index_b] ; //= '1' ;
                    }
                   
                    HijElement = 0;
                    Di_H_Dj(cellLength , KEsortedList, index_i, index_a, index_b, ib_spinDifferent, HijElement) ;
                    HijElement *= SIGN;
                    pSpawn = (delt* fabs(HijElement) ) / pGen ;
                    intpSpawn = floor(pSpawn) ;
                    remainderpSpawn = pSpawn - intpSpawn ;
                    metropolisRand = ((double) rand() / RAND_MAX) ; //(fastrand() % 100000 )/100000 ;

                    //std::cout << "pSpawn = " << pSpawn << std::endl;
                    //std::cout << "Excitation scheme = " << randChooseExcite << std::endl;
                    //std::cout << "pGen = " << pGen << std::endl;
                    //std::cout << "remainderpSpawn = " << remainderpSpawn << std::endl;
                    //std::cout << "metropolisRand = " << metropolisRand << std::endl;

                    /* SPAWN: YES or NO*/
                    if(remainderpSpawn >= metropolisRand){/* --> SPAWING EVENT SUCCESSFUL <--*/
                    //std::cout << "remainder > metropolis " << std::endl;

                        //concatStrings(alphaDetSTR, betaDetSTR, combinedAlphaBeta) ;
                        result = uniqueDetSet.insert( std::make_pair(alphaDet, betaDet) ) ;

                        if(result.second == true){/* New Determinant Found */
                        /*All 3 walker lists *MUST* each increase their size by 1,
                        hence the .push_back(0) for each .puch_back(1) 1*/
                            trueWalkerList.push_back(0) ;
                            // binaryToDecimal(alphaDetSTR, alphaDet);
                            // binaryToDecimal(betaDetSTR, betaDet);
                            alphaDets.push_back(alphaDet) ;
                            betaDets.push_back(betaDet) ;
                            if( (HijElement < 0) && (walkerNum > 0) ){
                                posWalkerList.push_back(1) ;
                                negWalkerList.push_back(0) ;
                            }
                            if( (HijElement > 0) && (walkerNum < 0) ){
                                posWalkerList.push_back(1) ;
                                negWalkerList.push_back(0) ;
                            }
                            if( (HijElement < 0) && (walkerNum < 0) ){
                                negWalkerList.push_back(1) ;
                                posWalkerList.push_back(0) ;
                            }
                            if( (HijElement > 0) && (walkerNum > 0) ){
                                negWalkerList.push_back(1) ;
                                posWalkerList.push_back(0) ;
                            }
                        }
                        if(result.second == false){/* Determinant ALREADY EXISTS */
                            iter = *result.first ;
                            /*Find position of determinant in list*/
                            position = INLgetPositionInList(iter, alphaDets, betaDets) ;
                            if( (HijElement < 0) && (walkerNum > 0) ){
                                posWalkerList[position] += 1 ;
                            }
                            if( (HijElement > 0) && (walkerNum < 0) ){
                                posWalkerList[position] += 1 ;
                            }
                            if( (HijElement < 0) && (walkerNum < 0) ){
                                negWalkerList[position] += 1 ;
                            }
                            if( (HijElement > 0) && (walkerNum > 0) ){
                                negWalkerList[position] += 1 ;
                            }
                        }
                    }

                    /*SPAWN AGAIN: YES or NO*/
                    if(intpSpawn >= 1){
                        //concatStrings(alphaDetSTR, betaDetSTR, combinedAlphaBeta) ;
                        result = uniqueDetSet.insert( std::make_pair(alphaDet, betaDet) ) ;
                        if(result.second == true){/* New Determinant Found */
                            trueWalkerList.push_back(0) ;
                            //binaryToDecimal(alphaDetSTR, alphaDet);
                            //binaryToDecimal(betaDetSTR, betaDet);
                            alphaDets.push_back(alphaDet) ;
                            betaDets.push_back(betaDet) ;
                            if( (HijElement < 0) && (walkerNum > 0) ){
                                posWalkerList.push_back(intpSpawn) ;
                                negWalkerList.push_back(0) ;
                            }
                            if( (HijElement > 0) && (walkerNum < 0) ){
                                posWalkerList.push_back(intpSpawn) ;
                                negWalkerList.push_back(0) ;
                            }
                            if( (HijElement < 0) && (walkerNum < 0) ){
                                negWalkerList.push_back(intpSpawn) ;
                                posWalkerList.push_back(0) ;
                            }
                            if( (HijElement > 0) && (walkerNum > 0) ){
                                negWalkerList.push_back(intpSpawn) ;
                                posWalkerList.push_back(0) ;
                            }
                        }/*END if result.second==true*/
                        if(result.second == false){/* Determinant ALREADY EXISTS */
                            iter = *result.first ; 
                            /*Find position of determinant in list*/
                            position = INLgetPositionInList(iter, alphaDets, betaDets) ;
                            if( (HijElement < 0) && (walkerNum > 0) ){
                                posWalkerList[position] += intpSpawn ;
                            }
                            if( (HijElement > 0) && (walkerNum < 0) ){
                                posWalkerList[position] += intpSpawn ;
                            }
                            if( (HijElement < 0) && (walkerNum < 0) ){
                                negWalkerList[position] += intpSpawn ;
                            }
                            if( (HijElement > 0) && (walkerNum > 0) ){
                                negWalkerList[position] += intpSpawn ;
                            }
                        }
                    }


                }/*End of IF(index_b != -1)*/
                else{/* COULD NOT find a k consevring orbital - SET Hij = 0 and carry on as usual!*/
                    HijElement = 0 ;
                    pSpawn = 0 ;
                    /*Do Nothing - we know with certainty that pSpawn = 0 since HijElement = 0*/
                }


            }/*End of inner walker FOR loop*/
        }/*End of IF loop*/
    } /*End of main FOR loop*/
}/* end of function SPAWN */



void DEATH_CLONE(const double cellLength, 
                std::vector<int>& trueWalkerList, 
                double (&KEsortedList)[ORB_SIZE][3], 
                double& SHIFT,
                std::vector<long int>& alphaDets, 
                std::vector<long int>& betaDets, 
                const double& HFEnergy){
    int walkerNum = 0;
    int INTnDeath = 0;
    int walkersToDie = 0;
    int walkersToClone = 0;
    int numDets = 0;
    double pDeath = 0.0;
    double nDeath = 0.0;
    double RandRound = 0.0;
    double HiiElement = 0;
    double KiiElement = 0;
    long int alphaDet;
    long int betaDet;
    //std::string alphaDet(ORB_SIZE, ' ') ;
    //std::string betaDet(ORB_SIZE, ' ') ;

    numDets = trueWalkerList.size() ;
    for(int det = 0; det<numDets; det++){
        walkerNum = trueWalkerList[det] ;
        if(walkerNum != 0){
            alphaDet = alphaDets[det];
            betaDet = betaDets[det];

            HiiElement = 0;
            Di_H_Di(cellLength, INTelectrons, alphaDet, betaDet, KEsortedList, HiiElement);
            KiiElement = HiiElement - HFEnergy ;

            pDeath = delt*(KiiElement - SHIFT) ;
            nDeath = fabs(pDeath*walkerNum) ;
            RandRound =  ((double) rand() / RAND_MAX) ; //(rand() % 1000000 )/1000000.0 ;
            INTnDeath = floor(nDeath + RandRound); 

            //std::cout << "Kii diagonal matrix = " << KiiElement << std::endl;
            //std::cout << "nDeath = " << nDeath << std::endl;
            //std::cout << " nDeath + RandRound = " << nDeath + RandRound << std::endl;
            //std::cout << "INTnDeath = " << INTnDeath << std::endl;

            if(pDeath > 0){
                walkersToDie = copysign(INTnDeath, walkerNum);
                trueWalkerList[det] -= walkersToDie;
            }
            if(pDeath < 0){
                walkersToClone = copysign(INTnDeath, walkerNum);
                trueWalkerList[det] += walkersToClone;
            }
        }
    }/* End of FOR loop */
}


void ANNIHILATION(int& step, 
                std::vector<int>& trueWalkerList, 
                std::vector<int>& posWalkerList, 
                std::vector<int>& negWalkerList,
                std::set< std::pair<long int, 
                long int> >& uniqueDeterminantSet, 
                std::vector<long int>& alphaDetsBinary, 
                std::vector<long int>& betaDetsBinary){

    int numDets = 0;
    bool prune = false;

    if(step%8==0){
        prune = true;
    }
    numDets = trueWalkerList.size();

    for(int det = 0; det<numDets; det++){
        trueWalkerList[det] += (posWalkerList[det] - negWalkerList[det]) ;
        /* Take this opportunity to clear temporary pos and neg walker lists */
        posWalkerList[det] = 0;
        negWalkerList[det] = 0;
        
    }

    if(prune==true){
        

        std::vector<int> TWL_copy;
        std::vector<long int> alphaDets_copy;
        std::vector<long int> betaDets_copy;
    
        TWL_copy.reserve(trueWalkerList.size());
        alphaDets_copy.reserve(alphaDetsBinary.size());
        betaDets_copy.reserve(betaDetsBinary.size());

        copy(trueWalkerList.begin(), trueWalkerList.end(), back_inserter(TWL_copy));
        copy(alphaDetsBinary.begin(), alphaDetsBinary.end(), back_inserter(alphaDets_copy));
        copy(betaDetsBinary.begin(), betaDetsBinary.end(), back_inserter(betaDets_copy));

        trueWalkerList.clear();
        alphaDetsBinary.clear();
        betaDetsBinary.clear();

        for(int el = 0; el<TWL_copy.size(); el++){
            if(TWL_copy[el] != 0){
                trueWalkerList.push_back(TWL_copy[el]);
                alphaDetsBinary.push_back(alphaDets_copy[el]);
                betaDetsBinary.push_back(betaDets_copy[el]);
            }
        }

        long int alphaBin;
        long int betaBin;
        // std::string alphaDet(ORB_SIZE, ' ');
        // std::string betaDet(ORB_SIZE, ' ');
        // std::string combinedAlphaBeta(ORB_SIZE*2 + 1, ' ');

        posWalkerList.resize(trueWalkerList.size()) ;
        negWalkerList.resize(trueWalkerList.size()) ;

        uniqueDeterminantSet.clear();
        for(int det = 0; det < alphaDetsBinary.size(); det++){
            alphaBin = alphaDetsBinary[det];
            betaBin = betaDetsBinary[det];
            //decimalToBinary(alphaBin, alphaDet);
            //decimalToBinary(betaBin, betaDet);
            //concatStrings(alphaDet, betaDet, combinedAlphaBeta);
            uniqueDeterminantSet.insert( std::make_pair(alphaBin, betaBin) );
        }
    }/*End Prune*/

}


/* - - - - - AUXILLIARY FUNCTIONS - - - - -*/
int totalWalkerNumber( std::vector<int>& trueWalkerList){
  int count = 0;
  int numDets = trueWalkerList.size() ;
  for(int det = 0; det < numDets; det++){
    //cout << "WALKER VAL = " << trueWalkerList[det] << endl;
    count += abs(trueWalkerList[det]);
  }
  return count;
}


/** variableShift function controls the shift in response to the population, and relates exactly to the equation:
 \f[
 S(\tau) = S(\tau - A\delta \tau) - \frac{\zeta}{A \delta \tau} ln\frac{N_w (\tau)}{N_w (\tau - A\delta \tau)} 
 \f]
 */
double variableShift( const double& delt, 
                      const int& AShift, 
                      int& intTimestep, 
                      const double& zeta, 
                      std::vector<double>& totalWalkerTracker, 
                      std::vector<double>& shiftTracker){

  double walkNumT = totalWalkerTracker[intTimestep] ;
  double walkNumTminusA = totalWalkerTracker[intTimestep - AShift] ;
  double shiftTminusA = shiftTracker[intTimestep - AShift] ;
  double newShift = 0 ;
  double fAshift = AShift;
  double logVal = log(walkNumT/walkNumTminusA) ;
  newShift = shiftTminusA - ( (zeta/(fAshift*delt)) * logVal ) ;
  return newShift ;
}

/** The Projector Energy gives an independent estimation of the energy from the shift, according to equation:
\f[
 E_{proj} =  \sum_{j \neq 0} \Big \langle D_j | H | D_0 \Big \rangle \frac{N_j (\tau)}{N_0 (\tau)} 
\f]
*/

double projectorEnergy(){
    return 1;
}





int main(void){
    srand(492831);

	
	const double cellVolume = (numElectrons*PI*4.0*rs*rs*rs)/3.0;
	const double cellLength = getCellLength(cellVolume); /* Return cell length in units of rs */


    for(int i = 0; i<ORB_SIZE; i++){
        pow2Array[i] = INLpow2(i);
    }

	

	int spinOrbitals;
    double kpoints[500][3];
    createPlaneWaveOrbitals(kpoints, Kc_CUTTOFF, spinOrbitals);

    std::cout << "Cell length = " << cellLength << '\n' << std::endl;
    std::cout << "Number spin orbitals = " << spinOrbitals << std::endl;
    std::cout << "Number of orbitals = " << ORB_SIZE  << std::endl;
    std::cout << "Number of Electrons = " << numElectrons << '\n' << std::endl;
    if(spinOrbitals/2 != ORB_SIZE){
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        std::cout << "ERROR: Orbital size is not correct: Alter in file < GLOBAL_ORBITAL.H > " << std::endl;
        std::cout << "---> SET < ORB_SIZE > TO BE : " << spinOrbitals/2 << std::endl;
        std::cout << "-----------------------------------------------------------------------" << std::endl;
        return 1;
    }

    double KEsortedKpoints[ORB_SIZE][3];
    /* Create a set of UNIQUE orbitals, in ascending Kinetic energy order */
    kPointsEnergySort(kpoints, KEsortedKpoints, ORB_SIZE);
    std::cout << "ORDERED LIST: " << std::endl;
    PRINTORBITALS(ORB_SIZE, KEsortedKpoints);


    /* -----> DEFINE IMPORTANT LISTS <----- */
    std::set< std::pair<long int, long int> > uniqueDeterminantSet ;
    std::pair< std::set< std::pair<long int, long int> >::iterator , bool > result;

    std::vector<long int> alphaDetsBinary;
    std::vector<long int> betaDetsBinary;

    std::vector<int> walkerList;
    std::vector<int> posChildList;
    std::vector<int> negChildList;

    std::vector<double> walkerNUMTracker;
    std::vector<double> SHIFTTracker;

    walkerList.push_back(initRefWalkers); /* ---> Put X walkers (e.g, X = 10) on HF determinent */
    posChildList.push_back(0); /*Makes sure that dhild lists are always same size as walker list*/
    negChildList.push_back(0);

    // ---> Fill with HF orbitals - lowest KE occupied!!
    long int zero = 0;
    std::string HF_STRING;
    decimalToBinary(zero, HF_STRING) ; 
    for(int el = 0; el < INTelectrons/2; el++){
        HF_STRING[ORB_SIZE-1-el] = '1';
    }
    
    //std::string combinedHF(ORB_SIZE*2 + 1, ' ');
    //concatStrings(HF_STRING, HF_STRING, combinedHF ) ;
    long int HF_binary;
    binaryToDecimal(HF_STRING, HF_binary);
    result = uniqueDeterminantSet.insert( std::make_pair(HF_binary, HF_binary) ) ; 
    
    if(result.second == true){
        std::cout << "INITIAL INPUT SUCCESSFUL" << std::endl;
        std::pair<long int, long int> iter = *result.first;
        std::cout << "HF Alpha + Beta DETERMINANT = " << iter.first << " " << iter.second << std::endl;
        alphaDetsBinary.push_back( HF_binary ) ;
        betaDetsBinary.push_back( HF_binary ) ;
    }
    double EnergyHFnonConst;
    Di_H_Di(cellLength, INTelectrons, HF_binary, HF_binary, KEsortedKpoints, EnergyHFnonConst ) ;
    const double HFEnergy = EnergyHFnonConst;
    std::cout << " < D_I | H | D_I > Element for HF = " << HFEnergy << '\n' << std::endl;

    //std::string large;
    //decimalToBinary(HF_binary - pow2(ORB_SIZE-1-56), large) ;
    //std::cout << large <<  std::endl;





    for(int i = 0; i<1; i++){
        int index_i = 0, index_a = 0, index_j = 0, index_b = 0;
        int sign = 0;
        excitationAlpha_iaBeta_jb(index_i, index_a, index_j, index_b, INTelectrons, ORB_SIZE, HF_binary, HF_binary, KEsortedKpoints, sign) ;
        if(index_b != -1){
            std::cout << "s_i --> s_a = " << index_i << " --> " << index_a << std::endl;
            std::cout << "i --> a = [" << KEsortedKpoints[index_i][0] << KEsortedKpoints[index_i][1] << KEsortedKpoints[index_i][2] ;
            std::cout << "] --> [" << KEsortedKpoints[index_a][0] << KEsortedKpoints[index_a][1] << KEsortedKpoints[index_a][2] <<"]"<< std::endl;
            std::cout << "s_j --> s_b = " << index_j << " --> " << index_b << std::endl;
            std::cout << "j --> b = [" << KEsortedKpoints[index_j][0] << KEsortedKpoints[index_j][1] << KEsortedKpoints[index_j][2] ;
            std::cout << "] --> [" << KEsortedKpoints[index_b][0] << KEsortedKpoints[index_b][1] << KEsortedKpoints[index_b][2] <<"]"<< std::endl;
            std::cout << "SIGN : " << sign << std::endl;
            std::cout << " " << std::endl;
        }
    }





    /* - - - - - - - - M A I N   L O O P - - - - - - - - - - - - */
    int beginAverageStep;
    double SHIFT = 0.1 ;
    double instantAverageShift = 0.0;
    double instantAvergeProjector = 0.0;
    double elapsedTime;  
    double currentProjectorEnergy = 0.0;
    double compareShiftProjector = 0.0;
    double compareMeanEnergies = 0.0;
    double referenceWalkers = 0;
    double currentPopulation = 0;
    double popToRefRatio = 0;

    clock_t start;
    clock_t end;
    bool shiftOn = false;

    std::ofstream shoulderplot;
    std::ofstream shiftPlot;
    shoulderplot.open ("SHOULDER_66S0_500000WC.txt");
    shiftPlot.open("SHIFT_66SO_50000WC.txt");
   
    for(int i = 0; i < numSteps; i++){
        
        currentPopulation = totalWalkerNumber(walkerList) ;
        walkerNUMTracker.push_back( currentPopulation );

        referenceWalkers = walkerList[0] ;
        popToRefRatio = currentPopulation / referenceWalkers ;

        std::cout << "-------------> STEP : " << i << std::endl;
        std::cout << "CURRENT POPULATION = " << walkerNUMTracker[i] << std::endl;
        std::cout << "REFERENCE WALKERS = " << referenceWalkers << std::endl;
        std::cout << "CURRENT NUMBER OF DETERMINANTS = " <<  uniqueDeterminantSet.size() << std::endl;
        std::cout << "CURRENT NUMBER OF ALPHAS = " <<  alphaDetsBinary.size() << std::endl;
        std::cout << "CURRENT NUMBER OF BETAS = " <<  betaDetsBinary.size() << '\n' << std::endl;

        /*TURN SHIFT ON NOW, *THEN* perform spawn-death/clone-annihilation*/

        if( (currentPopulation >= walkerCritical ) && !shiftOn ){
            shiftOn = true;
            beginAverageStep = i*2;
            std::cout << "THIS TEXT SHOULD APPEAR EXACTLY ONCE" << std::endl;
        }
        if(shiftOn && ( (i)%AShift == 0) ){ 
            SHIFT = variableShift(delt, AShift, i, zeta, walkerNUMTracker, SHIFTTracker) ;
        }
        SHIFTTracker.push_back(SHIFT) ;

        //projectorTracker[i] = currentProjectorEnergy;

        shoulderplot << currentPopulation << " " << popToRefRatio << std::endl ;
        shiftPlot << SHIFT << std::endl;


        start = clock();
        SPAWN(cellLength, walkerList, posChildList, negChildList, KEsortedKpoints, uniqueDeterminantSet, result, alphaDetsBinary, betaDetsBinary );
        DEATH_CLONE(cellLength, walkerList, KEsortedKpoints, SHIFT, alphaDetsBinary, betaDetsBinary, HFEnergy);
        ANNIHILATION(i, walkerList, posChildList, negChildList, uniqueDeterminantSet, alphaDetsBinary, betaDetsBinary);
        end = clock();
        elapsedTime = end-start;
        //std::cout << "Triple loop Time = " << elapsedTime << std::endl;
        std::cout<< "Elapsed time for triple loop: " << elapsedTime/CLOCKS_PER_SEC << std::endl;

    }

    for(int i = 0; i< walkerList.size(); i++ ){
        std::cout << "Walker on DET " << i << " = " << walkerList[i] << std::endl;
    }

    shoulderplot.close();

    //for (std::set<std::string>::iterator it = uniqueDeterminantSet.begin(); it != uniqueDeterminantSet.end(); it++) {
    //    std::string element = *it ;
    //    std::cout << element << std::endl;
   // }
  
    

    return 1;
}
















    /*

    int det_index = 0;
    int numDuplicates = 0;
    int numDetsCreated = 1;
    for(int run = 0; run < 1000; run++){

        int alpha_i = 0, alpha_a = 0;
        int beta_j = 0, beta_b = 0;
        std::string alphaDetStr = alphaDetsBinary[det_index] ;
        std::string betaDetstr =  betaDetsBinary[det_index] ;

        ////std::cout << "-------> Begin Run " << run << " <--------" << std::endl;

        excitationAlpha_iaBeta_jb(alpha_i, alpha_a, beta_j, beta_b, INTelectrons, ORB_SIZE, betaDetstr, alphaDetStr,  KEsortedKpoints );
        bool ib_SpinDifferent = true;
        if(beta_b != -1){
            numDetsCreated += 1;

            double hamltnij__ab =  Di_H_Dj(cellLength , KEsortedKpoints, alpha_i, alpha_a, beta_b, ib_SpinDifferent ) ;
            ////std::cout << " < D_I | H | D_J > Element =  < ij || ab > is : " << hamltnij__ab << '\n' << std::endl;
            //Create new determinant into list
            betaDetstr[ORB_SIZE-1-alpha_i] = '0' ;
            betaDetstr[ORB_SIZE-1-alpha_a] = '1' ;
            alphaDetStr[ORB_SIZE-1-beta_j] = '0' ;
            alphaDetStr[ORB_SIZE-1-beta_b] = '1' ;

            std::string combinedAlphaBeta(ORB_SIZE*2 + 1, ' ');
            concatStrings(alphaDetStr, betaDetstr, combinedAlphaBeta) ;
            result = uniqueDeterminantSet.insert( combinedAlphaBeta ) ;
            if(result.second == true){
                // INPUT SUCCESSFUL 
                alphaDetsBinary.push_back(alphaDetStr) ;
                betaDetsBinary.push_back(betaDetstr) ;
                det_index += 1 ;
                walkerList.push_back(1);
            }
            if(result.second == false){
                // INPUT FAILED 
                std::string iter = *result.first;
                int position = getPositionInList( iter, alphaDetsBinary, betaDetsBinary ) ;
                walkerList[position] += 1;
                numDuplicates += 1 ;
            }
        }

        int spin1_i = 0, spin1_a = 0 ;
        int spin1_j = 0, spin1_b = 0 ;
        std::string spinDetStr = alphaDetsBinary[det_index] ;

        excitationSameSpinij_ab(spin1_i, spin1_a, spin1_j, spin1_b, INTelectrons, ORB_SIZE, spinDetStr, KEsortedKpoints ) ;
        bool ib_SpinDiffAlpha = false;
        if(spin1_b != -1){
            numDetsCreated += 1;

            double hamltn_spinij__ab =  Di_H_Dj(cellLength , KEsortedKpoints, spin1_i, spin1_a, spin1_b, ib_SpinDiffAlpha ) ;
            ////std::cout << " < D_I | H | D_J > Element =  < ij || ab > is : " << hamltn_spinij__ab << std::endl;
            spinDetStr[ORB_SIZE-1- spin1_i ] = '0' ;
            spinDetStr[ORB_SIZE-1- spin1_j ] = '0' ;
            spinDetStr[ORB_SIZE-1- spin1_a ] = '1' ;
            spinDetStr[ORB_SIZE-1- spin1_b ] = '1' ;

            std::string combinedStrings(ORB_SIZE*2 + 1, ' ');
            concatStrings(spinDetStr, betaDetstr, combinedStrings) ;
            result = uniqueDeterminantSet.insert( combinedStrings) ;
            if(result.second == true){
                // INPUT SUCCESSFUL 
                alphaDetsBinary.push_back(spinDetStr) ;
                betaDetsBinary.push_back(betaDetstr) ;
                det_index += 1 ;
                walkerList.push_back(1);
            }
            if(result.second == false){
                // INPUT FAILED 
                std::string iter = *result.first;
                int position = getPositionInList( iter, alphaDetsBinary, betaDetsBinary ) ;
                walkerList[position] += 1;
                numDuplicates += 1 ;
            }
        }
    }




    std::cout << "Number of ALPHA determinants is now = " << alphaDetsBinary.size() << std::endl;
    std::cout << "Number of BETA determinants is now = " << betaDetsBinary.size() << std::endl;
    std::cout << " " << std::endl;
    int uniqueSetSize = uniqueDeterminantSet.size() ;
    std::cout << "Number of ALPHA determinants in set is now = " << uniqueSetSize << std::endl;
    std::cout << "First in Alpha unique SET =  " << *uniqueDeterminantSet.begin() << std::endl;
    std::cout << " " << std::endl;
    std::cout << "Number of Duplicates = " << numDuplicates << std::endl;
    std::cout << "Dets Generated in Total = " << numDetsCreated << std::endl;
    std::cout << "Dets Generated MINUS number of Duplicates = " << numDetsCreated - numDuplicates << std::endl;
    
    for (std::vector<int>::iterator it = walkerList.begin(); it != walkerList.end(); it++) {
        int element = *it ;
        std::cout << element << " ";
     }
    std::cout << " " << std::endl;
    */

