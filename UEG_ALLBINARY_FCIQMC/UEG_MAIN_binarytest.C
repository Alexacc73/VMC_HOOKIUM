/*! \mainpage FCIQMC HEAD PAGE
 * \section intro_sec Introduction
 *  ---> INCLUDED IN HEADER < GLOBAL_ORBITAL.H > is " const int ORB_SIZE = 125 " <--- 
    ---> NB - - NUMBER OF ORBITALS MUST BE KNOWN A PRIORI FOR "binaryManip.C* AND for KEsortedKpoints[size][3]
*/

#include "planeWaves.H"     /* Includes <iostream>  <set>  <math.h>  <stdlib.h>  <algorithm> */
#include "UEGHamiltonian.H" /* includes <bitset>  <climits>  <boost/math/special_functions/binomial.hpp>  <vector> */

#include <fstream>
#include <time.h>
#include <set>
#include <map>


/*
*-----> Parameters for initial setup of UEG <----- 
*/

/** A number, nobody knows what it means */
const double PI = 3.141592653589793;

/** rs controls density of the Electron Gas. 
* "rs" is the radius of the sphere whose volume is that of the cell divided by the number of electrons*/
const double rs = 0.5; 

/** Total number of electrons -> half allocated Alpha spin, the other half allocated Beta Spin.*/
const double numElectrons = 14; 
const int INTelectrons = numElectrons ;

/** Kc_CUTTOFF is the kinetic energy cutoff for the plane wave basis orbitals.
E.g, a cutoff of "2" will allow the orbital [4 0 0] but not [5 0 0]. Set cutoff = 2.4 for 57 Orbitals (114 Spin Orbitals) */
const double Kc_CUTTOFF = 2.4 ; 



/*
 *-----> Parameters for Main FCIQMC Algorithm <----- 
 */

/** delt is the Imaginary timestep for the propogation of the "walker" population */
const double delt = 0.00025 ;

/** Zeta is a damping parameter which controls the agressiveness of the "shift" in the variable shift mode of the algorithm */
const double zeta = 0.005 ;

/** AShift controls how frequently the shift is changed in response to the population in the variable shift mode (AShift = 1 means every step) */
const int AShift = 1 ;

/** Number of steps after which to terminate the algorithm*/
const int numSteps = 1000000;

/** After "walker critical" walkers have been spawned after a complete cycle (post annihilation) the variable shift mode is turned on */
const int walkerCritical = 300000;

/** initRefWalkers is the number of wlakers which are initially placed on the reference (i.e Hartree Fock) determinant to begin the spawning */
int initRefWalkers = 150;
long int pow2Array [ORB_SIZE];

/*
 *-----> OUTPUT FILES <----- 
 */
//const std::string FILE_shoulderPlot = "SHOULDER_114SO_rs0.5_WORKING.txt" ;
//const std::string FILE_shiftPlot = "SHIFT_114SO_rs0.5_WORKING.txt" ;

const std::string FILE_shoulderPlot = "SHOULDER_TEST.txt" ;
const std::string FILE_shiftPlot = "SHIFT_TEST.txt" ;




/**
Simply raises the input argument to the power 2, thus returning \f$ 2^x \f$
*/
inline constexpr std::uint64_t INLpow2 (std::uint64_t x)
{
    return std::uint64_t(1) << x;
}


/**
* This function purely resturns the position of a unique determinant within the list of Alpha and Beta determinants.
* It is necessary to know this, in order that we can add a new walker to the correct determinant (Our unique list is not
* sorted according to the walker number list, but the Alpha and Beta lists ARE sorted to associate with the walker number list)
*/
inline int INLgetPositionInList(std::pair<long int, long int>& uniqueDet, std::vector<long int>& alphaDetList, std::vector<long int>& betaDetList)
{
    int DETLIST_size = alphaDetList.size();
    int pos ;
    int j = 0;
    bool not_found = true ;

    long int alphaBIN = uniqueDet.first ;
    long int betaBIN = uniqueDet.second ;
    while( (j<DETLIST_size) && (not_found == true)  ){
        if( (alphaBIN == alphaDetList[j]) && (betaBIN == betaDetList[j]) ){
            pos = j ;
            not_found = false ;
        }
        j++ ;
    }
    if(not_found == true){
        std::cout << "OOPS POTISION FIND WENT WRONG!!!!!" << std::endl;
    }
    return pos ;
}


/**
* This function simply returns the number of bits, i.e, the number of ones in a binary number.
*/
inline size_t oneBitCount(long int n){
    std::bitset<sizeof(size_t) * CHAR_BIT> b(n);
    return b.count();
}


/**
* This inline function simply converts a base 10 integer (which represents the unique determinant)
* into its binary representation, but with type string. 
*/
inline  void INLdecimalToBinary(long int& decimal, std::string& binaryNum)
{
    binaryNum = std::bitset<ORB_SIZE>(decimal).to_string() ;
}

/**
* This inline function simply converts a string (of length defined at compile time) which represents the unique determinant
* into its decimal representation, type long int (up to 2^63).
*/
inline void INLbinaryToDecimal(std::string& binaryString, long int& decimal)
{
    decimal = std::bitset<ORB_SIZE>(binaryString).to_ulong() ;
}


/**
* This function is only used in conjunction with the Projector Energy routine,
* An very similar piece of code exists already within the Excitation operator
*/
inline int INLgetHijSign(int orb_i ,int orb_j ,int orb_a ,int orb_b)
{
    int index_i_excite;
    int index_j_excite;
    int sign = 0;
    int filledSO[INTelectrons/2] ;
    int tempFilledSO[INTelectrons/2] ;
    for(int i = 0; i<INTelectrons/2; i++){ /*Fill as a HF orbital, which we know it will be*/
        filledSO[i] = i;
        tempFilledSO[i] = i;
    }
    tempFilledSO[orb_i] = orb_a;
    tempFilledSO[orb_j] = orb_b;

    std::sort(tempFilledSO, tempFilledSO + INTelectrons/2 );
    for(int j = 0; j<INTelectrons/2; j++){
        if(tempFilledSO[j] == orb_a){
            index_i_excite = j ;
        }
        if(tempFilledSO[j] == orb_b){
            index_j_excite = j ;
        }
    }
    int idxSum = orb_i + orb_j + index_i_excite + index_j_excite;
    if(idxSum%2 == 0){// is EVEN!
        sign = 1;
    }
    else{
        sign = -1;
    }
    return sign;
}


/**
* This function is only used in conjunction with the Projector Energy routine,
* An very similar piece of code exists already within the Excitation operator
*/
inline int INLgetHijSignAB(int orb_i ,int orb_j ,int orb_a ,int orb_b)
{
    int sign;
    int alphai_excite;
    int betaj_excite;
    int filledSO[INTelectrons/2] ;
    int tempFilledSOA[INTelectrons/2] ;
    int tempFilledSOB[INTelectrons/2] ;
    for(int i = 0; i<INTelectrons/2; i++){ /*Fill as a HF orbital, which we know it will be*/
        filledSO[i] = i;
        tempFilledSOA[i] = i;
        tempFilledSOB[i] = i;
    }

    tempFilledSOA[orb_i] = orb_a;
    tempFilledSOB[orb_j] = orb_b;
    std::sort(tempFilledSOA, tempFilledSOA + INTelectrons/2) ;
    std::sort(tempFilledSOB, tempFilledSOB + INTelectrons/2) ;


    for(int j = 0; j<INTelectrons/2; j++){
        if(tempFilledSOA[j] == orb_a){
            alphai_excite = j;
        }
        if(tempFilledSOB[j] == orb_b){
            betaj_excite = j;
        }
    }
    int indexSum = orb_i + orb_j + alphai_excite + betaj_excite;
    if(indexSum%2 == 0){// is EVEN
        sign = 1;
    }
    else{
        sign = -1;
    }
    return sign;
}




/**This function represents first step in the population dynamics of the algorithm.
  * The spawning routine goes as follows, described verbatim from the original paper:
  * 
  * The spawning step: for each walker \f$ \alpha \f$ (located on
  *  \f$ D_{i \alpha} \f$), we select a (coupled) determinant \f$ D_j \f$ with
  * normalised probability pgen(j|iα) and attempt to
  * spawn a child there with probability
  \f[
    p_s(j|i_{\alpha}) = \frac{\delta \tau | K_{i_{\alpha}j} | }{p_{gen}j | i_{\alpha}} 
    (15)
   \f]
  * If a spawning event is successful, (i.e. if ps exceeds
  * a uniformly chosen random number between
  *0 and 1), then the sign of the child is determined
* by the sign of \f$ K_{i_{\alpha}j} \f$ and the sign of the parent: it
* is the same sign as the parent if \f$ K_{i_{\alpha}j} < 0 \f$ , and
* opposite to the parent otherwise. Our method to
* compute the generation probabilities \f$ p_{gen} \f$ is given
* in Appendix B. If \f$ p_s > 1 \f$ , then multiple copies of
* walkers are spawned on j (namely with probability
* 1, \f$ |ps| \f$ walkers are spawned, and with probability
* \f$ ps − |ps| \f$ an additional walker is spawned).
*/
void SPAWN( const double& cellLength, 
            std::vector<int>& trueWalkerList, 
            std::vector<int>& posWalkerList, 
            std::vector<int>& negWalkerList,
            double (&KEsortedList)[ORB_SIZE][3], 
            std::set< std::pair<long int, long int> >& uniqueDetSet, 
            std::pair< std::set< std::pair<long int, long int> >::iterator , bool >& result, 
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
    long int detChange = 0;
    long int alphaDet;
    long int betaDet;
    std::pair<long int, long int> iter ;
    bool ib_spinDifferent;

    int numDets = uniqueDetSet.size() ;
    for(int det = 0; det < numDets; det++){
        walkerNum = trueWalkerList[det] ;
        //  if(walkerNum != 0) {
        for(int num = 0; num < abs(walkerNum); num++){ /* abs(walkerNum) is VERY important */
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
                pGen = (2.0/( (numElectrons/2.0) *((numElectrons/2.0)-1.0)) ) * ( 2.0/(ORB_SIZE - (numElectrons/2)) ) ;
                pGen = pGen/4.0 ;
            }
            if(randChooseExcite == 1){/* i-->a AND j-->b are both BETA spin*/
                excitationSameSpinij_ab(index_i, index_a, index_j, index_b, INTelectrons, ORB_SIZE, betaDet, KEsortedList, SIGN);
                bool ib_spinDiffAlpha = false;
                pGen = (2.0/( (numElectrons/2.0) *((numElectrons/2.0)-1.0)) ) * ( 2.0/(ORB_SIZE - (numElectrons/2)) ) ;
                pGen = pGen/4.0 ;
            }
            if( (randChooseExcite == 2)  || (randChooseExcite == 3) ){ /* i-->a is ALPHA,  j-->b is BETA*/
                excitationAlpha_iaBeta_jb(index_i, index_a, index_j, index_b, INTelectrons, ORB_SIZE, alphaDet, betaDet, KEsortedList, SIGN);
                bool ib_spinDifferent = true;
                pGen = (1.0/((numElectrons/2.0)*(numElectrons/2.0)) ) * ( 1.0/(ORB_SIZE - (numElectrons/2)) ) ;
                pGen = pGen/2.0 ;
            }
    
            if(index_b != -1){ /*Finding a k consevring orbital was SUCCESSFUL*/
                /*1) Update newly excited determinants*/
                if(randChooseExcite==0){
                    detChange = -pow2Array[index_i] - pow2Array[index_j] + pow2Array[index_a] + pow2Array[index_b] ;
                    alphaDet += detChange;
                    //alphaDet -= pow2Array[  index_i ] ;//= '0' ;
                    //alphaDet += pow2Array[  index_a ] ;//= '1' ;
                }
                if(randChooseExcite==1){
                    detChange = -pow2Array[index_i] - pow2Array[index_j] + pow2Array[index_a] + pow2Array[index_b] ;
                    betaDet += detChange;
                    //betaDet -= pow2Array[ index_i ] ;// = '0' ;
                    //betaDet += pow2Array[ index_a ] ;// = '1' ;
                }
                if( (randChooseExcite == 2)  || (randChooseExcite == 3) ){
                    alphaDet -= pow2Array[ index_i] ; // = '0' ;
                    alphaDet += pow2Array[ index_a] ; // = '1' ;
                    betaDet  -= pow2Array[ index_j] ; // = '0' ;
                    betaDet  += pow2Array[ index_b] ; // = '1' ;
                }
                   
                HijElement = 0;
                Di_H_Dj(cellLength , KEsortedList, index_i, index_a, index_b, ib_spinDifferent, HijElement) ;
                HijElement *= SIGN ;
                pSpawn = (delt * fabs(HijElement) ) / pGen ;
                intpSpawn = floor(pSpawn) ;
                remainderpSpawn = pSpawn - intpSpawn ;
                metropolisRand = ((double) rand() / RAND_MAX) ; //(fastrand() % 100000 )/100000 ;

                /* SPAWN: YES or NO*/
                if(remainderpSpawn >= metropolisRand){/* --> SPAWING EVENT SUCCESSFUL <--*/
                    result = uniqueDetSet.insert( std::make_pair(alphaDet, betaDet) ) ;

                    if(result.second == true){/* New Determinant Found */
                    /*All 3 walker lists *MUST* each increase their size by 1,
                    hence the .push_back(0) for each .puch_back(1) 1*/
                        trueWalkerList.push_back(0) ;
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
                    result = uniqueDetSet.insert( std::make_pair(alphaDet, betaDet) ) ;

                    if(result.second == true){/* New Determinant Found */
                        trueWalkerList.push_back(0) ;
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
            } /*End of IF(index_b != -1)*/
                // else{/* COULD NOT find a k consevring orbital - SET Hij = 0 and carry on as usual!*/
                    /*Do Nothing - we know with certainty that pSpawn = 0 since HijElement = 0*/
                //}
        }/*End of inner walker FOR loop*/
        // end if }
    } /*End of main FOR loop*/
}/* end of function SPAWN */




/** This routine represents the second step in the population dynamics algorithm. The following description is verbatim
*from the original paper:
*
*The diagonal death/cloning step: for each (parent)
*walker compute
\f[
p_{d}(i_{\alpha}) = \delta \tau( K_{i_{\alpha} i_{\alpha}} - S )
(16)
\f]
* 
* If \f$ p_d > 0 \f$, the walker dies with probability \f$ p_d \f$,
* and if \f$ p_d < 0 \f$ the walker is cloned with probability
* \f$ |pd| \f$. The death event happens immediately,
* and such a parent does not participate in the following
* (annihilation) step to be described shortly.
* Cloning events are quite rare, and only occur for
* \f$ S > 0 \f$, and even then only on determinants for
* which \f$ \langle  D_i | K | D_i \rangle < S  \f$ . In simulations where we desire
* to grow the number of walkers rapidly, a positive
* value of S is adopted, and this can lead to
* cloning events. However, more often, the value of
* S is negative (as it tries to match the correlation
* energy), and in such cases there can be no cloning
* events at all.
*/
void DEATH_CLONE(const double& cellLength, 
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
            RandRound =  ((double) rand() / RAND_MAX) ; /** Generates a number between (0,1] */
            INTnDeath = floor(nDeath + RandRound); 
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


/** This represents the final (third) step in the population dynamics algorithm. The following description is verbatim
from the original paper:

* The annihilation step: In this (final) part of the
* algorithm, we run over all (newly-spawned, cloned
* and surviving parent) walkers, and annihilate pairs
* of walkers of opposite sign which are found to be on
* the same determinant. Each time an annihilation
* event occurs, the corresponding pair is removed
* from the list of walkers, and the total number of
* walkers \f$ N_w \f$ reduced by two.
*/
void ANNIHILATION(int& step, 
                std::vector<int>& trueWalkerList, 
                std::vector<int>& posWalkerList, 
                std::vector<int>& negWalkerList,
                std::set< std::pair<long int, long int> >& uniqueDeterminantSet, 
                std::vector<long int>& alphaDetsBinary, 
                std::vector<long int>& betaDetsBinary){
    long int alphaBin;
    long int betaBin;
    int numDets = 0;
    int NEWnumDets = 0;
    bool prune = false;
    if(step%10==1){
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
        /* Make copies of all the lists */
        std::vector<int> TWL_copy;
        std::vector<long int> alphaDets_copy;
        std::vector<long int> betaDets_copy;
      
        /* Move the values from the original lists to the copies */
        for(int num = 0; num < numDets; num++){
            TWL_copy.push_back(trueWalkerList[num]);
            alphaDets_copy.push_back(alphaDetsBinary[num]);
            betaDets_copy.push_back(betaDetsBinary[num]);
        }
        //copy(trueWalkerList.begin(), trueWalkerList.end(), back_inserter(TWL_copy));
        //copy(alphaDetsBinary.begin(), alphaDetsBinary.end(), back_inserter(alphaDets_copy));
        //copy(betaDetsBinary.begin(), betaDetsBinary.end(), back_inserter(betaDets_copy));
        /* Empty the old lists of their contents, ready to be refilled with the pruned ones */
        trueWalkerList.clear();
        alphaDetsBinary.clear();
        betaDetsBinary.clear();
        posWalkerList.clear();
        negWalkerList.clear();
        NEWnumDets = TWL_copy.size();
        /* CREATE NEW LISTS WITH ONLY DETERMINANTS WITH A NON-ZERO NUMBER OF WALKERS ASSOCIATED WITH THEM  */
        for(int el = 0; el<NEWnumDets; el++){
            if(TWL_copy[el] != 0){
                trueWalkerList.push_back(TWL_copy[el]);
                alphaDetsBinary.push_back(alphaDets_copy[el]);
                betaDetsBinary.push_back(betaDets_copy[el]);
                posWalkerList.push_back(0);
                negWalkerList.push_back(0);
            }
        }

        //posWalkerList.resize(trueWalkerList.size()) ;
        //negWalkerList.resize(trueWalkerList.size()) ;
        uniqueDeterminantSet.clear();

        for(int det = 0; det < alphaDetsBinary.size(); det++){
            alphaBin = alphaDetsBinary[det];
            betaBin = betaDetsBinary[det];
            uniqueDeterminantSet.insert( std::make_pair(alphaBin, betaBin) );
        }
    }/*End Prune*/

}


/* - - - - - AUXILLIARY FUNCTIONS - - - - -*/

/**
* Simply counts the current number of walkers associated with all determinants, irrespective of sign.
* Returns: \f$  N_w = \sum_{i} | N_i |  \f$
*/
int totalWalkerNumber( std::vector<int>& trueWalkerList){
  int count = 0;
  int numDets = trueWalkerList.size() ;
  for(int det = 0; det < numDets; det++){
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
 E_{proj} =  \sum_{J \neq 0} \Big \langle D_J | H | D_0 \Big \rangle \frac{N_J (\tau)}{N_0 (\tau)} 
\f]
*/
double projectorEnergy(const double& cellLength, 
                       long int& HF_BINARY,
                       std::vector<int>& trueWalkerList,
                       std::vector<long int>& alphaDetsBinary,
                       std::vector<long int>& betaDetsBinary,
                       double (&KEsortedList)[ORB_SIZE][3] ){
    int idx_i, idx_j, idx_a, idx_b;
    int alphaBits = 0, betaBits = 0, excitorCount = 0;
    int SIGN = 1;
    int idx;
    long int HFDet = HF_BINARY;
    long int alphaDetJ, betaDetJ, XORalpha, XORbeta;
    long int ijBin, abBin, iBin, jBin;
    long int aBin, bBin;
    bool ISdoubleExcitation = false;
    bool found_idx_i = false, found_idx_a = false;
    bool IB_spinDifferent = false;
    double HijElement = 0, walkerNumJ = 0, EProjSum = 0;
    int INTwalkerNum;
    std::string ijSTR(ORB_SIZE, ' '); 
    std::string abSTR(ORB_SIZE, ' ');
  
    int numDets = alphaDetsBinary.size();
    double refWalkerNum = trueWalkerList[0]; 
    long int one = 1;
    
    for(int det = 1; det<numDets; det++){ /*Begin at det = 1, since we do not care about < D_0 | H | D_0 >*/
        INTwalkerNum = trueWalkerList[det];
        walkerNumJ = trueWalkerList[det];
        if( INTwalkerNum != 0 ){

            
            ISdoubleExcitation = false;
            found_idx_i = false;
            found_idx_a = false;
            idx_i = 0;
            idx_j = 0;
            idx_a = 0;
            idx_b = 0;
            alphaDetJ = alphaDetsBinary[det];
            betaDetJ = betaDetsBinary[det];
            XORalpha = HFDet ^ alphaDetJ;
            XORbeta = HFDet ^ betaDetJ;
            alphaBits = oneBitCount(XORalpha);
            betaBits = oneBitCount(XORbeta);
            excitorCount = alphaBits + betaBits;
            if(excitorCount == 4){
                ISdoubleExcitation = true;
            }

            if(ISdoubleExcitation == true){ /* Double Excitation Found */
                if( (alphaBits == 4) ){ /* search for ij, ab, in ALPHA only!*/
                    bool IB_spinDifferent = false;
                    ijBin = HFDet & XORalpha;
                    abBin = alphaDetJ & XORalpha;
                    idx = 0;
                    for(int i = 0; i<ORB_SIZE; i++){
                        //idx = ORB_SIZE - i -1;
                        if( (found_idx_i == true) && ((ijBin & (one<<i)) != 0)  ){
                            idx_j = i;
                        }
                        if( (found_idx_i == false) && ((ijBin & (one<<i)) != 0) ){
                            idx_i = i;
                            found_idx_i = true;
                        }
                        if( (found_idx_a == true) &&  ((abBin & (one<<i)) != 0) ){
                            idx_b = i;
                        }
                        if( (found_idx_a == false) && ((abBin & (one<<i)) != 0) ){
                            idx_a = i;
                            found_idx_a = true;
                        }
                    }
                    SIGN = INLgetHijSign(idx_i, idx_j, idx_a, idx_b);
                }/*End ALPHA search for ij, ab*/

                if( (betaBits == 4) ){ /* search for ij, ab, in BETA only!*/
                    bool IB_spinDifferent = false;
                    ijBin = HFDet & XORbeta;
                    abBin = betaDetJ & XORbeta;
                    idx = 0;
                    for(int i = 0; i<ORB_SIZE; i++){
                        //idx = ORB_SIZE - i -1;
                        if( (found_idx_i == true) && ((ijBin & (one<<i)) != 0)  ){
                            idx_j = i;
                        }
                        if( (found_idx_i == false) && ((ijBin & (one<<i)) != 0) ){
                            idx_i = i;
                            found_idx_i = true;
                        }
                        if( (found_idx_a == true) &&  ((abBin & (one<<i)) != 0) ){
                            idx_b = i;
                        }
                        if( (found_idx_a == false) && ((abBin & (one<<i)) != 0) ){
                            idx_a = i;
                            found_idx_a = true;
                        }
                    }
                    SIGN = INLgetHijSign(idx_i, idx_j, idx_a, idx_b);
                }/*End ALPHA search for ij, ab*/

                if( (alphaBits == 2) || (betaBits == 2) ){ /* i->a are alpha, but i->b are beta*/
                    bool IB_spinDifferent = true;
                    iBin = HFDet & XORalpha;
                    jBin = HFDet & XORbeta;
                    aBin = alphaDetJ & XORalpha;
                    bBin = betaDetJ & XORbeta;
                    for(int i = 0; i<ORB_SIZE; i++){
                        if( iBin - (one<<i) == 0 ){
                            idx_i = i;
                        }
                        if( jBin - (one<<i) == 0 ){
                            idx_j = i;
                        }
                        if( aBin - (one<<i) == 0 ){
                            idx_a = i;
                        }
                        if( bBin - (one<<i) == 0 ){
                            idx_b = i;
                        }
                    }
                    SIGN = INLgetHijSignAB(idx_i, idx_j, idx_a, idx_b);
                }

                HijElement = 0;
                Di_H_Dj(cellLength, KEsortedList, idx_i, idx_a, idx_b, IB_spinDifferent, HijElement) ;
                HijElement *= SIGN;
                EProjSum += (HijElement * walkerNumJ ) ;
            }
        }
    }/* End of FOR loop */
    EProjSum *= (1.0/refWalkerNum) ;
    return EProjSum ;
}


/**
------------------------------------------------
M A I N  -  S T A R T S  -  H E R E 
------------------------------------------------
*/
int main(void){
    srand(492831);

    for(int i = 0; i<ORB_SIZE; i++){
        pow2Array[i] = INLpow2(i);
    }

	const double cellVolume = (numElectrons*PI*4.0*rs*rs*rs)/3.0;
	const double cellLength = getCellLength(cellVolume); /* Return cell length in units of rs */
	int spinOrbitals;
    double kpoints[1000][3];

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
    INLdecimalToBinary(zero, HF_STRING) ; 
    for(int el = 0; el < INTelectrons/2; el++){
        HF_STRING[ORB_SIZE-1-el] = '1';
    }
    
    long int HF_binary;
    INLbinaryToDecimal(HF_STRING, HF_binary);
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


    /* EXCIATATION TESTER 
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
    */





    /* - - - - - - - - - - - - M A I N   L O O P - - - - - - - - - - - - */
    int beginAverageStep;
    double SHIFT = 0 ;
    double PROJECTOR = 0;
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
    clock_t CLOCKSUM = 0;
    bool shiftOn = false;

    std::ofstream shoulderplot;
    std::ofstream shiftPlot;
    shoulderplot.open(FILE_shoulderPlot); /* File is defined at top of the code as const std::String*/
    shiftPlot.open(FILE_shiftPlot);  
   
    int PRINT_STEPS = 50;
    for(int i = 0; i < numSteps; i++){

        start = clock();
        SPAWN(cellLength, walkerList, posChildList, negChildList, KEsortedKpoints, uniqueDeterminantSet, result, alphaDetsBinary, betaDetsBinary );
        DEATH_CLONE(cellLength, walkerList, KEsortedKpoints, SHIFT, alphaDetsBinary, betaDetsBinary, HFEnergy);
        ANNIHILATION(i, walkerList, posChildList, negChildList, uniqueDeterminantSet, alphaDetsBinary, betaDetsBinary);
        end = clock();
        CLOCKSUM += (end-start);

        
        currentPopulation = totalWalkerNumber(walkerList) ;
        walkerNUMTracker.push_back( currentPopulation );

        referenceWalkers = walkerList[0] ;
        popToRefRatio = currentPopulation / referenceWalkers ;


        if(i%PRINT_STEPS ==0){
            std::cout << "-------------> STEP : " << i << std::endl;
            std::cout << "CURRENT POPULATION = " << walkerNUMTracker[i] << std::endl;
            std::cout << "REFERENCE WALKERS = " << referenceWalkers << std::endl;
            std::cout << "CURRENT NUMBER OF DETERMINANTS = " <<  uniqueDeterminantSet.size() << std::endl;
            std::cout << "CURRENT NUMBER OF ALPHAS = " <<  alphaDetsBinary.size() << std::endl;
            std::cout << "CURRENT NUMBER OF BETAS = " <<  betaDetsBinary.size() << '\n' << std::endl;
            std::cout << "CURRENT SHIFT ENERGY ESTIMATOR = " << SHIFT << std::endl;
        }


        /* TURN SHIFT ON NOW, *THEN* perform spawn-death/clone-annihilation */

        if( (currentPopulation >= walkerCritical ) && !shiftOn ){
            shiftOn = true;
            beginAverageStep = i*2;
            std::cout << "THIS TEXT SHOULD APPEAR EXACTLY ONCE" << std::endl;
        }
        if(shiftOn && ( (i+1)%AShift == 0) ){ 
            SHIFT = variableShift(delt, AShift, i, zeta, walkerNUMTracker, SHIFTTracker) ;
        }
        SHIFTTracker.push_back(SHIFT) ;
        PROJECTOR = projectorEnergy(cellLength, HF_binary, walkerList, alphaDetsBinary, betaDetsBinary, KEsortedKpoints);

        shoulderplot << currentPopulation << " " << popToRefRatio << std::endl ;
        shiftPlot << SHIFT << " " << PROJECTOR << std::endl;

        
        
        if(i%PRINT_STEPS == 0){
            elapsedTime = CLOCKSUM;
            std::cout<< "Elapsed time for 50 triple loop: " << elapsedTime/CLOCKS_PER_SEC << std::endl;
            CLOCKSUM = 0;
        }

    }
    shoulderplot.close();
    shiftPlot.close();

    /* PRINT OUT THE WALKER POPULATION PER DETERMINANT
    for(int i = 0; i< walkerList.size(); i++ ){
        std::cout << "Walker on DET " << i << " = " << walkerList[i] << std::endl;
    }
    */

    return 1;
}








