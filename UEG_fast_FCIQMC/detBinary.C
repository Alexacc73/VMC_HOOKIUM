#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <set>
#include <map>
#include <math.h>
#include <time.h>
#include <bitset>
#include <boost/math/special_functions/binomial.hpp>
#include <cmath>
#include <iomanip>
#include <climits>

//#include <random>

using namespace boost;
using namespace math;


//std::random_device rd; //seed
//std::mt19937 gen(rd()); //ruleset for rd(merzenne twister)
//std::uniform_int_distribution<> rng_dice(1, 6); ///rng2 range



double n = 19 ;
double r = 7 ;
int INTr = r;

std::string decimalToBinary(long int decimal){
	std::string binaryNum = std::bitset<57>(decimal).to_string() ;
	return binaryNum ;
}

long int binaryToDecimal(std::string binaryString){
	long int decimal = std::bitset<57>(binaryString).to_ulong();
	return decimal;
}

size_t oneBitCount(size_t n){
	std::bitset<sizeof(size_t) * CHAR_BIT> b(n);
    return b.count();
}


std::set< std::pair<long int, long int> > myDet;
std::set< std::string > mydetSTR;
std::set< std::pair<std::bitset<57>, std::bitset<57>> > myDetBit;


inline constexpr std::uint64_t pow2 (std::uint64_t i)
{
    return std::uint64_t(1) << i;
}


int main(void){
	std::cout << std::setprecision(100);
	//srand(492831) ;
	srand(time(NULL)) ;


	double RAND;
	double badrand;
	int       rand_int = 1000000;
	double rand_divide = 1000000;

    double RAND_sum = 0;
    double badrand_sum = 0 ;
    int steps = 5;
	for (int i=0; i<steps; i++){
		RAND = ((double) rand() / RAND_MAX) ;
		//badrand = ((double)(rand()% rand_int)) / rand_divide ;
		badrand = (rand() % 4 ) ;
		//std::cout << badrand << std::endl;

		RAND_sum += RAND;
		badrand_sum += badrand;
		//std::cout << "Good = " << RAND << std::endl;
		//std::cout << "Other = " << badrand << '\n' << std::endl;
	}
	//std::cout << "RAND mean = " << RAND_sum / ((double)steps) << std::endl;
	//std::cout << "badrand mean = " << badrand_sum / ((double)steps) << std::endl;

    long int num = 127;
    std::string binaryNum = decimalToBinary(num);
    std::cout << binaryNum << std::endl;
    for(int i = 0; i < binaryNum.size(); i++){
    	binaryNum[i] = '1';
    }
    binaryNum[1] = '1';
    long int det1 = binaryToDecimal(binaryNum);

    long int alphaHF = 127;
    long int betaHF = 127;
    std::string alphaSTR = decimalToBinary(alphaHF);

    //myDet.insert(std::make_pair(alphaHF, betaHF)); //1
    //myDet.insert(std::make_pair(179387, 131029));  //2
    //mydetSTR.insert( decimalToBinary(101029) );
    //mydetSTR.insert( decimalToBinary(127) );

    int cycle = 1000;
    long int insert = 0;
    std::bitset<57> INBIT (1);
    //myDetBit.insert(std::make_pair(1, 2));
    
    for(long int i = 0; i<cycle; i++){
    	std::bitset<57> INBIT (i);
    	mydetSTR.insert( decimalToBinary(insert) );
    	myDet.insert(std::make_pair(insert, insert));
    	insert += 1;
    }

    clock_t start1;
    clock_t end1;
    clock_t start2;
    clock_t end2;
    double elapsedTime1;
    double elapsedTime2;

    long int num1 = 10;
    //TIME THINGS!!!!!
    start1 = clock();
    for(int j = 0; j<cycle; j++){
    	mydetSTR.insert( decimalToBinary(num1) ) ;
    }
    end1 = clock();

    start2 = clock();
    for(int k = 0; k<cycle; k++){
    	myDet.insert(std::make_pair(num1, num1)) ;
    }
    end2 = clock();

    elapsedTime1 = end1-start1;
    elapsedTime2 = end2-start2 ;
    std::cout<< "Elapsed time for STRING compare: " << elapsedTime1/CLOCKS_PER_SEC << std::endl;
    std::cout << " " << std::endl;
    std::cout<< "Elapsed time for LONG INT compare: " << elapsedTime2/CLOCKS_PER_SEC << std::endl;
    std::cout << " " << std::endl;



    std::cout << "Det number long int = " << myDet.size() << std::endl;
    std::cout << "Det number string = " << mydetSTR.size() << std::endl;

    std::pair< std::set< std::pair<long int, long int> >::iterator , bool > result;
    result = myDet.insert(std::make_pair(179387, 131029));
    if(result.second == false){
    	std::cout << "False Bool = " << result.second << std::endl;
    }
    result = myDet.insert(std::make_pair(betaHF, 10287));
    std::pair<long int, long int> iter = *result.first;
    if(result.second == true){
    	std::cout << "ELement added successfully = " << iter.first << " " << iter.second << std::endl;
    }

    std::cout << "2^57 = " << pow2(63) << std::endl;

    alphaHF = 127;
    alphaHF +=  pow2(40);
    alphaHF += pow2(41);
    alphaHF -= pow2(1);
    alphaHF -= pow2(2);

    std::bitset<64> foo (127);
    std::bitset<64> one (1);
    foo.set(63,1);
    std::cout << foo << std::endl;
   
    return 1;

   

}

