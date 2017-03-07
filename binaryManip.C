#include "binaryManip.H"
#include "GLOBAL_ORBITAL.H"



void decimalToBinary(long int decimal, std::string& binaryNum){
	binaryNum = std::bitset<ORB_SIZE>(decimal).to_string() ;
	//return binaryNum ;
}

void binaryToDecimal(std::string binaryString, long int& decimal){
	decimal = std::bitset<ORB_SIZE>(binaryString).to_ulong() ;
	//return decimal ;
}

size_t oneBitCount(size_t number){
	std::bitset<sizeof(size_t) * CHAR_BIT> b(number) ;
    return b.count() ;
}

double binomialnCr(double n, double r){
	double nCr = boost::math::binomial_coefficient<double>(n, r) ; 
	return nCr ;
}

/*
void concatStrings(std::string& ALPHA, std::string& BETA, std::string& COMBINATION){
    for(int i = 0; i<ORB_SIZE; i++){
        COMBINATION[i] = ALPHA[i] ;
        COMBINATION[i+ORB_SIZE+1] = BETA[i] ;
    }
    COMBINATION[ORB_SIZE] = ' ' ;
}
*/

int getPositionInList(std::pair<long int, long int>& uniqueDet, std::vector<long int>& alphaDETLIST, std::vector<long int> betaDETLIST ){
    //std::string ALPHA(ORB_SIZE, ' ') ;
    //std::string BETA(ORB_SIZE, ' ') ;
    //for(int i = 0; i<ORB_SIZE; i++ ){
    //    ALPHA[i] = uniqueDet[i] ;
    //    BETA[i] = uniqueDet[i+ORB_SIZE+1] ;
    //}

    int DETLIST_size = alphaDETLIST.size();
    int position = 0;
    int j = 0;
    bool not_found = true ;

    long int alphaBIN = uniqueDet.first ;
    //binaryToDecimal(ALPHA, alphaBIN);
    long int betaBIN = uniqueDet.second ;
    //binaryToDecimal(BETA, betaBIN);
    while( (j<DETLIST_size) && (not_found == true)  ){
        if( (alphaBIN == alphaDETLIST[j]) && (betaBIN == betaDETLIST[j]) ){
            position = j ;
            not_found = false ;
        }
        j++ ;
    }
    if(not_found == true){
        position = -1;
        std::cout << "OOPS POTISION FIND WENT WRONG!!!!!" << std::endl;
    }
    return position ;
}