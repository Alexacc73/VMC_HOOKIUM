#include <stdlib.h>
#include <iostream>
#include <bitset>
#include <vector>
#include <climits>
#include <boost/math/special_functions/binomial.hpp>

void decimalToBinary(long int, std::string&);
void binaryToDecimal(std::string, long int&);
size_t oneBitCount(size_t);
double binomialnCr(double, double);
//void concatStrings(std::string& , std::string&, std::string& );
int getPositionInList(std::pair<long int, long int>& , std::vector<long int>& , std::vector<long int> );

