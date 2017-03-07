#include <iostream>
#include <set>
#include <math.h>
#include <stdlib.h>
#include <algorithm>

void createPlaneWaveOrbitals(double [][3], const double, int& );
void kPointsEnergySort(double[][3], double[][3], const int);
void PRINTORBITALS(const int, double [][3]);
const double getCellLength(const double);

