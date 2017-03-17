#include "planeWaves.H"




void createPlaneWaveOrbitals(double kPoints[][3], const double KECutoff, int& numSpinOrbitals){
	numSpinOrbitals = 0;
	int kval = 0 ;
	double kineticEnergy  = 0.0 ;
	int maxK = ceil(KECutoff*KECutoff) ;
	for(int n = -maxK; n < maxK+1; n++){
		for(int m = -maxK; m < maxK+1; m++){
			for(int l = -maxK; l < maxK+1; l++){
				kineticEnergy = pow( (n*n + m*m + l*l) , 0.5 ) ;
				if(kineticEnergy <= KECutoff){
					kPoints[kval][0] = n;
					kPoints[kval][1] = m;
					kPoints[kval][2] = l;
					numSpinOrbitals += 2;
					kval++;
				}
			}
		}
	}
}

void kPointsEnergySort(double kPoints[][3], double orderedkPoints[][3], const int kPointLength ){
	std::set<int> kineticList;
	for(int j = 0;  j < kPointLength; j++ ){
		double currentKinetic = kPoints[j][0]*kPoints[j][0] + kPoints[j][1]*kPoints[j][1] + kPoints[j][2]*kPoints[j][2];
		kineticList.insert(currentKinetic);
	}
	int k = 0;
	while(!kineticList.empty()){
		double currentKinetic =  *kineticList.begin() ;
		for(int i = 0; i<kPointLength; i++){
			double testKinetic = kPoints[i][0]*kPoints[i][0] + kPoints[i][1]*kPoints[i][1] + kPoints[i][2]*kPoints[i][2];
			if(testKinetic == currentKinetic){
				orderedkPoints[k][0] = kPoints[i][0] ;
				orderedkPoints[k][1] = kPoints[i][1] ;
				orderedkPoints[k][2] = kPoints[i][2] ;
				k++;
			}
		}
		kineticList.erase(kineticList.begin());
	}
}

void PRINTORBITALS(const int klistLength, double klist[][3] ){
	for(int i = 0; i < klistLength; i++){
		std::cout << "[" << klist[i][0] << klist[i][1] << klist[i][2] << "] " << std::flush;
	}
	std::cout << " " << std::endl;
}

const double getCellLength(const double cellVolume){
  const double length = pow(cellVolume, (1.0/3.0) );
  return length;
}