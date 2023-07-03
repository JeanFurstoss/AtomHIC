// AtomHic library files
#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"

using namespace std;

int main(int argc, char *argv[])
{
	unsigned int dim = 6.;
	vector<double> row(6,0.);
	vector<vector<double>> mat;
	vector<vector<double>> mat2;
	for(unsigned int i=0;i<dim;i++){
		mat.push_back(row);
		mat[i][i] = i*2+20;
		mat2.push_back(row);
		mat2[i][i] = 25+pow((double) i,2.);
	}
	for(unsigned int i=0;i<dim;i++){
		for(unsigned int j=i+1;j<dim;j++){
			mat[i][j] = j+i;
			mat[j][i] = j+i;
			mat2[i][j] = 2*j+4*i;
			mat2[j][i] = j+1*i;
		}
	}
	MathTools MT;
	MT.printMatVec(mat);
	double det;
	vector<vector<double>> inv;
	MT.invMat_LU(mat,inv,det);
	MT.printMatVec(inv);
	//if( argc < 3 ){
	//	cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename number_of_struct number_of_struct*(AtomIndexListFileName name_of_struct) CrystalType OutputFilename" << endl;
	//	cerr << "or: AnalyzeBicrystal_ARGS AtomicInputFilename number_of_struct number_of_struct*(x_lo x_hi y_lo y_hi z_lo z_hi name_of_struct) CrystalType OutputFilename" << endl; 
	//	return EXIT_FAILURE;
	//}
	//
	//InputFilename = argv[1];
	//istringstream iss_nbS(argv[2]);
	//iss_nbS >> nbStruct;
	return 0;
}
