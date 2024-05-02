// AtomHic library files
#include <AtomicSystem.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <cmath>
#include "MathTools.h"

using namespace std;

int main(int argc, char *argv[])
{
	// check the number of argument
	if( argc < 6 ){
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename GBNormalDirection CrystalType AtomicOutputFilename OrderParameterDensityFilename GaussianOrderParameterDensityFilename StatdataFilename" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename = argv[1];
	string NormalDir = argv[2];
	string DensityFilename = argv[3];
	string GaussDensityFilename = argv[4];
	string StatdataFilename = argv[5];
	AtomicSystem MySystem(InputFilename);
	unsigned int nbPts_i = 1000;// TODO maybe set a file reading sigma nbPts, lsph and rc for the different systems
	double Lz = 1.4825682700766561e+02;
        double sigma = 2.;
	unsigned int ind_DisoDens = MySystem.Compute1dDensity("Disorder", NormalDir, sigma, nbPts_i);
	double pos1_i = 55.;
	double pos2_i = 100.;
	unsigned int ind1_i = round(pos1_i*nbPts_i/Lz);
	unsigned int ind2_i = round(pos2_i*nbPts_i/Lz);
	cout << ind1_i << " " << ind2_i << endl;
	unsigned int maxpts = 150;
	double max1 = MySystem.getDensityProf(ind_DisoDens)[ind1_i*2];
	double max2 = MySystem.getDensityProf(ind_DisoDens)[ind2_i*2];
	double facred = 0.6;
	double low_1, low_2, hi_1, hi_2;
	bool found = false;
	for(unsigned int i=0;i<maxpts;i++){
		if( MySystem.getDensityProf(ind_DisoDens)[(ind1_i-i)*2] < max1*facred ){
			low_1 = MySystem.getDensityProf(ind_DisoDens)[(ind1_i-i)*2+1];
			found = true;
			break;
		}
	}
	if( !found ) cout << "low1 not found" << endl;
	found = false;
	for(unsigned int i=0;i<maxpts;i++){
		if( MySystem.getDensityProf(ind_DisoDens)[(ind1_i+i)*2] < max1*facred ){
			hi_1 = MySystem.getDensityProf(ind_DisoDens)[(ind1_i+i)*2+1];
			found = true;
			break;
		}
	}
	if( !found ) cout << "hi1 not found" << endl;
	found = false;
	for(unsigned int i=0;i<maxpts;i++){
		if( MySystem.getDensityProf(ind_DisoDens)[(ind2_i+i)*2] < max2*facred ){
			hi_2 = MySystem.getDensityProf(ind_DisoDens)[(ind2_i+i)*2+1];
			found = true;
			break;
		}
	}
	if( !found ) cout << "hi2 not found" << endl;
	found = false;
	for(unsigned int i=0;i<maxpts;i++){
		if( MySystem.getDensityProf(ind_DisoDens)[(ind2_i-i)*2] < max2*facred ){
			low_2 = MySystem.getDensityProf(ind_DisoDens)[(ind2_i-i)*2+1];
			found = true;
			break;
		}
	}
	if( !found ) cout << "low2 not found" << endl;
	MySystem.Print1dDensity(DensityFilename, "Disorder");
	ofstream writefile(StatdataFilename);
	writefile << hi_1-low_1 << " " << hi_2-low_2;
	writefile.close();
	return 0;
}
