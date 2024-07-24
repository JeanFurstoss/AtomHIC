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
	// check the number of argument
	if( argc < 8 ){
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename GBNormalDirection CrystalType AtomicOutputFilename OrderParameterDensityFilename GaussianOrderParameterDensityFilename StatdataFilename" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename = argv[1];
	string NormalDir = argv[2];
	string CrystalType = argv[3];
	string OutputFilename = argv[4];
	string DensityFilename = argv[5];
	string GaussDensityFilename = argv[6];
	string StatdataFilename = argv[7];
	Bicrystal MySystem(InputFilename, NormalDir, CrystalType);
	MySystem.printSystem_aux(OutputFilename, "Disorder");
	MySystem.Print1dDensity(DensityFilename, "GBProfile");
	MySystem.Print1dDensity(GaussDensityFilename, "GBProfile_Gauss");
	ofstream writefile(StatdataFilename);
	writefile << MySystem.getGBPos1() << " " << MySystem.getGBwidth1() << " " << MySystem.getExcessVol();
	writefile.close();
	return 0;
}
