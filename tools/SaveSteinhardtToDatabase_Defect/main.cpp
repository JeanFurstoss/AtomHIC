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
	if( argc < 2 ){
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename CrystalType filename" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename = argv[1];
	string CrystalType = argv[2];
	string filename = argv[3];
	Bicrystal MySystem(InputFilename, CrystalType);
	MySystem.get_CA()->SaveSteinhardtParamToDatabase_Defect(CrystalType,filename);
	return 0;
}
