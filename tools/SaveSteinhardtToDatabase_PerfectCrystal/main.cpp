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
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename CrystalType(The crystal type has to be defined in /data/Crystal/)" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename = argv[1];
	string CrystalType = argv[2];
	Bicrystal MySystem(InputFilename, CrystalType); // for the moment need to be a bicrystal object (to have the computeAuxiliary methods), it has no sense but its complicated to change (it should be an AtomicSystem instead but the AtomicSystem cannot have a ComputeAuxiliary)
	MySystem.get_CA()->SaveSteinhardtParamToDatabase_PerfectCrystal(CrystalType);
	return 0;
}
