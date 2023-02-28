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
	if( argc < 3 ){
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename CrystalName AtomicOutputFilename" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename = argv[1];
	string CrystalName = argv[2];
	string OutputFilename = argv[3];
	Bicrystal MySystem(InputFilename, CrystalName); // for the moment need to be a bicrystal object (to have the computeAuxiliary methods), it has no sense but its complicated to change (it should be an AtomicSystem instead but the AtomicSystem cannot have a ComputeAuxiliary)
	MySystem.setAux_vec(MySystem.get_CA()->StructuralAnalysis_Steinhardt(), MySystem.get_CA()->getNumberRefDef_Steinhardt(), "sAB");
	MySystem.printSystem_aux(OutputFilename, "sAB");
	return 0;
}
