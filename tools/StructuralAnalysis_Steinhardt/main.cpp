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
		cerr << "Usage: StructuralAnalysis_Steinhardt AtomicInputFilename CrystalName AtomicOutputFilename" << endl;
		cerr << "this executable performs structural analysis using the Steinhardt parameters in the files /data/Steinhardt/Crystal/*.dat and using Euclidian distances, this method is quite obsolete since the development of the SGMA method so prefer use StructuralAnalysis_Steinhardt_GMM executable" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename = argv[1];
	string CrystalName = argv[2];
	string OutputFilename = argv[3];
	Bicrystal MySystem(InputFilename, CrystalName); // for the moment need to be a bicrystal object (to have the computeAuxiliary methods), it has no sense but its complicated to change (it should be an AtomicSystem instead but the AtomicSystem cannot have a ComputeAuxiliary)
	MySystem.setAux_vec(MySystem.get_CA()->StructuralAnalysis_Steinhardt(), 2, "Struct");//MySystem.get_CA()->getNumberRefDef_Steinhardt(), "sAB");
	MySystem.printSystem_aux(OutputFilename, "Struct");
	return 0;
}
