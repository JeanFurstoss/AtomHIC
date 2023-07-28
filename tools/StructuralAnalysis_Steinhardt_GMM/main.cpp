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
	if( argc != 4 && argc != 5 ){
		cerr << "Usage: StructuralAnalysis_Steinhardt_GMM AtomicInputFilename CrystalName AtomicOutputFilename" << endl;
		cerr << "or: StructuralAnalysis_Steinhardt_GMM AtomicInputFilename CrystalName AtomicOutputFilename NameOfAuxiliary_SteinhardtParam" << endl;
		cerr << "in the first case, the Steinhardt parameters will be explicitely computed while in the later case, the Steinhardt parameters provided in the InputFilename will be used" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename = argv[1];
	string CrystalName = argv[2];
	string OutputFilename = argv[3];
	string aux_name;
	if( argc == 5 ) aux_name = argv[4];
	else aux_name = "none";
	Bicrystal MySystem(InputFilename, CrystalName); // for the moment need to be a bicrystal object (to have the computeAuxiliary methods), it has no sense but its complicated to change (it should be an AtomicSystem instead but the AtomicSystem cannot have a ComputeAuxiliary)
	MySystem.setAux_vec(MySystem.get_CA()->StructuralAnalysis_Steinhardt_GMM(aux_name), 2, "Struct");//MySystem.get_CA()->getNumberRefDef_Steinhardt(), "sAB");
	MySystem.printSystem_aux(OutputFilename, "Struct");
	return 0;
}
