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
	if( argc < 4 ){
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename l rcut AtomicOutputFilename" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename = argv[1];
	istringstream iss_l(argv[2]);
	int lsph;
	iss_l >> lsph;
	istringstream iss_r(argv[3]);
	double rcut;
	iss_r >> rcut;
	string OutputFilename = argv[4];
	Bicrystal MySystem(InputFilename); // for the moment need to be a bicrystal object (to have the computeAuxiliary methods), it has no sense but its complicated to change (it should be an AtomicSystem instead but the AtomicSystem cannot have a ComputeAuxiliary)
	MySystem.setAux_vec(MySystem.get_CA()->ComputeSteinhardtParameters_Mono(rcut,lsph), ((unsigned int) lsph + 1), "SteinhardtParams");
	MySystem.printSystem_aux(OutputFilename, "SteinhardtParams");
	return 0;
}
