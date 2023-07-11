// AtomHic library files
#include <AtomicSystem.h>
#include <ComputeAuxiliary.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"

using namespace std;

int main(int argc, char *argv[])
{
	string InputFilename, CrystalType, outfilename;
	if( argc == 4 ){
		InputFilename = argv[1];
		CrystalType = argv[2];
		outfilename = argv[3];
	}else if( argc == 3 ){
		InputFilename = argv[1];
		outfilename = argv[2];
	}else{
		cerr << "Usage: ./BondOriParam AtomicInputFilename CrystalType(The crystal type has to be defined in /data/Crystal/) outputFilename" << endl;
		cerr << "or: ./BondOriParam AtomicInputFilename outputFilename" << endl;
		cerr << "In the latter case, the crystal will be considered as a monosite crystal and the value of spherical harmonic degree and cutoff radius will be read in (/data/FixedParameters/FixedParameters.dat)" << endl; 
		cerr << "The first case should be used for multisite crystals for which the Steinhardt parameters have been stored in the AtomHic database using the SaveSteinhardtToDatabase_PerfectCrystal executable" << endl;
		return EXIT_FAILURE;
	}
	AtomicSystem MySystem(InputFilename);
	if( argc == 4 ) MySystem.setCrystal(CrystalType);
	ComputeAuxiliary CA(&MySystem);
	MySystem.setAux(CA.BondOrientationalParameter(), "Disorder");
	MySystem.printSystem_aux(outfilename, "Disorder");
	return 0;
}
