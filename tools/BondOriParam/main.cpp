// AtomHic library files
#include <AtomicSystem.h>
#include <ComputeAuxiliary.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

int main(int argc, char *argv[])
{
	cout << "Calculation running using " << omp_get_max_threads() << " threads" << endl;
	
	auto start = high_resolution_clock::now();

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
		cerr << "The first case should be used for multisite and non-centrosymmetric crystals for which the reference bond orientational parameters have been stored in the AtomHic database using the ./SaveNonCSCrystalBondOriParam executable" << endl;
		cerr << "More details on the computation of this order parameter can be found in Furstoss et al. (2024) Comp. Mat. Science." << endl;
		return EXIT_FAILURE;
	}
	AtomicSystem MySystem(InputFilename);
	if( argc == 4 ) MySystem.setCrystal(CrystalType);
	ComputeAuxiliary CA(&MySystem);
	MySystem.setAux(CA.BondOrientationalParameter(), "BondOriParam");
	MySystem.printSystem_aux(outfilename, "BondOriParam");
	

	auto end = high_resolution_clock::now();	
	auto duration = duration_cast<microseconds>(end - start);
	cout << "Resolution time = " << duration.count() << endl;
	return 0;
}
