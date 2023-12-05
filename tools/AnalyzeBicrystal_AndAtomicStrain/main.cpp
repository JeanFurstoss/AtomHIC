// AtomHic library files
#include <Bicrystal.h>
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
	// check the number of argument
	if( argc < 10 ){
		cerr << "Usage: AnalyzeBicrystal_AndAtomicStrain AtomicInputFilename GBNormalDirection CrystalType ReferenceAtomicStrainFilename cutoff_radius(AtomicStrain) AtomicOutputFilename OrderParameterDensityFilename GaussianOrderParameterDensityFilename StatdataFilename" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename = argv[1];
	string NormalDir = argv[2];
	string CrystalType = argv[3];
	string RefFilename = argv[4];
	double rc; 
	istringstream iss_rc(argv[5]);
	iss_rc >> rc;
	string OutputFilename = argv[6];
	string DensityFilename = argv[7];
	string GaussDensityFilename = argv[8];
	string StatdataFilename = argv[9];
	Bicrystal MySystem(InputFilename, NormalDir, CrystalType);

	AtomicSystem ReferenceSystem(RefFilename);
	ComputeAuxiliary CA(&ReferenceSystem);

	unsigned int nbAt = ReferenceSystem.getNbAtom();
	
	MySystem.setAux_vec(CA.Compute_AtomicStrain(MySystem,rc),8,"AtomicStrain");

	MySystem.printSystem_aux(OutputFilename, "Disorder AtomicStrain");
	//MySystem.Print1dDensity(DensityFilename, "GBProfile");
	MySystem.Print1dDensity(DensityFilename, "Disorder");
	MySystem.Print1dDensity(GaussDensityFilename, "GBProfile_Gauss");
	ofstream writefile(StatdataFilename);
	writefile << MySystem.getGBPos1() << " " << MySystem.getGBwidth1() << " " << MySystem.getExcessVol();
	writefile.close();
	return 0;
}
