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
	if( argc < 7 ){
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename GBNormalDirection CrystalType AtomicOutputFilename OrderParameterDensityFilename StatdataFilename" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename = argv[1];
	string NormalDir = argv[2];
	string CrystalType = argv[3];
	string OutputFilename = argv[4];
	string DensityFilename = argv[5];
	string StatdataFilename = argv[6];
	Bicrystal MySystem(InputFilename, NormalDir, CrystalType);
	MySystem.printSystem_aux(OutputFilename, "Disorder");
	MySystem.Print1dDensity(DensityFilename, "Disorder");
	ofstream writefile(StatdataFilename);
	writefile << "GBPosition GBWidth\n";
	writefile << MySystem.getGBPos1() << " " << MySystem.getGBwidth1();
	writefile.close();
	MathTools *MT = new MathTools;
	ofstream writefile_bis(StatdataFilename);
	writefile_bis << "Pos Distrib\n";
	//for(unsigned int i=0;i<100;i++)	writefile_bis << i*MySystem.getH3()[2]/99 << " " << MT->gaussian(i*MySystem.getH3()[2]/99, MySystem.getGBPos1(), MySystem.getGBwidth1()) << "\n";
	for(unsigned int i=0;i<100;i++)	writefile_bis << i*MySystem.getH3()[2]/99 << " " << MT->gaussian(i*MySystem.getH3()[2]/99, MySystem.getGBPos1(), 10.) << "\n";
	writefile_bis.close();
	delete MT;
	return 0;
}
