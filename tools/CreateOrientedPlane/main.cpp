// AtomHic library files
#include <AtomicSystem.h>
#include <Crystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "ComputeAuxiliary.h"
using namespace std;

int main(int argc, char *argv[])
{
	if( argc < 5 ){
		cerr << "Usage: CreateOrientedPlane h k l CrystalName(has to be defined in /data/Crystal/) OutputFilename" << endl;
		cerr << "This executable will create an orthogonal cell with a z-oriented (hkl) plane" << endl;
		cerr << "The numerical parameters used for this construction can be tuned in the /data/FixedParameters/FixedParameters.dat file" << endl;
		return EXIT_FAILURE;
	}
	int h, k ,l;
	istringstream iss_h(argv[1]);
	iss_h >> h;
	istringstream iss_k(argv[2]);
	iss_k >> k;
	istringstream iss_l(argv[3]);
	iss_l >> l;
	string crystalType=argv[4];
	string filename=argv[5];
	Crystal MyCrystal(crystalType);
	MyCrystal.RotateCrystal(h,k,l);
	MyCrystal.ConstructOrthogonalCell();
	MyCrystal.getOrientedSystem()->print_lmp(filename);
	return 0;
}
