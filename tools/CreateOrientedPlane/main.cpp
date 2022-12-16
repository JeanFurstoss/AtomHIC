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
	int h, k ,l;
	istringstream iss_h(argv[1]);
	iss_h >> h;
	istringstream iss_k(argv[2]);
	iss_k >> k;
	istringstream iss_l(argv[3]);
	iss_l >> l;
	string filename="Forsterite";
	Crystal MyCrystal(filename);
	MyCrystal.RotateCrystal(h,k,l);
	MyCrystal.ConstructOrthogonalCell();
	MyCrystal.getOrientedSystem()->print_lmp("test.lmp");
	return 0;
}
