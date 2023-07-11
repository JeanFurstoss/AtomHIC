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
	AtomicSystem MySystem("dump.xsf");
	MySystem.setCrystal("Forsterite");
	ComputeAuxiliary CA(&MySystem);
	MySystem.setAux(CA.BondOrientationalParameter(), "Disorder");
	MySystem.printSystem_aux("output.xsf", "Disorder");
	return 0;
}
