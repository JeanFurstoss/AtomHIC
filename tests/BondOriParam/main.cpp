// AtomHic library files
#include <AtomicSystem.h>
#include <ComputeAuxiliary.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

int main(int argc, char *argv[])
{
	AtomicSystem MySystem("dump.xsf");
	MySystem.setCrystal("Forsterite");
	ComputeAuxiliary CA(&MySystem);
	MySystem.setAux(CA.BondOrientationalParameter(), "BondOriParam");
	MySystem.printSystem_aux("output.xsf", "BondOriParam");
	return 0;
}
