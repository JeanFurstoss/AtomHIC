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
	Bicrystal MySystem("dump.lmp","Forsterite");
	unsigned int nbAt = MySystem.getNbAtom();
	unsigned int nbAt_type = MySystem.getCrystal()->getNbAtomType();
	double rc = 5.;
	int l_sph = 10;
	MySystem.setAux_vec(MySystem.get_CA()->ComputeSteinhardtParameters(rc,l_sph,"Mono","Mono"),(l_sph+1),"Q_mono_mono");
	MySystem.setAux_vec(MySystem.get_CA()->ComputeSteinhardtParameters(rc,l_sph,"Multi","Multi"),(l_sph+1),"Q_multi_multi");
	MySystem.printSystem_aux("output.xsf","Q_mono_mono Q_multi_multi");
	return 0;
}
