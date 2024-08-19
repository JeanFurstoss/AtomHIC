// AtomHic library files
#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include <SteinhardtDescriptors.h>

using namespace std;

int main(int argc, char *argv[])
{
	AtomicSystem MySystem("dump.lmp");
	double rc = 5.;
	int l_sph = 10;
	vector<string> Properties_1;
	Properties_1.push_back("STEINHARDT_MODE Full");
	Properties_1.push_back("CUTOFF_RADIUS 5.");
	Properties_1.push_back("STEINHARDT_STYLE Mono");
	Properties_1.push_back("AVE_STYLE Mono");
	Properties_1.push_back("NUMBER_OF_DIMENSION 10");
	vector<string> Properties_2;
	Properties_2.push_back("STEINHARDT_MODE Full");
	Properties_2.push_back("CUTOFF_RADIUS 5.");
	Properties_2.push_back("STEINHARDT_STYLE Multi");
	Properties_2.push_back("AVE_STYLE Multi");
	Properties_2.push_back("NUMBER_OF_DIMENSION 10");
	SteinhardtDescriptors St_1(&MySystem,Properties_1);
	SteinhardtDescriptors St_2(&MySystem,Properties_2);
	MySystem.setAux_vec(St_1.getDescriptors(),10,"Q_mono_mono");
	MySystem.setAux_vec(St_2.getDescriptors(),10,"Q_multi_multi");
	MySystem.printSystem_aux("output.xsf","Q_mono_mono Q_multi_multi");
	return 0;
}
