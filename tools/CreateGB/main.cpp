// AtomHic library files
#include <AtomicSystem.h>
#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "ComputeAuxiliary.h"

using namespace std;


int main(int argc, char *argv[])
{
	if( argc < 8 ){
		cerr << "Usage: CreateGB h_RotAxis k_RotAxis l_RotAxis RotAngle h_GBPlane k_GBPlane l_GBPlane CrystalName(has to be defined in /data/Crystal/)" << endl;
		return EXIT_FAILURE;
	}
	int h_a, k_a ,l_a, h_p, k_p, l_p;
	double theta;
	string crystalName;
	istringstream iss_ha(argv[1]);
	iss_ha >> h_a;
	istringstream iss_ka(argv[2]);
	iss_ka >> k_a;
	istringstream iss_la(argv[3]);
	iss_la >> l_a;
	istringstream iss_theta(argv[4]);
	iss_theta >> theta;
	istringstream iss_hp(argv[5]);
	iss_hp >> h_p;
	istringstream iss_kp(argv[6]);
	iss_kp >> k_p;
	istringstream iss_lp(argv[7]);
	iss_lp >> l_p;
	istringstream iss_cn(argv[8]);
	iss_cn >> crystalName;
	//Bicrystal MyGB(crystalName,h_a,k_a,l_a,theta,h_p,k_p,l_p);
	Bicrystal MyGB(crystalName,h_a,k_a,l_a,theta,h_p,k_p,l_p,false);
	MyGB.print_lmp("GB.lmp");
	MyGB.printCSL("CSL.lmp");
	MyGB.print_Grains();
	return 0;
}
