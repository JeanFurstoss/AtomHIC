// AtomHic library files
#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
using namespace std;

int main()
{
	string filename;
	cout << "Atomic file name : " ;
	cin >> filename;
	Bicrystal MySystem(filename, "Forsterite");
	MySystem.get_CA()->Compute_StrainTensor(1);
	MySystem.get_CA()->Compute_StrainTensor_invII();
	MySystem.setAux_vec(MySystem.get_CA()->get_StrainTensor(), 6, "Strain");
	MySystem.setAux(MySystem.get_CA()->get_StrainInvII(), "StrainInvII");
	cout << "output atomic file name : " ;
	cin >> filename;
	MySystem.printSystem_aux(filename, "Strain StrainInvII");
	return 0;
}
