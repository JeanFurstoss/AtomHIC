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
	Bicrystal MySystem(filename, "z", "Forsterite");
	MySystem.get_CA()->Compute_StrainTensor(1);
	//MySystem.setAux(MySystem.get_CA()->get_AtomSiteIndex(), "AtomSite");
	cout << "output atomic file name : " ;
	cin >> filename;
	MySystem.printSystem_aux(filename, "Disorder");
	//MySystem.printSystem_aux(filename, "Disorder AtomSite");
	//cout << "output GB profile file name : " ;
	//cin >> filename;
	//MySystem.Print1dDensity(filename, "GBProfile");
	//cout << "output estimated GB profile file name : " ;
	//cin >> filename;
	//MySystem.Print1dDensity(filename, "GBProfile_Gauss");
	//cout << "output density profile : " ;
	//cin >> filename;
	//MySystem.Print1dDensity(filename, "Mass");
	return 0;
}
