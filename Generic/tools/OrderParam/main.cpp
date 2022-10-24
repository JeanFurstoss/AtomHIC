// AtomHic library files
#include <AtomicSystem.h>
#include <Crystal.h>
#include <OrthorhombicCrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "ComputeAuxiliary.h"
using namespace std;

int main()
{
	string filename;
	cout << "Atomic file name : " ;
	cin >> filename;
	AtomicSystem MySystem(filename);
	//MySystem.setCrystal("Periclase");
	MySystem.setCrystal("Forsterite");
	ComputeAuxiliary CompAux(&MySystem);
	double rc=5.;
	int l;
	cout << "spherical harmonic degree : ";
	cin >> l;	
	MySystem.setAux(CompAux.BondOrientationalParameter(l,rc), "order");
	cout << "Output file name : ";
	cin >> filename;
	MySystem.printSystem_aux(filename,"order");
	return 0;
}
