// AtomHic library files
#include <AtomicSystem.h>
#include <Crystal.h>
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
	//cout << MySystem.getCrystal()->getA1()[0] << " " << MySystem.getCrystal()->getA1()[1] << " " << MySystem.getCrystal()->getA1()[2] << endl; 
	ComputeAuxiliary CompAux(&MySystem);
	double rc=5.;
	int l;
	cout << "spherical harmonic degree : ";
	cin >> l;	
	MySystem.setAux(CompAux.BondOrientationalParameter(l,rc), "order");
	//unsigned int *nbNeigh = new unsigned int[MySystem.getNbAtom()];
	//for(unsigned int i=0;i<MySystem.getNbAtom();i++) nbNeigh[i] = MySystem.getNeighbours()[i*(MySystem.getNbMaxN()+1)];
	//MySystem.setAux(nbNeigh, "nbNeigh");
	cout << "Output file name : ";
	cin >> filename;
	//MySystem.printSystem_aux(filename,"order nbNeigh");
	MySystem.printSystem_aux(filename,"order");
	//delete[] nbNeigh;
	return 0;
}
