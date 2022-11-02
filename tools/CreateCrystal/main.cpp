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
	string filename="Forsterite";
	Crystal MyCrystal(filename);
	cout << MyCrystal.getA1()[0] << " " << MyCrystal.getA1()[1] << " " << MyCrystal.getA1()[2] << endl; 
	cout << "Atomic file name : " ;
	cin >> filename;
	AtomicSystem MySystem(filename);
	double rc=5.;
	MySystem.searchNeighbours(rc);
	unsigned int nbMaxNeigh = MySystem.getNbMaxN();
	unsigned int nbAtom = MySystem.getNbAtom();
	unsigned int* Neigh=MySystem.getNeighbours();
	int* CLNeigh=MySystem.getCLNeighbours();
	unsigned int AtSampledId = 86;
	double *Aux = new double[nbAtom];
//	for(unsigned int i=0;i<Neigh[AtSampledId*(nbMaxNeigh+1)];i++) cout << Neigh[AtSampledId*(nbMaxNeigh+1)+i+1] << " " << CLNeigh[AtSampledId*(nbMaxNeigh*3)+i*3] << " " << CLNeigh[AtSampledId*(nbMaxNeigh*3)+i*3+1] << " " << CLNeigh[AtSampledId*(nbMaxNeigh*3)+i*3+2] << endl;
	for(unsigned int i=0;i<nbAtom;i++){
		for(unsigned int j=0;j<Neigh[AtSampledId*(nbMaxNeigh+1)];j++){
			if( i == Neigh[AtSampledId*(nbMaxNeigh+1)+j+1] ){
				Aux[i] = 1;
				break;
			}else{
				Aux[i] = 0;
			}
		}
	}
	MySystem.setAux(Aux, "neighbourOfAt");
	MySystem.setAux(Aux, "neighbourOfAt_bis");
	string aux="neighbourOfAt neighbourOfAt_bis";
	cout << "Output file name : ";
	cin >> filename;
	//MySystem.printSystem(filename);
	MySystem.printSystem_aux(filename,aux);
	delete[] Aux;
	return 0;
}
