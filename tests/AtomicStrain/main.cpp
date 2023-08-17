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
	AtomicSystem ReferenceSystem("AtStrainRef.lmp");
	ComputeAuxiliary CA(&ReferenceSystem);
	AtomicSystem AnalyzedSystem("Sheared.lmp");
	unsigned int nbAt = ReferenceSystem.getNbAtom();
	double *strains = new double[nbAt*8];
	AtomicSystem AtStrainToAddSystem("ToAdd.cfg");
	unsigned int size_s, ind_s;
	ind_s = AtStrainToAddSystem.getAuxIdAndSize("AtomicStrain",size_s);
	for(unsigned int i=0;i<nbAt;i++){
		for(unsigned int j=0;j<8;j++) strains[i*8+j] = AtStrainToAddSystem.getAux(ind_s)[i*8+j];
	}
	double *buffer_strains = CA.Compute_AtomicStrain(AnalyzedSystem,5);
	for(unsigned int i=0;i<nbAt;i++){
		for(unsigned int j=0;j<8;j++) buffer_strains[i*8+j] += strains[i*8+j];
	}
	AnalyzedSystem.setAux_vec(buffer_strains,8,"AtomicStrain");
	AnalyzedSystem.printSystem_aux("output.xsf","AtomicStrain");
	delete[] strains;

	return 0;
}
