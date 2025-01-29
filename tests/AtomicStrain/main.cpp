//**********************************************************************************
//*   AtomicStrain/main.cpp                                                        *
//**********************************************************************************
//* This file contains the implementation of the AtomicStrain test.                *
//**********************************************************************************
//* (C) Jan 2025 - Jean Furstoss                                                   *
//*     Université de Poitiers, Institut PPRIME                                    *
//*     UPR CNRS 3346, 86360 Chasseuneuil-du-Poitou, France                        *
//*     jean.furstoss@univ-poitiers.fr                                             *
//* Last modification: J. Furstoss - 28 Janv 2025                                  *
//**********************************************************************************
//* This program is free software: you can redistribute it and/or modify           *
//* it under the terms of the GNU General Public License as published by           *
//* the Free Software Foundation, either version 3 of the License, or              *
//* (at your option) any later version.                                            *
//*                                                                                *
//* This program is distributed in the hope that it will be useful,                *
//* but WITHOUT ANY WARRANTY; without even the implied warranty of                 *
//* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  *
//* GNU General Public License for more details.                                   *
//*                                                                                *
//* You should have received a copy of the GNU General Public License              *
//* along with this program.  If not, see <http://www.gnu.org/licenses/>.          *
//**********************************************************************************

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
