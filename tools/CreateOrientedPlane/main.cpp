//**********************************************************************************
//*   CreateOrientedPlane/main.cpp                                                 *
//**********************************************************************************
//* This file contains the implementation of the CreateOrientedPlane executable.   *
//* It allows to create an orthogonal atomic system containing a perfect crystal   *
//* z-oriented to a given plane							   *
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
//* What is still needed to do here:                                               *
//*	- 									   *
//**********************************************************************************

#include <AtomicSystem.h>
#include <Crystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "ComputeAuxiliary.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc < 5 ){
		cerr << "Usage: CreateOrientedPlane h k l CrystalName(has to be defined in /data/Crystal/) OutputFilename" << endl;
		cerr << "This executable will create an orthogonal cell with a z-oriented (hkl) plane" << endl;
		cerr << "The numerical parameters used for this construction can be tuned in the /data/FixedParameters/FixedParameters.dat file" << endl;
		return EXIT_FAILURE;
	}
	int h, k ,l;
	istringstream iss_h(argv[1]);
	iss_h >> h;
	istringstream iss_k(argv[2]);
	iss_k >> k;
	istringstream iss_l(argv[3]);
	iss_l >> l;
	string crystalType=argv[4];
	string filename=argv[5];
	Crystal MyCrystal(crystalType);
	MyCrystal.RotateCrystal(h,k,l);
	MyCrystal.ConstructOrthogonalCell();
	MyCrystal.getOrientedSystem()->print_lmp(filename);
	Dis.ExecutionTime();	
	return 0;
}
