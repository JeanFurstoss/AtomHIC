//**********************************************************************************
//*   MakeNeutralSurfaces/main.cpp                 		                   *
//**********************************************************************************
//* This file contains the implementation of the MakeNeutralSurfaces		   *
//* executable.							                   *
//* It allows to move atom from one surface to another to try make the surfaces as *
//* neutral as possible
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
//*	- work on Bicrystal class to have more clear outputs			   *
//**********************************************************************************

#include <AtomicSystem.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <Displays.h>

using namespace std;


int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc != 3 && argc != 4 ){
		cerr << "Usage: ./MakeNeutralSurfaces InputFilename (CrystalName) OutputFilename" << endl;
		cerr << "This executable tries to move ions from the bottom surface to the upper one to get the surface as neutral as possible" << endl;
		cerr << "The atomic system needs to define charges of the ions" << endl;
		cerr << "The considered surfaces are the one oriented normal to the z axis" << endl;
		cerr << "Providing the crystal name allows to apply the do not separe instructions readed from the crystal database (see /data/ExampleFiles/Crystal.ath for more informations)" << endl;
		return EXIT_FAILURE;
	}
	
	Dis.Printer_OnlyAuxProp();
	
	string InputFilename = argv[1];
	string OutputFilename, CrystalName;
	AtomicSystem MySys(InputFilename);
	if( argc == 3 )	OutputFilename = argv[2];
	else{
		CrystalName = argv[2];
		OutputFilename = argv[3];
		MySys.setCrystal(CrystalName);
		MySys.ComputeNotSepList();
	}
	MySys.MakeSurfaceNeutral();
	MySys.printSystem(OutputFilename);

	Dis.ExecutionTime();	
	return 0;
}
