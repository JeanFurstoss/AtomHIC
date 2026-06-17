//**********************************************************************************
//*   FitAndSaveGMM/main.cpp                                                       *
//**********************************************************************************
//* This file contains the implementation of the FitAndSaveGMM executable.         *
//* It allows to fit and save in the AtomHIC database a gaussian mixture model	   *
//* fitted on given descriptors and using the SGMA method			   *
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
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include "AtomicSystem.h"
#include "MathTools.h"
#include "Crystal.h"
#include "ComputeAuxiliary.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc != 5 && argc != 6 ){
		cerr << "Usage: ./SymmetrizeSurfaces InputFilename down/up z_cut (CrystalName) OutputFilename" << endl;
		cerr << "InputFilename should be an atomic system file" << endl;
		cerr << "down/up specify which surface should be duplicated" << endl;
		cerr << "zcut is the altitude at which to cut the system before pasting the surface" << endl;
		cerr << "if the crystal name is provided the stoichiometry will be verified before printing the system" << endl;
		Dis.Printer_FitAndSaveGMM();	
		return EXIT_FAILURE;
	}

	Dis.Printer_FitAndSaveGMM();	
	string InputFilename = argv[1];
	string surf2dup = argv[2];
	if( surf2dup != "down" && surf2dup != "up" ){
		cerr << "Second argument should be either \"down\" or \"up\"" << endl;
		exit(EXIT_FAILURE);
	}
	istringstream iss_zcut(argv[3]);
	double zcut;
	iss_zcut >> zcut;
	string CrystalName, OutputFilename;
	AtomicSystem AtSys(InputFilename);
	if( argc == 6 ){
		CrystalName = argv[4];
		OutputFilename = argv[5];
		AtSys.setCrystal(CrystalName);
	}else
		OutputFilename = argv[4];


	if( AtSys.SymmetrizeRelaxedSurfaces(surf2dup,zcut) ){
		if( ( argc == 6 && AtSys.IsSystemStoichiometric() ) || argc == 5 ) AtSys.printSystem(OutputFilename);
		else cout << "The resulting system is not stoichiometric, aborting" << endl;
	}else cout << "Cannot symmetrize relaxed surface" << endl;
	Dis.ExecutionTime();
	return 0;
}
