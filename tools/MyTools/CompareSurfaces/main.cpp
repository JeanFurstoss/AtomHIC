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
	if( argc < 2 ){
		cerr << "Usage: ./CompareSurfaces InputFilename1 up/down1 InputFilename2 up/down2" << endl;
		cerr << "InputFilename should be an atomic system file" << endl;
		cerr << "Bottom/Upper specify which surface should be duplicated" << endl;
		Dis.Printer_FitAndSaveGMM();	
		return EXIT_FAILURE;
	}

	Dis.Printer_FitAndSaveGMM();	
	string InputFilename1 = argv[1];
	string surf2dup1 = argv[2];
	if( surf2dup1 != "down" && surf2dup1 != "up" ){
		cerr << "Second argument should be either \"down\" or \"up\"" << endl;
		exit(EXIT_FAILURE);
	}
	AtomicSystem AtSys1(InputFilename1);
	
	string InputFilename2 = argv[3];
	string surf2dup2 = argv[4];
	if( surf2dup2 != "down" && surf2dup2 != "up" ){
		cerr << "Second argument should be either \"down\" or \"up\"" << endl;
		exit(EXIT_FAILURE);
	}
	AtomicSystem AtSys2(InputFilename2);

	ComputeAuxiliary CA;
	if( CA.AreSurfacesSame(&AtSys1,surf2dup1,&AtSys2,surf2dup2) ) cout << "Surfaces are the same" << endl;
	else cout << "Surfaces are not the same" << endl; 
	
	Dis.ExecutionTime();
	return 0;
}
