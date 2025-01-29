//**********************************************************************************
//*   Tool.cpp		                           				   *
//**********************************************************************************
//* This file contains example of utilization of main functions of the AtomHIC     *
//* library which can be used to construct a tailored executable.	   	   *							   *
//* For this, just create a directory in /tools/MyTools/, copy and modify a        *
//* CMakeLists.txt from another tool and write main.cpp code. Then recompile the   *
//* library and the produced executable will be present in /build/bin/MyTools/	   *
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
//*	- more documentation							   *
//**********************************************************************************

// AtomHIC includes
#include <AtomicSystem.h>
#include <Bicrystal.h>
#include <Crystal.h>
#include <ComputeAuxiliary.h>
#include <MathTools.h>
#include "MathTools.h"
#include "MyStructs.h"
#include "Descriptors.h"
#include "GaussianMixtureModel.h"
#include <Displays.h>
// other includes
#include <iostream>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <filesystem>
#include <dirent.h>

using namespace std;

int main(int argc, char *argv[])
{
	// Create a display object, this will store the current time to display the execution time using the ExecutionTime() method
	Displays Dis;
	// Display AtomHIC logo
	Dis.Logo();

	// Declare argument of the executable
	string InputFilename, OutputFilename;
	double rcut;

	// Verify if the provided number of argument is ok, if not print the description of the executable and the expected arguments
	if( argc == 6 ){
		// get the arguments
		InputFilename = argv[1];
		istringstream iss_rc(argv[2]);
		iss_rc >> rcut;
		OutputFilename = argv[3];
	}else{
		cerr << "Usage: ./MyTool InputFilename rcut OutputFilename" << endl;
		cerr << "TODO description" << endl;
		return EXIT_FAILURE;
	}
	
	// TODO

	// Display execution time
	Dis.ExecutionTime();	
	return 0;
}
