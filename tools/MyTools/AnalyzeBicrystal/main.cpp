//**********************************************************************************
//*   AnalyzeBicrystal/main.cpp                                                    *
//**********************************************************************************
//* This file contains the implementation of the AnalyzeBicrystal executable.      *
//* It allows to find the position, width and excess volume of a GB	 	   *
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
//*	- more documentation and test and put it in generic tools		   *
//**********************************************************************************

#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();

	// check the number of argument
	if( argc < 8 ){
		cerr << "Usage: AnalyzeBicrystal_ARGS AtomicInputFilename GBNormalDirection CrystalType AtomicOutputFilename OrderParameterDensityFilename GaussianOrderParameterDensityFilename StatdataFilename" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename = argv[1];
	string NormalDir = argv[2];
	string CrystalType = argv[3];
	string OutputFilename = argv[4];
	string DensityFilename = argv[5];
	string GaussDensityFilename = argv[6];
	string StatdataFilename = argv[7];
	Bicrystal MySystem(InputFilename, NormalDir, CrystalType);
	MySystem.printSystem_aux(OutputFilename, "Disorder");
	MySystem.Print1dDensity(DensityFilename, "GBProfile");
	MySystem.Print1dDensity(GaussDensityFilename, "GBProfile_Gauss");
	ofstream writefile(StatdataFilename);
	writefile << MySystem.getGBPos1() << " " << MySystem.getGBwidth1() << " " << MySystem.getExcessVol();
	writefile.close();
	Dis.ExecutionTime();	
	return 0;
}
