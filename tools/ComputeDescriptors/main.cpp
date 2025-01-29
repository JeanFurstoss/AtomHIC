//**********************************************************************************
//*   ComputeDescriptors/main.cpp                                                  *
//**********************************************************************************
//* This file contains the implementation of the ComputeDescriptors executable.    *
//* It allows to compute some atomic descriptors, print them and their properties  *
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
//*	- modify the code when other descriptors are implemented         	   *
//**********************************************************************************

#include <AtomicSystem.h>
#include <SteinhardtDescriptors.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc < 3 ){
		cerr << "Usage: ./ComputeDescriptors AtomicInputFilename NameOfDescriptor OutputFilename" << endl;
		cerr << "The descriptor properties will be read from /data/FixedParameters/FixedParameters.dat" << endl;
		cerr << "In addition to the output atomic file, this executable will generate a DescriptorProperties.ath file containing the descriptors properties and which can be used for fitting a ML for instance" << endl;
		return EXIT_FAILURE;
	}
	
	string InputFilename, OutputFilename, DescriptorName;
	InputFilename = argv[1];
	DescriptorName = argv[2];
	OutputFilename = argv[3];

	AtomicSystem MySystem(InputFilename);

	if( DescriptorName == "Steinhardt" ){
		// Compute the descriptor
		SteinhardtDescriptors MyDescriptors(&MySystem);
		// Set the auxiliary property
		MySystem.setAux_vec(MyDescriptors.getDescriptors(),MyDescriptors.getDim(),DescriptorName);
		// Print descriptor parameters
		ofstream writefile("DescriptorProperties.ath");
		MyDescriptors.printDescriptorsPropToDatabase(writefile);
		writefile.close();
	}else{ // other developped descriptors could be put here
		cerr << "The descriptor name does not correspond to a descriptor that AtomHIC can compute, aborting" << endl;
		return EXIT_FAILURE;
	}

	MySystem.printSystem_aux(OutputFilename,DescriptorName);
	Dis.ExecutionTime();	
	return 0;
}
