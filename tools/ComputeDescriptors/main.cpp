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

#include <AtomicSystemTrajectory.h>
#include <AtomicSystem.h>
#include <SteinhardtDescriptors.h>
#include <ACEDescriptors.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <chrono>
#include <Displays.h>
#include <vector>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc < 4 ){
		cerr << "Usage: ./ComputeDescriptors AtomicInputFilename NameOfDescriptor (yaml_input_filename) OutputFilename" << endl;
		cerr << "The descriptor properties will be read from /data/FixedParameters/FixedParameters.dat" << endl;
		cerr << "In addition to the output atomic file, this executable will generate a DescriptorProperties.ath file containing the descriptors properties and which can be used for fitting a ML for instance" << endl;
		cerr << "The yaml_input_filename shoudl be used only for ACE descriptors (it could be the name of the yaml file if it is in the working directory or the relative or absolute path to it if the file is not in the working directory)" << endl;
		cerr << "Available descriptors : " << endl;
		cerr << "\t - Steinhardt" << endl;
		cerr << "\t - ACE (uses PACE implementation Lysogorskiy et al. 2021)" << endl;
		return EXIT_FAILURE;
	}
	
	string InputFilename, OutputFilename, DescriptorName, yaml_input_filename;
	InputFilename = argv[1];
	DescriptorName = argv[2];
	if( argc == 4 ){
		OutputFilename = argv[3];
	}else{
		yaml_input_filename = argv[3];
		OutputFilename = argv[4];
	}

	if( DescriptorName == "Steinhardt" ) Dis.Printer_ComputeDescriptors_Steinhardt();
	else if( DescriptorName == "ACE" ) Dis.Printer_OnlyAuxProp();

	Descriptors *MyDescriptors;

	AtomicSystemTrajectory MyTraj;
	AtomicSystem *MySystem;
	
	unsigned int nbSys = 1;
	if( MyTraj.SearchIsTrajectory(InputFilename) ){
		MyTraj.setAtomicSystemList(InputFilename);
		nbSys = MyTraj.getNbSys();
	}

	for(unsigned int i=0;i<nbSys;i++){
		if( nbSys == 1 ) MySystem = new AtomicSystem(InputFilename);
		else{
			cout << "Treating timestep " << MyTraj.getTimestep(i) << endl;
			MySystem = MyTraj.getAtomicSystem(i);
		}

		if( DescriptorName == "Steinhardt" )
			MyDescriptors = new SteinhardtDescriptors(MySystem);
		else if( DescriptorName == "ACE" )
			MyDescriptors = new ACEDescriptors(MySystem,yaml_input_filename);
		else{ // other developped descriptors could be put here
			delete MyDescriptors;
			if( nbSys == 1 ) delete MySystem;
			cerr << "The descriptor name does not correspond to a descriptor that AtomHIC can compute, aborting" << endl;
			return EXIT_FAILURE;
		}
		// Set the auxiliary property
		MySystem->setAux_vec(MyDescriptors->getDescriptors(),MyDescriptors->getDim(),DescriptorName);
		cout << endl;

		// free memory
		if( i != nbSys-1 ) delete MyDescriptors;
	}

	// Print descriptor parameters
	ofstream writefile("DescriptorProperties.ath");
	MyDescriptors->printDescriptorsPropToDatabase(writefile);
	writefile.close();
	cout << "DescriptorProperties.ath file successfully writted" << endl;
	
	if( nbSys == 1 ){
		MySystem->printSystem_aux(OutputFilename,DescriptorName);
		delete MySystem;
	}else
		MyTraj.printSystem_aux(OutputFilename,DescriptorName);

	// free memory
	delete MyDescriptors;

	Dis.ExecutionTime();	
	return 0;
}
