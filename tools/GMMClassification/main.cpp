//**********************************************************************************
//*   GMMClassification/main.cpp                                                   *
//**********************************************************************************
//* This file contains the implementation of the GMMClassification executable.     *
//* It allows to classify unknown descriptor using a GMM saved in the AtomHIC	   *
//* machine-learning database (used for SGMA method)				   *
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
//*	- 									   *
//**********************************************************************************

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "Descriptors.h"
#include "SteinhardtDescriptors.h"
#include "ACEDescriptors.h"
#include "GaussianMixtureModel.h"
#include "AtomicSystem.h"
#include "AtomicSystemTrajectory.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc < 5 ){
		cerr << "Usage: ./GMMClassification InputFilename (NameOfDescriptor) NameOfDatabase OutputFilename OutputDescriptorsOrNot" << endl;
		cerr << "InputFilename could be an atomic system containing or not the descriptor values. If the atomic system contains the descriptors values the name of the auxiliary properties corresponding to the descriptors could be provided as the NameOfDescriptor argument" << endl;
		cerr << "InputFilename could also be a simple file containing the descriptors value" << endl;
		cerr << "The argument NameOfDatabse should be a directory in /data/MachineLearningModels/GaussianMixtureModel/ which could be generated using the ./FitAndSaveGMM executable" << endl;
		cerr << "The obtained output file contains the index of the labels (Struct[1]) and the maximum likelihood classifier (Struct[2])" << endl;
		cerr << "If OutputDescriptorsOrNot = 0 the descriptors wont be printed in the output file, if OutputDescriptorsOrNot = 1 the descriptors will be printed in addition with Struct[1] and Struct[2]" << endl;
		GaussianMixtureModel GMM;
		vector<string> basis = GMM.getAvailableDatabases();
		cerr << "Available GMM databases:" << endl;
		for(unsigned int i=0;i<basis.size();i++){
			cerr << basis[i] << endl;
		}
		return EXIT_FAILURE;
	}


	// read parameters
	unsigned int current_read_ind = 1;
	string InputFilename = argv[current_read_ind];
	current_read_ind++;
	
	string DescriptorName;
	if( argc == 6 ){
		DescriptorName = argv[current_read_ind];
		current_read_ind++;
	}

	string DatabaseFilename = argv[current_read_ind];
	current_read_ind++;

	string OutputFilename = argv[current_read_ind];
	current_read_ind++;

	unsigned int outdes_i;
	istringstream iss_o(argv[current_read_ind]);
	iss_o >> outdes_i;
	bool outdes = false;
	if( outdes_i == 1 ) outdes = true;

	// set system, nbSys = 0 => simple file (only descriptors), 1 => single atomic file, >1 atomic trajectory
	AtomicSystemTrajectory MyTraj;
	AtomicSystem *MySystem;

	unsigned int nbSys = 0;

	if( MyTraj.SearchIsTrajectory(InputFilename) ){
		MyTraj.setAtomicSystemList(InputFilename);
		nbSys = MyTraj.getNbSys();
	}else{
		MySystem = new AtomicSystem();
		if( MySystem->FilenameConstructor(InputFilename) ) nbSys = 1;
		else{
			delete MySystem;
			MySystem = nullptr;
		}
	}

	if( nbSys == 0 ) Dis.Printer_NoFixedParams();
	else Dis.Printer_OnlyAuxProp();

	// read GMM properties
	GaussianMixtureModel GMM;
	GMM.ReadModelParamFromDatabase(DatabaseFilename);
	string DescriptorNameGMM;
	string ftype = GMM.getFilteringType();
	if( nbSys > 0 && argc == 5 ){	
		size_t pos_name;
		string buffer_s;
		for(unsigned int s=0;s<GMM.getDescriptorProperties().size();s++){
			pos_name = GMM.getDescriptorProperties()[s].find("DESCRIPTOR_NAME");
			if( pos_name!=string::npos ){
				istringstream text(GMM.getDescriptorProperties()[s]);
				text >> buffer_s >> DescriptorNameGMM;
			}
		}
	}
	
	if( argc == 6 ){
		cout << endl;
		cout << "WARNING !!! From the database we expect that the provided descriptors (" << DescriptorName << ") are computed using the following properties: " << endl;
		for(unsigned int s=0;s<GMM.getDescriptorProperties().size();s++) cout << GMM.getDescriptorProperties()[s] << endl;
		cout << endl;
	}

	// set descriptors if file is not atomic file
	Descriptors *MyDescriptors;
	if( nbSys == 0 ){
		MyDescriptors = new Descriptors(InputFilename);
		nbSys = 1;
	}

	// loop over number of system
	for(unsigned int i=0;i<nbSys;i++){
		// set AtomicSystem in case of trajectory file
		if( nbSys > 1 ){
			cout << "Treating timestep " << MyTraj.getTimestep(i) << endl;
			MySystem = MyTraj.getAtomicSystem(i);
		}
		
		// compute descriptors for atomic files if the descriptor is not provided
		if( MySystem ){
			if( argc == 5 ){
				if( DescriptorNameGMM == "Steinhardt" )
					MyDescriptors = new SteinhardtDescriptors(MySystem,GMM.getDescriptorProperties());
				else if( DescriptorNameGMM == "ACE" )
					MyDescriptors = new ACEDescriptors(MySystem,GMM.getDescriptorProperties());
				else{ // other developped descriptors could be put here
					cerr << "The descriptor name does not correspond to a descriptor that AtomHIC can compute, aborting" << endl;
					return EXIT_FAILURE;
				}
			}else{
				MyDescriptors = new Descriptors(MySystem,DescriptorName,ftype);
				DescriptorNameGMM = DescriptorName;
			}
		}

		// perform classification
		GMM.setDescriptors(MyDescriptors);
		GMM.Classify();

		// set classificator and descriptors if asked
		if( MySystem ){
			MySystem->setAux_vec(GMM.getClassificator(),2,"Struct");
			if( outdes && argc == 5 ) MySystem->setAux_vec(MyDescriptors->getDescriptors(),MyDescriptors->getDim(),DescriptorNameGMM);
			// free descriptor memory
			delete MyDescriptors;
		}
	}

	if( MySystem ){
		string ForPrinting = "Struct ";
		if( outdes ) ForPrinting += DescriptorNameGMM;
		if( nbSys == 1 ){
			MySystem->printSystem_aux(OutputFilename,ForPrinting);
			delete MySystem;
		}else
			MyTraj.printSystem_aux(OutputFilename,ForPrinting);
	}else
		GMM.PrintClassifiedData(OutputFilename,outdes);

	Dis.ExecutionTime();	
	return 0;

}
