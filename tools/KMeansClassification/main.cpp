//**********************************************************************************
//*   KMeansClassification/main.cpp                                                *
//**********************************************************************************
//* This file contains the implementation of the KMeansClassification executable.  *
//* It allows to classify unknown descriptor using a KMeans saved in the AtomHIC   *
//* machine-learning database (used for SGMA method with GMM replaced by KMeans)   *
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
#include "KMeans.h"
#include "AtomicSystem.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc < 4 ){
		cerr << "Usage: ./KMeansClassification InputFilename (NameOfDescriptor) NameOfDatabase OutputFilename" << endl;
		cerr << "InputFilename could be an atomic system or a file directly containing the descriptor values, in the latter case NameOfDescriptor argument should not be used" << endl;
		cerr << "The argument NameOfDatabse should be a directory in /data/MachineLearningModels/KMeans/ which could be generated using the ./FitAndSaveKMeans executable" << endl;
		cerr << "In the case of an atomic system the obtained output file contains, the computed descriptors if the descriptors are not provided, the index of the labels (Struct[1])" << endl;
		cerr << "In the case of a simple file, the descriptors are simply printed again as well as Struct[1]" << endl;
		KMeans MyKM;
		vector<string> basis = MyKM.getAvailableDatabases();
		cerr << "Available KMeans databases:" << endl;
		for(unsigned int i=0;i<basis.size();i++){
			cerr << basis[i] << endl;
		}
		return EXIT_FAILURE;
	}
	
	Dis.Printer_NoFixedParams();
	
	if( argc == 4 ){	
		string InputFilename = argv[1];
		string DatabaseFilename = argv[2];
		string OutputFilename = argv[3];
		KMeans MyKM;
		MyKM.ReadModelParamFromDatabase(DatabaseFilename);
		AtomicSystem MySystem;
		if( MySystem.FilenameConstructor(InputFilename) ){ // in this case, the provided file is an atomic file without descriptors, thus compute the descriptors needed from the database
			string DescriptorName, buffer_s;
			size_t pos_name;
			for(unsigned int s=0;s<MyKM.getDescriptorProperties().size();s++){
				pos_name = MyKM.getDescriptorProperties()[s].find("DESCRIPTOR_NAME");
				if( pos_name!=string::npos ){
					istringstream text(MyKM.getDescriptorProperties()[s]);
					text >> buffer_s >> DescriptorName;
				}
			}
			if( DescriptorName == "Steinhardt" ){
				SteinhardtDescriptors MyDescriptors(&MySystem,MyKM.getDescriptorProperties());
				MyKM.setDescriptors(&MyDescriptors);
				MyKM.Classify();
				MySystem.setAux_vec(MyDescriptors.getDescriptors(),MyDescriptors.getDim(),"Steinhardt");
				MySystem.setAux_vec(MyKM.getClassificator(),2,"Struct");
				MySystem.printSystem_aux(OutputFilename,"Steinhardt Struct");
				return 0;
			}else{ // other developped descriptors could be put here
				cerr << "The descriptor name does not correspond to a descriptor that AtomHIC can compute, aborting" << endl;
				return EXIT_FAILURE;
			}
		}else{
			Descriptors MyDescriptors(InputFilename);
			MyKM.setDescriptors(&MyDescriptors);
			MyKM.Classify();
			MyKM.PrintClassifiedData(OutputFilename);
		}
	}else{
		string InputFilename = argv[1];
		string DescriptorName = argv[2];
		string DatabaseFilename = argv[3];
		string OutputFilename = argv[4];
		KMeans MyKM;
		MyKM.ReadModelParamFromDatabase(DatabaseFilename);
		string ftype = MyKM.getFilteringType();
		cout << endl;
		cout << "WARNING !!! From the database we expect that the provided descriptors (" << DescriptorName << ") are computed using the following properties: " << endl;
		for(unsigned int s=0;s<MyKM.getDescriptorProperties().size();s++) cout << MyKM.getDescriptorProperties()[s] << endl;
		cout << endl;
		AtomicSystem MySystem(InputFilename);
		Descriptors MyDescriptors(&MySystem,DescriptorName,ftype);
		MyKM.setDescriptors(&MyDescriptors);
		MyKM.Classify();
		MySystem.setAux_vec(MyKM.getClassificator(),2,"Struct");
		MySystem.printSystem_aux(OutputFilename,"Struct");
	}
	Dis.ExecutionTime();	
	return 0;
}
