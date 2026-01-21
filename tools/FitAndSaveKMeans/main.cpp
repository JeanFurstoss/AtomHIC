//**********************************************************************************
//*   FitAndSaveKMeans/main.cpp                                                    *
//**********************************************************************************
//* This file contains the implementation of the FitAndSaveKMeans executable.      *
//* It allows to fit and save in the AtomHIC database a KMeans model fitted on     *
//* given descriptors and using the SGMA method	(with GMM replaced by KMeans)	   *
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
#include "Descriptors.h"
#include "KMeans.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc < 6 ){
		cerr << "Usage: ./FitAndSaveKMeans InputDirectory (NameOfDescriptor) n_clust_min n_clust_max Strongly(1)/Softly(0)Labelled NameOfDatabase" << endl;
		cerr << "InputDirectory should contains subdirectories having the names of the labels and containing either: " << endl;
	        cerr << "\t- atomic files (in .cfg format such as printed by AtomHIC) where the name of the auxiliary property to use as descriptor for KMeans has to be provided in the NameOfDescriptor argument, in this case it is recommanded to provide the properties of the descriptor in a \"InputDirectory/DescriptorProperties.ath\" file (an example file could be found in /data/ExampleFiles/), the descriptors can be filtered (by type or element for instance) by setting the value of DESCRIPTORS_FILTERING_TYPE in FixedParameters.ath file in the working directory (if not default value will be used)" << endl;
		cerr << "\t- or files listing only the descriptors" << endl;
		cerr << "n_clust_min and n_clust_max are the range for searching the optimal KMeans number of cluster" << endl;
		cerr << "Softly labelled KMeans (0) means that the KMeans is fitted on all descriptors and labelled afterwards, meaning that eventually some labels migth not correspond to a KMeans cluster" << endl;
		cerr << "Strongly labelled KMeans (1) means that a KMeans is fitted for each label, which are gathered at the end to construct the final KMeans, meaning that each label are described by almost one KMeans cluster" << endl;
		cerr << "Once fitted and labelled, the obtained KMeans parameters will be stored in a database which could then be used to classify data using the ./KMeansClassification executable" << endl;
		Dis.Printer_FitAndSaveKMeans();
		return EXIT_FAILURE;
	}

	Dis.Printer_FitAndSaveKMeans();

	if( argc == 6 ){	
		string InputDir = argv[1];
		istringstream iss_nmin(argv[2]);
		unsigned int nmin;
		iss_nmin >> nmin;
		istringstream iss_nmax(argv[3]);
		unsigned int nmax;
		iss_nmax >> nmax;
		istringstream iss_method(argv[4]);
		unsigned int method;
		iss_method >> method;
		string DatabaseFilename = argv[5];
		Descriptors MyDescriptors(InputDir);
		KMeans MyKM;
		MyKM.setDescriptors(&MyDescriptors);
		if( method == 0 ) MyKM.TrainModel(nmin,nmax,true);
		else MyKM.TrainModel(nmin,nmax,false);
		MyKM.PrintToDatabase(DatabaseFilename);
	}else{
		string InputDir = argv[1];
		string DescriptorName = argv[2];
		istringstream iss_nmin(argv[3]);
		unsigned int nmin;
		iss_nmin >> nmin;
		istringstream iss_nmax(argv[4]);
		unsigned int nmax;
		iss_nmax >> nmax;
		istringstream iss_method(argv[5]);
		unsigned int method;
		iss_method >> method;
		string DatabaseFilename = argv[6];
		Descriptors MyDescriptors(InputDir,DescriptorName);
		KMeans MyKM;
		MyKM.setDescriptors(&MyDescriptors);
		if( method == 0 ) MyKM.TrainModel(nmin,nmax,true);
		else MyKM.TrainModel(nmin,nmax,false);
		MyKM.PrintToDatabase(DatabaseFilename);
	}
	Dis.ExecutionTime();	
	return 0;
}
