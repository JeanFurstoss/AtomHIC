//**********************************************************************************
//*   DBScanClustering/main.cpp                                                    *
//**********************************************************************************
//* This file contains the implementation of the DBScanClustering executable.      *
//* It allows to make a cluster analysis of atomic position using the density	   *
//* based clustering algorythm							   *
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
//*	- add possibility of filtering ions					   *
//*	- generalize to the clustering of descriptors and not only on position	   *
//**********************************************************************************

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "Descriptors.h"
#include "DBScan.h"
#include "AtomicSystem.h"
#include <AtomicSystemTrajectory.h>
#include <chrono>
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc != 3 && argc != 4 ){
		cerr << "Usage: ./DBScanClustering AtomicInputFilename (VariableToCluster) OutputFilename" << endl;
		cerr << "This executable performs a density based scan (DBScan) clustering of atoms with respect to a VariableToCluster, which could be either the atomic positions or a given auxiliary property (which could have arbitrary dimension)" << endl;
		cerr << "The program will return an ovito output file with an auxiliary property named ClusterId, with ClusterId[1] is the id of the cluster (-1 => means undefined point, 0 => noisy point, else the cluster id) and ClusterId[2] is the status of the atoms (i.e. 1 => core, 0 => outlier point or -10 => not treated due to filtering)" << endl;
		cerr << "Parameters of DBScan are read from FixedParameters.ath file if it exists, if not default values will be used" << endl;
		cerr << "If VariableToCluster is not provided, the variable to cluster will be the atomic position" << endl;
		cerr << "This executable does not consider periodic boundary conditions even when the variable to cluster is the atomic position" << endl;
		Dis.Printer_DBScanClustering();	
		return EXIT_FAILURE;
	}
	Dis.Printer_DBScanClustering();	
	string InputFilename = argv[1];
	string OutputFilename;
	string VarToClust = "Position";
	if( argc == 4 ){
		VarToClust = argv[2];
		OutputFilename = argv[3];
	}else OutputFilename = argv[2];

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

		Descriptors MyDes(MySystem,VarToClust);
		DBScan MyDB;
		MyDB.setDescriptors(&MyDes);
		MyDB.TrainModel();
		MySystem->setAux_vec(MyDB.getClassificator(),2,"ClusterId");
	}

	if( argc != 4 ) VarToClust = "";

	if( nbSys == 1 ){
		MySystem->printSystem_aux(OutputFilename,"ClusterId "+VarToClust);
		delete MySystem;
	}else
		MyTraj.printSystem_aux(OutputFilename,"ClusterId "+VarToClust);

	Dis.ExecutionTime();	
	return 0;
}
