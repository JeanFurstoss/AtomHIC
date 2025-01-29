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
#include <chrono>
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc != 3 ){
		cerr << "Usage: ./DBScanClustering AtomicInputFilename OutputFilename" << endl;
		cerr << "AtomicInputFilename should be an atomic system" << endl;
		cerr << "The program will return an ovito output file with the id of cluster (parameters of DBScan are read from FixedParameters.dat)" << endl;
		cerr << "For the moment ions cannot be filtered" << endl;
		return EXIT_FAILURE;
	}
	
	string InputFilename = argv[1];
	string OutputFilename = argv[2];
	AtomicSystem MySystem(InputFilename);
	unsigned int nbAt = MySystem.getNbAtom();
	double *aux = new double[nbAt];
	unsigned int Struct_ind, size_Struct;
	Struct_ind = MySystem.getAuxIdAndSize("Struct",size_Struct);
	for(unsigned int i=0;i<nbAt;i++){
		if( MySystem.getAux(Struct_ind)[i*size_Struct] == 0 || MySystem.getAux(Struct_ind)[i*size_Struct] == 5 ) aux[i] = 1;
		else aux[i] = 0;
	}
	MySystem.setAux(aux,"Crystal");
	Descriptors MyDes(&MySystem,"Position","Crystal");
	DBScan MyDB;
	MyDB.setDescriptors(&MyDes);
	MyDB.TrainModel(to_string(0.));
	MySystem.setAux_vec(MyDB.getClassificator(),2,"ClusterId");
	MySystem.printSystem_aux(OutputFilename,"ClusterId");
	Dis.ExecutionTime();	
	return 0;
}
