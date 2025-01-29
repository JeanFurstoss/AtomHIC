//**********************************************************************************
//*   ClusterAnalysis/main.cpp	          	                                   *
//**********************************************************************************
//* This file contains the implementation of the ClusterAnalysis executable.	   *
//* It is used to cluster atomic position using GMM				   *
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
//*	- make a generic ClusterAnalysis tool where we can choose the method and   *
//* the descriptors (cf. ../DBScanClustering/)					   *
//**********************************************************************************

#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "Descriptors.h"
#include "SteinhardtDescriptors.h"
#include "GaussianMixtureModel.h"
#include "AtomicSystem.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc < 5 ){
		cerr << "Usage: ./ClusterAnalysis AtomicInputFilename n_clust_min n_clust_max OutputFilename" << endl;
		cerr << "AtomicInputFilename should be an atomic system" << endl;
		cerr << "The executable returns two auxialiary properties Clust[1] containing the id of cluster which the ions belongs to and Clust[2] containing the MLC" << endl;
		cerr << "For the moment ions cannot be filtered" << endl;
		return EXIT_FAILURE;
	}
	
	string InputFilename = argv[1];
	unsigned int nclust_min, nclust_max;
	istringstream iss_min(argv[2]);
	iss_min >> nclust_min;
	istringstream iss_max(argv[3]);
	iss_max >> nclust_max;
	string OutputFilename = argv[4];
	AtomicSystem MySystem(InputFilename);
	double *pos = new double[MySystem.getNbAtom()*3];
	for(unsigned int i=0;i<MySystem.getNbAtom();i++){
		pos[i*3] = MySystem.getAtom(i).pos.x;
		pos[i*3+1] = MySystem.getAtom(i).pos.y;
		pos[i*3+2] = MySystem.getAtom(i).pos.z;
	}
	Descriptors MyDes(pos,MySystem.getNbAtom(),3);
	GaussianMixtureModel GMM;
	GMM.setDescriptors(&MyDes);
	GMM.fitOptimalGMM(nclust_min,nclust_max);
	GMM.Classify();
	MySystem.setAux_vec(GMM.getClassificator(),2,"Clust");
	MySystem.printSystem_aux(OutputFilename,"Clust");
	Dis.ExecutionTime();	
	delete[] pos;
}
