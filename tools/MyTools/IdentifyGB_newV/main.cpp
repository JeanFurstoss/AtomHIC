//**********************************************************************************
//*   IdentifyGrains/main.cpp		                                   	   *
//**********************************************************************************
//* This file contains the implementation of the IdentifyGrains executable.	   *
//* It is used to analyze polycrystaline deformation simulations		   *
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
#include "SteinhardtDescriptors.h"
#include "GaussianMixtureModel.h"
#include "AtomicSystem.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc != 5 ){
		cerr << "Usage: ./IdentifyGB_newV GrainIdInputFilename DBScanOutFilename rcut outputFilename" << endl;
		return EXIT_FAILURE;
	}
	

	string GrainIdInputFilename = argv[1];
	double rcut;
	istringstream iss_rcut(argv[3]);
	iss_rcut >> rcut;
	string DBScanFilename = argv[2];
	string OutputFilename = argv[4];

	AtomicSystem GISystem(GrainIdInputFilename);
	AtomicSystem CoreOutSystem(DBScanFilename);
	unsigned int CoreOut_Struct_ind, CoreOut_size_Struct;
	CoreOut_Struct_ind = CoreOutSystem.getAuxIdAndSize("ClusterId",CoreOut_size_Struct);
	unsigned int GI_Struct_ind, GI_size_Struct;
	GI_Struct_ind = GISystem.getAuxIdAndSize("grainId",GI_size_Struct);
	
	unsigned int nbAt = GISystem.getNbAtom();
	MathTools MT;

	// get index of ions to which we need to search neighbours and the index in which searching 
	double *auxGBId = new double[nbAt*2];
	vector<unsigned int> IndexToSearch;
	vector<unsigned int> IndexForSearch;
	cout << "Sorting ions" << endl;
	for(unsigned int i=0;i<nbAt;i++){
		auxGBId[i*2] = 0.;
		auxGBId[i*2+1] = 0.;
		if( CoreOutSystem.getAux(CoreOut_Struct_ind)[i*CoreOut_size_Struct] > 0. && CoreOutSystem.getAux(CoreOut_Struct_ind)[i*CoreOut_size_Struct+1] == 0. ) IndexForSearch.push_back(i);
		else if( CoreOutSystem.getAux(CoreOut_Struct_ind)[i*CoreOut_size_Struct] == -10. ) IndexToSearch.push_back(i);
	}
			
	unsigned long int nbMaxN = CoreOutSystem.searchNeighbours_restricted(rcut,IndexToSearch,IndexForSearch);
	cout << "Searching clostest grains" << endl;
	unsigned long int one_l(1), two_l(2), three_l(3);
	#pragma omp parallel for
	for(unsigned long int i=0;i<IndexToSearch.size();i++){
	       unsigned long int current_nbN = CoreOutSystem.getNeighbours()[i*(nbMaxN+one_l)];
       	       vector<double> dist_arr(2*current_nbN);
	       double xpos = CoreOutSystem.getWrappedPos(IndexToSearch[i]).x;
	       double ypos = CoreOutSystem.getWrappedPos(IndexToSearch[i]).y;
	       double zpos = CoreOutSystem.getWrappedPos(IndexToSearch[i]).z;
	       unsigned int first_gi_index = 0;
	       for(unsigned long int n=0;n<current_nbN;n++){
		       unsigned int current_index = CoreOutSystem.getNeighbours()[i*(nbMaxN+one_l)+one_l+n];
		       double dist = pow(CoreOutSystem.getWrappedPos(current_index).x+CoreOutSystem.getCLNeighbours(i*nbMaxN*three_l+n*three_l)*CoreOutSystem.getH1()[0]+CoreOutSystem.getCLNeighbours(i*nbMaxN*three_l+n*three_l+one_l)*CoreOutSystem.getH2()[0]+CoreOutSystem.getCLNeighbours(i*nbMaxN*three_l+n*three_l+two_l)*CoreOutSystem.getH3()[0]-xpos,2.); // CL !!!!
		       dist += pow(CoreOutSystem.getWrappedPos(current_index).y+CoreOutSystem.getCLNeighbours(i*nbMaxN*three_l+n*three_l)*CoreOutSystem.getH1()[1]+CoreOutSystem.getCLNeighbours(i*nbMaxN*three_l+n*three_l+one_l)*CoreOutSystem.getH2()[1]+CoreOutSystem.getCLNeighbours(i*nbMaxN*three_l+n*three_l+two_l)*CoreOutSystem.getH3()[1]-ypos,2.); // CL !!!!
		       dist += pow(CoreOutSystem.getWrappedPos(current_index).z+CoreOutSystem.getCLNeighbours(i*nbMaxN*three_l+n*three_l)*CoreOutSystem.getH1()[2]+CoreOutSystem.getCLNeighbours(i*nbMaxN*three_l+n*three_l+one_l)*CoreOutSystem.getH2()[2]+CoreOutSystem.getCLNeighbours(i*nbMaxN*three_l+n*three_l+two_l)*CoreOutSystem.getH3()[2]-zpos,2.); // CL !!!!
		       dist_arr[n*2] = dist;
		       dist_arr[n*2+1] = (double) current_index;
		       if( dist < dist_arr[first_gi_index*2] ) first_gi_index = n;
	       }
	       if( dist_arr.size() < 4 ){
		       cout << "Issue with number of neighbours" << endl;
		       continue;
	       }
	       //MT.sort(dist_arr,0,2,dist_arr);
		unsigned int first_gi = (unsigned int) GISystem.getAux(GI_Struct_ind)[((unsigned int) dist_arr[first_gi_index*2+1])*GI_size_Struct];
		if( first_gi > 30 ){
			cout << "ISSUE" << endl;
			cout << "neighbour index = " << dist_arr[first_gi_index*2+1] << endl;
			cout << "converted neighbour index = " << (unsigned int) dist_arr[first_gi_index*2+1] << endl;
			cout << "Grain index " << first_gi << endl;
		}
		auxGBId[IndexToSearch[i]*2] = first_gi;
		unsigned int second_gi_index = 0;
		double min_dist = 1.e6;
	       for(unsigned int n=0;n<current_nbN;n++){
		       unsigned int current_gi = (unsigned int) GISystem.getAux(GI_Struct_ind)[((unsigned int) dist_arr[n*2+1])*GI_size_Struct];
		       if( current_gi != first_gi ){
			       if( dist_arr[n*2] < min_dist ){
				       min_dist = dist_arr[n*2];
				       second_gi_index = n;
			       }
		       }
	       }
	       if( min_dist == 1e6 ) cout << "Issue when identifying second grain" << endl;
	       else auxGBId[IndexToSearch[i]*2+1] = (unsigned int) GISystem.getAux(GI_Struct_ind)[((unsigned int) dist_arr[second_gi_index*2+1])*GI_size_Struct];
	}



	GISystem.setAux_vec(auxGBId,2,"GBId");
	GISystem.setAux_vec(CoreOutSystem.getAux(CoreOut_Struct_ind),2,"clust");
	
	GISystem.printSystem_aux(OutputFilename,"GBId clust grainId");
	cout << "Done" << endl;


	delete[] auxGBId;
	Dis.ExecutionTime();	
	return 0;
}
