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
	if( argc != 4 && argc != 5 ){
		cerr << "Usage: ./UpdateGrainId AtomicInputFilename PreviousDumpFilename" << endl;
		cerr << "This will update the grain id of the current dump with the respect to previous one by compraing centroid positions" << endl;
		return EXIT_FAILURE;
	}
	

	string InputFilename = argv[1];
	string PrevInputFilename = argv[2];
	AtomicSystem MySystem(InputFilename);
	AtomicSystem PrevSystem(PrevInputFilename);
	unsigned int Struct_ind, size_Struct;
	Struct_ind = MySystem.getAuxIdAndSize("ClusterId",size_Struct);
	unsigned int prev_Struct_ind, prev_size_Struct;
	prev_Struct_ind = PrevSystem.getAuxIdAndSize("ClusterId",prev_size_Struct);
	
	unsigned int nbAt = MySystem.getNbAtom();

	// get the previous and current number of grains and store their centroids
	vector<unsigned int> prevGrainIndex;
	vector<vector<double>> prevGrainCentroids;
	vector<unsigned int> prevNbAtGrains;
	vector<unsigned int> GrainIndex;
	vector<vector<double>> GrainCentroids;
	vector<unsigned int> nbAtGrains;
	for(unsigned int i=0;i<nbAt;i++){
		bool already = false;
		int current_gi = (int) PrevSystem.getAux(prev_Struct_ind)[i*prev_size_Struct];
		if( current_gi > 0 ){
			unsigned int ci = (unsigned int) current_gi;
			for(unsigned int n=0;n<prevGrainIndex.size();n++){
				if( prevGrainIndex[n] == ci ){
					prevNbAtGrains[n]++;
					prevGrainCentroids[n][0] += PrevSystem.getAtom(i).pos.x;
					prevGrainCentroids[n][1] += PrevSystem.getAtom(i).pos.y;
					prevGrainCentroids[n][2] += PrevSystem.getAtom(i).pos.z;
					already = true;
					break;
				}
			}
			if( !already ){
				prevGrainIndex.push_back(ci);
				prevNbAtGrains.push_back(1);
				prevGrainCentroids.push_back(vector<double>());
				prevGrainCentroids[prevGrainCentroids.size()-1].push_back(PrevSystem.getAtom(i).pos.x);
				prevGrainCentroids[prevGrainCentroids.size()-1].push_back(PrevSystem.getAtom(i).pos.y);
				prevGrainCentroids[prevGrainCentroids.size()-1].push_back(PrevSystem.getAtom(i).pos.z);
			}
		}
		already = false;
		current_gi = (int) MySystem.getAux(Struct_ind)[i*size_Struct];
		if( current_gi > 0 ){
			unsigned int ci = (unsigned int) current_gi;
			for(unsigned int n=0;n<GrainIndex.size();n++){
				if( GrainIndex[n] == ci ){
					nbAtGrains[n]++;
					GrainCentroids[n][0] += MySystem.getAtom(i).pos.x;
					GrainCentroids[n][1] += MySystem.getAtom(i).pos.y;
					GrainCentroids[n][2] += MySystem.getAtom(i).pos.z;
					already = true;
					break;
				}
			}
			if( !already ){
				GrainIndex.push_back(ci);
				nbAtGrains.push_back(1);
				GrainCentroids.push_back(vector<double>());
				GrainCentroids[GrainCentroids.size()-1].push_back(MySystem.getAtom(i).pos.x);
				GrainCentroids[GrainCentroids.size()-1].push_back(MySystem.getAtom(i).pos.y);
				GrainCentroids[GrainCentroids.size()-1].push_back(MySystem.getAtom(i).pos.z);
			}
		}
	}


	unsigned int prev_nb_grains = prevGrainIndex.size();
	unsigned int nb_grains = GrainIndex.size();
	for(unsigned int n=0;n<prev_nb_grains;n++){
		for(unsigned int d=0;d<3;d++) prevGrainCentroids[n][d] /= prevNbAtGrains[n];
	}
	for(unsigned int n=0;n<nb_grains;n++){
		for(unsigned int d=0;d<3;d++) GrainCentroids[n][d] /= nbAtGrains[n];
	}
	
	unsigned int *corres_array = new unsigned int[nb_grains];
	double *min_dist = new double[nb_grains];
	double *dist = new double[prev_nb_grains];
	MathTools MT;
	
	unsigned int maxGrainIndex = MT.max_vec(GrainIndex);
	unsigned int *corres_array_int = new unsigned int[maxGrainIndex];
	for(unsigned int i=0;i<nb_grains;i++) corres_array_int[GrainIndex[i]] = i;
	
	// atribute grain id to the ones with the closest centroids
	for(unsigned int n1=0;n1<nb_grains;n1++){
		for(unsigned int n2=0;n2<prev_nb_grains;n2++){
			dist[n2] = 0.;
			for(unsigned int d=0;d<3;d++) dist[n2] += pow(GrainCentroids[n1][d]-prevGrainCentroids[n2][d],2);
		}
		unsigned int min_ind = MT.min_p_ind(dist,prev_nb_grains);
		corres_array[n1] = prevGrainIndex[min_ind];
		min_dist[n1] = dist[min_ind];
	}

	// if higher number of grains search the one with closest centroid
	if( nb_grains > prev_nb_grains ){
		unsigned int current_inc = 1;
		for(unsigned int n1=0;n1<nb_grains;n1++){
			bool only_one = true;
			vector<unsigned int> others;
			for(unsigned int n2=0;n2<nb_grains;n2++){
				if( n1 != n2 && corres_array[n1] == corres_array[n2] ){
					only_one = false;
					others.push_back(n2);
				}
			}
			if( !only_one ){
				others.push_back(n1);
				unsigned int min_ind = 0;
				for(unsigned int i=1;i<others.size();i++){
					if( min_dist[others[min_ind]] > min_dist[others[i]] ) min_ind = i;
				}
				for(unsigned int i=0;i<others.size();i++){
					if( i != min_ind ){
						corres_array[others[i]] = prev_nb_grains+current_inc;
						current_inc++;
					}
				}
			}
		}
	}





	double *auxGT = new double[nbAt];	
	
	for(unsigned int i=0;i<nbAt;i++){
		unsigned int current_gi = (int) MySystem.getAux(Struct_ind)[i*size_Struct];
		if( current_gi > 0 ){
			unsigned int ci = (unsigned int) current_gi;
			auxGT[i] = corres_array[corres_array_int[ci]];
		}else auxGT[i] = -1.;
	}


	MySystem.setAux_vec(auxGT,1,"grainId");

	string suf = "RelabedGrains_";
	string OutputFilename = suf+InputFilename;
	
	MySystem.printSystem_aux(OutputFilename,"grainId");

	delete[] corres_array;
	delete[] corres_array_int;
	delete[] auxGT;
	delete[] min_dist;
	delete[] dist;
	Dis.ExecutionTime();	
	return 0;
}
