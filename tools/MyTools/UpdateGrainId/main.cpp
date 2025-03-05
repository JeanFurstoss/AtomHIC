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
		cerr << "Usage: ./UpdateGrainId AtomicInputFilename_preffix current_timestep previousFilename structNameClust" << endl;
		cerr << "This will update the grain id of the current dump with the respect to previous one by compraing centroid positions" << endl;
		return EXIT_FAILURE;
	}
	

	string InputFilename_pref = argv[1];
	unsigned int timestep;
	istringstream iss_timestep(argv[2]);
	iss_timestep >> timestep;
	string PrevInputFilename = argv[3];
	string structnameclust = argv[4];
	string ext=".cfg";
	string InputFilename=InputFilename_pref+to_string(timestep)+ext;

	unsigned int type_uint_GT = 2; // type of the ions used to search the cell vectors of grains

	AtomicSystem MySystem(InputFilename);
	AtomicSystem PrevSystem(PrevInputFilename);
	unsigned int Struct_ind, size_Struct;
	Struct_ind = MySystem.getAuxIdAndSize("ClusterId",size_Struct);
	unsigned int prev_Struct_ind, prev_size_Struct;
	prev_Struct_ind = PrevSystem.getAuxIdAndSize(structnameclust,prev_size_Struct);
	
	unsigned int nbAt = MySystem.getNbAtom();

	// get the previous and current number of grains and store their centroids
	vector<unsigned int> prevGrainIndex;
	vector<vector<double>> prevGrainCentroids;
	vector<unsigned int> prevNbAtGrains;
	vector<unsigned int> GrainIndex;
	vector<vector<double>> GrainCentroids;
	vector<unsigned int> nbAtGrains;
	vector<vector<unsigned int>> IndexForSearchCV;

	double *auxCV = new double[nbAt];
	cout << "Computing grain centroids" << endl;
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
		auxCV[i] = 0;
		current_gi = (int) MySystem.getAux(Struct_ind)[i*size_Struct];
		if( current_gi > 0 ){
			unsigned int ci = (unsigned int) current_gi;
			for(unsigned int n=0;n<GrainIndex.size();n++){
				if( GrainIndex[n] == ci ){
					nbAtGrains[n]++;
					GrainCentroids[n][0] += MySystem.getAtom(i).pos.x;
					GrainCentroids[n][1] += MySystem.getAtom(i).pos.y;
					GrainCentroids[n][2] += MySystem.getAtom(i).pos.z;
					if( MySystem.getAtom(i).type_uint == type_uint_GT ) IndexForSearchCV[n].push_back(i);
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
				IndexForSearchCV.push_back(vector<unsigned int>());
				if( MySystem.getAtom(i).type_uint == type_uint_GT ) IndexForSearchCV[IndexForSearchCV.size()-1].push_back(i);
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
	

	cout << "Changing grain id correspondance" << endl;	
	unsigned int *corres_array = new unsigned int[nb_grains];
	double *min_dist = new double[nb_grains];
	double *dist = new double[prev_nb_grains];
	MathTools MT;
	
	unsigned int maxGrainIndex = MT.max_vec(GrainIndex);
	unsigned int *corres_array_int = new unsigned int[maxGrainIndex+1];
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
	// attribute new grain id
	for(unsigned int i=0;i<nbAt;i++){
		int current_gi = (int) MySystem.getAux(Struct_ind)[i*size_Struct];
		if( current_gi > 0 ){
			unsigned int ci = (unsigned int) current_gi;
			auxGT[i] = corres_array[corres_array_int[ci]];
		}else auxGT[i] = -1.;
	}

	cout << "Done" << endl;
	MySystem.setAux_vec(auxGT,1,"grainId");

	// search cell vectors
	// get cell vectors from crystal
	Crystal MyCrystal("Forsterite");
	const double *CV = MyCrystal.getALength();
	double *CV_squared = new double[3];
	double rc = MT.max_p(CV,3);
	rc *= 1.25;
	for(unsigned int d=0;d<3;d++){
		CV_squared[d] = pow(CV[d],2);
	}
	double tol_length = 2.;
	double tol_sp = 1.e-1;
	vector<unsigned int> ori_ion(1);
	vector<vector<double>> FinalCellVecs;
	vector<unsigned int> CVTag;
	for(unsigned int n=0;n<nb_grains;n++){
		// use the ion closest to the centroid of grain to search cell vectors
		vector<double> dist;
		for(unsigned int i=0;i<IndexForSearchCV[n].size();i++){
			dist.push_back(0.);
			dist[i] += pow(MySystem.getAtom(IndexForSearchCV[n][i]).pos.x-GrainCentroids[n][0],2);
			dist[i] += pow(MySystem.getAtom(IndexForSearchCV[n][i]).pos.y-GrainCentroids[n][1],2);
			dist[i] += pow(MySystem.getAtom(IndexForSearchCV[n][i]).pos.z-GrainCentroids[n][2],2);
		}
		ori_ion[0] = IndexForSearchCV[n][MT.min(dist)];
		auxCV[ori_ion[0]] = corres_array[corres_array_int[GrainIndex[n]]]*10+1;
		bool cv_found = false;
		unsigned int count(0), max_count(100);
		while( !cv_found && count < max_count ){
			unsigned int nbMaxN = MySystem.searchNeighbours_restricted(rc,ori_ion,IndexForSearchCV[n]);
			unsigned int nbN = MySystem.getNeighbours()[0];
			double xpos = MySystem.getAtom(ori_ion[0]).pos.x;
			double ypos = MySystem.getAtom(ori_ion[0]).pos.y;
			double zpos = MySystem.getAtom(ori_ion[0]).pos.z;
			vector<vector<double>> dist_cv;
			for(unsigned int d=0;d<3;d++) dist_cv.push_back(vector<double>());
			for(unsigned int i=0;i<nbN;i++){
				double current_dist = 0.;
				current_dist += pow(xpos-MySystem.getAtom(MySystem.getNeighbours()[i+1]).pos.x, 2);
				current_dist += pow(ypos-MySystem.getAtom(MySystem.getNeighbours()[i+1]).pos.y, 2);
				current_dist += pow(zpos-MySystem.getAtom(MySystem.getNeighbours()[i+1]).pos.z, 2);
				for(unsigned int d=0;d<3;d++){
					if( fabs(current_dist-CV_squared[d]) < tol_length ){
						dist_cv[d].push_back(current_dist);
						dist_cv[d].push_back(i);
					}
				}
			}
			for(unsigned int d=0;d<3;d++){
				if( dist_cv[d].size() == 0 ) cout << "Issue when searching cell vectors" << endl;
				MT.sort(dist_cv[d],0,2,dist_cv[d]);
			}
			vector<vector<double>> good_combi;
			for(unsigned ind_a=0;ind_a<dist_cv[0].size()/2;ind_a++){
				unsigned int current_a_index = (unsigned int) dist_cv[0][ind_a*2+1];
				for(unsigned ind_b=0;ind_b<dist_cv[1].size()/2;ind_b++){
					unsigned int current_b_index = (unsigned int) dist_cv[1][ind_b*2+1];
					for(unsigned ind_c=0;ind_c<dist_cv[2].size()/2;ind_c++){
						unsigned int current_c_index = (unsigned int) dist_cv[2][ind_c*2+1];
						double sp_a(0.),sp_b(0.),sp_c(0.);
						double a_dist(0.),b_dist(0.),c_dist(0.);
						a_dist += pow(xpos-MySystem.getAtom(MySystem.getNeighbours()[current_a_index+1]).pos.x,2);
						a_dist += pow(ypos-MySystem.getAtom(MySystem.getNeighbours()[current_a_index+1]).pos.y,2);
						a_dist += pow(zpos-MySystem.getAtom(MySystem.getNeighbours()[current_a_index+1]).pos.z,2);
						a_dist = sqrt(a_dist);
						b_dist += pow(xpos-MySystem.getAtom(MySystem.getNeighbours()[current_b_index+1]).pos.x,2);
						b_dist += pow(ypos-MySystem.getAtom(MySystem.getNeighbours()[current_b_index+1]).pos.y,2);
						b_dist += pow(zpos-MySystem.getAtom(MySystem.getNeighbours()[current_b_index+1]).pos.z,2);
						b_dist = sqrt(b_dist);
						c_dist += pow(xpos-MySystem.getAtom(MySystem.getNeighbours()[current_c_index+1]).pos.x,2);
						c_dist += pow(ypos-MySystem.getAtom(MySystem.getNeighbours()[current_c_index+1]).pos.y,2);
						c_dist += pow(zpos-MySystem.getAtom(MySystem.getNeighbours()[current_c_index+1]).pos.z,2);
						c_dist = sqrt(c_dist);
						sp_a += (xpos-MySystem.getAtom(MySystem.getNeighbours()[current_c_index+1]).pos.x)*(xpos-MySystem.getAtom(MySystem.getNeighbours()[current_b_index+1]).pos.x);
						sp_a += (ypos-MySystem.getAtom(MySystem.getNeighbours()[current_c_index+1]).pos.y)*(ypos-MySystem.getAtom(MySystem.getNeighbours()[current_b_index+1]).pos.y);
						sp_a += (zpos-MySystem.getAtom(MySystem.getNeighbours()[current_c_index+1]).pos.z)*(zpos-MySystem.getAtom(MySystem.getNeighbours()[current_b_index+1]).pos.z);
						sp_b += (xpos-MySystem.getAtom(MySystem.getNeighbours()[current_a_index+1]).pos.x)*(xpos-MySystem.getAtom(MySystem.getNeighbours()[current_c_index+1]).pos.x);
						sp_b += (ypos-MySystem.getAtom(MySystem.getNeighbours()[current_a_index+1]).pos.y)*(ypos-MySystem.getAtom(MySystem.getNeighbours()[current_c_index+1]).pos.y);
						sp_b += (zpos-MySystem.getAtom(MySystem.getNeighbours()[current_a_index+1]).pos.z)*(zpos-MySystem.getAtom(MySystem.getNeighbours()[current_c_index+1]).pos.z);
						sp_c += (xpos-MySystem.getAtom(MySystem.getNeighbours()[current_a_index+1]).pos.x)*(xpos-MySystem.getAtom(MySystem.getNeighbours()[current_b_index+1]).pos.x);
						sp_c += (ypos-MySystem.getAtom(MySystem.getNeighbours()[current_a_index+1]).pos.y)*(ypos-MySystem.getAtom(MySystem.getNeighbours()[current_b_index+1]).pos.y);
						sp_c += (zpos-MySystem.getAtom(MySystem.getNeighbours()[current_a_index+1]).pos.z)*(zpos-MySystem.getAtom(MySystem.getNeighbours()[current_b_index+1]).pos.z);
						sp_a = fabs(sp_a/(b_dist*c_dist));
						sp_b = fabs(sp_b/(a_dist*c_dist));
						sp_c = fabs(sp_c/(a_dist*b_dist));
						if( sp_a < tol_sp && sp_b < tol_sp && sp_c < tol_sp ){
							cv_found = true;
							unsigned int current_gi_good = corres_array[corres_array_int[GrainIndex[n]]];
							auxCV[MySystem.getNeighbours()[current_a_index+1]] = current_gi_good*10+2;
							auxCV[MySystem.getNeighbours()[current_b_index+1]] = current_gi_good*10+3;
							auxCV[MySystem.getNeighbours()[current_c_index+1]] = current_gi_good*10+4;
							CVTag.push_back(current_gi_good);
							FinalCellVecs.push_back(vector<double>());
							FinalCellVecs[FinalCellVecs.size()-1].push_back(xpos-MySystem.getAtom(MySystem.getNeighbours()[current_a_index+1]).pos.x);
							FinalCellVecs[FinalCellVecs.size()-1].push_back(ypos-MySystem.getAtom(MySystem.getNeighbours()[current_a_index+1]).pos.y);
							FinalCellVecs[FinalCellVecs.size()-1].push_back(zpos-MySystem.getAtom(MySystem.getNeighbours()[current_a_index+1]).pos.z);
							FinalCellVecs[FinalCellVecs.size()-1].push_back(xpos-MySystem.getAtom(MySystem.getNeighbours()[current_b_index+1]).pos.x);
							FinalCellVecs[FinalCellVecs.size()-1].push_back(ypos-MySystem.getAtom(MySystem.getNeighbours()[current_b_index+1]).pos.y);
							FinalCellVecs[FinalCellVecs.size()-1].push_back(zpos-MySystem.getAtom(MySystem.getNeighbours()[current_b_index+1]).pos.z);
							FinalCellVecs[FinalCellVecs.size()-1].push_back(xpos-MySystem.getAtom(MySystem.getNeighbours()[current_c_index+1]).pos.x);
							FinalCellVecs[FinalCellVecs.size()-1].push_back(ypos-MySystem.getAtom(MySystem.getNeighbours()[current_c_index+1]).pos.y);
							FinalCellVecs[FinalCellVecs.size()-1].push_back(zpos-MySystem.getAtom(MySystem.getNeighbours()[current_c_index+1]).pos.z);
							break;
						}
					}
					if( cv_found ) break;
				}
				if( cv_found ) break;
			}
			if( !cv_found ){
				cout << "We do not find the cell vectors" << endl;
				cout << "changing origin ion" << endl;
				count++;
				auxCV[ori_ion[0]] = 0.;
				ori_ion[0] = IndexForSearchCV[n][rand()%(IndexForSearchCV[n].size())];
				auxCV[ori_ion[0]] = corres_array[corres_array_int[GrainIndex[n]]]*10+1;
			}
		}
		if( !cv_found ) cout << "We finally not managed to find cell vectors" << endl;
		else cout << "CELL VEC FIND !" << endl;
	}
	MySystem.setAux_vec(auxCV,1,"GrainTag");
	
	string suf = "RelabeledGrains_";
	string OutputFilename = suf+to_string(timestep)+ext;
	
	MySystem.printSystem_aux(OutputFilename,"grainId GrainTag");

	cout << "Writing GrainCellVec.dat file.." << endl;
	ofstream writefile("GrainCellVec_"+to_string(timestep)+".dat");
	for(unsigned int n=0;n<CVTag.size();n++){
		for(unsigned int d1=0;d1<3;d1++){
			writefile << CVTag[n] << " " << d1 << " ";
			for(unsigned int d2=0;d2<3;d2++) writefile << FinalCellVecs[n][d1*3+d2] << " ";
			writefile << endl;
		}
	}
	writefile.close();
	cout << "Done" << endl;


	delete[] corres_array;
	delete[] corres_array_int;
	delete[] auxGT;
	delete[] auxCV;
	delete[] min_dist;
	delete[] dist;
	delete[] CV_squared;
	Dis.ExecutionTime();	
	return 0;
}
