//**********************************************************************************
//*   AtomicStrain/main.cpp                                                        *
//**********************************************************************************
//* This file contains the implementation of the AtomicStrain executable.          *
//* It allows to compute the atomic strain of given dump files relative to a 	   *
//* reference dump file								   *
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
//*	- make a generic function in AtomicSystem to adjust stoichiometry and put  *
//* this executable in generic tools						   *
//**********************************************************************************


#include <AtomicSystem.h>
#include <Bicrystal.h>
#include <Crystal.h>
#include <ComputeAuxiliary.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include "MyStructs.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	
	string InputFilename, OutputFilename, crystalName, AuxAtVolName;
	if( argc < 7 ){
		cerr << "Usage: ./GetStress2DMap InputFilename CrystalName AuxAtVolumeName NumberOfStressComp StressCompName1 StressCompName2 .. OutputNameFile" << endl;
		return EXIT_FAILURE;
	}
	InputFilename = argv[1];
	crystalName = argv[2];
	AuxAtVolName = argv[3];
	unsigned int nbaux;
	istringstream iss_nbaux(argv[4]);
	iss_nbaux >> nbaux;
	unsigned int current_read_ind = 5;
	vector<string> aux2treat;
	for(unsigned int i=0;i<nbaux;i++){
		aux2treat.push_back(argv[current_read_ind]);
		current_read_ind++;
	}
	OutputFilename = argv[current_read_ind];

	AtomicSystem MySystem(InputFilename);
	MySystem.setCrystal(crystalName);
	bool DDNS = false;
	vector<int> *DNS;
	unsigned int *Neigh;
	int *CLNeigh;
	unsigned nbN;
	cout << "Compute not sep list" << endl;
	MySystem.ComputeNotSepList();
	cout << "done" << endl;
	DNS = MySystem.getNotSepTag();
	
	cout << "Compute neighbours.." << endl;
	nbN = MySystem.searchNeighbours(7.);
	cout << "done" << endl;
	Neigh = MySystem.getNeighbours();
	CLNeigh = MySystem.getCLNeighbours();
	
	cout << "Getting aux" << endl;
	unsigned int buffer;
	unsigned int ind_v = MySystem.getAuxIdAndSize(AuxAtVolName,buffer);
	double *Vol = MySystem.getAux(ind_v);
	unsigned int fullauxsize = 0;
	vector<unsigned int> aux_sizes(nbaux);
	vector<double*> auxes(nbaux);
	unsigned int count_a;
	for(unsigned int i=0;i<nbaux;i++){
		unsigned int ind = MySystem.getAuxIdAndSize(aux2treat[i],aux_sizes[i]);
		fullauxsize += aux_sizes[i];
		auxes[i] = MySystem.getAux(ind);
	}
	const unsigned int nbAt = MySystem.getNbAtom();
	cout << "done" << endl;
	vector<double> coords;
	vector<double> toprint;
	unsigned int count;
	//unsigned int *Account = new unsigned int[nbAt];
	//for(unsigned int i=0;i<nbAt;i++) Account[i] = 0;
	for(unsigned int i=0;i<nbAt;i++){
		if( DNS[i][0] <= 0 ) continue;
		else{
			if( DNS[i][0] != 6 ){ // TODO to generalized
				//cout << "WARNING" << endl;
				continue;
			}
			//Account[i] += 1;
			count = 1;
			count += DNS[i][0];
			double mean_y = MySystem.getAtom(i).pos.y;
			double mean_z = MySystem.getAtom(i).pos.z;
			double cur_Vol = Vol[i];
			vector<double> curaux;
			for(unsigned int a=0;a<nbaux;a++){
				for(unsigned int d=0;d<aux_sizes[a];d++) curaux.push_back(auxes[a][i*aux_sizes[a]+d]);
			}
			// acount for the DNS list
			for(unsigned int n=0;n<DNS[i][0];n++){
				//Account[DNS[i][n+1]] += 1;
				cur_Vol += Vol[DNS[i][n+1]];
				mean_y += MySystem.getAtom(DNS[i][n+1]).pos.y;
				mean_z += MySystem.getAtom(DNS[i][n+1]).pos.z;
				count_a = 0;
				for(unsigned int a=0;a<nbaux;a++){
					for(unsigned int d=0;d<aux_sizes[a];d++){
						curaux[count_a] += auxes[a][i*aux_sizes[a]+d];
						count_a++;
					}
				}
			}
			// account for neighbours
			for(unsigned int j=0;j<Neigh[i*(nbN+1)];j++){
				unsigned int id=Neigh[i*(nbN+1)+1+j];
				if( DNS[id][0] != 6 ) continue;
				mean_y += MySystem.getAtom(id).pos.y;
				mean_z += MySystem.getAtom(id).pos.z;
				count++;
				cur_Vol += Vol[id];
				count_a = 0;
				for(unsigned int a=0;a<nbaux;a++){
					for(unsigned int d=0;d<aux_sizes[a];d++){
						curaux[count_a] += auxes[a][id*aux_sizes[a]+d];
						count_a++;
					}
				}
				for(unsigned int n=0;n<DNS[id][0];n++){
					cur_Vol += Vol[DNS[id][n+1]];
					mean_y += MySystem.getAtom(DNS[id][n+1]).pos.y;
					mean_z += MySystem.getAtom(DNS[id][n+1]).pos.z;
					count++;
					count_a = 0;
					for(unsigned int a=0;a<nbaux;a++){
						for(unsigned int d=0;d<aux_sizes[a];d++){
							curaux[count_a] += auxes[a][id*aux_sizes[a]+d];
							count_a++;
						}
					}
				}
			}

			coords.push_back(mean_y / count);
			coords.push_back(mean_z / count);
			count_a = 0;
			for(unsigned int a=0;a<nbaux;a++){
				for(unsigned int d=0;d<aux_sizes[a];d++){
					toprint.push_back(curaux[count_a] / cur_Vol);
					count_a++;
				}
			}
		}
	}
	//MySystem.setAux(Account,"aux");
	//MySystem.printSystem_aux("Verif.cfg","aux");
	//cout << "done, NbAt treated =" << count << ", trueNbAt= " << nbAt << endl;
	ofstream file(OutputFilename);
	file << "Y Z";
	for(unsigned int a=0;a<nbaux;a++){
		if( aux_sizes[a] > 1 ){
			for(unsigned int d=0;d<aux_sizes[a];d++) file << " " << aux2treat[a] << "[" << d+1 << "]";
		}else file << " " << aux2treat[a];
	}
	file << endl;
	unsigned int nb2print = coords.size()/2;
	for(unsigned int i=0;i<nb2print;i++){
		file << coords[i*2] << " " << coords[i*2+1];
		unsigned int count_a = 0;
		for(unsigned int a=0;a<nbaux;a++){
			for(unsigned int d=0;d<aux_sizes[a];d++){
				file << " " << toprint[i*fullauxsize+count_a];
				count_a++;
			}
		}
		file << endl;
	}
	file.close();

	Dis.ExecutionTime();	
	return 0;
}
