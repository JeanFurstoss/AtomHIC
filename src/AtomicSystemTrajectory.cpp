//**********************************************************************************
//*   AtomicSystemTrajectory.cpp                                                   *
//**********************************************************************************
//* This file contains the implementation of the AtomicSystemTrajectory class      *
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


#include "AtomicSystemTrajectory.h"
#include "AtomHicConfig.h"
#include <cmath>
#include <complex>
#include <omp.h>
#include <iomanip>

using namespace std; 

AtomicSystemTrajectory::AtomicSystemTrajectory(){
}

bool AtomicSystemTrajectory::SearchIsTrajectory(const string &_input_filename){
	timesteps.clear();
	ifstream file(_input_filename, ios::in);
	unsigned int current_nbSys(0);
	vector<unsigned int> ReadOk;
	if(file){
		vector<unsigned int> line_dt;
		unsigned int count(0);
		unsigned int buffer_1;
		size_t pos_dt, pos_nbat, pos_box, pos_list;
		string line;
		while(file){
			getline(file,line);
			// find timestep
			pos_dt=line.find("TIMESTEP");
			if( pos_dt!=string::npos ){
				line_dt.push_back(count);
				ReadOk.push_back(1);
				current_nbSys++;
			}
			if( current_nbSys != 0 ){
				if( count == line_dt[current_nbSys-1]+1 ){
					istringstream text(line);
					text >> buffer_1;
				       	timesteps.push_back(buffer_1);
				}
				pos_nbat=line.find("NUMBER OF ATOMS");
				if( pos_nbat!=string::npos &&  count == line_dt[current_nbSys-1]+2 ) ReadOk[current_nbSys-1]++;
				pos_box=line.find("BOX BOUNDS");
				if( pos_box!=string::npos &&  count == line_dt[current_nbSys-1]+4 ) ReadOk[current_nbSys-1]++;
				pos_list=line.find("ITEM: ATOMS");
				if( pos_list!=string::npos &&  count == line_dt[current_nbSys-1]+8 ) ReadOk[current_nbSys-1]++;
			}
			count++;
		}
		file.close();
	}else{
		cerr << "The file " << _input_filename << " cannot be openned" << endl;
		exit(EXIT_FAILURE);
	}

	if( current_nbSys < 2 ) return false;
	nbSys = current_nbSys;
	vector<unsigned int> NotOk;
	for(unsigned int i=0;i<nbSys;i++)
		if( ReadOk[i] != 4 ) NotOk.push_back(timesteps[i]);

	if( NotOk.size() != 0 ){
		cerr << "Issue when reading file " << _input_filename << ", for timesteps ";
		for(unsigned int i=0;i<NotOk.size();i++) cerr << NotOk[i] << " ";
		cerr << endl;
		exit(EXIT_FAILURE);
	}
	IsTraj = true;	
	return true;
}

void AtomicSystemTrajectory::setAtomicSystemList(const string &_input_filename){
	for(unsigned int i=0;i<AtSysList.size();i++)
		delete AtSysList[i];
	AtSysList.clear();
	AtSysList = vector<AtomicSystem*>(nbSys);
	if( IsTraj ){
		cout << "Reading " << _input_filename << " trajectory file .." << endl;
		for(unsigned int i=0;i<nbSys;i++)
			AtSysList[i] = new AtomicSystem(_input_filename,(int) timesteps[i]);
		cout << "Done !" << endl << endl;
	}
}

void AtomicSystemTrajectory::printSystem_aux(const std::string& filename, const std::string& AuxId){
	if( IsTraj ){
		unsigned int atListSize = AtSysList.size();
		AtSysList[0]->printSystem_aux(filename,AuxId);
		for(unsigned int i=1;i<atListSize;i++)
			AtSysList[i]->printSystem_aux(filename,AuxId,true);
	}
}

AtomicSystemTrajectory::~AtomicSystemTrajectory(){
	for(unsigned int i=0;i<AtSysList.size();i++)
		delete AtSysList[i];
	AtSysList.clear();
	for(unsigned int i=0;i<BicrystalList.size();i++)
		delete BicrystalList[i];
	BicrystalList.clear();
}
