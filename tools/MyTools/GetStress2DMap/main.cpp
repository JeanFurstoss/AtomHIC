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
	
	string InputFilename, OutputFilename, crystalName;
	if( argc >= 3 ){
		InputFilename = argv[1];
		crystalName = argv[2];
		OutputFilename = argv[3];
	}else{
		cerr << "Usage: ./AdjustStoichiometry InputFilename OutputFilename AuxName2Print_1 AuxName2Print_2.." << endl;
		cerr << "Adjsut the stoichiometry of Mg2SiO4 system by removing Si, Mg and O ions, if we have to remove Si ions we assume that SiO4 tetrahedra have not been separated" << endl;
		cerr << "TODO AuxName2Print implementation" << endl;
		return EXIT_FAILURE;
	}
	AtomicSystem MySystem(InputFilename);
	MySystem.setCrystal(crystalName);
	cout << "Compute not sep list" << endl;
	MySystem.ComputeNotSepList();
	cout << "done" << endl;
	vector<int> *DNS = MySystem.getNotSepTag();
	unsigned int buffer;
	cout << "Getting aux" << endl;
	unsigned int ind_p = MySystem.getAuxIdAndSize("Pressure",buffer);
	double *Press = MySystem.getAux(ind_p);
	unsigned int ind_v = MySystem.getAuxIdAndSize("AtomicVolume",buffer);
	double *Vol = MySystem.getAux(ind_v);
	unsigned int ind_tau = MySystem.getAuxIdAndSize("ShearStress",buffer);
	double *Tau = MySystem.getAux(ind_tau);
	const unsigned int nbAt = MySystem.getNbAtom();
	cout << "done" << endl;
	vector<double> coords;
	vector<double> PAndS;
	cout << "Compute pressure and shear stress" << endl;
	for(unsigned int i=0;i<nbAt;i++){
		if( DNS[i][0] <= 0 ) continue;
		else{
			double mean_y = MySystem.getAtom(i).pos.y;
			double mean_z = MySystem.getAtom(i).pos.z;
			double cur_Vol = Vol[i];
			double cur_Press = Press[i];
			double cur_Tau = Tau[i];
			for(unsigned int n=0;n<DNS[i][0];n++){
				cur_Vol += Vol[DNS[i][n+1]];
				cur_Press += Press[DNS[i][n+1]];
				cur_Tau += Tau[DNS[i][n+1]];
				mean_y += MySystem.getAtom(DNS[i][n+1]).pos.y;
				mean_z += MySystem.getAtom(DNS[i][n+1]).pos.z;
			}
			coords.push_back(mean_y / (DNS[i][0]+1));
			coords.push_back(mean_z / (DNS[i][0]+1));
			PAndS.push_back(cur_Press / cur_Vol);
			PAndS.push_back(cur_Tau / cur_Vol);
		}
	}
	cout << "done" << endl;

	ofstream file(OutputFilename);
	file << "Y Z Press(bar) Tau(bar)" << endl;
	for(unsigned int i=0;i<coords.size()/2;i++) file << coords[i*2] << " " << coords[i*2+1] << " " << PAndS[i*2] << " " << PAndS[i*2+1] << endl;
	file.close();
	Dis.ExecutionTime();	
	return 0;
}
