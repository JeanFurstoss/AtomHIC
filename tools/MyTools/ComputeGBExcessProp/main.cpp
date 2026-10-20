//**********************************************************************************
//*   AnalyzeBicrystal_AndAtomicStrain/main.cpp                                    *
//**********************************************************************************
//* This file contains the implementation of the AnalyzeBicrystal_AndAtomicStrain  *
//* executable.									   *
//* It allows to find the position, width and excess volume of a GB and compute    *
//* the atomic strain								   *
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

#include <Bicrystal.h>
#include <AtomicSystem.h>
#include <ComputeAuxiliary.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();

	// check the number of argument
	if( argc < 4 ){
		cerr << "Usage: ComputeGBExcessProp AtomicInputFilename CrystalName GBNormalDirection NumberOfProperties NameOfProp1 NameOfProp2 .. NameOfPropN" << endl;
		return EXIT_FAILURE;
	}
	string InputFilename = argv[1];
	string NormalDir = argv[3];
	string CN = argv[2];
	unsigned int nbaux;
	istringstream iss_rc(argv[4]);
	iss_rc >> nbaux;
	unsigned int current_read_ind = 5;
	vector<string> aux2treat;
	for(unsigned int i=0;i<nbaux;i++){
		aux2treat.push_back(argv[current_read_ind]);
		current_read_ind++;
	}

	Bicrystal MySystem(InputFilename, NormalDir, CN);
	MySystem.Print1dDensity("ExcessVolume.dat", "GBProfile");
	for(unsigned int i=0;i<nbaux;i++){
		MySystem.ComputeExcessProp(aux2treat[i]);
		MySystem.Print1dDensity("Density_"+aux2treat[i]+".dat", aux2treat[i]);
	}
	unsigned int nbGB = MySystem.getNbGB();
	ofstream writefile("ExcessProperties.dat");
	if( nbGB == 2 ) writefile << "\t GB1 \t GB2" << endl;
	writefile << "Position\t" << MySystem.getGBPos1();
        if( nbGB == 2 ) writefile << "\t" << MySystem.getGBPos2();
        writefile << endl;
	writefile << "Width\t" << "\t" << MySystem.getGBwidth1();
        if( nbGB == 2 ) writefile << "\t" << MySystem.getGBwidth2();
        writefile << endl;
	writefile << "ExcessVolume\t" << MySystem.getExcessVol1();
        if( nbGB == 2 ) writefile << "\t" << MySystem.getExcessVol2();
        writefile << endl;
	vector<string> PropName = MySystem.getGBExcessPropName();
	vector<double> EP = MySystem.getGBExcessProp();
	for(unsigned int i=0;i<nbaux;i++){
		for(unsigned int j=0;j<PropName.size();j++){
			if( PropName[j] == aux2treat[i] ){
				writefile << "Excess_" << aux2treat[i] << "\t" << EP[j*nbGB];
				if( nbGB == 2 ) writefile << "\t" << EP[j*2+1];
				writefile << endl;
				break;
			}
		}
	}
	writefile.close();
	cout << "ExcessProperties.dat file successfully writted" << endl;
	Dis.ExecutionTime();	
	return 0;
}
