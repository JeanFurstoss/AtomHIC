//**********************************************************************************
//*   CreateGB/main.cpp                                                            *
//**********************************************************************************
//* This file contains the implementation of the CreateGB executable.              *
//* It allows to create an atomic system with a GB with a given misorientation and *
//* GB plane									   *
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
//*	- work on Bicrystal class to have more clear outputs			   *
//**********************************************************************************

#include <AtomicSystem.h>
#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "ComputeAuxiliary.h"
#include <Displays.h>

using namespace std;


int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc < 9 ){
		cerr << "Usage: CreateGB h_RotAxis k_RotAxis l_RotAxis RotAngle(in degree) h_GBPlane k_GBPlane l_GBPlane CrystalName(has to be defined in /data/Crystal/) Rationalize" << endl;
		cerr << "Rationalize can be either 0 or 1" << endl;
		cerr << "1 => rationalize the GB (i.e. search the closest CSL GB to the provided parameters, in this case the parameters for CSL calculation in /data/FixedParameters/FixedParameters.dat can be important)" << endl;
		cerr << "0 => do not rationalize the GB, in this case the parameters for constructing crystal and bicrystal in /data/FixedParameters/FixedParameters.dat can be important" << endl;
		cerr << "The program will return 4 dump files containing the GB system, the CSL lattice and the two grains (the two latters can be used to change shift between crystals)" << endl;
		return EXIT_FAILURE;
	}
	int h_a, k_a ,l_a, h_p, k_p, l_p;
	double theta;
	unsigned int rat;
	string crystalName;
	istringstream iss_ha(argv[1]);
	iss_ha >> h_a;
	istringstream iss_ka(argv[2]);
	iss_ka >> k_a;
	istringstream iss_la(argv[3]);
	iss_la >> l_a;
	istringstream iss_theta(argv[4]);
	iss_theta >> theta;
	theta *= M_PI/180.;
	istringstream iss_hp(argv[5]);
	iss_hp >> h_p;
	istringstream iss_kp(argv[6]);
	iss_kp >> k_p;
	istringstream iss_lp(argv[7]);
	iss_lp >> l_p;
	istringstream iss_cn(argv[8]);
	iss_cn >> crystalName;
	istringstream iss_rat(argv[9]);
	iss_rat >> rat;
	bool rat_b;
	if( rat == 0 ) rat_b = false;
	else rat_b = true;
	Bicrystal MyGB(crystalName,h_a,k_a,l_a,theta,h_p,k_p,l_p,rat_b);
	MyGB.print_lmp("GB.lmp");
	MyGB.printCSL("CSL.lmp");
	MyGB.print_Grains();
	Dis.ExecutionTime();	
	return 0;
}
