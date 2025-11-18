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
	if( argc != 13 && argc != 15 ){
		cerr << "Usage: SampleGBComplexionCSL h_RotAxis k_RotAxis (i_RotAxis) l_RotAxis RotAngle(in degree) h_GBPlane k_GBPlane (i_GBPlane) l_GBPlane CrystalName(has to be defined in /data/Crystal/) Rationalize nx ny nz" << endl;
		cerr << "The i Miller indexes (for rotation axis and GB plane) should only be used if the crystal is hexagonal" << endl;
		cerr << "Rationalize can be either 0 or 1" << endl;
		cerr << "1 => rationalize the GB (i.e. search the closest CSL GB to the provided parameters, in this case the fixed parameters for CSL calculation can be important, they are read from FixedParameters.ath file if exist if not defaults values are used)" << endl;
		cerr << "0 => do not rationalize the GB" << endl;
		cerr << "nx ny nz specify the number of sampling points along each of the three CSL lattice vectors" << endl;
		cerr << "The program will generate nx x ny x nz configurations with different relative shifts between grains." << endl;
		cerr << "Each configuration is saved as a separate dump file (e.g., GB_CSL_Shift_0_0_0.lmp, etc)." << endl;
		cerr << "In addition the program will also return 3 dump files containing the CSL lattice and the two grains" << endl;
		return EXIT_FAILURE;
	}
	Dis.Printer_CreateGB();
	int h_a, k_a ,l_a, h_p, k_p, l_p;
	double theta;
	unsigned int rat;
	unsigned int nx, ny, nz;
	bool rat_b;
	string crystalName;
	if( argc == 13 ){
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
		if( rat == 0 ) rat_b = false;
		else rat_b = true;
		// 
		istringstream iss_nx(argv[10]);
		iss_nx >> nx;
		istringstream iss_ny(argv[11]); 
		iss_ny >> ny;
		istringstream iss_nz(argv[12]); 
		iss_nz >> nz;
		//
	}else{
		int i_a, i_p;
		istringstream iss_ha(argv[1]);
		iss_ha >> h_a;
		istringstream iss_ka(argv[2]);
		iss_ka >> k_a;
		istringstream iss_ia(argv[3]);
		iss_ia >> i_a;
		if( i_a != (-h_a-k_a) ) cout << "Warning, i index of rotation axis is different than -h-k, i will be considered as equal to " << -h_a-k_a << endl;
		istringstream iss_la(argv[4]);
		iss_la >> l_a;
		istringstream iss_theta(argv[5]);
		iss_theta >> theta;
		theta *= M_PI/180.;
		istringstream iss_hp(argv[6]);
		iss_hp >> h_p;
		istringstream iss_kp(argv[7]);
		iss_kp >> k_p;
		istringstream iss_ip(argv[8]);
		iss_ip >> i_p;
		if( i_p != (-h_p-k_p) ) cout << "Warning, i index of GB plane is different than -h-k, i will be considered as equal to " << -h_p-k_p << endl;
		istringstream iss_lp(argv[9]);
		iss_lp >> l_p;
		istringstream iss_cn(argv[10]);
		iss_cn >> crystalName;
		istringstream iss_rat(argv[11]);
		iss_rat >> rat;
		if( rat == 0 ) rat_b = false;
		else rat_b = true;
		istringstream iss_nx(argv[12]);
		iss_nx >> nx;
		istringstream iss_ny(argv[13]); 
		iss_ny >> ny;
		istringstream iss_nz(argv[14]); 
		iss_nz >> nz;
	}
	Bicrystal MyGB(crystalName,h_a,k_a,l_a,theta,h_p,k_p,l_p,rat_b);
	// 
	MyGB.ShiftGrainsAlongCSL(nx, ny, nz);

	MyGB.print_lmp("GB.lmp");
	MyGB.printCSL("CSL.lmp");
	MyGB.print_Grains();
	Dis.ExecutionTime();	
	return 0;
}
