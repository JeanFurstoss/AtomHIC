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
		cerr << "Usage: SampleGBComplexionDSC h_RotAxis k_RotAxis (i_RotAxis) l_RotAxis RotAngle(in degree) h_GBPlane k_GBPlane (i_GBPlane) l_GBPlane CrystalName(has to be defined in /data/Crystal/) Rationalize nx ny nz" << endl;
		cerr << "The i Miller indexes (for rotation axis and GB plane) should only be used if the crystal is hexagonal" << endl;
		cerr << "Rationalize can be either 0 or 1" << endl;
		cerr << "1 => rationalize the GB (i.e. search the closest CSL GB to the provided parameters, in this case the parameters for CSL calculation in /data/FixedParameters/FixedParameters.dat can be important)" << endl;
		cerr << "0 => do not rationalize the GB, in this case the parameters for constructing crystal and bicrystal in /data/FixedParameters/FixedParameters.dat can be important" << endl;
		cerr << "       In that case, the grain construction uses fixed input parameters." << endl;
		//
	    cerr << "nx, ny, nz specify the number of shifts (sampling points) along each of the three DSC lattice vectors." << endl;
    	cerr << "  These values control how many configurations of the GB are generated with different relative vectors." << endl;
    	cerr << "  For each shift, Grain 1 is translated by a linear combination of the DSC basis vectors, wraped into the box,and pasted on top of Grain 2 along the z-direction (normal to the GB plane)." << endl;
		cerr << "  Each configuration is saved as a separate dump file (e.g., GB_DSC_Shift_0_0_0.lmp, etc)." << endl;
		//
		cerr << "In addition the program will also return 3 dump files containing the CSL lattice and the two grains (the two latters can be used to change shift between crystals)" << endl;
		return EXIT_FAILURE;
	}
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
		istringstream iss_lp(argv[9]);
		iss_lp >> l_p;
		istringstream iss_cn(argv[10]);
		iss_cn >> crystalName;
		istringstream iss_rat(argv[11]);
		iss_rat >> rat;
		istringstream iss_nx(argv[12]);
		iss_nx >> nx;
		istringstream iss_ny(argv[13]); 
		iss_ny >> ny;
		istringstream iss_nz(argv[14]); 
		iss_nz >> nz;
		if( rat == 0 ) rat_b = false;
		else rat_b = true;
	}
	Bicrystal MyGB(crystalName,h_a,k_a,l_a,theta,h_p,k_p,l_p,rat_b);
	// 
    //MyGB.ShiftGrainsAlongDSC(2, 2, 1);  // e.g.
	MyGB.ShiftGrainsAlongDSC(nx, ny, nz);

	MyGB.print_lmp("GB.lmp");
	MyGB.printCSL("CSL.lmp");
	MyGB.print_Grains();
	//
	if (MyGB.getIsDSC_Basis()) {
		MyGB.printDSC();
	}
	Dis.ExecutionTime();	
	return 0;
}
