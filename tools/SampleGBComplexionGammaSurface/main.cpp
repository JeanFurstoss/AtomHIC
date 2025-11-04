//**********************************************************************************
//*   SampleGBComplexionGammaSurface/main.cpp                                      *
//**********************************************************************************
//* This file contains the implementation of the SampleGBComplexionGammaSurface    *
//* executable.							                   *
//* It allows to create atomic systems containing a GB with a given misorientation *
//* and GB plane with shift between the two crystals distributed along the unit    *
//* cell plane of the GB							   *
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
		cerr << "Usage: SampleGBComplexionGammaSurface h_RotAxis k_RotAxis (i_RotAxis) l_RotAxis RotAngle(in degree) h_GBPlane k_GBPlane (i_GBPlane) l_GBPlane CrystalName(has to be defined in /data/Crystal/) Rationalize nx ny VacuumOrNot" << endl;
		cerr << "The i Miller indexes (for rotation axis and GB plane) should only be used if the crystal is hexagonal" << endl;
		cerr << "Rationalize can be either 0 or 1" << endl;
		cerr << "1 => rationalize the GB (i.e. search the closest CSL GB to the provided parameters, in this case the parameters for CSL calculation in /data/FixedParameters/FixedParameters.dat can be important)" << endl;
		cerr << "0 => do not rationalize the GB, in this case the parameters for constructing crystal and bicrystal in /data/FixedParameters/FixedParameters.dat can be important" << endl;
		cerr << "       In that case, the grain construction uses fixed input parameters." << endl;
		//
	    cerr << "nx, ny specify the number of shifts (sampling points) along each of the two vectors in the GB plane defining the unit cell of the GB (gamma-surface)" << endl;
    	cerr << "  These values control how many configurations of the GB are generated with different relative vectors." << endl;
	cerr << "VacuumOrNot could be 0 or 1, 0 => no vaccum is added to above and bellow the GB (i.e. the periodic BC will lead to 2 GBs in the system), 1 => add vacuum and center the GB in the middle of the cell" << endl;
		cerr << "  Each configuration is saved as a separate dump file (e.g., GB_Shift_0_0.lmp, etc), also containing the values of the applied shift" << endl;
		//
		cerr << "In addition the program will also return 3 dump files containing the CSL lattice and the two grains (the two latters can be used to change shift between crystals), if vacuum = 1 Grain1_vacuum.lmp and Grain2_vacuum.lmp will be also generated" << endl;
		return EXIT_FAILURE;
	}
	int h_a, k_a ,l_a, h_p, k_p, l_p;
	double theta;
	unsigned int rat, vac;
	unsigned int nx, ny;
	bool rat_b, vac_b;
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
		istringstream iss_vac(argv[12]); 
		iss_vac >> vac;
		if( vac == 0 ) vac_b = false;
		else vac_b = true;
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
		istringstream iss_vac(argv[14]); 
		iss_vac >> vac;
		if( rat == 0 ) rat_b = false;
		else rat_b = true;
		if( vac == 0 ) vac_b = false;
		else vac_b = true;
	}
	Bicrystal MyGB(crystalName,h_a,k_a,l_a,theta,h_p,k_p,l_p,rat_b);
	// 
	MyGB.ShiftGrainsAlongUCInPlane(nx, ny, vac_b);

	MyGB.print_lmp("GB.lmp");
	MyGB.printCSL("CSL.lmp");
	MyGB.print_Grains(vac_b);

	Dis.ExecutionTime();	
	return 0;
}
