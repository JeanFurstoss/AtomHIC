//**********************************************************************************
//*   SampleGBShiftDoF/main.cpp			                                   *
//**********************************************************************************
//* This file contains the implementation of the SampleGBShiftDoF		   *
//* executable.							                   *
//* It allows to create atomic systems containing a GB with a given misorientation *
//* and GB plane with shift between the two crystals distributed along the unit    *
//* cell plane of the GB (gamma-surface), or along DSC basis, or along CSL basis   *
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
#include <Displays.h>

using namespace std;

void ExecMsg(){
	cerr << "Usage: SampleGBShiftDoF SamplingMode DiscretizationMode nx ny (nz) h_RotAxis k_RotAxis (i_RotAxis) l_RotAxis RotAngle(in degree) h_GBPlane k_GBPlane (i_GBPlane) l_GBPlane CrystalName(has to be defined in /data/Crystal/) Rationalize VacuumOrNot" << endl;
	cerr << "This executable creates atomic systems containing a GB with a given misorientation and GB plane and varying the translational (i.e. shift) degree of freedom" << endl << endl;
	cerr << "SamplingMode can be \"InPlane\", \"DSC\" or \"CSL\":" << endl;
	cerr << "\t - InPlane => the shift between the two crystals distributed along the unit cell plane of the GB (gamma-surface like sampling)" << endl;
	cerr << "\t - DSC => the shift between the two crystals distributed along DSC basis (3D sampling)" << endl;
	cerr << "\t - CSL => the shift between the two crystals distributed along CSL basis (3D sampling)" << endl;
	cerr << "For the two last modes (DSC or CSL) the nz parameter must be provided (as the sampling is performed in 3D)" << endl << endl;
	cerr << "DiscretizationMode can be \"NbPts\" or \"Distance\":" << endl;
	cerr << "\t - NbPts => nx, ny (nz) specify the number of shifts (sampling points) along each of the two (or three) sampled vectors" << endl;
	cerr << "\t - Density =>  nx, ny (nz) specify the distances between two shifts in each direction (then the number of shifts will be computed from these distances)" << endl << endl;
	cerr << "The i Miller indexes (for rotation axis and GB plane) should only be used if the crystal is hexagonal" << endl;
	cerr << "Rationalize can be either 0 or 1" << endl;
	cerr << "1 => rationalize the GB (i.e. search the closest CSL GB to the provided parameters, in this case the fixed parameters for CSL calculation can be important, they are read from FixedParameters.ath file if exist if not defaults values are used)" << endl;
	cerr << "0 => do not rationalize the GB" << endl;
	cerr << "VacuumOrNot could be 0 or 1, 0 => no vaccum is added to above and bellow the GB (i.e. the periodic BC will lead to 2 GBs in the system), 1 => add vacuum and center the GB in the middle of the cell" << endl;
	cerr << "Each configuration is saved as a separate dump file (e.g., GB_Shift_0_0.lmp or GB_DSC_Shift_0_0_0.lmp or GB_CSL_Shift_0_0_0.lmp, etc), also containing the values of the applied shift" << endl;
	cerr << "In addition the program will also return 3 dump files containing the CSL lattice and the two grains (the two latters can be used to change shift between crystals), if vacuum = 1 Grain1_vacuum.lmp and Grain2_vacuum.lmp will be also generated" << endl;
	cerr << "When SamplingMode is InPlane, this program enforce the fixed parameter MIN_BOX_ASIDE to a value of 2A to be sure that the created GB has only one unit cell" << endl;
	exit(EXIT_FAILURE);
}

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc < 3 ) ExecMsg();
	string SampMode = argv[1];
        string DisMode = argv[2];
	if( SampMode != "InPlane" && SampMode != "DSC" && SampMode != "CSL" ){
		cout << "\t *** Error : Unknown SamplingMode !!! ***" << endl << endl;
		ExecMsg();
	}
	if( DisMode != "NbPts" && DisMode != "Distance" ){
		cout << "\t *** Error : Unknown DiscretizationMode !!! ***" << endl << endl;
		ExecMsg();
	}
	if( !( SampMode == "InPlane" && ( argc == 15 || argc == 17 ) ) && !( SampMode != "InPlane" && ( argc == 16 || argc == 18 ) ) ){
		cout << "\t *** Error : Wrong number of arguments !!! ***" << endl << endl;
		ExecMsg();
	}
	if( SampMode == "InPlane" ) Dis.Printer_SampleGB_GammaSurface();
	else Dis.Printer_CreateGB();
	int h_a, k_a ,l_a, h_p, k_p, l_p, i_a, i_p;
	double theta;
	unsigned int rat, vac;
	unsigned int nx, ny, nz;
	double dx, dy, dz;
	bool rat_b, vac_b;
	bool isDist = false;
	string crystalName;
	unsigned int current_readind = 3;
	if( DisMode == "NbPts" ){
		istringstream iss_nx(argv[current_readind]);
		iss_nx >> nx;
		current_readind++;
		istringstream iss_ny(argv[current_readind]); 
		iss_ny >> ny;
		current_readind++;
		if( SampMode != "InPlane" ){
			istringstream iss_nz(argv[current_readind]); 
			iss_nz >> nz;
			current_readind++;
		}
	}else{
		istringstream iss_dx(argv[current_readind]);
		iss_dx >> dx;
		current_readind++;
		istringstream iss_dy(argv[current_readind]); 
		iss_dy >> dy;
		current_readind++;
		if( SampMode != "InPlane" ){
			istringstream iss_dz(argv[current_readind]); 
			iss_dz >> dz;
			current_readind++;
		}
	}

	istringstream iss_ha(argv[current_readind]);
	iss_ha >> h_a;
	current_readind++;
	istringstream iss_ka(argv[current_readind]);
	iss_ka >> k_a;
	current_readind++;
	if( argc == 17 || argc == 18 ){
		istringstream iss_ia(argv[current_readind]);
		iss_ia >> i_a;
		current_readind++;
		if( i_a != (-h_a-k_a) ) cout << "Warning, i index of rotation axis is different than -h-k, i will be considered as equal to " << -h_a-k_a << endl;
	}
	istringstream iss_la(argv[current_readind]);
	iss_la >> l_a;
	current_readind++;
	istringstream iss_theta(argv[current_readind]);
	iss_theta >> theta;
	current_readind++;
	theta *= M_PI/180.;
	istringstream iss_hp(argv[current_readind]);
	iss_hp >> h_p;
	current_readind++;
	istringstream iss_kp(argv[current_readind]);
	iss_kp >> k_p;
	current_readind++;
	if( argc == 17 || argc == 18 ){
		istringstream iss_ip(argv[current_readind]);
		iss_ip >> i_p;
		current_readind++;
		if( i_p != (-h_p-k_p) ) cout << "Warning, i index of GB plane is different than -h-k, i will be considered as equal to " << -h_p-k_p << endl;
	}
	istringstream iss_lp(argv[current_readind]);
	iss_lp >> l_p;
	current_readind++;
	istringstream iss_cn(argv[current_readind]);
	iss_cn >> crystalName;
	current_readind++;
	istringstream iss_rat(argv[current_readind]);
	iss_rat >> rat;
	current_readind++;
	if( rat == 0 ) rat_b = false;
	else rat_b = true;
	istringstream iss_vac(argv[current_readind]); 
	iss_vac >> vac;
	current_readind++;
	if( vac == 0 ) vac_b = false;
	else vac_b = true;
	
	vector<string> Properties;
	
	// set minimum aside to very small value to be sure that we only have the UC of the GB
	if( SampMode == "InPlane" ) Properties.push_back("MIN_BOX_ASIDE 2.");

	Bicrystal MyGB(crystalName,h_a,k_a,l_a,theta,h_p,k_p,l_p,rat_b,Properties);
	//
	if( DisMode == "Distance" ){
		double Lx = MyGB.getH1()[0];
		double Ly = MyGB.getH2()[1];
		double Lz = MyGB.getH3()[2];
		nx = round(Lx/dx);
		ny = round(Ly/dy);
		nz = round(Lz/dz);
		if( nx < 1 ){
			cout << "Warning from the provided dx value nx is null => nx has been set equal to 1" << endl;
			nx = 1;
		}
		if( ny < 1 ){
			cout << "Warning from the provided dy value ny is null => ny has been set equal to 1" << endl;
			ny = 1;
		}
		if( nz < 1 ){
			cout << "Warning from the provided dz value nz is null => nz has been set equal to 1" << endl;
			nz = 1;
		}
	}
		
	if( SampMode == "InPlane" ) MyGB.ShiftGrainsAlongUCInPlane(nx, ny, vac_b);
	else if( SampMode == "DSC" ){
		MyGB.printDSC();
		MyGB.ShiftGrainsAlongDSC(nx, ny, nz, vac_b); 
	}else if( SampMode == "CSL" ) MyGB.ShiftGrainsAlongCSL(nx, ny, nz, vac_b); 

	MyGB.print_lmp("GB.lmp");
	MyGB.printCSL("CSL.lmp");
	MyGB.print_Grains(vac_b);

	Dis.ExecutionTime();	
	return 0;
}
