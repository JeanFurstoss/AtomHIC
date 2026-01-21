//**********************************************************************************
//*   CreateOrientedPlane/main.cpp                                                 *
//**********************************************************************************
//* This file contains the implementation of the CreateOrientedPlane executable.   *
//* It allows to create an orthogonal atomic system containing a perfect crystal   *
//* z-oriented to a given plane							   *
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

#include <AtomicSystem.h>
#include <Crystal.h>
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
	if( argc != 6 && argc != 7 && argc != 9 && argc != 10 ){
		cerr << "Usage: CreateOrientedPlane h k (i) l (h_x k_x (i_x) l_x) CrystalName(has to be defined in /data/Crystal/) OutputFilename" << endl;
		cerr << "This executable will create an orthogonal cell with a z-oriented (hkl) plane" << endl;
		cerr << "If h_x, k_x (i_x) and l_x parameters are provided, the (h_x k_x l_x) crystallographic plane will be oriented along the x axis. In this case the (hkl) and (hxkxlx) planes should be orthogonal" << endl;
		cerr << "If h_x k_x and l_x are not provided (or if the two planes are not orthogonal) the crystal orientation along the x axis will be automatically computed" << endl;
		cerr << "Index i should be only used for hexagonal crystals" << endl;
		Dis.Printer_CreateOrientedPlane();
		return EXIT_FAILURE;
	}
	Dis.Printer_CreateOrientedPlane();
	int h, k ,l, i;
	int h_x(0), k_x(0), l_x(0), i_x;
	string crystalType, filename;
	unsigned int current_readind = 1;
	// reading parameters
	istringstream iss_h(argv[current_readind]);
	iss_h >> h;
	current_readind++;
	istringstream iss_k(argv[current_readind]);
	iss_k >> k;
	current_readind++;
	if( argc == 7 || argc == 10 ){
		istringstream iss_i(argv[current_readind]);
		iss_i >> i;
		current_readind++;
		if( i != (-h-k) ) cout << "Warning, i index for plane along z axis is different than -h-k, i will be considered as equal to " << -h-k << endl;
	}
	istringstream iss_l(argv[current_readind]);
	iss_l >> l;
	current_readind++;
	if( argc == 9 || argc == 10 ){
		istringstream iss_h_x(argv[current_readind]);
		iss_h_x >> h_x;
		current_readind++;
		istringstream iss_k_x(argv[current_readind]);
		iss_k_x >> k_x;
		current_readind++;
		if( argc == 10 ){
			istringstream iss_i_x(argv[current_readind]);
			iss_i_x >> i_x;
			current_readind++;
			if( i_x != (-h_x-k_x) ) cout << "Warning, i index for plane along x axis is different than -h-k, i will be considered as equal to " << -h_x-k_x << endl;
		}
		istringstream iss_l_x(argv[current_readind]);
		iss_l_x >> l_x;
		current_readind++;
	}
	crystalType=argv[current_readind];
	iss_l >> l;
	current_readind++;
	filename=argv[current_readind];
	iss_l >> l;
	current_readind++;

	string i_str = "";
	if( argc == 7 ) i_str = to_string(-h-k)+"_";
	cout << "Creating orthogonal cell of " << crystalType << " crystal, with (" << h << "_" << k << "_" << i_str << l << ") plane normal aligned with z direction" << endl;
	string i_x_str = "";
	if( argc == 10 ){
		i_x_str = to_string(-h_x-k_x)+"_";
		cout << "and with (" << h_x << "_" << k_x << "_" << i_x_str << l_x << ") plane normal aligned with x direction" << endl;
	}
	cout << endl;
	Crystal MyCrystal(crystalType);
	MyCrystal.RotateCrystal(h,k,l,h_x,k_x,l_x);
	MyCrystal.ConstructOrthogonalCell();
	cout << "Building orthogonal cell required deforming the crystal such that:" << endl;
	vector<vector<string>> arr_element = {
		{"X","Y","Z","Shear"},
		{to_string(MyCrystal.GetCrystalDef()[0])+" %",to_string(MyCrystal.GetCrystalDef()[1])+" %",to_string(MyCrystal.GetCrystalDef()[2])+" %",to_string(MyCrystal.GetCrystalDef()[3])+" %"}
	};
	vector<vector<unsigned int>> arr_fusion = {{1,1,1,1},{1,1,1,1}};
	Dis.DisplayArray(arr_element,arr_fusion);
	cout << "If these values are too large it may cause issues during the relaxation of the cell" << endl;
	cout << "These values can be reduced by decreasing TOL_ORTHO_BOX and TOL_ORTHO_BOX_Z in data/FixedParameters/FixedParameters.dat" << endl;
	cout << endl;

	MyCrystal.ComputeOrthogonalPlanesAndDirections();
	Dis.DisplayOrthogonalCell(&MyCrystal);
	MyCrystal.getOrientedSystem()->print_lmp(filename);
	Dis.ExecutionTime();	
	return 0;
}
