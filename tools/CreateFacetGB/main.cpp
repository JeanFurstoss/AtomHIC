//**********************************************************************************
//*   CreateFacetGB/main.cpp                                                       *
//**********************************************************************************
//* This file contains the implementation of the CreateFacetGB executable.         *
//* It allows to create an atomic system containing a faceted GB with a given  	   *
//* misorientation and GB facet planes						   *
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
#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "ComputeAuxiliary.h"
#include <vector>
#include <Displays.h>

using namespace std;

// TODO haxegonal case
int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc != 19 && argc != 16 ){
		cerr << "Usage: CreateFacetGB h_RotAxis k_RotAxis l_RotAxis RotAngle(in degree) h_GBPlane k_GBPlane l_GBPlane (hx_Plane kx_Plane lx_Plane) h_facet1 k_facet1 l_facet1 h_facet2 k_facet2 l_facet2 Lfacet CrystalName" << endl;
		cerr << "This executable creates a faceted GB with two type of facets by providing the geometrical parameters of the GB (rotation axis and angle and GB plane normal)" << endl;
		cerr << "h_GBPlane, k_GBPlane, l_GBPlane are the Miller indices of the GB plane in the crystal reference frame of the lower grain" << endl;
		cerr << "The Miller indices of the facet planes (h_facet1, k_facet1, ..) are expressed in the crystal reference frame of the lower grain" << endl;
		cerr << "Lfacet is the desired length of the first facet (the second facet length is automatically computed in order to respect periodic boundary conditions (PBC))" << endl;
		cerr << "If the obtained facet length is too far from the provided one, the parameters MAX_VAR_FACET_LENGTh and MAX_DUP_FACET can be tuned" << endl;
		cerr << "CrystalName correspond to the crystal in which the system will be created, this crystal should be defined in /data/Crystal/" << endl;
		cerr << "If h_x, k_x and l_x parameters are provided, the (h_x k_x l_x) crystallographic plane of the lower grain will be oriented along the x axis. In this case the (hkl) and (hxkxlx) planes should be orthogonal" << endl;
		cerr << "If h_x k_x and l_x are not provided (or if the two planes are not orthogonal) the crystal orientation along the x axis will be automatically computed" << endl << endl;
		cerr << "If the provided facet types do not allows to create a system respecting PBC, the program will propose to compute the possible facet combination for the given GB" << endl;
		cerr << "The corresponding Miller indices of the GB and facet plane in the reference frame of the upper grain will be automatically computed by the program (if it does not work the MAX_HKL_SEARCH parameter can be tuned)" << endl;
		Dis.Printer_CreateFacetGB();
		return EXIT_FAILURE;
	}
	Dis.Printer_CreateFacetGB();
	int h_a, k_a ,l_a, h_p, k_p, l_p, hf1, kf1, lf1, hf2, kf2, lf2, hx_p(0), kx_p(0), lx_p(0);
	double theta, Lfacet;
	string crystalName;
	unsigned int current_readind = 1;
	istringstream iss_ha(argv[current_readind]);
	iss_ha >> h_a;
	current_readind++;
	istringstream iss_ka(argv[current_readind]);
	iss_ka >> k_a;
	current_readind++;
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
	istringstream iss_lp(argv[current_readind]);
	iss_lp >> l_p;
	current_readind++;
	if( argc == 19 ){
		istringstream iss_hpx(argv[current_readind]);
		iss_hpx >> hx_p;
		current_readind++;
		istringstream iss_kpx(argv[current_readind]);
		iss_kpx >> kx_p;
		current_readind++;
		istringstream iss_lpx(argv[current_readind]);
		iss_lpx >> lx_p;
		current_readind++;
	}
	istringstream iss_hf1(argv[current_readind]);
	iss_hf1 >> hf1;
	current_readind++;
	istringstream iss_kf1(argv[current_readind]);
	iss_kf1 >> kf1;
	current_readind++;
	istringstream iss_lf1(argv[current_readind]);
	iss_lf1 >> lf1;
	current_readind++;
	istringstream iss_hf2(argv[current_readind]);
	iss_hf2 >> hf2;
	current_readind++;
	istringstream iss_kf2(argv[current_readind]);
	iss_kf2 >> kf2;
	current_readind++;
	istringstream iss_lf2(argv[current_readind]);
	iss_lf2 >> lf2;
	current_readind++;
	istringstream iss_L(argv[current_readind]);
	iss_L >> Lfacet;
	current_readind++;
	istringstream iss_cn(argv[current_readind]);
	iss_cn >> crystalName;
	current_readind++;
	vector<int> FacetType;
	FacetType.push_back(hf1);
	FacetType.push_back(kf1);
	FacetType.push_back(lf1);
	FacetType.push_back(hf2);
	FacetType.push_back(kf2);
	FacetType.push_back(lf2);
	vector<string> Properties;
	Bicrystal MyGB(crystalName,h_a,k_a,l_a,theta,h_p,k_p,l_p, FacetType, Lfacet, Properties, hx_p, kx_p, lx_p);
	MyGB.print_lmp("GB.lmp");
	MyGB.print_Grains();
	Dis.ExecutionTime();	
	return 0;
}
