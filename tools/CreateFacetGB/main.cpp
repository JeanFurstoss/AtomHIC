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

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	if( argc != 19 && argc != 16 ){
		cerr << "Usage: CreateGB h_RotAxis k_RotAxis l_RotAxis RotAngle(in degree) h_GBPlane k_GBPlane l_GBPlane hx_GBPlane kx_GBPlane lx_GBPlane h_facet1 k_facet1 l_facet1 h_facet2 k_facet2 l_facet2 Nfacet lCrystalName" << endl;
		cerr << "This executable creates a faceted GB with two type of facets which have to have a null x normal plane component" << endl;
		cerr << "warning here h_facet1(2) k_facet1(2) l_facet1(2) represent directions (and not plane normals)" << endl;
		cerr << "Nfacet is an integer used to control the length of the facets" << endl;
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
