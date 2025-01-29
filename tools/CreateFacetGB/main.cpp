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
	if( argc < 14 ){
		cerr << "Usage: CreateGB h_RotAxis k_RotAxis l_RotAxis RotAngle(in degree) h_GBPlane k_GBPlane l_GBPlane h_facet1 k_facet1 l_facet1 h_facet2 k_facet2 l_facet2 Nfacet lCrystalName" << endl;
		cerr << "This executable creates a faceted GB with two type of facets which have to have a null x normal plane component" << endl;
		cerr << "warning here h_facet1(2) k_facet1(2) l_facet1(2) represent directions (and not plane normals)" << endl;
		cerr << "Nfacet is an integer used to control the length of the facets" << endl;
		return EXIT_FAILURE;
	}
	int h_a, k_a ,l_a, h_p, k_p, l_p, hf1, kf1, lf1, hf2, kf2, lf2;
	unsigned int N;
	double theta;
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
	istringstream iss_hf1(argv[8]);
	iss_hf1 >> hf1;
	istringstream iss_kf1(argv[9]);
	iss_kf1 >> kf1;
	istringstream iss_lf1(argv[10]);
	iss_lf1 >> lf1;
	istringstream iss_hf2(argv[11]);
	iss_hf2 >> hf2;
	istringstream iss_kf2(argv[12]);
	iss_kf2 >> kf2;
	istringstream iss_lf2(argv[13]);
	iss_lf2 >> lf2;
	istringstream iss_N(argv[14]);
	iss_N >> N;
	istringstream iss_cn(argv[15]);
	iss_cn >> crystalName;
	vector<int> FacetType;
	FacetType.push_back(hf1);
	FacetType.push_back(kf1);
	FacetType.push_back(lf1);
	FacetType.push_back(hf2);
	FacetType.push_back(kf2);
	FacetType.push_back(lf2);
	Bicrystal MyGB(crystalName,h_a,k_a,l_a,theta,h_p,k_p,l_p, FacetType, N);
	MyGB.print_lmp("GB.lmp");
	MyGB.print_Grains();
	Dis.ExecutionTime();	
	return 0;
}
