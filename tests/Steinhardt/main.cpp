//**********************************************************************************
//*   Steinhardt/main.cpp                                                          *
//**********************************************************************************
//* This file contains the implementation of the Steinhardt test.                  *
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

#include <Bicrystal.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include "MathTools.h"
#include <SteinhardtDescriptors.h>

using namespace std;

int main(int argc, char *argv[])
{
	AtomicSystem MySystem("dump.lmp");
	double rc = 5.;
	int l_sph = 10;
	vector<string> Properties_1;
	Properties_1.push_back("STEINHARDT_MODE Full");
	Properties_1.push_back("CUTOFF_RADIUS 5.");
	Properties_1.push_back("STEINHARDT_STYLE Mono");
	Properties_1.push_back("AVE_STYLE Mono");
	Properties_1.push_back("NUMBER_OF_DIMENSION 10");
	vector<string> Properties_2;
	Properties_2.push_back("STEINHARDT_MODE Full");
	Properties_2.push_back("CUTOFF_RADIUS 5.");
	Properties_2.push_back("STEINHARDT_STYLE Multi");
	Properties_2.push_back("AVE_STYLE Multi");
	Properties_2.push_back("NUMBER_OF_DIMENSION 10");
	SteinhardtDescriptors St_1(&MySystem,Properties_1);
	SteinhardtDescriptors St_2(&MySystem,Properties_2);
	MySystem.setAux_vec(St_1.getDescriptors(),10,"Q_mono_mono");
	MySystem.setAux_vec(St_2.getDescriptors(),10,"Q_multi_multi");
	MySystem.printSystem_aux("output.xsf","Q_mono_mono Q_multi_multi");
	return 0;
}
