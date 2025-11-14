//**********************************************************************************
//*   CreateGBs/main.cpp                                                           *
//**********************************************************************************
//* This file contains the implementation of the CreateGBs test.                   *
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
#include <Displays.h>

using namespace std;

int main(int argc, char *argv[])
{
	Displays Dis;
	Dis.Logo();
	vector<string> Properties;
	Properties.push_back("TOL_ORTHO_BOX 0.5");
	Properties.push_back("TOL_ORTHO_BOX_Z 0.5");
	Properties.push_back("MIN_BOX_HEIGHT 15.");
	Properties.push_back("MIN_BOX_ASIDE 3.");
	Properties.push_back("NB_MAX_LC 50");
	Properties.push_back("MOTIF_SHIFT 0 0 0");
	Properties.push_back("MAX_MISFIT 0.02");
	Properties.push_back("MAX_DUP 100");
	Properties.push_back("GB_SPACE 2.");
	Properties.push_back("");
	Bicrystal MyGB1("Forsterite",1,0,0,-60.8*M_PI/180.,0,1,1,0,Properties);
	MyGB1.print_lmp("GB1.lmp");
	Dis.ExecutionTime();	
	Bicrystal MyGB2("Ice",0,0,1,141.787*M_PI/180.,1,3,0,0,Properties);
	MyGB2.print_lmp("GB2.lmp");
	Dis.ExecutionTime();	
	return 0;
}
