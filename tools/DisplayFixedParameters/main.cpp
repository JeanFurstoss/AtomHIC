//**********************************************************************************
//*   DisplayFixedParameters/main.cpp                                              *
//**********************************************************************************
//* This file contains the implementation of the DisplayFixedParameters executable *
//* It allows to display all the fixed parameters used in AtomHIC by reading the   *
//* parameters and their description in /data/FixedParameters/FixedParameters.ath  *
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
	Dis.PrintFixedParameters();
	Dis.ExecutionTime();	
	return 0;
}
