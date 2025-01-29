//**********************************************************************************
//*   Displays.h                                                                   *
//**********************************************************************************
//* This file contains the declaration of the Displays class used to display stuff *
//* during the execution of AtomHIC
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
//*	- display of the resume of bicrystal construction			   * 						                           *
//*	- 						                           *
//**********************************************************************************

#ifndef DISPLAYS_H
#define DISPLAYS_H

#include "AtomHicExport.h"
#include <chrono>

class ATOMHIC_EXPORT Displays {
private:
	std::chrono::high_resolution_clock::time_point time_beg;
public:
	// constructors
	Displays();
	// methods
	void Logo();
	void ExecutionTime();	
	// destructor
	~Displays(){};
	
};

#endif
