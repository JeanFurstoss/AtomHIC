//**********************************************************************************
//*   AtomicSystemTrajectory.h                                                     *
//**********************************************************************************
//* This file contains the declaration of the AtomicSystemTrajectory class         *
//* It allows to read and treat atomic system trajectories such as .lammpstrj files*
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
//*	- 						                           *
//**********************************************************************************


#ifndef ATOMICSYSTEMTRAJECTORY_H
#define ATOMICSYSTEMTRAJECTORY_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "Bicrystal.h"
#include "MathTools.h"
#include <string>
#include <vector>

class ATOMHIC_EXPORT AtomicSystemTrajectory {
private:
	// Properties
	std::string input_filename;
	bool IsTraj = false;
	std::vector<AtomicSystem*> AtSysList;
	std::vector<Bicrystal*> BicrystalList;
	unsigned int nbSys = 0;
	std::vector<unsigned int> timesteps;
public:
	// constructors
	AtomicSystemTrajectory();
	// methods
	bool SearchIsTrajectory(const std::string &_input_filename);
	void setAtomicSystemList(const std::string &_input_filename);
	void printSystem_aux(const std::string& filename, const std::string& AuxId);
	// getters
	unsigned int getNbSys(){ return nbSys; }
	unsigned int getTimestep(const unsigned int &i){ return timesteps[i]; }
	AtomicSystem* getAtomicSystem(const unsigned int &i){ return AtSysList[i]; } 
	// destructor
	~AtomicSystemTrajectory();
	
};

#endif
