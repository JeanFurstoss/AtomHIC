//**********************************************************************************
//*   SteinhardtDescriptors.h                                                                     *
//**********************************************************************************
//* This file contains the declaration of the SteinhardtDescriptors class (herited *
//* from Descriptors) and allowing to compute Steinhardt parameters                *
//* This class mainly contains the methods for computing the different versions of *
//* the Steinhardt parameters (can be set from FixedParameters):							   *
//*	- the original ones (Steinhardt et al., 1983, Phys. Rev. B)		   *
//*	- the averaged ones (Lechner and Dellago, 2008, J. Chem. Phys.) with       *
//* different type of averaging							   *
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


#ifndef STEINHARDTDESCRIPTORS_H
#define STEINHARDTDESCRIPTORS_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "Descriptors.h"

class ATOMHIC_EXPORT SteinhardtDescriptors : public Descriptors {
private:
	// variables for calculation
	unsigned int nbNMax; // maximum number of neighbours of the AtomicSystem (used to access to the different arrays)
	unsigned int *Malpha;// array containing the index of neighbours of the same species (or same site in case of multisite crystal) with the first line corresponding to the number of neighbours, i.e. Malpha[i*(nbNMax+1)] = nb of neighbour of atom i, Malpha[i*(nbNMax+1)+j+1] = id of the jth neighbour of atom i
	double *Calpha; // normalization factor, only used for OneL mode (BondOriParam) 
	std::complex<double> *Qlm; // complex array containing the spherical harmonic for the different modes Qlm[i*(l_sph*2+1)*(l_sph+1)+l*(l_sph*2+1)+m] gives the spherical harmonic for atom i and degree l and m (a bit different for OneL mode see function)
	unsigned int lsph2; // (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1; // l_sph*2+1;

	// Variables for calculation and printing
	double zeronum = 1e-8;
	const int bar_length = 30;
	double prog;
	unsigned int count_t;
	
	// Properties of Steinhardt descriptors (FixedParameters)
	double rc = 5.; // cutoff radius for neighbor research
	int l_sph = 10; // harmonic degree
	std::string SteinhardtStyle = "Mono"; // multi, mono
	std::string AverageStyle = "Multi"; // multi mono
	std::string mode = "Full"; // OneL Full
public:
	// constructors
	SteinhardtDescriptors(AtomicSystem *_MySystem);
	SteinhardtDescriptors(AtomicSystem *_MySystem, std::vector<std::string> _Properties);
	// methods
	void readProperties(std::vector<std::string> _Properties);
	void ComputeDescriptors();
	void setProperties();
	void readFixedParams();
	void InitializeArrays();
	void ComputeBondOriParam();
	void ComputeSteinhardtParameters_Mono();
	void ComputeSteinhardtParameters_Multi();
	void AverageSteinhardtParameters_Mono();
	void AverageSteinhardtParameters_Multi();
	void printDescriptorsPropToDatabase(std::ofstream &writefile);
	// getters
	int get_l_sph(){ return l_sph; }
	double get_rc(){ return rc; }
	std::string getMode(){ return mode; }
	// destructor
	~SteinhardtDescriptors();
	
};

#endif
