//**********************************************************************************
//*   ACEDescriptors.h                                                                     *
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


#ifndef ACEDESCRIPTORS_H
#define ACEDESCRIPTORS_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "Descriptors.h"
#include "ML-PACE/ace/ace_b_basis.h"

class ATOMHIC_EXPORT ACEDescriptors : public Descriptors {
private:
	// Properties of ACE descriptors
	double rc; // cutoff radius for neighbor research
	int l_sph; // harmonic degree
	std::string yaml_filename;
	ACEBBasisSet *MyACEBBase;
	    Array2D<DOUBLE_TYPE> A_rank1 = Array2D<DOUBLE_TYPE>(
            "A_rank1"); ///< 2D-array for storing A's for rank=1, shape: A(mu_j,n)
    Array4DLM<ACEComplex> A = Array4DLM<ACEComplex>(
            "A"); ///< 4D array with (l,m) last indices  for storing A's for rank>1: A(mu_j, n, l, m)
    double Y00 = 1.;
	//BBasisConfiguration MyBBasisConf;
	// from what I understood, we need to provide:
	// 0. mu the number of chemical species 
	// 1. R the rank of the expansion (R-body interaction)
	// 2. lmax the maximum angular momentum (degree of the spherical harmonics)
	// 3. nmax the number of radial functions
	// => this should gives N basis functions (descriptors) with N = ??
public:
	// constructors
	ACEDescriptors(AtomicSystem *_MySystem, std::string _yaml_filename);
	// methods
	void ComputeDescriptors();

	// destructor
	~ACEDescriptors();
	
};

#endif
