//**********************************************************************************
//*   ComputeAuxiliary.h                                                           *
//**********************************************************************************
//* This file contains the declaration of the ComputeAuxiliary class which is used *
//* to compute some auxiliary properties of an AtomicSystem. For the moment, the   *
//* possible auxiliary properties are:						   *
//*	- a bond orientational parameter (Furstoss et al., 2024, Comp. Mat. Sc.)   *
//*	- the crystallographic site index (based on the BondOriParam)		   *
//*	- the affine transformation matrix, D2Min and atomic strain (based on the  *
//* work of (Shimizu et al., 2007, Mater. Trans.))				   *
//*	- the strain tensor based either or atom numerotation of crystallographic  *
//* site
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
//*	-                                                                          *
//**********************************************************************************

#ifndef COMPUTEAUXILIARY_H
#define COMPUTEAUXILIARY_H

#include "SteinhardtDescriptors.h"
#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>

class ATOMHIC_EXPORT ComputeAuxiliary {
protected:
	AtomicSystem *_MySystem;
	MathTools *MT;
	SteinhardtDescriptors *StDes;
	bool IsSteinhardtDescriptor = false;
	unsigned int *Atom_SiteIndex; // determined using comparison of Steinhardt param for l=l_sph (using the database)
	bool IsAtomSiteIndex = false;
	double *StrainTensor;
	double *Strain_invII;

	bool IsStrainTensor = false;
	bool IsStrainInvII = false;
	// For atomic strain and D2Min
	double *Ji; // affine transformation matrix
	double *d0; // reference state variable for affine transfo matrix
	double *current_d; // current delta for affine transfo matrix
	double *Vi_inv; // reference state variable for affine transfo matrix
	double *AtomicStrain;
	double *D2Min;
	bool Reference_AtomicStrain_Computed = false;
	bool IsAtomicStrain = false;
	bool IsJi = false;
	bool IsD2Min = false;
	// Parameters to read
	std::string FixedParam_Filename = "Fixed_Parameters.dat";
	double tolSites;
public:
	// constructors
	ComputeAuxiliary(){};
	ComputeAuxiliary(AtomicSystem *_MySystem): _MySystem(_MySystem){
		this->MT = new MathTools;
		read_params();
	};
	// getters
	// methods
	void ComputeAtomSiteIndex();
	double* BondOrientationalParameter();
	double* Compute_AffineTransfoMatrix(AtomicSystem &ReferenceSystem, double rc); // Ji matrix as defined in ovito (i.e. Shimizu, Ogata, Li: Mater. Trans. 48 (2007), 2923)
	double* Compute_AtomicStrain(AtomicSystem &ReferenceSystem, double rc); // atomic strain as defined in ovito (i.e. Shimizu, Ogata, Li: Mater. Trans. 48 (2007), 2923)
	double* Compute_D2Min(AtomicSystem &ReferenceSystem, double rc); // D2Min as defined in Delbecq et al. 2023
	double* Compute_StrainTensor();
	double* Compute_StrainTensor(unsigned int FromNum);
	double* Compute_StrainTensor_invII();
	unsigned int *get_AtomSiteIndex(){ return this->Atom_SiteIndex; }
	double *get_StrainTensor(){ return this->StrainTensor; }
	double *get_StrainInvII(){ return this->Strain_invII; }
	void read_params();
	// destructor
	~ComputeAuxiliary();
	
};

#endif
