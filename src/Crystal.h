//**********************************************************************************
//*   Crystal.h                                                                    *
//**********************************************************************************
//* This file contains the declaration of the Crystal class which is used for:     *
//*	- compute reciproqual lattice
//*	- contructing orthogonal oriented AtomicSystem (with a given plane)        *
//*	- read the crystal database including the crystallography (cubic,          *
//* orthorhombic..), the atomic coordinates, types, crystallographic sites and     *
//* reference bond orientational parameters, the "do not separe" list (used for    *
//* dont separating some ions at surface)
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
//*	- extend methodology to hexagonal, monoclinic and triclinic systems (for   *
//* the latters we have to test if the current version work)                       *
//*	-                                                                          *
//**********************************************************************************

#ifndef CRYSTAL_H
#define CRYSTAL_H

#include "AtomHicExport.h"
#include "MyStructs.h"
#include <string>
#include "MathTools.h"

class AtomicSystem;

class ATOMHIC_EXPORT Crystal {
private:
	std::string name;
	std::string crystallo; // Cubic, Hexagonal, Tetragonal, Orthorhombic, Monoclinic or Triclinic
	std::string path2database;
	unsigned int nbAtom;
	unsigned int nbAtomType;
	unsigned int MaxAtomType = 15;
	std::string *AtomType; // the different atom species present in the system
	unsigned int *AtomType_uint;
	unsigned int *AtomSite;
	unsigned int *NbAtomSite; // array such as NbAtomSite[i] = number of site for atom type AtomType[i]
	bool IsCharge = false;
	bool IsMolId = false; // is molecule id (atom_style full in LAMMPS)
	bool IsBond = false; // bond between 2 atoms (e.g. O-H bonds in ice with TIP4P potential)
	unsigned int nbBonds;
	unsigned int *MolId; // MolId[i] => index of the molecule of atom i
	unsigned int *Bonds; // Bonds[i*2] => index-1 of the atom 1 involved in bond i, Bond[i*2+1] => index-1 of atom 2 involved in bond i
	unsigned int *BondType; // BondType[i] type of the bond i
	unsigned int nbBondType;
	bool IsAngle = false; // angle between 3 atoms (e.g. H-O-H angle in ice with TIP4P potential)
	unsigned int nbAngles;
	unsigned int nbAngleType;
	unsigned int *Angles; // Angles[i*3], Angles[i*3+1], Angles[i*3+2], indexes-1 atom 1, 2 and 3 involved in angle i (atom 2 is the center of the angle)
	unsigned int *AngleType;
	double *AtomMass;
	double *AtomCharge;
	// Cell vectors
	double *a1;
	double *a2;
	double *a3;
	double *alength;
	// reciproqual lattice
	double *a1_star;
	double *a2_star;
	double *a3_star;
	double V; // elementary volum
	int *OrthogonalPlanes;
	int *OrthogonalDirs;
	Atom *Motif;
	bool IsOrientedSystem = false;
	bool IsMultisite = false;
	AtomicSystem *OrientedSystem;
	MathTools *MT;
	std::vector<std::vector<unsigned int>> DoNotSep; // contain the number neighbour which should not be separe from a given atom type <=> do not separe atom type DoNotSep[i][0] from its DoNotSep[i][1] first neighbors of atom type DoNotSep[i][2]
	std::vector<std::vector<int>> NotSepList; // contain the id of neighbors which should not be separe from the given atom, NotSepList[i][j*4] = id of neighbor j which should not be separe from atom i, NotSepList[i][j*4+k] = multiplicative coefficient with ak(a1, a2, a3) for this atom to be the neighbot
	bool IsDoNotSep = false;
	double *rot_mat_total; // rotation matrix used to pass from the database oriented crystal to the current orientation
	double *TiltTrans_xyz; // transformation matrix used to construct oriented systems
        unsigned int *Stoichiometry;
	bool IsReferenceBondOriParam = false;
	std::vector<double> *ReferenceBondOriParam; // ReferenceBondOriParam[t][s] = reference bond orientational parameter (used for searching atom site and computing order parameter of multisite and non centrosymmetric crystal) of atom type t and crystallographic site s
	std::vector<std::string> BondOriParamProperties;
	// Parameters to read
	std::string FixedParam_Filename = "Fixed_Parameters.dat";
	int CLsearch = 150;
	double TolOrthoBox = 1.;
	double TolOrthoBoxZ = 1.;
	double MinBoxHeight = 10.;
	double MinBoxAside = 10.;
	double *crystal_def; // contains the deformation applied to the crystal for creating orthogonal box (crystal_def[0,1,2] = compression/dilatation in x,y,z direction, crystal_def[3] = shear invariant)
	double shift_x = 0.;
	double shift_y = 0.;
	double shift_z = 0.;
public:
	// constructors
	Crystal(){};
	Crystal(const std::string& crystalName); // constructor searching in the database the different parameter for the construction
	// getters
	const std::string& getName(){ return this->name; }
	const std::string& getCrystallo(){ return this->crystallo; }
	Atom* getMotif(){ return this->Motif; } 
	double* getA1(){ return this->a1; };
	double* getA2(){ return this->a2; };
	double* getA3(){ return this->a3; };
	const double* getA1_star(){ return this->a1_star; };
	const double* getA2_star(){ return this->a2_star; };
	const double* getA3_star(){ return this->a3_star; };
	const double* getALength(){ return this->alength; };
	const int* getOrthogonalPlanes(){ return this->OrthogonalPlanes; };
	const int* getOrthogonalDirs(){ return this->OrthogonalDirs; };
	const double getVol(){ return this->V; };
	unsigned int getNbAtom(){ return this->nbAtom; }
	unsigned int getNbAtomType(){ return this->nbAtomType; }
	unsigned int getNbBondType(){ return this->nbBondType; }
	unsigned int getNbAngleType(){ return this->nbAngleType; }
	unsigned int getNbBond(){ return this->nbBonds; }
	unsigned int* getBonds(){ return this->Bonds; }
	unsigned int getBondType(unsigned int& i){ return this->BondType[i]; }
	unsigned int getAngleType(unsigned int& i){ return this->AngleType[i]; }
	unsigned int getNbAngle(){ return this->nbAngles; }
	unsigned int* getAngles(){ return this->Angles; }
	unsigned int getAtomSite(const unsigned int i){ return this->AtomSite[i]; }
	const bool getIsMultisite(){ return this->IsMultisite; }
	unsigned int getNbAtomSite(const unsigned int typeuint){ return this->NbAtomSite[typeuint-1]; }
	unsigned int *getNbAtomSite(){ return this->NbAtomSite; }
	const double* getTiltTrans(){ return this->TiltTrans_xyz; }
	const std::string getAtomType(const unsigned int typeuint){ return this->AtomType[typeuint-1]; }
	std::string* getAtomType(){ return this->AtomType; }
	const double getAtomMass(const unsigned int typeuint){ return this->AtomMass[typeuint-1]; }
	double* getAtomMass(){ return this->AtomMass; }
	const double getAtomCharge(const unsigned int typeuint){ return this->AtomCharge[typeuint-1]; }
	double* getAtomCharge(){ return this->AtomCharge; }
	const bool getIsCharge(){ return this->IsCharge; }
	const bool getIsMolId(){ return this->IsMolId; }
	const unsigned int getMolId(unsigned int &i){ return this->MolId[i]; }
	const bool getIsBond(){ return this->IsBond; }
	const bool getIsAngle(){ return this->IsAngle; }
	AtomicSystem* getOrientedSystem(){ return this->OrientedSystem; } 
	const double* getRotMat(){ return this->rot_mat_total; };
	bool getIsDoNotSep(){ return this->IsDoNotSep; }
	unsigned int getNotSepList_size(const unsigned int i){ return this->NotSepList[i].size()/4; }
	int getNotSepList(const unsigned int i, const unsigned int j){ return this->NotSepList[i][j]; }
	std::vector<std::vector<unsigned int>> getDoNotSep(){ return this->DoNotSep; }
	unsigned int *getStoich(){ return this->Stoichiometry; }
	std::string getDatabasePath(std::string crystalName);
	std::vector<std::string> getBondOriParamProperties(){ return BondOriParamProperties; }
	std::vector<double> *getReferenceBondOriParam(){ return ReferenceBondOriParam; }
	bool getIsReferenceBondOriParam(){ return IsReferenceBondOriParam; }
	double *GetCrystalDef(){ return crystal_def; }
	// methods
	void read_params();
	void ReadProperties(std::vector<std::string> Properties);
	void read_database();
	void ComputeCrystalDef();
	void computeReciproqual();
	void ComputeOrthogonalPlanesAndDirections();
	void RotateCrystal(const int& h_p, const int& k_p, const int& l_p);
	void RotateCrystal(const double *RotMat);
	void ConstructNotSepList();
	void ConstructNotSepListFromMolId();
	void ConstructOrthogonalCell();
	void computeStoich();
	void ChangeTypes(unsigned int *CorresArray);
	// destructor
	~Crystal();
	
};

#include "AtomicSystem.h"

#endif
