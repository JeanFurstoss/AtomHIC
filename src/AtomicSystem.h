//**********************************************************************************
//*   AtomicSystem.h                                                               *
//**********************************************************************************
//* This file contains the declaration of the AtomicSystem class which one of the  *
//* basis class of AtomHIC. This class is mainly used to read/print dump files,    *
//* computing neighbours (warning the variable for accessing neighbours should 	   *
//* really be an unsigned long int variable), etc.                                                     *
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
//* 	- change the Atom struct because we don't need it and it imply to do some  *
//* copies in other classes							   *
//*	- think about one common neighbour research in all code (for the moment    *
//* there is an other one in the Descriptor class)                                 *
//*	- add non periodic boundary conditions (could be link to the above because *
//* the neighbour search in Descriptors does not have PBC)                         *
//*	- put security factor for neighbour research in FixedParameters            *
//*	-                                                                          *
//**********************************************************************************

#ifndef ATOMICSYSTEM_H
#define ATOMICSYSTEM_H

#include "AtomHicExport.h"
#include "MathTools.h"
#include "MyStructs.h"
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "Crystal.h"

class Displays;

class ATOMHIC_EXPORT AtomicSystem {
protected:
	unsigned int nbAtom;
	Displays *Dis;
	unsigned int nbAtomType;
	unsigned int MaxAtomType=100;
	std::string *AtomType; // the different atom species present in the system
	double *AtomMass;
	double *AtomCharge;
	Atom *AtomList = nullptr; // List of the atom belonging to the system
	bool IsAtomListMine = true;
	Position *WrappedPos;
	bool IsWrappedPos = false;
	unsigned long int nbMaxN; // Maximum number of neighbours (computed as function of cuttof radius)
	unsigned int *Neighbours; // List of neigbors, structured as 2D array [ nbAtom x (nbMaxNeighbours+1) ] where the first value of each line corresponds to the number of neighbours for the ith atom, i.e. Neighbours[i*(nbMaxN+1)] = nb of neighbour of atom i, Neighbours[i*(nbMaxN+1)+j+1] = id of the jth neighbour of atom i
	int *CLNeighbours; // List of coefficient (NclX,Ncly,Nclz) applied to the atom position to be a neighbour [ nbAtom x (nbMaxNeighbours*3) ] where the first value of each line corresponds to the number of neighbours for the ith atom, i.e. CLNeighbours[i*nbMaxN*3+j*3] = integer applied to H1 to the jth neighbour of atom i for being its neighbour, CLNeighbours[i*nbMaxN*3+j*3+1] = integer applied to H2 to the jth neighbour of atom i for being its neighbour, CLNeighbours[i*nbMaxN*3+j*3+2] = integer applied to H3 to the jth neighbour of atom i for being its neighbour
	bool IsNeighbours = false;
	std::vector<int> *NotSepTag; // array containing the atom which should not be separated: NotSepTag[i][0] = -1 => the atom i is stored because of the NotSep list, 0 => this atom is not related to the NotSep list, n => the atom has stored n neighbour because of the NotSep list, in the latter case: NotSepTag[i][j+1] return the index of the jth atom which has been stored with the ith one due to the NotSepList
	bool IsNotSepTag = false;
	double *H1 = nullptr; // cell vectors of the system
	double *H2 = nullptr;
	double *H3 = nullptr;
	bool IsCellVecMine = true;
	double *G1 = nullptr; // inverse cell vectors of the system
	double *G2 = nullptr;
	double *G3 = nullptr;
	bool IsG = false;
	bool IsElem = false;
	bool IsVel = false;
	double timestep = 0.; // timestep related to the atomic system (cfg lammps file)
	bool IsCharge = false;
	bool IsMolId = false; // is molecule id (atom_style full in LAMMPS)
	bool IsMolIdMine = true; // is molecule id (atom_style full in LAMMPS)
	bool IsBond = false; // bond between 2 atoms (e.g. O-H bonds in ice with TIP4P potential)
	bool IsBondMine = true; // bond between 2 atoms (e.g. O-H bonds in ice with TIP4P potential)
	unsigned int nbBonds = 0;
	unsigned int *MolId; // MolId[i] => index of the molecule of atom i
	unsigned int *Bonds; // Bonds[i*2] => index-1 of the atom 1 involved in bond i, Bond[i*2+1] => index-1 of atom 2 involved in bond i
	unsigned int *BondType; // BondType[i] type of the bond i
	unsigned int nbBondType = 0;
	bool IsAngle = false; // angle between 3 atoms (e.g. H-O-H angle in ice with TIP4P potential)
	bool IsAngleMine = true; // angle between 3 atoms (e.g. H-O-H angle in ice with TIP4P potential)
	unsigned int nbAngles = 0;
	unsigned int nbAngleType = 0;
	unsigned int *Angles; // Angles[i*3], Angles[i*3+1], Angles[i*3+2], indexes-1 atom 1, 2 and 3 involved in angle i (atom 2 is the center of the angle)
	unsigned int *AngleType;
	bool IsPeriodicArr = false; // .data file of LAMMPS write_data containing info on how apply BC for bonds/angles to be effective
	int *PeriodicArr; // PeriodicArr[i*3+d] = times to apply the PC in d cell vector for atom i to be at the write place for bonds and angles to be preserved
	bool IsTilted = false;
	bool IsCrystalDefined = false;
	bool IsCrystalMine = false;
	bool FilenameConstructed = false;
	bool AtomListConstructed = false;
	Crystal *_MyCrystal;	
	std::vector<double*> Aux; // auxiliary atom properties => Aux[a][i*Aux_size[a]+d] = dimension d (over Aux_size[a]) of the auxiliary property with name Aux_name[a] of atom i
	std::vector<unsigned int> Aux_size; // size of auxiliary atom properties
	bool IsSetAux = false;
	std::vector<std::string> Aux_name; // auxiliary atom properties
	MathTools *MT;
	std::vector<double*> density_prof; // density_prof[i][j*2] = density of auxiliary property i at sampled point j, density_prof[i][j*2+1] = coordinate (in density_name[i*2+1] direction) of the sample point j
	std::vector<int> density_nbPts; // density_nbPts[i] = number of sampled point
	std::vector<std::string*> density_name; // auxiliary atom properties => density_name[i*2] = auxiliary property used to compute density, density_name[i*2+1] = direction along which the density has been computed
	std::string File_Heading; // head of lmp printed file
	// Parameters read from Fixed_Parameter.dat file
	std::string FixedParam_Filename = "Fixed_Parameters.dat";
	double r_cut_n;
	double current_rc_neigh;
	int l_sph_st;
public:
	AtomicSystem();
	AtomicSystem(Crystal *_MyCrystal, double xhi, double yhi, double zhi, std::vector<int> cl_box); // construct atomic system from crystal and cell size
	AtomicSystem(const std::string& filename); // construct AtomicSystem by reading file
	AtomicSystem(AtomicSystem *AtSys, unsigned int &nbSys, std::string &dir); // construct AtomicSystem with multiple AtomicSystems by merging along a given direction (to be used for replacing bicrystal construction in addition with the new duplicate method)
	AtomicSystem(Atom *AtomList, unsigned int nbAtom, Crystal *_MyCrystal, double *H1, double *H2, double *H3); // construct AtomicSystem giving AtomList and cell vectors 
	AtomicSystem(Atom *AtomList, unsigned int nbAtom, Crystal *_MyCrystal, double *H1, double *H2, double *H3, unsigned int *MolId); // construct AtomicSystem giving AtomList and cell vectors 
	AtomicSystem(Atom *AtomList, unsigned int nbAtom, Crystal *_MyCrystal, double *H1, double *H2, double *H3, unsigned int *MolId, unsigned int nbBonds, unsigned int nbBondType, unsigned int *Bonds, unsigned int *BondType); // construct AtomicSystem giving AtomList and cell vectors 
	AtomicSystem(Atom *AtomList, unsigned int nbAtom, Crystal *_MyCrystal, double *H1, double *H2, double *H3, unsigned int *MolId, unsigned int nbBonds, unsigned int nbBondType, unsigned int *Bonds, unsigned int *BondType, unsigned int nbAngles, unsigned int nbAngleType, unsigned int *Angles, unsigned int *AngleType); // construct AtomicSystem giving AtomList and cell vectors 
	bool FilenameConstructor(const std::string& filename);
	void AtomListConstructor(Atom *AtomList, unsigned int nbAtom, Crystal *_MyCrystal, double *H1, double *H2, double *H3); // construct AtomicSystem giving AtomList and cell vectors
	// getters
	std::string getAtomType(const unsigned int i){ return this->AtomType[i]; };
	double getAtomMass(const unsigned int i){ return this->AtomMass[i]; };
	double getAtomCharge(const unsigned int i){ return this->AtomCharge[i]; };
	unsigned int getNbAtomType(){ return this->nbAtomType; };
	unsigned int getNbAtom(){ return this->nbAtom; }
	unsigned int getNbBonds(){ return this->nbBonds; }
	unsigned int getNbBondType(){ return this->nbBondType; }
	unsigned int getNbAngles(){ return this->nbAngles; }
	unsigned int getNbAngleType(){ return this->nbAngleType; }
	unsigned int *getBonds(){ return this->Bonds; }
	unsigned int *getMolId(){ return this->MolId; }
	unsigned int getMolId(unsigned int &i){ return this->MolId[i]; }
	unsigned int getBondType(unsigned int &i){ return this->BondType[i]; }
	unsigned int *getAngles(){ return this->Angles; }
	unsigned int getAngleType(unsigned int &i){ return this->AngleType[i]; }
	Atom getAtom(const unsigned int Id){ return this->AtomList[Id]; }
	double* getAux(const unsigned int AuxId){ return this->Aux[AuxId]; }
	unsigned int getNbAux(){ return this->Aux.size(); }
	unsigned int getAux_size(unsigned int &i){ return this->Aux_size[i]; }
	std::string getAux_name(unsigned int &i){ return this->Aux_name[i]; }
	double* getDensityProf(const unsigned int DensityId){ return this->density_prof[DensityId]; }
	unsigned int getAuxIdAndSize(std::string auxname, unsigned int &size); 
	double get_current_rc(){ return this->current_rc_neigh; }
	Position getWrappedPos(const unsigned int AuxId){ if( !IsWrappedPos ) computeWrap(); return this->WrappedPos[AuxId]; }
	unsigned int* getNeighbours(){ return this->Neighbours; }
	unsigned int getNeighbours(const unsigned long int Id){ return this->Neighbours[Id]; } // Warning Id should really be a long variable !!!
	int* getCLNeighbours(){ return this->CLNeighbours; }
	std::vector<int>* getNotSepTag(){ return this->NotSepTag; }
	int getCLNeighbours(const unsigned long int Id){ return this->CLNeighbours[Id]; } // Warning Id should really be a long variable !!!
	unsigned int getNbMaxN(){ return this->nbMaxN; }
	bool getIsNeighbours(){ return this->IsNeighbours; }
	bool getIsCrystalDefined(){ return this->IsCrystalDefined; }
	bool getIsElem(){ return this->IsElem; }
	bool getIsCharge(){ return this->IsCharge; }
	bool getIsMolId(){ return this->IsMolId; }
	bool getIsBond(){ return this->IsBond; }
	bool getIsAngle(){ return this->IsAngle; }
	bool getIsVel(){ return this->IsVel; }
	bool getIsSetAux(){ return this->IsSetAux; }
	Crystal* getCrystal(){ return this->_MyCrystal; }
	double* getH1(){ return this->H1; }
	double* getH2(){ return this->H2; }
	double* getH3(){ return this->H3; }
	double* getG1(){ return this->G1; }
	double* getG2(){ return this->G2; }
	double* getG3(){ return this->G3; }
	double get_rcut(){ return this->r_cut_n; }
	int get_lsph(){ return this->l_sph_st; }
	// setters
	void setAux(const double* aux, const std::string& AuxName);
	void setAux_vec(const double* aux, const unsigned int size, const std::string& AuxName);
	void modifyAux_vec(const double* aux, const std::string& AuxName);
	void setAux(const unsigned int* aux, const std::string& AuxName);
	void setAux(const int* aux, const std::string& AuxName);
	void setCrystal(Crystal* MyCrystal);
	void setCrystal(const std::string& CrystalName);
	void set_File_Heading(const std::string& Heading){ this->File_Heading = Heading; }
	// methods
	void MakeSurfaceNeutral(); // try to have neutral surfaces (considering only z-oriented surface)
	double ComputeAverageDistance();
	void ComputeNotSepList();
	void UpdateTypes2Crystal();
	void read_params_atsys();
	void computeInverseCellVec();
	bool ReadAtomicFile(const std::string& filename);
	bool read_lmp_file(const std::string& filename);
	bool read_cfg_file(const std::string& filename);
	bool read_other_cfg(const std::string& filename);
	void printSystem(const std::string& filename); // print a file containing atom type, charge, masse and position
	void print_lmp(const std::string& filename);
	void print_cfg(const std::string& filename);
	void printSystem_aux(const std::string& filename, const std::string& AuxId);
	std::string SystemCharacteristics(bool cfg=false);
	unsigned int searchNeighbours(const double& rc); // return nbNMax which is crucial for findings neighbours from the list
	unsigned int searchNeighbours_restricted(const double& rc, const std::vector<unsigned int> & IndexToSearch, const std::vector<unsigned int> & IndexForSearch); // return nbNMax which is crucial for findings neighbours from the list
	void computeWrap();
	unsigned int Compute1dDensity(std::string auxname, std::string dir, double sigma, unsigned int nbPts); // compute and store the 1D density profile of a given auxiliary property, the metho return the index of the given density
	void Print1dDensity(std::string filename, std::string auxname);
	//applyshift
	void ApplyShift(const double &shift_x, const double &shift_y, const double &shift_z);
	void duplicate(const unsigned int &nx, const unsigned int &ny, const unsigned int &nz);
	Position getWrappedPosition(unsigned int &i) const { return this->WrappedPos[i]; }
	void ComputePeriodicArr();
	Atom& getAtomRef(unsigned int &i) { return this->AtomList[i]; }
	//
	//void ApplyShift(const double &shift_x, const double &shift_y, const double &shift_z);
	//Position getWrappedPosition(unsigned int i) const { return this->WrappedPos[i];}
	//Atom& getAtomRef(unsigned int i) { return this->AtomList[i];}
	
	//
	void deleteNeighList(){
		if( this->IsNeighbours ){
			delete[] this->Neighbours;
			delete[] this->CLNeighbours;
		}
	}
	std::vector<unsigned int> selectAtomInBox(const double x_lo,const double x_hi,const double y_lo,const double y_hi,const double z_lo,const double z_hi);
	~AtomicSystem();
};

#include "Displays.h"

#endif
