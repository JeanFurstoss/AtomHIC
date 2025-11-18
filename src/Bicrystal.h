//**********************************************************************************
//*   Bicrystal.h                                                                  *
//**********************************************************************************
//* This file contains the declaration of the Bicrystal class which is mainly used *
//* for constructing (given the five macroscopic degrees of freedom) and analyzing *
//* (search GB position, width, excess volume) bycristalline systems.              *
//* This class has numerous bicrystallophic tools such as an iterative algorythm   *
//* allowing to find (near) CSL and DSC lattices even for non cubic crystals, this *
//* part is based on the work of Bonnet and Rolland 1975 and allow also for        *
//* rationalize a general GB							   *
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
//*	- more generic definition of the number of mesh points and gaussian width  *
//* used for searching GB position and width                                       *
//*	- refine algorythm for rationalizing GB                                    *
//*	- implement non iterative method for searching CSL of cubic crystals       *
//*	- better outputs                                                           *
//*	-                                                                          *
//**********************************************************************************

#ifndef BICRYSTAL_H
#define BICRYSTAL_H

#include "AtomHicExport.h"
#include "Crystal.h"
#include "AtomicSystem.h"
#include "ComputeAuxiliary.h"
#include <string>
#include <iomanip>

class ATOMHIC_EXPORT Bicrystal : public AtomicSystem {
protected:
	std::string NormalDir;
	double Ldir; // box length in the GB normal direction
	double ExcessVol;
	bool IsVacuum;
	bool IsCentered; // in the case where there is vacuum, is the system centered
	double VacuumLo; // low bound coordinate of the vacuum
	double VacuumHi; // high bound coordinate of the vacuum
	double GBPos1; // position of the first GB
	double GBPos2; // position of the second GB (if there is no vacuum)
	double GBwidth1; // width of the first GB
	double GBwidth2; // width of the second GB (if there is no vacuum)
	double MaxPos; // maximum pos of atoms
	double MinPos; // min --
	int h_a, k_a, l_a; // miller indices of rotation axis
	int h_p, k_p, l_p; // miller indices of GB plane
        double theta; // misorientation angle
	double SystemLength;
	ComputeAuxiliary *CA;
	bool IsCA = false;
	double prefac_test;
	bool IsMassDensity;
	Crystal *_MyCrystal2;
	AtomicSystem *Grain1;
	AtomicSystem *Grain2;
	Atom *AtomList_G1;
	double *H1_G1, *H2_G1, *H3_G1;
	Atom *AtomList_G2;
	double *H1_G2, *H2_G2, *H3_G2;
	double *RotGrain1ToGrain2;
	double *RotCartToGrain2;
	bool IsRotMatDefine = false;
	bool AreGrainsDefined = false;
	bool IsCrystal2 = false;
	double xl1, xl2, yl1, yl2;
	std::vector<double*> CSL;
	std::vector<double*> NodesG1;
	std::vector<double*> NodesG2;
	double *CSL_Basis;
	double sigma; // double to have an information of how good is the basis (normally unsigned int but if double, the computation of CSL is not very precise)
	bool IsCSL = false;
	bool IsCSL_Basis = false;
	unsigned int dupX1, dupX2, dupY1, dupY2; // number of duplicate in X and Y dir for crystal 1 and 2 for constructing the bicrystal
	double Mx1, Mx2, My1, My2; // misfit applied to crystal 1 and 2 in dir X and Y for constructing the bicrystal
	//
	double* DSC_Basis = nullptr; //does not point to anything at creation
    bool IsDSC_Basis = false; //whether DSC_Basis is valid or properly initialized
    	// for full/bonds/angles atoms
	unsigned int *MolId_G1, *MolId_G2;
        unsigned int nbBonds_G1, nbBonds_G2, nbBondType_G1, nbBondType_G2;
        unsigned int *Bonds_G1, *Bonds_G2;
	unsigned int *BondType_G1, *BondType_G2;
        unsigned int nbAngles_G1, nbAngles_G2, nbAngleType_G1, nbAngleType_G2;
        unsigned int *Angles_G1, *Angles_G2;
	unsigned int *AngleType_G1, *AngleType_G2;
	//
	// Parameters read from Fixed_Parameter.dat file
	std::string FixedParam_Filename = "FixedParameters.ath";
	bool FullGrains = true;
	double theta_max_rot_ax_rat = 1.7e-2;
	unsigned int MaxHKL_rot_angle_rat = 75;
	double tol_dist_rot_angle_rat = 5e-2;
	double tolpos_known_CSL = 5e-1;
	double tol_CSL_integer = 1e-2;
	double tolAlignment_CSL = 0.1;
        double MaxMisfit = 0.02;
	unsigned int MaxDup = 100;
	double GBspace = 2.;

public:
	// constructors
	Bicrystal(){};
	Bicrystal(const std::string& filename);
	Bicrystal(const std::string& filename, const std::string CrystalName);
	Bicrystal(const std::string& filename, const std::string NormalDir, const std::string CrystalName);
	Bicrystal(const std::string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p, bool rationalize, std::vector<std::string> Properties=std::vector<std::string>());// Constructor for bicrystal with plane GB with given misorientation and GB plane
	Bicrystal(const std::string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p, std::vector<int> FacetsType, unsigned int N_facet);// Constructor for bicrystal with facetted GB with given misorientation and GB plane and facet type
	Bicrystal(const std::string& crystalName, int h_a, int k_a, int l_a, double theta); // constructor using only misorientation freedom degree (bicrystallo analyzis)
	// getters
	ComputeAuxiliary *get_CA(){ return this->CA; }
	double getGBPos1(){ return this->GBPos1; }
	double getGBwidth1(){ return this->GBwidth1; }
	double getPrefac1(){ return this->prefac_test; }
	double getExcessVol() { return this->ExcessVol; }
	Crystal* getCrystal2(){ return this->_MyCrystal2; }
	double getxl1(){ return this->xl1; }
	double getyl1(){ return this->yl1; }
	double getxl2(){ return this->xl2; }
	double getyl2(){ return this->yl2; }
	double getSigma(){ return this->sigma; }
	double* getCSL_Basis(){ return this->CSL_Basis; }
	//
	double* getDSC_Basis() { return this->DSC_Basis; }
    bool getIsDSC_Basis() { return this->IsDSC_Basis; }
    void printDSC();
	// Accesses the DSC_Basis data structure and checks its validity, with a method to display its contents.
	// methods
	void read_params();
	void print_Grains(bool vacuum = false);
	void searchGBPos();
	void ComputeExcessVolume();
        //void searchCSL(int h_a, int k_a, int l_a, double theta, int *CSL_vec, unsigned int verbose=0);
        bool searchCSL(double *rot_ax, double theta, int *CSL_vec, unsigned int verbose=0);
	void generateCSL();
	void solve_DSC(const int *u, const unsigned int L, const double *B, double *DSC_Base, double tol);
	void setOrientedCrystals(const std::string& crystalName, bool rationalize, std::vector<std::string> Properties=std::vector<std::string>()); // method initializing 2 crystals with a given misorientation relationship and plane
	double RationalizeOri(int h_a, int k_a, int l_a, double theta, double *rot_ax, int *CSL_vec);// return the rotation angle corresponding to the closest rational GB and a known CSL vector due to this rationalization
	void searchGBSize(const int h_p_func, const int k_p_func, const int l_p_func);
	void printCSL(const std::string filename);
	//
	void ShiftGrainsAlongDSC(unsigned int n1, unsigned int n2, unsigned int n3);
	void ShiftGrainsAlongCSL(unsigned int n1, unsigned int n2, unsigned int n3);
	void ShiftGrainsAlongUCInPlane(unsigned int n1, unsigned int n2, bool vacuum=false); 
	void PasteGrains(AtomicSystem* grain1, AtomicSystem* grain2);
	void ReadProperties(std::vector<std::string> Properties);
	//changed name SampleGBComplexionDSC
	// destructor
	~Bicrystal();
};

#endif
