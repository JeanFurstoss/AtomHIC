//**********************************************************************************
//*   Bicrystal.h                                                                  *
//**********************************************************************************
//* This file contains the declaration of the Bicrystal class which is mainly used *
//* for constructing (given the five macroscopic degrees of freedom) and analyzing *
//* (search GB position, width, excess volume) bycristalline systems.              *
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
//*	- better outputs                                                           *
//*	-                                                                          *
//**********************************************************************************

#ifndef BICRYSTAL_H
#define BICRYSTAL_H

#include "AtomHicExport.h"
#include "Crystal.h"
#include "AtomicSystem.h"
#include "BicrystalloTools.h"
#include "ComputeAuxiliary.h"
#include <string>
#include <iomanip>

class ATOMHIC_EXPORT Bicrystal : public AtomicSystem {
protected:
	BicrystalloTools *BT = nullptr;
	std::string NormalDir;
	double Ldir; // box length in the GB normal direction
	double ExcessVol1 = 0.;
	double ExcessVol2 = 0.;
	std::vector<double> ExcessProp;
	std::vector<std::string> ExcessPropName;
	unsigned int nbGB;
	bool IsVacuum;
	bool IsCentered; // in the case where there is vacuum, is the system centered
	double VacuumLo; // low bound coordinate of the vacuum
	double VacuumHi; // high bound coordinate of the vacuum
	double GBPos1 = 0.; // position of the first GB
	double GBPos2 = 0.; // position of the second GB (if there is no vacuum)
	double GBwidth1 = -1.; // width of the first GB
	double GBwidth2 = -1.; // width of the second GB (if there is no vacuum)
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
	double *RotCartToGrain2 = nullptr;
	bool IsRotMatDefine = false;
	bool AreGrainsDefined = false;
	bool IsCrystal2 = false;
	double xl1, xl2, yl1, yl2, zl1, zl2;
	double *CSL_Basis = nullptr;
	double *DSC_Basis = nullptr;
	double sigma; // double to have an information of how good is the basis (normally unsigned int but if double, the computation of CSL is not very precise)
	unsigned int dupX1, dupX2, dupY1, dupY2; // number of duplicate in X and Y dir for crystal 1 and 2 for constructing the bicrystal
	double Mx1, Mx2, My1, My2; // misfit applied to crystal 1 and 2 in dir X and Y for constructing the bicrystal
	//
    	// for full/bonds/angles atoms
	unsigned int *MolId_G1, *MolId_G2;
        unsigned int nbBonds_G1, nbBonds_G2, nbBondType_G1, nbBondType_G2;
        unsigned int *Bonds_G1, *Bonds_G2;
	unsigned int *BondType_G1, *BondType_G2;
        unsigned int nbAngles_G1, nbAngles_G2, nbAngleType_G1, nbAngleType_G2;
        unsigned int *Angles_G1, *Angles_G2;
	unsigned int *AngleType_G1, *AngleType_G2;
	double fac_vacuum = 1.5;
	//
	// Parameters read from Fixed_Parameter.dat file
	std::string FixedParam_Filename = "FixedParameters.ath";
	unsigned int dupMaxFac = 10;
	double maxVarL = 0.3;
	double lowgrain_shift_x = 0.;
	double lowgrain_shift_y = 0.;
        double lowgrain_shift_z = 0.;
	double upgrain_shift_x = 0.;
	double upgrain_shift_y = 0.;
        double upgrain_shift_z = 0.;
	bool FullGrains = true;
        double MaxMisfit = 0.02;
	unsigned int MaxDup = 100;
	double GBspace = 2.;
public:
	// constructors
	Bicrystal(){};
	Bicrystal(const std::string& filename);
	Bicrystal(const std::string& filename, const std::string CrystalName);
	Bicrystal(const std::string& filename, const std::string NormalDir, const std::string CrystalName);
	Bicrystal(const std::string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p, bool rationalize, std::vector<std::string> Properties=std::vector<std::string>(), int h_px=0, int k_px=0, int l_px=0, bool vacuum=false);// Constructor for bicrystal with plane GB with given misorientation and GB plane
	Bicrystal(const std::string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p, std::vector<int> FacetsType, double L_facet, std::vector<std::string> Properties=std::vector<std::string>(), int h_px=0, int k_px=0, int l_px=0);// Constructor for bicrystal with facetted GB with given misorientation and GB plane and facet type
	Bicrystal(const std::string& crystalName, int h_a, int k_a, int l_a, double theta); // constructor using only misorientation freedom degree (bicrystallo analyzis)
	// getters
	ComputeAuxiliary *get_CA(){ return this->CA; }
	unsigned int getNbGB(){ return nbGB; }
	std::vector<double> getGBExcessProp(){ return ExcessProp; }
	std::vector<std::string> getGBExcessPropName(){ return ExcessPropName; }
	double getGBPos1(){ return this->GBPos1; }
	double getGBwidth1(){ return this->GBwidth1; }
	double getExcessVol1() { return this->ExcessVol1; }
	double getGBPos2(){ return this->GBPos2; }
	double getGBwidth2(){ return this->GBwidth2; }
	double getExcessVol2() { return this->ExcessVol2; }
	double getPrefac1(){ return this->prefac_test; }
	Crystal* getCrystal2(){ return this->_MyCrystal2; }
	double getxl1(){ return this->xl1; }
	double getyl1(){ return this->yl1; }
	double getxl2(){ return this->xl2; }
	double getyl2(){ return this->yl2; }
	double getDupX1(){ return this->dupX1; }
	double getDupX2(){ return this->dupX2; }
	double getDupY1(){ return this->dupY1; }
	double getDupY2(){ return this->dupY2; }
	double getMx1(){ return this->Mx1; }
	double getMx2(){ return this->Mx2; }
	double getMy1(){ return this->My1; }
	double getMy2(){ return this->My2; }
	double *getRotMatG2(){ return this->RotCartToGrain2; }
	double getSigma(){ return this->sigma; }
	double* getCSL_Basis(){ return this->CSL_Basis; }
	//
	double* getDSC_Basis() { return this->DSC_Basis; }
	// methods
	void read_params();
	void print_Grains(bool vacuum = false, std::string filenameG1="Grain1.lmp", std::string filenameG2="Grain2.lmp");
	void searchGBPos();
	void ComputeExcessVolume();
	void ComputeExcessProp(std::string PropName);
	void setOrientedCrystals(const std::string& crystalName, bool rationalize, std::vector<std::string> Properties=std::vector<std::string>(), int h_px=0, int k_px=0, int l_px=0); // method initializing 2 crystals with a given misorientation relationship and plane
	//
	void ShiftGrainsAlongDSC(unsigned int n1, unsigned int n2, unsigned int n3, bool vacuum=false);
	void ShiftGrainsAlongCSL(unsigned int n1, unsigned int n2, unsigned int n3, bool vacuum=false);
	void ShiftGrainsAlongUCInPlane(unsigned int n1, unsigned int n2, bool vacuum=false); 
	void PasteGrains(AtomicSystem* grain1, AtomicSystem* grain2);
	void ReadProperties(std::vector<std::string> Properties);
	void ShowPossibleFacets(unsigned int max_hkl_u=2);
	void InitCSL_DSC();
	//changed name SampleGBComplexionDSC
	// destructor
	~Bicrystal();
};

#endif
