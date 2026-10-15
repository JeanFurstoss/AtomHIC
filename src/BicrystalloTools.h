//**********************************************************************************
//*   BicrystalloTools.h                                                           *
//**********************************************************************************
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
//*	- refine algorythm for rationalizing GB                                    *
//*	- implement non iterative method for searching CSL of cubic crystals       *
//*	-                                                                          *
//**********************************************************************************

#ifndef BICRYSTALLOTOOLS_H
#define BICRYSTALLOTOOLS_H

#include "AtomHicExport.h"
#include "MathTools.h"
#include "Crystal.h"
#include "Displays.h"
#include <string>
#include <vector>
#include <iomanip>

class ATOMHIC_EXPORT BicrystalloTools {
protected:
	
	MathTools *MT;
	Displays Dis;
	
	double *current_a1 = nullptr;
	double *current_a2 = nullptr;
	double *current_CSL = nullptr;
	double *current_DSC = nullptr;
	int current_sigma;
	Crystal *MyCSL_Crystal = nullptr;

	// Parameters read from Fixed_Parameter.dat file
	std::string FixedParam_Filename = "FixedParameters.ath";
	double theta_max_rot_ax_rat = 1.7e-2;
	unsigned int MaxHKL_rot_angle_rat = 75;
	double tol_dist_rot_angle_rat = 5e-2;
	double SigmaMax = 500;

public:
	// constructors
	BicrystalloTools();
	// CSL and DSC tools
	// the following methods return sigma value and store the compute CSL, DSC, lattices and sigma value
	int ComputeCSL_DSC(Crystal *Cryst1, Crystal *Cryst2, double *CSL_Basis, double *DSC_Basis, unsigned int verbose);
	int ComputeCSL_DSC(Crystal *Cryst, double *rot_mat, double *CSL_Basis, double *DSC_Basis, unsigned int verbose);
	int ComputeCSL_DSC(Crystal *Cryst, int h_a, int k_a, int l_a, double theta, double *CSL_Basis, double *DSC_Basis, unsigned int verbose);
	int ComputeCSL_DSC(Crystal *Cryst, double *rot_ax_func, double theta_func, double *CSL_Basis, double *DSC_Basis, unsigned int verbose);
	int ComputeCSL_DSC(double *a1, double *a2, double *CSL_Basis, double *DSC_Basis, unsigned int verbose);
	// the following method return succes of the CSL/DSC research
	bool searchCSL(Crystal *Cryst1, Crystal *Cryst2, double *CSL_Basis, unsigned int verbose); // search CSL of two crystals
	bool searchCSL(Crystal *Cryst, double *rot_mat, double *CSL_Basis, unsigned int verbose); // search CSL of a crystal with a given rotation matrix
	bool searchCSL(Crystal *Cryst, int h_a, int k_a, int l_a, double theta, double *CSL_Basis, unsigned int verbose); // search CSL of a crystal given [hkl] rotation axis and angle
	bool searchCSL(Crystal *Cryst, double *rot_ax_func, double theta_func, double *CSL_Basis, unsigned int verbose); // search CSL of a crystal given rotation axis (in cartesian frame) and angle
	bool searchCSL(double *a1, double *a2, double *CSL_Basis, unsigned int verbose); // search CSL primitive lattice based on the work of Xie et al. 2020 and Bonnet and Rolland 1975
	bool searchDSC(Crystal *Cryst1, Crystal *Cryst2, double *DSC_Basis, unsigned int verbose); // search DSC of two crystals
	bool searchDSC(Crystal *Cryst, double *rot_mat, double *DSC_Basis, unsigned int verbose); // search DSC of a crystal with a given rotation matrix
	bool searchDSC(Crystal *Cryst, int h_a, int k_a, int l_a, double theta, double *DSC_Basis, unsigned int verbose); // search DSC of a crystal given [hkl] rotation axis and angle
	bool searchDSC(Crystal *Cryst, double *rot_ax_func, double theta_func, double *DSC_Basis, unsigned int verbose); // search DSC of a crystal given rotation axis (in cartesian frame) and angle
	bool searchDSC(double *a1, double *a2, double *DSC_Basis, unsigned int verbose); // search DSC primitive lattice assuming that DSC is reciprocal of the CSL of the reciproqual lattices (Grimmer's relationship)
											 //
	//double getSigma(){ return this->sigma; }
	//double* getCSL_Basis(){ return this->CSL_Basis; }
	////
	//double* getDSC_Basis() { return this->DSC_Basis; }
	//bool getIsDSC_Basis() { return this->IsDSC_Basis; }
	//void printDSC();
	// methods
	void read_params();
	void ReadProperties(std::vector<std::string> Properties);
	void SaveLattices(double *a1_s, double *a2_s);
	void SaveCSL_DSC(double *csl, double *dsc);
	void PrintCurrentCSL_DSC();
	int EndComputeCSL_DSC(int sigma_csl, int sigma_dsc, double *CSL_Basis, double *DSC_Basis);
	bool Check_CurrentBasis(); // Check if volum of CSL and DSC is consistent with Sigma value and lattice volume
	void generateCSL();
	double RationalizeOri(Crystal *Cryst, int h_a_func, int k_a_func, int l_a_func, double theta_func, int &hp_near, int &kp_near, int &lp_near, int *CSL_vec);
	//double RationalizeOri(int h_a, int k_a, int l_a, double theta, double *rot_ax, int *CSL_vec);// return the rotation angle corresponding to the closest rational GB and a known CSL vector due to this rationalization
	//void printCSL(const std::string filename);
	// destructor
	~BicrystalloTools();
};

#endif
