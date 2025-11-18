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
#include <vector>
#include <string>
#include "Crystal.h"

class ATOMHIC_EXPORT Displays {
private:
	std::chrono::high_resolution_clock::time_point time_beg;
	// all fixed parameters
	std::vector<std::string> FixedParameters;
		// for crystal 
	double TolOrthoBox = 0.25;
	double TolOrthoBoxZ = 0.25;
	double MinBoxHeight = 20.;
	double MinBoxAside = 3.;
	unsigned int CLsearch = 150;
	double shift_x = 0.;
	double shift_y = 0.;
	double shift_z = 0.;
		// bicrystal
        double MaxMisfit = 0.02;
	unsigned int MaxDup = 100;
	double GBspace = 2.;
	unsigned int FullGrains = 1;
		// for CSL/DSC
	double theta_max_rot_ax_rat = 1.7e-2;
	unsigned int MaxHKL_rot_angle_rat = 75;
	double tol_dist_rot_angle_rat = 5e-2;
	double tolpos_known_CSL = 5e-1;
	double tol_CSL_integer = 1e-2;
	double tolAlignment_CSL = 0.1;
       		// for ML and descriptors
	double RatioTestTrain = 0.1;
	std::string FilteringType = "none";
		// Steinhardt descriptors
	std::string mode = "Full";
	double rc = 5.;
	int l_sph = 10;
	std::string SteinhardtStyle = "Mono";
	std::string AverageStyle = "none";
		// for GMM
	double tol_Lkh_EM = 1e-4;
	unsigned int MaxIter_EM = 100;
	unsigned int nbInit_GMM = 10;
	std::string InitMethod_GMM = "KMEANSPP";
	unsigned int nb_bic_increase_GMM = 1;
	double fac_elbow_GMM = 0.1;
	std::string after_elbow_choice_GMM = "Max";
		// for KMeans
	double tol_KM = 1e-5;
	unsigned int MaxIter_KM = 500;
	unsigned int nbInit_KM = 10;
	unsigned int nb_sil_increase_KM = 1;
	double fac_elbow_KM = 0.1;
	std::string after_elbow_choice_KM = "Max";
		// For DBScan
	std::string eps_meth = "AUTO";
	std::string minPts_meth = "AUTO";
	unsigned int nbClustMax_DBScan = 200;
	
public:
	// constructors
	Displays();
	// methods
	void Logo();
	void ExecutionTime();
	void DisplayArray(const std::vector<std::vector<std::string>>& elements, const std::vector<std::vector<unsigned int>>& fusion);
	void DisplayArray(const std::vector<std::vector<std::string>>& elements, const std::vector<std::vector<unsigned int>>& fusion, std::ofstream& filetoprint);
	std::string center(const std::string& text, int width);
	void DisplayGB(Crystal *Crystal1, Crystal *Crystal2);
	void DisplayOrthogonalCell(Crystal *Crystal);
	void ProgressBar(unsigned int &Final, unsigned int &current);
	// printers of the different fixed parameters of executables
	void PrintHeadExecMsg();
	void PrintFootExecMsg();
	void ReadFixedParams(std::string filename); 
	void ReadCurrentFixedParams(); // read the FixedParameters.ath file if it exists in the working dir and store the values
	void ReadDatabaseFixedParams(); // read the FixedParameters.ath file in /data/FixedParameters/ 
	void PrintFixedParameters();
	void Printer_AnalyzeCSL();
	void Printer_NoFixedParams();
	void Printer_BondOriParam();
	void Printer_ComputeDescriptors_Steinhardt();
	void Printer_CreateGB();
	void Printer_CreateOrientedPlane();
	void Printer_DBScanClustering();
	void Printer_FitAndSaveGMM();
	void Printer_FitAndSaveKMeans();
	void Printer_SampleGB_GammaSurface();
	// destructor
	~Displays(){};
	
};

#endif
