//**********************************************************************************
//*   Displays.cpp                                                                   *
//**********************************************************************************
//* This file contains the implementation of the Displays class                      *
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


#include "Displays.h"
#include "AtomHicConfig.h"
#include <iostream>
#include <stdlib.h>
#include <omp.h>
#include <iomanip>

using namespace std; 
using namespace std::chrono;

Displays::Displays(){
	time_beg = high_resolution_clock::now();	
}

void Displays::Logo(){

    cout << R"(
       █████╗ ████████╗ ██████╗ ███╗   ███╗██╗  ██╗██╗ ██████╗
      ██╔══██╗╚══██╔══╝██╔═══██╗████╗ ████║██║  ██║██║██╔════╝
      ███████║   ██║   ██║   ██║██╔████╔██║███████║██║██║
      ██╔══██║   ██║   ██║   ██║██║╚██╔╝██║██╔══██║██║██║
      ██║  ██║   ██║   ╚██████╔╝██║ ╚═╝ ██║██║  ██║██║╚██████╗
      ╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝ ╚═════╝
    )";
    cout << "                  (C) 2025 Jean Furstoss " << endl;

    cout << "              =====================================" << endl;
    cout << endl;
    cout << "Calculation running using " << omp_get_max_threads() << " threads" << endl;
    cout << endl;
}

void Displays::ExecutionTime(){
	std::chrono::high_resolution_clock::time_point end = high_resolution_clock::now();
	auto duration = duration_cast<milliseconds>(end - time_beg);
	cout << "Execution time : " << duration.count()/1.e3 << " s" << endl;
}

string Displays::center(const string& text, int width) {
	int padding = width - text.length();
	int leftPadding = padding / 2;
	int rightPadding = padding - leftPadding;
	return string(leftPadding, ' ') + text + string(rightPadding, ' ');
}

void Displays::DisplayArray(const vector<vector<string>>& elements, const vector<vector<unsigned int>>& fusion) {
	unsigned int nbRow = fusion.size();
	unsigned int nbCols = 0;
	for(unsigned int i=0;i<fusion[0].size();i++) nbCols += fusion[0][i];
	for(unsigned int i=1;i<nbRow;i++){
	        unsigned int buf = 0;
		for(unsigned int ii=0;ii<fusion[i].size();ii++) buf += fusion[i][ii];
	    	if( buf != nbCols ){
	    		cerr << "The number of columns is not consistent for printing array, aborting" << endl;
	    		exit(EXIT_FAILURE);
	    	}
	}
	
	// compute column width
	unsigned int width = 0;
	for(unsigned int i=0;i<elements.size();i++){
		for(unsigned int j=0;j<elements[i].size();j++){
	    		if( elements[i][j].length() > width ) width = elements[i][j].length();
	    }
	}
	width += 4;
	
	string line = "+";
	for(unsigned int i=0;i<nbCols;i++) line += string(width, '-') + "+";
	
	cout << line << endl;
	for (unsigned int i = 0; i < nbRow; i++) {
		cout << "|";
		for(unsigned int j=0;j<fusion[i].size();j++){
			unsigned int colspan = fusion[i][j];
			unsigned int mergedWidth = width*colspan+(colspan-1);
			cout << center(elements[i][j], mergedWidth) << "|";
		}
		cout << endl << line << endl;
	}
}

void Displays::DisplayArray(const vector<vector<string>>& elements, const vector<vector<unsigned int>>& fusion, ofstream& filetoprint) {
	unsigned int nbRow = fusion.size();
	unsigned int nbCols = 0;
	for(unsigned int i=0;i<fusion[0].size();i++) nbCols += fusion[0][i];
	for(unsigned int i=1;i<nbRow;i++){
	        unsigned int buf = 0;
		for(unsigned int ii=0;ii<fusion[i].size();ii++) buf += fusion[i][ii];
	    	if( buf != nbCols ){
	    		cerr << "The number of columns is not consistent for printing array, aborting" << endl;
	    		exit(EXIT_FAILURE);
	    	}
	}
	
	// compute column width
	unsigned int width = 0;
	for(unsigned int i=0;i<elements.size();i++){
		for(unsigned int j=0;j<elements[i].size();j++){
	    		if( elements[i][j].length() > width ) width = elements[i][j].length();
	    }
	}
	width += 4;
	
	string line = "+";
	for(unsigned int i=0;i<nbCols;i++) line += string(width, '-') + "+";
	
	filetoprint << line << endl;
	for (unsigned int i = 0; i < nbRow; i++) {
		filetoprint << "|";
		for(unsigned int j=0;j<fusion[i].size();j++){
			unsigned int colspan = fusion[i][j];
			unsigned int mergedWidth = width*colspan+(colspan-1);
			filetoprint << center(elements[i][j], mergedWidth) << "|";
		}
		filetoprint << endl << line << endl;
	}
}

void Displays::DisplayOrthogonalCell(Crystal *Crystal){
	string *Dirs1 = new string[3];
	string *Planes1 = new string[3];
	unsigned int width = 0;
	string *i_1_d = new string[3];
	string *i_1_p = new string[3];
	for(unsigned int i=0;i<3;i++){
		i_1_d[i] = "";
		i_1_p[i] = "";
	}
	if( Crystal->getCrystallo() == "Hexagonal" ){
		for(unsigned int i=0;i<3;i++){
			i_1_d[i] = "_"+to_string(-Crystal->getOrthogonalDirs()[i*3]-Crystal->getOrthogonalDirs()[i*3+1]);
			i_1_p[i] = "_"+to_string(-Crystal->getOrthogonalPlanes()[i*3]-Crystal->getOrthogonalPlanes()[i*3+1]);
		}
	}

	for(unsigned int i=0;i<3;i++){
		Dirs1[i] = "["+to_string(Crystal->getOrthogonalDirs()[i*3]);
		Planes1[i] = "("+to_string(Crystal->getOrthogonalPlanes()[i*3]);
		
		Dirs1[i] += "_"+to_string(Crystal->getOrthogonalDirs()[i*3+1]);
		Planes1[i] += "_"+to_string(Crystal->getOrthogonalPlanes()[i*3+1]);
		
		Dirs1[i] += i_1_d[i];
		Planes1[i] += i_1_p[i];
		
		Dirs1[i] += "_"+to_string(Crystal->getOrthogonalDirs()[i*3+2]);
		Planes1[i] += "_"+to_string(Crystal->getOrthogonalPlanes()[i*3+2]);
		
		Dirs1[i] += "]";
		Planes1[i] += ")";
		if( Dirs1[i].length() > width ) width = Dirs1[i].length();
		if( Planes1[i].length() > width ) width = Planes1[i].length();
	}

	width += 12;
	string line = string(width+2, '-');
	string side = "|";
	side += string(width, ' ')+'|';
	cout << "The system is constructed following :" << endl << endl;
	// Display grain 1
	cout << "        " << "   /" << string(width, ' ') << "/" << endl;
	cout << "        " << "  /" << string(width, ' ') << "/" << endl;
	cout << "        " << " /" << center(Planes1[2], width) << "/ " << Planes1[0] << endl;
	cout << "        " << line << endl;
	cout << "        " << "|" << center(Planes1[1], width) << "|" << endl;
	cout << "        " << side << endl;
	cout << "        " << side << endl;
	cout << "        " << "|" << center("", width) << "|        ^ z   " << endl;
	cout << "        " << side << "        |     " << endl;
	cout << "        " << side << "        ---> x" << endl;
	cout << "        " << "|     ^ " << Dirs1[2] << string(width-7-Dirs1[2].length(), ' ') << "|     y /   " << endl;
	cout << "        " << "|     |" << string(width-6,' ') << "|      v      " << endl;
	cout << "        " << "|     ---> " << Dirs1[0] << string(width-10-Dirs1[0].length(), ' ') << "|" << endl;
	cout << "        " << "|    /" << string(width-5,' ') << "|" << endl;
	cout << "        " << "|   v " << Dirs1[1] << string(width-5-Dirs1[1].length(), ' ') << "|  /" << endl;
	cout << "        " << side << " /" << endl;
	cout << "        " << side << "/" << endl;
	cout << "        " << line << endl;


	delete[] Dirs1;
	delete[] Planes1;
	delete[] i_1_d;
	delete[] i_1_p;

}

void Displays::DisplayGB(Crystal *Crystal1, Crystal *Crystal2){
	string *Dirs1 = new string[3];
	string *Dirs2 = new string[3];
	string *Planes1 = new string[3];
	string *Planes2 = new string[3];
	unsigned int width = 0;
	string *i_1_d = new string[3];
	string *i_2_d = new string[3];
	string *i_1_p = new string[3];
	string *i_2_p = new string[3];
	for(unsigned int i=0;i<3;i++){
		i_1_d[i] = "";
		i_2_d[i] = "";
		i_1_p[i] = "";
		i_2_p[i] = "";
	}
	if( Crystal1->getCrystallo() == "Hexagonal" ){
		for(unsigned int i=0;i<3;i++){
			i_1_d[i] = "_"+to_string(-Crystal1->getOrthogonalDirs()[i*3]-Crystal1->getOrthogonalDirs()[i*3+1]);
			i_1_p[i] = "_"+to_string(-Crystal1->getOrthogonalPlanes()[i*3]-Crystal1->getOrthogonalPlanes()[i*3+1]);
		}
	}
	if( Crystal2->getCrystallo() == "Hexagonal" ){
		for(unsigned int i=0;i<3;i++){
			i_2_d[i] = "_"+to_string(-Crystal2->getOrthogonalDirs()[i*3]-Crystal2->getOrthogonalDirs()[i*3+1]);
			i_2_p[i] = "_"+to_string(-Crystal2->getOrthogonalPlanes()[i*3]-Crystal2->getOrthogonalPlanes()[i*3+1]);
		}
	}
	
	for(unsigned int i=0;i<3;i++){
		Dirs1[i] = "["+to_string(Crystal1->getOrthogonalDirs()[i*3]);
		Planes1[i] = "("+to_string(Crystal1->getOrthogonalPlanes()[i*3]);
		Dirs2[i] = "["+to_string(Crystal2->getOrthogonalDirs()[i*3]);
		Planes2[i] = "("+to_string(Crystal2->getOrthogonalPlanes()[i*3]);
		
		Dirs1[i] += "_"+to_string(Crystal1->getOrthogonalDirs()[i*3+1]);
		Planes1[i] += "_"+to_string(Crystal1->getOrthogonalPlanes()[i*3+1]);
		Dirs2[i] += "_"+to_string(Crystal2->getOrthogonalDirs()[i*3+1]);
		Planes2[i] += "_"+to_string(Crystal2->getOrthogonalPlanes()[i*3+1]);
		
		Dirs1[i] += i_1_d[i];
		Planes1[i] += i_1_p[i];
		Dirs2[i] += i_2_d[i];
		Planes2[i] += i_2_p[i];
		
		Dirs1[i] += "_"+to_string(Crystal1->getOrthogonalDirs()[i*3+2]);
		Planes1[i] += "_"+to_string(Crystal1->getOrthogonalPlanes()[i*3+2]);
		Dirs2[i] += "_"+to_string(Crystal2->getOrthogonalDirs()[i*3+2]);
		Planes2[i] += "_"+to_string(Crystal2->getOrthogonalPlanes()[i*3+2]);
		
		Dirs1[i] += "]";
		Planes1[i] += ")";
		Dirs2[i] += "]";
		Planes2[i] += ")";
		if( Dirs1[i].length() > width ) width = Dirs1[i].length();
		if( Dirs2[i].length() > width ) width = Dirs2[i].length();
		if( Planes1[i].length() > width ) width = Planes1[i].length();
		if( Planes2[i].length() > width ) width = Planes2[i].length();
	}

	width += 12;
	string line = string(width+2, '-');
	string side = "|";
	side += string(width, ' ')+'|';
	cout << "The GB is constructed following :" << endl << endl;
	// Display grain 1
	cout << "        " << "   /" << string(width, ' ') << "/" << endl;
	cout << "        " << "  /" << string(width, ' ') << "/" << endl;
	cout << "        " << " /" << center(Planes1[2], width) << "/ " << Planes1[0] << endl;
	cout << "        " << line << endl;
	cout << "        " << "|" << center(Planes1[1], width) << "|" << endl;
	cout << "        " << side << endl;
	cout << "        " << side << endl;
	cout << "        " << "|" << center("Grain 1", width) << "|        ^ z   " << endl;
	cout << "        " << side << "        |     " << endl;
	cout << "        " << side << "        ---> x" << endl;
	cout << "        " << "|     ^ " << Dirs1[2] << string(width-7-Dirs1[2].length(), ' ') << "|     y /   " << endl;
	cout << "        " << "|     |" << string(width-6,' ') << "|      v      " << endl;
	cout << "        " << "|     ---> " << Dirs1[0] << string(width-10-Dirs1[0].length(), ' ') << "|" << endl;
	cout << "        " << "|    /" << string(width-5,' ') << "|" << endl;
	cout << "        " << "|   v " << Dirs1[1] << string(width-5-Dirs1[1].length(), ' ') << "|  /" << endl;
	cout << "        " << side << " /" << endl;
	cout << "        " << side << "/" << endl;
	cout << "        " << line << endl;

	// Display grain 2
	cout << "        " << "   /" << string(width, ' ') << "/" << endl;
	cout << "        " << "  /" << string(width, ' ') << "/" << endl;
	cout << "        " << " /" << center(Planes2[2], width) << "/ " << Planes2[0] << endl;
	cout << "        " << line << endl;
	cout << "        " << "|" << center(Planes2[1], width) << "|" << endl;
	cout << "        " << side << endl;
	cout << "        " << side << endl;
	cout << "        " << "|" << center("Grain 2", width) << "|" << endl;
	cout << "        " << side << endl;
	cout << "        " << side << endl;
	cout << "        " << "|     ^ " << Dirs2[2] << string(width-7-Dirs2[2].length(), ' ') << "|" << endl;
	cout << "        " << "|     |" << string(width-6,' ') << "|" << endl;
	cout << "        " << "|     ---> " << Dirs2[0] << string(width-10-Dirs2[0].length(), ' ') << "|" << endl;
	cout << "        " << "|    /" << string(width-5,' ') << "|" << endl;
	cout << "        " << "|   v " << Dirs2[1] << string(width-5-Dirs2[1].length(), ' ') << "|  /" << endl;
	cout << "        " << side << " /" << endl;
	cout << "        " << side << "/" << endl;
	cout << "        " << line << endl;

	delete[] Dirs1;
	delete[] Dirs2;
	delete[] Planes1;
	delete[] Planes2;
	delete[] i_1_d;
	delete[] i_2_d;
	delete[] i_1_p;
	delete[] i_2_p;
}


void Displays::ProgressBar(unsigned int &Final, unsigned int &current){
	const int bar_length = 30;
	if( current+1 == Final ){
		cout << endl;
		return;
	}
	double prog = double(current+1)/double(Final);
	cout << "\r[" << string(floor(bar_length*prog),'X') << string(ceil(bar_length*(1-prog)),'-') << "] " << setprecision(3) << 100*prog << "%";
}

void Displays::ReadFixedParams(string filename){
	ifstream file(filename, ios::in);
	string buffer_s, line;
	// crystal
	size_t pos_tolOrthoBox, pos_tolOrthoBoxZ, pos_minBoxHeight, pos_minBoxAside, pos_shift, pos_nblc;
	// bicrystals and CSL/DSC
	size_t pos_thetamax, pos_MaxHKL, pos_toldist, pos_tolpos_kC, pos_tolCSLint, pos_tolAlign, pos_MaxMisfit, pos_GBSpace, pos_MaxDup, pos_FullGrains;
	// descriptors and ML
	size_t pos_filter, pos_rattrain;
	// Steinhardt descriptors
	size_t pos_rcut, pos_lsph, pos_StStyle, pos_AveStyle, pos_mode;
	// GMM
	size_t pos_elfac_GMM, pos_nb_bic_inc, pos_after_el_GMM;
	size_t pos_tol_GMM, pos_maxIter_GMM, pos_nbInit_GMM, pos_InitMethod_GMM;
	// KMeans
	size_t pos_elfac_KM, pos_nb_sil_inc, pos_after_el_KM;
	size_t pos_tol_KM, pos_maxIter_KM, pos_nbInit_KM;
	// DBScan
	size_t pos_eps, pos_minPts, pos_nbClustMax;
	if(file){
		while(file){
			getline(file,line);
			// crystal
			pos_tolOrthoBox=line.find("TOL_ORTHO_BOX ");
			if(pos_tolOrthoBox!=string::npos){
				FixedParameters.push_back("\t *** Parameters related to crystal (construction of orthogonal simulation cells): ***");
				FixedParameters.push_back("");
				istringstream text(line);
				text >> buffer_s >> this->TolOrthoBox;
				FixedParameters.push_back(line);
			}
			pos_tolOrthoBoxZ=line.find("TOL_ORTHO_BOX_Z ");
			if(pos_tolOrthoBoxZ!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->TolOrthoBoxZ;
				FixedParameters.push_back(line);
			}
			pos_minBoxHeight=line.find("MIN_BOX_HEIGHT");
			if(pos_minBoxHeight!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->MinBoxHeight;
				FixedParameters.push_back(line);
			}
			pos_minBoxAside=line.find("MIN_BOX_ASIDE");
			if(pos_minBoxAside!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->MinBoxAside;
				FixedParameters.push_back(line);
			}
			pos_shift=line.find("MOTIF_SHIFT");
			if(pos_shift!=string::npos){
				istringstream text(line);
				text >> buffer_s >> shift_x >> shift_y >> shift_z;
				FixedParameters.push_back(line);
			}
			pos_nblc=line.find("NB_MAX_LC");
			if(pos_nblc!=string::npos){
				istringstream text(line);
				text >> buffer_s >> CLsearch;
				FixedParameters.push_back(line);
			}
			// Bicrystal and CSL/DSC
			pos_thetamax=line.find("THETA_MAX_ROT_ANGLE_RAT ");
			if(pos_thetamax!=string::npos){
				FixedParameters.push_back("");
				FixedParameters.push_back("\t *** Parameters related to CSL/DSC calculations and GB construction: ***");
				FixedParameters.push_back("");
				istringstream text(line);
				text >> buffer_s >> this->theta_max_rot_ax_rat;
				FixedParameters.push_back(line);
			}
			pos_MaxHKL=line.find("MAX_HKL_ROT_ANGLE ");
			if(pos_MaxHKL!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->MaxHKL_rot_angle_rat;
				FixedParameters.push_back(line);
			}
			pos_toldist=line.find("TOL_DIST_RAT_ANGLE ");
			if(pos_toldist!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->tol_dist_rot_angle_rat;
				FixedParameters.push_back(line);
			}
			pos_tolpos_kC=line.find("TOL_POS_KNOWN_CSL_VEC ");
			if(pos_tolpos_kC!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->tolpos_known_CSL;
				FixedParameters.push_back(line);
			}
			pos_tolCSLint=line.find("TOL_CSL_INTEGER ");
			if(pos_tolCSLint!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->tol_CSL_integer;
				FixedParameters.push_back(line);
			}
			pos_tolAlign=line.find("TOL_ALIGNMENT_CSL ");
			if(pos_tolAlign!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->tolAlignment_CSL;
				FixedParameters.push_back(line);
			}
			pos_GBSpace=line.find("GB_SPACE ");
			if(pos_GBSpace!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->GBspace;
				FixedParameters.push_back(line);
			}
			pos_MaxMisfit=line.find("MAX_MISFIT ");
			if(pos_MaxMisfit!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->MaxMisfit;
				FixedParameters.push_back(line);
			}
			pos_MaxDup=line.find("MAX_DUP ");
			if(pos_MaxDup!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->MaxDup;
				FixedParameters.push_back(line);
			}
			pos_FullGrains=line.find("FULL_GRAINS ");
			if(pos_FullGrains!=string::npos){
				istringstream text(line);
				text >> buffer_s >> FullGrains;
				FixedParameters.push_back(line);
			}
			// Descriptors and ML
			pos_filter=line.find("DESCRIPTORS_FILTERING_TYPE ");
			if(pos_filter!=string::npos){
				FixedParameters.push_back("");
				FixedParameters.push_back("\t *** Parameters related to descriptors and machine learning (ML) models: ***");
				FixedParameters.push_back("");
				istringstream text(line);
				text >> buffer_s >> FilteringType;
				FixedParameters.push_back(line);
			}
			pos_rattrain=line.find("ML_RATIO_TEST_TRAIN");
			if(pos_rattrain!=string::npos){
				istringstream text(line);
				text >> buffer_s >> RatioTestTrain;
				FixedParameters.push_back(line);
			}
			// Steinhardt descriptors
			pos_mode=line.find("STEINHARDT_DESCRIPTORS_MODE ");
			if(pos_mode!=string::npos){
				FixedParameters.push_back("");
				FixedParameters.push_back("\t *** Parameters related to Steinhardt descriptors: ***");
				FixedParameters.push_back("");
				istringstream text(line);
				text >> buffer_s >> mode;
				FixedParameters.push_back(line);
			}
			pos_rcut=line.find("STEINHARDT_DESCRIPTORS_RCUT ");
			if(pos_rcut!=string::npos){
				istringstream text(line);
				text >> buffer_s >> rc;
				FixedParameters.push_back(line);
			}
			pos_lsph=line.find("STEINHARDT_DESCRIPTORS_L_SPH ");
			if(pos_lsph!=string::npos){
				istringstream text(line);
				text >> buffer_s >> l_sph;
				FixedParameters.push_back(line);
			}
			pos_StStyle=line.find("STEINHARDT_DESCRIPTORS_STYLE ");
			if(pos_StStyle!=string::npos){
				istringstream text(line);
				text >> buffer_s >> SteinhardtStyle;
				FixedParameters.push_back(line);
			}
			pos_AveStyle=line.find("STEINHARDT_DESCRIPTORS_AVERAGE_STYLE ");
			if(pos_AveStyle!=string::npos){
				istringstream text(line);
				text >> buffer_s >> AverageStyle;
				FixedParameters.push_back(line);
			}
			// GMM
			pos_elfac_GMM=line.find("GMM_ELBOW_FACTOR");
			if(pos_elfac_GMM!=string::npos){
				FixedParameters.push_back("");
				FixedParameters.push_back("\t *** Parameters related to the ML Gaussian Mixture Model (GMM): ***");
				FixedParameters.push_back("");
				istringstream text(line);
				text >> buffer_s >> fac_elbow_GMM;
				FixedParameters.push_back(line);
			}
			pos_nb_bic_inc=line.find("GMM_NB_BIC_INCREASE_FOR_MIN");
			if(pos_nb_bic_inc!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nb_bic_increase_GMM;
				FixedParameters.push_back(line);
			}
			pos_after_el_GMM=line.find("GMM_NO_BIC_MIN_NO_ELBOW_CHOICE");
			if(pos_after_el_GMM!=string::npos){
				istringstream text(line);
				text >> buffer_s >> after_elbow_choice_GMM;
				FixedParameters.push_back(line);
			}
			pos_tol_GMM=line.find("GMM_TOL_LKH_EM");
			if(pos_tol_GMM!=string::npos){
				istringstream text(line);
				text >> buffer_s >> tol_Lkh_EM;
				FixedParameters.push_back(line);
			}
			pos_maxIter_GMM=line.find("GMM_MAX_ITER_EM");
			if(pos_maxIter_GMM!=string::npos){
				istringstream text(line);
				text >> buffer_s >> MaxIter_EM;
				FixedParameters.push_back(line);
			}
			pos_nbInit_GMM=line.find("GMM_NB_INIT");
			if(pos_nbInit_GMM!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nbInit_GMM;
				FixedParameters.push_back(line);
			}
			pos_InitMethod_GMM=line.find("GMM_INIT_METHOD");
			if(pos_InitMethod_GMM!=string::npos){
				istringstream text(line);
				text >> buffer_s >> InitMethod_GMM;
				FixedParameters.push_back(line);
			}
			// KMeans
			pos_elfac_KM=line.find("KMEANS_ELBOW_FACTOR");
			if(pos_elfac_KM!=string::npos){
				FixedParameters.push_back("");
				FixedParameters.push_back("\t *** Parameters related to KMeans ML model: ***");
				FixedParameters.push_back("");
				istringstream text(line);
				text >> buffer_s >> fac_elbow_KM;
				FixedParameters.push_back(line);
			}
			pos_nb_sil_inc=line.find("KMEANS_NB_SIL_INCREASE_FOR_MIN");
			if(pos_nb_sil_inc!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nb_sil_increase_KM;
				FixedParameters.push_back(line);
			}
			pos_after_el_KM=line.find("KMEANS_NO_BIC_MIN_NO_ELBOW_CHOICE");
			if(pos_after_el_KM!=string::npos){
				istringstream text(line);
				text >> buffer_s >> after_elbow_choice_KM;
				FixedParameters.push_back(line);
			}
			pos_tol_KM=line.find("KMEANS_TOL");
			if(pos_tol_KM!=string::npos){
				istringstream text(line);
				text >> buffer_s >> tol_KM;
				FixedParameters.push_back(line);
			}
			pos_maxIter_KM=line.find("KMEANS_MAX_ITER");
			if(pos_maxIter_KM!=string::npos){
				istringstream text(line);
				text >> buffer_s >> MaxIter_KM;
				FixedParameters.push_back(line);
			}
			pos_nbInit_KM=line.find("KMEANS_NB_INIT");
			if(pos_nbInit_KM!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nbInit_KM;
				FixedParameters.push_back(line);
			}
			// DBScan
			pos_nbClustMax=line.find("DBSCAN_NB_MAX_CLUSTER");
			if(pos_nbClustMax!=string::npos){
				FixedParameters.push_back("");
				FixedParameters.push_back("\t *** Parameters related to Density Based Scan model: ***");
				FixedParameters.push_back("");
				istringstream text(line);
				text >> buffer_s >> nbClustMax_DBScan;
				FixedParameters.push_back(line);
			}
			pos_eps=line.find("DBSCAN_EPS");
			if(pos_eps!=string::npos){
				istringstream text(line);
				text >> buffer_s >> eps_meth;
				FixedParameters.push_back(line);
			}
			pos_minPts=line.find("DBSCAN_MINPTS");
			if(pos_minPts!=string::npos){
				istringstream text(line);
				text >> buffer_s >> minPts_meth;
				FixedParameters.push_back(line);
			}
		}
	}
}

void Displays::ReadCurrentFixedParams(){
	string filename = "FixedParameters.ath";
	ReadFixedParams(filename);
}

void Displays::ReadDatabaseFixedParams(){
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+"FixedParameters.ath";
	ReadFixedParams(filename);
}

void Displays::PrintFixedParameters(){
	ReadDatabaseFixedParams();
	cout << "List of fixed parameters used by AtomHIC:" << endl << endl;
	for(unsigned int i=0;i<FixedParameters.size();i++) cout << FixedParameters[i] << endl;
	cout << endl;
}

void Displays::PrintHeadExecMsg(){
	cout << "List of fixed parameter used by this executable:" << endl;
}

void Displays::PrintFootExecMsg(){
	cout << "Description of parameters can be displayed using the DisplayFixedParameters executable or within the /data/FixedParameters/FixedParameters.ath file" << endl;
	cout << "The values of these parameters can be set to other values by creating a FixedParameters.ath file in the working directory" << endl << endl;
}

void Displays::Printer_AnalyzeCSL(){
	ReadCurrentFixedParams();
	PrintHeadExecMsg();
	cout << "\t THETA_MAX_ROT_ANGLE_RAT " << theta_max_rot_ax_rat << endl;
	cout << "\t MAX_HKL_ROT_ANGLE " << MaxHKL_rot_angle_rat << endl;
	cout << "\t TOL_DIST_RAT_ANGLE " << tol_dist_rot_angle_rat << endl;
	cout << "\t TOL_POS_KNOWN_CSL_VEC " << tolpos_known_CSL << endl;
	cout << "\t TOL_CSL_INTEGER " << tol_CSL_integer << endl;
	cout << "\t TOL_ALIGNMENT_CSL " << tolAlignment_CSL << endl;
	PrintFootExecMsg();
}

void Displays::Printer_NoFixedParams(){
	cout << "No fixed parameters are used for this executable" << endl << endl;
}
	
void Displays::Printer_BondOriParam(){
	ReadCurrentFixedParams();
	PrintHeadExecMsg();
	cout << "\t STEINHARDT_DESCRIPTORS_L_SPH " << l_sph << endl;
	cout << "\t STEINHARDT_DESCRIPTORS_STYLE " << SteinhardtStyle << endl;
	PrintFootExecMsg();
}

void Displays::Printer_ComputeDescriptors_Steinhardt(){
	ReadCurrentFixedParams();
	PrintHeadExecMsg();
	cout << "\t STEINHARDT_DESCRIPTORS_MODE " << mode << endl;
	cout << "\t STEINHARDT_DESCRIPTORS_RCUT " << rc << endl;
	cout << "\t STEINHARDT_DESCRIPTORS_L_SPH " << l_sph << endl;
	cout << "\t STEINHARDT_DESCRIPTORS_STYLE " << SteinhardtStyle << endl;
	cout << "\t STEINHARDT_DESCRIPTORS_AVERAGE " << AverageStyle << endl;
	PrintFootExecMsg();
}

void Displays::Printer_CreateGB(){
	ReadCurrentFixedParams();
	PrintHeadExecMsg();
	cout << "\t TOL_ORTHO_BOX " << TolOrthoBox << endl;
	cout << "\t TOL_ORTHO_BOX_Z " << TolOrthoBoxZ << endl;
	cout << "\t MIN_BOX_HEIGHT " << MinBoxHeight << endl;
	cout << "\t MIN_BOX_ASIDE " << MinBoxAside << endl;
	cout << "\t NB_MAX_LC " << CLsearch << endl;
	cout << "\t MOTIF_SHIFT " << shift_x << " " << shift_y << " " << shift_z << endl;
	cout << "\t MAX_MISFIT " << MaxMisfit << endl;
	cout << "\t MAX_DUP " << MaxDup << endl;
	cout << "\t GB_SPACE " << GBspace << endl;
	cout << "\t FULL_GRAINS " << FullGrains << endl;
	cout << "\t THETA_MAX_ROT_ANGLE_RAT " << theta_max_rot_ax_rat << endl;
	cout << "\t MAX_HKL_ROT_ANGLE " << MaxHKL_rot_angle_rat << endl;
	cout << "\t TOL_DIST_RAT_ANGLE " << tol_dist_rot_angle_rat << endl;
	cout << "\t TOL_POS_KNOWN_CSL_VEC " << tolpos_known_CSL << endl;
	cout << "\t TOL_CSL_INTEGER " << tol_CSL_integer << endl;
	cout << "\t TOL_ALIGNMENT_CSL " << tolAlignment_CSL << endl;
	PrintFootExecMsg();
}

void Displays::Printer_CreateOrientedPlane(){
	ReadCurrentFixedParams();
	PrintHeadExecMsg();
	cout << "\t TOL_ORTHO_BOX " << TolOrthoBox << endl;
	cout << "\t TOL_ORTHO_BOX_Z " << TolOrthoBoxZ << endl;
	cout << "\t MIN_BOX_HEIGHT " << MinBoxHeight << endl;
	cout << "\t MIN_BOX_ASIDE " << MinBoxAside << endl;
	cout << "\t NB_MAX_LC " << CLsearch << endl;
	cout << "\t MOTIF_SHIFT " << shift_x << " " << shift_y << " " << shift_z << endl;
	PrintFootExecMsg();
}

void Displays::Printer_DBScanClustering(){
	ReadCurrentFixedParams();
	PrintHeadExecMsg();
	cout << "\t DBSCAN_NB_MAX_CLUSTER " << nbClustMax_DBScan << endl;
	cout << "\t DBSCAN_EPS " << eps_meth << endl;
	cout << "\t DBSCAN_MINPTS " << minPts_meth << endl;
	PrintFootExecMsg();
}

void Displays::Printer_FitAndSaveGMM(){
	ReadCurrentFixedParams();
	PrintHeadExecMsg();
	cout << "\t DESCRIPTORS_FILTERING_TYPE " << FilteringType << endl;
	cout << "\t ML_RATIO_TEST_TRAIN " << RatioTestTrain << endl;
	cout << "\t GMM_TOL_LKH_EM " << tol_Lkh_EM << endl;
	cout << "\t GMM_MAX_ITER_EM " << MaxIter_EM << endl;
	cout << "\t GMM_NB_INIT " << nbInit_GMM << endl;
	cout << "\t GMM_INIT_METHOD " << InitMethod_GMM << endl;
	cout << "\t GMM_NB_BIC_INCREASE_FOR_MIN " << nb_bic_increase_GMM << endl;
	cout << "\t GMM_ELBOW_FACTOR " << fac_elbow_GMM << endl;
	cout << "\t GMM_NO_BIC_MIN_NO_ELBOW_CHOICE " << after_elbow_choice_GMM << endl;
	PrintFootExecMsg();
}

void Displays::Printer_FitAndSaveKMeans(){
	ReadCurrentFixedParams();
	PrintHeadExecMsg();
	cout << "\t DESCRIPTORS_FILTERING_TYPE " << FilteringType << endl;
	cout << "\t ML_RATIO_TEST_TRAIN " << RatioTestTrain << endl;
	cout << "\t KMEANS_TOL " << tol_KM << endl;
	cout << "\t KMEANS_MAX_ITER " << MaxIter_KM << endl;
	cout << "\t KMEANS_NB_INIT " << nbInit_KM << endl;
	cout << "\t KMEANS_NB_SIL_INCREASE_FOR_MIN " << nb_sil_increase_KM << endl;
	cout << "\t KMEANS_ELBOW_FACTOR " << fac_elbow_KM << endl;
	cout << "\t KMEANS_NO_SIL_MIN_NO_ELBOW_CHOICE " << after_elbow_choice_KM << endl;
	PrintFootExecMsg();
}

void Displays::Printer_SampleGB_GammaSurface(){
	ReadCurrentFixedParams();
	PrintHeadExecMsg();
	cout << "\t TOL_ORTHO_BOX " << TolOrthoBox << endl;
	cout << "\t TOL_ORTHO_BOX_Z " << TolOrthoBoxZ << endl;
	cout << "\t MIN_BOX_HEIGHT " << MinBoxHeight << endl;
	cout << "\t NB_MAX_LC " << CLsearch << endl;
	cout << "\t MOTIF_SHIFT " << shift_x << " " << shift_y << " " << shift_z << endl;
	cout << "\t MAX_MISFIT " << MaxMisfit << endl;
	cout << "\t MAX_DUP " << MaxDup << endl;
	cout << "\t GB_SPACE " << GBspace << endl;
	cout << "\t FULL_GRAINS " << FullGrains << endl;
	cout << "\t THETA_MAX_ROT_ANGLE_RAT " << theta_max_rot_ax_rat << endl;
	cout << "\t MAX_HKL_ROT_ANGLE " << MaxHKL_rot_angle_rat << endl;
	cout << "\t TOL_DIST_RAT_ANGLE " << tol_dist_rot_angle_rat << endl;
	cout << "\t TOL_POS_KNOWN_CSL_VEC " << tolpos_known_CSL << endl;
	cout << "\t TOL_CSL_INTEGER " << tol_CSL_integer << endl;
	cout << "\t TOL_ALIGNMENT_CSL " << tolAlignment_CSL << endl;
	PrintFootExecMsg();
}


