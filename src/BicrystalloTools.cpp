//**********************************************************************************
//*   BicrystalloTools.cpp                                                         *
//**********************************************************************************
//* This file contains the implementation of the BicrystalloTools class            *
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

#include "AtomHicConfig.h"
#include <iostream>
#include <vector>
#include "BicrystalloTools.h"
#include <omp.h>
#include <MyStructs.h>

using namespace std;

BicrystalloTools::BicrystalloTools(){
	this->MT = new MathTools;
	read_params();
	current_a1 = new double[9];
	current_a2 = new double[9];
	current_CSL = new double[9];
	current_DSC = new double[9];
	for(unsigned int i=0;i<9;i++){
		current_a1[i] = 0.;
		current_a2[i] = 0.;
		current_CSL[i] = 0.;
		current_DSC[i] = 0.;
	}
}

int BicrystalloTools::ComputeCSL_DSC(Crystal *Cryst1, Crystal *Cryst2, double *CSL_Basis, double *DSC_Basis, unsigned int verbose){
	searchCSL(Cryst1,Cryst2,CSL_Basis,verbose);
	int sigma_csl = current_sigma;
	searchDSC(Cryst1,Cryst2,DSC_Basis,verbose);
	int sigma_dsc = current_sigma;
	current_sigma = EndComputeCSL_DSC(sigma_csl,sigma_dsc,CSL_Basis,DSC_Basis);
	Check_CurrentBasis();
	return current_sigma;
}

int BicrystalloTools::ComputeCSL_DSC(Crystal *Cryst, double *rot_mat, double *CSL_Basis, double *DSC_Basis, unsigned int verbose){
	searchCSL(Cryst,rot_mat,CSL_Basis,verbose);
	int sigma_csl = current_sigma;
	searchDSC(Cryst,rot_mat,DSC_Basis,verbose);
	int sigma_dsc = current_sigma;
	current_sigma = EndComputeCSL_DSC(sigma_csl,sigma_dsc,CSL_Basis,DSC_Basis);
	Check_CurrentBasis();
	return current_sigma;
}

int BicrystalloTools::ComputeCSL_DSC(Crystal *Cryst, int h_a, int k_a, int l_a, double theta, double *CSL_Basis, double *DSC_Basis, unsigned int verbose){
	searchCSL(Cryst,h_a,k_a,l_a,theta,CSL_Basis,verbose);
	int sigma_csl = current_sigma;
	searchDSC(Cryst,h_a,k_a,l_a,theta,DSC_Basis,verbose);
	int sigma_dsc = current_sigma;
	current_sigma = EndComputeCSL_DSC(sigma_csl,sigma_dsc,CSL_Basis,DSC_Basis);
	Check_CurrentBasis();
	return current_sigma;
}

int BicrystalloTools::ComputeCSL_DSC(Crystal *Cryst, double *rot_ax_func, double theta_func, double *CSL_Basis, double *DSC_Basis, unsigned int verbose){
	searchCSL(Cryst,rot_ax_func,theta_func,CSL_Basis,verbose);
	int sigma_csl = current_sigma;
	searchDSC(Cryst,rot_ax_func,theta_func,DSC_Basis,verbose);
	int sigma_dsc = current_sigma;
	current_sigma = EndComputeCSL_DSC(sigma_csl,sigma_dsc,CSL_Basis,DSC_Basis);
	Check_CurrentBasis();
	return current_sigma;
}

int BicrystalloTools::ComputeCSL_DSC(double *a1, double *a2, double *CSL_Basis, double *DSC_Basis, unsigned int verbose){
	searchCSL(a1,a2,CSL_Basis,verbose);
	int sigma_csl = current_sigma;
	searchDSC(a1,a2,DSC_Basis,verbose);
	int sigma_dsc = current_sigma;
	current_sigma = EndComputeCSL_DSC(sigma_csl,sigma_dsc,CSL_Basis,DSC_Basis);
	Check_CurrentBasis();
	return current_sigma;
}

// Searching CSL methods
bool BicrystalloTools::searchCSL(Crystal *Cryst1, Crystal *Cryst2, double *CSL_Basis, unsigned int verbose){
	double *a1 = new double[9];
	double *a2 = new double[9];
	for(unsigned int i=0;i<3;i++){
        	a1[i*3] = Cryst1->getA1()[i];
        	a1[i*3+1] = Cryst1->getA2()[i];
        	a1[i*3+2] = Cryst1->getA3()[i];
        	a2[i*3] = Cryst2->getA1()[i];
        	a2[i*3+1] = Cryst2->getA2()[i];
        	a2[i*3+2] = Cryst2->getA3()[i];
	}
	bool success = searchCSL(a1,a2,CSL_Basis,verbose);
	SaveLattices(a1,a2);
	delete[] a1;
	delete[] a2;
	return success;
}

bool BicrystalloTools::searchCSL(Crystal *Cryst, double *rot_mat, double *CSL_Basis, unsigned int verbose){
	double *a1 = new double[9];
	double *a2 = new double[9];
	for(unsigned int i=0;i<3;i++){
        	a1[i*3] = Cryst->getA1()[i];
        	a1[i*3+1] = Cryst->getA2()[i];
        	a1[i*3+2] = Cryst->getA3()[i];
	}
        MT->MatDotMat(rot_mat,a1,a2);
	bool success = searchCSL(a1,a2,CSL_Basis,verbose);
	SaveLattices(a1,a2);
	delete[] a1;
	delete[] a2;
	return success;
}

bool BicrystalloTools::searchCSL(Crystal *Cryst, int h_a, int k_a, int l_a, double theta, double *CSL_Basis, unsigned int verbose){
	double *a1 = new double[9];
	double *a2 = new double[9];
	double *rot_mat = new double[9];
	for(unsigned int i=0;i<3;i++){
        	a1[i*3] = Cryst->getA1()[i];
        	a1[i*3+1] = Cryst->getA2()[i];
        	a1[i*3+2] = Cryst->getA3()[i];
	}
	Cryst->ComputeRotMat(h_a,k_a,l_a,theta,rot_mat);
        MT->MatDotMat(rot_mat,a1,a2);
	bool success = searchCSL(a1,a2,CSL_Basis,verbose);
	SaveLattices(a1,a2);
	delete[] a1;
	delete[] a2;
	delete[] rot_mat;
	return success;
}

bool BicrystalloTools::searchCSL(Crystal *Cryst, double *rot_ax_func, double theta_func, double *CSL_Basis, unsigned int verbose){
	double *a1 = new double[9];
	double *a2 = new double[9];
	double *rot_mat = new double[9];
	for(unsigned int i=0;i<3;i++){
        	a1[i*3] = Cryst->getA1()[i];
        	a1[i*3+1] = Cryst->getA2()[i];
        	a1[i*3+2] = Cryst->getA3()[i];
	}
        MT->Vec2rotMat(rot_ax_func,theta_func,rot_mat);
        MT->MatDotMat(rot_mat,a1,a2);
	bool success = searchCSL(a1,a2,CSL_Basis,verbose);
	SaveLattices(a1,a2);
	delete[] a1;
	delete[] a2;
	delete[] rot_mat;
	return success;
}

// Searching DSC methods
bool BicrystalloTools::searchDSC(Crystal *Cryst1, Crystal *Cryst2, double *DSC_Basis, unsigned int verbose){
	double *a1 = new double[9];
	double *a2 = new double[9];
	for(unsigned int i=0;i<3;i++){
        	a1[i*3] = Cryst1->getA1()[i];
        	a1[i*3+1] = Cryst1->getA2()[i];
        	a1[i*3+2] = Cryst1->getA3()[i];
        	a2[i*3] = Cryst2->getA1()[i];
        	a2[i*3+1] = Cryst2->getA2()[i];
        	a2[i*3+2] = Cryst2->getA3()[i];
	}
	bool success = searchDSC(a1,a2,DSC_Basis,verbose);
	SaveLattices(a1,a2);
	delete[] a1;
	delete[] a2;
	return success;
}

bool BicrystalloTools::searchDSC(Crystal *Cryst, double *rot_mat, double *DSC_Basis, unsigned int verbose){
	double *a1 = new double[9];
	double *a2 = new double[9];
	for(unsigned int i=0;i<3;i++){
        	a1[i*3] = Cryst->getA1()[i];
        	a1[i*3+1] = Cryst->getA2()[i];
        	a1[i*3+2] = Cryst->getA3()[i];
	}
        MT->MatDotMat(rot_mat,a1,a2);
	bool success = searchDSC(a1,a2,DSC_Basis,verbose);
	SaveLattices(a1,a2);
	delete[] a1;
	delete[] a2;
	return success;
}

bool BicrystalloTools::searchDSC(Crystal *Cryst, int h_a, int k_a, int l_a, double theta, double *DSC_Basis, unsigned int verbose){
	double *a1 = new double[9];
	double *a2 = new double[9];
	double *rot_mat = new double[9];
	for(unsigned int i=0;i<3;i++){
        	a1[i*3] = Cryst->getA1()[i];
        	a1[i*3+1] = Cryst->getA2()[i];
        	a1[i*3+2] = Cryst->getA3()[i];
	}
	Cryst->ComputeRotMat(h_a,k_a,l_a,theta,rot_mat);
        MT->MatDotMat(rot_mat,a1,a2);
	bool success = searchDSC(a1,a2,DSC_Basis,verbose);
	SaveLattices(a1,a2);
	delete[] a1;
	delete[] a2;
	delete[] rot_mat;
	return success;
}

bool BicrystalloTools::searchDSC(Crystal *Cryst, double *rot_ax_func, double theta_func, double *DSC_Basis, unsigned int verbose){
	double *a1 = new double[9];
	double *a2 = new double[9];
	double *rot_mat = new double[9];
	for(unsigned int i=0;i<3;i++){
        	a1[i*3] = Cryst->getA1()[i];
        	a1[i*3+1] = Cryst->getA2()[i];
        	a1[i*3+2] = Cryst->getA3()[i];
	}
        MT->Vec2rotMat(rot_ax_func,theta_func,rot_mat);
        MT->MatDotMat(rot_mat,a1,a2);
	bool success = searchDSC(a1,a2,DSC_Basis,verbose);
	SaveLattices(a1,a2);
	delete[] a1;
	delete[] a2;
	delete[] rot_mat;
	return success;
}

bool BicrystalloTools::searchDSC(double *a1, double *a2, double *DSC_Basis, unsigned int verbose){
	double *a1_s = new double[9];
	double *a2_s = new double[9];
	MT->ComputeReciproqual(a1,a1_s);
	MT->ComputeReciproqual(a2,a2_s);
	bool success = searchCSL(a1_s,a2_s,DSC_Basis,verbose);
	SaveLattices(a1,a2);
	MT->ComputeReciproqual(DSC_Basis,DSC_Basis);
	delete[] a1_s;
	delete[] a2_s;
	return success;
}

// search CSL primitive lattice based on the work of Xie et al. 2020 and Bonnet and Rolland 1975 (iterative brut force algo)
bool BicrystalloTools::searchCSL(double *a1, double *a2, double *CSL_Basis, unsigned int verbose){
	double tolpos_known_CSL = 5e-1; // these were initially in FixedParameters but they are too difficult to understand
	double tol_CSL_integer = 1e-2;
	double tolAlignment_CSL = 0.1;
	
	double *a1_inv = new double[9]; // Basis vector of lattice 1 
	double *a2_a = new double[9]; // auxiliary lattice for CSL base computation 
	double *U = new double[9]; 
	double *U_inv = new double[9]; 
	double *U_1 = new double[9]; 
	double *Backup_U1 = new double[9]; 
	double *U_1_temp = new double[9]; 
	double *U_2 = new double[9]; 
	double *Ei = new double[9]; 
	double *Ei_inv = new double[9]; 
	double *Fi = new double[9]; 
	double *Fi_inv = new double[9]; 
	//
	double *DSC_Basis = new double[9];
	double *DSC_inv = new double[9]; 
	double *U_a = new double[9];
	double *searchVec = new double[3];
	double *CSL_Basis2 = new double[9]; 
	double *CSL_Basis_temp = new double[9]; 
	double *buffer_mat = new double[9]; 
	int *k = new int[3];
	double *rot_ax_func = new double[3];
	double *R = new double[9];
	MT->invert3x3(a1,a1_inv);
	MT->MatDotMat(a2,a1_inv,R);
	rot_ax_func[0] = R[7]-R[5];
	rot_ax_func[1] = R[2]-R[6];
	rot_ax_func[2] = R[3]-R[1];
	double NormRotAx = 0.;
	for(unsigned int i=0;i<3;i++)
		NormRotAx += pow(rot_ax_func[i],2.);
	NormRotAx = sqrt(NormRotAx);
	for(unsigned int i=0;i<3;i++) rot_ax_func[i] /= NormRotAx;
	double *Known_CSL = new double[3];
	double *buffer_vec = new double[3];
	double tmp(0.), tmp_1(0.), tmp_2(0.), BasePrec(0.);
	for(unsigned int i=0;i<3;i++)
		Known_CSL[i] = 0.;
	unsigned int L, M, N;
	bool BaseFound = false;
	bool IsFindIntVec = false;
	bool BaseAligned = false;
	bool FirstBase = true;
	double tol_DSC = 1e-6;
	double tol_IntVec, sp;
	MT->invert3x3(a1,a1_inv);
	MT->MatDotMat(a1_inv,a2,U);
	MT->invert3x3(U,U_inv);
	// Compute the symmetric of U with respect to second diag
	MT->dia_sym_mtx(U,U_a);
	// compute the auxiliary lattice
	MT->MatDotMat(a1,U_a,a2_a);
	double exp = 5.;
	bool OneBaseFound = false;
	bool success;
	// loop on tol_IntVec (increasing tolerance) until the primitive lattice vector of CSL fall within the box
	while( exp > 0. && ( ( !BaseFound && !OneBaseFound ) || ( BaseFound && OneBaseFound ) ) ){
		tol_IntVec = pow(10.,-exp);
		if( verbose == 2 ) cout << "Trying to find CSL base using a tolerance of : " << tol_IntVec << endl;
		// search if a2_a_1 could be the basic vector for computing DSC lattice
		for(unsigned int i=0;i<3;i++) searchVec[i] = a2_a[i*3];
		MT->MatDotRawVec(a1_inv,searchVec,searchVec);
		L = MT->find_integer_vector(searchVec,tol_IntVec,SigmaMax,k,IsFindIntVec);
		if( !IsFindIntVec ){
			if( verbose == 2 ) cout << "Fail 1 ! Increasing tolerance !" << endl;
			if( exp > 1.5 )	exp -= .5;
			else exp -= .1;
			continue;
		}
		//
		MT->solve_DSC(k,L,a1,Ei,tol_DSC);
		// search if this Basis is also good for a2_a_2
		for(unsigned int i=0;i<3;i++) searchVec[i] = a2_a[i*3+1];
		MT->invert3x3(Ei,Ei_inv);
		MT->MatDotRawVec(Ei_inv,searchVec,searchVec);
		L = MT->find_integer_vector(searchVec,tol_IntVec,SigmaMax,k,IsFindIntVec);
		if( !IsFindIntVec ){
			if( verbose == 2 ) cout << "Fail 2 ! Increasing tolerance !" << endl;
			if( exp > 1.5 )	exp -= .5;
			else exp -= .1;
			continue;
		}
		if( L == 1 ){
			// Yes, test if the Basis is also good for a2_a_3
			for(unsigned int i=0;i<3;i++) searchVec[i] = a2_a[i*3+2];
			MT->MatDotRawVec(Ei_inv,searchVec,searchVec);
			M = MT->find_integer_vector(searchVec,tol_IntVec,SigmaMax,k,IsFindIntVec);
			if( !IsFindIntVec ){
				if( verbose == 2 ) cout << "Fail 3 ! Increasing tolerance !" << endl;
				if( exp > 1.5 )	exp -= .5;
				else exp -= .1;
				continue;
			}
			if( M != 1 ) MT->solve_DSC(k,M,Ei,DSC_Basis,tol_DSC); // No, use this vector for DSC Basis
			else for(unsigned int i=0;i<9;i++) DSC_Basis[i] = Ei[i];
		}else{
			// No, use this vector for DSC Basis
			MT->solve_DSC(k,L,Ei,Fi,tol_DSC);
			// search if this Basis is also good for a2_a_2
			for(unsigned int i=0;i<3;i++) searchVec[i] = a2_a[i*3+2];
			MT->invert3x3(Fi,Fi_inv);
			MT->MatDotRawVec(Fi_inv,searchVec,searchVec);
			M = MT->find_integer_vector(searchVec,tol_IntVec,SigmaMax,k,IsFindIntVec);
			if( !IsFindIntVec ){
				if( verbose == 2 ) cout << "Fail 4 ! Increasing tolerance !" << endl;
				if( exp > 1.5 )	exp -= .5;
				else exp -= .1;
				continue;
			}
			if( M != 1 ){
				MT->solve_DSC(k,M,Fi,DSC_Basis,tol_DSC); // No, use this vector for DSC Basis
			}
			else for(unsigned int i=0;i<9;i++) DSC_Basis[i] = Fi[i];
		}
		MT->get_right_hand(DSC_Basis,DSC_Basis);
		MT->LLL(DSC_Basis,DSC_Basis);
		MT->invert3x3(DSC_Basis,DSC_inv);
		MT->MatDotMat(DSC_inv,a2_a,U_1_temp);
		MT->dia_sym_mtx(U_1_temp,U_1_temp);
		MT->MatDotMat(U_inv,U_1_temp,U_2);
		MT->MatDotMat(a1,U_1_temp,CSL_Basis_temp);
		MT->MatDotMat(a2,U_2,CSL_Basis2);
		MT->LLL(CSL_Basis_temp,CSL_Basis_temp);
		MT->get_right_hand(CSL_Basis_temp,CSL_Basis_temp);
		MT->LLL(CSL_Basis2,CSL_Basis2);
		MT->get_right_hand(CSL_Basis2,CSL_Basis2);
		if( CSL_Basis_temp[0] == CSL_Basis_temp[0] ){
			// keep the base which permits to get the known CSL vector and with the higher tolIntVec and which gives a sigma value higher than 1
			this->current_sigma = round(fabs(MT->det(U_1_temp)));
			if( verbose == 1 ){
				cout << "We have found CSL Basis corresponding to Sigma = " << this->current_sigma << endl;
				MT->printMat(U_1_temp);
			}
			if( OneBaseFound && this->current_sigma < 1.5 ){
				BaseFound = true;
				break;
			}
			BaseFound = false;
			MT->invert3x3(CSL_Basis_temp,buffer_mat);
			MT->MatDotRawVec(buffer_mat,Known_CSL,buffer_vec);
			tmp = 0.;
			tmp_1 = 0.;
			tmp_2 = 0.;
			for(unsigned int i=0;i<9;i++) tmp_2 += fabs(U_1_temp[i]-round(U_1_temp[i]));
			for(unsigned int i=0;i<3;i++){
				tmp += fabs(buffer_vec[i]-round(buffer_vec[i]));
				tmp_1 += fabs(round(buffer_vec[i]));
			}
			if( ( tmp < tolpos_known_CSL && tmp_1 > tolpos_known_CSL && tmp_2 < tol_CSL_integer ) || tmp_2 < tol_DSC ){
				BaseFound = true;
				OneBaseFound = true;
				if( verbose == 1 ) cout << "Good base found !" << endl;
				for(unsigned int ii=0;ii<9;ii++){
					CSL_Basis[ii] = CSL_Basis_temp[ii];
					U_1[ii] = U_1_temp[ii];
				}
			}
			// one of the CSL Basis vector should be aligned with the rotation axis => verify that the backup base satisfy this criterium
			for(unsigned int v=0;v<3;v++){
				sp = 0.;
				NormRotAx = 0.;
				BaseAligned = false;
				for(unsigned int i=0;i<3;i++){
					sp += CSL_Basis_temp[i*3+v]*rot_ax_func[i];
					NormRotAx += pow(CSL_Basis_temp[i*3+v],2.);
				}
				if( fabs(fabs((sp/sqrt(NormRotAx)))-1.) < tolAlignment_CSL ){
					BaseAligned = true;
					break;
				}
			}
			if( FirstBase && this->current_sigma < SigmaMax && BaseAligned ){
				BasePrec = tmp_2;
				for(unsigned int i=0;i<9;i++) Backup_U1[i] = U_1_temp[i];
				FirstBase = false;
			}else if( tmp_2 < BasePrec && this->current_sigma < SigmaMax && BaseAligned ){
				BasePrec = tmp_2;
				for(unsigned int i=0;i<9;i++) Backup_U1[i] = U_1_temp[i];
			}
			if( !BaseFound && OneBaseFound ){
				BaseFound = true;
				break;
			}
		}else{
			if( verbose == 2 ) cout << "Fail ! Increasing tolerance !" << endl;
			if( exp > 1.5 )	exp -= .5;
			else exp -= .1;
		}
		if( ( !BaseFound && !OneBaseFound ) || ( BaseFound && OneBaseFound ) ){
			if( exp > 1.5 )	exp -= .5;
			else exp -= .1;
		}
	}
	if( BaseFound ){
		this->current_sigma = round(fabs(MT->det(U_1)));
		//for(unsigned int i=0;i<9;i++){// still something TODO for sharing the CSL to the two lattices (for the moment CSL coincide with grain 2) 
		//	CSL_Basis[i] += CSL_Basis2[i];
		//	CSL_Basis[i] /= 2.;
		//}
		//cout << "Success ! We find a CSL corresponding to sigma =" << this->sigma << ", using a tolerance of " << tol_IntVec << endl;
		//cout << "This misorientation corresponds to a GB with sigma = " << this->current_sigma << endl;
		success=true;
	}else if( !FirstBase ){
		MT->MatDotMat(a1,Backup_U1,CSL_Basis);
		this->current_sigma = round(fabs(MT->det(Backup_U1)));
		double dev_CSL = 0.;
		for(unsigned int i=0;i<9;i++) dev_CSL += fabs(round(Backup_U1[i]/1.)-(Backup_U1[i]/1.));
		//MT->printMat(Backup_U1);
		cout << "We find a CSL Basis corresponding to sigma = " << this->current_sigma << endl;
		cout << "Nevertheless the deviation from a perfect CSL lattice is : " << dev_CSL  << " (sum of differences from integer of all elements of U)" << endl;
		success=true;
	}else{
		cout << "Fail to find CSL Basis" << endl;
		success=false;
	}
	delete[] a1_inv;
	delete[] a2_a;
	delete[] R;
	delete[] U;
	delete[] U_inv;
	delete[] U_1;
	delete[] Backup_U1;
	delete[] U_1_temp;
	delete[] U_2;
	delete[] Ei;
	delete[] Ei_inv;
	delete[] Fi;
	delete[] Fi_inv;
	delete[] DSC_Basis;
	delete[] DSC_inv;
	delete[] U_a;
	delete[] searchVec;
	delete[] CSL_Basis2;
	delete[] CSL_Basis_temp;
	delete[] buffer_mat;
	delete[] k;
	delete[] rot_ax_func;
	delete[] Known_CSL;
	delete[] buffer_vec;
	return success;
}

void BicrystalloTools::generateCSL(){
	// Construct motif of CSL, DSC, lattice 1 and 2 points to create a crystal object
	vector<Atom> Motif_v;
	string *AtomType = new string[4];
	unsigned int *AtomType_uint = new unsigned int[4];
	for(unsigned int i=0;i<4;i++) AtomType_uint[i] = i+1;
	// CSL => only one point at the origin
	Motif_v.push_back({1,0,0,0});
	AtomType[0] = "CSL";
	AtomType[1] = "DSC";
	AtomType[2] = "Lattice 1";
	AtomType[3] = "Lattice 2";
	// lattice 1 and 2
	unsigned int max_cl = 100;
	int cl = current_sigma;
	if( cl > max_cl ) cl = max_cl;
	double *cur_pos = new double[3];
	double *zero_vec = new double[3];
	for(unsigned int d=0;d<3;d++) zero_vec[d] = 0.;
	for(int i1=-cl;i1<cl+1;i1++){
		for(int i2=-cl;i2<cl+1;i2++){
			for(int i3=-cl;i3<cl+1;i3++){
				for(unsigned int d=0;d<3;d++) cur_pos[d] = current_a1[3*d]*i1 + current_a1[3*d+1]*i2 + current_a1[3*d+2]*i3;
				if( MT->IsInParallelepiped(current_CSL,zero_vec,cur_pos) ) Motif_v.push_back({3,cur_pos[0],cur_pos[1],cur_pos[2]});
				for(unsigned int d=0;d<3;d++) cur_pos[d] = current_a2[3*d]*i1 + current_a2[3*d+1]*i2 + current_a2[3*d+2]*i3;
				if( MT->IsInParallelepiped(current_CSL,zero_vec,cur_pos) ) Motif_v.push_back({4,cur_pos[0],cur_pos[1],cur_pos[2]});
			}
		}
	}
	// dsc
	cl = current_sigma*current_sigma;
	if( cl > max_cl ) cl = max_cl;
	for(int i1=-cl;i1<cl+1;i1++){
		for(int i2=-cl;i2<cl+1;i2++){
			for(int i3=-cl;i3<cl+1;i3++){
				for(unsigned int d=0;d<3;d++) cur_pos[d] = current_DSC[3*d]*i1 + current_DSC[3*d+1]*i2 + current_DSC[3*d+2]*i3;
				if( MT->IsInParallelepiped(current_CSL,zero_vec,cur_pos) ) Motif_v.push_back({2,cur_pos[0],cur_pos[1],cur_pos[2]});
			}
		}
	}
	unsigned int nbAtom = Motif_v.size();
	Atom *Motif = new Atom[nbAtom];
	for(unsigned int n=0;n<nbAtom;n++) Motif[n] = Motif_v[n];
	if( MyCSL_Crystal ) delete MyCSL_Crystal;
	MyCSL_Crystal = new Crystal(current_CSL,Motif,nbAtom,4,AtomType,AtomType_uint);
	MyCSL_Crystal->ConstructOrthogonalCell(true);
	MyCSL_Crystal->getOrientedSystem()->print_lmp("CSL.lmp");
	delete[] AtomType;
	delete[] AtomType_uint;
	delete[] cur_pos;
	delete[] zero_vec;
	delete[] Motif;
}

//void Bicrystal::printCSL(const std::string filename){
//	if( !IsCSL ){
//		cout << "CSL has not been generated, cannot print it" << endl;
//		return;
//	}
//	ofstream writefile(filename);
//	writefile << " # File generated using AtomHic\n";
//	writefile << this->File_Heading;
//	// TEST
//	double l1 = _MyCrystal->getOrientedSystem()->getH1()[0];
//	double l2 = _MyCrystal->getOrientedSystem()->getH2()[1];
//	double l3 = _MyCrystal->getOrientedSystem()->getH3()[2];
//        writefile << "\n\t" << NodesG1.size()+NodesG2.size()+CSL.size() << "\tatoms\n\t" << 3 << "\tatom types\n\n\t0.000000000000\t" << l1 << "\txlo xhi\n\t0.000000000000\t" << l2 << "\tylo yhi\n\t0.000000000000\t" << l3 << "\tzlo zhi\n";
//	// END
//        //writefile << "\n\t" << NodesG1.size()+NodesG2.size()+CSL.size() << "\tatoms\n\t" << 3 << "\tatom types\n\n\t0.000000000000\t" << H1[0] << "\txlo xhi\n\t0.000000000000\t" << H2[1] << "\tylo yhi\n\t0.000000000000\t" << H3[2] << "\tzlo zhi\n";
//	if( this->IsTilted ) writefile << "\t" << H2[0] << "\t" << H3[0] << "\t" << H3[1] << "\txy xz yz\n";
//	writefile << "\nMasses\n\n";
//	writefile << "\t" << 1 << "\t0\t# Lattice1\n";
//	writefile << "\t" << 2 << "\t0\t# Lattice2\n";
//	writefile << "\t" << 3 << "\t0\t# CSL\n";
//	writefile << "\nAtoms # atomic\n\n";
//
//	for(unsigned int i=0;i<NodesG1.size();i++) writefile << i << "\t" << 1 << "\t" << NodesG1[i][0] << "\t" << NodesG1[i][1] << "\t" << NodesG1[i][2] << "\n";
//	for(unsigned int i=0;i<NodesG2.size();i++) writefile << NodesG1.size()+1+i << "\t" << 2 << "\t" << NodesG2[i][0] << "\t" << NodesG2[i][1] << "\t" << NodesG2[i][2] << "\n";
//	for(unsigned int i=0;i<CSL.size();i++) writefile << NodesG1.size()+NodesG2.size()+1+i << "\t" << 3 << "\t" << CSL[i][0] << "\t" << CSL[i][1] << "\t" << CSL[i][2] << "\n";
//	writefile.close();
//	cout << "File " << filename << " successfully writted !" << endl;
//}
//
//// search the rotation axis and angle corresponding to the closest rational GB // TODO warning will not work for cubic primitive cell
double BicrystalloTools::RationalizeOri(Crystal *Cryst, int h_a_func, int k_a_func, int l_a_func, double theta_func, int &hp_near, int &kp_near, int &lp_near, int *CSL_vec){
	cout << "Rationalizing the rotation axis and angle" << endl;
	double sp, norm1, norm2, theta_test;
	vector<double> theta_list;
	vector<int> hkl_list;
	double tolScalarProd = 1e-3;
	double *vec1 = new double[3];
	double *vec2 = new double[3];
	double *vec3 = new double[3];
	double *vec4 = new double[3];
	double *vec5 = new double[3];
	double *rot_ax_func = new double[3];
	bool DoNotRatAxis = false;
	int TrueMaxCL = 150; // TODO parameters	
	// rationalize the rotation axis by searching the closest (with quasi the same orientation) normal to a crystal plane (with relatively low integer Miller indices)
	int arr[3] = {abs(h_a_func),abs(k_a_func),abs(l_a_func)};
	int MaxHKL_Norm;
	if( MT->max(arr,3) < 3 && ( ( h_a_func == 0 && k_a_func == 0 ) || ( h_a_func == 0 && l_a_func == 0 ) || ( k_a_func == 0 && l_a_func == 0 ) ) ) MaxHKL_Norm = MT->max(arr,3)*5;
	else if( MT->max(arr,3) < 3 ) MaxHKL_Norm = MT->max(arr,3)*100;
	else if( MT->max(arr,3) < 10 ) MaxHKL_Norm = MT->max(arr,3)*50;
	else MaxHKL_Norm = MT->max(arr,3)*10;
	//cout << MaxHKL_Norm << endl;
	if( ( h_a_func == 0 && k_a_func == 0 ) || ( h_a_func == 0 && l_a_func == 0 ) || ( k_a_func == 0 && l_a_func == 0 ) ) DoNotRatAxis = true;
	if( DoNotRatAxis || Cryst->getCrystallo() == "Cubic" || Cryst->getCrystallo() == "BCC" || Cryst->getCrystallo() == "FCC" ){// for cubic crystals or simple rot axis, the scalar product between plane and direction vector is expressed simply
		hp_near = h_a_func;
		kp_near = k_a_func;
		lp_near = l_a_func;
	}else{
		norm1 = 0.;
		for(unsigned int i=0;i<3;i++){
			vec1[i] = Cryst->getA1()[i]*h_a_func + Cryst->getA2()[i]*k_a_func + Cryst->getA3()[i]*l_a_func;
			norm1 += pow(vec1[i],2.);
		}
		norm1 = sqrt(norm1);
		for(int i=-MaxHKL_Norm;i<MaxHKL_Norm+1;i++){
			for(int j=-MaxHKL_Norm;j<MaxHKL_Norm+1;j++){
				for(int k=-MaxHKL_Norm;k<MaxHKL_Norm+1;k++){
					sp = 0.;
					norm2 = 0.;
					for(unsigned int ii=0;ii<3;ii++){
						vec2[ii] = Cryst->getA1_star()[ii]*i + Cryst->getA2_star()[ii]*j + Cryst->getA3_star()[ii]*k;
						sp += vec2[ii]*vec1[ii];
						norm2 += pow(vec2[ii],2.);
					}
					theta_test = acos(sp/(norm1*sqrt(norm2)));
					if( theta_test < this->theta_max_rot_ax_rat ){
						theta_list.push_back(theta_test);
						hkl_list.push_back(i);
						hkl_list.push_back(j);
						hkl_list.push_back(k);
					}
				}
			}
		}
		unsigned int ind_thetalist;
		if( theta_list.size() != 0 ){
			ind_thetalist = MT->min(theta_list);
			hp_near = hkl_list[ind_thetalist*3];
			kp_near = hkl_list[ind_thetalist*3+1];
			lp_near = hkl_list[ind_thetalist*3+2];
			cout << "The rotation axis has been rationalized to be the normal of the (" << hp_near << "_" << kp_near << "_" << lp_near <<") plane" << endl;
		}else{
			cerr << "We don't have found a close normal plane vector to the rotation axis for rationalizing the GB, aborting calculation" << endl;
			exit(EXIT_FAILURE);
		}
	}
	
	// compute the new rotation axis
	Cryst->ComputeVec(hp_near,kp_near,lp_near,rot_ax_func,false);
	
	// Now hp_near, kp_near and lp_near are plane indices which permits to simplify scalar product and the searching of vector normal to rotation axis	
	// get cell parameter from crystal
	double a_c = Cryst->getALength()[0];
	double b_c = Cryst->getALength()[1];
	double c_c = Cryst->getALength()[2];
	vector<double> Normal1;
	vector<double> Normal2;
	vector<int> Normal_ind;
	unsigned int ind_norm;
	int *int_vec = new int[3];
	unsigned int *uint_vec = new unsigned int[3];
	double norm;
	bool ToStore;
	// search the two smallest normal direction to the rotation axis with integer Miller indices
	// loop for the first one
	for(int i=-MaxHKL_Norm;i<MaxHKL_Norm+1;i++){
		for(int j=-MaxHKL_Norm;j<MaxHKL_Norm+1;j++){
			for(int k=-MaxHKL_Norm;k<MaxHKL_Norm+1;k++){
				if( ( abs(i*hp_near+j*kp_near+k*lp_near) < tolScalarProd ) && ( i != 0 || j != 0 || k!= 0 ) ){
					ToStore = true; // do not store colinear vectors
					norm = sqrt(pow(a_c*i,2.)+pow(b_c*j,2.)+pow(c_c*k,2.));
					for(unsigned int n=0;n<Normal1.size()/4;n++){
						sp = ( i*Normal1[n*4+1]*pow(a_c,2.) + j*Normal1[n*4+2]*pow(b_c,2.) + k*Normal1[n*4+3]*pow(c_c,2.) ) / ( norm*Normal1[n*4] );
						if( fabs(sp-1.) < 1e-3 ){
							ToStore = false;
							break;
						}
					}
					if( ToStore ){
						int_vec[0] = i;
						int_vec[1] = j;
						int_vec[2] = k;
						MT->reduce_vec(int_vec,int_vec);
						Normal1.push_back(sqrt(pow(a_c*int_vec[0],2.)+pow(b_c*int_vec[1],2.)+pow(c_c*int_vec[2],2.)));
						Normal1.push_back(int_vec[0]);
						Normal1.push_back(int_vec[1]);
						Normal1.push_back(int_vec[2]);
					}
				}
			}
		}
	}
	if( Normal1.size() == 0 ){
		cerr << "We don't find normal direction to the roation axis" << endl;
		exit(EXIT_FAILURE);
	}
	double DeltaTheta_orth = M_PI/12.;
	double cos_lo = cos((M_PI/2.)+DeltaTheta_orth);
	double cos_hi = cos((M_PI/2.)-DeltaTheta_orth);
        for(unsigned int i=0;i<Normal1.size()/4;i++){
        	for(unsigned int j=0;j<Normal1.size()/4;j++){
			sp = ( Normal1[i*4+1]*Normal1[j*4+1]*pow(a_c,2.) + Normal1[i*4+2]*Normal1[j*4+2]*pow(b_c,2.) + Normal1[i*4+3]*Normal1[j*4+3]*pow(c_c,2.) ) / ( Normal1[i*4]*Normal1[j*4] );
		        if( sp > cos_lo && sp < cos_hi ){
				Normal2.push_back(0.); // criteria for choosing this => small Miller indices (not for the moment + angle between vectors near 90°)
				for(unsigned int ii=0;ii<3;ii++) Normal2[Normal2.size()-1] += abs(Normal1[i*4+ii+1]) + abs(Normal1[j*4+ii+1]);
				Normal_ind.push_back(round(Normal1[i*4+1]));
				Normal_ind.push_back(round(Normal1[i*4+2]));
				Normal_ind.push_back(round(Normal1[i*4+3]));
				Normal_ind.push_back(round(Normal1[j*4+1]));
				Normal_ind.push_back(round(Normal1[j*4+2]));
				Normal_ind.push_back(round(Normal1[j*4+3]));
			}
		}
	}	

	if( Normal2.size() == 0 ){
		cerr << "We don't find normal direction to the rotation axis" << endl;
		exit(EXIT_FAILURE);
	}
	// get the smallest ones
	ind_norm = MT->min(Normal2);
	// change or not the sign of the second vector to have a direct Basis
	Cryst->ComputeVec(Normal_ind[ind_norm*6],Normal_ind[ind_norm*6+1],Normal_ind[ind_norm*6+2],vec1,true,false);
	Cryst->ComputeVec(Normal_ind[ind_norm*6+3],Normal_ind[ind_norm*6+4],Normal_ind[ind_norm*6+5],vec2,true,false);
	MT->crossProd(vec1,vec2,vec3);
	sp = 0.;
	for(unsigned int i=0;i<3;i++) sp += rot_ax_func[i]*vec3[i];
	if( sp < 0 ) for(unsigned int i=0;i<3;i++) Normal_ind[ind_norm*6+3+i] *= -1;
	double theta_star;
	unsigned int ind_rattheta;
	int max_i, max_j;
	for(unsigned int i=0;i<3;i++) arr[i] = abs(Normal_ind[ind_norm*6+3+i]);
	max_j = MT->max(arr,3);
	for(unsigned int i=0;i<3;i++) arr[i] = abs(Normal_ind[ind_norm*6+i]);
	max_i = MT->max(arr,3);
	int MaxHKL_i, MaxHKL_j;
	if( max_i >= max_j ){
		MaxHKL_i = this->MaxHKL_rot_angle_rat;
		MaxHKL_j = this->MaxHKL_rot_angle_rat*max_i/max_j;
		if( MaxHKL_j > 100 ) MaxHKL_j = 100;
	}else{
		MaxHKL_j = this->MaxHKL_rot_angle_rat;
		MaxHKL_i = this->MaxHKL_rot_angle_rat*max_j/max_i;
		if( MaxHKL_i > 100 ) MaxHKL_i = 100;
	}
	// generate lattices of the two vectors for the two orientations
	vector<double> Lat1; // contains distance and x y z coordinates
	vector<double> DeltaTheta_temp;
	vector<double> temp_vec;
	vector<int> indices;
	double *RotMat = new double[9];
	MT->Vec2rotMat(rot_ax_func,theta_func,RotMat);
	for(int i=-MaxHKL_i;i<MaxHKL_i+1;i++){
		for(int j=-MaxHKL_j;j<MaxHKL_j+1;j++){
			norm = 0.;
			for(unsigned int ii=0;ii<3;ii++){
				vec3[ii] = vec1[ii]*i + vec2[ii]*j;
				norm += pow(vec3[ii],2.);
			}
			if( norm != 0. ){
				Lat1.push_back(sqrt(norm));
				for(unsigned int ii=0;ii<3;ii++) Lat1.push_back(vec3[ii]);
				indices.push_back(i);
				indices.push_back(j);
			}
		}
	}
	// search the ones with the same distance
	for(unsigned int i=0;i<Lat1.size()/4;i++){
		temp_vec.clear();
		for(unsigned int j=0;j<Lat1.size()/4;j++){
			if( fabs(Lat1[i*4]-Lat1[j*4]) < this->tol_dist_rot_angle_rat ){
				// rotate the vector to give the corresponding crystal 2 vector
				for(unsigned int ii=0;ii<3;ii++){
					vec3[ii] = Lat1[i*4+ii+1];
					vec4[ii] = Lat1[j*4+ii+1];
				}
				MT->MatDotRawVec(RotMat,vec4,vec4);
				// compare the angles between the two vector and store the difference
				sp = 0.;
				for(unsigned int ii=0;ii<3;ii++) sp += vec3[ii]*vec4[ii];
				theta_star = acos(sp/(Lat1[i*4]*Lat1[j*4]));
				// use cross product for the sign of the angle
				MT->crossProd(vec3,vec4,vec5);
				sp = 0.;
				for(unsigned int ii=0;ii<3;ii++) sp += vec5[ii]*rot_ax_func[ii];
				if( sp <= 0 ) temp_vec.push_back(theta_star);
				else temp_vec.push_back(-theta_star);
			}
		}
		if( temp_vec.size() != 0 ) DeltaTheta_temp.push_back(MT->min_vec_abs(temp_vec));
		else DeltaTheta_temp.push_back(M_PI);
		DeltaTheta_temp.push_back(fabs(indices[i*2])+fabs(indices[i*2+1]));
		DeltaTheta_temp.push_back(indices[i*2]);
		DeltaTheta_temp.push_back(indices[i*2+1]);
	}
	MT->sort_abs(DeltaTheta_temp,0,4,DeltaTheta_temp);
	double min = DeltaTheta_temp[0];
	vector<double> DeltaTheta;
	vector<double> true_DeltaTheta;
	double tolSelect = 2.5e-2; // in rads
	indices.clear();
	DeltaTheta.push_back(DeltaTheta_temp[1]);
	true_DeltaTheta.push_back(DeltaTheta_temp[0]);
	indices.push_back(DeltaTheta_temp[2]);
	indices.push_back(DeltaTheta_temp[3]);
	for(unsigned int i=0;i<(DeltaTheta_temp.size()/4)-1;i++){
		if( fabs(fabs(DeltaTheta_temp[i*4])-fabs(min)) < tolSelect ){
			DeltaTheta.push_back(DeltaTheta_temp[(i+1)*4+1]);
			true_DeltaTheta.push_back(DeltaTheta_temp[(i+1)*4]);
			indices.push_back(DeltaTheta_temp[(i+1)*4+2]);
			indices.push_back(DeltaTheta_temp[(i+1)*4+3]);
		}
	}
	ind_rattheta = MT->min(DeltaTheta);
	if( fabs(true_DeltaTheta[ind_rattheta]+theta_func) < 1e-3 ){
		cerr << "We don't have found rational angle" << endl;
		exit(EXIT_FAILURE);
	}	
	
	CSL_vec[0] = Normal_ind[ind_norm*6]  *indices[ind_rattheta*2] + Normal_ind[ind_norm*6+3]*indices[ind_rattheta*2+1];
	CSL_vec[1] = Normal_ind[ind_norm*6+1]*indices[ind_rattheta*2] + Normal_ind[ind_norm*6+4]*indices[ind_rattheta*2+1];
	CSL_vec[2] = Normal_ind[ind_norm*6+2]*indices[ind_rattheta*2] + Normal_ind[ind_norm*6+5]*indices[ind_rattheta*2+1];
	
	// search the GB plane corresponding to a STGB for this rotation angle
	vector<double> normSTGB_true;
	vector<int> MilSTGB_true;
	vector<double> normSTGB_near;
	vector<int> MilSTGB_near;
	for(unsigned int i=0;i<3;i++) uint_vec[i] = abs(CSL_vec[i]);
	if( MT->max_p(uint_vec,3) < 100 ){
		int CL_STGB = (int) MT->max_p(uint_vec,3)*2;
		for(int i=-CL_STGB;i<CL_STGB;i++){
			for(int j=-CL_STGB;j<CL_STGB;j++){
				for(int k=-CL_STGB;k<CL_STGB;k++){
					Cryst->ComputeVec(i,j,k,vec1,false,false);
					sp = 0.;
					for(unsigned int ii=0;ii<3;ii++) sp += vec1[ii]*rot_ax_func[ii];
					if( abs(i*CSL_vec[0]+j*CSL_vec[1]+k*CSL_vec[2]) < tolScalarProd && fabs(sp) < tolScalarProd && ( i != 0 || j != 0 || k!= 0 ) ){
						normSTGB_true.push_back(pow(i,2.)+pow(j,2.)+pow(k,2.));
						MilSTGB_true.push_back(i);
						MilSTGB_true.push_back(j);
						MilSTGB_true.push_back(k);
					}else if( abs(i*CSL_vec[0]+j*CSL_vec[1]+k*CSL_vec[2]) < tolScalarProd && abs(i*h_a_func+j*k_a_func+k*l_a_func) < tolScalarProd && ( i != 0 || j != 0 || k!= 0 ) ){
						normSTGB_near.push_back(pow(i,2.)+pow(j,2.)+pow(k,2.));
						MilSTGB_near.push_back(i);
						MilSTGB_near.push_back(j);
						MilSTGB_near.push_back(k);
					}
				}
			}
		}
	}

	
	cout << "We can rationalize this GB by taking a rotation axis corresponding to the normal of the (" << hp_near << "_" << kp_near << "_" << lp_near << ") plane (instead of direction [" << h_a_func << "_" << k_a_func << "_" << l_a_func << "]) and a rotation angle of " << 180.*(theta_func+true_DeltaTheta[ind_rattheta])/M_PI << " (instead of " << theta_func*180./M_PI << ")" << endl;
	if( MT->max_p(uint_vec,3) < 100 ){
		if( MilSTGB_true.size() != 0 ){
			unsigned int ind_STGB = MT->min(normSTGB_true);
			cout << "For the rotation axis (" << hp_near << kp_near << lp_near << ") and this angle, the corresponding symmetrical tilt GB have a (" << MilSTGB_true[ind_STGB*3] << "_" << MilSTGB_true[ind_STGB*3+1] << "_" << MilSTGB_true[ind_STGB*3+2] << ") GB plane" << endl;
		}else if( MilSTGB_near.size() != 0 ){
			unsigned int ind_STGB = MT->min(normSTGB_near);
			cout << "For the rotation axis [" << h_a_func << k_a_func << l_a_func << "] (not the one rationalized) and this angle, the corresponding symmetrical tilt GB have a (" << MilSTGB_near[ind_STGB*3] << "_" << MilSTGB_near[ind_STGB*3+1] << "_" << MilSTGB_near[ind_STGB*3+2] << ") GB plane" << endl;
		}else{
			cout << "We don't have found a corresponding STGB !" << endl;
		}
	}else{
		cout << "We don't have search the corresponding STGB of this misorientation (too long to compute)" << endl;
	}

	delete[] vec1;
	delete[] vec2;
	delete[] vec3;
	delete[] vec4;
	delete[] vec5;
	delete[] int_vec;
	delete[] uint_vec;
	delete[] RotMat;
	delete[] rot_ax_func;
	return theta_func+true_DeltaTheta[ind_rattheta];
}

//// Based on the CSL Basis, we search the size of the GB plane (in x and y dir) for which there is a CSL site // TODO warning with primitive cubic cell
//void Bicrystal::searchGBSize(const int h_p_func, const int k_p_func, const int l_p_func){
//	if( !this->IsCrystalDefined ){
//		cerr << "The crystal is not defined, we cannot compute the GB size" << endl;
//		exit(EXIT_FAILURE);
//	}
//	if( !IsCSL_Basis ){
//		cout << "The CSL Basis have not been computed, need to perform the calculation !" << endl;
//		//searchCSL(this->h_a,this->k_a,this->l_a,this->theta,CSL_vec,0);
//	}
//	// search two orthogonal vectors of the CSL belonging to the GB plane
//	double tolScalarProd = 1e-3;
//	int MaxCL_CSL = 100;
//	double *vec1 = new double[3];
//	double *vec2 = new double[3];
//	double sp;
//	// get the normal direction to the GB plane
//	for(unsigned int i=0;i<3;i++) vec1[i] = h_p_func*_MyCrystal->getA1_star()[0] + k_p_func*_MyCrystal->getA2_star()[i] + l_p_func*_MyCrystal->getA3_star()[i];
//	// search a first normal vector
//	for(int i=-MaxCL_CSL;i<MaxCL_CSL+1;i++){
//		for(int j=-MaxCL_CSL;j<MaxCL_CSL+1;j++){
//			for(int k=-MaxCL_CSL;k<MaxCL_CSL+1;k++){
//				for(unsigned int ii=0;ii<3;ii++) vec2[ii] = i*this->CSL_Basis[ii*3] + j*this->CSL_Basis[ii*3+1] + k*this->CSL_Basis[ii*3+2];
//				sp = 0;
//				for(unsigned int ii=0;ii<3;ii++) sp += vec1[ii]*vec2[ii];
//				if( fabs(sp) < tolScalarProd && ( i != 0 || j != 0 || k != 0 ) ){
//
//				}
//			}
//		}
//	}
//	delete[] vec1;
//	delete[] vec2;
//	// TODO continue this function
//}

void BicrystalloTools::read_params(){
	//string fp;
	//#ifdef FIXEDPARAMETERS
	//fp = FIXEDPARAMETERS;
	//#endif
	//string backslash="/";
	//string filename=fp+backslash+FixedParam_Filename;
	//ifstream file(filename, ios::in);
	ifstream file(FixedParam_Filename, ios::in);
	size_t pos_thetamax, pos_MaxHKL, pos_toldist, pos_sigmamax;
	string buffer_s, line;
	if(file){
		while(file){
			getline(file,line);
			pos_thetamax=line.find("THETA_MAX_ROT_ANGLE_RAT ");
			if(pos_thetamax!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->theta_max_rot_ax_rat;
			}
			pos_MaxHKL=line.find("MAX_HKL_ROT_ANGLE ");
			if(pos_MaxHKL!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->MaxHKL_rot_angle_rat;
			}
			pos_toldist=line.find("TOL_DIST_RAT_ANGLE ");
			if(pos_toldist!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->tol_dist_rot_angle_rat;
			}
			pos_sigmamax=line.find("SIGMA_MAX ");
			if(pos_sigmamax!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->SigmaMax;
			}
		}
		file.close();
	}//else{
	//	cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
	//	exit(EXIT_FAILURE);
	//}
}

void BicrystalloTools::ReadProperties(vector<string> Properties){
	size_t pos_thetamax, pos_MaxHKL, pos_toldist, pos_sigmamax;
	string buffer_s, line;
	for(unsigned int i=0;i<Properties.size();i++){
		pos_thetamax=Properties[i].find("THETA_MAX_ROT_ANGLE_RAT ");
		if(pos_thetamax!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->theta_max_rot_ax_rat;
		}
		pos_MaxHKL=Properties[i].find("MAX_HKL_ROT_ANGLE ");
		if(pos_MaxHKL!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->MaxHKL_rot_angle_rat;
		}
		pos_toldist=Properties[i].find("TOL_DIST_RAT_ANGLE ");
		if(pos_toldist!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->tol_dist_rot_angle_rat;
		}
		pos_sigmamax=Properties[i].find("SIGMA_MAX ");
		if(pos_sigmamax!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->SigmaMax;
		}
	}
}

void BicrystalloTools::SaveLattices(double *a1_s, double *a2_s){
	for(unsigned int i=0;i<9;i++){
		current_a1[i] = a1_s[i];
		current_a2[i] = a2_s[i];
	}
}

void BicrystalloTools::SaveCSL_DSC(double *csl, double *dsc){
	for(unsigned int i=0;i<9;i++){
		current_CSL[i] = csl[i];
		current_DSC[i] = dsc[i];
	}
}

void BicrystalloTools::PrintCurrentCSL_DSC(){
	vector<vector<unsigned int>> arr_fusion = {{1,1,1,1},{1,1,1,1},{1,1,1,1},{1,1,1,1}};
	cout << endl << "The current lattices have a Sigma = " << current_sigma << endl;
	cout << "The CSL and DSC basis in the Cartesian frame are :" << endl;
	vector<vector<string>> arr_element = {
		{"Basis vector","X","Y","Z"},
		{"CSL[1]",to_string(current_CSL[0]),to_string(current_CSL[3]),to_string(current_CSL[6])},
		{"CSL[2]",to_string(current_CSL[1]),to_string(current_CSL[4]),to_string(current_CSL[7])},
		{"CSL[3]",to_string(current_CSL[2]),to_string(current_CSL[5]),to_string(current_CSL[8])}
	};
	Dis.DisplayArray(arr_element,arr_fusion);
	cout << endl;
	arr_element = {
	    {"Basis vector", "X", "Y", "Z"},
	    {"DSC[1]", to_string(this->current_DSC[0]), to_string(this->current_DSC[3]), to_string(this->current_DSC[6])},
	    {"DSC[2]", to_string(this->current_DSC[1]), to_string(this->current_DSC[4]), to_string(this->current_DSC[7])},
	    {"DSC[3]", to_string(this->current_DSC[2]), to_string(this->current_DSC[5]), to_string(this->current_DSC[8])}
	};
	Dis.DisplayArray(arr_element, arr_fusion);
	cout << endl;
}

int BicrystalloTools::EndComputeCSL_DSC(int sigma_csl, int sigma_dsc, double *CSL_Basis, double *DSC_Basis){
	if( abs(sigma_csl-sigma_dsc) > 1 ){
		cout << "Warning, the sigma values obtained through CSL and DSC research are different (" << sigma_csl << " and " << sigma_dsc << "), try to rationalize this GB to fix the issue" << endl;
		current_sigma = (sigma_csl+sigma_dsc)/2;
	}else if( sigma_csl != sigma_dsc ){
		if( sigma_csl % 2 == 0 ) current_sigma = sigma_dsc;
		else current_sigma = sigma_csl;
	}else if( sigma_csl % 2 == 0 ) cout << "Warning the sigma value is even, which could reflect something weird" << endl;
	SaveCSL_DSC(CSL_Basis,DSC_Basis);
	PrintCurrentCSL_DSC();
	return current_sigma;
}

bool BicrystalloTools::Check_CurrentBasis(){
	double V_csl = fabs(MT->det(current_CSL));
	double V_dsc = fabs(MT->det(current_DSC));
	double V_lat1 = fabs(MT->det(current_a1));
	double tol_vol = 1.e-5;
	bool cons = true;
	double diff = fabs(V_csl-(current_sigma*V_lat1));
	if( diff > tol_vol ){
		cout << "Warning V_CSL (" << to_string(V_csl) << ") != Sigma * V_lattice ( " << to_string(current_sigma) << " * " << to_string(V_lat1) << " = " << to_string(current_sigma*V_lat1) << ") which may means that either the computed CSL is not primitive or something wrong happened" << endl;
		cons = false;
	}
	diff = fabs((V_dsc*current_sigma)-V_lat1);
	if( diff > tol_vol ){
		cout << "Warning V_DSC (" << to_string(V_dsc) << ") != V_lattice / Sigma ( " << to_string(V_lat1) << " / " << to_string(current_sigma) << " = " << to_string(V_lat1/current_sigma) << ") which may means that either the computed DSC is not primitive or something wrong happened" << endl;
		cons = false;
	}
	return cons;
}

BicrystalloTools::~BicrystalloTools(){
	delete[] current_a1;
	delete[] current_a2;
	delete[] current_CSL;
	delete[] current_DSC;
	delete MT;
	if( MyCSL_Crystal ) delete MyCSL_Crystal;
}
