//**********************************************************************************
//*   Crystal.cpp                                                                  *
//**********************************************************************************
//* This file contains the implementation of the Crystal class                     *
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


#include "Crystal.h"
#include "MyStructs.h"
#include "AtomHicConfig.h"
#include <iostream>
#include <dirent.h>
#include <vector>
#include <fstream>
#include <sstream>
#include "MathTools.h"
#include <cmath>

using namespace std;

Crystal::Crystal(const string& crystalName){
	read_params();
	this->a1 = new double[3];
	this->a2 = new double[3];
	this->a3 = new double[3];
	this->a1_star = new double[3];
	this->a2_star = new double[3];
	this->a3_star = new double[3];
	this->MT = new MathTools;
	this->IsCharge = false;
	// initialize these matrix as identity
	this->rot_mat_total = new double[9];
	this->rot_mat_total[0] = 1.;
	this->rot_mat_total[1] = 0.;
	this->rot_mat_total[2] = 0.;
	this->rot_mat_total[3] = 0.;
	this->rot_mat_total[4] = 1.;
	this->rot_mat_total[5] = 0.;
	this->rot_mat_total[6] = 0.;
	this->rot_mat_total[7] = 0.;
	this->rot_mat_total[8] = 1.;
	this->TiltTrans_xyz = new double[9];
	this->TiltTrans_xyz[0] = 1.;
	this->TiltTrans_xyz[1] = 0.;
	this->TiltTrans_xyz[2] = 0.;
	this->TiltTrans_xyz[3] = 0.;
	this->TiltTrans_xyz[4] = 1.;
	this->TiltTrans_xyz[5] = 0.;
	this->TiltTrans_xyz[6] = 0.;
	this->TiltTrans_xyz[7] = 0.;
	this->TiltTrans_xyz[8] = 1.;
	this->path2database = getDatabasePath(crystalName);
	read_database();
	this->alength = new double[3];
	this->alength[0] = sqrt(pow(this->a1[0],2.)+pow(this->a1[1],2.)+pow(this->a1[2],2.));
	this->alength[1] = sqrt(pow(this->a2[0],2.)+pow(this->a2[1],2.)+pow(this->a2[2],2.));
	this->alength[2] = sqrt(pow(this->a3[0],2.)+pow(this->a3[1],2.)+pow(this->a3[2],2.));
	if( this->IsDoNotSep == true ) ConstructNotSepList();
	computeReciproqual();
	computeStoich();
	for(unsigned int i=0;i<this->nbAtomType;i++){
		if( this->NbAtomSite[i] > 1 ){
			this->IsMultisite = true;
			break;
		}
	}
}

string Crystal::getDatabasePath(string crystalName){
	// Read the crystal database
	char *database_env = getenv("CRYSTAL_DATABASE");
	string database, returned_path;
	unsigned int crystal_index;
	if (database_env) {
		database = database_env;
	} else {
		#ifdef CRYSTAL_DATABASE
		database = CRYSTAL_DATABASE;
		#endif
	}
	if( database.empty() ){
		cerr << "Warning database environment for crystal is empty" << endl;
		exit(EXIT_FAILURE);
	} else {
		DIR *dir;
		struct dirent *diread;
		vector<string> AvailableCrystals;
		const char *env = database.c_str();
		string buffer_s;
		size_t pos;
		if( (dir = opendir(env) ) != nullptr ){
			while( (diread = readdir(dir)) != nullptr ){
				buffer_s = diread->d_name;
				pos = buffer_s.find(".ath");
				if(pos!=string::npos) AvailableCrystals.push_back(buffer_s.erase(buffer_s.size()-4));
			}
			closedir(dir);
		}else{
			perror("opendir");
		}

		// search into the database if crystalName is present
		bool crystalok = false;
		for(unsigned int i=0;i<AvailableCrystals.size();i++){
			if( crystalName == AvailableCrystals[i] ){
				crystal_index = i;
				this->name = crystalName; 
				crystalok = true;
				break;
			}
		}
		if( !crystalok ){
			cerr << "The crystal \"" << crystalName << "\" does not exist in the crystal database, please create \"" << crystalName << ".ath\" in the /data/Crystal/ directory following the example in /data/ExampleFiles/Crystal.ath or use an other Crystal constructor" << endl;
			cerr << "List of available crystals :" << endl;
			for(unsigned int i=0;i<AvailableCrystals.size();i++) cerr << AvailableCrystals[i] << endl;
			exit(EXIT_FAILURE);
		}else{
			returned_path = database+"/"+AvailableCrystals[crystal_index]+".ath";
		}
	}
	return returned_path;

}

// construct a z oriented (hkl) unit cell plane
void Crystal::RotateCrystal(const int& h_p, const int& k_p, const int& l_p){
	int arr[3] = {abs(h_p),abs(k_p),abs(l_p)};
	int maxhkl = MT->max(arr,3)*15;
	int max_hkl_cryst = 100;
	if( maxhkl > max_hkl_cryst ) maxhkl = max_hkl_cryst;
	double tolScalarProd = 1e-2;
	vector<double> buffer_vector_d;
	vector<int> buffer_vector_i;
	double *normalDir = new double[3];
	double *rot_mat = new double[9];
	double *rot_axis = new double[3];
	double *test_vec = new double[3];
	double *orthoDir = new double[3];
	double xbox,ybox,zbox; // box dimension
	cout << maxhkl << endl;
	// direction normal to the wanted plane
	normalDir[0] = h_p*(this->a1_star[0]) + k_p*(this->a2_star[0]) + l_p*(this->a3_star[0]);	
	normalDir[1] = h_p*(this->a1_star[1]) + k_p*(this->a2_star[1]) + l_p*(this->a3_star[1]);	
	normalDir[2] = h_p*(this->a1_star[2]) + k_p*(this->a2_star[2]) + l_p*(this->a3_star[2]);
	// search the smallest orthogonal direction with integer Miller indices to the wanted one => to be aligned with the cartesian x axis
	for(int i=-maxhkl;i<maxhkl;i++){
		for(int j=-maxhkl;j<maxhkl;j++){
			for(int k=-maxhkl;k<maxhkl;k++){
				if( ( abs(i*h_p+j*k_p+k*l_p) < tolScalarProd ) && ( i != 0 || j != 0 || k!= 0 ) ){
					buffer_vector_d.push_back(pow(i,2.)+pow(j,2.)+pow(k,2.));
					buffer_vector_i.push_back(i);
					buffer_vector_i.push_back(j);
					buffer_vector_i.push_back(k);
				}
			}
		}
	}
	unsigned int minhkl = MT->min(buffer_vector_d);
	orthoDir[0] = (buffer_vector_i[minhkl*3]*(this->a1[0])) + (buffer_vector_i[minhkl*3+1]*(this->a2[0])) + (buffer_vector_i[minhkl*3+2]*(this->a3[0]));
	orthoDir[1] = (buffer_vector_i[minhkl*3]*(this->a1[1])) + (buffer_vector_i[minhkl*3+1]*(this->a2[1])) + (buffer_vector_i[minhkl*3+2]*(this->a3[1]));
	orthoDir[2] = (buffer_vector_i[minhkl*3]*(this->a1[2])) + (buffer_vector_i[minhkl*3+1]*(this->a2[2])) + (buffer_vector_i[minhkl*3+2]*(this->a3[2]));
	// compute the rotation angles needed to align the wanted direction with the cartesian z axis and have the x cartesian axis corresponding to the smallest crystallographic direction possible
	// for rotation around the z axis to align the normal dir with x axis
	double theta_z, theta_y;
	if( (h_p == 0) && (k_p == 0) ) theta_z = 0.;
 	else if( k_p >= 0 ) theta_z = -acos(normalDir[0]/sqrt(pow(normalDir[0],2.)+pow(normalDir[1],2.)));
	else if( k_p < 0 ) theta_z = acos(normalDir[0]/sqrt(pow(normalDir[0],2.)+pow(normalDir[1],2.))); 
	rot_axis[0] = 0;
	rot_axis[1] = 0;
	rot_axis[2] = 1;
       	MT->Vec2rotMat(rot_axis,theta_z,this->rot_mat_total);
	for(unsigned int i=0;i<3;i++) test_vec[i] = normalDir[i];
	MT->MatDotRawVec(this->rot_mat_total,test_vec,test_vec);
	// second rotation around y to align the normal dir with the z axis // I have remplaced normalDir by test_vec in the two following lines
	if( test_vec[0] >= 0 ) theta_y = -acos(test_vec[2]/sqrt(pow(test_vec[0],2.)+pow(test_vec[1],2.)+pow(test_vec[2],2.)));	
	else theta_y = acos(test_vec[2]/sqrt(pow(test_vec[0],2.)+pow(test_vec[1],2.)+pow(test_vec[2],2.)));	
	rot_axis[0] = 0;
	rot_axis[1] = 1;
	rot_axis[2] = 0;
	MT->Vec2rotMat(rot_axis,theta_y,rot_mat);
	MT->MatDotMat(rot_mat,this->rot_mat_total,this->rot_mat_total);
	// last rotation about z axis to align the smallest orthogonal vector with the x direction
	MT->MatDotRawVec(this->rot_mat_total,orthoDir,orthoDir);
	if( orthoDir[1] <= 0. ) theta_z = acos(orthoDir[0]/sqrt(pow(orthoDir[0],2.)+pow(orthoDir[1],2.)+pow(orthoDir[2],2.)));
	else theta_z = -acos(orthoDir[0]/sqrt(pow(orthoDir[0],2.)+pow(orthoDir[1],2.)+pow(orthoDir[2],2.)));
	rot_axis[0] = 0;
	rot_axis[1] = 0;
	rot_axis[2] = 1;
       	MT->Vec2rotMat(rot_axis,theta_z,rot_mat);
	// last rotation for the orthoDir vector
	MT->MatDotRawVec(rot_mat,orthoDir,orthoDir);
	// compute the total rotation matrix
	MT->MatDotMat(rot_mat,this->rot_mat_total,this->rot_mat_total);
	// full rotation of the normal dir
	MT->MatDotRawVec(this->rot_mat_total,normalDir,normalDir);
	// verify is the rotation has been well achieved
	if( fabs(normalDir[0]) > tolScalarProd || fabs(normalDir[1]) > tolScalarProd || fabs(orthoDir[1]) > tolScalarProd || fabs(orthoDir[2]) > tolScalarProd ){
		cout << "The crystal basis cannot be aligned with the cartesian one, aborting calculation" << endl;
		return;
	}else cout << "We are constructing a system with the (" << h_p << k_p << l_p << ") plane parallel with (xy) plane and the [" << buffer_vector_i[minhkl*3] << buffer_vector_i[minhkl*3+1] << buffer_vector_i[minhkl*3+2] << "] direction aligned with the x axis" << endl;

	// rotate lattice vectors
	MT->MatDotRawVec(this->rot_mat_total,this->a1,this->a1);
	MT->MatDotRawVec(this->rot_mat_total,this->a2,this->a2);
	MT->MatDotRawVec(this->rot_mat_total,this->a3,this->a3);
	// recompute reciproqual space as cell parameters may have changed
	computeReciproqual();

	delete[] orthoDir;
	delete[] rot_mat;
	delete[] rot_axis;
	delete[] test_vec;
	delete[] normalDir;
}

void Crystal::RotateCrystal(const double *RotMat){
	double xbox, ybox, zbox;
	for(unsigned int i=0;i<9;i++) this->rot_mat_total[i] = RotMat[i];
	// rotate lattice vectors
	MT->MatDotRawVec(RotMat,this->a1,this->a1);
	MT->MatDotRawVec(RotMat,this->a2,this->a2);
	MT->MatDotRawVec(RotMat,this->a3,this->a3);
	// recompute reciproqual space as cell parameters may have changed
	computeReciproqual();
}


void Crystal::ConstructNotSepList(){
	// construct full neighbor list in a large cutoff radius
	double arr[3] = {this->alength[0], this->alength[1], this->alength[2]};
	double rc_squared = pow(this->MT->max(arr,3),2.);
	double d_squared, xpos, ypos, zpos;
	int cl = 4;
	vector<vector<double>> buffer_vector_d;
	for(unsigned int i=0;i<this->nbAtom;i++){
		this->NotSepList.push_back(vector<int> ());
		buffer_vector_d.push_back(vector<double> ());
		for(int xcl=-cl;xcl<cl+1;xcl++){
			for(int ycl=-cl;ycl<cl+1;ycl++){
				for(int zcl=-cl;zcl<cl+1;zcl++){
					for(unsigned int j=0;j<this->nbAtom;j++){
						xpos = this->Motif[j].pos.x + this->a1[0]*xcl+this->a2[0]*ycl+this->a3[0]*zcl;
						ypos = this->Motif[j].pos.y + this->a1[1]*xcl+this->a2[1]*ycl+this->a3[1]*zcl;
						zpos = this->Motif[j].pos.z + this->a1[2]*xcl+this->a2[2]*ycl+this->a3[2]*zcl;
						d_squared = pow(xpos-this->Motif[i].pos.x,2.) + pow(ypos-this->Motif[i].pos.y,2.) + pow(zpos-this->Motif[i].pos.z,2.); 
						if( d_squared < rc_squared && d_squared > 1e-6 ){
							buffer_vector_d[i].push_back(d_squared);
							buffer_vector_d[i].push_back(j);
							buffer_vector_d[i].push_back(xcl);
							buffer_vector_d[i].push_back(ycl);
							buffer_vector_d[i].push_back(zcl);
						}
					}
				}
			}
		}
	}
	unsigned int count;
	bool Continue;
	for(unsigned int i=0;i<this->nbAtom;i++){
		this->MT->sort(buffer_vector_d[i], 0, 5, buffer_vector_d[i]);
		for(unsigned int j=0;j<this->DoNotSep.size();j++){
			if( DoNotSep[j][0] == this->Motif[i].type_uint ){
				count = 0;
				for(unsigned int n=0;n<buffer_vector_d[i].size()/5;n++){
					if( count == this->DoNotSep[j][1] ) break;
					else if( this->Motif[(int) buffer_vector_d[i][n*5+1]].type_uint == DoNotSep[j][2] ){
						// verify that this atom has been stored before
						Continue = false;
						for(unsigned int s_1=0;s_1<NotSepList.size();s_1++){
							for(unsigned int s_2=0;s_2<NotSepList[s_1].size()/4;s_2++){
								// Case where all the motif is not necessarily get
								if( round(buffer_vector_d[i][n*5+1]) == NotSepList[s_1][s_2*4] && round(buffer_vector_d[i][n*5+2]) == NotSepList[s_1][s_2*4+1] && round(buffer_vector_d[i][n*5+3]) == NotSepList[s_1][s_2*4+2] && round(buffer_vector_d[i][n*5+4]) == NotSepList[s_1][s_2*4+3] ){
								// Case to have the all motif
								//if( round(buffer_vector_d[i][n*5+1]) == NotSepList[s_1][s_2*4] ){
									Continue = true;
									break;
								}
							}
							if( Continue ) break;
						}
						if( Continue ) continue;
						count += 1;
						NotSepList[i].push_back(round(buffer_vector_d[i][n*5+1]));
						NotSepList[i].push_back(round(buffer_vector_d[i][n*5+2]));
						NotSepList[i].push_back(round(buffer_vector_d[i][n*5+3]));
						NotSepList[i].push_back(round(buffer_vector_d[i][n*5+4]));
					}
				}
			}
		}
	}
}

void Crystal::ConstructOrthogonalCell(){
	int CLsearch = 150;
	// search the smallest orthogonal box for this orientation by testing linear combination of a1 a2 a3 giving cell vectors aligned with cartesian axis
	// expect for the x axis which is already aligned with a crystallographic direction
	vector<int> cl_box;
	vector<int> xh_list;
	vector<int> yh_list;
	vector<int> zh_list;
	vector<double> buffer_vector_d_x, buffer_vector_d_y, buffer_vector_d_z;
	double absX, absY, absZ;
	for(int i=-CLsearch;i<CLsearch+1;i++){
		for(int j=-CLsearch;j<CLsearch+1;j++){
			for(int k=-CLsearch;k<CLsearch+1;k++){
				absX = fabs(i*this->a1[0]+j*this->a2[0]+k*this->a3[0]);
				absY = fabs(i*this->a1[1]+j*this->a2[1]+k*this->a3[1]);
				absZ = fabs(i*this->a1[2]+j*this->a2[2]+k*this->a3[2]);
				if( absY < TolOrthoBox && absZ < TolOrthoBox && i*this->a1[0]+j*this->a2[0]+k*this->a3[0] > MinBoxAside ){
					buffer_vector_d_x.push_back(absX);
					xh_list.push_back(i);
					xh_list.push_back(j);
					xh_list.push_back(k);
				}
				if( absX < TolOrthoBoxZ && absY < TolOrthoBoxZ && i*this->a1[2]+j*this->a2[2]+k*this->a3[2] > MinBoxHeight ){
					buffer_vector_d_z.push_back(absZ);
					zh_list.push_back(i);
					zh_list.push_back(j);
					zh_list.push_back(k);
				}
				if( absX < TolOrthoBox && absZ < TolOrthoBox && i*this->a1[1]+j*this->a2[1]+k*this->a3[1] > MinBoxAside ){
					buffer_vector_d_y.push_back(absY);
					yh_list.push_back(i);
					yh_list.push_back(j);
					yh_list.push_back(k);
				}
			}
		}
	}
	if( xh_list.size() == 0 || yh_list.size() == 0 || zh_list.size() == 0 ){
		cerr << "There is no linear combination of crystal vectors giving an orthogonal cell, consider increase the tolerance or the number of CL used for the research, aborting computation" << endl;
		exit(EXIT_FAILURE);
	}

	unsigned int ind_x, ind_y, ind_z;
	ind_x = MT->min(buffer_vector_d_x);
	ind_y = MT->min(buffer_vector_d_y);
	ind_z = MT->min(buffer_vector_d_z);
	// fill the cl_box array
	for(unsigned int i=0;i<3;i++) cl_box.push_back(xh_list[ind_x*3+i]);
	for(unsigned int i=0;i<3;i++) cl_box.push_back(yh_list[ind_x*3+i]);
	for(unsigned int i=0;i<3;i++) cl_box.push_back(zh_list[ind_x*3+i]);
	double *xh = new double[3];
	double *yh = new double[3];
	double *zh = new double[3];
	for(unsigned int i=0;i<3;i++){
		xh[i] = this->a1[i]*xh_list[ind_x*3]+this->a2[i]*xh_list[ind_x*3+1]+this->a3[i]*xh_list[ind_x*3+2];
		yh[i] = this->a1[i]*yh_list[ind_y*3]+this->a2[i]*yh_list[ind_y*3+1]+this->a3[i]*yh_list[ind_y*3+2];
		zh[i] = this->a1[i]*zh_list[ind_z*3]+this->a2[i]*zh_list[ind_z*3+1]+this->a3[i]*zh_list[ind_z*3+2];
	}

	double xbox = sqrt( pow(xh[0],2.) + pow(xh[1],2.) + pow(xh[2],2.));
	double ybox = sqrt( pow(yh[0],2.) + pow(yh[1],2.) + pow(yh[2],2.));
	double zbox = sqrt( pow(zh[0],2.) + pow(zh[1],2.) + pow(zh[2],2.));

	MT->computeTiltTrans(xh,yh,zh,this->TiltTrans_xyz);

	// apply this transformation to the crystal vectors
	MT->MatDotRawVec(TiltTrans_xyz,this->a1,this->a1);
	MT->MatDotRawVec(TiltTrans_xyz,this->a2,this->a2);
	MT->MatDotRawVec(TiltTrans_xyz,this->a3,this->a3);
	// verify that the system has now an orthogonal cell
	for(unsigned int i=0;i<3;i++){
		xh[i] = this->a1[i]*xh_list[ind_x*3]+this->a2[i]*xh_list[ind_x*3+1]+this->a3[i]*xh_list[ind_x*3+2];
		yh[i] = this->a1[i]*yh_list[ind_y*3]+this->a2[i]*yh_list[ind_y*3+1]+this->a3[i]*yh_list[ind_y*3+2];
		zh[i] = this->a1[i]*zh_list[ind_z*3]+this->a2[i]*zh_list[ind_z*3+1]+this->a3[i]*zh_list[ind_z*3+2];
	}
	double tolBase = 1e-10;
	if( fabs(xh[1]) > tolBase || fabs(xh[2]) > tolBase || fabs(yh[0]) > tolBase || fabs(yh[2]) > tolBase || fabs(zh[0]) > tolBase || fabs(zh[1]) > tolBase || fabs(xh[0]-xbox) > tolBase || fabs(yh[1]-ybox) > tolBase || fabs(zh[2]-zbox) > tolBase ){
		cout << "Warning there were issues during the construction of orthogonal cell" << endl;
	}
	// recompute reciproqual space as cell parameters may have changed
	computeReciproqual();
	// compute the total transformation matrix
	double *TotalTrans = new double[9];
	MT->MatDotMat(TiltTrans_xyz,this->rot_mat_total,TotalTrans);
	// shift the motif
	for(unsigned int i=0;i<this->nbAtom;i++){
		Motif[i].pos.x += shift_x;
		Motif[i].pos.y += shift_y;
		Motif[i].pos.z += shift_z;
	}
	// rotate the motif
	for(unsigned int i=0;i<this->nbAtom;i++) MT->MatDotAt(TotalTrans,Motif[i],Motif[i]);	
	
	// construct the AtomicSystem
	this->OrientedSystem = new AtomicSystem(this, xbox, ybox, zbox, cl_box);
	// rescale the crystal vectors as they have been modified by the function
	double *invTrans = new double[9];
	MT->invert3x3(this->TiltTrans_xyz,invTrans);
	//MT->MatDotRawVec(invTrans,this->a1,this->a1);
	//MT->MatDotRawVec(invTrans,this->a2,this->a2);
	//MT->MatDotRawVec(invTrans,this->a3,this->a3);
	// recompute reciproqual space as cell parameters may have changed
	computeReciproqual();

	this->IsOrientedSystem = true;
	delete[] invTrans;
	delete[] TotalTrans;
	delete[] xh;
	delete[] yh;
	delete[] zh;
}

void Crystal::computeReciproqual(){
	// compute cell volume
	double *mixedProd = new double[3];
	MT->mixedProd(this->a2,this->a3,mixedProd);
	this->V = MT->dotProd(this->a1,mixedProd);
	// compute reciproqual vectors
	MT->mixedProd(this->a2,this->a3,this->a1_star);
	MT->mixedProd(this->a3,this->a1,this->a2_star);
	MT->mixedProd(this->a1,this->a2,this->a3_star);
	for(unsigned int i=0;i<3;i++){
		this->a1_star[i] /= this->V;
		this->a2_star[i] /= this->V;
		this->a3_star[i] /= this->V;
	}
	delete[] mixedProd;
}

void Crystal::read_database(){
	ifstream file(this->path2database, ios::in);
	size_t pos_at, pos_x, pos_y, pos_z, pos_attype, pos_Mass, pos_At, pos_Crystal, pos_tilt, pos_DNS, pos_bondori;
	unsigned int line_Mass(1000), line_At(1000), buffer_uint, buffer_uint_1, buffer_uint_2, count(0), line_bondori(1000), nbref2read(0);
	double buffer_1, buffer_2, buffer_3, buffer_4;
	string buffer_s, buffer_s_1, buffer_s_2, line;
	if(file){
		while(file){
			getline(file,line);
			// find crystallogaphy
			pos_Crystal=line.find("CRYSTAL");
			if(pos_Crystal!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1;
				this->crystallo = buffer_s_1;
			}
			// find number of atom
			pos_at=line.find("atoms");
			if(pos_at!=string::npos){
				istringstream text(line);
				text >> buffer_uint;
				this->nbAtom = buffer_uint;
				Motif = new Atom[this->nbAtom];
				this->AtomSite = new unsigned int[this->nbAtom];
			}

			// find H1 vector
			pos_x=line.find("xlo xhi");
			if(pos_x!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				this->a1[0] = buffer_2-buffer_1;
			}

			// find H2 vector
			pos_y=line.find("ylo yhi");
			if(pos_y!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				this->a2[1] = buffer_2-buffer_1;
			}

			// find H3 vector
			pos_z=line.find("zlo zhi");
			if(pos_z!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				this->a3[2] = buffer_2-buffer_1;
			}
			// search tilt component for non cubic crystal
			if( this->crystallo == "Hexagonal"){
				pos_tilt=line.find("xy xz yz");
				if(pos_tilt!=string::npos){
					istringstream text(line);
					text >> buffer_1 >> buffer_2 >> buffer_3;
					this->a2[0] = buffer_1;
				}
			}

			// find nb atom type
			pos_attype=line.find("atom types");
			if(pos_attype!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_s >> buffer_s_1;
				this->nbAtomType = buffer_1;
				AtomType = new string[this->MaxAtomType];
				AtomType_uint = new unsigned int[this->MaxAtomType];
				NbAtomSite = new unsigned int[this->MaxAtomType];
				AtomMass = new double[this->MaxAtomType];
				AtomCharge = new double[this->MaxAtomType];
			}

			// read DoNotSep array
			pos_DNS = line.find("DONOTSEPARE");
			if(pos_DNS!=string::npos){
				this->IsDoNotSep = true;
				istringstream text(line);
				text >> buffer_s >> buffer_uint >> buffer_uint_1 >> buffer_uint_2;
				this->DoNotSep.push_back(vector<unsigned int> ());
				this->DoNotSep[this->DoNotSep.size()-1].push_back(buffer_uint);
				this->DoNotSep[this->DoNotSep.size()-1].push_back(buffer_uint_1);
				this->DoNotSep[this->DoNotSep.size()-1].push_back(buffer_uint_2);
			}

			// get lines where are the keywords Masses and Atoms to get atom type masses and positions
			pos_Mass=line.find("Masses");
			if(pos_Mass!=string::npos) line_Mass = count;
			if( count > line_Mass+1 && count < line_Mass+2+this->nbAtomType ){
				istringstream text(line);
				text >> buffer_uint >> buffer_1 >> buffer_s_1 >> buffer_s >> buffer_uint_1;
				this->AtomMass[buffer_uint-1] = buffer_1;
				this->AtomType[buffer_uint-1] = buffer_s;
				this->AtomType_uint[buffer_uint-1] = buffer_uint;
				this->NbAtomSite[buffer_uint-1] = buffer_uint_1;
			}
			pos_At=line.find("Atoms");
			if(pos_At!=string::npos){
				istringstream text(line);
				text >> buffer_s_1 >> buffer_s_2 >> buffer_s;
				if( buffer_s == "charge" ) this->IsCharge = true;
			       	line_At = count;
			}
			if( count > line_At+1 && count < line_At+nbAtom+2 ){
				istringstream text(line);
				if( this->IsCharge ){
					text >> buffer_uint >> buffer_uint_1 >> buffer_1 >> buffer_2 >> buffer_3 >> buffer_4 >> buffer_uint_2;
					this->Motif[buffer_uint-1].pos.x = buffer_2;
					this->Motif[buffer_uint-1].pos.y = buffer_3;
					this->Motif[buffer_uint-1].pos.z = buffer_4;
					this->Motif[buffer_uint-1].type_uint = buffer_uint_1;
					this->AtomSite[buffer_uint-1] = buffer_uint_2-1;
					this->AtomCharge[buffer_uint_1-1] = buffer_1;
				}else{
					text >> buffer_uint >> buffer_uint_1 >> buffer_2 >> buffer_3 >> buffer_4 >> buffer_uint_2;
					this->Motif[buffer_uint-1].pos.x = buffer_2;
					this->Motif[buffer_uint-1].pos.y = buffer_3;
					this->Motif[buffer_uint-1].pos.z = buffer_4;
					this->Motif[buffer_uint-1].type_uint = buffer_uint_1;
					this->AtomSite[buffer_uint-1] = buffer_uint_2-1;
				}
			}
			pos_bondori=line.find("REFERENCE_BOND_ORIENTATIONAL_PARAMETERS");
			if(pos_bondori!=string::npos){
				line_bondori = count;
				IsReferenceBondOriParam = true;
				ReferenceBondOriParam = new vector<double>[MaxAtomType];
				for(unsigned int t=0;t<nbAtomType;t++) nbref2read += NbAtomSite[t];
			}
			if( count > line_bondori && count < line_bondori+4 ) BondOriParamProperties.push_back(line);
			if( count > line_bondori+3 && count < line_bondori+4+nbref2read ){
				istringstream text(line);
				unsigned int type, site;
				double bondori;
				text >> type >> site >> bondori;
				//if( type != ReferenceBondOriParam.size() ) ReferenceBondOriParam.push_back(vector<double>());
				ReferenceBondOriParam[type-1].push_back(bondori);
			}
			count += 1;
		}
		if( !this->IsCharge ) for(unsigned int i=0;i<this->nbAtomType;i++) AtomCharge[i] = 0.;
		file.close();
		this->a1[1] = 0;
		this->a1[2] = 0;
		this->a2[2] = 0;
		this->a3[0] = 0;
		this->a3[1] = 0;
		if( this->crystallo == "Cubic" || this->crystallo == "Orthorhombic" ){
			this->a2[0] = 0;
		}
	}else{
		cout << "The file " << this->path2database << " cannot be openned" << endl;
	}
}

void Crystal::computeStoich(){
	this->Stoichiometry = new unsigned int[this->MaxAtomType];
	for(unsigned int i=0;i<this->MaxAtomType;i++) this->Stoichiometry[i] = 0;
	for(unsigned int i=0;i<this->nbAtom;i++){
		for(unsigned int t=0;t<this->nbAtomType;t++){
			if( this->Motif[i].type_uint == this->AtomType_uint[t] ){
				this->Stoichiometry[t] += 1;
				break;
			}
		}
	}
}

void Crystal::read_params(){
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+FixedParam_Filename;
	ifstream file(filename, ios::in);
	size_t pos_tolOrthoBox, pos_tolOrthoBoxZ, pos_minBoxHeight, pos_minBoxAside, pos_shift;
	string buffer_s, line;
	if(file){
		while(file){
			getline(file,line);
			pos_tolOrthoBox=line.find("TOL_ORTHO_BOX ");
			if(pos_tolOrthoBox!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->TolOrthoBox;
			}
			pos_tolOrthoBoxZ=line.find("TOL_ORTHO_BOX_Z ");
			if(pos_tolOrthoBoxZ!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->TolOrthoBoxZ;
			}
			pos_minBoxHeight=line.find("MIN_BOX_HEIGHT");
			if(pos_minBoxHeight!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->MinBoxHeight;
			}
			pos_minBoxAside=line.find("MIN_BOX_ASIDE");
			if(pos_minBoxAside!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->MinBoxAside;
			}
			pos_shift=line.find("MOTIF_SHIFT");
			if(pos_shift!=string::npos){
				istringstream text(line);
				text >> buffer_s >> shift_x >> shift_y >> shift_z;
			}

		}
	}else{
		cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
		exit(EXIT_FAILURE);
	}
}

void Crystal::ChangeTypes(unsigned int *CorresArray){
	// Variables to modify:
	// AtomType AtomType_uint NbAtomSite AtomMass AtomCharge, Motif, DoNotSep, Stoichiometry, ReferenceBondOriParam
	// Update Motif
	for(unsigned int n=0;n<nbAtom;n++){
		unsigned int old_type = Motif[n].type_uint;
		Motif[n].type_uint = CorresArray[old_type-1]+1;
	}
	// Update DoNotSep list
	for(unsigned int d=0;d<DoNotSep.size();d++){
		unsigned int old_type1 = DoNotSep[d][0];
		unsigned int old_type2 = DoNotSep[d][2];
		DoNotSep[d][0] = CorresArray[old_type1-1]+1;
		DoNotSep[d][2] = CorresArray[old_type2-1]+1;
	}
	// Copy the data
	string *AtomType_tmp = new string[MaxAtomType];	
	unsigned int *AtomType_uint_tmp = new unsigned int[MaxAtomType];
	unsigned int *NbAtomSite_tmp = new unsigned int[MaxAtomType];
	double *AtomMass_tmp = new double[MaxAtomType];
	double *AtomCharge_tmp = new double[MaxAtomType];
	vector<double> *ReferenceBondOriParam_tmp = new vector<double>[MaxAtomType];
	unsigned int *Stoichiometry_tmp = new unsigned int[this->MaxAtomType];
	for(unsigned int n=0;n<MaxAtomType;n++){
		AtomType_tmp[n] = AtomType[n];
		AtomType[n] = "";
		AtomType_uint_tmp[n] = AtomType_uint[n];
		AtomType_uint[n] = 0;
		NbAtomSite_tmp[n] = NbAtomSite[n];
		NbAtomSite[n] = 0;
		AtomMass_tmp[n] = AtomMass[n];
		AtomMass[n] = 0.;
		AtomCharge_tmp[n] = AtomCharge[n];
		AtomCharge[n] = 0.;
		Stoichiometry_tmp[n] = Stoichiometry[n];
		Stoichiometry[n] = 0;
		for(unsigned int b=0;b<ReferenceBondOriParam[n].size();b++) ReferenceBondOriParam_tmp[n].push_back(ReferenceBondOriParam[n][b]);
		ReferenceBondOriParam[n].clear();
	}
	// Update the data
	for(unsigned int n=0;n<nbAtomType;n++){
		AtomType[n] = AtomType_tmp[CorresArray[n]];
		AtomType_uint[n] = AtomType_uint_tmp[CorresArray[n]];
		NbAtomSite[n] = NbAtomSite_tmp[CorresArray[n]];
		AtomMass[n] = AtomMass_tmp[CorresArray[n]];
		AtomCharge[n] = AtomCharge_tmp[CorresArray[n]];
		Stoichiometry[n] = Stoichiometry_tmp[CorresArray[n]];
		for(unsigned int b=0;b<ReferenceBondOriParam_tmp[CorresArray[n]].size();b++) ReferenceBondOriParam[n].push_back(ReferenceBondOriParam_tmp[CorresArray[n]][b]);
	}
	delete[] AtomType_tmp;	
	delete[] AtomType_uint_tmp;
	delete[] NbAtomSite_tmp;
	delete[] AtomMass_tmp;
	delete[] AtomCharge_tmp;
	for(unsigned int t=0;t<MaxAtomType;t++) ReferenceBondOriParam_tmp[t].clear();
	delete[] ReferenceBondOriParam_tmp;
	delete[] Stoichiometry_tmp;

}

Crystal::~Crystal(){
	delete[] a1;
	delete[] a2;
	delete[] a3;
	delete[] a1_star;
	delete[] a2_star;
	delete[] a3_star;
	delete MT;
	delete[] AtomType;
	delete[] AtomType_uint;
	delete[] AtomMass;
	delete[] NbAtomSite;
	delete[] this->Motif;
	delete[] rot_mat_total;
	delete[] TiltTrans_xyz;
	delete[] alength;
	delete[] Stoichiometry;
	delete[] AtomSite;
	delete[] AtomCharge;
	if( IsOrientedSystem ) delete OrientedSystem;
	if( IsReferenceBondOriParam ){
		for(unsigned int t=0;t<MaxAtomType;t++) ReferenceBondOriParam[t].clear();
		delete[] ReferenceBondOriParam;
	}
}
