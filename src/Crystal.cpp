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
#include <unordered_set>
#include <algorithm>

using namespace std;

Crystal::Crystal(double *lattice_params, Atom *Motif_s, unsigned int nbAtom, unsigned int nbAtomType, string *AtomType_s, unsigned int *AtomType_uint_s, string name, string crystallo): nbAtom(nbAtom), nbAtomType(nbAtomType), name(name), crystallo(crystallo){
	InitVars();
	for(unsigned int i=0;i<3;i++){
		a1[i] = lattice_params[i*3];
		a2[i] = lattice_params[i*3+1];
		a3[i] = lattice_params[i*3+2];
	}
	this->alength[0] = sqrt(pow(this->a1[0],2.)+pow(this->a1[1],2.)+pow(this->a1[2],2.));
	this->alength[1] = sqrt(pow(this->a2[0],2.)+pow(this->a2[1],2.)+pow(this->a2[2],2.));
	this->alength[2] = sqrt(pow(this->a3[0],2.)+pow(this->a3[1],2.)+pow(this->a3[2],2.));
	Motif = new Atom[nbAtom];
	for(unsigned int i=0;i<nbAtom;i++) Motif[i] = Motif_s[i];
	AtomType = new string[nbAtomType];
	AtomType_uint = new unsigned int[nbAtomType];
	AtomMass = new double[nbAtomType];
	for(unsigned int i=0;i<nbAtomType;i++){
		AtomType[i] = AtomType_s[i];
		AtomType_uint[i] = AtomType_uint_s[i];
		AtomMass[i] = 0.;
	}
	computeStoich();
	computeReciproqual();
	IsMultisite = false;
	NbAtomSite = new unsigned int[nbAtomType];
	for(unsigned int i=0;i<nbAtomType;i++) NbAtomSite[i] = 1;
}

Crystal::Crystal(const string& crystalName){
	InitVars();
	this->path2database = getDatabasePath(crystalName);
	read_database();
	this->alength[0] = sqrt(pow(this->a1[0],2.)+pow(this->a1[1],2.)+pow(this->a1[2],2.));
	this->alength[1] = sqrt(pow(this->a2[0],2.)+pow(this->a2[1],2.)+pow(this->a2[2],2.));
	this->alength[2] = sqrt(pow(this->a3[0],2.)+pow(this->a3[1],2.)+pow(this->a3[2],2.));
	computeStoich();
	if( this->IsDoNotSep == true ) ConstructNotSepList();
	if( this->IsMolId == true ) ConstructNotSepListFromMolId();
	computeReciproqual();
	for(unsigned int i=0;i<this->nbAtomType;i++){
		if( this->NbAtomSite[i] > 1 ){
			this->IsMultisite = true;
			break;
		}
	}
}

void Crystal::InitVars(){
	read_params();
	crystal_def = new double[4];
	this->a1 = new double[3];
	this->a2 = new double[3];
	this->a3 = new double[3];
	this->a1_star = new double[3];
	this->a2_star = new double[3];
	this->a3_star = new double[3];
	this->OrthogonalPlanes = new int[9];
	this->OrthogonalDirs = new int[9];
	this->MT = new MathTools;
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
	this->alength = new double[3];
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
void Crystal::RotateCrystal(const int& h_p, const int& k_p, const int& l_p, const int& h_px, const int& k_px, const int& l_px){
	bool full_constrained = true;
	int h_p_s, k_p_s, l_p_s, h_px_s(0), k_px_s(0), l_px_s(0);
	if( h_px == 0 && k_px == 0 && l_px == 0 ) full_constrained = false;
	else{
		if( h_p*h_px + k_p*k_px + l_p*l_px != 0 ){
			cout << "The provided planes for orienting the crystal (i.e. (" << h_p << "_" << k_p << "_" << l_p << ") oriented along z axis and (" << h_px << "_" << k_px << "_" << l_px << ") oriented along x axis) are not normal to each others, an orthogonal cell cannot then be constructed" << endl;
			cout << "The crystal plane aligned with the x axis will thus be computed automatically" << endl;
			full_constrained = false;
		}else{
			h_px_s = h_px;
			k_px_s = k_px;
			l_px_s = l_px;
			if( IsPrimitive ) Convert2Primitive_plane(h_px_s,k_px_s,l_px_s);
		}
	}
	h_p_s = h_p;
	k_p_s = k_p;
	l_p_s = l_p;
	if( IsPrimitive ) Convert2Primitive_plane(h_p_s,k_p_s,l_p_s);
	int arr[3] = {abs(h_p_s),abs(k_p_s),abs(l_p_s)};
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
	// direction normal to the wanted plane
	normalDir[0] = h_p_s*(this->a1_star[0]) + k_p_s*(this->a2_star[0]) + l_p_s*(this->a3_star[0]);	
	normalDir[1] = h_p_s*(this->a1_star[1]) + k_p_s*(this->a2_star[1]) + l_p_s*(this->a3_star[1]);	
	normalDir[2] = h_p_s*(this->a1_star[2]) + k_p_s*(this->a2_star[2]) + l_p_s*(this->a3_star[2]);
	if( !full_constrained ){
		// search the smallest orthogonal direction with integer Miller indices to the wanted one => to be aligned with the cartesian x axis
		for(int i=-maxhkl;i<maxhkl;i++){
			for(int j=-maxhkl;j<maxhkl;j++){
				for(int k=-maxhkl;k<maxhkl;k++){
					if( ( abs(i*h_p_s+j*k_p_s+k*l_p_s) < tolScalarProd ) && ( i != 0 || j != 0 || k!= 0 ) ){
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
	}else{
		orthoDir[0] = h_px_s*(this->a1_star[0]) + k_px_s*(this->a2_star[0]) + l_px_s*(this->a3_star[0]);	
		orthoDir[1] = h_px_s*(this->a1_star[1]) + k_px_s*(this->a2_star[1]) + l_px_s*(this->a3_star[1]);	
		orthoDir[2] = h_px_s*(this->a1_star[2]) + k_px_s*(this->a2_star[2]) + l_px_s*(this->a3_star[2]);
	}
	// compute the rotation angles needed to align the wanted direction with the cartesian z axis and have the x cartesian axis corresponding to the smallest crystallographic direction possible
	// for rotation around the z axis to align the normal dir with x axis
	double theta_z, theta_y;
	if( (h_p == 0) && (k_p == 0) ) theta_z = 0.;
 	else if( normalDir[1] >= 0 ) theta_z = -acos(normalDir[0]/sqrt(pow(normalDir[0],2.)+pow(normalDir[1],2.)));
	else if( normalDir[1] < 0 ) theta_z = acos(normalDir[0]/sqrt(pow(normalDir[0],2.)+pow(normalDir[1],2.))); 
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
	}//else cout << "We are constructing a system with the (" << h_p << k_p << l_p << ") plane parallel with (xy) plane and the [" << buffer_vector_i[minhkl*3] << buffer_vector_i[minhkl*3+1] << buffer_vector_i[minhkl*3+2] << "] direction aligned with the x axis" << endl;
	cout << "FIRST GRAIN : " << endl;
	MT->printMat(rot_mat_total);
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

void Crystal::ComputeRotMat(int h, int k, int l, double theta, double *rotmat, bool IsDirOrPlane){
	double *rot_ax = new double[3];
	double *buffer = new double[9];
	ComputeVec(h,k,l,rot_ax,IsDirOrPlane);
	MT->Vec2rotMat(rot_ax,theta,buffer);
	// add the second rotation to be in the cartesian frame
	MT->MatDotMat(buffer,rot_mat_total,rotmat);
	delete[] rot_ax;
	delete[] buffer;
}

void Crystal::ComputeVec(int h, int k, int l, double *Vec, bool IsDirOrPlane, bool Normalize){
	int h_s, k_s, l_s;
	h_s = h;
	k_s = k;
	l_s = l;
	if( !IsDirOrPlane ) computeReciproqual();
	if( IsPrimitive ){
		if( IsDirOrPlane ) Convert2Primitive_dir(h_s,k_s,l_s);
		else Convert2Primitive_plane(h_s,k_s,l_s);
	}	
	double NormRotAx = 0.;
	for(unsigned int i=0;i<3;i++){
		if( IsDirOrPlane ) Vec[i] = a1[i]*h_s + a2[i]*k_s + a3[i]*l_s;
		else Vec[i] = a1_star[i]*h_s + a2_star[i]*k_s + a3_star[i]*l_s;
		if( Normalize ) NormRotAx += pow(Vec[i],2.);
	}
	if( Normalize ){
		NormRotAx = sqrt(NormRotAx);
		for(unsigned int i=0;i<3;i++) Vec[i] /= NormRotAx;
	}
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
	if( IsMolId ){
		cout << "Warning the DONOTSEPARE tag is used in addition with molecule id (atom_style full) in Crystal " << name << ", the molecule has the priority over the DONOTSEPARE tag which is thus ignored";
		return;
	}
	// verify that the provided DoNotSepare tag allows to respect the stoichiometry of the crystal (for instance 2 4 3 and 2 3 1 does not work for Mg2SiO4 because it will lead to Mg3SiO4 stoich)
	for(unsigned int t=0;t<this->nbAtomType;t++){
		unsigned int count = 0;
		for(unsigned int j=0;j<this->DoNotSep.size();j++){
			if( AtomType_uint[t] == DoNotSep[j][2] ){
				unsigned int current_t2;
				bool t2found = false;
				for(unsigned int t2=0;t2<this->nbAtomType;t2++){
					if( AtomType_uint[t2] == DoNotSep[j][0] ){
						current_t2 = t2;
						t2found = true;
						break;
					}
				}
				if( !t2found ){
					cerr << "Issue when searching atom type for constructing DoNotSep list of the current crystal" << endl;
					exit(EXIT_FAILURE);
				}
				count += Stoichiometry[current_t2]*DoNotSep[j][1];
			}
		}
		if( count > Stoichiometry[t] ){
			cerr << "The provided DoNotSepare informations in crystal \"" << name << "\" cannot be used as they imply to break the stoichiometry of the crystal" << endl;
		        cerr << "Please correct the DONOTSEPARE tag in /data/Crystals/" << name << ".ath file" << endl;
			cerr << "The current issue is coming from instruction for storing neighboring " << AtomType[t] << " atoms" << endl;
			cerr << "Aborting calculation" << endl;
			exit(EXIT_FAILURE);
		}
	}
		
	// construct full neighbor list in a large cutoff radius
	double arr[3] = {this->alength[0], this->alength[1], this->alength[2]};
	double rc_squared = pow(this->MT->max(arr,3),2.);
	double d_squared, xpos, ypos, zpos;
	int cl = 4;
	vector<vector<double>> buffer_vector_d;
	for(unsigned int i=0;i<this->nbAtom;i++){
		this->NotSepList.push_back(vector<int> ());
		buffer_vector_d.push_back(vector<double> ());
		bool Continue = true;
		for(unsigned int j=0;j<this->DoNotSep.size();j++){
			if( DoNotSep[j][0] == this->Motif[i].type_uint ){
				Continue = false;
				break;
			}
		}
		if( Continue ) continue;
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

	for(unsigned int i=0;i<this->nbAtom;i++)
		this->MT->sort(buffer_vector_d[i], 0, 5, buffer_vector_d[i]);

	// We should use an ion of the motif only once (independament of the CL) in order to generate more compact structures
	// search if the current crystal is degenerated (i.e. if it exists multiple possible combination of N first neighbours having the same distances) and construct the PossibleShuffling array
	bool degenerated = false;
	vector<vector<unsigned int>> DegeneratedAtom(nbAtom); // DegeneratedAtom[i][j] = index of the jth atom which could be associated with the neighbour having index i
	unsigned int count;
	bool Continue;
	for(unsigned int i=0;i<this->nbAtom;i++){
		for(unsigned int j=0;j<this->DoNotSep.size();j++){
			if( DoNotSep[j][0] == this->Motif[i].type_uint ){
				count = 0;
				for(unsigned int n=0;n<buffer_vector_d[i].size()/5;n++){
					if( count == this->DoNotSep[j][1] ) break;
					else if( this->Motif[(int) buffer_vector_d[i][n*5+1]].type_uint == DoNotSep[j][2] ){
						// verify that this atom has not been stored before
						Continue = false;
						for(unsigned int s_1=0;s_1<NotSepList.size();s_1++){
							for(unsigned int s_2=0;s_2<NotSepList[s_1].size()/4;s_2++){
								if( fabs(round(buffer_vector_d[i][n*5+1])-NotSepList[s_1][s_2*4]) < 1e-9 ){
									// This atom has already been stored before => degenerated case 
									degenerated = true;
									unsigned int current_neighid = (unsigned int) buffer_vector_d[i][n*5+1];
									bool already_1 = false;
									bool already_2 = false;
									for(unsigned int ne=0;ne<DegeneratedAtom[current_neighid].size();ne++){
										if( DegeneratedAtom[current_neighid][ne] == s_1 ) already_1 = true;
										if( DegeneratedAtom[current_neighid][ne] == i ) already_2 = true;
									}
									if( !already_1 ) DegeneratedAtom[current_neighid].push_back(s_1);
									if( !already_2 ) DegeneratedAtom[current_neighid].push_back(i);
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
				if( count != DoNotSep[j][1] ) cout << "Warning the crystal DoNotSep list has not been successfully computed !" << endl;
			}
		}
	}
	if( degenerated ){
		// Compute the possible combination of degenerated neighbours
		vector<vector<unsigned int>> PossibleCombs; // PossibleCombs[i][k*2] = id of neighbour for combination i and kth choice, [k*2+1] = id of atom
    		vector<unsigned int> current;
    		unordered_set<unsigned int> usedId2;
    		MT->GenerateCombinationsFromList(DegeneratedAtom, 0, current, usedId2, PossibleCombs);
    		
		double current_d_score;
		double min_d_score = 1.e6;
		unsigned int opt_comb;
		for(unsigned int ne1=0;ne1<PossibleCombs.size();ne1++){
			// initialize array and d_score
			vector<unsigned int> AlreadyStored(nbAtom,0);
			for(unsigned int i=0;i<nbAtom;i++) NotSepList[i].clear();
			current_d_score = 0.;
			// Assign neighbours to atoms for the current combination
			for(unsigned int ne2=0;ne2<PossibleCombs[ne1].size()/2;ne2++){
				unsigned int current_n_ind = PossibleCombs[ne1][ne2*2];
				unsigned int current_ind = PossibleCombs[ne1][ne2*2+1];
				AlreadyStored[current_n_ind] = 1;
				unsigned int current_index_in_buf;
				for(unsigned int i=0;i<buffer_vector_d[current_ind].size()/5;i++){
					if( (unsigned int) buffer_vector_d[current_ind][i*5+1] == current_n_ind ){
						current_index_in_buf = i;
						break;
					}
				}
				current_d_score += buffer_vector_d[current_ind][current_index_in_buf*5];
				NotSepList[current_ind].push_back(round(buffer_vector_d[current_ind][current_index_in_buf*5+1]));
				NotSepList[current_ind].push_back(round(buffer_vector_d[current_ind][current_index_in_buf*5+2]));
				NotSepList[current_ind].push_back(round(buffer_vector_d[current_ind][current_index_in_buf*5+3]));
				NotSepList[current_ind].push_back(round(buffer_vector_d[current_ind][current_index_in_buf*5+4]));
			}
			// Assign the rest of neighbours/atoms
			for(unsigned int i=0;i<this->nbAtom;i++){
				for(unsigned int j=0;j<this->DoNotSep.size();j++){
					if( DoNotSep[j][0] == this->Motif[i].type_uint ){
						count = 0;
						// search number of neighbour with this typr already stored
						for(unsigned int n=0;n<NotSepList[i].size()/4;n++){
							unsigned int current_neighid = (unsigned int) NotSepList[i][n*4];
							if( this->Motif[current_neighid].type_uint == DoNotSep[j][2] ) count++;
						}
						if( count == this->DoNotSep[j][1] ) continue;
						if( count > this->DoNotSep[j][1] ){
							cerr << "Warning something went wrong when computing DoNotSep list in crystal" << endl;
							exit(EXIT_FAILURE);
						}
						for(unsigned int n=0;n<buffer_vector_d[i].size()/5;n++){
							unsigned int current_neighid = (unsigned int) buffer_vector_d[i][n*5+1];
							if( count == this->DoNotSep[j][1] ) break;
							else if( this->Motif[current_neighid].type_uint == DoNotSep[j][2] && AlreadyStored[current_neighid] == 0 ){
								AlreadyStored[current_neighid] = 1;
								current_d_score += buffer_vector_d[i][n*5];
								NotSepList[i].push_back(round(buffer_vector_d[i][n*5+1]));
								NotSepList[i].push_back(round(buffer_vector_d[i][n*5+2]));
								NotSepList[i].push_back(round(buffer_vector_d[i][n*5+3]));
								NotSepList[i].push_back(round(buffer_vector_d[i][n*5+4]));
								count++;
							}
						}
						if( count != DoNotSep[j][1] ) cout << "Warning the crystal DoNotSep list has not been successfully computed !" << endl;
					}
				}
			}
			if( current_d_score < min_d_score ){
				opt_comb = ne1;
				min_d_score = current_d_score;
			}
		}
		
		// Store DoNotSep list with the optimal combination
		vector<unsigned int> AlreadyStored(nbAtom,0);
		for(unsigned int i=0;i<nbAtom;i++) NotSepList[i].clear();
		current_d_score = 0.;
		// Assign neighbours to atoms for the current combination
		for(unsigned int ne2=0;ne2<PossibleCombs[opt_comb].size()/2;ne2++){
			unsigned int current_n_ind = PossibleCombs[opt_comb][ne2*2];
			unsigned int current_ind = PossibleCombs[opt_comb][ne2*2+1];
			AlreadyStored[current_n_ind] = 1;
			unsigned int current_index_in_buf;
			for(unsigned int i=0;i<buffer_vector_d[current_ind].size()/5;i++){
				if( (unsigned int) buffer_vector_d[current_ind][i*5+1] == current_n_ind ){
					current_index_in_buf = i;
					break;
				}
			}
			current_d_score += buffer_vector_d[current_ind][current_index_in_buf*5];
			NotSepList[current_ind].push_back(round(buffer_vector_d[current_ind][current_index_in_buf*5+1]));
			NotSepList[current_ind].push_back(round(buffer_vector_d[current_ind][current_index_in_buf*5+2]));
			NotSepList[current_ind].push_back(round(buffer_vector_d[current_ind][current_index_in_buf*5+3]));
			NotSepList[current_ind].push_back(round(buffer_vector_d[current_ind][current_index_in_buf*5+4]));
		}
		// Assign the rest of neighbours/atoms
		for(unsigned int i=0;i<this->nbAtom;i++){
			for(unsigned int j=0;j<this->DoNotSep.size();j++){
				if( DoNotSep[j][0] == this->Motif[i].type_uint ){
					count = 0;
					// search number of neighbour with this typr already stored
					for(unsigned int n=0;n<NotSepList[i].size()/4;n++){
						unsigned int current_neighid = (unsigned int) NotSepList[i][n*4];
						if( this->Motif[current_neighid].type_uint == DoNotSep[j][2] ) count++;
					}
					if( count == this->DoNotSep[j][1] ) continue;
					if( count > this->DoNotSep[j][1] ){
						cerr << "Warning something went wrong when computing DoNotSep list in crystal" << endl;
						exit(EXIT_FAILURE);
					}
					for(unsigned int n=0;n<buffer_vector_d[i].size()/5;n++){
						unsigned int current_neighid = (unsigned int) buffer_vector_d[i][n*5+1];
						if( count == this->DoNotSep[j][1] ) break;
						else if( this->Motif[current_neighid].type_uint == DoNotSep[j][2] && AlreadyStored[current_neighid] == 0 ){
							AlreadyStored[current_neighid] = 1;
							current_d_score += buffer_vector_d[i][n*5];
							NotSepList[i].push_back(round(buffer_vector_d[i][n*5+1]));
							NotSepList[i].push_back(round(buffer_vector_d[i][n*5+2]));
							NotSepList[i].push_back(round(buffer_vector_d[i][n*5+3]));
							NotSepList[i].push_back(round(buffer_vector_d[i][n*5+4]));
							count++;
						}
					}
					if( count != DoNotSep[j][1] ) cout << "Warning the crystal DoNotSep list has not been successfully computed !" << endl;
				}
			}
		}
		if( fabs(current_d_score-min_d_score) > 1e-9 ) cout << "Warning something went wrong" << endl;


	}
}

void Crystal::ConstructNotSepListFromMolId(){
	this->IsDoNotSep = true;
	vector<unsigned int> index_notused(nbAtom);
	for(unsigned int i=0;i<this->nbAtom;i++) index_notused[i] = i;
	for(unsigned int i=0;i<this->nbAtom;i++){
		this->NotSepList.push_back(vector<int> ());
		vector<unsigned int> index_torm;
		for(unsigned int j=0;j<index_notused.size();j++){
			if( index_notused[j] == i ){
				index_notused.erase(index_notused.begin()+j);
				break;
			}
		}
		for(unsigned int j=0;j<index_notused.size();j++){
			if( MolId[i] == MolId[index_notused[j]] ){
				NotSepList[i].push_back(index_notused[j]);
				index_torm.push_back(j);
				for(unsigned int d=0;d<3;d++) NotSepList[i].push_back(0); // no boundary conditions to apply
			}
		}
		for(unsigned int j=0;j<index_torm.size();j++) index_notused.erase(index_notused.begin()+index_torm[index_torm.size()-j-1]);
	}

}

void Crystal::ConstructOrthogonalCell(bool AllowTilt){
	// search the smallest orthogonal box for this orientation by testing linear combination of a1 a2 a3 giving cell vectors aligned with cartesian axis
	// expect for the x axis which is already aligned with a crystallographic direction
	vector<int> cl_box;
	vector<int> xh_list;
	vector<int> yh_list;
	vector<int> zh_list;
	vector<double> buffer_vector_d_x, buffer_vector_d_y, buffer_vector_d_z;
	double absX, absY, absZ, vecL;
	// search smallest lattice param
	double min_length = MT->min_p(alength,3)*0.9;	
	for(int i=-CLsearch;i<CLsearch+1;i++){ // TODO change this loop with a while to exit when vectors are found with a maximum number of search
		for(int j=-CLsearch;j<CLsearch+1;j++){
			for(int k=-CLsearch;k<CLsearch+1;k++){
				absX = i*this->a1[0]+j*this->a2[0]+k*this->a3[0];
				absY = i*this->a1[1]+j*this->a2[1]+k*this->a3[1];
				absZ = i*this->a1[2]+j*this->a2[2]+k*this->a3[2];
				vecL = absX*absX + absY*absY + absZ*absZ;
				absX = fabs(absX);
				absY = fabs(absY);
				absZ = fabs(absZ);
				if( absY < TolOrthoBox && absZ < TolOrthoBox && absX > min_length ){ // TODO change tolerance for being a deformation more than a distance
					buffer_vector_d_x.push_back(vecL);
					xh_list.push_back(i);
					xh_list.push_back(j);
					xh_list.push_back(k);
				}
				bool cond = absZ < TolOrthoBox;
				if( !AllowTilt ) cond = cond && absX < TolOrthoBox;
				if( cond && absY > min_length ){
					buffer_vector_d_y.push_back(vecL);
					yh_list.push_back(i);
					yh_list.push_back(j);
					yh_list.push_back(k);
				}
				cond = true;
				if( !AllowTilt ) cond = absX < TolOrthoBoxZ && absY < TolOrthoBoxZ;
				if( cond && absZ > min_length ){
					buffer_vector_d_z.push_back(vecL);
					zh_list.push_back(i);
					zh_list.push_back(j);
					zh_list.push_back(k);
				}
			}
		}
	}
	if( xh_list.size() == 0 || yh_list.size() == 0 || zh_list.size() == 0 ){
		cerr << "There is no linear combination of crystal vectors giving an orthogonal cell, consider increase the tolerances (TOL_ORTHO_BOX or TOL_ORTHO_BOX_Z) or the number of linear combinations used for the research (NB_MAX_LC) in a FixedParameters.ath file, aborting computation" << endl;
		exit(EXIT_FAILURE);
	}

	unsigned int ind_x, ind_y, ind_z;
	ind_x = MT->min(buffer_vector_d_x);
	ind_y = MT->min(buffer_vector_d_y);
	ind_z = MT->min(buffer_vector_d_z);
	// fill the cl_box array
	int arr1[3] = {abs(xh_list[ind_x*3]), abs(yh_list[ind_y*3]), abs(zh_list[ind_z*3])};
	cl_box.push_back(MT->max(arr1,3));
	int arr2[3] = {abs(xh_list[ind_x*3+1]), abs(yh_list[ind_y*3+1]), abs(zh_list[ind_z*3+1])};
	cl_box.push_back(MT->max(arr2,3));
	int arr3[3] = {abs(xh_list[ind_x*3+2]), abs(yh_list[ind_y*3+2]), abs(zh_list[ind_z*3+2])};
	cl_box.push_back(MT->max(arr3,3));
	double *xh = new double[3];
	double *yh = new double[3];
	double *zh = new double[3];
	for(unsigned int i=0;i<3;i++){
		xh[i] = this->a1[i]*xh_list[ind_x*3]+this->a2[i]*xh_list[ind_x*3+1]+this->a3[i]*xh_list[ind_x*3+2];
		yh[i] = this->a1[i]*yh_list[ind_y*3]+this->a2[i]*yh_list[ind_y*3+1]+this->a3[i]*yh_list[ind_y*3+2];
		zh[i] = this->a1[i]*zh_list[ind_z*3]+this->a2[i]*zh_list[ind_z*3+1]+this->a3[i]*zh_list[ind_z*3+2];
	}
	if( xh[0] < 0. ){
		for(unsigned int i=0;i<3;i++){
			xh_list[ind_x*3+i] *= -1;
			xh[i] *= -1.;
		}
	}
	if( yh[1] < 0. ){
		for(unsigned int i=0;i<3;i++){
			yh_list[ind_y*3+i] *= -1;
			yh[i] *= -1.;
		}
	}
	if( zh[2] < 0. ){
		for(unsigned int i=0;i<3;i++){
			zh_list[ind_z*3+i] *= -1;
			zh[i] *= -1.;
		}
	}
	// add conditions required by lammps for triclinic boxes TODO
	if( AllowTilt ){
		if( fabs(yh[0]) > xh[0]/2. || fabs(zh[0]) > xh[0]/2. || fabs(zh[1]) > yh[1]/2. ) cout << "Warning, the box is too skewed for being computed using LAMMPS.." << endl;
	}

	double xbox = sqrt( pow(xh[0],2.) + pow(xh[1],2.) + pow(xh[2],2.));
	double ybox = sqrt( pow(yh[0],2.) + pow(yh[1],2.) + pow(yh[2],2.));
	double zbox = sqrt( pow(zh[0],2.) + pow(zh[1],2.) + pow(zh[2],2.));

	MT->computeTiltTrans(xh,yh,zh,this->TiltTrans_xyz, AllowTilt);

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
	if( AllowTilt ){
		if( fabs(xh[1]) > tolBase || fabs(xh[2]) > tolBase || fabs(yh[2]) > tolBase ) cout << "Warning there were issues during the construction of orthogonal cell" << endl;
	}else{
		if( fabs(xh[1]) > tolBase || fabs(xh[2]) > tolBase || fabs(yh[0]) > tolBase || fabs(yh[2]) > tolBase || fabs(zh[0]) > tolBase || fabs(zh[1]) > tolBase || fabs(xh[0]-xbox) > tolBase || fabs(yh[1]-ybox) > tolBase || fabs(zh[2]-zbox) > tolBase ) cout << "Warning there were issues during the construction of orthogonal cell" << endl;
	}
	// recompute reciproqual space as cell parameters may have changed
	computeReciproqual();
	// compute the total transformation matrix
	double *TotalTrans = new double[9];
	MT->MatDotMat(TiltTrans_xyz,this->rot_mat_total,TotalTrans);
	// shift the motif
	ShiftMotif(shift_x,shift_y,shift_z);
	// rotate the motif
	for(unsigned int i=0;i<this->nbAtom;i++) MT->MatDotAt(TotalTrans,Motif[i],Motif[i]);	
	// construct the AtomicSystem
	this->OrientedSystem = new AtomicSystem(this, xh, yh, zh, cl_box);
	// rescale the crystal vectors as they have been modified by the function
	double *invTrans = new double[9];
	MT->invert3x3(this->TiltTrans_xyz,invTrans);
	//MT->MatDotRawVec(invTrans,this->a1,this->a1);
	//MT->MatDotRawVec(invTrans,this->a2,this->a2);
	//MT->MatDotRawVec(invTrans,this->a3,this->a3);
	// recompute reciproqual space as cell parameters may have changed
	//computeReciproqual();
	ComputeCrystalDef();
	this->IsOrientedSystem = true;
	delete[] invTrans;
	delete[] TotalTrans;
	delete[] xh;
	delete[] yh;
	delete[] zh;
}

void Crystal::ComputeCrystalDef(){
	crystal_def[0] = (TiltTrans_xyz[0]-1.)*100.;
	crystal_def[1] = (TiltTrans_xyz[4]-1.)*100.;
	crystal_def[2] = (TiltTrans_xyz[8]-1.)*100.;
	crystal_def[3] = pow(TiltTrans_xyz[1],2.)+pow(TiltTrans_xyz[2],2.)+pow(TiltTrans_xyz[5],2.);
	crystal_def[3] += (pow(TiltTrans_xyz[4]-TiltTrans_xyz[8],2.)+pow(TiltTrans_xyz[0]-TiltTrans_xyz[8],2.)+pow(TiltTrans_xyz[0]-TiltTrans_xyz[4],2.))/8.;
	crystal_def[3] = sqrt(crystal_def[3])*100.;
}

void Crystal::computeReciproqual(){
	// compute cell volume
	double *mixedProd = new double[3];
	MT->mixedProd(this->a2,this->a3,mixedProd);
	this->V = fabs(MT->dotProd(this->a1,mixedProd));
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

void Crystal::SearchCrystallographicPlanesAndDirections(double *Dir, int &h, int &k, int &l, int &u, int &v, int &w){
	computeReciproqual();
	// search the crystalligraphic plane and direction aligned with a given direction
	// initialize indices to zero (if not found the plane and dir will be zero)
	h = 0; k = 0; l = 0; u = 0; v = 0; w = 0;
	double *DirNorm = new double[3];
	double norm = 0.;
	for(unsigned int i=0;i<3;i++){
		norm += Dir[i]*Dir[i];
		DirNorm[i] = Dir[i];
	}
	for(unsigned int i=0;i<3;i++) DirNorm[i] /= sqrt(norm);

	vector<int> diff_list;
	vector<int> diff_list_p;
	vector<double> buffer_vector_d, buffer_vector_d_p;
	double sp, X, Y, Z;
	double tol_angle = cos(2.*M_PI/180.);
	for(int i=-max_hkl_search;i<max_hkl_search+1;i++){
		for(int j=-max_hkl_search;j<max_hkl_search+1;j++){
			for(int k=-max_hkl_search;k<max_hkl_search+1;k++){
				if( i != 0 || j !=0 || k != 0 ){
					// directions
					X = (i*this->a1[0]+j*this->a2[0]+k*this->a3[0]);
					Y = (i*this->a1[1]+j*this->a2[1]+k*this->a3[1]);
					Z = (i*this->a1[2]+j*this->a2[2]+k*this->a3[2]);
					//sp = fabs(acos(( X*DirNorm[0] + Y*DirNorm[1] + Z*DirNorm[2] ) / sqrt(X*X + Y*Y + Z*Z)));
					sp = ( X*DirNorm[0] + Y*DirNorm[1] + Z*DirNorm[2] ) / sqrt(X*X + Y*Y + Z*Z);
					if( sp > tol_angle ){
						buffer_vector_d.push_back(sp);
						diff_list.push_back(i);
						diff_list.push_back(j);
						diff_list.push_back(k);
					}
					// planes
					X = i*this->a1_star[0]+j*this->a2_star[0]+k*this->a3_star[0];
					Y = i*this->a1_star[1]+j*this->a2_star[1]+k*this->a3_star[1];
					Z = i*this->a1_star[2]+j*this->a2_star[2]+k*this->a3_star[2];
					//sp = fabs(acos(( X*DirNorm[0] + Y*DirNorm[1] + Z*DirNorm[2] ) / sqrt(X*X + Y*Y + Z*Z)));
					sp = ( X*DirNorm[0] + Y*DirNorm[1] + Z*DirNorm[2] ) / sqrt(X*X + Y*Y + Z*Z);
					if( sp > tol_angle ){
						buffer_vector_d_p.push_back(sp);
						diff_list_p.push_back(i);
						diff_list_p.push_back(j);
						diff_list_p.push_back(k);
					}
				}
			}
		}
	}
	if( diff_list.size() == 0 ) cout << "We cannot find crystallographic directions aligned with the provided direction (increase MAX_HKL_SEARCH parameter in a FixedParameters.ath file)" << endl;
	if( diff_list_p.size() == 0 ) cout << "We cannot find crystallographic planes aligned with the provided direction (increase MAX_HKL_SEARCH parameter in a FixedParameters.ath file)" << endl;

	if( diff_list.size() != 0 && diff_list_p.size() != 0 ){
		unsigned int ind;
		int *buffer_vec = new int[3];
		// directions
		ind = MT->max(buffer_vector_d);
		for(unsigned int j=0;j<3;j++) buffer_vec[j] = diff_list[ind*3+j];
		MT->reduce_vec(buffer_vec,buffer_vec);
		if( buffer_vec[0] < 0 ) for(unsigned int j=0;j<3;j++) buffer_vec[j] *= -1;
		else if( buffer_vec[0] == 0 ){
			if( buffer_vec[1] < 0 ) for(unsigned int j=0;j<3;j++) buffer_vec[j] *= -1;
		        else if( buffer_vec[1] == 0 && buffer_vec[2] < 0 ) for(unsigned int j=0;j<3;j++) buffer_vec[j] *= -1;
		}
		h = buffer_vec[0];
		k = buffer_vec[1];
		l = buffer_vec[2];
		// planes
		ind = MT->max(buffer_vector_d_p);
		for(unsigned int j=0;j<3;j++) buffer_vec[j] = diff_list_p[ind*3+j];
		MT->reduce_vec(buffer_vec,buffer_vec);
		if( buffer_vec[0] < 0 ) for(unsigned int j=0;j<3;j++) buffer_vec[j] *= -1;
		else if( buffer_vec[0] == 0 ){
			if( buffer_vec[1] < 0 ) for(unsigned int j=0;j<3;j++) buffer_vec[j] *= -1;
		        else if( buffer_vec[1] == 0 && buffer_vec[2] < 0 ) for(unsigned int j=0;j<3;j++) buffer_vec[j] *= -1;
		}
		u = buffer_vec[0];
		v = buffer_vec[1];
		w = buffer_vec[2];
		delete[] buffer_vec;
	}
	if( IsPrimitive ){
		Convert2NonPrimitive_dir(h,k,l);
		Convert2NonPrimitive_plane(u,v,w);
	}
	delete[] DirNorm;
}

void Crystal::Convert2NonPrimitive_dir(int &h, int &k, int &l){
	int h_s, k_s, l_s;
	h_s = h;
	k_s = k;
	l_s = l;
	if( crystallo == "FCC" ){
		h = 2 * k_s + 2 * h_s;
		k = 2 * k_s + 2 * l_s;
		l = 2 * h_s + 2 * l_s;
		int gcd = MT->gcd3(h,k,l);
		h /= gcd;
		k /= gcd;
		l /= gcd;
	}else if( crystallo == "BCC" ){
		cout << "BCC_conv <=> BCC_prim to be implemented TODO" << endl;
	}
}

void Crystal::Convert2Primitive_dir(int &h, int &k, int &l){
	int h_s, k_s, l_s;
	h_s = h;
	k_s = k;
	l_s = l;
	if( crystallo == "FCC" ){
		h = h_s - k_s + l_s;
		k = h_s + k_s - l_s;
		l = -h_s + k_s + l_s;
		int gcd = MT->gcd3(h,k,l);
		h /= gcd;
		k /= gcd;
		l /= gcd;
	}else if( crystallo == "BCC" ){
		cout << "BCC_conv <=> BCC_prim to be implemented TODO" << endl;
	}
}

void Crystal::Convert2NonPrimitive_plane(int &h, int &k, int &l){
	int h_s, k_s, l_s;
	h_s = h;
	k_s = k;
	l_s = l;
	if( crystallo == "FCC" ){
		h = h_s + k_s - l_s;
		k = - h_s + k_s + l_s;
		l = h_s - k_s + l_s;
		int gcd = MT->gcd3(h,k,l);
		h /= gcd;
		k /= gcd;
		l /= gcd;
	}else if( crystallo == "BCC" ){
		cout << "BCC_conv <=> BCC_prim to be implemented TODO" << endl;
	}
}

void Crystal::Convert2Primitive_plane(int &h, int &k, int &l){
	int h_s, k_s, l_s;
	h_s = h;
	k_s = k;
	l_s = l;
	if( crystallo == "FCC" ){
		h = 2 * h_s + 2 * l_s;
		k = 2 * h_s + 2 * k_s;
		l = 2 * k_s + 2 * l_s;
		int gcd = MT->gcd3(h,k,l);
		h /= gcd;
		k /= gcd;
		l /= gcd;
	}else if( crystallo == "BCC" ){
		cout << "BCC_conv <=> BCC_prim to be implemented TODO" << endl;
	}
}
void Crystal::ComputeOrthogonalPlanesAndDirections(){
	// search which crystallographic planes are aligned with the x, y, and z planes and which crystallographic directions are aligned with the x, y and z directions
	// initialize indices to zero (if not found the plane and dir will be zero)
	for(unsigned int i=0;i<9;i++){
		OrthogonalPlanes[i] = 0;		
		OrthogonalDirs[i] = 0;
	}

	computeReciproqual();
	vector<int> xh_list;
	vector<int> yh_list;
	vector<int> zh_list;
	vector<int> xh_list_p;
	vector<int> yh_list_p;
	vector<int> zh_list_p;
	vector<double> buffer_vector_d_x, buffer_vector_d_y, buffer_vector_d_z, buffer_vector_d_x_p, buffer_vector_d_y_p, buffer_vector_d_z_p;
	double absX, absY, absZ, X, Y, Z;
	for(int i=-max_hkl_search;i<max_hkl_search+1;i++){
		for(int j=-max_hkl_search;j<max_hkl_search+1;j++){
			for(int k=-max_hkl_search;k<max_hkl_search+1;k++){
				if( i != 0 || j !=0 || k != 0 ){
					// directions
					X = i*this->a1[0]+j*this->a2[0]+k*this->a3[0];
					Y = i*this->a1[1]+j*this->a2[1]+k*this->a3[1];
					Z = i*this->a1[2]+j*this->a2[2]+k*this->a3[2];
					absX = fabs(X);
					absY = fabs(Y);
					absZ = fabs(Z);
					if( absY < TolOrthoBox && absZ < TolOrthoBoxZ && X > 0 ){
						buffer_vector_d_x.push_back(absY+absZ);
						xh_list.push_back(i);
						xh_list.push_back(j);
						xh_list.push_back(k);
					}
					if( absX < TolOrthoBox && absY < TolOrthoBox && Z > 0 ){
						buffer_vector_d_z.push_back(absX+absY);
						zh_list.push_back(i);
						zh_list.push_back(j);
						zh_list.push_back(k);
					}
					if( absX < TolOrthoBox && absZ < TolOrthoBoxZ && Y > 0 ){
						buffer_vector_d_y.push_back(absX+absZ);
						yh_list.push_back(i);
						yh_list.push_back(j);
						yh_list.push_back(k);
					}
					// planes
					X = i*this->a1_star[0]+j*this->a2_star[0]+k*this->a3_star[0];
					Y = i*this->a1_star[1]+j*this->a2_star[1]+k*this->a3_star[1];
					Z = i*this->a1_star[2]+j*this->a2_star[2]+k*this->a3_star[2];
					absX = fabs(X);
					absY = fabs(Y);
					absZ = fabs(Z);
					if( absY < TolOrthoBox && absZ < TolOrthoBoxZ && X > 0 ){
						buffer_vector_d_x_p.push_back(absY+absZ);
						xh_list_p.push_back(i);
						xh_list_p.push_back(j);
						xh_list_p.push_back(k);
					}
					if( absX < TolOrthoBox && absY < TolOrthoBox && Z > 0 ){
						buffer_vector_d_z_p.push_back(absX+absY);
						zh_list_p.push_back(i);
						zh_list_p.push_back(j);
						zh_list_p.push_back(k);
					}
					if( absX < TolOrthoBox && absZ < TolOrthoBoxZ && Y > 0 ){
						buffer_vector_d_y_p.push_back(absX+absZ);
						yh_list_p.push_back(i);
						yh_list_p.push_back(j);
						yh_list_p.push_back(k);
					}
				}
			}
		}
	}
	int *buffer_vec = new int[3];
	unsigned int ind;
	if( xh_list.size() == 0 ) cout << "We cannot find crystallographic direction aligned with x cartesian axis (increase MAX_HKL_SEARCH parameter in a FixedParameters.ath file)" << endl;
	else{
		ind = MT->min(buffer_vector_d_x);
		for(unsigned int j=0;j<3;j++) buffer_vec[j] = xh_list[ind*3+j];
		MT->reduce_vec(buffer_vec,buffer_vec);
		for(unsigned int j=0;j<3;j++) OrthogonalDirs[j] = buffer_vec[j];
	}
	if( yh_list.size() == 0 ) cout << "We cannot find crystallographic direction aligned with y cartesian axis (increase MAX_HKL_SEARCH parameter in a FixedParameters.ath file)" << endl;
	else{
		ind = MT->min(buffer_vector_d_y);
		for(unsigned int j=0;j<3;j++) buffer_vec[j] = yh_list[ind*3+j];
		MT->reduce_vec(buffer_vec,buffer_vec);
		for(unsigned int j=0;j<3;j++) OrthogonalDirs[j+3] = buffer_vec[j];
	}
	if( zh_list.size() == 0 ) cout << "We cannot find crystallographic direction aligned with z cartesian axis (increase MAX_HKL_SEARCH parameter in a FixedParameters.ath file)" << endl;
	else{
		ind = MT->min(buffer_vector_d_z);
		for(unsigned int j=0;j<3;j++) buffer_vec[j] = zh_list[ind*3+j];
		MT->reduce_vec(buffer_vec,buffer_vec);
		for(unsigned int j=0;j<3;j++) OrthogonalDirs[j+6] = buffer_vec[j];
	}
	if( xh_list_p.size() == 0 ) cout << "We cannot find crystallographic plane aligned with x cartesian axis (increase MAX_HKL_SEARCH parameter in a FixedParameters.ath file)" << endl;
	else{
		ind = MT->min(buffer_vector_d_x_p);
		for(unsigned int j=0;j<3;j++) buffer_vec[j] = xh_list_p[ind*3+j];
		MT->reduce_vec(buffer_vec,buffer_vec);
		for(unsigned int j=0;j<3;j++) OrthogonalPlanes[j] = buffer_vec[j];
	}
	if( yh_list_p.size() == 0 ) cout << "We cannot find crystallographic plane aligned with y cartesian axis (increase MAX_HKL_SEARCH parameter in a FixedParameters.ath file)" << endl;
	else{
		ind = MT->min(buffer_vector_d_y_p);
		for(unsigned int j=0;j<3;j++) buffer_vec[j] = yh_list_p[ind*3+j];
		MT->reduce_vec(buffer_vec,buffer_vec);
		for(unsigned int j=0;j<3;j++) OrthogonalPlanes[j+3] = buffer_vec[j];
	}
	if( zh_list_p.size() == 0 ) cout << "We cannot find crystallographic plane aligned with z cartesian axis (increase MAX_HKL_SEARCH parameter in a FixedParameters.ath file)" << endl;
	else{
		ind = MT->min(buffer_vector_d_z_p);
		for(unsigned int j=0;j<3;j++) buffer_vec[j] = zh_list_p[ind*3+j];
		MT->reduce_vec(buffer_vec,buffer_vec);
		for(unsigned int j=0;j<3;j++) OrthogonalPlanes[j+6] = buffer_vec[j];
	}
	if( IsPrimitive ){
		for(unsigned int i=0;i<3;i++){
			Convert2NonPrimitive_plane(OrthogonalPlanes[i*3], OrthogonalPlanes[i*3+1], OrthogonalPlanes[i*3+2]);
			Convert2NonPrimitive_dir(OrthogonalDirs[i*3], OrthogonalDirs[i*3+1], OrthogonalDirs[i*3+2]);
		}
	}
	delete[] buffer_vec;
}

void Crystal::read_database(){
	ifstream file(this->path2database, ios::in);
	size_t pos_at, pos_x, pos_y, pos_z, pos_attype, pos_Mass, pos_At, pos_Crystal, pos_tilt, pos_DNS, pos_bondori, pos_bond, pos_bondtype, pos_angle, pos_angletype, pos_Bond, pos_Angle;
	unsigned int h_uint = 1e8;
	unsigned int line_Mass(h_uint), line_At(h_uint), buffer_uint, buffer_uint_1, buffer_uint_2, buffer_uint_3, buffer_uint_4, count(0), line_bondori(h_uint), nbref2read(0), buffer_molid, line_Bond(h_uint), line_Angle(h_uint);
	double buffer_1, buffer_2, buffer_3, buffer_4;
	string buffer_s, buffer_s_1, buffer_s_2, line;
	for(unsigned int i=0;i<3;i++){
		this->a1[i] = 0.;
		this->a2[i] = 0.;
		this->a3[i] = 0.;
	}
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
			pos_bond=line.find("bonds");
			if(pos_bond!=string::npos){
				istringstream text(line);
				text >> buffer_uint;
				this->IsBond = true;
				this->nbBonds = buffer_uint;
				Bonds = new unsigned int[nbBonds*2];
				BondType = new unsigned int[nbBonds];
			}
			pos_bondtype=line.find("bond type");
			if(pos_bondtype!=string::npos){
				istringstream text(line);
				text >> buffer_uint;
				this->nbBondType = buffer_uint;
			}
			pos_angle=line.find("angles");
			if(pos_angle!=string::npos){
				istringstream text(line);
				text >> buffer_uint;
				this->IsAngle = true;
				this->nbAngles = buffer_uint;
				Angles = new unsigned int[nbAngles*3];
				AngleType = new unsigned int[nbAngles];
			}
			pos_angletype=line.find("angle type");
			if(pos_angletype!=string::npos){
				istringstream text(line);
				text >> buffer_uint;
				this->nbAngleType = buffer_uint;
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
			pos_tilt=line.find("xy xz yz");
			if(pos_tilt!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2 >> buffer_3;
				this->a2[0] = buffer_1;
				this->a3[0] = buffer_2;
				this->a3[1] = buffer_3;
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
				string buffer_site;
				text >> buffer_s >> buffer_uint >> buffer_uint_1 >> buffer_site;
				this->DoNotSep.push_back(vector<unsigned int> ());
				this->DoNotSep[this->DoNotSep.size()-1].push_back(buffer_uint);
				this->DoNotSep[this->DoNotSep.size()-1].push_back(buffer_uint_1);
				size_t pos_und = buffer_site.find("_");
				if( pos_und!=string::npos ){
					replace(buffer_site.begin(), buffer_site.end(), '_', ' ');
					istringstream text2(buffer_site);
					text2 >> buffer_uint_2;
					this->DoNotSep[this->DoNotSep.size()-1].push_back(buffer_uint_2);
					text2 >> buffer_uint_2;
					this->DoNotSep[this->DoNotSep.size()-1].push_back(buffer_uint_2);
				}else{
					istringstream text2(buffer_site);
					text2 >> buffer_uint_2;
					this->DoNotSep[this->DoNotSep.size()-1].push_back(buffer_uint_2);
				}
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
				else if( buffer_s == "full" ){
					this->IsMolId = true;
					MolId = new unsigned int[nbAtom];
					this->IsCharge = true;
				}
			       	line_At = count;
			}
			if( count > line_At+1 && count < line_At+nbAtom+2 ){
				istringstream text(line);
				text >> buffer_uint;
				if( this->IsMolId ){
					text >> buffer_molid;
					MolId[buffer_uint-1] = buffer_molid;
				}
				text >> buffer_uint_1;
				if( this->IsCharge ){
					text >> buffer_1;
					this->AtomCharge[buffer_uint_1-1] = buffer_1;
				}
				text >> buffer_2 >> buffer_3 >> buffer_4 >> buffer_uint_2;
				this->Motif[buffer_uint-1].pos.x = buffer_2;
				this->Motif[buffer_uint-1].pos.y = buffer_3;
				this->Motif[buffer_uint-1].pos.z = buffer_4;
				this->Motif[buffer_uint-1].type_uint = buffer_uint_1;
				this->AtomSite[buffer_uint-1] = buffer_uint_2-1;
			}
			// Read bonds
			if( IsBond ){
				pos_Bond=line.find("Bonds");
				if(pos_Bond!=string::npos) line_Bond = count;
				if( count > line_Bond+1 && count < line_Bond+nbBonds+2 ){
					istringstream text(line);
					text >> buffer_uint >> buffer_uint_1 >> buffer_uint_2 >> buffer_uint_3;
					BondType[buffer_uint-1] = buffer_uint_1;
					Bonds[(buffer_uint-1)*2] = buffer_uint_2;
					Bonds[(buffer_uint-1)*2+1] = buffer_uint_3;
				}
			}
			// Read angles
			if( IsAngle ){
				pos_Angle=line.find("Angles");
				if(pos_Angle!=string::npos) line_Angle = count;
				if( count > line_Angle+1 && count < line_Angle+nbAngles+2 ){
					istringstream text(line);
					text >> buffer_uint >> buffer_uint_1 >> buffer_uint_2 >> buffer_uint_3 >> buffer_uint_4;
					AngleType[buffer_uint-1] = buffer_uint_1;
					Angles[(buffer_uint-1)*3] = buffer_uint_2;
					Angles[(buffer_uint-1)*3+1] = buffer_uint_3;
					Angles[(buffer_uint-1)*3+2] = buffer_uint_4;
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
	}else{
		cout << "The file " << this->path2database << " cannot be openned" << endl;
	}
	// check if the provided crystal is a primitive cell
	// for the moment only consider monoatomic cubic cells:
	// - check if there is only one atom and if cell vectors are not normals
	if( ( crystallo == "FCC" || crystallo == "BCC" ) && nbAtom != 1 ){
		double sp1(0.), sp2(0.), sp3(0.);
		double tol_sp = 1e-3;
		for(unsigned int i=0;i<3;i++){
			sp1 += a1[i]*a2[i];
			sp2 += a1[i]*a3[i];
			sp3 += a3[i]*a2[i];
		}
		if( fabs(sp1) < tol_sp && fabs(sp2) < tol_sp && fabs(sp3) < tol_sp ){
			cout << "Warning, the provided cubic crystal is not a primitive cell, which may lead to calculation of non-primitive lattices related to this crystal (e.g. CSL, DSC)" << endl;
			IsPrimitive = false;
		}
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

double Crystal::ComputeD_hkl(const int& h, const int& k, const int& l){
	int h_s, k_s, l_s;
	h_s = h;
	k_s = k;
	l_s = l;
	if( IsPrimitive ) Convert2Primitive_plane(h_s,k_s,l_s);
	double dist = 0.;
	for(unsigned int i=0;i<3;i++)
		dist += ((h_s*a1_star[i]) + (k_s*a2_star[i]) + (l_s*a3_star[i])) * ((h_s*a1_star[i]) + (k_s*a2_star[i]) + (l_s*a3_star[i]));
	dist = sqrt(dist);
	return 1./dist;
}

void Crystal::read_params(){
	//string fp;
	//#ifdef FIXEDPARAMETERS
	//fp = FIXEDPARAMETERS;
	//#endif
	//string backslash="/";
	//string filename=fp+backslash+FixedParam_Filename;
	//ifstream file(filename, ios::in);
	ifstream file(FixedParam_Filename, ios::in);
	size_t pos_tolOrthoBox, pos_tolOrthoBoxZ, pos_shift, pos_nblc, pos_max_hkl_search;
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
			pos_shift=line.find("MOTIF_SHIFT");
			if(pos_shift!=string::npos){
				istringstream text(line);
				text >> buffer_s >> shift_x >> shift_y >> shift_z;
			}
			pos_nblc=line.find("NB_MAX_LC");
			if(pos_nblc!=string::npos){
				istringstream text(line);
				text >> buffer_s >> CLsearch;
			}
			pos_max_hkl_search=line.find("MAX_HKL_SEARCH");
			if(pos_max_hkl_search!=string::npos){
				istringstream text(line);
				text >> buffer_s >> max_hkl_search;
			}

		}
		file.close();
	}//else{
	//	cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
	//	exit(EXIT_FAILURE);
	//}
}

void Crystal::ReadProperties(vector<string> Properties){
	size_t pos_tolOrthoBox, pos_tolOrthoBoxZ, pos_shift, pos_nblc, pos_max_hkl_search;
	string buffer_s;
	for(unsigned int i=0;i<Properties.size();i++){
		pos_tolOrthoBox=Properties[i].find("TOL_ORTHO_BOX ");
		if(pos_tolOrthoBox!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->TolOrthoBox;
		}
		pos_tolOrthoBoxZ=Properties[i].find("TOL_ORTHO_BOX_Z ");
		if(pos_tolOrthoBoxZ!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->TolOrthoBoxZ;
		}
		pos_shift=Properties[i].find("MOTIF_SHIFT");
		if(pos_shift!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> shift_x >> shift_y >> shift_z;
		}
		pos_nblc=Properties[i].find("NB_MAX_LC");
		if(pos_nblc!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> CLsearch;
		}
		pos_max_hkl_search=Properties[i].find("MAX_HKL_SEARCH");
		if(pos_max_hkl_search!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> max_hkl_search;
		}

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

void Crystal::ShiftMotif(const double &_shift_x, const double &_shift_y, const double &_shift_z){
	for(unsigned int i=0;i<this->nbAtom;i++){
		Motif[i].pos.x += _shift_x;
		Motif[i].pos.y += _shift_y;
		Motif[i].pos.z += _shift_z;
	}
}

Crystal::~Crystal(){
	delete[] a1;
	delete[] a2;
	delete[] a3;
	delete[] a1_star;
	delete[] a2_star;
	delete[] a3_star;
	delete[] OrthogonalPlanes;
	delete[] OrthogonalDirs;
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
	if( AtomSite ) delete[] AtomSite;
	if( AtomCharge ) delete[] AtomCharge;
	delete[] crystal_def;
	if( IsBond ){
		delete[] Bonds;
		delete[] BondType;
	}
	if( IsAngle ){
		delete[] Angles;
		delete[] AngleType;
	}
	if( IsMolId ) delete[] MolId;
	if( IsOrientedSystem ) delete OrientedSystem;
	if( IsReferenceBondOriParam ){
		for(unsigned int t=0;t<MaxAtomType;t++) ReferenceBondOriParam[t].clear();
		delete[] ReferenceBondOriParam;
	}
}
