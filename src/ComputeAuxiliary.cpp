//**********************************************************************************
//*   ComputeAuxiliary.cpp                                                         *
//**********************************************************************************
//* This file contains the implementation of the ComputeAuxiliary class            *
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


#include "ComputeAuxiliary.h"
#include "AtomHicConfig.h"
#include "Crystal.h"
#include <iostream>
#include <vector>
#include <iomanip>
#include <omp.h>
#include <cstring>
#include <sstream>
#include <fstream>
#include <dirent.h>

using namespace std;

// Compute bond orientational parameter, based on the work of Steinhardt, P. J. et al. 1983, modified by Chua et al. 2010 and modified by me for accounting for multisite crystals 
double* ComputeAuxiliary::BondOrientationalParameter(){
	vector<string> DesProp;
	DesProp.push_back("STEINHARDT_MODE OneL");
	if( _MySystem->getIsCrystalDefined() ){
		if( _MySystem->getCrystal()->getIsMultisite() || _MySystem->getCrystal()->getIsReferenceBondOriParam() ){
			if( !_MySystem->getCrystal()->getIsReferenceBondOriParam() ){
				cerr << "The reference bond orientational parameters have not been computed for this crystal" << endl;
				cerr << "Have a look on the SaveNonCSCrystalBondOriParam executable to compute them" << endl;
				cerr << "Aborting" << endl;
				exit(EXIT_FAILURE);
			}
			if( !IsSteinhardtDescriptor ){
				StDes = new SteinhardtDescriptors(_MySystem,_MySystem->getCrystal()->getBondOriParamProperties());
				IsSteinhardtDescriptor = true;
			}
			ComputeAtomSiteIndex();
			for(unsigned int i=0;i<_MySystem->getNbAtom();i++){
				unsigned int current_type = _MySystem->getAtom(i).type_uint;
				StDes->getDescriptors()[i] /= _MySystem->getCrystal()->getReferenceBondOriParam()[current_type-1][Atom_SiteIndex[i]];
				if( StDes->getDescriptors()[i] > 1. ) StDes->getDescriptors()[i] = 1.;
				StDes->getDescriptors()[i] = 1. - StDes->getDescriptors()[i];
			}
			return StDes->getDescriptors();
		}
		double rcut = 0.; // set rcut to be the average of cell parameters length times security factor
		for(unsigned int i=0;i<3;i++) rcut += _MySystem->getCrystal()->getALength()[i];
		rcut /= 3.;
		rcut *= 1.25;
		DesProp.push_back("CUTOFF_RADIUS "+to_string(rcut));
	}
	if( !IsSteinhardtDescriptor ){
		StDes = new SteinhardtDescriptors(_MySystem,DesProp);
		IsSteinhardtDescriptor = true;
	}else{
		StDes->readProperties(DesProp);
		StDes->InitializeArrays();
		StDes->ComputeDescriptors();
	}
	for(unsigned int i=0;i<_MySystem->getNbAtom();i++)
		StDes->getDescriptors()[i] = 1. - StDes->getDescriptors()[i];
	return StDes->getDescriptors();
}

void ComputeAuxiliary::ComputeAtomSiteIndex(){
	if( _MySystem->getCrystal()->getIsMultisite() ){
		if( !_MySystem->getCrystal()->getIsReferenceBondOriParam() ){
			cerr << "The reference bond orientational parameters have not been computed for this crystal" << endl;
			cerr << "Have a look on the SaveNonCSCrystalBondOriParam executable to compute them" << endl;
			cerr << "Aborting" << endl;
			exit(EXIT_FAILURE);
		}
	}//else{
	//	cerr << "The crystal is not multisite, we then cannot compute atom site index, aborting" << endl;
	//	exit(EXIT_FAILURE);
	//}
	if( !IsSteinhardtDescriptor ){
		StDes = new SteinhardtDescriptors(_MySystem,_MySystem->getCrystal()->getBondOriParamProperties());
		IsSteinhardtDescriptor = true;
	}
	if( !IsAtomSiteIndex ){
		Atom_SiteIndex = new unsigned int[_MySystem->getNbAtom()];
		IsAtomSiteIndex = true;
	}
	unsigned int NbMaxSite = MT->max_p(_MySystem->getCrystal()->getNbAtomSite(),_MySystem->getCrystal()->getNbAtomType());
	double *diff2ref = new double[NbMaxSite];
// TODO here warning with type and element (no need if we do this in AtomicSystem::setCrystal()
	for(unsigned int i=0;i<_MySystem->getNbAtom();i++){
		unsigned int current_type = _MySystem->getAtom(i).type_uint;
		unsigned int current_nsite = _MySystem->getCrystal()->getNbAtomSite(current_type);
		for(unsigned int s=0;s<current_nsite;s++) diff2ref[s] = fabs(StDes->getDescriptors()[i]-_MySystem->getCrystal()->getReferenceBondOriParam()[current_type-1][s]);
		Atom_SiteIndex[i] = MT->min_p_ind(diff2ref,current_nsite);
	}
	delete[] diff2ref;
}

double* ComputeAuxiliary::Compute_StrainTensor(){
	if( !_MySystem->getIsCrystalDefined() ){
		cerr << "The crystal need to be defined for the strain tensor calculation" << endl;
		exit(EXIT_FAILURE);
	}
	if( !IsAtomSiteIndex ){
		cerr << "The bond orientationnal parameter need to have been computed before computing strain tensor" << endl;
		exit(EXIT_FAILURE);
	}
	const unsigned int nbAt = _MySystem->getNbAtom();
	unsigned int id;
	this->StrainTensor = new double[nbAt*6]; // the xx, yy, zz, xy, xz, yz strain component for each atom
	vector<vector<unsigned int>> RestrictedN; // neighbour list restricted to atom of same site
	vector<vector<double>> vec_RestrictedN; // vector and distance to the restricted neighbours 
	this->IsStrainTensor = true;
	double xp, yp ,zp, xpos, ypos, zpos;
	double Fac_HighestLength = 1.5;
	double rc = MT->max_p(_MySystem->getCrystal()->getALength(),3)*Fac_HighestLength;
	_MySystem->searchNeighbours(rc);
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	
	for(unsigned int i=0;i<nbAt;i++){
		RestrictedN.push_back(vector<unsigned int>());
		vec_RestrictedN.push_back(vector<double>());
		xpos = _MySystem->getWrappedPos(i).x;
		ypos = _MySystem->getWrappedPos(i).y;
		zpos = _MySystem->getWrappedPos(i).z;
		for(unsigned int j=0;j<_MySystem->getNeighbours(i*(nbNMax+1));j++){
			id = _MySystem->getNeighbours(i*(nbNMax+1)+j+1);
			if( Atom_SiteIndex[i] == Atom_SiteIndex[id] ){
				RestrictedN[i].push_back(id);
				xp = _MySystem->getWrappedPos(id).x+_MySystem->getCLNeighbours(i*nbNMax*3+j*3)*_MySystem->getH1()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+1)*_MySystem->getH2()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+2)*_MySystem->getH3()[0]-xpos;
				yp = _MySystem->getWrappedPos(id).y+_MySystem->getCLNeighbours(i*nbNMax*3+j*3)*_MySystem->getH1()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+1)*_MySystem->getH2()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+2)*_MySystem->getH3()[1]-ypos;
				zp = _MySystem->getWrappedPos(id).z+_MySystem->getCLNeighbours(i*nbNMax*3+j*3)*_MySystem->getH1()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+1)*_MySystem->getH2()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+2)*_MySystem->getH3()[2]-zpos;
				vec_RestrictedN[i].push_back(xp);
				vec_RestrictedN[i].push_back(yp);
				vec_RestrictedN[i].push_back(zp);
				vec_RestrictedN[i].push_back(sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)));
			}
		}
		if( RestrictedN.size() < 3 ){
			cerr << "We don't have found enough neighbours of the same site to compute the strain tensor" << endl;
			exit(EXIT_FAILURE);
		}
	}

	// search vectors closest to the cell parameters

	
	return this->StrainTensor;
}

double* ComputeAuxiliary::Compute_StrainTensor(unsigned int FromNum){
	if( FromNum == 1 ) cout << "Computing strain tensor based on atom numerotation" << endl;
	if( FromNum == 0 ){
		cout << "Computing strain tensor based on atom type" << endl;
		return Compute_StrainTensor();
	}
	const unsigned int nbAt = _MySystem->getNbAtom();
	this->StrainTensor = new double[nbAt*6]; // the xx, yy, zz, xy, xz, yz strain component for each atom
	vector<vector<unsigned int>> RestrictedN; // neighbour list restricted to atom of same site
	vector<vector<double>> vec_RestrictedN; // vector and distance to the restricted neighbours 
	this->IsStrainTensor = true;
	double xp, yp ,zp, xpos, ypos, zpos;
	double Fac_HighestLength = 1.5;
	double rc = MT->max_p(_MySystem->getCrystal()->getALength(),3)*Fac_HighestLength;
	unsigned int nbAt_Cryst = _MySystem->getCrystal()->getNbAtom();
	_MySystem->searchNeighbours(rc);
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	unsigned int id;
	if( FromNum == 1 ){	
		for(unsigned int i=0;i<nbAt;i++){
			RestrictedN.push_back(vector<unsigned int>());
			vec_RestrictedN.push_back(vector<double>());
			xpos = _MySystem->getWrappedPos(i).x;
			ypos = _MySystem->getWrappedPos(i).y;
			zpos = _MySystem->getWrappedPos(i).z;
			for(unsigned int j=0;j<_MySystem->getNeighbours(i*(nbNMax+1));j++){
				id = _MySystem->getNeighbours(i*(nbNMax+1)+j+1);
				if( abs((int) (i-id))%nbAt_Cryst == 0 ){
					RestrictedN[i].push_back(id); // TODO maybe no need of this array
					xp = _MySystem->getWrappedPos(id).x+_MySystem->getCLNeighbours(i*nbNMax*3+j*3)*_MySystem->getH1()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+1)*_MySystem->getH2()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+2)*_MySystem->getH3()[0]-xpos;
					yp = _MySystem->getWrappedPos(id).y+_MySystem->getCLNeighbours(i*nbNMax*3+j*3)*_MySystem->getH1()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+1)*_MySystem->getH2()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+2)*_MySystem->getH3()[1]-ypos;
					zp = _MySystem->getWrappedPos(id).z+_MySystem->getCLNeighbours(i*nbNMax*3+j*3)*_MySystem->getH1()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+1)*_MySystem->getH2()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+2)*_MySystem->getH3()[2]-zpos;
					vec_RestrictedN[i].push_back(xp);
					vec_RestrictedN[i].push_back(yp);
					vec_RestrictedN[i].push_back(zp);
					vec_RestrictedN[i].push_back(sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)));
				}
			}
			if( RestrictedN[i].size() < 3 ){
				cerr << "We don't have found enough neighbours of the same site to compute the strain tensor" << endl;
				exit(EXIT_FAILURE);
			}
		}
	}else{
		for(unsigned int i=0;i<nbAt;i++){
			RestrictedN.push_back(vector<unsigned int>());
			vec_RestrictedN.push_back(vector<double>());
			xpos = _MySystem->getWrappedPos(i).x;
			ypos = _MySystem->getWrappedPos(i).y;
			zpos = _MySystem->getWrappedPos(i).z;
			for(unsigned int j=0;j<_MySystem->getNeighbours(i*(nbNMax+1));j++){
				id = _MySystem->getNeighbours(i*(nbNMax+1)+j+1);
				if( _MySystem->getAtom(i).type_uint == _MySystem->getAtom(id).type_uint ){
					RestrictedN[i].push_back(id); // TODO maybe no need of this array
					xp = _MySystem->getWrappedPos(id).x+_MySystem->getCLNeighbours(i*nbNMax*3+j*3)*_MySystem->getH1()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+1)*_MySystem->getH2()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+2)*_MySystem->getH3()[0]-xpos;
					yp = _MySystem->getWrappedPos(id).y+_MySystem->getCLNeighbours(i*nbNMax*3+j*3)*_MySystem->getH1()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+1)*_MySystem->getH2()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+2)*_MySystem->getH3()[1]-ypos;
					zp = _MySystem->getWrappedPos(id).z+_MySystem->getCLNeighbours(i*nbNMax*3+j*3)*_MySystem->getH1()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+1)*_MySystem->getH2()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+2)*_MySystem->getH3()[2]-zpos;
					vec_RestrictedN[i].push_back(xp);
					vec_RestrictedN[i].push_back(yp);
					vec_RestrictedN[i].push_back(zp);
					vec_RestrictedN[i].push_back(sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)));
				}
			}
			if( RestrictedN[i].size() < 3 ){
				cerr << "We don't have found enough neighbours of the same type to compute the strain tensor" << endl;
				exit(EXIT_FAILURE);
			}
		}
	}

	// search vectors closest to cell vectors
	vector<double> Diff_CellParam;
	// finally we should have an array containing the three vectors closest to the cell vectors
	double *cell_vec = new double[9];
	double *rot_mat = new double[9];
	double *buffer_mat = new double[9];
	double *ref_crystal = new double[9];
	double *rot_ax = new double[3];
	unsigned int first_ind, second_ind, third_ind;
	double tol_ZeroScalarProd_init = 5e-2; //TODO test this value
	double tol_ZeroScalarProd;
	double sp, norm, theta;
	double *crossProd = new double[3];
	bool Found;
	for(unsigned int i=0;i<9;i++) ref_crystal[i] = 0.;
	for(unsigned int i=0;i<3;i++) ref_crystal[i*3+i] = _MySystem->getCrystal()->getALength()[i];
	MT->invert3x3(ref_crystal,ref_crystal);
	for(unsigned int i=0;i<nbAt;i++){
		Diff_CellParam.clear();
		for(unsigned int j=0;j<RestrictedN[i].size();j++){
			for(unsigned int k=0;k<3;k++){
				Diff_CellParam.push_back(fabs(vec_RestrictedN[i][j*4+3]-_MySystem->getCrystal()->getALength()[k]));
				Diff_CellParam.push_back(k);
				Diff_CellParam.push_back(j);
			}
		}
		MT->sort(Diff_CellParam,0,3,Diff_CellParam);
		// get the first vector (the one with the norm closest to one of cell vec)
		first_ind = round(Diff_CellParam[1]);
		norm = 0.;
		for(unsigned int dim=0;dim<3;dim++){
			cell_vec[(dim*3)+first_ind] = vec_RestrictedN[i][round(Diff_CellParam[2])*4+dim];
			norm += pow(cell_vec[(dim*3)+first_ind],2.);
		}
		norm = sqrt(norm);
		// get the second vector (the one with the norm closest to the second cell vec and close normal to the first vec)
		if( first_ind == 2 ) second_ind = 0;
		else second_ind = first_ind+1;
		Found = false;
		tol_ZeroScalarProd = tol_ZeroScalarProd_init;
		while( !Found ){
			for(unsigned int j=1;j<Diff_CellParam.size()/3;j++){
				if( round(Diff_CellParam[j*3+1]) == second_ind ){
					sp = 0.;
					for(unsigned int dim=0;dim<3;dim++) sp += cell_vec[(dim*3)+first_ind]*vec_RestrictedN[i][round(Diff_CellParam[j*3+2])*4+dim];
					sp /= norm*vec_RestrictedN[i][round(Diff_CellParam[j*3+2])*4+3];
					if( fabs(sp) < tol_ZeroScalarProd ){
						for(unsigned int dim=0;dim<3;dim++) cell_vec[(dim*3)+second_ind] = vec_RestrictedN[i][round(Diff_CellParam[j*3+2])*4+dim];
						Found = true;
						break;
					}
				}
			}
			if( Found ) break;
			else tol_ZeroScalarProd *= 2;
		}
		// get the third vector (the one with the norm closest to the third cell vec and close normal to the first and second vec and forming a direct basis)
		// in practise we search the vector closest to the third cell vector with thei angle closest to the cross product of the first two
		crossProd[0] = cell_vec[3+first_ind]*cell_vec[6+second_ind] - cell_vec[6+first_ind]*cell_vec[3+second_ind];
		crossProd[1] = cell_vec[6+first_ind]*cell_vec[0+second_ind] - cell_vec[0+first_ind]*cell_vec[6+second_ind];
		crossProd[2] = cell_vec[0+first_ind]*cell_vec[3+second_ind] - cell_vec[3+first_ind]*cell_vec[0+second_ind];
		norm = 0.;
		for(unsigned int dim=0;dim<3;dim++) norm += pow(crossProd[dim],2.);
		norm = sqrt(norm);
		if( second_ind == 2 ) third_ind = 0;
		else third_ind = second_ind+1;
		Found = false;
		tol_ZeroScalarProd = tol_ZeroScalarProd_init;
		while( !Found ){
			for(unsigned int j=1;j<Diff_CellParam.size()/3;j++){
				if( round(Diff_CellParam[j*3+1]) == third_ind ){
					sp = 0.;
					for(unsigned int dim=0;dim<3;dim++) sp += crossProd[dim]*vec_RestrictedN[i][round(Diff_CellParam[j*3+2])*4+dim];
					sp /= norm*vec_RestrictedN[i][round(Diff_CellParam[j*3+2])*4+3];
					if( fabs(sp-1.) < tol_ZeroScalarProd ){
						for(unsigned int dim=0;dim<3;dim++) cell_vec[(dim*3)+third_ind] = vec_RestrictedN[i][round(Diff_CellParam[j*3+2])*4+dim];
						Found = true;
						break;
					}
				}
			}
			if( Found ) break;
			else tol_ZeroScalarProd *= 2;
		}
		// rotate the cell vectors for first vector to be aligned with x axis and second vector to be in (x,y) plane
		// first rot around z axis
		if( i == 1 ){
			cout << "break" << endl;
		}
		norm = sqrt(pow(cell_vec[3],2.)+pow(cell_vec[0],2.));
		if( norm > 1e-6 ){
			theta = asin(cell_vec[3]/norm);
			if( cell_vec[0] < 0 ){
				if( cell_vec[3] < 0 ) theta = -theta-M_PI;
				else theta = M_PI-theta;
			}
			rot_ax[0] = 0.;
			rot_ax[1] = 0.;
			rot_ax[2] = 1.;
			MT->Vec2rotMat(rot_ax,-theta,rot_mat);
			MT->MatDotMat(rot_mat,cell_vec,cell_vec);
		}
		// second rot around y axis
		norm = sqrt(pow(cell_vec[6],2.)+pow(cell_vec[0],2.));
		if( norm > 1e-6 ){
			theta = asin(cell_vec[6]/norm);
			if( cell_vec[0] < 0 ){
				if( cell_vec[6] < 0 ) theta = -theta-M_PI;
				else theta = M_PI-theta;
			}
			rot_ax[0] = 0.;
			rot_ax[1] = 1.;
			rot_ax[2] = 0.;
			MT->Vec2rotMat(rot_ax,theta,rot_mat);
			MT->MatDotMat(rot_mat,cell_vec,cell_vec);
		}
		// last rot about x axis
		norm = sqrt(pow(cell_vec[7],2.)+pow(cell_vec[4],2.));
		if( norm > 1e-6 ){
			theta = asin(cell_vec[7]/norm);
			if( cell_vec[4] < 0 ){
				if( cell_vec[7] < 0 ) theta = -theta-M_PI;
				else theta = M_PI-theta;
			}
			rot_ax[0] = 1.;
			rot_ax[1] = 0.;
			rot_ax[2] = 0.;
			MT->Vec2rotMat(rot_ax,-theta,rot_mat);
			MT->MatDotMat(rot_mat,cell_vec,cell_vec);
		}
		MT->MatDotMat(ref_crystal,cell_vec,buffer_mat);
		this->StrainTensor[i*6+0] = buffer_mat[0];
		this->StrainTensor[i*6+1] = buffer_mat[4];
		this->StrainTensor[i*6+2] = buffer_mat[8];
		this->StrainTensor[i*6+3] = buffer_mat[1];
		this->StrainTensor[i*6+4] = buffer_mat[2];
		this->StrainTensor[i*6+5] = buffer_mat[5];
	}
	delete[] cell_vec;
	delete[] rot_mat;
	delete[] buffer_mat;
	delete[] ref_crystal;
	delete[] rot_ax;
	delete[] crossProd;
	return this->StrainTensor;
}	
	
double* ComputeAuxiliary::Compute_StrainTensor_invII(){
	if( !IsStrainTensor ) Compute_StrainTensor(0);
	const unsigned int nbAt = _MySystem->getNbAtom();
	this->Strain_invII = new double[nbAt];
	this->IsStrainInvII = true;
	for(unsigned int i=0;i<nbAt;i++) this->Strain_invII[i] = this->StrainTensor[i*6+0]*this->StrainTensor[i*6+1] + this->StrainTensor[i*6+1]*this->StrainTensor[i*6+2] + this->StrainTensor[i*6+0]*this->StrainTensor[i*6+2] + this->StrainTensor[i*6+3]*this->StrainTensor[i*6+3] + this->StrainTensor[i*6+4]*this->StrainTensor[i*6+4] + this->StrainTensor[i*6+5]*this->StrainTensor[i*6+5];
       return this->Strain_invII;	
}


void ComputeAuxiliary::read_params(){ // nothing to read for the moment
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+FixedParam_Filename;
	ifstream file(filename, ios::in);
	size_t pos_tolSite;
	string buffer_s, line;
	if(file){
		while(file){
			getline(file,line);
			pos_tolSite=line.find("TOL_SITES ");
			if(pos_tolSite!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->tolSites;
			}
		}
	}
	//else{
	//	cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
	//	exit(EXIT_FAILURE);
	//}
}

// atomic strain as defined in ovito (i.e. Shimizu, Ogata, Li: Mater. Trans. 48 (2007), 2923)
double* ComputeAuxiliary::Compute_AffineTransfoMatrix(AtomicSystem &AnalyzedSystem, double rc){
	if( !_MySystem->getIsNeighbours() || _MySystem->get_current_rc() != rc ){
		_MySystem->searchNeighbours(rc);
	}

	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	const unsigned int nbAt_ref = AnalyzedSystem.getNbAtom();

	double xpos, ypos, zpos, xpos_r, ypos_r, zpos_r, id;
	double *buffer_mat = new double[9];

	if( nbAt != nbAt_ref ){
		cout << "The number of atom in the reference system and system to compute atomic strain is not the same !" << endl;
		return 0;
	}
	
	// Compute d0 and Vi if not already computed
	if( !this->Reference_AtomicStrain_Computed ){
		this->Vi_inv = new double[nbAt*9];
		this->d0 = new double[nbAt*nbNMax*3];
		
		for(unsigned int i=0;i<nbAt;i++){
			xpos_r = _MySystem->getWrappedPos(i).x;
			ypos_r = _MySystem->getWrappedPos(i).y;
			zpos_r = _MySystem->getWrappedPos(i).z;
			for(unsigned int j=0;j<9;j++) buffer_mat[j] = 0.;
			for(unsigned int n=0;n<_MySystem->getNeighbours(i*(nbNMax+1));n++){
				id = _MySystem->getNeighbours(i*(nbNMax+1)+n+1);
				this->d0[i*nbNMax*3+n*3] = _MySystem->getWrappedPos(id).x+_MySystem->getCLNeighbours(i*nbNMax*3+n*3)*_MySystem->getH1()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+1)*_MySystem->getH2()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+2)*_MySystem->getH3()[0]-xpos_r;
				this->d0[i*nbNMax*3+n*3+1] = _MySystem->getWrappedPos(id).y+_MySystem->getCLNeighbours(i*nbNMax*3+n*3)*_MySystem->getH1()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+1)*_MySystem->getH2()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+2)*_MySystem->getH3()[1]-ypos_r;
				this->d0[i*nbNMax*3+n*3+2] = _MySystem->getWrappedPos(id).z+_MySystem->getCLNeighbours(i*nbNMax*3+n*3)*_MySystem->getH1()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+1)*_MySystem->getH2()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+2)*_MySystem->getH3()[2]-zpos_r;
				for(unsigned int d1=0;d1<3;d1++){
					for(unsigned int d2=0;d2<3;d2++) buffer_mat[d1*3+d2] += this->d0[i*nbNMax*3+n*3+d1]*this->d0[i*nbNMax*3+n*3+d2];
				}
			}
			MT->invert3x3(buffer_mat,buffer_mat);
			for(unsigned int d=0;d<9;d++) Vi_inv[i*9+d] = buffer_mat[d];
		}
		this->Reference_AtomicStrain_Computed = true;
	}

	double *W = new double[9];
	double *buffer_vec = new double[3];
	double *delta = new double[nbAt*3];
	if( !this->IsJi ){
		this->Ji = new double[nbAt*9];
		this->current_d = new double[nbAt*nbNMax*3];
		this->IsJi = true;
	}

	for(unsigned int i=0;i<nbAt;i++){
		for(unsigned int d=0;d<3;d++) delta[i*3+d] = 0.;
		xpos = AnalyzedSystem.getWrappedPos(i).x;
		ypos = AnalyzedSystem.getWrappedPos(i).y;
		zpos = AnalyzedSystem.getWrappedPos(i).z;
		xpos_r = _MySystem->getWrappedPos(i).x;
		ypos_r = _MySystem->getWrappedPos(i).y;
		zpos_r = _MySystem->getWrappedPos(i).z;
		if( (xpos-xpos_r) > AnalyzedSystem.getH1()[0]/2. ) for(unsigned int d=0;d<3;d++) delta[i*3+d] -= AnalyzedSystem.getH1()[d];
		else if( (xpos-xpos_r) < -AnalyzedSystem.getH1()[0]/2. ) for(unsigned int d=0;d<3;d++) delta[i*3+d] += AnalyzedSystem.getH1()[d];
		if( (ypos-ypos_r) > AnalyzedSystem.getH2()[1]/2. ) for(unsigned int d=0;d<3;d++) delta[i*3+d] -= AnalyzedSystem.getH2()[d];
		else if( (ypos-ypos_r) < -AnalyzedSystem.getH2()[1]/2. ) for(unsigned int d=0;d<3;d++) delta[i*3+d] += AnalyzedSystem.getH2()[d];
		if( (zpos-zpos_r) > AnalyzedSystem.getH3()[2]/2. ) for(unsigned int d=0;d<3;d++) delta[i*3+d] -= AnalyzedSystem.getH3()[d];
		else if( (zpos-zpos_r) < -AnalyzedSystem.getH3()[2]/2. ) for(unsigned int d=0;d<3;d++) delta[i*3+d] += AnalyzedSystem.getH3()[d];
	}

	for(unsigned int i=0;i<nbAt;i++){
		xpos = AnalyzedSystem.getWrappedPos(i).x;
		ypos = AnalyzedSystem.getWrappedPos(i).y;
		zpos = AnalyzedSystem.getWrappedPos(i).z;
		for(unsigned int j=0;j<9;j++){
			W[j] = 0.;
			buffer_mat[j] = this->Vi_inv[i*9+j];
		}
		// Compute d and W
		for(unsigned int n=0;n<_MySystem->getNeighbours(i*(nbNMax+1));n++){
			id = _MySystem->getNeighbours(i*(nbNMax+1)+n+1);
			 this->current_d[i*nbNMax*3+n*3] = AnalyzedSystem.getWrappedPos(id).x+_MySystem->getCLNeighbours(i*nbNMax*3+n*3)*AnalyzedSystem.getH1()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+1)*AnalyzedSystem.getH2()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+2)*AnalyzedSystem.getH3()[0]-xpos-delta[i*3];
			 this->current_d[i*nbNMax*3+n*3+1] = AnalyzedSystem.getWrappedPos(id).y+_MySystem->getCLNeighbours(i*nbNMax*3+n*3)*AnalyzedSystem.getH1()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+1)*AnalyzedSystem.getH2()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+2)*AnalyzedSystem.getH3()[1]-ypos-delta[i*3+1];
			 this->current_d[i*nbNMax*3+n*3+2] = AnalyzedSystem.getWrappedPos(id).z+_MySystem->getCLNeighbours(i*nbNMax*3+n*3)*AnalyzedSystem.getH1()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+1)*AnalyzedSystem.getH2()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+2)*AnalyzedSystem.getH3()[2]-zpos-delta[i*3+2];
			if( this->current_d[i*nbNMax*3+n*3] > AnalyzedSystem.getH1()[0]/2. ) for(unsigned int d=0;d<3;d++) this->current_d[i*nbNMax*3+n*3+d] -= AnalyzedSystem.getH1()[d];
			else if( this->current_d[i*nbNMax*3+n*3] < -AnalyzedSystem.getH1()[0]/2. ) for(unsigned int d=0;d<3;d++) this->current_d[i*nbNMax*3+n*3+d] += AnalyzedSystem.getH1()[d];
			if( this->current_d[i*nbNMax*3+n*3+1] > AnalyzedSystem.getH2()[1]/2. ) for(unsigned int d=0;d<3;d++) this->current_d[i*nbNMax*3+n*3+d] -= AnalyzedSystem.getH2()[d];
			else if( this->current_d[i*nbNMax*3+n*3+1] < -AnalyzedSystem.getH2()[1]/2. ) for(unsigned int d=0;d<3;d++) this->current_d[i*nbNMax*3+n*3+d] += AnalyzedSystem.getH2()[d];
			if( this->current_d[i*nbNMax*3+n*3+2] > AnalyzedSystem.getH3()[2]/2. ) for(unsigned int d=0;d<3;d++) this->current_d[i*nbNMax*3+n*3+d] -= AnalyzedSystem.getH3()[d];
			else if( this->current_d[i*nbNMax*3+n*3+2] < -AnalyzedSystem.getH3()[2]/2. ) for(unsigned int d=0;d<3;d++) this->current_d[i*nbNMax*3+n*3+d] += AnalyzedSystem.getH3()[d];
			for(unsigned int d1=0;d1<3;d1++){
				for(unsigned int d2=0;d2<3;d2++) W[d1*3+d2] += this->d0[i*nbNMax*3+n*3+d1]*(this->current_d[i*nbNMax*3+n*3+d2]);
			}
		}
		// Compute deformation gradient tensor
		MT->MatDotMat(buffer_mat,W,W);
		for(unsigned int d=0;d<9;d++) this->Ji[i*9+d] = W[d];
	}
	
	delete[] W;
	delete[] buffer_vec;
	delete[] buffer_mat;
	delete[] delta;
	
	return this->Ji;
}
		
double* ComputeAuxiliary::Compute_AtomicStrain(AtomicSystem &AnalyzedSystem, double rc){

	Compute_AffineTransfoMatrix(AnalyzedSystem,rc);
	
	const unsigned int nbAt = _MySystem->getNbAtom();
	
	double *buffer_mat = new double[9];
	double *buffer_mat_1 = new double[9];
	
	if( !this->IsAtomicStrain ){
		this->AtomicStrain = new double[nbAt*8]; // sorted as eps_xx,eps_yy,eps_zz,_eps_xy,eps_xz,eps_yz,shear_inv,hydro_inv
		this->IsAtomicStrain = true;
	}

	for(unsigned int i=0;i<nbAt;i++){
		for(unsigned int d=0;d<9;d++) buffer_mat_1[d] = this->Ji[i*9+d];
		// Compute Lagrangian strain tensor (.5*(J.Jt-Id))	
		MT->MatDotTMat(buffer_mat_1,buffer_mat_1,buffer_mat);

		for(unsigned int d=0;d<3;d++){
			buffer_mat[d*3+d] -= 1.;
			buffer_mat[d*3+d] *= .5; //originally all strain component must be divided by two but regarding analytical tests we don't do that here
		}
		// sort the tensor 
		AtomicStrain[i*8] = buffer_mat[0];
		AtomicStrain[i*8+1] = buffer_mat[4];
		AtomicStrain[i*8+2] = buffer_mat[8];
		AtomicStrain[i*8+3] = buffer_mat[1];
		AtomicStrain[i*8+4] = buffer_mat[2];
		AtomicStrain[i*8+5] = buffer_mat[5];
		// compute shear invariant
		AtomicStrain[i*8+6] = 0.;
		for(unsigned int d1=0;d1<3;d1++) AtomicStrain[i*8+6] += pow(AtomicStrain[i*8+3+d1],2.);
		AtomicStrain[i*8+6] += (pow(buffer_mat[4]-buffer_mat[8],2.)+pow(buffer_mat[0]-buffer_mat[8],2.)+pow(buffer_mat[0]-buffer_mat[4],2.))/8.;
		AtomicStrain[i*8+6] = sqrt(AtomicStrain[i*8+6]);
		// compute hydrostatic invariant

		AtomicStrain[i*8+7] = 0.;
		for(unsigned int d1=0;d1<3;d1++) AtomicStrain[i*8+7] += buffer_mat[d1*3+d1];
		AtomicStrain[i*8+7] /= 3.;
	}
	
	delete[] buffer_mat;
	delete[] buffer_mat_1;

	return AtomicStrain;
}

// D2Min as defined in Delbecq et al. 2023
double* ComputeAuxiliary::Compute_D2Min(AtomicSystem &AnalyzedSystem, double rc){ 

	Compute_AffineTransfoMatrix(AnalyzedSystem,rc);
	
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	if( !this->IsD2Min ){
		this->D2Min = new double[nbAt];
		IsD2Min = true;
	}

	double *buffer_vec = new double[3];
	double *buffer_mat = new double[9];

	for(unsigned int i=0;i<nbAt;i++){
		this->D2Min[i] = 0.;
		for(unsigned int d=0;d<9;d++) buffer_mat[d] = this->Ji[i*9+d];
		for(unsigned int n=0;n<_MySystem->getNeighbours(i*(nbNMax+1));n++){
			// Compute Ji*delta0
			for(unsigned int d1=0;d1<3;d1++) buffer_vec[d1] = this->d0[i*nbNMax*3+n*3+d1];
			MT->MatDotRawVec(buffer_mat,buffer_vec,buffer_vec);
			// Compute sum( (delta - Ji*delta0)^2
			for(unsigned int d=0;d<3;d++) this->D2Min[i] += pow(this->current_d[i*nbNMax*3+n*3+d]-buffer_vec[d],2.);
		}
		// Normalize D2Min with respect to the number of neighbours
		// this->D2Min[i] /= _MySystem->getNeighbours(i*(nbNMax+1));
	}

	delete[] buffer_vec;
	delete[] buffer_mat;

	return this->D2Min;
}

bool ComputeAuxiliary::AreSurfacesSame(AtomicSystem *AtSys1, std::string updown1, AtomicSystem *AtSys2, std::string updown2, bool symm){
	if( ( updown1 != "down" && updown1 != "up" ) || ( updown2 != "down" && updown2 != "up" ) ){
		cout << "The surface to symmetrize should be either \"down\" or \"up\"" << endl;
		return false;		
	}
	
	bool reverse = false;
	if( ( updown1 == "down" && updown2 == "up" ) || ( updown1 == "up" && updown2 == "down" ) ) reverse = true;

	double slab_width = 10.;
	// search min and max z pos of the system
	double min_z1 = std::numeric_limits<double>::max();
	double max_z1 = -std::numeric_limits<double>::max();
	double min_z2 = std::numeric_limits<double>::max();
	double max_z2 = -std::numeric_limits<double>::max();
	double *refpos1 = new double[3];
	double *refpos2 = new double[3];
	unsigned int nbAtom1 = AtSys1->getNbAtom();
	unsigned int nbAtom2 = AtSys2->getNbAtom();
	for(unsigned int i=0;i<nbAtom1;i++){
		if( AtSys1->getAtom(i).pos.z < min_z1 ){
			min_z1 = AtSys1->getAtom(i).pos.z;
			if( updown1 == "down" ){
				refpos1[0] = AtSys1->getAtom(i).pos.x;
				refpos1[1] = AtSys1->getAtom(i).pos.y;
				refpos1[2] = AtSys1->getAtom(i).pos.z;
			}
		}
		if( AtSys1->getAtom(i).pos.z > max_z1 ){
			max_z1 = AtSys1->getAtom(i).pos.z;
			if( updown1 == "up" ){
				refpos1[0] = AtSys1->getAtom(i).pos.x;
				refpos1[1] = AtSys1->getAtom(i).pos.y;
				refpos1[2] = AtSys1->getAtom(i).pos.z;
			}
		}
	}
	for(unsigned int i=0;i<nbAtom2;i++){
		if( AtSys2->getAtom(i).pos.z < min_z2 ){
			min_z2 = AtSys2->getAtom(i).pos.z;
			if( updown2 == "down" ){
				refpos2[0] = AtSys2->getAtom(i).pos.x;
				refpos2[1] = AtSys2->getAtom(i).pos.y;
				refpos2[2] = AtSys2->getAtom(i).pos.z;
			}
		}
		if( AtSys2->getAtom(i).pos.z > max_z2 ){
			max_z2 = AtSys2->getAtom(i).pos.z;
			if( updown2 == "up" ){
				refpos2[0] = AtSys2->getAtom(i).pos.x;
				refpos2[1] = AtSys2->getAtom(i).pos.y;
				refpos2[2] = AtSys2->getAtom(i).pos.z;
			}
		}
	}

	vector<Atom> SubSys1, SubSys2;
	for(unsigned int i=0;i<nbAtom1;i++){
		if( updown1 == "down" && AtSys1->getAtom(i).pos.z < min_z1+slab_width ) SubSys1.push_back(AtSys1->getAtom(i));
		else if( updown1 == "up" && AtSys1->getAtom(i).pos.z > max_z1-slab_width ) SubSys1.push_back(AtSys1->getAtom(i));
	}
	for(unsigned int i=0;i<nbAtom2;i++){
		if( updown2 == "down" && AtSys2->getAtom(i).pos.z < min_z2+slab_width ) SubSys2.push_back(AtSys2->getAtom(i));
		else if( updown2 == "up" && AtSys2->getAtom(i).pos.z > max_z2-slab_width ) SubSys2.push_back(AtSys2->getAtom(i));
	}
	
	double minmax;
	if( reverse ){
		for(unsigned int i=0;i<SubSys1.size();i++) SubSys1[i].pos.z *= -1.;
		if( !symm ){ 
			double *vec = new double[3];
			vec[0] = 0.;
			vec[1] = 0.;
			vec[2] = 1.;
			double theta = M_PI;
			double *rotmat = new double[9];
			MT->Vec2rotMat(vec,theta,rotmat);
			for(unsigned int i=0;i<SubSys1.size();i++) MT->MatDotAt(rotmat,SubSys1[i],SubSys1[i]);
			delete[] vec;
			delete[] rotmat;
		}
		// recompute refpos1
		if( updown1 == "down" ) minmax = -std::numeric_limits<double>::max();
		else minmax = std::numeric_limits<double>::max();
		for(unsigned int i=0;i<SubSys1.size();i++){
			if( ( updown1 == "down" && SubSys1[i].pos.z > minmax ) || ( updown1 == "up" && SubSys1[i].pos.z < minmax ) ){
				minmax = SubSys1[i].pos.z;
				refpos1[0] = SubSys1[i].pos.x;
				refpos1[1] = SubSys1[i].pos.y;
				refpos1[2] = SubSys1[i].pos.z;
			}
		}
	}else if( updown1 == "up" ) minmax = max_z1;
	else minmax = min_z1;

	bool debug = false;
	if( debug ){
		Atom *AtList = new Atom[SubSys1.size()+SubSys2.size()];
		double *aux2print = new double[SubSys1.size()+SubSys2.size()];
		for(unsigned int i=0;i<SubSys1.size();i++){
			AtList[i] = SubSys1[i];
			aux2print[i] = 0.;
		}
		for(unsigned int i=0;i<SubSys2.size();i++){
			AtList[i+SubSys1.size()] = SubSys2[i];
			aux2print[i+SubSys1.size()] = 1.;
		}
		Crystal *MyC = new Crystal("Forsterite");
		AtomicSystem *ToPrint = new AtomicSystem(AtList,SubSys1.size()+SubSys2.size(),MyC,AtSys1->getH1(),AtSys1->getH2(),AtSys1->getH3());
		ToPrint->setAux(aux2print,"aux");
		ToPrint->printSystem_aux("Debug.cfg","aux");
	}
	
	// shift sys1 to be in contact with sys2
	for(unsigned int i=0;i<SubSys1.size();i++){
		SubSys1[i].pos.x += refpos2[0] - refpos1[0];
		SubSys1[i].pos.y += refpos2[1] - refpos1[1];
		SubSys1[i].pos.z += refpos2[2] - refpos1[2];
	}
	delete[] refpos1;
	delete[] refpos2;

	// search if 1 coincide with 2 with bc
	for(unsigned int i=0;i<SubSys1.size();i++){
		if( ( ( updown1 == "down" && reverse ) || (updown1 == "up" && !reverse ) ) && SubSys1[i].pos.z < (minmax-(0.75*slab_width)) ) continue;
		else if( ( ( updown1 == "up" && reverse ) || (updown1 == "down" && !reverse ) ) && SubSys1[i].pos.z > (minmax+(0.75*slab_width)) ) continue;
		else{
			double sigma = 0.25;
			double min_sp = 1.e-1;
			double cur_at_sp = 0.;
			double xr = SubSys1[i].pos.x;
			double yr = SubSys1[i].pos.y;
			double zr = SubSys1[i].pos.z;
			for(unsigned int j=0;j<SubSys2.size();j++){
				for(int bx=-1;bx<2;bx++){
					for(int by=-1;by<2;by++){
						double xw = SubSys2[j].pos.x + AtSys2->getH1()[0]*bx + AtSys2->getH2()[0]*by;
						double yw = SubSys2[j].pos.y + AtSys2->getH1()[1]*bx + AtSys2->getH2()[1]*by;
						double zw = SubSys2[j].pos.z;
						cur_at_sp += MT->gaussian(xw, yw, zw, xr, yr, zr, sigma);
					}
				}
			}
			if( cur_at_sp < min_sp ) return false;
		}
	}
	// search if 2 coincide with 1 with bc
	for(unsigned int i=0;i<SubSys2.size();i++){
		if( updown2 == "up" && SubSys2[i].pos.z < (max_z2-(0.75*slab_width)) ) continue;
		else if( updown1 == "down" && SubSys2[i].pos.z > (min_z2+(0.75*slab_width)) ) continue;
		else{
			double sigma = 0.25;
			double min_sp = 1.e-1;
			double cur_at_sp = 0.;
			double xr = SubSys2[i].pos.x;
			double yr = SubSys2[i].pos.y;
			double zr = SubSys2[i].pos.z;
			for(unsigned int j=0;j<SubSys1.size();j++){
				for(int bx=-1;bx<2;bx++){
					for(int by=-1;by<2;by++){
						double xw = SubSys1[j].pos.x + AtSys1->getH1()[0]*bx + AtSys1->getH2()[0]*by;
						double yw = SubSys1[j].pos.y + AtSys1->getH1()[1]*bx + AtSys1->getH2()[1]*by;
						double zw = SubSys1[j].pos.z;
						cur_at_sp += MT->gaussian(xw, yw, zw, xr, yr, zr, sigma);
					}
				}
			}
			if( cur_at_sp < min_sp ) return false;
		}
	}
	return true;
}


ComputeAuxiliary::~ComputeAuxiliary(){
	delete MT;
	if( this->Reference_AtomicStrain_Computed ){
		delete[] Vi_inv;
		delete[] d0;
		delete[] current_d;
	}
	if( this->IsJi ) delete[] Ji;
	if( this->IsAtomicStrain ) delete[] AtomicStrain;
	if( this->IsD2Min ) delete[] D2Min;
	if( IsStrainTensor ) delete[] StrainTensor;
	if( IsStrainInvII ) delete[] Strain_invII;
	if( IsSteinhardtDescriptor ){
		//StDes->~SteinhardtDescriptor();
		delete StDes;
	}
	if( IsAtomSiteIndex ) delete[] Atom_SiteIndex;
}

