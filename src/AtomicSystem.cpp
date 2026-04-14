//**********************************************************************************
//*   AtomicSystem.cpp                                                             *
//**********************************************************************************
//* This file contains the implementation of the AtomicSystem class                *
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

#include "AtomicSystem.h"
#include "AtomHicConfig.h"
#include "Descriptors.h"
#include "KMeans.h"
#include "Crystal.h"
#include <cmath>
#include <sstream>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <iomanip>
#include <chrono>
#include <random>

using namespace std;

AtomicSystem::AtomicSystem(){
	Dis = new Displays;
}

AtomicSystem::AtomicSystem(Crystal *_MyCrystal, double xhi, double yhi, double zhi, vector<int> cl_box):_MyCrystal(_MyCrystal){
	Dis = new Displays;
	read_params_atsys();
	this->IsCrystalDefined = true;
	// Compute the number of atom and verify it gives an integer number
	double nbAt_d = _MyCrystal->getNbAtom()*xhi*yhi*zhi/_MyCrystal->getVol();
	if( fabs(nbAt_d-round(nbAt_d)) > 1e-3 ) cout << "Warning the cell dimension does not correspond to integer number of atom" << endl;
	this->nbAtom = round(nbAt_d);
	// initialize pointers
	this->AtomList = new Atom[this->nbAtom];
	this->H1 = new double[3]; 
	this->H2 = new double[3]; 
	this->H3 = new double[3]; 
	this->H1[0] = xhi;
	this->H1[1] = 0.;
	this->H1[2] = 0.;
	this->H2[0] = 0.;
	this->H2[1] = yhi;
	this->H2[2] = 0.;
	this->H3[0] = 0.;
	this->H3[1] = 0.;
	this->H3[2] = zhi;
	this->MT = new MathTools;
	this->nbAtomType = _MyCrystal->getNbAtomType();
	this->AtomType = _MyCrystal->getAtomType();
	this->AtomMass = _MyCrystal->getAtomMass();
	this->AtomCharge = _MyCrystal->getAtomCharge();
	this->IsCharge = _MyCrystal->getIsCharge();
	this->IsMolId = _MyCrystal->getIsMolId();
	if( IsMolId ) MolId = new unsigned int[nbAtom];
	this->IsBond = _MyCrystal->getIsBond();
	if( this->IsBond ){
		this->nbBondType = _MyCrystal->getNbBondType();
		this->nbBonds = round(_MyCrystal->getNbBond()*xhi*yhi*zhi/_MyCrystal->getVol());
		Bonds = new unsigned int[this->nbBonds*2];
		BondType = new unsigned int[this->nbBonds];
	}
	this->IsAngle = _MyCrystal->getIsAngle();
	if( this->IsAngle ){
		this->nbAngleType = _MyCrystal->getNbAngleType();
		this->nbAngles = round(_MyCrystal->getNbAngle()*xhi*yhi*zhi/_MyCrystal->getVol());
		Angles = new unsigned int[this->nbAngles*3];
		AngleType = new unsigned int[this->nbAngles];
	}
	this->IsTilted = false;
	this->G1 = new double[3];
	this->G2 = new double[3];
	this->G3 = new double[3];
	this->IsG = true;
	computeInverseCellVec();
	// compute the atomic positions by testing linear combination of unit cell crystal
	double cla1, cla2, cla3;
       	int arr[3];
	if( _MyCrystal->getCrystallo() == "Cubic" || _MyCrystal->getCrystallo() == "Orthorhombic" || _MyCrystal->getCrystallo() == "Tetragonal" ){
		if( fabs(_MyCrystal->getA1()[0]) >1e-3 ) arr[0] = round(xhi/fabs(_MyCrystal->getA1()[0]));
		else arr[0] = 0;
		if( fabs(_MyCrystal->getA1()[1]) >1e-3 ) arr[1] = round(yhi/fabs(_MyCrystal->getA1()[1]));
		else arr[1] = 0;
		if( fabs(_MyCrystal->getA1()[2]) >1e-3 ) arr[2] = round(zhi/fabs(_MyCrystal->getA1()[2]));
		else arr[2] = 0;
		cla1 = MT->max(arr,3)+2;
		if( fabs(_MyCrystal->getA2()[0]) >1e-3 ) arr[0] = round(xhi/fabs(_MyCrystal->getA2()[0]));
		else arr[0] = 0;
		if( fabs(_MyCrystal->getA2()[1]) >1e-3 ) arr[1] = round(yhi/fabs(_MyCrystal->getA2()[1]));
		else arr[1] = 0;
		if( fabs(_MyCrystal->getA2()[2]) >1e-3 ) arr[2] = round(zhi/fabs(_MyCrystal->getA2()[2]));
		else arr[2] = 0;
		cla2 = MT->max(arr,3)+2;
		if( fabs(_MyCrystal->getA3()[0]) >1e-3 ) arr[0] = round(xhi/fabs(_MyCrystal->getA3()[0]));
		else arr[0] = 0;
		if( fabs(_MyCrystal->getA3()[1]) >1e-3 ) arr[1] = round(yhi/fabs(_MyCrystal->getA3()[1]));
		else arr[1] = 0;
		if( fabs(_MyCrystal->getA3()[2]) >1e-3 ) arr[2] = round(zhi/fabs(_MyCrystal->getA3()[2]));
		else arr[2] = 0;
		cla3 = MT->max(arr,3)+2;
	}else{
		arr[0] = round(xhi/fabs(_MyCrystal->getA1()[0] + _MyCrystal->getA2()[0] + _MyCrystal->getA3()[0]));
		arr[1] = round(yhi/fabs(_MyCrystal->getA1()[1] + _MyCrystal->getA2()[1] + _MyCrystal->getA3()[1]));
		arr[2] = round(zhi/fabs(_MyCrystal->getA1()[2] + _MyCrystal->getA2()[2] + _MyCrystal->getA3()[2]));
		cla1 = MT->max(arr,3)+2;
		cla2 = cla1;
		cla3 = cla1;
	}

	double delta_x, delta_y, delta_z, xpos, ypos, zpos;
	unsigned int count;
	double tolIonPos = 1e-9;
	bool break_comp = false;
	bool box_fill = false;
	unsigned int count_box_fill = 0;
	unsigned int max_count_box_fill = 10;
	while( !box_fill && count_box_fill < max_count_box_fill ){
		count = 0;
		if( !_MyCrystal->getIsDoNotSep() ){
			// crystals where there is no constraints on atom separation
			for(int i=-cla1;i<cla1+1;i++){
				for(int j=-cla2;j<cla2+1;j++){
					for(int k=-cla3;k<cla3+1;k++){
						delta_x = i*_MyCrystal->getA1()[0]+j*_MyCrystal->getA2()[0]+k*_MyCrystal->getA3()[0];
						delta_y = i*_MyCrystal->getA1()[1]+j*_MyCrystal->getA2()[1]+k*_MyCrystal->getA3()[1];
						delta_z = i*_MyCrystal->getA1()[2]+j*_MyCrystal->getA2()[2]+k*_MyCrystal->getA3()[2];
						for(unsigned int n=0;n<_MyCrystal->getNbAtom();n++){
							xpos = _MyCrystal->getMotif()[n].pos.x + delta_x;
							ypos = _MyCrystal->getMotif()[n].pos.y + delta_y;
							zpos = _MyCrystal->getMotif()[n].pos.z + delta_z;
							if( xpos > -tolIonPos && xpos < xhi-tolIonPos && ypos > -tolIonPos && ypos < yhi-tolIonPos && zpos > -tolIonPos && zpos < zhi-tolIonPos ){
								if( count == this->nbAtom ){
									break_comp = true;
									break;
								}
								this->AtomList[count] = _MyCrystal->getMotif()[n];
								this->AtomList[count].pos.x = xpos;
								this->AtomList[count].pos.y = ypos;
								this->AtomList[count].pos.z = zpos;
								count += 1;
							}
						}
						if( break_comp) break;
					}
					if( break_comp) break;
				}
				if( break_comp) break;
			}
		}else{
			// crystals for which some atoms should not be separeted
			// create array of already stored atoms
			vector<vector<double>> AlreadyStored;
			for(unsigned int i=0;i<_MyCrystal->getNbAtom();i++) AlreadyStored.push_back(vector<double> ());
			// array of atom types which may be stored without being in the box
			vector<unsigned int> AtType_stored;
			for(unsigned int i=0;i<_MyCrystal->getDoNotSep().size();i++) AtType_stored.push_back(_MyCrystal->getDoNotSep()[i][2]);
			bool ToStore, ToTreat;
			unsigned int id_s;
			int cl_pos_n_i, cl_pos_n_j, cl_pos_n_k, id_notsep;
			double xpos_n, ypos_n, zpos_n, xdiff, ydiff, zdiff, xpos_n_w, ypos_n_w, zpos_n_w;
			double tol = 1e-2;
			// initialize the NotSepTag array
			this->NotSepTag = new vector<int>[this->nbAtom];
			this->IsNotSepTag = true;
			unsigned int count_mol = 1;
			unsigned int count_bonds = 0;
			unsigned int count_angles = 0;
			unsigned int complete_bonds;
			unsigned int complete_angles;
			unsigned int current_bond;
			unsigned int current_angle;
			vector<unsigned int> b_prevs;
			vector<unsigned int> a_prevs;
			vector<unsigned int> a_compl;
			unsigned int nbBondCryst = _MyCrystal->getNbBond();
			unsigned int nbAngleCryst = _MyCrystal->getNbAngle();
			// first loop to store atoms which should not be seperated
			for(int i=-cla1;i<cla1+1;i++){
				for(int j=-cla2;j<cla2+1;j++){
					for(int k=-cla3;k<cla3+1;k++){
						delta_x = i*_MyCrystal->getA1()[0]+j*_MyCrystal->getA2()[0]+k*_MyCrystal->getA3()[0];
						delta_y = i*_MyCrystal->getA1()[1]+j*_MyCrystal->getA2()[1]+k*_MyCrystal->getA3()[1];
						delta_z = i*_MyCrystal->getA1()[2]+j*_MyCrystal->getA2()[2]+k*_MyCrystal->getA3()[2];
						for(unsigned int n=0;n<_MyCrystal->getNbAtom();n++){
							ToTreat = true;
							for(unsigned int sep=0;sep<AtType_stored.size();sep++){
								if( _MyCrystal->getMotif()[n].type_uint == AtType_stored[sep] ){
									ToTreat = false;
									break;
								}
							}
							if( IsMolId &&  _MyCrystal->getNotSepList_size(n) == 0 ) ToTreat = false;
							if( ToTreat ){
								xpos = _MyCrystal->getMotif()[n].pos.x + delta_x;
								ypos = _MyCrystal->getMotif()[n].pos.y + delta_y;
								zpos = _MyCrystal->getMotif()[n].pos.z + delta_z;
								if( xpos > -tolIonPos && xpos < xhi-tolIonPos && ypos > -tolIonPos && ypos < yhi-tolIonPos && zpos > -tolIonPos && zpos < zhi-tolIonPos ){
									if( count == this->nbAtom ){
										break_comp = true;
										break;
									}
									this->AtomList[count] = _MyCrystal->getMotif()[n];
									this->AtomList[count].pos.x = xpos;
									this->AtomList[count].pos.y = ypos;
									this->AtomList[count].pos.z = zpos;
									this->NotSepTag[count].push_back(0);
									if( IsMolId ) MolId[count] = count_mol;
									id_notsep = count;
									// store bonds and angle
									if( IsBond ){
										current_bond = 0;
										b_prevs.clear();
										for(unsigned int b=0;b<nbBondCryst;b++){ // an atom could be involved in multiple bonds, do not break the loop
											for(unsigned int b_i=0;b_i<2;b_i++){
												if( _MyCrystal->getBonds()[(b*2)+b_i]-1 == n ){ // +1 ?
													Bonds[((count_bonds+current_bond)*2)+b_i] = count+1;
													BondType[count_bonds+current_bond] = _MyCrystal->getBondType(b);
													b_prevs.push_back(b);
													current_bond++;
												}
											}
										}
										complete_bonds = 0;
									}
									if( IsAngle ){
										current_angle = 0;
										a_prevs.clear();
										a_compl.clear();
										for(unsigned int a=0;a<nbAngleCryst;a++){ // an atom could be involved in multiple bonds, do not break the loop
											for(unsigned int a_i=0;a_i<3;a_i++){
												if( _MyCrystal->getAngles()[(a*3)+a_i]-1 == n ){ // +1 ?
													Angles[((count_angles+current_angle)*3)+a_i] = count+1;
													AngleType[count_angles+current_angle] = _MyCrystal->getAngleType(a);
													a_prevs.push_back(a);
													a_compl.push_back(1);
													current_angle++;
												}
											}
										}
										complete_angles = 0;
									}
									count += 1;
									// store all its neighboring atoms which should not be separeted from him
									for(unsigned int s=0;s<_MyCrystal->getNotSepList_size(n);s++){
										id_s = _MyCrystal->getNotSepList(n,s*4);
										xpos_n = _MyCrystal->getMotif()[id_s].pos.x + (i+_MyCrystal->getNotSepList(n,s*4+1))*_MyCrystal->getA1()[0] + (j+_MyCrystal->getNotSepList(n,s*4+2))*_MyCrystal->getA2()[0] + (k+_MyCrystal->getNotSepList(n,s*4+3))*_MyCrystal->getA3()[0];
										ypos_n = _MyCrystal->getMotif()[id_s].pos.y + (i+_MyCrystal->getNotSepList(n,s*4+1))*_MyCrystal->getA1()[1] + (j+_MyCrystal->getNotSepList(n,s*4+2))*_MyCrystal->getA2()[1] + (k+_MyCrystal->getNotSepList(n,s*4+3))*_MyCrystal->getA3()[1];
										zpos_n = _MyCrystal->getMotif()[id_s].pos.z + (i+_MyCrystal->getNotSepList(n,s*4+1))*_MyCrystal->getA1()[2] + (j+_MyCrystal->getNotSepList(n,s*4+2))*_MyCrystal->getA2()[2] + (k+_MyCrystal->getNotSepList(n,s*4+3))*_MyCrystal->getA3()[2];
										if( xpos_n < -tolIonPos ) xpos_n_w = xpos_n+xhi;
										else if(xpos_n > xhi-tolIonPos ) xpos_n_w = xpos_n-xhi;
										else xpos_n_w = xpos_n;
										if( ypos_n < -tolIonPos ) ypos_n_w = ypos_n+yhi;
										else if(ypos_n > yhi-tolIonPos ) ypos_n_w = ypos_n-yhi;
										else ypos_n_w = ypos_n;
										if( zpos_n < -tolIonPos ) zpos_n_w = zpos_n+zhi;
										else if(zpos_n > zhi-tolIonPos ) zpos_n_w = zpos_n-zhi;
										else zpos_n_w = zpos_n;
										ToStore = true;
										for(unsigned int t=0;t<AlreadyStored[id_s].size()/3;t++){
											xdiff = fabs(xpos_n_w-AlreadyStored[id_s][t*3]);
											ydiff = fabs(ypos_n_w-AlreadyStored[id_s][t*3+1]);
											zdiff = fabs(zpos_n_w-AlreadyStored[id_s][t*3+2]);
											if( xdiff < tol && ydiff < tol && zdiff < tol ){
												ToStore = false;
												break;
											}
										}
										if( !ToStore ) continue;
										if( count == this->nbAtom ){
											break_comp = true;
											break;
										}
										this->AtomList[count] = _MyCrystal->getMotif()[id_s];
										AlreadyStored[id_s].push_back(xpos_n_w);
										AlreadyStored[id_s].push_back(ypos_n_w);
										AlreadyStored[id_s].push_back(zpos_n_w);
										this->AtomList[count].pos.x = xpos_n;
										this->AtomList[count].pos.y = ypos_n;
										this->AtomList[count].pos.z = zpos_n;
										this->NotSepTag[count].push_back(-1-id_notsep);
										this->NotSepTag[id_notsep].push_back(count);
										this->NotSepTag[id_notsep][0] += 1;
										if( IsMolId ) MolId[count] = count_mol;
										// Store bonds and angles
										if( IsBond ){
											for(unsigned int b=0;b<nbBondCryst;b++){ // an atom could be involved in multiple bonds, do not break the loop
												for(unsigned int b_i=0;b_i<2;b_i++){
													if( _MyCrystal->getBonds()[(b*2)+b_i]-1 == id_s ){ // +1 ?
														bool prev_found = false;
														for(unsigned int b_prev=0;b_prev<current_bond;b_prev++){
															if( b_prevs[b_prev] == b ){
																if( count_bonds+b_prev >= nbBonds ){
																	cout << "WARNING: Exceeding allocated number of bonds (" << nbBonds << ") !" << endl;
																	prev_found = true;
																	break;
																}
																Bonds[((count_bonds+b_prev)*2)+b_i] = count+1;
																complete_bonds++;
																prev_found = true;
															}
														}
														if( !prev_found ){
															Bonds[((count_bonds+current_bond)*2)+b_i] = count+1;
															BondType[count_bonds+current_bond] = _MyCrystal->getBondType(b);
															b_prevs.push_back(b);
															current_bond++;
														}
													}
												}
											}
										}
										if( IsAngle ){
											for(unsigned int a=0;a<nbAngleCryst;a++){ // an atom could be involved in multiple bonds, do not break the loop
												for(unsigned int a_i=0;a_i<3;a_i++){
													if( _MyCrystal->getAngles()[(a*3)+a_i]-1 == id_s ){ // +1 ?
														bool prev_found = false;
														for(unsigned int a_prev=0;a_prev<current_angle;a_prev++){
															if( a_prevs[a_prev] == a ){
																if( count_angles+a_prev >= nbAngles ){
																	cout << "WARNING: Exceeding allocated number of angles (" << nbAngles << ") !" << endl;
																	prev_found = true;
																	break;
																}

																Angles[((count_angles+a_prev)*3)+a_i] = count+1;
																if( a_compl[a_prev] == 1 ) a_compl[a_prev]++;
																else if( a_compl[a_prev] == 2 ) complete_angles++;
																else cout << "WARNING: More than 3 atoms are stored in 1 angle !" << endl;
																prev_found = true;
															}
														}
														if( !prev_found ){
															Angles[((count_angles+current_angle)*3)+a_i] = count+1;
															AngleType[count_angles+current_angle] = _MyCrystal->getAngleType(a);
															a_prevs.push_back(a);
															a_compl.push_back(1);
															current_angle++;
														}
													}
												}
											}
										}

										count += 1;
									} // end loop on neighboring atom to do not sep
									count_bonds += complete_bonds;
									count_angles += complete_angles;
									count_mol++;
								} // end if atom is inside box limits
							} // end if ToTreat DoNotSep or molecule
							if( break_comp) break;
						} // end loop on Crystal motif
						if( break_comp) break;
					} // end loop on CL_z
					if( break_comp) break;
				} // end loop on CL_y
				if( break_comp) break;
			} // end loop on CL_x
			// second loop to store atoms which may have not been stored
			if( count < this->nbAtom-1 ){
				for(int i=-cla1;i<cla1+1;i++){
					for(int j=-cla2;j<cla2+1;j++){
						for(int k=-cla3;k<cla3+1;k++){
							delta_x = i*_MyCrystal->getA1()[0]+j*_MyCrystal->getA2()[0]+k*_MyCrystal->getA3()[0];
							delta_y = i*_MyCrystal->getA1()[1]+j*_MyCrystal->getA2()[1]+k*_MyCrystal->getA3()[1];
							delta_z = i*_MyCrystal->getA1()[2]+j*_MyCrystal->getA2()[2]+k*_MyCrystal->getA3()[2];
							for(unsigned int n=0;n<_MyCrystal->getNbAtom();n++){
								ToTreat = false;
								for(unsigned int sep=0;sep<AtType_stored.size();sep++){
									if( _MyCrystal->getMotif()[n].type_uint == AtType_stored[sep] ){
										ToTreat = true;
										break;
									}
								}
								if( ToTreat ){
									xpos = _MyCrystal->getMotif()[n].pos.x + delta_x;
									ypos = _MyCrystal->getMotif()[n].pos.y + delta_y;
									zpos = _MyCrystal->getMotif()[n].pos.z + delta_z;
									if( xpos > -tolIonPos && xpos < xhi-tolIonPos && ypos > -tolIonPos && ypos < yhi-tolIonPos && zpos > -tolIonPos && zpos < zhi-tolIonPos ){
									// search if this atom type is one of atom which could be stored without being inside the box
										ToStore = true;
										for(unsigned int t=0;t<AlreadyStored[n].size()/3;t++){ //TODO case for which ion is in (0,0,0)
											xdiff = fabs(xpos-AlreadyStored[n][t*3]);
											ydiff = fabs(ypos-AlreadyStored[n][t*3+1]);
											zdiff = fabs(zpos-AlreadyStored[n][t*3+2]);
											if( ( xdiff < tol && ydiff < tol && zdiff < tol ) ){
												ToStore = false;
												break;
											}
										}
										// if not store it
										if( ToStore ){
											if( count == this->nbAtom ){
												break_comp = true;
												break;
											}
											this->AtomList[count] = _MyCrystal->getMotif()[n];
											this->AtomList[count].pos.x = xpos;
											this->AtomList[count].pos.y = ypos;
											this->AtomList[count].pos.z = zpos;
											this->NotSepTag[count].push_back(0);
											count += 1;
										}
									}
								}
							}
						if( break_comp) break;
						}
					if( break_comp) break;
					}
				if( break_comp) break;
				}
			}
		} // END DoNotSep case
		if( count < this->nbAtom ){
			count_box_fill++;
			cla1 *= 2;
			cla2 *= 2;
			cla3 *= 2;
		}else box_fill = true;
	} // end while
	if( count < this->nbAtom ){
		cout << "WARGNING :The number of atoms after constructing orthogonal cell is lower than expected !" << endl;
	}
}
// construct AtomicSystem giving AtomList, number of atom and cell vectors 
AtomicSystem::AtomicSystem(Atom *AtomList, unsigned int nbAtom, Crystal *_MyCrystal, double *H1, double *H2, double *H3){
	Dis = new Displays;
	AtomListConstructor(AtomList,nbAtom,_MyCrystal,H1,H2,H3);
}

AtomicSystem::AtomicSystem(Atom *AtomList, unsigned int nbAtom, Crystal *_MyCrystal, double *H1, double *H2, double *H3, unsigned int *MolId){
	Dis = new Displays;
	AtomListConstructor(AtomList,nbAtom,_MyCrystal,H1,H2,H3);
	this->MolId = MolId;
	IsMolIdMine = false;
	IsMolId = true;
}

AtomicSystem::AtomicSystem(Atom *AtomList, unsigned int nbAtom, Crystal *_MyCrystal, double *H1, double *H2, double *H3, unsigned int *MolId, unsigned int nbBonds, unsigned int nbBondType, unsigned int *Bonds, unsigned int *BondType): nbBonds(nbBonds), nbBondType(nbBondType){
	Dis = new Displays;
	AtomListConstructor(AtomList,nbAtom,_MyCrystal,H1,H2,H3);
	this->MolId = MolId;
	IsMolIdMine = false;
	IsMolId = true;
	this->Bonds = Bonds;
	this->BondType = BondType;
	IsBondMine = false;
	IsBond = true;
}

AtomicSystem::AtomicSystem(Atom *AtomList, unsigned int nbAtom, Crystal *_MyCrystal, double *H1, double *H2, double *H3, unsigned int *MolId, unsigned int nbBonds, unsigned int nbBondType, unsigned int *Bonds, unsigned int *BondType, unsigned int nbAngles, unsigned int nbAngleType, unsigned int *Angles, unsigned int *AngleType): nbBonds(nbBonds), nbBondType(nbBondType), nbAngles(nbAngles), nbAngleType(nbAngleType){
	Dis = new Displays;
	AtomListConstructor(AtomList,nbAtom,_MyCrystal,H1,H2,H3);
	this->MolId = MolId;
	IsMolIdMine = false;
	IsMolId = true;
	this->Bonds = Bonds;
	this->BondType = BondType;
	IsBondMine = false;
	IsBond = true;
	this->Angles = Angles;
	this->AngleType = AngleType;
	IsAngleMine = false;
	IsAngle = true;
}


void AtomicSystem::AtomListConstructor(Atom *AtomList, unsigned int nbAtom, Crystal *_MyCrystal, double *H1, double *H2, double *H3){
	if( this->AtomListConstructed ){
		//delete[] H1;
		//delete[] H2;
		//delete[] H3;
		delete[] G1;
		delete[] G2;
		delete[] G3;
		MT->~MathTools();
	}
	this->AtomList = AtomList;
	this->nbAtom = nbAtom;
	this->_MyCrystal = _MyCrystal;
	this->H1 = H1;
	this->H2 = H2;
	this->H3 = H3;
	read_params_atsys();
	this->IsWrappedPos = false;
	this->IsAtomListMine = false;
	IsElem = true;
	this->IsCellVecMine = false;
	IsCrystalDefined = true;
	this->nbAtomType = _MyCrystal->getNbAtomType();
	this->AtomType = _MyCrystal->getAtomType();
	this->AtomMass = _MyCrystal->getAtomMass();
	this->AtomCharge = _MyCrystal->getAtomCharge();
	this->IsCharge = _MyCrystal->getIsCharge();
	if( H2[0] != 0 || H3[1] != 0 || H3[2] != 0 ) this->IsTilted = true;
	else this->IsTilted = false;
	this->MT = new MathTools;
	this->G1 = new double[3];
	this->G2 = new double[3];
	this->G3 = new double[3];
	this->IsG = true;
	computeInverseCellVec();
	this->AtomListConstructed = true;
}

// Constructor using the filename of an atomic file (.cfg, .lmp, .xsf,..)
AtomicSystem::AtomicSystem(const string& filename, const int &_timestep){
	Dis = new Displays;
	if( !FilenameConstructor(filename, _timestep) ){
		cerr << "Maybe try to change its format using atomsk or generate it with AtomHic (supported format lmp, cfg, lammpstrj)" << endl;
		exit(EXIT_FAILURE);
	}
}

bool AtomicSystem::FilenameConstructor(const string& filename, const int &_timestep){
	if( FilenameConstructed ){
		delete[] AtomType;
		delete[] AtomMass;
		delete[] AtomCharge;
		delete[] H1;
		delete[] H2;
		delete[] H3;
		delete[] G1;
		delete[] G2;
		delete[] G3;
		MT->~MathTools();
		delete[] AtomList; 
	}
	read_params_atsys();
	this->AreTypeMassChargeMine = true;
	this->AtomType = new string[this->MaxAtomType];
	this->AtomMass = new double[this->MaxAtomType];
	this->AtomCharge = new double[this->MaxAtomType];
	for(unsigned int i=0;i<MaxAtomType;i++){
		AtomType[i] = "";
		AtomMass[i] = 0.;
		AtomCharge[i] = 0.;
	}
	this->H1 = new double[3];
	this->H2 = new double[3];
	this->H3 = new double[3];
	for(unsigned int i=0;i<3;i++){
		this->H1[i] = 0.;
		this->H2[i] = 0.;
		this->H3[i] = 0.;
	}
	this->IsCharge=false;
	this->IsTilted=false;
	this->nbAtomType = 0;
	this->MT = new MathTools;
	this->G1 = new double[3];
	this->G2 = new double[3];
	this->G3 = new double[3];
	this->IsG = true;

	if( !ReadAtomicFile(filename,_timestep) ){
		cout << "The atomic file readers cannot read the given file." << endl;
		IsAtomListMine = false;
		return false;
	}
	computeInverseCellVec();
	computeWrap();
	this->FilenameConstructed = true;
	return true;
}

AtomicSystem::AtomicSystem(AtomicSystem *AtSys, unsigned int &nbSys, std::string &dir){
	Dis = new Displays;
	nbAtom = 0;
	nbAtomType = 0;
	nbBonds = 0;
	nbBondType = 0;
	nbAngles = 0;
	nbAngleType = 0;
	this->AreTypeMassChargeMine = true;
	this->AtomType = new string[this->MaxAtomType];
	this->AtomMass = new double[this->MaxAtomType];
	this->AtomCharge = new double[this->MaxAtomType];
	for(unsigned int i=0;i<MaxAtomType;i++){
		AtomType[i] = "";
		AtomMass[i] = 0.;
		AtomCharge[i] = 0.;
	}
	// get properties of the system
	IsElem = AtSys[0].getIsElem();
	vector<unsigned int> *corres_array_elem;
	IsCharge = AtSys[0].getIsCharge();
	IsMolId = AtSys[0].getIsMolId();
	IsBond = AtSys[0].getIsBond();
	IsAngle = AtSys[0].getIsAngle();
	IsVel = AtSys[0].getIsVel();
	if( !IsElem ){
		cout << "The elements are not defined in the AtomicSystems to merge, we will consider only types (integer based) when merging the systems" << endl;
		for(unsigned int i=0;i<AtSys[0].getNbAtomType();i++){
			AtomMass[i] = AtSys[0].getAtomMass(i);
			if( IsCharge ) AtomCharge[i] = AtSys[0].getAtomCharge(i);
		}
	}else corres_array_elem = new vector<unsigned int>[nbSys];
	IsSetAux = AtSys[0].getIsSetAux();
	vector<string> Aux_name_temp;
	vector<unsigned int> Aux_size_temp;
	vector<unsigned int> *corres_array_aux;
	if( IsSetAux ){
		unsigned int beg=0;
		if( IsVel ){
			Aux_name_temp.push_back("Velocities");
			Aux_size_temp.push_back(3);
			beg = 1;
		}
		for(unsigned int n=beg;n<AtSys[0].getNbAux();n++){
			Aux_name_temp.push_back(AtSys[0].getAux_name(n));
			Aux_size_temp.push_back(AtSys[0].getAux_size(n));
		}
	}
	
	this->H1 = new double[3];
	this->H2 = new double[3];
	this->H3 = new double[3];
	for(unsigned int i=0;i<3;i++){
		this->H1[i] = AtSys[0].getH1()[i];
		this->H2[i] = AtSys[0].getH2()[i];
		this->H3[i] = AtSys[0].getH3()[i];
	}
	if( dir != "x" && dir != "y" && dir != "z" ){
		cerr << "Merging direction should be x, y, or z, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	
	// first loop on number of system to determine the global attribute of the merged system (nbAtom, nbAtomType, nbAux/Bond/Angle..)
	for(unsigned int i=0;i<nbSys;i++){
		nbAtom += AtSys[i].getNbAtom();
		if( i != 0 ){
			if( dir == "x" ){
				for(unsigned int d=0;d<3;d++){
					H1[d] += AtSys[i].getH1()[d];
					if( H2[d] != AtSys[i].getH2()[d] || H3[d] != AtSys[i].getH3()[d] ) cout << "Warning, the cell vectors in the direction other than the merging one are different between the systems to merge, the ones of the first system provided will be kept for the final merged system" << endl;
				}
			}else if( dir == "y" ){
				for(unsigned int d=0;d<3;d++){
					H2[d] += AtSys[i].getH2()[d];
					if( H1[d] != AtSys[i].getH1()[d] || H3[d] != AtSys[i].getH3()[d] ) cout << "Warning, the cell vectors in the direction other than the merging one are different between the systems to merge, the ones of the first system provided will be kept for the final merged system" << endl;
				}
			}else{
				for(unsigned int d=0;d<3;d++){
					H3[d] += AtSys[i].getH3()[d];
					if( H2[d] != AtSys[i].getH2()[d] || H1[d] != AtSys[i].getH1()[d] ) cout << "Warning, the cell vectors in the direction other than the merging one are different between the systems to merge, the ones of the first system provided will be kept for the final merged system" << endl;
				}
			}
			// manage aux prop
			if( IsSetAux ){
				unsigned int beg = 0;
				if( IsVel ){
				       	if( !AtSys[i].getIsVel() ){
						cout << "All the systems have not velocities defined, in the merged system the velocities will then not be defined" << endl;
						Aux_name_temp.erase(Aux_name_temp.begin());
						Aux_size_temp.erase(Aux_size_temp.begin());
						IsVel = false;
					}else beg = 1;
				}
			        if( !AtSys[i].getIsSetAux() ){
					cout << "All the systems have not auxiliary properties defined, in the merged system the auxiliary properties will then not be defined" << endl;
					IsSetAux = false;
				}else{
					vector<unsigned int> indextodel;
					for(unsigned int n=beg;n<Aux_name_temp.size();n++){
						bool already = false;
						for(unsigned int n2=0;n2<AtSys[i].getNbAux();n2++){
							if( Aux_name_temp[n] == AtSys[i].getAux_name(n2) && Aux_size_temp[n] == AtSys[i].getAux_size(n2) ){
								already = true;
								break;
							}
						}
						if( !already ){
							cout << "The auxiliary property " << Aux_name_temp[n] << " is not defined in all system to merge, it will then not be defined in the merged system" << endl;
							indextodel.push_back(n);
						}
					}
					for(unsigned int n=indextodel.size()-1;n==0;n--){
						Aux_name_temp.erase(Aux_name_temp.begin()+indextodel[n]);
						Aux_size_temp.erase(Aux_size_temp.begin()+indextodel[n]);
					}
				}
				if( Aux_name_temp.size() == 0 ) IsSetAux = false;
			}
		} // end if i != 0

		if( IsElem ){
			if( !AtSys[i].getIsElem() ){
				cout << "Warning, the elements are not defined in all systems, this may cause mislabelling of elements in the merged system" << endl;
				IsElem = false;
			}else{
				for(unsigned int n=0;n<AtSys[i].getNbAtomType();n++){
					if( !AtSys[i].getIsElem() ) cout << "Warning the elements are not defined in all the AtomicSystems to merge, it may lead to mislabelling of atom type" << endl;
					if( IsCharge && !AtSys[i].getIsCharge() )  cout << "Warning the charges are not defined in all the AtomicSystems to merge, we will then consider a charge of 0 when it is not defined" << endl;
					bool already = false;
					for(unsigned int n2=0;n2<nbAtomType;n2++){
						if( AtomType[n2] == AtSys[i].getAtomType(n) ){
							corres_array_elem[i].push_back(n2+1);
							already = true;
							break;
						}
					}
					if( !already ){
						AtomType[nbAtomType] = AtSys[i].getAtomType(n);
						AtomMass[nbAtomType] = AtSys[i].getAtomMass(n);
						if( IsCharge ){
							if( !AtSys[i].getIsCharge() ) AtomCharge[nbAtomType] = 0.;
							else AtomCharge[nbAtomType] = AtSys[i].getAtomCharge(n);
						}
						nbAtomType++;
						if( nbAtomType >= MaxAtomType ){
						       cerr << "The number of atom type is higher than the maximum number of atom type allowed (" << MaxAtomType << "), aborting" << endl;
					       		exit(EXIT_FAILURE);
						}		
						corres_array_elem[i].push_back(nbAtomType);
					}
				}
			}
		}else if( AtSys[i].getNbAtomType() > nbAtomType ){
			for(unsigned int n=nbAtomType;n<AtSys[i].getNbAtomType();n++){
				AtomMass[n] = AtSys[i].getAtomMass(n);
				AtomCharge[n] = AtSys[i].getAtomCharge(n);
			}
			nbAtomType = AtSys[i].getNbAtomType();
		}
		
		if( IsMolId && !AtSys[i].getIsMolId() ){
			cout << "All the systems have not molecule id defined, in the merged system the molecule id will then not be defined" << endl;
			IsMolId = false;
			
		}
		if( IsBond ){
			if( !AtSys[i].getIsBond() ){
				cout << "All the systems have not bonds defined, in the merged system the bonds will then not be defined" << endl;
				IsBond = false;
			}else{
				nbBonds += AtSys[i].getNbBonds();
				if( AtSys[i].getNbBondType() > nbBondType ) nbBondType = AtSys[i].getNbBondType();
			}
		}
		if( IsAngle ){
		        if( !AtSys[i].getIsAngle() ){
				cout << "All the systems have not angles defined, in the merged system the angles will then not be defined" << endl;
				IsAngle = false;
			}else{
				nbAngles += AtSys[i].getNbAngles();
				if( AtSys[i].getNbAngleType() > nbAngleType ) nbAngleType = AtSys[i].getNbAngleType();
			}

		}
	} // end first loop on nbSys
	
	// Initialize pointers that depends on nbAtom/Bond/Angle/Aux
	AtomList = new Atom[nbAtom];
	if( IsMolId ) MolId = new unsigned int[nbAtom];
	if( IsBond ){
		Bonds = new unsigned int[nbBonds*2];
		BondType = new unsigned int[nbBonds];
	}
	if( IsAngle ){
		Angles = new unsigned int[nbAngles*3];
		AngleType = new unsigned int[nbAngles];
	}
	if( IsSetAux ){
		corres_array_aux = new vector<unsigned int>[nbSys];
		unsigned int buffer;
		for(unsigned int i=0;i<Aux_name_temp.size();i++){
			for(unsigned int n=0;n<nbSys;n++) corres_array_aux[n].push_back(AtSys[n].getAuxIdAndSize(Aux_name_temp[i],buffer));
			Aux_name.push_back(Aux_name_temp[i]);
			Aux_size.push_back(Aux_size_temp[i]);
			Aux.push_back(new double[nbAtom*Aux_size[i]]);
		}
	}

	// Merge the systems
	unsigned int current_nbAt = 0;
	unsigned int current_nbBonds = 0;
	unsigned int current_nbAngles = 0;
	unsigned int prev_nbAt = 0;
	unsigned int prev_nbMol = 0;
	double *current_H = new double[3];
	for(unsigned int d=0;d<3;d++) current_H[d] = 0.;
	for(unsigned int i=0;i<nbSys;i++){
		// Atoms, Molecule Id and Aux
		for(unsigned int n=0;n<AtSys[i].getNbAtom();n++){
			if( current_nbAt >= nbAtom ){
				cerr << "The number of stored atoms exceed the computed number of atom, aborting" << endl;
				exit(EXIT_FAILURE);
			}
			AtomList[current_nbAt] = AtSys[i].getAtom(n);
			AtomList[current_nbAt].pos.x += current_H[0];
			AtomList[current_nbAt].pos.y += current_H[1];
			AtomList[current_nbAt].pos.z += current_H[2];
			unsigned int curtype = AtomList[current_nbAt].type_uint;
			if( IsElem ) AtomList[current_nbAt].type_uint = corres_array_elem[i][curtype-1];
			if( IsMolId ) MolId[current_nbAt] = AtSys[i].getMolId(n) + prev_nbMol;
			if( IsSetAux ){
				for(unsigned int a=0;a<Aux.size();a++){
					for(unsigned int ad=0;ad<Aux_size[a];ad++) Aux[a][current_nbAt*Aux_size[a]+ad] = AtSys[i].getAux(corres_array_aux[i][a])[n];
				}
			}
			current_nbAt++;
		}
		if( dir == "x" ) for(unsigned int d=0;d<3;d++) current_H[d] += AtSys[i].getH1()[d];
		else if( dir == "y" ) for(unsigned int d=0;d<3;d++) current_H[d] += AtSys[i].getH2()[d];
		else if( dir == "z" ) for(unsigned int d=0;d<3;d++) current_H[d] += AtSys[i].getH3()[d];
		// Bonds
		if( IsBond ){
			for(unsigned int n=0;n<AtSys[i].getNbBonds();n++){
				if( current_nbBonds >= nbBonds ){
					cerr << "The number of stored bonds exceed the computed number of bond, aborting" << endl;
					exit(EXIT_FAILURE);
				}
				for(unsigned int d=0;d<2;d++) Bonds[current_nbBonds*2+d] = AtSys[i].getBonds()[n*2+d] + prev_nbAt;
				BondType[current_nbBonds] = AtSys[i].getBondType(n);
				current_nbBonds++;
			}
		}
		// Angles
		if( IsAngle ){
			for(unsigned int n=0;n<AtSys[i].getNbAngles();n++){
				if( current_nbAngles >= nbAngles ){
					cerr << "The number of stored angles exceed the computed number of angle, aborting" << endl;
					exit(EXIT_FAILURE);
				}
				for(unsigned int d=0;d<3;d++) Angles[current_nbAngles*3+d] = AtSys[i].getAngles()[n*3+d] + prev_nbAt;
				AngleType[current_nbAngles] = AtSys[i].getAngleType(n);
				current_nbAngles++;
			}
		}

		prev_nbAt += AtSys[i].getNbAtom();
		if( IsMolId ) prev_nbMol += MT->max_p(MolId,current_nbAt);
	}

	this->MT = new MathTools;
	this->G1 = new double[3];
	this->G2 = new double[3];
	this->G3 = new double[3];
	this->IsG = true;
	computeInverseCellVec();

	if( IsElem ){
		for(unsigned int i=0;i<nbSys;i++) corres_array_elem[i].clear();
		delete[] corres_array_elem;
	}
	if( IsSetAux ){
		for(unsigned int i=0;i<nbSys;i++) corres_array_aux[i].clear();
		delete[] corres_array_aux;
	}
	delete[] current_H;

}

bool AtomicSystem::ReadAtomicFile(const string &filename, const int &_timestep){
	string ext=filename.substr(filename.find_last_of(".") + 1);
	if( _timestep == -1 ) cout << "Reading " << filename << " file..";
	else cout << "Timestep " << _timestep;
	if( ext == "lmp" ){
		if( !this->read_lmp_file(filename) ){
			if( !this->read_cfg_file(filename,_timestep) ){
				if( !this->read_other_cfg(filename) ){
					return false;
				}else{
					cout << " done ! (" << SystemCharacteristics() << ")" << endl;
					return true;
				}
			}else{
				cout << " done ! (" << SystemCharacteristics() << ")" << endl;
				return true;
			}
		}else{
			cout << " done ! (" << SystemCharacteristics() << ")" << endl;
				return true;
		}
	}else{
		if( !this->read_cfg_file(filename,_timestep) ){
			if( !this->read_other_cfg(filename) ){
				if( !this->read_lmp_file(filename) ){
					return false;
				}else{
					cout << " done ! (" << SystemCharacteristics() << ")" << endl;
					return true;
				}
			}else{
				cout << " done ! (" << SystemCharacteristics() << ")" << endl;
				return true;
			}
		}else{
			cout << " done ! (" << SystemCharacteristics() << ")" << endl;
			return true;
		}
	}
}

string AtomicSystem::SystemCharacteristics(bool cfg){
	string msg=to_string(nbAtom);
	if( IsMolId ) msg += " full";
	else if( IsCharge ) msg += " charged";
	msg += " atoms";
	if( IsVel ) msg += " and velocities";
	if( !cfg ){
		if( IsBond ) msg += ", "+to_string(nbBonds)+" bonds";
		if( IsAngle ) msg += ", "+to_string(nbAngles)+" angles";
	}
	return msg;
}

void AtomicSystem::computeInverseCellVec(){
	if( !this->IsG ){
		this->G1 = new double[3];
		this->G2 = new double[3];
		this->G3 = new double[3];
		this->IsG = true;
	}
	double det = this->H1[0]*this->H2[1]*this->H3[2];
	this->G1[0] = this->H2[1]*this->H3[2]/det;
	this->G1[1] = 0.;
	this->G1[2] = 0.;
	this->G2[0] = -(this->H2[0]*this->H3[2])/det;
	this->G2[1] = this->H1[0]*this->H3[2]/det;
	this->G2[2] = -(this->H1[0]*this->H2[2])/det;
	this->G3[0] = ((this->H2[0]*this->H3[1])-(this->H2[1]*this->H3[0]))/det;
	this->G3[1] = -(this->H3[1]*this->H1[0])/det;
	this->G3[2] = this->H2[1]*this->H1[0]/det;
}

unsigned int AtomicSystem::Compute2dDensity(std::string auxname, std::string dir, double sigma, unsigned int nbPtsDir1, unsigned int nbPtsDir2){
	if( !IsWrappedPos ) computeWrap();
	bool already_stored = false;
	unsigned int ind_dens = 0;
	for(unsigned int i=0;i<density_prof_2D.size();i++){
		if( auxname == density_name_2D[i][0] && dir == density_name_2D[i][1] ){
		       ind_dens = i;
		       already_stored = true;
	       	       break;
		}
	}
	if( auxname != "Mass" ){
		cerr << "Computing 2D density with property different than mass is not implemented for the moment" << endl;
		exit(EXIT_FAILURE);
	}
	unsigned int nbPts = nbPtsDir1*nbPtsDir2;
	if( !already_stored ){
		density_prof_2D.push_back(new double[nbPts*3]); // density_prof_2D[i*3] = x coordinate of point i [i*3+1] = y coord [i*3+2] value of density at point i
		ind_dens = density_prof_2D.size()-1;
	}
	int indexaux=-1;
	for(unsigned int i=0;i<this->Aux_name.size();i++){
		if( auxname == Aux_name[i] ){
			indexaux=i;
			break;
		}
	}
	if( indexaux == -1 ){
		if( auxname == "Mass" ){ // first not clean implemantation
			cout << "Computing 2d mass density" << endl;
			if( this->AtomMass[0] == 0. ){
				cout << "Atomic mass is not defined, setting all masses to unity" << endl;
				for(unsigned int i=0;i<nbAtomType;i++) this->AtomMass[i] = 1.;
			}
			if( dir == "x" ){ // Dir1 => y, Dir2 => z
				#pragma omp parallel for
				for(unsigned int i1=0;i1<nbPtsDir1;i1++){
					for(unsigned int i2=0;i2<nbPtsDir2;i2++){
						unsigned int ind1 = (i1*nbPtsDir2+i2)*3;
						unsigned int ind2 = (i1*nbPtsDir2+i2)*3+1;
						unsigned int ind3 = (i1*nbPtsDir2+i2)*3+2;
						this->density_prof_2D[ind_dens][ind1] = this->H2[1]*i1/(nbPtsDir1-1.); // consider only not tilted box for the moment (TODO)
						this->density_prof_2D[ind_dens][ind2] = this->H3[2]*i2/(nbPtsDir2-1.); // consider only not tilted box for the moment (TODO)
						this->density_prof_2D[ind_dens][ind3] = 0.;
						for(unsigned int j=0;j<this->nbAtom;j++){
							// consider boundary condition
							for(int bc1=-1;bc1<=1;bc1++){
								for(int bc2=-1;bc2<=1;bc2++){
									double coord1 = this->WrappedPos[j].y+(this->H2[1]*bc1);
									double coord2 = this->WrappedPos[j].z+(this->H3[2]*bc2);
									this->density_prof_2D[ind_dens][ind3] += MT->gaussian(this->density_prof_2D[ind_dens][ind1], this->density_prof_2D[ind_dens][ind2], coord1, coord2, sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
								}
							}
						}
					}
				}
			}else if( dir == "y" ){ // Dir1 => z, Dir2 => x
				#pragma omp parallel for
				for(unsigned int i1=0;i1<nbPtsDir1;i1++){
					for(unsigned int i2=0;i2<nbPtsDir2;i2++){
						unsigned int ind1 = (i1*nbPtsDir2+i2)*3;
						unsigned int ind2 = (i1*nbPtsDir2+i2)*3+1;
						unsigned int ind3 = (i1*nbPtsDir2+i2)*3+2;
						this->density_prof_2D[ind_dens][ind1] = this->H3[2]*i1/(nbPtsDir1-1.); // consider only not tilted box for the moment (TODO)
						this->density_prof_2D[ind_dens][ind2] = this->H1[0]*i2/(nbPtsDir2-1.); // consider only not tilted box for the moment (TODO)
						this->density_prof_2D[ind_dens][ind3] = 0;
						for(unsigned int j=0;j<this->nbAtom;j++){
							// consider boundary condition
							for(int bc1=-1;bc1<=1;bc1++){
								for(int bc2=-1;bc2<=1;bc2++){
									double coord1 = this->WrappedPos[j].z+(this->H3[2]*bc1);
									double coord2 = this->WrappedPos[j].x+(this->H1[0]*bc2);
									this->density_prof_2D[ind_dens][ind3] += MT->gaussian(this->density_prof_2D[ind_dens][ind1], this->density_prof_2D[ind_dens][ind2], coord1, coord2, sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
								}
							}
						}
					}
				}
			}else if( dir == "z" ){ // Dir1 => x, Dir2 => y
				#pragma omp parallel for
				for(unsigned int i1=0;i1<nbPtsDir1;i1++){
					for(unsigned int i2=0;i2<nbPtsDir2;i2++){
						unsigned int ind1 = (i1*nbPtsDir2+i2)*3;
						unsigned int ind2 = (i1*nbPtsDir2+i2)*3+1;
						unsigned int ind3 = (i1*nbPtsDir2+i2)*3+2;
						this->density_prof_2D[ind_dens][ind1] = this->H1[0]*i1/(nbPtsDir1-1.); // consider only not tilted box for the moment (TODO)
						this->density_prof_2D[ind_dens][ind2] = this->H2[1]*i2/(nbPtsDir2-1.); // consider only not tilted box for the moment (TODO)
						this->density_prof_2D[ind_dens][ind3] = 0;
						for(unsigned int j=0;j<this->nbAtom;j++){
							// consider boundary condition
							for(int bc1=-1;bc1<=1;bc1++){
								for(int bc2=-1;bc2<=1;bc2++){
									double coord1 = this->WrappedPos[j].x+(this->H1[0]*bc1);
									double coord2 = this->WrappedPos[j].y+(this->H2[1]*bc2);
									this->density_prof_2D[ind_dens][ind3] += MT->gaussian(this->density_prof_2D[ind_dens][ind1], this->density_prof_2D[ind_dens][ind2], coord1, coord2, sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
								}
							}
						}
					}
				}
			}
			cout << "Done !" << endl;
		}
	}
	if( already_stored ){
		this->density_nbPts_2D[ind_dens*2] = nbPtsDir1;
		this->density_nbPts_2D[ind_dens*2+1] = nbPtsDir2;
	}else{
		this->density_name_2D.push_back(new string[2]);
		this->density_name_2D[this->density_name_2D.size()-1][0] = auxname;
		this->density_name_2D[this->density_name_2D.size()-1][1] = dir;
		this->density_nbPts_2D.push_back(nbPtsDir1);
		this->density_nbPts_2D.push_back(nbPtsDir2);
	}

	return ind_dens;
}

unsigned int AtomicSystem::Compute1dDensity(std::string auxname, std::string dir, double sigma, unsigned int nbPts){
	if( !IsWrappedPos ) computeWrap();
	// search if this density has already been computed or not
	bool already_stored = false;
	unsigned int ind_dens = 0;
	for(unsigned int i=0;i<density_prof.size();i++){
		if( auxname == density_name[i][0] && dir == density_name[i][1] ){
		       ind_dens = i;
		       already_stored = true;
	       	       break;
		}
	}
	if( !already_stored ){
		density_prof.push_back(new double[nbPts*2]);
		ind_dens = density_prof.size()-1;
	}
	int indexaux=-1;
	for(unsigned int i=0;i<this->Aux_name.size();i++){
		if( auxname == Aux_name[i] ){
			indexaux=i;
			break;
		}
	}
	if( indexaux == -1 ){
		if( auxname == "Mass" ){
			if( dir == "x" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[ind_dens][i*2] = 0;
					this->density_prof[ind_dens][i*2+1] = this->H1[0]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++){
						this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].x, sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
						// consider boundary conditions
						this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].x+this->H1[0], sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
						this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].x-this->H1[0], sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
					}
				}
			}else if( dir == "y" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[ind_dens][i*2] = 0;
					this->density_prof[ind_dens][i*2+1] = this->H2[1]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++){
						this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].y, sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
						// BC
						this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].y+this->H2[1], sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
						this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].y-this->H2[1], sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
					}
				}
			}else if( dir == "z" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[ind_dens][i*2] = 0;
					this->density_prof[ind_dens][i*2+1] = this->H3[2]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].z, sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
				}
			}else{
			        cout << "The provided direction \"" << dir << "\" for density computation has not been recognize" << endl;
				return 0;
			}
		}else if( auxname == "Charge" ){
			if( dir == "x" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[ind_dens][i*2] = 0;
					this->density_prof[ind_dens][i*2+1] = this->H1[0]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].x, sigma)*(this->AtomCharge[this->AtomList[j].type_uint-1]);
				}
			}else if( dir == "y" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[ind_dens][i*2] = 0;
					this->density_prof[ind_dens][i*2+1] = this->H2[1]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].y, sigma)*(this->AtomCharge[this->AtomList[j].type_uint-1]);
				}
			}else if( dir == "z" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[ind_dens][i*2] = 0;
					this->density_prof[ind_dens][i*2+1] = this->H3[2]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].z, sigma)*(this->AtomCharge[this->AtomList[j].type_uint-1]);
				}
			}else{
			        cout << "The provided direction \"" << dir << "\" for density computation has not been recognize" << endl;
				return 0;
			}
		}else if ( auxname == "Atomic" ){
			if( dir == "x" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[ind_dens][i*2] = 0;
					this->density_prof[ind_dens][i*2+1] = this->H1[0]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].x, sigma);
				}
			}else if( dir == "y" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[ind_dens][i*2] = 0;
					this->density_prof[ind_dens][i*2+1] = this->H2[1]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++){
						this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].y, sigma);
						//BC
						this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].y+this->H2[1], sigma);
						this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].y-this->H2[1], sigma);
					}
				}
			}else if( dir == "z" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[ind_dens][i*2] = 0;
					this->density_prof[ind_dens][i*2+1] = this->H3[2]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].z, sigma);
				}
			}else{
			        cout << "The provided direction \"" << dir << "\" for density computation has not been recognize (directions available : x y z)" << endl;
				return 0;
			}
		}else{
			cout << "The property used to compute density does not exist (properties available : Mass, Charge, Atomic";
		        for(unsigned int i=0;i<Aux_name.size();i++) cout << ", " << Aux_name[i];
	       	        cout << ")" << endl;
			return 0;
		}
	}else{
		if( dir == "x" ){
			for(unsigned int i=0;i<nbPts;i++){
				this->density_prof[ind_dens][i*2] = 0;
				this->density_prof[ind_dens][i*2+1] = this->H1[0]*i/(nbPts-1.);
				for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].x, sigma)*(this->Aux[indexaux][j]);
			}
		}else if( dir == "y" ){
			for(unsigned int i=0;i<nbPts;i++){
				this->density_prof[ind_dens][i*2] = 0;
				this->density_prof[ind_dens][i*2+1] = this->H2[1]*i/(nbPts-1.);
				for(unsigned int j=0;j<this->nbAtom;j++){
					this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].y, sigma)*(this->Aux[indexaux][j]);
					//BC
					this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].y+this->H2[1], sigma)*(this->Aux[indexaux][j]);
					this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].y-this->H2[1], sigma)*(this->Aux[indexaux][j]);
				}
			}
		}else if( dir == "z" ){
			for(unsigned int i=0;i<nbPts;i++){
				this->density_prof[ind_dens][i*2] = 0;
				this->density_prof[ind_dens][i*2+1] = this->H3[2]*i/(nbPts-1.);
				for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[ind_dens][i*2] += MT->gaussian(this->density_prof[ind_dens][i*2+1], this->WrappedPos[j].z, sigma)*(this->Aux[indexaux][j]);
			}
		}else{
		        cout << "The provided direction \"" << dir << "\" for density computation has not been recognize (directions available : x y z)" << endl;
			return 0;
		}
	}
	if( already_stored ) this->density_nbPts[ind_dens] = nbPts;
	else{
		this->density_name.push_back(new string[2]);
		this->density_name[this->density_name.size()-1][0] = auxname;
		this->density_name[this->density_name.size()-1][1] = dir;
		this->density_nbPts.push_back(nbPts);
	}
	return ind_dens;
}

void AtomicSystem::Print1dDensity(string filename, string auxname){
	int indexaux = -1;
	for(unsigned int i=0;i<this->density_name.size();i++){
		if( this->density_name[i][0] == auxname ){
			indexaux = i;
			break;
		}
	}
	if( indexaux == -1 ) cout << "The density property to print does not exist" << endl;
	ofstream writefile(filename);
	writefile << this->density_name[indexaux][1] << " " << this->density_name[indexaux][0] << "Density" << endl;
	for(unsigned int i=0;i<this->density_nbPts[indexaux];i++) writefile << this->density_prof[indexaux][i*2+1] << " " << this->density_prof[indexaux][i*2] << endl; 
	writefile.close();
	cout << filename << " file writted" << endl;
}

void AtomicSystem::Print2dDensity(string filename, string auxname){
	int indexaux = -1;
	for(unsigned int i=0;i<this->density_name_2D.size();i++){
		if( this->density_name_2D[i][0] == auxname ){
			indexaux = i;
			break;
		}
	}
	if( indexaux == -1 ) cout << "The density property to print does not exist" << endl;
	ofstream writefile(filename);
	writefile << "#" << this->density_name_2D[indexaux][1] << " " << this->density_name_2D[indexaux][0] << "Density" << endl;
	unsigned int count = 0;
	for(unsigned int i=0;i<this->density_nbPts_2D[indexaux*2];i++){
		for(unsigned int j=0;j<this->density_nbPts_2D[indexaux*2+1];j++){
			writefile << this->density_prof_2D[indexaux][count*3] << " " << this->density_prof_2D[indexaux][count*3+1] << " " << this->density_prof_2D[indexaux][count*3+2] << endl; 
			count++;
		}
	}
	writefile.close();
	cout << filename << " file writted" << endl;
}


void AtomicSystem::setCrystal(Crystal* MyCrystal){
	this->_MyCrystal = MyCrystal;
	this->IsCrystalDefined = true;
	if( FilenameConstructed && nbAtomType > 1 ) UpdateTypes2Crystal();
}

void AtomicSystem::setCrystal(const std::string& CrystalName){
	cout << "Setting crystal : " << CrystalName << endl;
	this->_MyCrystal = new Crystal(CrystalName);
	this->IsCrystalDefined = true;
	this->IsCrystalMine = true;
	if( FilenameConstructed && nbAtomType > 1 ) UpdateTypes2Crystal();
}

void AtomicSystem::ComputeNotSepList(){
	if( !IsCrystalDefined ){
		cout << "Warning, the crystal is not defined, we cannot compute the DoNotSepare list" << endl;
		return;
	}else if( !_MyCrystal->getIsDoNotSep() ){
		cout << "The crystal does not have a DoNotSepare list, we then cannot compute the list for the atomic system" << endl;
		return;
	}else if( IsNotSepTag ){
		//delete[] NotSepTag;
	       //cout << "The DoNotSepare list is already computed for this system, it will then not be computed" << endl;
	       return;
	}
	this->NotSepTag = new vector<int>[this->nbAtom];
	for(unsigned int i=0;i<nbAtom;i++) NotSepTag[i].push_back(0);
	this->IsNotSepTag = true;
	double rcut = MT->max_p(_MyCrystal->getALength(),3);
	rcut *= 1.2;
	searchNeighbours(rcut);
	unsigned int DNS_size = _MyCrystal->getDoNotSep().size();
	//#pragma omp parallel for
	for(unsigned int i=0;i<nbAtom;i++){
		double xp, yp ,zp, xpos, ypos, zpos;
		bool ToTreat = false;
		vector<unsigned int> type2search, nbNeigh;
		for(unsigned int j=0;j<DNS_size;j++){
			if( AtomList[i].type_uint == _MyCrystal->getDoNotSep()[j][0] ){
				ToTreat = true;
				nbNeigh.push_back(_MyCrystal->getDoNotSep()[j][1]);
				type2search.push_back(_MyCrystal->getDoNotSep()[j][2]);
			}
		}
		if( !ToTreat ) continue;
		xpos = WrappedPos[i].x;
		ypos = WrappedPos[i].y;
		zpos = WrappedPos[i].z;
		for(unsigned int n=0;n<nbNeigh.size();n++){
			vector<double> ToSort;
			for(unsigned int j=0;j<Neighbours[i*(nbMaxN+1)];j++){
				unsigned int id = Neighbours[i*(nbMaxN+1)+j+1];
				if( AtomList[id].type_uint == type2search[n] ){
					long unsigned int id1 = i*nbMaxN*3+j*3;
					long unsigned int id2 = id1+1;
					long unsigned int id3 = id1+2;
					xp = WrappedPos[id].x+CLNeighbours[id1]*H1[0]+CLNeighbours[id2]*H2[0]+CLNeighbours[id3]*H3[0]-xpos;
					yp = WrappedPos[id].y+CLNeighbours[id1]*H1[1]+CLNeighbours[id2]*H2[1]+CLNeighbours[id3]*H3[1]-ypos;
					zp = WrappedPos[id].z+CLNeighbours[id1]*H1[2]+CLNeighbours[id2]*H2[2]+CLNeighbours[id3]*H3[2]-zpos;
					ToSort.push_back(xp*xp+yp*yp+zp*zp);
					ToSort.push_back(id);
				}
			}
			MT->sort(ToSort,0,2,ToSort);
			unsigned int currentneigh = ToSort.size()/2;
			unsigned int j=0;
			unsigned int nbneighstored = 0;
			while( nbneighstored < nbNeigh[n] ){
				if( j == currentneigh ){
					cout << "Warning, not enough ions have been found to construct the DoNotSepare list, the cutoff should be increased" << endl;
					break;
				}
				if( NotSepTag[(unsigned int) (ToSort[j*2+1])][0] != -1 ){
					NotSepTag[i][0]++;
					NotSepTag[i].push_back((unsigned int) (ToSort[j*2+1]));
					NotSepTag[(unsigned int) (ToSort[j*2+1])][0] = -1-i;
					nbneighstored++;
				}
				j++;
			}
		}
	}

}

// This function modifies the ordering of type and type_uint of _MyCrystal for consistency
void AtomicSystem::UpdateTypes2Crystal(){
	if( !IsCrystalDefined ){
		cerr << "The crystal is not defined, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	if( !IsElem ){
		vector<unsigned int> type_t2e;
		vector<string> element_t2e;
		unsigned int buffer_i;
		string buffer_s, line;
		cout << "The element names are not provided in the atomic system" << endl;
		ifstream file_e2t("Type2Element.ath", ifstream::in);
		if( file_e2t ){
			cout << "Reading the correspondance between type and element in Type2Element.ath" << endl;
			while(getline(file_e2t,line)){
				istringstream text(line);
				text >> buffer_i >> buffer_s;
				type_t2e.push_back(buffer_i);
				element_t2e.push_back(buffer_s);
			}
			nbAtomType = type_t2e.size();
			for(unsigned int t=0;t<nbAtomType;t++) AtomType[type_t2e[t]-1] = element_t2e[t];
			IsElem = true;
		}else{
			cout << "You can provide a Type2Element.ath file giving the correspondance between type and element in the working directory (an example is present in /data/ExampleFiles/)" << endl;
			cout << "As this file has not been found, we will simply based correspondance with crystal database on type (integer-based) which could lead to dramatic confusion depending on what you are doing" << endl;
		}
	}
	bool ok = true;
	if( nbAtomType == 0 || nbAtomType < _MyCrystal->getNbAtomType() ){ 
		cerr << "We cannot update atomic system atom types with crystal" << endl;
	       exit(EXIT_FAILURE);
	}
	// we want to allows considering ion types which are not in the crystal (e.g. solute ions)
	if( IsElem ){
		for(unsigned int t=0;t<_MyCrystal->getNbAtomType();t++){
			if( AtomType[t] != _MyCrystal->getAtomType()[t] ){
				ok = false;
				break;
			}
		}
		if( ok ) return;
		else{
			// Compute the correspondance array
			unsigned int *CorresArray = new unsigned int[_MyCrystal->getNbAtomType()];
			unsigned int nb_type_found = 0;
			for(unsigned int tc=0;tc<_MyCrystal->getNbAtomType();tc++){
				for(unsigned int ta=0;ta<nbAtomType;ta++){
					if( AtomType[ta] == _MyCrystal->getAtomType()[tc] ){
						CorresArray[tc] = ta;
						nb_type_found++;
						break;
					}
				}
			}
			if( nb_type_found != _MyCrystal->getNbAtomType() ){
				cout << "Warning, it seems that the provided crystal does not correspond to the atomic system" << endl;
				return;
			}
			_MyCrystal->ChangeTypes(CorresArray);
			delete[] CorresArray;
		}
	}
}

void AtomicSystem::computeWrap(){
	//
	// 
    if(!this->IsWrappedPos){
	WrappedPos = new Position[this->nbAtom]; 
	IsWrappedPos = true;
    }
	computeInverseCellVec();
	double x,y,z;
	for(unsigned int i=0;i<this->nbAtom;i++){
		// compute reduced coordinates
		x = this->AtomList[i].pos.x*G1[0]+this->AtomList[i].pos.y*G2[0]+this->AtomList[i].pos.z*G3[0];
		y = this->AtomList[i].pos.x*G1[1]+this->AtomList[i].pos.y*G2[1]+this->AtomList[i].pos.z*G3[1];
		z = this->AtomList[i].pos.x*G1[2]+this->AtomList[i].pos.y*G2[2]+this->AtomList[i].pos.z*G3[2];
		if( x >= 1. || x < 0. ) x = x-floor(x);
		if( y >= 1. || y < 0. ) y = y-floor(y);
		if( z >= 1. || z < 0. ) z = z-floor(z);
		// cartesian coordinates
		this->WrappedPos[i].x = x*H1[0]+y*H2[0]+z*H3[0];
		this->WrappedPos[i].y = x*H1[1]+y*H2[1]+z*H3[1];
		this->WrappedPos[i].z = x*H1[2]+y*H2[2]+z*H3[2];
	}
}

// Searching neighbours using cell list algorithm
unsigned int AtomicSystem::searchNeighbours(const double& rc){
	cout << "Performing neighbour research.." << endl;
	//auto start = chrono::high_resolution_clock::now();
	if( !IsWrappedPos ) computeWrap();
	// construct the cells
	// set maximum neighbour by estimating the maximum number of atom in a rc x rc x rc cube
	double rc_squared = pow(rc,2.);
	double zeronum = 1e-6;
	unsigned int nbCellX, nbCellY, nbCellZ;
	int NeighCellX, NeighCellY, NeighCellZ;
	const int bar_length = 30;
	double CellSizeX, CellSizeY, CellSizeZ, d_squared;
        nbCellX = floor(this->H1[0]/rc);
        if( nbCellX == 0 ){
                nbCellX = 1;
                CellSizeX = this->H1[0];
                NeighCellX = ceil(rc/this->H1[0]);
        }else{
                CellSizeX = this->H1[0]/nbCellX;
                NeighCellX = 1;
        }
        nbCellY = floor(this->H2[1]/rc);
        if( nbCellY == 0 ){
                nbCellY = 1;
                CellSizeY = this->H2[1];
                NeighCellY = ceil(rc/this->H2[1]);
        }else{
                CellSizeY = this->H2[1]/nbCellY;
                NeighCellY = 1;
        }
        nbCellZ = floor(this->H3[2]/rc);
        if( nbCellZ == 0 ){
                nbCellZ = 1;
                CellSizeZ = this->H3[2];
                NeighCellZ = ceil(rc/this->H3[2]);
        }else{
                CellSizeZ = this->H3[2]/nbCellZ;
                NeighCellZ = 1;
        }
	vector<vector<unsigned int>> Cells(nbCellX*nbCellY*nbCellZ);
	if( this->IsTilted ){
		// compute plane equations (as everywhere the only tilts considered are xy, yz, xz)
		double H1H3_z;
	        if(fabs(this->H3[2]) > 1e-6) H1H3_z = -(this->H3[1]/this->H3[2]);
		else H1H3_z = 0.;
		double H2H3_y = -(this->H2[0]/this->H2[1]); 
		double H2H3_z = ((this->H2[0]*this->H3[1]/this->H2[1])-this->H3[0])/this->H3[2];
		double planeH1H3, planeH2H3;
		#pragma omp parallel for private(planeH1H3,planeH2H3)
		for(unsigned int i=0;i<nbCellX;i++){
        	        for(unsigned int j=0;j<nbCellY;j++){
        	                for(unsigned int k=0;k<nbCellZ;k++){
        	                        for(unsigned at=0;at<this->nbAtom;at++){
						planeH1H3 = this->WrappedPos[at].y+this->WrappedPos[at].z*H1H3_z;	
						planeH2H3 = this->WrappedPos[at].x+this->WrappedPos[at].y*H2H3_y+this->WrappedPos[at].z*H2H3_z;	
        	                        	if( (i == nbCellX-1) && (j == nbCellY-1) && (k == nbCellZ-1) ){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<=(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<=(j+1)*CellSizeY) &&  (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<=(k+1)*CellSizeZ) ) Cells[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                        	}else if( (i == nbCellX-1) && (j == nbCellY-1) ){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<=(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<=(j+1)*CellSizeY) &&  (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<(k+1)*CellSizeZ) ) Cells[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                        	}else if( (j == nbCellY-1) && (k == nbCellZ-1) ){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<=(j+1)*CellSizeY) &&  (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<=(k+1)*CellSizeZ) ) Cells[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                        	}else if( (i == nbCellX-1) && (k == nbCellZ-1) ){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<=(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<(j+1)*CellSizeY) &&  (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<=(k+1)*CellSizeZ) ) Cells[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                        	}else if(i == nbCellX-1){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<=(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<(j+1)*CellSizeY) &&  (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<(k+1)*CellSizeZ) ) Cells[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                        	}else if(j == nbCellY-1){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<=(j+1)*CellSizeY) &&  (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<(k+1)*CellSizeZ) ) Cells[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                        	}else if(k == nbCellZ-1){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<(j+1)*CellSizeY) &&  (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<=(k+1)*CellSizeZ) ) Cells[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
						}else{
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<(j+1)*CellSizeY) &&  (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<(k+1)*CellSizeZ) ) Cells[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
						}
					}
				}
			}
		}
	}else{
		#pragma omp parallel for
		for(unsigned int i=0;i<nbCellX;i++){
        	        for(unsigned int j=0;j<nbCellY;j++){
        	                for(unsigned int k=0;k<nbCellZ;k++){
        	                        if( (i == nbCellX-1) && (j == nbCellY-1) && (k == nbCellZ-1) ){
        	                                for(unsigned at=0;at<this->nbAtom;at++){
        	                                        if( (this->WrappedPos[at].x>=i*CellSizeX) && (this->WrappedPos[at].x<=(i+1)*CellSizeX) && (this->WrappedPos[at].y>=j*CellSizeY) && (this->WrappedPos[at].y<=(j+1)*CellSizeY) && (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<=(k+1)*CellSizeZ) ) Cells[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                                }
        	                        }else if( (i == nbCellX-1) && (j == nbCellY-1) ){
        	                                for(unsigned at=0;at<this->nbAtom;at++){
        	                                        if( (this->WrappedPos[at].x>=i*CellSizeX) && (this->WrappedPos[at].x<=(i+1)*CellSizeX) && (this->WrappedPos[at].y>=j*CellSizeY) && (this->WrappedPos[at].y<=(j+1)*CellSizeY) && (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<(k+1)*CellSizeZ) ) Cells[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                                }
        	                        }else if( (j == nbCellY-1) && (k == nbCellZ-1) ){
        	                                for(unsigned at=0;at<this->nbAtom;at++){
        	                                        if( (this->WrappedPos[at].x>=i*CellSizeX) && (this->WrappedPos[at].x<(i+1)*CellSizeX) && (this->WrappedPos[at].y>=j*CellSizeY) && (this->WrappedPos[at].y<=(j+1)*CellSizeY) && (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<=(k+1)*CellSizeZ) ) Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
        	                                }
        	                        }else if( (i == nbCellX-1) && (k == nbCellZ-1) ){
        	                                for(unsigned at=0;at<this->nbAtom;at++){
        	                                        if( (this->WrappedPos[at].x>=i*CellSizeX) && (this->WrappedPos[at].x<=(i+1)*CellSizeX) && (this->WrappedPos[at].y>=j*CellSizeY) && (this->WrappedPos[at].y<(j+1)*CellSizeY) && (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<=(k+1)*CellSizeZ) ) Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
        	                                }
        	                        }else if(i == nbCellX-1){
        	                                for(unsigned at=0;at<this->nbAtom;at++){
        	                                        if( (this->WrappedPos[at].x>=i*CellSizeX) && (this->WrappedPos[at].x<=(i+1)*CellSizeX) && (this->WrappedPos[at].y>=j*CellSizeY) && (this->WrappedPos[at].y<(j+1)*CellSizeY) && (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<(k+1)*CellSizeZ) ) Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
        	                                }
        	                        }else if(j == nbCellY-1){
        	                                for(unsigned at=0;at<this->nbAtom;at++){
        	                                        if( (this->WrappedPos[at].x>=i*CellSizeX) && (this->WrappedPos[at].x<(i+1)*CellSizeX) && (this->WrappedPos[at].y>=j*CellSizeY) && (this->WrappedPos[at].y<=(j+1)*CellSizeY) && (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<(k+1)*CellSizeZ) ) Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
        	                                }
        	                        }else if(k == nbCellZ-1){
        	                                for(unsigned at=0;at<this->nbAtom;at++){
        	                                        if( (this->WrappedPos[at].x>=i*CellSizeX) && (this->WrappedPos[at].x<(i+1)*CellSizeX) && (this->WrappedPos[at].y>=j*CellSizeY) && (this->WrappedPos[at].y<(j+1)*CellSizeY) && (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<=(k+1)*CellSizeZ) ) Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
        	                                }
        	                        }else{
        	                                for(unsigned at=0;at<this->nbAtom;at++){
        	                                        if( (this->WrappedPos[at].x>=i*CellSizeX) && (this->WrappedPos[at].x<(i+1)*CellSizeX) && (this->WrappedPos[at].y>=j*CellSizeY) && (this->WrappedPos[at].y<(j+1)*CellSizeY) && (this->WrappedPos[at].z>=k*CellSizeZ) && (this->WrappedPos[at].z<(k+1)*CellSizeZ) ) Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
        	                                }
        	                        }
        	                }
        	        }
        	}
	}
	this->nbMaxN = Cells[0].size();
	unsigned int nbAt_test = Cells[0].size();
	for(unsigned int i=1;i<nbCellX*nbCellY*nbCellZ;i++){
		nbAt_test += Cells[i].size();
		if( Cells[i].size() > this->nbMaxN ) this->nbMaxN = Cells[i].size();
	}
	if( nbAt_test != this->nbAtom ) cout << "We miss atoms during cell list" << endl;
	this->nbMaxN *= (int) (2.*4.*M_PI*pow(rc,3.)/(3.*CellSizeX*CellSizeY*CellSizeZ)); // 2. is a security factor TODO put it in FixedParameters
	if( this->IsNeighbours ){
		delete[] this->Neighbours;
		delete[] this->CLNeighbours;
	}
	this->Neighbours = new unsigned int[(this->nbMaxN+1)*this->nbAtom];
	this->CLNeighbours = new int[(this->nbMaxN*3)*this->nbAtom]; // contain the periodic condition (Nclx, Ncly, Nclz) applied for atom to be a neighbour
	// Perform neighbour research
	double xpos,ypos,zpos;
	int ibx, jby, kbz;
	int Nclx, Ncly, Nclz;
	double prog=0.;
	unsigned int countN = 0;
	unsigned int currentId, currentId2;
	//cout << "Performing neighbour research" << endl;
	//cout << "\r[" << string(bar_length*prog,'X') << string(bar_length*(1-prog),'-') << "] " << setprecision(3) << 100*prog << "%";
	#pragma omp parallel for private(countN,currentId,xpos,ypos,zpos,Nclx,Ncly,Nclz,ibx,jby,kbz,currentId2,d_squared)
	for(unsigned int i=0;i<nbCellX;i++){
		for(unsigned int j=0;j<nbCellY;j++){
			for(unsigned int k=0;k<nbCellZ;k++){
				//prog = double(i*nbCellZ*nbCellY+j*nbCellZ+k)/double(nbCellX*nbCellY*nbCellZ);
				//cout << "\r[" << string(floor(bar_length*prog),'X') << string(ceil(bar_length*(1-prog)),'-') << "] " << setprecision(3) << 100*prog << "%";
				for(unsigned int at1 = 0; at1<Cells[i*nbCellZ*nbCellY+j*nbCellZ+k].size(); at1++){
				        countN = 0;
					currentId = Cells[i*nbCellZ*nbCellY+j*nbCellZ+k][at1];
					this->Neighbours[currentId*(this->nbMaxN+1)] = 0; // initialize to zero the neighbour counters
					xpos = this->WrappedPos[currentId].x;
					ypos = this->WrappedPos[currentId].y;
					zpos = this->WrappedPos[currentId].z;
					for(int bx=-NeighCellX;bx<NeighCellX+1;bx++){
						for(int by=-NeighCellY;by<NeighCellY+1;by++){
							for(int bz=-NeighCellZ;bz<NeighCellZ+1;bz++){
								// Search using periodic BC which cell use in case of border cell and what CL applied to ion pos
								// for a non border cell =>
								Nclx = 0;
								Ncly = 0;
								Nclz = 0;
								ibx = bx+i;
								jby = by+j;
								kbz = bz+k;
								// border cells
								if( (int) i >= ((int) nbCellX)-NeighCellX && bx >= 1 ){
									ibx = 0;
									Nclx = bx;
								}else if( (int) i <= NeighCellX-1 && bx <= -1 ){
									ibx = nbCellX-1;
									Nclx = bx;
								}
								if( (int) j >= ((int) nbCellY)-NeighCellY && by >= 1 ){
									jby = 0;
									Ncly = by;
								}else if( (int) j <= NeighCellY-1  && by <= -1 ){
									jby = nbCellY-1;
									Ncly = by;
								}
								if( (int) k >= ((int) nbCellZ)-NeighCellZ && bz >= 1 ){
									kbz = 0;
									Nclz = bz;
								}else if( (int) k <= NeighCellZ-1 && bz <= -1 ){
									kbz = nbCellZ-1;
									Nclz = bz;
								}
								for(unsigned int at2=0;at2<Cells[ibx*nbCellZ*nbCellY+jby*nbCellZ+kbz].size();at2++){
									currentId2 = Cells[ibx*nbCellY*nbCellZ+jby*nbCellZ+kbz][at2];
									d_squared = pow(this->WrappedPos[currentId2].x-xpos+Nclx*H1[0]+Ncly*H2[0]+Nclz*H3[0],2.)+pow(this->WrappedPos[currentId2].y-ypos+Nclx*H1[1]+Ncly*H2[1]+Nclz*H3[1],2.)+pow(this->WrappedPos[currentId2].z-zpos+Nclx*H1[2]+Ncly*H2[2]+Nclz*H3[2],2.);
									if( d_squared > zeronum && d_squared < rc_squared ){
										long unsigned int id1 = currentId*this->nbMaxN*3+countN*3;
										long unsigned int id2 = id1+1;
										long unsigned int id3 = id1+2;
										this->Neighbours[currentId*(this->nbMaxN+1)] += 1; // add this neighbours to the neighbour count
										this->Neighbours[currentId*(this->nbMaxN+1)+countN+1] = currentId2; // add this neighbours to the neighbour list 
										this->CLNeighbours[id1] = Nclx; // store the cl used for this neighbour
										this->CLNeighbours[id2] = Ncly; // store the cl used for this neighbour
										this->CLNeighbours[id3] = Nclz; // store the cl used for this neighbour
										countN += 1;
									}
								}
							}
						}
						if( countN >= this->nbMaxN ) cout << "Warning the number of found neighbour exceed the maximum number of neighbour allowed" << endl;
					}
				}
			}
		}
	}
	IsNeighbours = true;
	this->current_rc_neigh = rc;
	cout << "Done !" << endl;
        //auto end = chrono::high_resolution_clock::now();
        //auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
        //cout << "Execution time : " << duration.count() << endl;
	vector<vector<unsigned int>>().swap(Cells); // free the memory allocated for cells
	return this->nbMaxN;
}

// Searching neighbours using cell list algorithm restricted for
unsigned int AtomicSystem::searchNeighbours_restricted(const double& rc, const vector<unsigned int> & IndexToSearch, const vector<unsigned int> & IndexForSearch){
	unsigned int nbToSearch = IndexToSearch.size();
	unsigned int nbForSearch = IndexForSearch.size();
	//cout << "Performing neighbour research" << endl;
	//auto start = chrono::high_resolution_clock::now();
	if( !IsWrappedPos ) computeWrap();
	// construct the cells
	// set maximum neighbour by estimating the maximum number of atom in a rc x rc x rc cube
	double rc_squared = pow(rc,2.);
	double zeronum = 1e-6;
	unsigned int nbCellX, nbCellY, nbCellZ;
	int NeighCellX, NeighCellY, NeighCellZ;
	const int bar_length = 30;
	double CellSizeX, CellSizeY, CellSizeZ, d_squared;
        nbCellX = floor(this->H1[0]/rc);
        if( nbCellX == 0 ){
                nbCellX = 1;
                CellSizeX = this->H1[0];
                NeighCellX = ceil(rc/this->H1[0]);
        }else{
                CellSizeX = this->H1[0]/nbCellX;
                NeighCellX = 1;
        }
        nbCellY = floor(this->H2[1]/rc);
        if( nbCellY == 0 ){
                nbCellY = 1;
                CellSizeY = this->H2[1];
                NeighCellY = ceil(rc/this->H2[1]);
        }else{
                CellSizeY = this->H2[1]/nbCellY;
                NeighCellY = 1;
        }
        nbCellZ = floor(this->H3[2]/rc);
        if( nbCellZ == 0 ){
                nbCellZ = 1;
                CellSizeZ = this->H3[2];
                NeighCellZ = ceil(rc/this->H3[2]);
        }else{
                CellSizeZ = this->H3[2]/nbCellZ;
                NeighCellZ = 1;
        }
	vector<vector<unsigned int>> Cells_ToSearch(nbCellX*nbCellY*nbCellZ);
	vector<vector<unsigned int>> Cells_ForSearch(nbCellX*nbCellY*nbCellZ);
	if( this->IsTilted ){
		// compute plane equations (as everywhere the only tilts considered are xy, yz, xz)
		double H1H3_z;
	        if(fabs(this->H3[2]) > 1e-6) H1H3_z = -(this->H3[1]/this->H3[2]);
		else H1H3_z = 0.;
		double H2H3_y = -(this->H2[0]/this->H2[1]); 
		double H2H3_z = ((this->H2[0]*this->H3[1]/this->H2[1])-this->H3[0])/this->H3[2];
		double planeH1H3, planeH2H3;
		//#pragma omp parallel for private(planeH1H3,planeH2H3)
		for(unsigned int i=0;i<nbCellX;i++){
        	        for(unsigned int j=0;j<nbCellY;j++){
        	                for(unsigned int k=0;k<nbCellZ;k++){
        	                        for(unsigned at=0;at<nbToSearch;at++){
						planeH1H3 = this->WrappedPos[IndexToSearch[at]].y+this->WrappedPos[IndexToSearch[at]].z*H1H3_z;	
						planeH2H3 = this->WrappedPos[IndexToSearch[at]].x+this->WrappedPos[IndexToSearch[at]].y*H2H3_y+this->WrappedPos[IndexToSearch[at]].z*H2H3_z;	
        	                        	if( (i == nbCellX-1) && (j == nbCellY-1) && (k == nbCellZ-1) ){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<=(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<=(j+1)*CellSizeY) &&  (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                        	}else if( (i == nbCellX-1) && (j == nbCellY-1) ){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<=(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<=(j+1)*CellSizeY) &&  (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                        	}else if( (j == nbCellY-1) && (k == nbCellZ-1) ){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<=(j+1)*CellSizeY) &&  (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                        	}else if( (i == nbCellX-1) && (k == nbCellZ-1) ){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<=(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<(j+1)*CellSizeY) &&  (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                        	}else if(i == nbCellX-1){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<=(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<(j+1)*CellSizeY) &&  (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                        	}else if(j == nbCellY-1){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<=(j+1)*CellSizeY) &&  (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                        	}else if(k == nbCellZ-1){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<(j+1)*CellSizeY) &&  (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
						}else{
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<(j+1)*CellSizeY) &&  (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
						}
					}
        	                        for(unsigned at=0;at<nbForSearch;at++){
						planeH1H3 = this->WrappedPos[IndexForSearch[at]].y+this->WrappedPos[IndexForSearch[at]].z*H1H3_z;	
						planeH2H3 = this->WrappedPos[IndexForSearch[at]].x+this->WrappedPos[IndexForSearch[at]].y*H2H3_y+this->WrappedPos[IndexForSearch[at]].z*H2H3_z;	
        	                        	if( (i == nbCellX-1) && (j == nbCellY-1) && (k == nbCellZ-1) ){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<=(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<=(j+1)*CellSizeY) &&  (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                        	}else if( (i == nbCellX-1) && (j == nbCellY-1) ){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<=(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<=(j+1)*CellSizeY) &&  (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                        	}else if( (j == nbCellY-1) && (k == nbCellZ-1) ){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<=(j+1)*CellSizeY) &&  (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                        	}else if( (i == nbCellX-1) && (k == nbCellZ-1) ){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<=(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<(j+1)*CellSizeY) &&  (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                        	}else if(i == nbCellX-1){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<=(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<(j+1)*CellSizeY) &&  (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                        	}else if(j == nbCellY-1){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<=(j+1)*CellSizeY) &&  (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                        	}else if(k == nbCellZ-1){
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<(j+1)*CellSizeY) &&  (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(IndexForSearch[at]);
						}else{
        	                                        if( (planeH2H3>=i*CellSizeX) && (planeH2H3<(i+1)*CellSizeX) && (planeH1H3>=j*CellSizeY) && (planeH1H3<(j+1)*CellSizeY) &&  (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(IndexForSearch[at]);
						}
					}
				}
			}
		}
	}else{
		//#pragma omp parallel for
		for(unsigned int i=0;i<nbCellX;i++){
        	        for(unsigned int j=0;j<nbCellY;j++){
        	                for(unsigned int k=0;k<nbCellZ;k++){
        	                        if( (i == nbCellX-1) && (j == nbCellY-1) && (k == nbCellZ-1) ){
        	                                for(unsigned at=0;at<nbToSearch;at++){
        	                                        if( (this->WrappedPos[IndexToSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].x<=(i+1)*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].y<=(j+1)*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                                }
        	                                for(unsigned at=0;at<nbForSearch;at++){
        	                                        if( (this->WrappedPos[IndexForSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].x<=(i+1)*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].y<=(j+1)*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                                }
        	                        }else if( (i == nbCellX-1) && (j == nbCellY-1) ){
        	                                for(unsigned at=0;at<nbToSearch;at++){
        	                                        if( (this->WrappedPos[IndexToSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].x<=(i+1)*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].y<=(j+1)*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(at);
        	                                }
        	                                for(unsigned at=0;at<nbForSearch;at++){
        	                                        if( (this->WrappedPos[IndexForSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].x<=(i+1)*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].y<=(j+1)*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellY*nbCellZ+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                                }
        	                        }else if( (j == nbCellY-1) && (k == nbCellZ-1) ){
        	                                for(unsigned at=0;at<nbToSearch;at++){
        	                                        if( (this->WrappedPos[IndexToSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].x<(i+1)*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].y<=(j+1)*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
        	                                }
        	                                for(unsigned at=0;at<nbForSearch;at++){
        	                                        if( (this->WrappedPos[IndexForSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].x<(i+1)*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].y<=(j+1)*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                                }
        	                        }else if( (i == nbCellX-1) && (k == nbCellZ-1) ){
        	                                for(unsigned at=0;at<nbToSearch;at++){
        	                                        if( (this->WrappedPos[IndexToSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].x<=(i+1)*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].y<(j+1)*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
        	                                }
        	                                for(unsigned at=0;at<nbForSearch;at++){
        	                                        if( (this->WrappedPos[IndexForSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].x<=(i+1)*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].y<(j+1)*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                                }
        	                        }else if(i == nbCellX-1){
        	                                for(unsigned at=0;at<nbToSearch;at++){
        	                                        if( (this->WrappedPos[IndexToSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].x<=(i+1)*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].y<(j+1)*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
        	                                }
        	                                for(unsigned at=0;at<nbForSearch;at++){
        	                                        if( (this->WrappedPos[IndexForSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].x<=(i+1)*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].y<(j+1)*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                                }
        	                        }else if(j == nbCellY-1){
        	                                for(unsigned at=0;at<nbToSearch;at++){
        	                                        if( (this->WrappedPos[IndexToSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].x<(i+1)*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].y<=(j+1)*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
        	                                }
        	                                for(unsigned at=0;at<nbForSearch;at++){
        	                                        if( (this->WrappedPos[IndexForSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].x<(i+1)*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].y<=(j+1)*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                                }
        	                        }else if(k == nbCellZ-1){
        	                                for(unsigned at=0;at<nbToSearch;at++){
        	                                        if( (this->WrappedPos[IndexToSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].x<(i+1)*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].y<(j+1)*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
        	                                }
        	                                for(unsigned at=0;at<nbForSearch;at++){
        	                                        if( (this->WrappedPos[IndexForSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].x<(i+1)*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].y<(j+1)*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<=(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                                }
        	                        }else{
        	                                for(unsigned at=0;at<nbToSearch;at++){
        	                                        if( (this->WrappedPos[IndexToSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].x<(i+1)*CellSizeX) && (this->WrappedPos[IndexToSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].y<(j+1)*CellSizeY) && (this->WrappedPos[IndexToSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexToSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ToSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(at);
        	                                }
        	                                for(unsigned at=0;at<nbForSearch;at++){
        	                                        if( (this->WrappedPos[IndexForSearch[at]].x>=i*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].x<(i+1)*CellSizeX) && (this->WrappedPos[IndexForSearch[at]].y>=j*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].y<(j+1)*CellSizeY) && (this->WrappedPos[IndexForSearch[at]].z>=k*CellSizeZ) && (this->WrappedPos[IndexForSearch[at]].z<(k+1)*CellSizeZ) ) Cells_ForSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].push_back(IndexForSearch[at]);
        	                                }
        	                        }
        	                }
        	        }
        	}
	}
	//cout << "Cells done"<< endl;
	this->nbMaxN = Cells_ForSearch[0].size();
	unsigned int nbAt_test = Cells_ForSearch[0].size();
	unsigned int nbAt_test_ToSearch = Cells_ToSearch[0].size();
	for(unsigned int i=1;i<nbCellX*nbCellY*nbCellZ;i++){
		nbAt_test += Cells_ForSearch[i].size();
		nbAt_test_ToSearch += Cells_ToSearch[i].size();
		if( Cells_ForSearch[i].size() > this->nbMaxN ) this->nbMaxN = Cells_ForSearch[i].size();
	}
	if( nbAt_test != nbForSearch ) cout << "We miss atoms during cell list" << endl;
	if( nbAt_test_ToSearch != nbToSearch ) cout << "We miss atoms during cell list" << endl;
	this->nbMaxN *= (int) (2.*4.*M_PI*pow(rc,3.)/(3.*CellSizeX*CellSizeY*CellSizeZ)); // 2. is a security factor
	if( this->IsNeighbours ){
		delete[] this->Neighbours;
		delete[] this->CLNeighbours;
	} // TODO maybe issue here if we delete the var we may need to reclare them ?
	unsigned long int SizeNeigh = ((unsigned long int) (this->nbMaxN) + 1)*((unsigned long int) nbToSearch);
	unsigned long int SizeCLNeigh = ((unsigned long int) (this->nbMaxN))*((unsigned long int) nbToSearch)*((unsigned long int) 3);
	this->Neighbours = new unsigned int[SizeNeigh];
	this->CLNeighbours = new int[SizeCLNeigh]; // contain the periodic condition (Nclx, Ncly, Nclz) applied for atom to be a neighbour
	// Perform neighbour research
	double xpos,ypos,zpos;
	int ibx, jby, kbz;
	int Nclx, Ncly, Nclz;
	double prog=0.;
	unsigned long int countN = 0;
	unsigned long int currentId, currentId2;
	unsigned long int one_l = 1;
	unsigned long int two_l = 1;
	unsigned long int three_l = 3;
	//cout << "done ! number of neighbour max = " << this->nbMaxN << endl;
	//cout << "Performing neighbour research" << endl;
	//cout << "\r[" << string(bar_length*prog,'X') << string(bar_length*(1-prog),'-') << "] " << setprecision(3) << 100*prog << "%";
	//#pragma omp parallel for private(countN,currentId,xpos,ypos,zpos,Nclx,Ncly,Nclz,ibx,jby,kbz,currentId2,d_squared)
	for(unsigned int i=0;i<nbCellX;i++){
		for(unsigned int j=0;j<nbCellY;j++){
			for(unsigned int k=0;k<nbCellZ;k++){
				//prog = double(i*nbCellZ*nbCellY+j*nbCellZ+k)/double(nbCellX*nbCellY*nbCellZ);
				//cout << "\r[" << string(floor(bar_length*prog),'X') << string(ceil(bar_length*(1-prog)),'-') << "] " << setprecision(3) << 100*prog << "%";
				for(unsigned int at1 = 0; at1<Cells_ToSearch[i*nbCellZ*nbCellY+j*nbCellZ+k].size(); at1++){
					bool Over = false;
				        countN = 0;
					currentId = Cells_ToSearch[i*nbCellZ*nbCellY+j*nbCellZ+k][at1];
	if( currentId > nbToSearch ) cout << "Issue HERE" << endl;
					this->Neighbours[currentId*(this->nbMaxN+1)] = 0; // initialize to zero the neighbour counters
					xpos = this->WrappedPos[IndexToSearch[currentId]].x;
					ypos = this->WrappedPos[IndexToSearch[currentId]].y;
					zpos = this->WrappedPos[IndexToSearch[currentId]].z;
					for(int bx=-NeighCellX;bx<NeighCellX+1;bx++){
						for(int by=-NeighCellY;by<NeighCellY+1;by++){
							for(int bz=-NeighCellZ;bz<NeighCellZ+1;bz++){
								// Search using periodic BC which cell use in case of border cell and what CL applied to ion pos
								// for a non border cell =>
								Nclx = 0;
								Ncly = 0;
								Nclz = 0;
								ibx = bx+i;
								jby = by+j;
								kbz = bz+k;
								// border cells
								if( (int) i >= ((int) nbCellX)-NeighCellX && bx >= 1 ){
									ibx = 0;
									Nclx = bx;
								}else if( (int) i <= NeighCellX-1 && bx <= -1 ){
									ibx = nbCellX-1;
									Nclx = bx;
								}
								if( (int) j >= ((int) nbCellY)-NeighCellY && by >= 1 ){
									jby = 0;
									Ncly = by;
								}else if( (int) j <= NeighCellY-1  && by <= -1 ){
									jby = nbCellY-1;
									Ncly = by;
								}
								if( (int) k >= ((int) nbCellZ)-NeighCellZ && bz >= 1 ){
									kbz = 0;
									Nclz = bz;
								}else if( (int) k <= NeighCellZ-1 && bz <= -1 ){
									kbz = nbCellZ-1;
									Nclz = bz;
								}
								for(unsigned int at2=0;at2<Cells_ForSearch[ibx*nbCellZ*nbCellY+jby*nbCellZ+kbz].size();at2++){
									currentId2 = Cells_ForSearch[ibx*nbCellY*nbCellZ+jby*nbCellZ+kbz][at2];
									d_squared = pow(this->WrappedPos[currentId2].x-xpos+Nclx*H1[0]+Ncly*H2[0]+Nclz*H3[0],2.)+pow(this->WrappedPos[currentId2].y-ypos+Nclx*H1[1]+Ncly*H2[1]+Nclz*H3[1],2.)+pow(this->WrappedPos[currentId2].z-zpos+Nclx*H1[2]+Ncly*H2[2]+Nclz*H3[2],2.);
									if( d_squared > zeronum && d_squared < rc_squared ){
										this->Neighbours[currentId*(this->nbMaxN+one_l)] += 1; // add this neighbours to the neighbour count
										this->Neighbours[currentId*(this->nbMaxN+one_l)+countN+one_l] = currentId2; // add this neighbours to the neighbour list 
										this->CLNeighbours[(this->nbMaxN)*currentId*three_l+countN*three_l] = Nclx; // store the cl used for this neighbour
										this->CLNeighbours[(this->nbMaxN)*currentId*three_l+countN*three_l+one_l] = Ncly; // store the cl used for this neighbour
										this->CLNeighbours[(this->nbMaxN)*currentId*three_l+countN*three_l+two_l] = Nclz; // store the cl used for this neighbour
										countN += 1;
										if( countN >= this->nbMaxN ){
											Over = true;
											cout << "Warning the number of found neighbour exceed the maximum number of neighbour allowed" << endl;
											break;
										}
									}
								}
								if( Over ) break;
							}
							if( Over ) break;
						}
						if( Over ) break;
					}
				}
			}
		}
	}
	IsNeighbours = true;
	this->current_rc_neigh = rc;
	//cout << " done !" << endl;
        //auto end = chrono::high_resolution_clock::now();
        //auto duration = chrono::duration_cast<chrono::microseconds>(end - start);
        //cout << "Execution time : " << duration.count() << endl;
	vector<vector<unsigned int>>().swap(Cells_ToSearch); // free memory allocated for cells
	vector<vector<unsigned int>>().swap(Cells_ForSearch);
	return this->nbMaxN;
}


void AtomicSystem::setAux(const double* aux, const string& AuxName){
	IsSetAux = true;
	for(unsigned int i=0;i<this->Aux_name.size();i++){
		if( AuxName == this->Aux_name[i] ){
			cout << "Auxiliary propery already stored, aborting" << endl;
			return;
		}
	}
	this->Aux.push_back(new double[this->nbAtom]);
	this->Aux_name.push_back(AuxName);
	this->Aux_size.push_back(1);
	for(unsigned int i=0;i<this->nbAtom;i++) Aux[Aux.size()-1][i] = aux[i];
}

void AtomicSystem::setAux_vec(const double* aux, const unsigned int size, const string& AuxName){
	IsSetAux = true;
	for(unsigned int i=0;i<this->Aux_name.size();i++){
		if( AuxName == this->Aux_name[i] ){
			cout << "Auxiliary propery already stored, aborting" << endl;
			return;
		}
	}
	this->Aux.push_back(new double[this->nbAtom*size]);
	this->Aux_name.push_back(AuxName);
	this->Aux_size.push_back(size);
	for(unsigned int i=0;i<this->nbAtom;i++){
		for(unsigned int j=0;j<size;j++){
			Aux[Aux.size()-1][i*size+j] = aux[i*size+j];
		}
	}
}

void AtomicSystem::modifyAux_vec(const double* aux, const string& AuxName){
	unsigned int aux_ind,size;
	bool ok = false;
	for(unsigned int i=0;i<this->Aux_name.size();i++){
		if( AuxName == this->Aux_name[i] ){
			aux_ind = i;
			ok = true;
			break;
		}
	}
	if( ok ){
		size = this->Aux_size[aux_ind];
		for(unsigned int i=0;i<this->nbAtom;i++){
			for(unsigned int j=0;j<size;j++){
				Aux[aux_ind][i*size+j] = aux[i*size+j];
			}
		}
	}else{
		cout << "Auxiliary property to modify does not exist, aborting" << endl;
	}
}

void AtomicSystem::modifyAux_vec(const unsigned int* aux, const string& AuxName){
	unsigned int aux_ind,size;
	bool ok = false;
	for(unsigned int i=0;i<this->Aux_name.size();i++){
		if( AuxName == this->Aux_name[i] ){
			aux_ind = i;
			ok = true;
			break;
		}
	}
	if( ok ){
		size = this->Aux_size[aux_ind];
		for(unsigned int i=0;i<this->nbAtom;i++){
			for(unsigned int j=0;j<size;j++){
				Aux[aux_ind][i*size+j] = (double) aux[i*size+j];
			}
		}
	}else{
		cout << "Auxiliary property to modify does not exist, aborting" << endl;
	}
}
void AtomicSystem::setAux(const int* aux, const string& AuxName){
	IsSetAux = true;
	this->Aux.push_back(new double[this->nbAtom]);
	this->Aux_name.push_back(AuxName);
	this->Aux_size.push_back(1);
	for(unsigned int i=0;i<this->nbAtom;i++) Aux[Aux.size()-1][i] = (double) aux[i];
}

void AtomicSystem::setAux(const unsigned int* aux, const string& AuxName){
	IsSetAux = true;
	this->Aux.push_back(new double[this->nbAtom]);
	this->Aux_name.push_back(AuxName);
	this->Aux_size.push_back(1);
	for(unsigned int i=0;i<this->nbAtom;i++) Aux[Aux.size()-1][i] = (double) aux[i];
}

void AtomicSystem::printSystem(const string& filename){
	string ext=filename.substr(filename.find_last_of(".") + 1);
	if( ext == "lmp" ){
		this->print_lmp(filename);
	}else if( ext == "cfg" || ext == "xsf" ){
		this->print_cfg(filename);
	}else{
		cerr << "The extension \"" << ext << "\" of the output file is unknown (possible extension : lmp, cfg(or xsf)), aborting the print" << endl;
	}
}

bool AtomicSystem::read_lmp_file(const string& filename){
	ifstream file(filename, ios::in);
	size_t pos_at, pos_x, pos_y, pos_z, pos_tilt, pos_attype, pos_Mass, pos_At, pos_vel, pos_bond, pos_bondtype, pos_angle, pos_angletype, pos_Bond, pos_Angle;
	unsigned int h_uint = 1e8;
	unsigned int line_Mass(h_uint), line_At(h_uint), line_vel(h_uint), buffer_uint, buffer_uint_1, buffer_uint_2, buffer_uint_3, buffer_uint_4, count(0), buffer_molid, line_Bond(h_uint), line_Angle(h_uint);
	double buffer_1, buffer_2, buffer_3, buffer_4;
	double xlo,xhi,ylo,yhi,zlo,zhi;
	string buffer_s, buffer_s_1, buffer_s_2, line;
	unsigned int ReadOk(0), NbAtRead(0), NbVelRead(0);
        unsigned int count_at(0), nbFields(0);
	if(file){
		while(getline(file,line)){
			// find number of atom
			pos_at=line.find("atoms");
			if(pos_at!=string::npos){
				istringstream text(line);
				text >> buffer_uint >> buffer_s;
				if( buffer_s == "atoms" ){
					this->nbAtom = buffer_uint;
					AtomList = new Atom[this->nbAtom];
					ReadOk++;
				}
			}
			// manage bonds and angles
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
				xhi = buffer_2;
				xlo = buffer_1;
				ReadOk++;
			}

			// find H2 vector
			pos_y=line.find("ylo yhi");
			if(pos_y!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				yhi = buffer_2;
				ylo = buffer_1;
				ReadOk++;
			}

			// find H3 vector
			pos_z=line.find("zlo zhi");
			if(pos_z!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				zhi = buffer_2;
				zlo = buffer_1;
				ReadOk++;
			}

			// find tilts
			pos_tilt=line.find("xy xz yz");
			if(pos_tilt!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2 >> buffer_3;
				this->H2[0] = buffer_1;
				this->H3[0] = buffer_2;
				this->H3[1] = buffer_3;
				this->IsTilted = true;
			}

			// find nb atom type
			pos_attype=line.find("atom types");
			if(pos_attype!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_s >> buffer_s_1;
				if( buffer_s == "atom" && buffer_s_1 == "types" ) this->nbAtomType = buffer_1;
			}

			// get lines where are the keywords Masses and Atoms to get atom type masses and positions
			pos_Mass=line.find("Masses");
			if(pos_Mass!=string::npos){
				istringstream text(line);
				text >> buffer_s;
				if( buffer_s == "Masses" ) line_Mass = count;
			}
			if( count > line_Mass+1 && count < line_Mass+2+this->nbAtomType ){
				istringstream text(line);
				text >> buffer_uint >> buffer_1 >> buffer_s_1 >> buffer_s;
				this->AtomMass[buffer_uint-1] = buffer_1;
				if( buffer_s_1 == "#" ){
					this->AtomType[buffer_uint-1] = buffer_s;
					this->IsElem = true;
				}
			}
			pos_At=line.find("Atoms");
			if(pos_At!=string::npos){
				istringstream text(line);
				text >> buffer_s_1 >> buffer_s_2 >> buffer_s;
				if( buffer_s_1 == "Atoms" ){
					if( buffer_s == "charge" ) this->IsCharge = true;
					else if( buffer_s == "full" ){
						this->IsMolId = true;
						MolId = new unsigned int[this->nbAtom];
						this->IsCharge = true;
					}
			       		line_At = count;
					ReadOk++;
				}
			}
			if( count == line_At+2 ){
				istringstream text(line);
				while( text >> buffer_s ) nbFields++;
				if( nbFields == 10 ){
					IsPeriodicArr = true;
					PeriodicArr = new int[nbAtom*3];
				}
			}	

			if( count > line_At+1 && count < line_At+2+this->nbAtom ){
				istringstream text(line);
				text >> buffer_uint;
				if( buffer_uint > this->nbAtom ){
					cerr << "An index of atom (" << buffer_uint << ") is higher than the total number of atom (" << this->nbAtom << "), aborting" << endl;
					exit(EXIT_FAILURE);
				}
				if( this->IsMolId ){
					text >> buffer_molid;
					MolId[buffer_uint-1] = buffer_molid;
				}
				text >> buffer_uint_1;
				if( this->IsCharge ){
					text >> buffer_1;
					this->AtomCharge[buffer_uint_1-1] = buffer_1;
				}
				text >> buffer_2 >> buffer_3 >> buffer_4;
				this->AtomList[buffer_uint-1].pos.x = buffer_2;
				this->AtomList[buffer_uint-1].pos.y = buffer_3;
				this->AtomList[buffer_uint-1].pos.z = buffer_4;
				this->AtomList[buffer_uint-1].type_uint = buffer_uint_1;
				if( IsPeriodicArr ) for(unsigned int d=0;d<3;d++) text >> PeriodicArr[(buffer_uint-1)*3+d];//buffer_int_1 >> buffer_int_2 >> buffer_int_2;
				NbAtRead++;
			}
			// Read bonds
			if( IsBond ){
				pos_Bond=line.find("Bonds");
				if(pos_Bond!=string::npos) line_Bond = count;
				if( count > line_Bond+1 && count < line_Bond+nbBonds+2 ){
					istringstream text(line);
					text >> buffer_uint >> buffer_uint_1 >> buffer_uint_2 >> buffer_uint_3;
					if( buffer_uint > nbBonds ){
						cerr << "An index of bond (" << buffer_uint << ") is higher than the total number of bond (" << nbBonds << "), aborting" << endl;
						exit(EXIT_FAILURE);
					}
					if( buffer_uint_2 > this->nbAtom || buffer_uint_3 > this->nbAtom ){
						cerr << "An index of atom in bond (" << buffer_uint << ") is higher than the total number of atom (" << this->nbAtom << "), aborting" << endl;
						exit(EXIT_FAILURE);
					}
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
					if( buffer_uint > nbAngles ){
						cerr << "An index of angle (" << buffer_uint << ") is higher than the total number of angle (" << nbAngles << "), aborting" << endl;
						exit(EXIT_FAILURE);
					}
					if( buffer_uint_2 > this->nbAtom || buffer_uint_3 > this->nbAtom || buffer_uint_4 > this->nbAtom ){
						cerr << "An index of atom in angle (" << buffer_uint << ") is higher than the total number of atom (" << this->nbAtom << "), aborting" << endl;
						exit(EXIT_FAILURE);
					}
					AngleType[buffer_uint-1] = buffer_uint_1;
					Angles[(buffer_uint-1)*3] = buffer_uint_2;
					Angles[(buffer_uint-1)*3+1] = buffer_uint_3;
					Angles[(buffer_uint-1)*3+2] = buffer_uint_4;
				}
			}

			pos_vel=line.find("Velocities");
			if(pos_vel!=string::npos){
				istringstream text(line);
				text >> buffer_s;
				this->Aux_size.push_back(3);
				this->Aux_name.push_back("Velocities");
				this->Aux.push_back(nullptr);
				unsigned int mys = Aux.size();
				Aux[mys-1] = new double[3*(this->nbAtom)];
				unsigned int totnum = nbAtom*3;
				for(unsigned int i=0;i<totnum;i++) Aux[mys-1][i] = 0.;
			       	line_vel = count;
				IsVel = true;
				IsSetAux = true;
			}
			count += 1;
			if( !file ) break;
		}
		file.close();
		// compute the cell vectors
		this->H1[0] = xhi-xlo;
		this->H2[1] = yhi-ylo;
		this->H3[2] = zhi-zlo;
		if( IsPeriodicArr ){
			for(unsigned int i=0;i<nbAtom;i++){
				this->AtomList[i].pos.x += H1[0]*PeriodicArr[i*3] + H2[0]*PeriodicArr[i*3+1] + H3[0]*PeriodicArr[i*3+2];
				this->AtomList[i].pos.y += H1[1]*PeriodicArr[i*3] + H2[1]*PeriodicArr[i*3+1] + H3[1]*PeriodicArr[i*3+2];
				this->AtomList[i].pos.z += H1[2]*PeriodicArr[i*3] + H2[2]*PeriodicArr[i*3+1] + H3[2]*PeriodicArr[i*3+2];
			}
		}
		ifstream file2;
		file2.open(filename.c_str(), ifstream::in);
		count = 0;
		if(file2){
			while(getline(file2,line)){
				if( IsVel && count > line_vel+1 && count < line_vel+2+this->nbAtom ){
					istringstream text2(line);
					text2 >> buffer_uint_1 >> buffer_1 >> buffer_2 >> buffer_3;	
					if( buffer_uint_1 > this->nbAtom ){
						cerr << "An index of velocity (" << buffer_uint_1 << ") is higher than the total number of atom (" << this->nbAtom << "), aborting" << endl;
						exit(EXIT_FAILURE);
					}
					this->Aux[0][(buffer_uint_1-1)*3] = buffer_1;
					this->Aux[0][(buffer_uint_1-1)*3+1] = buffer_2;
					this->Aux[0][(buffer_uint_1-1)*3+2] = buffer_3;
					NbVelRead++;
				}
				count++;
			}
		}
		file2.close();
	}else{
		cout << "The file " << filename << " cannot be openned" << endl;
		exit(EXIT_FAILURE);
		return false;
	}
	if( ReadOk == 5 ){
		if( NbAtRead == this->nbAtom ){
			if( IsVel && NbVelRead == this->nbAtom ) return true;
			else{
				if( IsVel ) cout << "warning ! the number of atomic velocities provided does not correspond to the number of atom" << endl;
				return true;
			}
		}else{
			cout << "the number of atom provided does not correspond to the true number of atom" << endl;
			exit(EXIT_FAILURE);
			return false;
		}
	}else{
		return false;
	}
}

bool AtomicSystem::read_other_cfg(const string& filename){
	ifstream file(filename, ios::in);
	if(file){
		unsigned int buffer_uint, count(0), nbAux, line_aux(1000), indId, count_at(0), current_type_uint, count_aux, current_ind;
		size_t pos_At, pos_H1_x, pos_H1_y, pos_H1_z, pos_H2_x, pos_H2_y, pos_H2_z, pos_H3_x, pos_H3_y, pos_H3_z, pos_nbAux;
		double buffer_1, buffer_2, buffer_3, buffer_4, buffer_5;
		bool IsId = false, TypeStored;
		string buffer_s, buffer_s_1, buffer_s_2, buffer_s_3, line;
		unsigned int ReadOk(0);
		this->nbAtomType = 0;
		while(file){
			getline(file,line);
			// find number of atom 
			pos_At=line.find("Number of particles");
			if(pos_At!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1 >> buffer_s_2 >> buffer_s_3 >> buffer_uint;
			       	this->nbAtom = buffer_uint;
				AtomList = new Atom[this->nbAtom];
				ReadOk++;
			}
			// find cell vectors
			pos_H1_x=line.find("H0(1,1)");
			if(pos_H1_x!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1 >> buffer_1;
				this->H1[0] = buffer_1;
				ReadOk++;
			}
			pos_H1_y=line.find("H0(1,2)");
			if(pos_H1_y!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1 >> buffer_1;
				this->H1[1] = buffer_1;
				ReadOk++;
			}
			pos_H1_z=line.find("H0(1,3)");
			if(pos_H1_z!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1 >> buffer_1;
				this->H1[2] = buffer_1;
				ReadOk++;
			}
			pos_H2_x=line.find("H0(2,1)");
			if(pos_H2_x!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1 >> buffer_1;
				this->H2[0] = buffer_1;
				ReadOk++;
			}
			pos_H2_y=line.find("H0(2,2)");
			if(pos_H2_y!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1 >> buffer_1;
				this->H2[1] = buffer_1;
				ReadOk++;
			}
			pos_H2_z=line.find("H0(2,3)");
			if(pos_H2_z!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1 >> buffer_1;
				this->H2[2] = buffer_1;
				ReadOk++;
			}
			pos_H3_x=line.find("H0(3,1)");
			if(pos_H3_x!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1 >> buffer_1;
				this->H3[0] = buffer_1;
				ReadOk++;
			}
			pos_H3_y=line.find("H0(3,2)");
			if(pos_H3_y!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1 >> buffer_1;
				this->H3[1] = buffer_1;
				ReadOk++;
			}
			pos_H3_z=line.find("H0(3,3)");
			if(pos_H3_z!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1 >> buffer_1;
				this->H3[2] = buffer_1;
				ReadOk++;
			}
			pos_nbAux=line.find("entry_count");
			if(pos_nbAux!=string::npos){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1 >> nbAux;
				nbAux -= 3; // the three atomic coordinates
				line_aux = count;
				if( nbAux > 0 ) this->IsSetAux = true;
				ReadOk++;
			}
			if( ReadOk == 11 && count > line_aux && count <= line_aux+nbAux ){
				istringstream text(line);
				text >> buffer_s >> buffer_s_1 >> buffer_s_2;
				if( buffer_s_2 == "id" ){ // TODO : same with charge
					IsId = true;
					indId = Aux_name.size();
					this->Aux_size.push_back(1); //TODO for vector aux
				}else{
					this->Aux_name.push_back(buffer_s_2);
					this->Aux.push_back(new double[this->nbAtom]);
					this->Aux_size.push_back(1); //TODO for vector aux
				}
			}
			if( ReadOk == 11 && count > line_aux+nbAux ){
				istringstream text2(line);
				int nbCol = 0;
				do{
					string sub;
					text2 >> sub;
					if( sub.length() )
						++nbCol;
				}
				while ( text2 );
				if( nbCol > 1 ){
					istringstream text2(line);
					text2 >> buffer_1 >> buffer_2 >> buffer_3;
					count_aux = 0;
					for(unsigned int i=0;i<nbAux;i++){
						if( IsId ){
							if( i == indId ) text2 >> current_ind;
							else{
								//text2 >> Aux[count_aux][current_ind-1];
								text2 >> Aux[count_aux][count_at];
								count_aux += 1;
							}
						}else text2 >> Aux[i][count_at];
					}
					//if( IsId ){
					//	this->AtomList[current_ind-1].pos.x = buffer_1*H1[0]+buffer_2*H2[0]+buffer_3*H3[0];
					//	this->AtomList[current_ind-1].pos.y = buffer_1*H1[1]+buffer_2*H2[1]+buffer_3*H3[1];
					//	this->AtomList[current_ind-1].pos.z = buffer_1*H1[2]+buffer_2*H2[2]+buffer_3*H3[2];
					//	this->AtomList[current_ind-1].type_uint = current_type_uint;
					//}else{
						this->AtomList[count_at].pos.x = buffer_1*H1[0]+buffer_2*H2[0]+buffer_3*H3[0];
						this->AtomList[count_at].pos.y = buffer_1*H1[1]+buffer_2*H2[1]+buffer_3*H3[1];
						this->AtomList[count_at].pos.z = buffer_1*H1[2]+buffer_2*H2[2]+buffer_3*H3[2];
						this->AtomList[count_at].type_uint = current_type_uint;
					//}
					count_at += 1;
				} else {
					istringstream text(line);
					text >> buffer_1;
					getline(file,line);
					istringstream text1(line);
					text1 >> buffer_s;
					TypeStored = false;
					for(unsigned int i=0;i<this->nbAtomType;i++){
						if( buffer_s == this->AtomType[i] ){
							TypeStored = true;
							current_type_uint = i+1;
							break;
						}
					}
					if( !TypeStored ){
						this->AtomType[nbAtomType] = buffer_s;
						this->AtomMass[nbAtomType] = buffer_1;
						this->nbAtomType += 1;
						current_type_uint = this->nbAtomType;
					}
				}
			}
			count += 1;
		}
		if( H2[0] != 0 || H3[0] != 0 || H3[1] != 0 ) IsTilted = true;
		if( ReadOk == 11 ){
			if( nbAtom == count_at ){
				return true;
			}else{
				cout << "The number of atomic data does not correspond to the number of atom" << endl;
				exit(EXIT_FAILURE);
				return false;
			}
		}else{
			return false;
		}
	}else{
		cout << "The file " << filename << " cannot be openned" << endl;
		exit(EXIT_FAILURE);
		return false;
	}
}

bool AtomicSystem::read_cfg_file(const string& filename,const int &_timestep){
	ifstream file(filename, ios::in);
	unsigned int ReadOk(0);
	unsigned int NbAtRead = 0;
	if(file){
		unsigned int hi_number_uint = std::numeric_limits<unsigned int>::max() - 5e8; // 5e8 is now the maximum number of atoms
		unsigned int line_dt(hi_number_uint), line_At(hi_number_uint), line_H(hi_number_uint), line_at(hi_number_uint), final_line(hi_number_uint), buffer_uint, buffer_uint_1, count_H(0), count(0), nbAux(0), aux_count;
		size_t pos_dt, pos_At, pos_H, pos_charge, pos_at, pos_aux_vec, pos_elem, pos_typeuint, pos_id;
		double buffer_1, buffer_2, buffer_3, buffer_4, buffer_5;
		double xlo,xhi,ylo,yhi,zlo,zhi;
		string buffer_s, buffer_s_1, buffer_s_2, line, aux_name;
		bool other_cfg = false, IsReducedCoords = false;
		bool IsType_uint = false;
		bool IsId = false;
		bool GoodTimestep = true;
		if( _timestep != -1 ) GoodTimestep = false;
		unsigned int nbAux_norm=3;
		vector<string> befAuxNames;
		nbAtomType = 0;
		while(file){
			getline(file,line);
			// find timestep
			pos_dt=line.find("TIMESTEP");
			if( pos_dt!=string::npos ) line_dt = count;
			if( count == line_dt+1 ){
				istringstream text(line);
				text >> buffer_1;
			       	this->timestep = buffer_1;
				if( !GoodTimestep && _timestep == (int) timestep ) GoodTimestep = true;
			}

			// find number of atom
			pos_At=line.find("NUMBER OF ATOMS");
			if( pos_At!=string::npos ) line_At = count;
			if( count == line_At+1 && GoodTimestep ){
				istringstream text(line);
				text >> buffer_uint;
			       	this->nbAtom = buffer_uint;
				AtomList = new Atom[this->nbAtom];
				ReadOk++;
			}

			// find box vectors
			pos_H=line.find("BOX BOUNDS");
			if( pos_H!=string::npos && GoodTimestep ){
				line_H = count;
				istringstream text2(line);
				unsigned int nbCol = 0;
				do{
					string sub;
					text2 >> sub;
					if( sub.length() ) ++nbCol;
				}while( text2 );
				if( nbCol > 6 )	this->IsTilted = true;
				ReadOk++;
			}
			if( IsTilted && count > line_H && count < line_H+4 && GoodTimestep ){
				istringstream text(line);
				text >> buffer_1 >> buffer_2 >> buffer_3;
				if( count_H == 0 ){
					xlo = buffer_1;
					xhi = buffer_2;
					this->H2[0] = buffer_3;
				}else if( count_H == 1 ){
					ylo = buffer_1;
					yhi = buffer_2;
					this->H3[0] = buffer_3;
				}else if( count_H == 2 ){
					zlo = buffer_1;
					zhi = buffer_2;
					this->H3[2] = buffer_2-buffer_1;
					this->H3[1] = buffer_3;
				double arr[4] = {0.,this->H2[0],this->H3[0],this->H2[0]+this->H3[0]};
				this->H1[0] = xhi-xlo+this->MT->min(arr,4)-this->MT->max(arr,4);
				double arr_2[2] = {0.,this->H3[1]};
				this->H2[1] = yhi-ylo+this->MT->min(arr_2,2)-this->MT->max(arr_2,2);
				}
				count_H += 1;
			}else if( !IsTilted && count > line_H && count < line_H+4 && GoodTimestep ){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				if( count_H == 0 ){
					xlo = buffer_1;
					xhi = buffer_2;
					this->H1[0] = buffer_2-buffer_1;
				}else if( count_H == 1 ){
					ylo = buffer_1;
					yhi = buffer_2;
					this->H2[1] = buffer_2-buffer_1;
				}else if( count_H == 2 ){
					zlo = buffer_1;
					zhi = buffer_2;
					this->H3[2] = buffer_2-buffer_1;
				}
				count_H += 1;
			}

			// search if atom are charged
			pos_charge=line.find(" q");
			if( pos_charge!=string::npos && GoodTimestep ){
				this->IsCharge = true;
				nbAux_norm += 1;
			}
			pos_elem=line.find(" element");
			if( pos_elem!=string::npos && GoodTimestep ){
				IsElem = true;
				nbAux_norm += 1;
			}
			pos_typeuint=line.find(" type");
			if( pos_typeuint!=string::npos && GoodTimestep ){
				IsType_uint = true;
				nbAux_norm += 1;
			}
			pos_id=line.find(" id");
			if( pos_id!=string::npos && GoodTimestep ){
				IsId = true;
				nbAux_norm += 1;
			}
			// find and get atom positions
			pos_at=line.find("ITEM: ATOMS");
			if( pos_at!=string::npos && GoodTimestep ){
				ReadOk++;
				istringstream text(line);
				text >> buffer_s >> buffer_s;
				while(text >> buffer_s){
					nbAux += 1;
					if( nbAux > nbAux_norm ){
					       this->IsSetAux = true;
					       // search if auxiliary property is vector (contains "[" follwed by an integer number)
					       bool isauxvec = false;
					       pos_aux_vec=buffer_s.find("[");
					       if(pos_aux_vec!=string::npos){
						       bool isint = true;
						       for(size_t ch=pos_aux_vec+1;ch<buffer_s.size()-1;ch++){
							       if( !isdigit(buffer_s[ch]) ){
								       isint = false;
								       break;
							       }
						       }
						       if( isint ){
							       isauxvec = true;
							       aux_name = buffer_s.substr(0,pos_aux_vec);
						       }
					       }
					       if( isauxvec ){
						       bool already_stored = false;
						       unsigned int aux_ind;
						       for(unsigned int av=0;av<Aux_name.size();av++){
							       if( aux_name == Aux_name[av] ){
								       already_stored = true;
								       aux_ind = av;
								       break;
							       }
						       }
						       if( already_stored ) Aux_size[aux_ind] += 1;
						       else{
							       Aux_name.push_back(aux_name);
							       Aux_size.push_back(1);
						       }
					       }else{
						       Aux_name.push_back(buffer_s);
						       Aux_size.push_back(1);
					       }
					}else{
						befAuxNames.push_back(buffer_s);		
					}		
				}
				for( unsigned int bs=0;bs<befAuxNames.size();bs++){
					if( befAuxNames[bs] == "xs" || befAuxNames[bs] == "ys" || befAuxNames[bs] == "zs" ){
						IsReducedCoords = true;
						break;
					}
				}
				for(unsigned int au=0;au<Aux_name.size();au++) this->Aux.push_back(new double[this->nbAtom*Aux_size[au]]);
				line_at = count;
				nbAux -= nbAux_norm;
				final_line = line_at+1+this->nbAtom;
			}
			if( count > line_at  && count < final_line && GoodTimestep ){
				istringstream text(line);
				for(unsigned int bs=0;bs<befAuxNames.size();bs++){
					if( befAuxNames[bs] == "id" ) text >> buffer_uint;
					else if( befAuxNames[bs] == "type" ) text >> buffer_uint_1;
					else if( befAuxNames[bs] == "element" ) text >> buffer_s;
					else if( befAuxNames[bs] == "xs" || befAuxNames[bs] == "x" || befAuxNames[bs] == "xu" ) text >> buffer_1;
					else if( befAuxNames[bs] == "ys" || befAuxNames[bs] == "y" || befAuxNames[bs] == "yu" ) text >> buffer_2;
					else if( befAuxNames[bs] == "zs" || befAuxNames[bs] == "z" || befAuxNames[bs] == "zu" ) text >> buffer_3;
					else if( befAuxNames[bs] == "q" ) text >> buffer_4;
				}
				if( IsReducedCoords ){
					this->AtomList[NbAtRead].pos.x = buffer_1*H1[0]+buffer_2*H2[0]+buffer_3*H3[0];
					this->AtomList[NbAtRead].pos.y = buffer_1*H1[1]+buffer_2*H2[1]+buffer_3*H3[1];
					this->AtomList[NbAtRead].pos.z = buffer_1*H1[2]+buffer_2*H2[2]+buffer_3*H3[2];
				}else{
					this->AtomList[NbAtRead].pos.x = buffer_1;
					this->AtomList[NbAtRead].pos.y = buffer_2;
					this->AtomList[NbAtRead].pos.z = buffer_3;
				}
				if( !IsType_uint ){
					if( IsElem ){
						bool elem_already_stored = false;
						unsigned int current_typeuint = 0;
						for(unsigned int i=0;i<nbAtomType;i++){
							if( buffer_s == AtomType[i] ){
								elem_already_stored = true;
								current_typeuint = i;
								break;
							}
						}
						if( !elem_already_stored ){
							if( nbAtomType >= MaxAtomType-1 ){
								cerr << "Number of atom type exceed maximum number of ato type allowed" << endl;
								exit(EXIT_FAILURE);
							}
							AtomType[nbAtomType] = buffer_s;
							nbAtomType += 1;
							this->AtomList[NbAtRead].type_uint = nbAtomType;
						}else this->AtomList[NbAtRead].type_uint = current_typeuint+1;
					}else{
						this->AtomList[NbAtRead].type_uint = 1;
					}
				}else{
					this->AtomList[NbAtRead].type_uint = buffer_uint_1;
					if( IsElem ){
						if( buffer_uint_1 >= MaxAtomType ){
							cerr << "Number of atom type exceed maximum number of ato type allowed" << endl;
							exit(EXIT_FAILURE);
						}
						this->AtomType[buffer_uint_1-1] = buffer_s;
					}
				}
				if( IsCharge ) this->AtomCharge[AtomList[NbAtRead].type_uint-1] = buffer_4;
				if(nbAux>0){
					for(unsigned int au=0;au<Aux_name.size();au++){
						for(unsigned int au_v=0;au_v<Aux_size[au];au_v++) text >> Aux[au][NbAtRead*Aux_size[au]+au_v];
					}
				}
				++NbAtRead;
			}
			count += 1;
			if( count >= final_line && GoodTimestep ) break;
		}
		file.close();
		if( ReadOk != 3 ) return false;
		// search the number of atom type
		vector<unsigned int> AtTypeDif;
		bool already;
		for(unsigned int i=0;i<nbAtom;i++){
			already = false;
			for(unsigned int k=0;k<AtTypeDif.size();k++){
				if( AtomList[i].type_uint == AtTypeDif[k] ){
					already = true;
					break;
				}
			}
			if( !already ) AtTypeDif.push_back(AtomList[i].type_uint); // TODO here do better thing with type_uint
		}
		nbAtomType = AtTypeDif.size();

                // read MASSES database to get the masses of ions
                //vector<string> element;
                //vector<double> masses;
                //char *database_env = getenv("MASSES_DATABASE");
                //string database;
                //if (database_env) {
                //      database = database_env;
                //} else {
                //      #ifdef MASSES_DATABASE
                //      database = MASSES_DATABASE;
                //      #endif
                //}
                //if( database.empty() ) cout << "Warning database environment for masses is empty" << endl;
                //else{
                //      string dataname = database + "/Masses.txt";
                //      ifstream filedata(dataname, ios::in);
                //      if( filedata ){
                //              while( filedata ){
                //                      filedata >> buffer_s >> buffer_1;
                //                      if( !filedata ) break;
                //                      element.push_back(buffer_s);
                //                      masses.push_back(buffer_1);
                //              }
                //              filedata.close();
                //      }
                //}

                //bool element_find;
                //for(unsigned int j=0;j<this->nbAtomType;j++){
                //      element_find= false;
                //      for(unsigned int i=0;i<element.size();i++){
                //              if( this->AtomType[j] == element[i] ){
                //                      this->AtomMass[j] = masses[i];
                //                      element_find = true;
                //                      break;
                //              }
                //      }
                //      if( !element_find ) cout << "The mass for element " << this->AtomType[j] << " has not been found in the masses databse" << endl;
                //}

	// end read cfg (xsf) file
	}else{
		cout << "The file " << filename << " cannot be openned" << endl;
		exit(EXIT_FAILURE);
	}
	if( ReadOk == 3 ){
		if( (NbAtRead) == nbAtom ){
			return true;
		}else{
			cout << "the number of atom provided does not correspond to the true number of atom" << endl;
			exit(EXIT_FAILURE);
			return false;
		}
	}else{
		return false;
	}
}

// Compute periodic array (for atom_style full, allowing to get the right bond and angle relationship) assuming only that ions outside the box have non-null component of the periodic array
void AtomicSystem::ComputePeriodicArr(){
	if( !IsPeriodicArr ){
		PeriodicArr = new int[nbAtom*3];
		IsPeriodicArr = true;
	}
	// Initialize to zero
	unsigned int loop=nbAtom*3;
	for(unsigned int i=0;i<loop;i++) PeriodicArr[i] = 0;
	
	// compute reduced coordinates
	double x,y,z;
	for(unsigned int i=0;i<this->nbAtom;i++){
		x = this->AtomList[i].pos.x*G1[0]+this->AtomList[i].pos.y*G2[0]+this->AtomList[i].pos.z*G3[0];
		PeriodicArr[i*3] = -floor(x);
		y = this->AtomList[i].pos.x*G1[1]+this->AtomList[i].pos.y*G2[1]+this->AtomList[i].pos.z*G3[1];
		PeriodicArr[i*3+1] = -floor(y);
		z = this->AtomList[i].pos.x*G1[2]+this->AtomList[i].pos.y*G2[2]+this->AtomList[i].pos.z*G3[2];
		PeriodicArr[i*3+2] = -floor(z);
	}
}

// TODO : change format for having more prec
void AtomicSystem::print_lmp(const string& filename){
	ofstream writefile(filename);
	writefile << " # File generated using AtomHic\n";
	writefile << this->File_Heading;
        writefile << "\n\t" << this->nbAtom << "\tatoms\n\t";
	if( IsBond ) ComputePeriodicArr();
	if( IsBond ) writefile << nbBonds << "\tbonds\n\t";
	if( IsAngle ) writefile << nbAngles << "\tangles\n\t";
        writefile << this->nbAtomType << "\tatom types\n\t";
	if( IsBond ) writefile << nbBondType << "\tbond types\n\t";
	if( IsAngle ) writefile << nbAngleType << "\tangle types\n\t";
	writefile << "\n\t0.000000000000\t" << this->H1[0] << "\txlo xhi\n\t0.000000000000\t" << H2[1] << "\tylo yhi\n\t0.000000000000\t" << H3[2] << "\tzlo zhi\n";
	if( this->IsTilted ) writefile << "\t" << H2[0] << "\t" << H3[0] << "\t" << H3[1] << "\txy xz yz\n";
	writefile << "\nMasses\n\n";
	for(unsigned int i=0;i<this->nbAtomType;i++) writefile << "\t" << i+1 << "\t" << this->AtomMass[i] << "\t# " << this->AtomType[i] << "\n";
	if( IsMolId ){
		writefile << "\nAtoms # full\n\n";
		for(unsigned int i=0;i<this->nbAtom;i++) writefile << i+1 << "\t" << this->MolId[i] << "\t" << this->AtomList[i].type_uint << "\t" << this->AtomCharge[this->AtomList[i].type_uint-1] << "\t" << this->AtomList[i].pos.x << "\t" << this->AtomList[i].pos.y << "\t" << this->AtomList[i].pos.z << "\t" << this->PeriodicArr[i*3] << "\t" << this->PeriodicArr[i*3+1] << "\t" << this->PeriodicArr[i*3+2] << "\n"; 
	}else if( IsCharge ){
		writefile << "\nAtoms # charge\n\n";
		for(unsigned int i=0;i<this->nbAtom;i++) writefile << i+1 << "\t" << this->AtomList[i].type_uint << "\t" << this->AtomCharge[this->AtomList[i].type_uint-1] << "\t" << this->AtomList[i].pos.x << "\t" << this->AtomList[i].pos.y << "\t" << this->AtomList[i].pos.z << "\n"; 
	}else{
		writefile << "\nAtoms # atomic\n\n";
		for(unsigned int i=0;i<this->nbAtom;i++) writefile << i+1 << "\t" << this->AtomList[i].type_uint << "\t" << this->AtomList[i].pos.x << "\t" << this->AtomList[i].pos.y << "\t" << this->AtomList[i].pos.z << "\n"; 
	}
	if( IsVel ){
		writefile << "\nVelocities\n\n";
		for(unsigned int i=0;i<nbAtom;i++) writefile << i+1 << "\t" << Aux[0][(i*3)] << "\t" << Aux[0][(i*3)+1] << "\t" << Aux[0][(i*3)+2] << "\n";
	}
	if( IsBond ){
		writefile << "\nBonds\n\n";
		for(unsigned int i=0;i<nbBonds;i++) writefile << i+1 << "\t" << BondType[i] << "\t" << Bonds[i*2] << "\t" << Bonds[i*2+1] << "\n";
	}
	if( IsAngle ){
		writefile << "\nAngles\n\n";
		for(unsigned int i=0;i<nbAngles;i++) writefile << i+1 << "\t" << AngleType[i] << "\t" << Angles[i*3] << "\t" << Angles[i*3+1] << "\t" << Angles[i*3+2] << "\n";
	}

	writefile.close();
	cout << "File " << filename << " successfully writted ! (" << SystemCharacteristics() << ")" << endl;
}

void AtomicSystem::print_cfg(const string& filename){
	ofstream writefile(filename);
	writefile << "ITEM: TIMESTEP\n" << (int) this->timestep << "\nITEM: NUMBER OF ATOMS\n" << this->nbAtom << "\nITEM: ";
	if( IsTilted ){
		// compute the cell vectors
		double arr[4] = {0.,this->H2[0],this->H3[0],this->H2[0]+this->H3[0]};
		double arr_2[2] = {0.,this->H3[1]};
		// TODO issues here
	       	writefile << "BOX BOUNDS xy xz yz pp pp pp\n" << this->MT->min(arr,4) << "\t" << this->H1[0]+this->MT->max(arr,4) << "\t" << H2[0] << "\n" << this->MT->min(arr_2,2) << "\t" << this->H2[1]+this->MT->max(arr_2,2) << "\t" << H3[0] << "\n0\t" << H3[2] << "\t" << H3[1] << "\n";
	}
	else writefile << "BOX BOUNDS pp pp pp\n" << "0\t" << H1[0] << "\n0\t" << H2[1] << "\n0\t" << H3[2] << "\n";
	if( IsMolId ){
		writefile << "ITEM: ATOMS id mol type element xu yu zu q\n";
		for(unsigned int i=0;i<this->nbAtom;i++) writefile << i+1 << " " << this->MolId[i] << " " << this->AtomList[i].type_uint << " " << this->AtomType[this->AtomList[i].type_uint-1] << " " << this->AtomList[i].pos.x << " " << this->AtomList[i].pos.y << " " << this->AtomList[i].pos.z << " " << this->AtomCharge[this->AtomList[i].type_uint-1] << "\n";
	}else if( IsElem ){
		writefile << "ITEM: ATOMS id type element xu yu zu";
		if( IsCharge ) writefile << " q";
		writefile << "\n";
		for(unsigned int i=0;i<this->nbAtom;i++){
			writefile << i+1 << " " << this->AtomList[i].type_uint << " " << this->AtomType[this->AtomList[i].type_uint-1] << " " << this->AtomList[i].pos.x << " " << this->AtomList[i].pos.y << " " << this->AtomList[i].pos.z;
		       if( IsCharge ) writefile	<< " " << this->AtomCharge[this->AtomList[i].type_uint-1];
		       writefile << "\n";
		}
	}else{
		writefile << "ITEM: ATOMS id type xu yu zu";
		if( IsCharge ) writefile << " q";
		writefile << "\n";
		for(unsigned int i=0;i<this->nbAtom;i++){
			writefile << i+1 << " " << this->AtomList[i].type_uint << " " << this->AtomList[i].pos.x << " " << this->AtomList[i].pos.y << " " << this->AtomList[i].pos.z;
		       if( IsCharge ) writefile	<< " " << this->AtomCharge[this->AtomList[i].type_uint-1];
		       writefile << "\n";
		}
	}
	writefile.close();
	cout << "File " << filename << " successfully writted ! (" << SystemCharacteristics(true) << ")" << endl;
}

void AtomicSystem::printSystem_aux(const string& filename, const string& AuxName, bool append){
	vector<string> AuxNames_v;
	string buffer;
	string AuxName_a = AuxProp2Print+AuxName;
	istringstream names(AuxName_a);
	while( getline(names, buffer, ' ') ){
		AuxNames_v.push_back(buffer);
	}
	vector<unsigned int> AuxId;
	for(unsigned int j=0;j<AuxNames_v.size();j++){
		bool aux_find = false;
		for(unsigned int i=0;i<this->Aux_name.size();i++){
			if(AuxNames_v[j] == Aux_name[i]){
				AuxId.push_back(i);
				aux_find = true;
				break;
			}
		}
		if( aux_find == false ) cout << "The auxiliary property \"" << AuxNames_v[j] << "\" has not been found" << endl;
	}
	
	ofstream writefile;
	if( append ) writefile.open(filename, ios_base::app);
	else writefile.open(filename);

	writefile << "ITEM: TIMESTEP\n" << (int) this->timestep << "\nITEM: NUMBER OF ATOMS\n" << this->nbAtom << "\nITEM: ";
	if( IsTilted ){
		// compute the cell vectors
		double arr[4] = {0.,this->H2[0],this->H3[0],this->H2[0]+this->H3[0]};
		double arr_2[2] = {0.,this->H3[1]};
	       	writefile << "BOX BOUNDS xy xz yz pp pp pp\n" << this->MT->min(arr,4) << "\t" << this->H1[0]+this->MT->max(arr,4) << "\t" << H2[0] << "\n" << this->MT->min(arr_2,2) << "\t" << this->H2[1]+this->MT->max(arr_2,2) << "\t" << H3[0] << "\n0\t" << H3[2] << "\t" << H3[1] << "\n";
	}
	else writefile << "BOX BOUNDS pp pp pp\n" << "0\t" << H1[0] << "\n0\t" << H2[1] << "\n0\t" << H3[2] << "\n";
	if( IsMolId ){
		writefile << "ITEM: ATOMS id mol type element xu yu zu";
		if( this->IsCharge ) writefile << " q";
        	for(unsigned int i=0;i<AuxId.size();i++){
			if( Aux_size[AuxId[i]] == 1 ) writefile << " " << this->Aux_name[AuxId[i]];
			else for(unsigned int j=0;j<Aux_size[AuxId[i]];j++) writefile << " " << this->Aux_name[AuxId[i]] << "[" << j+1 << "]";
		}
		writefile << "\n";
		if( this->IsCharge ){
			for(unsigned int i=0;i<this->nbAtom;i++){
				writefile << i+1 << " " << this->MolId[i] << " " << this->AtomList[i].type_uint << " " << this->AtomType[this->AtomList[i].type_uint-1] << " " << this->AtomList[i].pos.x << " " << this->AtomList[i].pos.y << " " << this->AtomList[i].pos.z << " " << this->AtomCharge[this->AtomList[i].type_uint-1];
        			for(unsigned int j=0;j<AuxId.size();j++){
					if( Aux_size[AuxId[j]] == 1 ) writefile << " " << this->Aux[AuxId[j]][i];
					else for(unsigned int k=0;k<Aux_size[AuxId[j]];k++) writefile << " " << this->Aux[AuxId[j]][i*Aux_size[AuxId[j]]+k];
				}
				writefile << "\n";
			}
		}else{
			for(unsigned int i=0;i<this->nbAtom;i++){
				writefile << i+1 << " " << this->AtomList[i].type_uint << " " << this->AtomType[this->AtomList[i].type_uint-1] << " " << this->AtomList[i].pos.x << " " << this->AtomList[i].pos.y << " " << this->AtomList[i].pos.z;
        			for(unsigned int j=0;j<AuxId.size();j++){
					if( Aux_size[AuxId[j]] == 1 ) writefile << " " << this->Aux[AuxId[j]][i];
					else for(unsigned int k=0;k<Aux_size[AuxId[j]];k++) writefile << " " << this->Aux[AuxId[j]][i*Aux_size[AuxId[j]]+k];
				}
				writefile << "\n";
			}
		}
	}else if( IsElem ){
		writefile << "ITEM: ATOMS id type element xu yu zu";
		if( this->IsCharge ) writefile << " q";
        	for(unsigned int i=0;i<AuxId.size();i++){
			if( Aux_size[AuxId[i]] == 1 ) writefile << " " << this->Aux_name[AuxId[i]];
			else for(unsigned int j=0;j<Aux_size[AuxId[i]];j++) writefile << " " << this->Aux_name[AuxId[i]] << "[" << j+1 << "]";
		}
		writefile << "\n";
		if( this->IsCharge ){
			for(unsigned int i=0;i<this->nbAtom;i++){
				writefile << i+1 << " " << this->AtomList[i].type_uint << " " << this->AtomType[this->AtomList[i].type_uint-1] << " " << this->AtomList[i].pos.x << " " << this->AtomList[i].pos.y << " " << this->AtomList[i].pos.z << " " << this->AtomCharge[this->AtomList[i].type_uint-1];
        			for(unsigned int j=0;j<AuxId.size();j++){
					if( Aux_size[AuxId[j]] == 1 ) writefile << " " << this->Aux[AuxId[j]][i];
					else for(unsigned int k=0;k<Aux_size[AuxId[j]];k++) writefile << " " << this->Aux[AuxId[j]][i*Aux_size[AuxId[j]]+k];
				}
				writefile << "\n";
			}
		}else{
			for(unsigned int i=0;i<this->nbAtom;i++){
				writefile << i+1 << " " << this->AtomList[i].type_uint << " " << this->AtomType[this->AtomList[i].type_uint-1] << " " << this->AtomList[i].pos.x << " " << this->AtomList[i].pos.y << " " << this->AtomList[i].pos.z;
        			for(unsigned int j=0;j<AuxId.size();j++){
					if( Aux_size[AuxId[j]] == 1 ) writefile << " " << this->Aux[AuxId[j]][i];
					else for(unsigned int k=0;k<Aux_size[AuxId[j]];k++) writefile << " " << this->Aux[AuxId[j]][i*Aux_size[AuxId[j]]+k];
				}
				writefile << "\n";
			}
		}
	}else{
		writefile << "ITEM: ATOMS id type xu yu zu";
		if( this->IsCharge ) writefile << " q";
        	for(unsigned int i=0;i<AuxId.size();i++){
			if( Aux_size[AuxId[i]] == 1 ) writefile << " " << this->Aux_name[AuxId[i]];
			else for(unsigned int j=0;j<Aux_size[AuxId[i]];j++) writefile << " " << this->Aux_name[AuxId[i]] << "[" << j+1 << "]";
		}
		writefile << "\n";
		if( this->IsCharge ){
			for(unsigned int i=0;i<this->nbAtom;i++){
				writefile << i+1 << " " << this->AtomList[i].type_uint << " " << this->AtomList[i].pos.x << " " << this->AtomList[i].pos.y << " " << this->AtomList[i].pos.z << " " << this->AtomCharge[this->AtomList[i].type_uint-1];
        			for(unsigned int j=0;j<AuxId.size();j++){
					if( Aux_size[AuxId[j]] == 1 ) writefile << " " << this->Aux[AuxId[j]][i];
					else for(unsigned int k=0;k<Aux_size[AuxId[j]];k++) writefile << " " << this->Aux[AuxId[j]][i*Aux_size[AuxId[j]]+k];
				}
				writefile << "\n";
			}
		}else{
			for(unsigned int i=0;i<this->nbAtom;i++){
				writefile << i+1 << " " << this->AtomList[i].type_uint << " " << this->AtomList[i].pos.x << " " << this->AtomList[i].pos.y << " " << this->AtomList[i].pos.z;
        			for(unsigned int j=0;j<AuxId.size();j++){
					if( Aux_size[AuxId[j]] == 1 ) writefile << " " << this->Aux[AuxId[j]][i];
					else for(unsigned int k=0;k<Aux_size[AuxId[j]];k++) writefile << " " << this->Aux[AuxId[j]][i*Aux_size[AuxId[j]]+k];
				}
				writefile << "\n";
			}
		}
	}
	writefile.close();
	if( append ) cout << "timestep " << timestep;
	else cout << "File " << filename;
        cout << " successfully writted ! (" << SystemCharacteristics(true);
	cout << ", auxiliary properties :";
        for(unsigned int i=0;i<AuxId.size();i++) cout << " " << this->Aux_name[AuxId[i]] << " - " << Aux_size[AuxId[i]] << " dim -";
	cout << ")" << endl;
}

void AtomicSystem::read_params_atsys(){
	//string fp;
	//#ifdef FIXEDPARAMETERS
	//fp = FIXEDPARAMETERS;
	//#endif
	//string backslash="/";
	//string filename=fp+backslash+FixedParam_Filename;
	//ifstream file(filename, ios::in);
	ifstream file(FixedParam_Filename, ios::in);
	size_t pos_auxprop;
	string buffer_s, line;
	if(file){
		while(file){
			getline(file,line);
			pos_auxprop=line.find("AUX_PROPERTIES_TO_PRINT ");
			if(pos_auxprop!=string::npos){
				istringstream text(line);
				text >> buffer_s;
				while( text >> buffer_s ){
					if( buffer_s == "//" ) break;
					AuxProp2Print += buffer_s+" ";
				}
			}
		}
	}
	//else{
	//	cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
	//	exit(EXIT_FAILURE);
	//}
}

vector<unsigned int> AtomicSystem::selectAtomInBox(const double x_lo,const double x_hi,const double y_lo,const double y_hi,const double z_lo,const double z_hi, bool OnWrap){
	if( !IsWrappedPos && OnWrap ) computeWrap();
	vector<unsigned int> AtList;
	if( OnWrap ){
		for(unsigned int i=0;i<this->nbAtom;i++)
			if( this->WrappedPos[i].x > x_lo && this->WrappedPos[i].x < x_hi && this->WrappedPos[i].y > y_lo && this->WrappedPos[i].y < y_hi && this->WrappedPos[i].z > z_lo && this->WrappedPos[i].z < z_hi ){ AtList.push_back(i); }
	}else{
		for(unsigned int i=0;i<this->nbAtom;i++)
			if( this->AtomList[i].pos.x > x_lo && this->AtomList[i].pos.x < x_hi && this->AtomList[i].pos.y > y_lo && this->AtomList[i].pos.y < y_hi && this->AtomList[i].pos.z > z_lo && this->AtomList[i].pos.z < z_hi ){ AtList.push_back(i); }
	}
	
	return AtList;
}

unsigned int AtomicSystem::getAuxIdAndSize(std::string auxname, unsigned int &size){
	unsigned int ind;
	bool found = false;
	for(unsigned int i=0;i<Aux.size();i++){
		if( Aux_name[i] == auxname ){
		 	ind = i;
			size = Aux_size[i];
			found = true;
			break;
		}
	}
	if( !found ){
		cout << "Warning, the auxiliary property : \"" << auxname << "\" has not been found !" << endl;
		ind = 0;
		size = 0;
	}
	return ind; 
}

void AtomicSystem::ApplyShift(const double& shift_x, const double& shift_y, const double& shift_z){
	for(unsigned int i = 0; i < this->nbAtom; i++){
		this->AtomList[i].pos.x += shift_x;
		this->AtomList[i].pos.y += shift_y;
		this->AtomList[i].pos.z += shift_z;
	}
}

void AtomicSystem::duplicate(const unsigned int& nx, const unsigned int& ny, const unsigned int& nz){
	if( nx == 0 || ny == 0 || nz == 0 ){
		cerr << "Cannot duplicate if one value is equal to zero, aborting" << endl;
		exit(EXIT_FAILURE);
	}
	// TODO manage when IsAtomListMine is false and replace in bicrystal constructor with this method
	if( !IsAtomListMine ){
		cerr << "The AtomList array does not belong to the AtomicSystem object, cannot duplicate (to be implemented), aborting" << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Duplicating the system (" << nx << "," << ny << "," << nz << ")" << endl;
	
	unsigned int new_nbAtom = this->nbAtom*nx*ny*nz;
	Atom *AtomList_temp = new Atom[nbAtom];
	for(unsigned int i=0;i<nbAtom;i++) AtomList_temp[i] = AtomList[i];
	delete[] AtomList;
	AtomList = new Atom[new_nbAtom];
	unsigned int new_nbBonds, new_nbAngles, nbMol;
	unsigned int *MolId_temp, *Bonds_temp, *BondType_temp, *Angles_temp, *AngleType_temp;
	vector<double*> Aux_temp;
	if( IsMolId ){
		MolId_temp = new unsigned int[nbAtom];
		for(unsigned int i=0;i<nbAtom;i++) MolId_temp[i] = MolId[i];
		nbMol = MT->max_p(MolId,nbAtom);
		delete[] MolId;
		MolId = new unsigned int[new_nbAtom];
	}
	if( IsBond ){
		new_nbBonds = nbBonds*nx*ny*nz;
		Bonds_temp = new unsigned int[nbBonds*2];
		BondType_temp = new unsigned int[nbBonds];
		for(unsigned int i=0;i<nbBonds;i++){
			for(unsigned int d=0;d<2;d++) Bonds_temp[i*2+d] = Bonds[i*2+d];
			BondType_temp[i] = BondType[i];
		}
		delete[] Bonds;
		delete[] BondType;
		Bonds = new unsigned int[new_nbBonds*2];
		BondType = new unsigned int[new_nbBonds];
	}
	if( IsAngle ){
		new_nbAngles = nbAngles*nx*ny*nz;
		Angles_temp = new unsigned int[nbAngles*3];
		AngleType_temp = new unsigned int[nbAngles];
		for(unsigned int i=0;i<nbAngles;i++){
			for(unsigned int d=0;d<3;d++) Angles_temp[i*3+d] = Angles[i*3+d];
			AngleType_temp[i] = AngleType[i];
		}
		delete[] Angles;
		delete[] AngleType;
		Angles = new unsigned int[new_nbAngles*3];
		AngleType = new unsigned int[new_nbAngles];
	}
	if( IsSetAux ){
		for(unsigned int i=0;i<Aux.size();i++){
			Aux_temp.push_back(new double[(nbAtom)*Aux_size[i]]);
			for(unsigned int n=0;n<nbAtom;n++){
				for(unsigned int d=0;d<Aux_size[i];d++) Aux_temp[i][n*Aux_size[i]+d] = Aux[i][n*Aux_size[i]+d];
			}
		}
		for(unsigned int i=0;i<Aux.size();i++){
			delete[] Aux[i];
			Aux[i] = new double[(new_nbAtom)*Aux_size[i]];
		}
	}
	if( IsPeriodicArr ){
		delete[] PeriodicArr;
		PeriodicArr = new int[new_nbAtom*3];
	}
	
	unsigned int ind, ind_a;
	for(unsigned int i=0;i<nx;i++){
		for(unsigned int j=0;j<ny;j++){
			for(unsigned int k=0;k<nz;k++){
				ind = i*ny*nz*nbAtom + j*nz*nbAtom + k*nbAtom;
				for(unsigned int n=0;n<nbAtom;n++){
					ind_a = ind + n;
					AtomList[ind_a] = AtomList_temp[n];
					AtomList[ind_a].pos.x += i*H1[0] + j*H2[0] + k*H3[0];
					AtomList[ind_a].pos.y += i*H1[1] + j*H2[1] + k*H3[1];
					AtomList[ind_a].pos.z += i*H1[2] + j*H2[2] + k*H3[2];
					if( IsMolId ) MolId[ind_a] = MolId_temp[n] + i*ny*nz*nbMol + j*nz*nbMol + k*nbMol;
					if( IsSetAux ){
						for(unsigned int a=0;a<Aux_temp.size();a++){
							for(unsigned int d=0;d<Aux_size[a];d++) Aux[a][ind_a*Aux_size[a]+d] = Aux_temp[a][n*Aux_size[a]+d];
						}
					}
				}
				if( IsBond ){
					ind = i*ny*nz*nbBonds+j*nz*nbBonds+k*nbBonds;
					ind_a = i*ny*nz*nbAtom + j*nz*nbAtom + k*nbAtom;
					for(unsigned int n=0;n<nbBonds;n++){
						unsigned int c_ind_n = ind+n;
						for(unsigned int d=0;d<2;d++) Bonds[c_ind_n*2+d] = Bonds_temp[n*2+d] + ind_a;
						BondType[c_ind_n] = BondType_temp[n];
					}
				}
				if( IsAngle ){
					ind = i*ny*nz*nbAngles+j*nz*nbAngles+k*nbAngles;
					ind_a = i*ny*nz*nbAtom + j*nz*nbAtom + k*nbAtom;
					for(unsigned int n=0;n<nbAngles;n++){
						unsigned int c_ind_n = ind+n;
						for(unsigned int d=0;d<3;d++) Angles[c_ind_n*3+d] = Angles_temp[n*3+d] + ind_a;
						AngleType[c_ind_n] = AngleType_temp[n];
					}
				}
			}
		}
	}
	nbAtom = new_nbAtom;
	delete[] AtomList_temp;
	if( IsMolId ) delete[] MolId_temp;
	if( IsBond ){
		nbBonds = new_nbBonds;
		delete[] Bonds_temp;
		delete[] BondType_temp;
	}
	if( IsAngle ){
		nbAngles = new_nbAngles;
		delete[] Angles_temp;
		delete[] AngleType_temp;
	}
	if( IsSetAux ) for(unsigned int i=0;i<Aux_temp.size();i++) delete[] Aux_temp[i];
	
	for(unsigned int i=0;i<3;i++){
		this->H1[i] *= nx;
		this->H2[i] *= ny;
		this->H3[i] *= nz;
	}
	computeInverseCellVec();

	if( this->IsWrappedPos ){
		delete[] WrappedPos;
		this->IsWrappedPos = false;
	}
}

// The two next functions have not been succesfully implemented one day may be pushed further
void AtomicSystem::MakeSurfaceNeutral_3dBased(string ext){
	if( IsCharge == false ) return;
	bool Possible = false; // at least we should have two atom types having charge with opposite sign
	bool pos(false),neg(false);
	double zero_num = 1e-8;
	if( !_MyCrystal->getIsDoNotSep() ){
		for(unsigned int n=0;n<nbAtomType;n++){
			if( AtomCharge[n] > zero_num ) pos = true;
			else if( AtomCharge[n] < -zero_num ) neg = true;
			if( pos && neg ){
				Possible = true;
				break;
			}
		}
	}else{
		ComputeNotSepList();
		vector<vector<unsigned int>> DNS = _MyCrystal->getDoNotSep();
		unsigned int *Stoich = _MyCrystal->getStoich();
		for(unsigned int n=0;n<nbAtomType;n++){
			double block_charge = 0.; 
			for(unsigned int i=0;i<DNS.size();i++){
				if( DNS[i][0]-1 == n ){
					if( i == 0 ) block_charge = AtomCharge[DNS[i][0]-1];
					block_charge += AtomCharge[DNS[i][2]-1]*DNS[i][1];
					if( block_charge > zero_num ) pos = true;
					else if( block_charge < -zero_num ) neg = true;
					if( Stoich[DNS[i][2]-1] > Stoich[DNS[i][0]-1]*DNS[i][1] ){
						if( AtomCharge[DNS[i][2]-1] > zero_num ) pos = true;
						else if( AtomCharge[DNS[i][2]-1] < -zero_num ) neg = true;
					}
					if( pos && neg ){
						Possible = true;
						break;
					}
				}
			}
		}
		for(unsigned int n=0;n<nbAtomType;n++){
			bool treat = true;
			for(unsigned int i=0;i<DNS.size();i++){
				if( DNS[i][0]-1 == n ){
					treat = false;
					break;
				}
				if( DNS[i][2]-1 == n ){
					treat = false;
					break;
				}
			}
			if( treat ){
				if( AtomCharge[n] > zero_num ) pos = true;
				else if( AtomCharge[n] < -zero_num ) neg = true;
				if( pos && neg ){
					Possible = true;
					break;
				}
			}
		}
	}	

	if( !Possible ){
		cout << "It is not possible to construct charge neutral surfaces as there are not two types of ions (or blocks of ions due to the DoNotSep crystal tag) having charge with opposite signs" << endl;
		return;
	}

	cout << "Trying to make the surface charge neutral" << endl;
	// Compute charge density at each ion position
	double fac_sigma_gauss = 1.5;
	double *ChargeDens = new double[nbAtom];
	double sigma_gauss = MT->max_p(_MyCrystal->getALength(),3);
	sigma_gauss *= fac_sigma_gauss;
	int bcx = 2;
	int bcy = 2;
	int bcz = 2;
	#pragma omp parallel for
	for(unsigned int i=0;i<nbAtom;i++){
		ChargeDens[i] = 0.;
		for(unsigned int j=0;j<nbAtom;j++){ // think about restrict to neighbour in cutoff
			for(int bc1=-bcx;bc1<=bcx;bc1++){
				double coord1 = this->AtomList[j].pos.x+(this->H1[0]*bc1);
				for(int bc2=-bcy;bc2<=bcy;bc2++){
					double coord2 = this->AtomList[j].pos.y+(this->H2[1]*bc2);
					for(int bc3=-bcz;bc3<=bcz;bc3++){
						double coord3 = this->AtomList[j].pos.z+(this->H3[2]*bc3);
						ChargeDens[i] += MT->gaussian(AtomList[i].pos.x, AtomList[i].pos.y, AtomList[i].pos.z, coord1, coord2, coord3, sigma_gauss)*AtomCharge[AtomList[j].type_uint-1];
					}
				}
			}
		}
	}
	setAux(ChargeDens,"ChargeDens");

	// remove ions in areas of high charge density with the same sign of the ion and exclude its neighbours for this pass
	double rcut = MT->min_p(_MyCrystal->getALength(),3);
	//searchNeighbours(rcut);
	sigma_gauss /= fac_sigma_gauss;
	searchNeighbours(sigma_gauss);
	sigma_gauss *= fac_sigma_gauss;
	double max_charge = 7.e-5;
	unsigned int *removed = new unsigned int[nbAtom];
	double *classif = new double[nbAtom];
	double *current_charge = new double[nbAtom];
	unsigned int *dontreat = new unsigned int[nbAtom];
	double total_charge = 0.;
	vector<unsigned int> ions2remove;
	for(unsigned int i=0;i<nbAtom;i++){
		dontreat[i] = 0;
		removed[i] = 0;
		classif[i] = 1. / (double) Neighbours[i*(nbMaxN+1)];
		//classif[i] = fabs( ChargeDens[i] / ((double) Neighbours[i*(nbMaxN+1)] * (double) Neighbours[i*(nbMaxN+1)] * (double) Neighbours[i*(nbMaxN+1)]) );
		//classif[i] = fabs( ChargeDens[i] / ((double) Neighbours[i*(nbMaxN+1)]) );
		current_charge[i] = AtomCharge[AtomList[i].type_uint-1];
		total_charge += current_charge[i];
	}

	double mean_q, stdev_q;
	unsigned int current_nbAtoms;
	current_nbAtoms = nbAtom - ions2remove.size();
	mean_q = 0.;
	stdev_q = 0.;
	for(unsigned int i=0;i<nbAtom;i++) mean_q += ChargeDens[i];
	mean_q /= current_nbAtoms;
	for(unsigned int n=0;n<ions2remove.size();n++) ChargeDens[ions2remove[n]] = mean_q;
	for(unsigned int i=0;i<nbAtom;i++)
		//stdev_q += (mean_q-ChargeDens[i])*(mean_q-ChargeDens[i]);
		stdev_q += fabs(ChargeDens[i]);
	stdev_q /= current_nbAtoms;
	for(unsigned int n=0;n<ions2remove.size();n++) ChargeDens[ions2remove[n]] = 0.;
	bool stop = false;
	unsigned int iter = 0;
	unsigned int max_iter = nbAtom;
	ofstream file("distribQ_0.dat");
	for(unsigned int i=0;i<nbAtom;i++) file << ChargeDens[i] << endl;
	file.close();
	file.open("metric");
	setAux(classif,"classbef");
	setAux(classif,"class");
	setAux(removed,"rm");

	vector<unsigned int> select_ind; 
	vector<double> select_stdev;
	select_ind.push_back(ions2remove.size());
	select_stdev.push_back(stdev_q); // TODO use sum(abs(chargedens))/nbAt as loss fct
	
			printSystem_aux("Charges"+ext+"_0.cfg","ChargeDens class rm");
	//while( !stop && iter < max_iter ){ // the loop should be stoped when all ions having charge density higher than max_charge have been removed or assigned to dontreat = 1
	while( iter < max_iter ){ // the loop should be stoped when all ions having charge density higher than max_charge have been removed or assigned to dontreat = 1
		cout << "ITER " << iter << " over " << max_iter << " (" << ions2remove.size() << " removed), stdevQ =" << stdev_q << endl;
		iter++;
		stop = false;
		unsigned int current_ind = MT->max_p_ind(classif,nbAtom);
		vector<unsigned int> current_ions2remove;
		if( IsNotSepTag ){
			if( NotSepTag[current_ind][0] == 0 ) current_ions2remove.push_back(current_ind);
			else if( NotSepTag[current_ind][0] < 0 ){
				unsigned int neigh_ind = (unsigned int) 1-NotSepTag[current_ind][0];
				current_ions2remove.push_back(neigh_ind);	
				for(unsigned int n=0;n<NotSepTag[neigh_ind][0];n++) current_ions2remove.push_back(NotSepTag[neigh_ind][n+1]);
			}else{
				current_ions2remove.push_back(current_ind);
				for(unsigned int n=0;n<NotSepTag[current_ind][0];n++) current_ions2remove.push_back(NotSepTag[current_ind][n+1]);
			}
		}else current_ions2remove.push_back(current_ind);
		double block_charge = 0.;
		for(unsigned int n=0;n<current_ions2remove.size();n++) block_charge += AtomCharge[AtomList[current_ions2remove[n]].type_uint-1];

		if( fabs(ChargeDens[current_ind]) >  max_charge && signbit(ChargeDens[current_ind]) == signbit(block_charge) ){
			for(unsigned int n=0;n<current_ions2remove.size();n++){
				// remove contribution of this ion to neighbour list
				//for(unsigned int j=0;j<Neighbours[current_ions2remove[n]*(nbMaxN+1)];j++){
				//	unsigned int current_nid = Neighbours[current_ions2remove[n]*(nbMaxN+1)+1+j];
				//	unsigned int index2rm = 0;
				//	bool found = false;
				//	for(unsigned int j2=0;j2<Neighbours[current_nid*(nbMaxN+1)];j2++){
				//		if( Neighbours[current_nid*(nbMaxN+1)+1+j2] == current_ions2remove[n] ){
				//			index2rm = j2;
				//			found = true;
				//			break;
				//		}
				//	}
				//	if(found){
				//		for(unsigned int j2=index2rm+1;j2<Neighbours[current_nid*(nbMaxN+1)];j2++) Neighbours[current_nid*(nbMaxN+1)+j2] = Neighbours[current_nid*(nbMaxN+1)+1+j2];
				//		Neighbours[current_nid*(nbMaxN+1)]--;
				//	}else cout << "NOTFOUND!!" << endl;
				//}
				// remove contribution of this ion to charge density
				for(unsigned int i=0;i<nbAtom;i++){
					bool treat = true;
					for(unsigned int j=0;j<ions2remove.size();j++){
						if( i == ions2remove[j] ){
							treat = false;
							break;
						}
					}
					if( !treat ) continue;
					for(int bc1=-bcx;bc1<=bcx;bc1++){
						double coord1 = this->AtomList[current_ind].pos.x+(this->H1[0]*bc1);
						for(int bc2=-bcx;bc2<=bcx;bc2++){
							double coord2 = this->AtomList[current_ind].pos.y+(this->H2[1]*bc2);
							for(int bc3=-bcx;bc3<=bcx;bc3++){
								double coord3 = this->AtomList[current_ind].pos.z+(this->H3[2]*bc3);
								ChargeDens[i] -= MT->gaussian(AtomList[i].pos.x, AtomList[i].pos.y, AtomList[i].pos.z, coord1, coord2, coord3, sigma_gauss)*AtomCharge[AtomList[current_ions2remove[n]].type_uint-1];
							}
						}
					}
					// update classif
					classif[i] = 1. / (double) Neighbours[i*(nbMaxN+1)];
					//classif[i] = fabs( ChargeDens[i] / ((double) Neighbours[i*(nbMaxN+1)] * (double) Neighbours[i*(nbMaxN+1)] * (double) Neighbours[i*(nbMaxN+1)]) );
					//classif[i] = fabs( ChargeDens[i] / ((double) Neighbours[i*(nbMaxN+1)] ) );
				}
				classif[current_ions2remove[n]] = 0.;
				removed[current_ions2remove[n]] += 1;
				total_charge -= current_charge[current_ions2remove[n]]; // trouver une loss global pour juger qualite du truc courant et arreter la fct quand le loss est pas mal et quand total_charge = 0
				current_charge[current_ions2remove[n]] = 0.;
				ChargeDens[current_ions2remove[n]] = 0.;
				ions2remove.push_back(current_ions2remove[n]);
			}
			if( fabs(total_charge) < zero_num ){
				// Compute standard dev of chagre density (used as metric)
				current_nbAtoms = nbAtom - ions2remove.size();
				mean_q = 0.;
				stdev_q = 0.;
				for(unsigned int i=0;i<nbAtom;i++) mean_q += ChargeDens[i];
				mean_q /= current_nbAtoms;
				for(unsigned int n=0;n<ions2remove.size();n++) ChargeDens[ions2remove[n]] = mean_q;
				for(unsigned int i=0;i<nbAtom;i++)
					//stdev_q += (mean_q-ChargeDens[i])*(mean_q-ChargeDens[i]);
					stdev_q += fabs(ChargeDens[i]);
				stdev_q /= current_nbAtoms;
				for(unsigned int n=0;n<ions2remove.size();n++) ChargeDens[ions2remove[n]] = 0.;
				select_ind.push_back(ions2remove.size());
				select_stdev.push_back(stdev_q);
			}

			string fname = "Charges"+ext+"_"+to_string(ions2remove.size())+".cfg";
			modifyAux_vec(ChargeDens,"ChargeDens");
			modifyAux_vec(classif,"class");
			modifyAux_vec(removed,"rm");
			printSystem_aux(fname,"ChargeDens class rm");
		}else
			for(unsigned int n=0;n<current_ions2remove.size();n++) classif[current_ions2remove[n]] = 0.; // doublement pas bon ici => essayer de les remettre en jeu
	}
	if( select_ind.size() == 0 ){
		cout << "Neutral system has not been found, aborting making neutral surfaces" << endl;
		return;
	}else{
		unsigned int opt_ind = MT->min(select_stdev);
		cout << "Number of possible neutral systems : " << select_ind.size() << endl;
		cout << "OPTIND corresponds to Nrm = " << select_ind[opt_ind] << " over " << ions2remove.size() << endl;
		unsigned int init_size = ions2remove.size();
		for(unsigned int i=select_ind[opt_ind];i<init_size;i++) ions2remove.pop_back(); // to test
	}
	RemoveAtoms(ions2remove);

	file.close();
	file.open("distribQ_1.dat");
	for(unsigned int i=0;i<nbAtom;i++) file << ChargeDens[i] << endl;
	file.close();
	cout << "nb ions rem = " << ions2remove.size() << ", total charge = " << total_charge << endl;
	setAux(classif,"class");
	setAux(removed,"rm");
	setAux(ChargeDens,"ChargeDensAfter");
	printSystem_aux("TESTCharge.cfg","ChargeDensBefore classbef class ChargeDensAfter rm");

	delete[] ChargeDens;
	delete[] removed;
	delete[] classif;
	delete[] current_charge;
	delete[] dontreat;
}
//
void AtomicSystem::MakeSurfaceNeutral_3dBased_bis(string ext){
	if( IsCharge == false ) return;
	bool Possible = false; // at least we should have two atom types having charge with opposite sign
	bool pos(false),neg(false);
	double zero_num = 1e-8;
	if( !_MyCrystal->getIsDoNotSep() ){
		for(unsigned int n=0;n<nbAtomType;n++){
			if( AtomCharge[n] > zero_num ) pos = true;
			else if( AtomCharge[n] < -zero_num ) neg = true;
			if( pos && neg ){
				Possible = true;
				break;
			}
		}
	}else{
		ComputeNotSepList();
		vector<vector<unsigned int>> DNS = _MyCrystal->getDoNotSep();
		unsigned int *Stoich = _MyCrystal->getStoich();
		for(unsigned int n=0;n<nbAtomType;n++){
			double block_charge = 0.; 
			for(unsigned int i=0;i<DNS.size();i++){
				if( DNS[i][0]-1 == n ){
					if( i == 0 ) block_charge = AtomCharge[DNS[i][0]-1];
					block_charge += AtomCharge[DNS[i][2]-1]*DNS[i][1];
					if( block_charge > zero_num ) pos = true;
					else if( block_charge < -zero_num ) neg = true;
					if( Stoich[DNS[i][2]-1] > Stoich[DNS[i][0]-1]*DNS[i][1] ){
						if( AtomCharge[DNS[i][2]-1] > zero_num ) pos = true;
						else if( AtomCharge[DNS[i][2]-1] < -zero_num ) neg = true;
					}
					if( pos && neg ){
						Possible = true;
						break;
					}
				}
			}
		}
		for(unsigned int n=0;n<nbAtomType;n++){
			bool treat = true;
			for(unsigned int i=0;i<DNS.size();i++){
				if( DNS[i][0]-1 == n ){
					treat = false;
					break;
				}
				if( DNS[i][2]-1 == n ){
					treat = false;
					break;
				}
			}
			if( treat ){
				if( AtomCharge[n] > zero_num ) pos = true;
				else if( AtomCharge[n] < -zero_num ) neg = true;
				if( pos && neg ){
					Possible = true;
					break;
				}
			}
		}
	}	

	if( !Possible ){
		cout << "It is not possible to construct charge neutral surfaces as there are not two types of ions (or blocks of ions due to the DoNotSep crystal tag) having charge with opposite signs" << endl;
		return;
	}

	cout << "Trying to make the surface charge neutral" << endl;
	
	// search surface ions and compute the initial total_charge
	double rcut = MT->max_p(_MyCrystal->getALength(),3);
	searchNeighbours(rcut);
	double *classif = new double[nbAtom];
	for(unsigned int i=0;i<nbAtom;i++)
		classif[i] = 1. / (double) Neighbours[i*(nbMaxN+1)];
	Descriptors *Des = new Descriptors(classif,nbAtom,1);
	unsigned int nb_clust_KM = 2;
	KMeans *KM = new KMeans();
	KM->setDescriptors(Des);
	KM->TrainModel(nb_clust_KM);
	unsigned int clust_ind = 0;
	for(unsigned int i=1;i<nb_clust_KM;i++)
		if( KM->getOptimalCentroid(i)[0] > KM->getOptimalCentroid(clust_ind)[0] ) clust_ind = i;
	vector<unsigned int> ions2treat;
	double total_charge = 0.;
	//ofstream filep("test");
	//for(unsigned int i=0;i<nbAtom;i++){
	//	filep << classif[i] << endl;
	//	total_charge += AtomCharge[AtomList[i].type_uint-1];
	//	if( classif[i] < 0.99*KM->getOptimalCentroid(clust_ind)[0] ) classif[i] = 0.;
	//	else{
	//		classif[i] = 1.;
	//		ions2treat.push_back(i);
	//	}
	//}
	//filep.close();
	//cout << KM->getOptimalCentroid(0)[0] << " " << KM->getOptimalCentroid(0)[0] << " " << KM->getOptimalCentroid(clust_ind)[0] << endl;
	for(unsigned int i=0;i<nbAtom;i++){
		classif[i] = 0.;
		total_charge += AtomCharge[AtomList[i].type_uint-1];
		if( KM->getKMeansTools(0)->getData2Cluster()[i] == clust_ind ){
			classif[i] = 1.;
			ions2treat.push_back(i);
		}
	}
	unsigned int nbIons2Treat = ions2treat.size();
	cout << nbIons2Treat << " surface ions will be used to try to get neutral surfaces" << endl;

	// Compute charge density at each ions to treat position
	double fac_sigma_gauss = 1.5;
	double *ChargeDens = new double[nbIons2Treat];
	double sigma_gauss = rcut;
	sigma_gauss *= fac_sigma_gauss;
	int bcx = 2; // to be computed depending on box size
	int bcy = 2;
	int bcz = 2;
	#pragma omp parallel for
	for(unsigned int i=0;i<nbIons2Treat;i++){
		ChargeDens[i] = 0.;
		for(unsigned int j=0;j<nbAtom;j++){
			for(int bc1=-bcx;bc1<=bcx;bc1++){
				double coord1 = this->AtomList[j].pos.x+(this->H1[0]*bc1);
				for(int bc2=-bcy;bc2<=bcy;bc2++){
					double coord2 = this->AtomList[j].pos.y+(this->H2[1]*bc2);
					for(int bc3=-bcz;bc3<=bcz;bc3++){
						double coord3 = this->AtomList[j].pos.z+(this->H3[2]*bc3);
						ChargeDens[i] += MT->gaussian(AtomList[ions2treat[i]].pos.x, AtomList[ions2treat[i]].pos.y, AtomList[ions2treat[i]].pos.z, coord1, coord2, coord3, sigma_gauss)*AtomCharge[AtomList[j].type_uint-1];
					}
				}
			}
		}
	}
	for(unsigned int i=0;i<nbIons2Treat;i++) classif[ions2treat[i]] = ChargeDens[i];

	setAux(classif,"class");
	printSystem_aux("TEST.cfg","class");

	// compute initial loss only if the system is charge neutral
	bool initialy_neutral = false; 
	double initial_loss;
	if( fabs(total_charge) < zero_num ){
		initialy_neutral = true;
		initial_loss = 0.;
		for(unsigned int i=0;i<nbIons2Treat;i++) initial_loss += fabs(ChargeDens[i]);
		initial_loss /= nbIons2Treat;
	}



	mt19937 gen(random_device{}());
	//mt19937 gen;
	//gen.seed(42);
	
	vector<vector<unsigned int>> ions2remove; // ions2remove[r*nbRand+n] = ions 2 remove for rcut = rcut[r] and n random
	vector<vector<unsigned int>> select_ind; // select_ind[r*nbRand+n][i] = index from which to remove element in ions2remove[r*nbRand+n] array to have the loss select_loss[r*..][i]
	vector<vector<double>> select_loss;
	unsigned int exp_prob = 2;
	unsigned int nb_rcut = 20;
	double max_facrcut = 10;
	double min_facrcut = 0.1;
	unsigned int nbRand = 20;
	unsigned int total_nbTry = nb_rcut*nbRand;
	// loop or fac_rcut should begin here
	// compute restricted neighbours
	for(unsigned int nrc=0;nrc<nb_rcut;nrc++){
	double fac_rcut = min_facrcut + (max_facrcut-min_facrcut)*nrc/(nb_rcut-1.);
	//cout << "facrcut = " << fac_rcut << endl;
	double current_rcut = rcut*fac_rcut;
	searchNeighbours_restricted(current_rcut,ions2treat,ions2treat);
	// loop on nbRand should begin here
	for(unsigned int nr=0;nr<nbRand;nr++){

		// construct a first set of possible ion to remove
		vector<vector<unsigned int>> poss2rm;
		vector<double> probs;

        	for(unsigned int i=0;i<nbIons2Treat;i++){
			unsigned int current_ind = ions2treat[i];
			bool already = false;
			for(unsigned int j1=0;j1<poss2rm.size();j1++){
				for(unsigned int j2=0;j2<poss2rm[j1].size();j2++){
					if( poss2rm[j1][j2] == current_ind ){
						already = true;
						break;
					}
				}
				if( already ) break;
			}
			if( already ) continue;
			vector<unsigned int> current_ions2remove;
			if( IsNotSepTag ){
				if( NotSepTag[current_ind][0] == 0 ) current_ions2remove.push_back(current_ind);
				else if( NotSepTag[current_ind][0] < 0 ){
					unsigned int neigh_ind = (unsigned int) -1-NotSepTag[current_ind][0];
					current_ions2remove.push_back(current_ind);
					current_ions2remove.push_back(neigh_ind);	
					for(unsigned int n=0;n<NotSepTag[neigh_ind][0];n++)
						if( current_ind != NotSepTag[neigh_ind][n+1] ) current_ions2remove.push_back(NotSepTag[neigh_ind][n+1]);
				}else{
					current_ions2remove.push_back(current_ind);
					for(unsigned int n=0;n<NotSepTag[current_ind][0];n++) current_ions2remove.push_back(NotSepTag[current_ind][n+1]);
				}
			}else current_ions2remove.push_back(current_ind);
			double block_charge = 0.;
			for(unsigned int n=0;n<current_ions2remove.size();n++) block_charge += AtomCharge[AtomList[current_ions2remove[n]].type_uint-1];

			//if( signbit(ChargeDens[i]) == signbit(block_charge) ){
			if( ions2treat[i] == 15911 ){
				cout << "CELUI LA!!" << endl;
				cout << (signbit(ChargeDens[i]) == signbit(block_charge)) << " " << (signbit(block_charge) == signbit(total_charge))  << " " << (fabs(block_charge) - fabs(total_charge) < zero_num) << endl;
				cout << block_charge << " " << total_charge << endl;
			}
			if( signbit(ChargeDens[i]) == signbit(block_charge) && signbit(block_charge) == signbit(total_charge)  && fabs(block_charge) - fabs(total_charge) < zero_num ){
				cout << "STOCKEEEE" << endl;
				poss2rm.push_back(current_ions2remove);
				probs.push_back(pow(fabs(ChargeDens[i]),exp_prob)); // maybe bias more the probability by taking ChargeDens^2,3,4? 
			}
		}


		//cout << "NRAND =" << nr << endl;
	ions2remove.push_back(vector<unsigned int>());
	select_ind.push_back(vector<unsigned int>());
	select_loss.push_back(vector<double>());
	// copy initial data, we need ChargeDensity
	
	vector<double> current_ChargeDens(nbIons2Treat);
	vector<unsigned int> current_ions2treat(nbIons2Treat);
	for(unsigned int i=0;i<nbIons2Treat;i++){
		current_ChargeDens[i] = ChargeDens[i];
		current_ions2treat[i] = ions2treat[i];
	}
	unsigned int cur_i2m = ions2remove.size()-1;

	unsigned int current_nbIons2Treat = current_ChargeDens.size();
	double current_total_charge = total_charge;
	unsigned int nbreset = 0;
	unsigned int maxnbreset = 10;
	while( probs.size() != 0 && fabs(current_total_charge) > zero_num ){ // this loop is also broke if current_nbIons2Treat == 0
	// randomly select ions to remove from poss2rm using ChargeDensity as probability
		unsigned int nb_continue = 0;
		while( probs.size() != 0 && fabs(current_total_charge) > zero_num ){
    			discrete_distribution<size_t> d{probs.begin(), probs.end()};
			unsigned int ind2rm = d(gen);

			// help to got neutral system
			double block_charge = 0.;	
			if( fabs(current_total_charge) > zero_num ){
				for(unsigned int i=0;i<poss2rm[ind2rm].size();i++)
					block_charge += AtomCharge[AtomList[poss2rm[ind2rm][i]].type_uint-1];
				if( signbit(current_total_charge) != signbit(block_charge) ){
					nb_continue++;
					if( nb_continue < poss2rm.size() )
						continue;
					else
						probs.clear();
						poss2rm.clear();
						nbreset++;
						cout << "resseting " << nbreset << endl;
						break;
				}
			}
			nb_continue = 0;

			// remove these ions from ions2treat and ChargeDens lists
			for(unsigned int i=0;i<poss2rm[ind2rm].size();i++){
				// add these ions to the ion2remove list
				ions2remove[cur_i2m].push_back(poss2rm[ind2rm][i]);
				// remove it from the current_ChargeDens and current_ions2treat arrays
				bool found = false;
				unsigned int curind2rm;
				for(unsigned int n=0;n<current_ions2treat.size();n++){
					if( poss2rm[ind2rm][i] == current_ions2treat[n] ){
						curind2rm = n;
						found = true;
						break;
					}
				}
				if( found ){
					current_ions2treat.erase(current_ions2treat.begin()+curind2rm);
					current_ChargeDens.erase(current_ChargeDens.begin()+curind2rm);
				}
			}
			current_nbIons2Treat = current_ions2treat.size();

			// remove contribution of these ions to charge density and total_charge
			for(unsigned int i=0;i<poss2rm[ind2rm].size();i++){
				current_total_charge -= AtomCharge[AtomList[poss2rm[ind2rm][i]].type_uint-1];
				for(unsigned int n=0;n<current_ChargeDens.size();n++){
					for(int bc1=-bcx;bc1<=bcx;bc1++){
						double coord1 = this->AtomList[poss2rm[ind2rm][i]].pos.x+(this->H1[0]*bc1);
						for(int bc2=-bcx;bc2<=bcx;bc2++){
							double coord2 = this->AtomList[poss2rm[ind2rm][i]].pos.y+(this->H2[1]*bc2);
							for(int bc3=-bcx;bc3<=bcx;bc3++){
								double coord3 = this->AtomList[poss2rm[ind2rm][i]].pos.z+(this->H3[2]*bc3);
								current_ChargeDens[n] -= MT->gaussian(AtomList[current_ions2treat[n]].pos.x, AtomList[current_ions2treat[n]].pos.y, AtomList[current_ions2treat[n]].pos.z, coord1, coord2, coord3, sigma_gauss)*AtomCharge[AtomList[poss2rm[ind2rm][i]].type_uint-1];
							}
						}
					}
				}
			}

			// recompute loss and store indexes if total charge = 0
			if( fabs(current_total_charge) < zero_num ){
				select_ind[cur_i2m].push_back(ions2remove[cur_i2m].size());
				select_loss[cur_i2m].push_back(0.);
				for(unsigned int n=0;n<current_ChargeDens.size();n++) select_loss[cur_i2m][select_loss[cur_i2m].size()-1] += fabs(current_ChargeDens[n]);
				select_loss[cur_i2m][select_loss[cur_i2m].size()-1] /= current_nbIons2Treat;
			}

			// remove ions and neighboring ions from the probs and poss2rm list
			unsigned int ind_in_N;
		        bool found = false;
			for(unsigned int n=0;n<ions2treat.size();n++){
				if( poss2rm[ind2rm][0] == ions2treat[n] ){
					found = true;
					ind_in_N = n;
					break;
				}
			}
			if( !found ) cout << "OTHER ISSUE IN NEIGH!!" << endl;
			poss2rm.erase(poss2rm.begin()+ind2rm);
			probs.erase(probs.begin()+ind2rm);
			for(unsigned int n=0;n<Neighbours[ind_in_N*(nbMaxN+1)];n++){
				unsigned int Nind = Neighbours[ind_in_N*(nbMaxN+1)+1+n];
				unsigned int n2rm;
				found = false;
				for(unsigned int j1=0;j1<poss2rm.size();j1++){
					for(unsigned int j2=0;j2<poss2rm[j1].size();j2++){
						if( Nind == poss2rm[j1][j2] ){
							n2rm = j1;
							found = true;
							break;
						}
					}
					if( found ) break;
				}
				if( found ){
					poss2rm.erase(poss2rm.begin()+n2rm);
					probs.erase(probs.begin()+n2rm);
				}
			}
			// maybe update probs here ?
			for(unsigned int i=0;i<probs.size();i++){
				bool found = false;
				for(unsigned int j=0;j<current_ions2treat.size();j++){
					if( poss2rm[i][0] == current_ions2treat[j] ){
						found = true;
						probs[i] = pow(fabs(current_ChargeDens[j]),exp_prob);
						break;
					}
				}
			}
		} // end while probs.size() != 0
		// compute new poss2rm and probs lists if current_nbIons2Treat != 0 
		if( current_nbIons2Treat == 0 || nbreset == maxnbreset ) break;
        	for(unsigned int i=0;i<current_nbIons2Treat;i++){
			unsigned int current_ind = current_ions2treat[i];
			bool already = false;
			for(unsigned int j1=0;j1<poss2rm.size();j1++){
				for(unsigned int j2=0;j2<poss2rm[j1].size();j2++){
					if( poss2rm[j1][j2] == current_ind ){
						already = true;
						break;
					}
				}
				if( already ) break;
			}
			if( already ) continue;
			vector<unsigned int> current_ions2remove;
			if( IsNotSepTag ){
				if( NotSepTag[current_ind][0] == 0 ) current_ions2remove.push_back(current_ind);
				else if( NotSepTag[current_ind][0] < 0 ){
					unsigned int neigh_ind = (unsigned int) -1-NotSepTag[current_ind][0];
					current_ions2remove.push_back(current_ind);
					current_ions2remove.push_back(neigh_ind);	
					for(unsigned int n=0;n<NotSepTag[neigh_ind][0];n++)
						if( current_ind != NotSepTag[neigh_ind][n+1] ) current_ions2remove.push_back(NotSepTag[neigh_ind][n+1]);
				}else{
					current_ions2remove.push_back(current_ind);
					for(unsigned int n=0;n<NotSepTag[current_ind][0];n++) current_ions2remove.push_back(NotSepTag[current_ind][n+1]);
				}
			}else current_ions2remove.push_back(current_ind);
			double block_charge = 0.;
			for(unsigned int n=0;n<current_ions2remove.size();n++) block_charge += AtomCharge[AtomList[current_ions2remove[n]].type_uint-1];

			if( signbit(ChargeDens[i]) == signbit(block_charge) && signbit(block_charge) == signbit(total_charge)  && fabs(block_charge) - fabs(total_charge) < zero_num ){
				poss2rm.push_back(current_ions2remove);
				probs.push_back(pow(fabs(current_ChargeDens[i]),exp_prob)); // maybe bias more the probability by taking ChargeDens^2,3,4? 
			}
		}
	} // end while probs.size() != 0
	  cout << "At the end of try " << cur_i2m << ", still " << current_nbIons2Treat << " ions 2 treat" << endl;
	}}
	// end loops on rcut and nrand
	
	// choose the optimal system
	bool neutral_found = false;
	for(unsigned int n=0;n<total_nbTry;n++){
		if( select_ind[n].size() != 0 ){
			neutral_found = true;
			break;
		}
	}
	
	unsigned int opt_ind_glob, opt_ind;
	if( initialy_neutral && !neutral_found ){
		cout << "Only initial system is neutral, keeping it" << endl;
		return;
	}else if( !initialy_neutral && !neutral_found ){
		cout << "Cannot find neutral system" << endl;
		return;
	}else{
		double min_loss = std::numeric_limits<double>::max();
		ofstream file("lossfct.dat");
		unsigned int glob_count = 0;
		for(unsigned int n1=0;n1<total_nbTry;n1++){
			for(unsigned int n2=0;n2<select_loss[n1].size();n2++){
				file << glob_count << " " << select_loss[n1][n2] << endl;
				glob_count++;
				if( select_loss[n1][n2] < min_loss ){
					opt_ind_glob = n1;
					opt_ind = n2;
					min_loss = select_loss[n1][n2];
				}
			}
		}
		file.close();
		if( initialy_neutral && initial_loss < min_loss )
			cout << "The initial system is better than all investigated systems, keeping it intact" << endl;
		else{
			unsigned int init_size = ions2remove[opt_ind_glob].size();
			cout << "A better system have been found by removing " << select_ind[opt_ind_glob][opt_ind] << " over " << ions2remove[opt_ind_glob].size() << endl;
			for(unsigned int i=select_ind[opt_ind_glob][opt_ind];i<init_size;i++) ions2remove[opt_ind_glob].pop_back(); // to test
			RemoveAtoms(ions2remove[opt_ind_glob]);
		}
	}

	//file.close();
	//file.open("distribQ_1.dat");
	//for(unsigned int i=0;i<nbAtom;i++) file << ChargeDens[i] << endl;
	//file.close();
	//cout << "nb ions rem = " << ions2remove.size() << ", total charge = " << total_charge << endl;
	//setAux(classif,"class");
	//setAux(removed,"rm");
	//setAux(ChargeDens,"ChargeDensAfter");
	//printSystem_aux("TESTCharge.cfg","ChargeDensBefore classbef class ChargeDensAfter rm");
	//
	delete[] classif;
	delete Des;
	delete KM;
	delete[] ChargeDens;
}

void AtomicSystem::RemoveAtoms(vector<unsigned int> index2rm, bool DNS){
	//if( !IsAtomListMine ){
	//	cout << "The atom list does not belong to the atomic system, cannot delete atoms" << endl;
	//	return;
	//}
	if( IsWrappedPos ){
		delete[] WrappedPos;
		WrappedPos = nullptr;
		IsWrappedPos = false;
	}
	if( IsNeighbours ){
		delete[] Neighbours;
		Neighbours = nullptr;
		delete[] CLNeighbours;
		CLNeighbours = nullptr;
		IsNeighbours = false;
	}
	if( DNS && !_MyCrystal->getIsDoNotSep() ) DNS = false;
	else if( DNS && !IsNotSepTag ){
		cout << "Computing not sep list.." << endl;
		ComputeNotSepList();
	}


	//unsigned int newNbAtom = nbAtom-index2rm.size();
	//Atom *tempAtomList = new Atom[newNbAtom];
	//vector<double*> tempAux;
	//for(unsigned int i=0;i<Aux_size.size();i++) tempAux.push_back(new double[newNbAtom*Aux_size[i]]);
	vector<Atom> tempAtomList;
	vector<vector<double>> tempAux(Aux_size.size(),vector<double>());

	// copy vector
	vector<unsigned int> current_index2rm;
	for(unsigned int i=0;i<index2rm.size();i++) current_index2rm.push_back(index2rm[i]);

	// adjust the list to account for the do not separe list
	if( DNS ){
		for(unsigned int i=0;i<index2rm.size();i++){
			if( NotSepTag[index2rm[i]][0] == 0 ) continue;
			else if( NotSepTag[index2rm[i]][0] < 0 ){
				bool toadd = true;
				unsigned int neigh_ind = (unsigned int) -1-NotSepTag[index2rm[i]][0];
				for(unsigned int j=0;j<index2rm.size();j++){
					if( index2rm[j] == neigh_ind ){
						toadd = false;
						break;
					}
				}
				if( toadd ) current_index2rm.push_back(neigh_ind);
				for(unsigned int n=0;n<NotSepTag[neigh_ind][0];n++){
					if( index2rm[i] == NotSepTag[neigh_ind][n+1] ) continue;
					else{
						toadd = true;
						for(unsigned int j=0;j<index2rm.size();j++){
							if( index2rm[j] == NotSepTag[neigh_ind][n+1] ){
								toadd = false;
								break;
							}
						}
						if( toadd ) current_index2rm.push_back(NotSepTag[neigh_ind][n+1]);
					}
				}
			}else{
				for(unsigned int n=0;n<NotSepTag[index2rm[i]][0];n++){
					bool toadd = true;
					for(unsigned int j=0;j<index2rm.size();j++){
						if( index2rm[j] == NotSepTag[index2rm[i]][n+1] ){
							toadd = false;
							break;
						}
					}
					if( toadd ) current_index2rm.push_back(NotSepTag[index2rm[i]][n+1]);
				}
			}
		}
	}
	unsigned int count=0;
	for(unsigned int i=0;i<nbAtom;i++){
		bool treat = true;
		unsigned int j;
		for(j=0;j<current_index2rm.size();j++){
			if( current_index2rm[j] == i ){
				treat = false;
				break;
			}
		}
		if( treat ){
			tempAtomList.push_back(AtomList[i]);
			for(unsigned int a=0;a<Aux_size.size();a++){
				for(unsigned int d=0;d<Aux_size[a];d++) tempAux[a].push_back(Aux[a][i*Aux_size[a]+d]);
			}
			count++;
		}else current_index2rm.erase(current_index2rm.begin()+j);
	}

	delete[] AtomList;
	for(unsigned int a=0;a<Aux_size.size();a++) delete[] Aux[a];
	Aux.clear();
	
	nbAtom = count;
	AtomList = new Atom[nbAtom];
	for(unsigned int i=0;i<Aux_size.size();i++) Aux.push_back(new double[nbAtom*Aux_size[i]]);
	for(unsigned int i=0;i<nbAtom;i++){
		AtomList[i] = tempAtomList[i];
		for(unsigned int a=0;a<Aux_size.size();a++){
			for(unsigned int d=0;d<Aux_size[a];d++) Aux[a][i*Aux_size[a]+d] = tempAux[a][i*Aux_size[a]+d];
		}
	}

	if( IsNotSepTag ){
		delete[] NotSepTag;
		NotSepTag = nullptr;
		IsNotSepTag = false;
	}

	//delete[] tempAtomList;
	//for(unsigned int a=0;a<Aux_size.size();a++) delete[] tempAux[a];
	//tempAux.clear();
		
	// TODO
	//unsigned int *tempMolId, *tempBonds, *tempAngles, *tempBondType, *tempAngleType;
	//int *tempPeriodicArr;
	//
	
}

void AtomicSystem::MakeSystemNeutral(){
	if( IsCharge == false ) return;
	bool Possible = false; // at least we should have two atom types having charge with opposite sign
	bool pos(false),neg(false);
	double zero_num = 1e-8;
	if( !_MyCrystal->getIsDoNotSep() ){
		for(unsigned int n=0;n<nbAtomType;n++){
			if( AtomCharge[n] > zero_num ) pos = true;
			else if( AtomCharge[n] < -zero_num ) neg = true;
			if( pos && neg ){
				Possible = true;
				break;
			}
		}
	}else{
		ComputeNotSepList();
		vector<vector<unsigned int>> DNS = _MyCrystal->getDoNotSep();
		unsigned int *Stoich = _MyCrystal->getStoich();
		for(unsigned int n=0;n<nbAtomType;n++){
			double block_charge = 0.; 
			for(unsigned int i=0;i<DNS.size();i++){
				if( DNS[i][0]-1 == n ){
					if( i == 0 ) block_charge = AtomCharge[DNS[i][0]-1];
					block_charge += AtomCharge[DNS[i][2]-1]*DNS[i][1];
					if( block_charge > zero_num ) pos = true;
					else if( block_charge < -zero_num ) neg = true;
					if( Stoich[DNS[i][2]-1] > Stoich[DNS[i][0]-1]*DNS[i][1] ){
						if( AtomCharge[DNS[i][2]-1] > zero_num ) pos = true;
						else if( AtomCharge[DNS[i][2]-1] < -zero_num ) neg = true;
					}
					if( pos && neg ){
						Possible = true;
						break;
					}
				}
			}
		}
		for(unsigned int n=0;n<nbAtomType;n++){
			bool treat = true;
			for(unsigned int i=0;i<DNS.size();i++){
				if( DNS[i][0]-1 == n ){
					treat = false;
					break;
				}
				if( DNS[i][2]-1 == n ){
					treat = false;
					break;
				}
			}
			if( treat ){
				if( AtomCharge[n] > zero_num ) pos = true;
				else if( AtomCharge[n] < -zero_num ) neg = true;
				if( pos && neg ){
					Possible = true;
					break;
				}
			}
		}
	}	

	if( !Possible ){
		cout << "It is not possible to construct charge neutral surfaces as there are not two types of ions (or blocks of ions due to the DoNotSep crystal tag) having charge with opposite signs" << endl;
		return;
	}
	
	vector<unsigned int> ions2treat;
	double total_charge = 0.;
	for(unsigned int i=0;i<nbAtom;i++) total_charge += AtomCharge[AtomList[i].type_uint-1];
	
	// Compute charge density at each ions position
	double fac_sigma_gauss = 1.5;
	double *ChargeDens = new double[nbAtom];
	double sigma_gauss = MT->max_p(_MyCrystal->getALength(),3);
	sigma_gauss *= fac_sigma_gauss;
	int bcx = 2; // to be computed depending on box size
	int bcy = 2;
	int bcz = 2;
	double initial_metric = 0.;
	#pragma omp parallel for
	for(unsigned int i=0;i<nbAtom;i++){
		ChargeDens[i] = 0.;
		for(unsigned int j=0;j<nbAtom;j++){
			for(int bc1=-bcx;bc1<=bcx;bc1++){
				double coord1 = this->AtomList[j].pos.x+(this->H1[0]*bc1);
				for(int bc2=-bcy;bc2<=bcy;bc2++){
					double coord2 = this->AtomList[j].pos.y+(this->H2[1]*bc2);
					for(int bc3=-bcz;bc3<=bcz;bc3++){
						double coord3 = this->AtomList[j].pos.z+(this->H3[2]*bc3);
						ChargeDens[i] += MT->gaussian(AtomList[ions2treat[i]].pos.x, AtomList[ions2treat[i]].pos.y, AtomList[ions2treat[i]].pos.z, coord1, coord2, coord3, sigma_gauss)*AtomCharge[AtomList[j].type_uint-1];
					}
				}
			}
		}
		initial_metric += fabs(ChargeDens[i]);
	}
	initial_metric /= nbAtom;

	//mt19937 gen(random_device{}());
	mt19937 gen;
	gen.seed(42);

	unsigned int iter_totC = 0;
	unsigned int maxIter_totC = 100;

	//while( fabs(total_charge) > 1e-8 && iter_totC < maxIter_totC ){
	//	vector<vector<unsigned int>> poss2rm;
	//	vector<double> probs;
	//	// search ions having ChargeDens of same sign than both block_charge and total_charge
	//	for(unsigned int i=0;i<nbAtom;i++){
	//		vector<unsigned int> current_ions2remove;
	//		if( IsNotSepTag ){
	//			if( NotSepTag[i][0] == 0 ) current_ions2remove.push_back(i);
	//			else if( NotSepTag[i][0] < 0 ){
	//				unsigned int neigh_ind = (unsigned int) -1-NotSepTag[i][0];
	//				current_ions2remove.push_back(i);
	//				current_ions2remove.push_back(neigh_ind);	
	//				for(unsigned int n=0;n<NotSepTag[neigh_ind][0];n++)
	//					if( i != NotSepTag[neigh_ind][n+1] ) current_ions2remove.push_back(NotSepTag[neigh_ind][n+1]);
	//			}else{
	//				current_ions2remove.push_back(i);
	//				for(unsigned int n=0;n<NotSepTag[i][0];n++) current_ions2remove.push_back(NotSepTag[i][n+1]);
	//			}
	//		}else current_ions2remove.push_back(i);
	//		double block_charge = 0.;
	//		for(unsigned int n=0;n<current_ions2remove.size();n++) block_charge += AtomCharge[AtomList[current_ions2remove[n]].type_uint-1];

	//		if( signbit(ChargeDens[i]) == signbit(block_charge) && signbit(block_charge) == signbit(total_charge)  && fabs(block_charge) <= fabs(total_charge) ){
	//			poss2rm.push_back(current_ions2remove);
	//			double mean_charge_dens = 0.;
	//			for(unsigned int n=0;n<current_ions2remove.size();n++) mean_charge_dens += fabs(ChargeDens[current_ions2remove[n]]);
	//			mean_charge_dens /= current_ions2remove.size();
	//			probs.push_back(mean_charge_dens);
	//		}
	//	}

    	//	discrete_distribution<size_t> d{probs.begin(), probs.end()};
	//	unsigned int ind2rm = d(gen);

	//	// remove these ions from ions2treat and ChargeDens lists
	//	for(unsigned int i=0;i<poss2rm[ind2rm].size();i++){
	//		// add these ions to the ion2remove list
	//		ions2remove[cur_i2m].push_back(poss2rm[ind2rm][i]);
	//		// remove it from the current_ChargeDens and current_ions2treat arrays
	//		bool found = false;
	//		unsigned int curind2rm;
	//		for(unsigned int n=0;n<current_ions2treat.size();n++){
	//			if( poss2rm[ind2rm][i] == current_ions2treat[n] ){
	//				curind2rm = n;
	//				found = true;
	//				break;
	//			}
	//		}
	//		if( found ){
	//			current_ions2treat.erase(current_ions2treat.begin()+curind2rm);
	//			current_ChargeDens.erase(current_ChargeDens.begin()+curind2rm);
	//		}
	//	}
	//	current_nbIons2Treat = current_ions2treat.size();

	//} // end while total_charge != 0

}

void AtomicSystem::MakeSurfaceNeutral(vector<int> Oris, vector<double> shift, vector<double> Misfit, vector<double> PlaneNormal, double VertPlaneSlope, vector<double> VertPlaneEq, vector<double> FullPlaneEq, double shift_z_box, string outputfilename){
	double zero = 0.;
	unsigned int nbSample_S_z = 5;
	double max_s_z = 4.;
	vector<vector<unsigned int>> opt_index2rm(nbSample_S_z,vector<unsigned int>());
	#pragma omp parallel for
	for(unsigned int s_z_i=0;s_z_i<nbSample_S_z;s_z_i++){
	double s_z = max_s_z*(double) s_z_i/((double) nbSample_S_z - 1.);
	Crystal OriCryst_f1(_MyCrystal->getName());
	OriCryst_f1.RotateCrystal(Oris[0], Oris[1], Oris[2], Oris[3], Oris[4], Oris[5]);
	OriCryst_f1.ShiftMotif(shift[0], shift[1], shift[2]);
	OriCryst_f1.ConstructOrthogonalCell();
	AtomicSystem *OriSys_f1 = OriCryst_f1.getOrientedSystem();
	double *invTiltTrans = new double[9];
	double *invRotMat = new double[9];
	MT->invert3x3(OriCryst_f1.getTiltTrans(),invTiltTrans);
	MT->invert3x3(OriCryst_f1.getRotMat(),invRotMat);
	OriSys_f1->ApplyShift(zero, zero, s_z);
	double z_shift_OriSys = OriSys_f1->MakeSurfaceNeutral(false);
	if( outputfilename != "" ) OriSys_f1->print_lmp(outputfilename);
	// keep only a small part of the system near to the surface to limit computational cost
	double select_width = 10.;
	double zero_num = 1.e-8;
	bool xpar = false; // TODO
	double hi_nb = 1.e8;
	// test whether we should keep the bottom or the top surface to be aligned with the normal plane
	double *zaxis = new double[3];
	zaxis[0] = 0.;
	zaxis[1] = 0.;
	zaxis[2] = 1.;
	MT->MatDotRawVec(invTiltTrans,zaxis,zaxis);
	MT->MatDotRawVec(invRotMat,zaxis,zaxis);
	MT->MatDotRawVec(_MyCrystal->getRotMat(),zaxis,zaxis);
	MT->MatDotRawVec(_MyCrystal->getTiltTrans(),zaxis,zaxis);
        double sp_z = 0.;
	for(unsigned int i=0;i<3;i++) sp_z += zaxis[i]*PlaneNormal[i];
	vector<unsigned int> index2rm;
	if( sp_z > 0 ) index2rm = OriSys_f1->selectAtomInBox(-hi_nb, hi_nb, -hi_nb, hi_nb, -hi_nb, OriSys_f1->getH3()[2]-select_width, false);
	else index2rm = OriSys_f1->selectAtomInBox(-hi_nb, hi_nb, -hi_nb, hi_nb, select_width, hi_nb, false);
	OriSys_f1->RemoveAtoms(index2rm);
	unsigned int nbAt_OriSys = OriSys_f1->getNbAtom();
	// remove z shift from the MakeSurfaceNeutral function
	OriSys_f1->ApplyShift(zero, zero, z_shift_OriSys-s_z);
	// ** Set the curent grain def and misfit to the OriSys ** //
	// 1. remove the crystal def to the ori sys
	for(unsigned int i=0;i<nbAt_OriSys;i++) MT->MatDotAt(invTiltTrans, OriSys_f1->getAtomRef(i), OriSys_f1->getAtomRef(i));
	for(unsigned int i=0;i<nbAt_OriSys;i++) MT->MatDotAt(invRotMat, OriSys_f1->getAtomRef(i), OriSys_f1->getAtomRef(i));
	// 2. Rotate the Ori sys for the plane normal to be aligned with the facet normal 
	for(unsigned int i=0;i<nbAt_OriSys;i++) MT->MatDotAt(_MyCrystal->getRotMat(), OriSys_f1->getAtomRef(i), OriSys_f1->getAtomRef(i));
	// 3. Apply grain1 def
	for(unsigned int i=0;i<nbAt_OriSys;i++) MT->MatDotAt(_MyCrystal->getTiltTrans(), OriSys_f1->getAtomRef(i), OriSys_f1->getAtomRef(i));
	// 4. Apply misfit
	double *mean_pos_orisys = new double[3];
	for(unsigned int i=0;i<3;i++) mean_pos_orisys[i] = 0.;
	for(unsigned int i=0;i<nbAt_OriSys;i++){
		OriSys_f1->getAtomRef(i).pos.x *= Misfit[0];
		OriSys_f1->getAtomRef(i).pos.y *= Misfit[1];
		OriSys_f1->getAtomRef(i).pos.z *= Misfit[2];
		mean_pos_orisys[0] += OriSys_f1->getAtomRef(i).pos.x;
		mean_pos_orisys[1] += OriSys_f1->getAtomRef(i).pos.y;
		mean_pos_orisys[2] += OriSys_f1->getAtomRef(i).pos.z;
	}
	for(unsigned int i=0;i<3;i++) mean_pos_orisys[i] /= nbAt_OriSys;

	// ** Make coincide the ori system with the atomic system ** //
	// Set a subsystem of the atomic system
	double para_width = 2.*OriSys_f1->getH2()[1];
	double *mean_pos_subsys = new double[3];
	for(unsigned int i=0;i<3;i++) mean_pos_subsys[i] = 0.;

	double *diff_orisubsys = new double[3];
	double *at_sp = new double[nbAt_OriSys];
	double *current_shift = new double[3];
	double *swap_shift = new double[3];

	// loop on number of plane 
	for(unsigned int f=0;f<FullPlaneEq.size();f++){
	//for(unsigned int f=0;f<1;f++){
		vector<Atom> SubSys;
		for(unsigned int i=0;i<3;i++) mean_pos_subsys[i] = 0.;
		for(unsigned int i=0;i<nbAtom;i++){
			double xpos = AtomList[i].pos.x;
			double ypos = AtomList[i].pos.y;
			double zpos = AtomList[i].pos.z;
			double buffer_v;
			if( xpar ) buffer_v = xpos;
			else buffer_v = ypos - VertPlaneSlope*xpos;
			if( ( buffer_v > VertPlaneEq[f*2] ) && ( buffer_v < VertPlaneEq[2*f+1] ) ){
				// Then test if atom is bellow the real plane
				double buffer_f = xpos*PlaneNormal[0] + ypos*PlaneNormal[1] + (zpos-shift_z_box)*PlaneNormal[2];
				if( buffer_f > FullPlaneEq[f]-(select_width*2.) ){
					SubSys.push_back(AtomList[i]);
					mean_pos_subsys[0] += AtomList[i].pos.x;
					mean_pos_subsys[1] += AtomList[i].pos.y;
					mean_pos_subsys[2] += AtomList[i].pos.z;
				}
			}
		}

		// shift orisys using cell vec to be inside the subsys
		for(unsigned int i=0;i<3;i++){
			mean_pos_subsys[i] /= SubSys.size();
			diff_orisubsys[i] = mean_pos_subsys[i] - mean_pos_orisys[i];
		}
		int max_hkl_search = 150;
		vector<unsigned int> search_hkl;
		for(unsigned int i=0;i<3;i++){
			if( fabs(_MyCrystal->getA1()[i]) > zero_num ) search_hkl.push_back(round(fabs(diff_orisubsys[i]/_MyCrystal->getA1()[i])));
			if( fabs(_MyCrystal->getA2()[i]) > zero_num ) search_hkl.push_back(round(fabs(diff_orisubsys[i]/_MyCrystal->getA2()[i])));
			if( fabs(_MyCrystal->getA3()[i]) > zero_num ) search_hkl.push_back(round(fabs(diff_orisubsys[i]/_MyCrystal->getA3()[i])));
		}
		int hkl_search = (int) MT->max_vec(search_hkl)+3;
		if( hkl_search > max_hkl_search ) hkl_search = max_hkl_search;
		vector<double> optshift(3);
		double min_dist = std::numeric_limits<double>::max();
		for(int i=-hkl_search;i<=hkl_search;i++){
			for(int j=-hkl_search;j<=hkl_search;j++){
				for(int k=-hkl_search;k<=hkl_search;k++){
					double cur_dist = 0.;
					for(unsigned int d=0;d<3;d++){
						double cur_vec = _MyCrystal->getA1()[d]*i + _MyCrystal->getA2()[d]*j + _MyCrystal->getA3()[d]*k;
						cur_dist += (cur_vec-diff_orisubsys[d])*(cur_vec-diff_orisubsys[d]); 
					}
					if( cur_dist < min_dist ){
						for(unsigned int d=0;d<3;d++) optshift[d] = _MyCrystal->getA1()[d]*i + _MyCrystal->getA2()[d]*j + _MyCrystal->getA3()[d]*k;
						min_dist = cur_dist;
					}
				}
			}
		}
		OriSys_f1->ApplyShift(optshift[0],optshift[1],optshift[2]);
		for(unsigned int i=0;i<3;i++) mean_pos_orisys[i] += optshift[i];


		double sigma = 0.25;
		double tol = 1.e-3;
		double init_sp_sum = 0.;
		double min_sp = 1.e-1;
		bool toprint = false;
		for(unsigned int i=0;i<nbAt_OriSys;i++){
			at_sp[i] = 0.;
			for(unsigned int j=0;j<SubSys.size();j++) at_sp[i] += MT->gaussian(OriSys_f1->getAtomRef(i).pos.x, OriSys_f1->getAtomRef(i).pos.y, OriSys_f1->getAtomRef(i).pos.z, SubSys[j].pos.x, SubSys[j].pos.y, SubSys[j].pos.z, sigma);
			if( ( i != 0 && fabs(at_sp[i]-at_sp[i-1]) > tol ) || at_sp[i] < min_sp ){
				cout << "WARNING the oriented system seems to not coincide with grain 1" << endl;
				toprint = true;
			}
			init_sp_sum += at_sp[i];
		}
		if( toprint ){
			cout << diff_orisubsys[0] << " "<< diff_orisubsys[1] << " "<< diff_orisubsys[2] << endl; 
			// BEGIN  PRINTING
			//OriSys->setAux(at_sp,"foraux");
			//OriSys->printSystem_aux("SP.cfg","foraux");
			unsigned int sizesys = SubSys.size()+nbAt_OriSys;
			//unsigned int sizesys = nbAtom1+nbAt_OriSys;
			Atom *AtList_temp = new Atom[sizesys];
			double *foraux = new double[sizesys];
			for(unsigned int i=0;i<SubSys.size();i++){
				AtList_temp[i] = SubSys[i];
				foraux[i] = 0.;
			}
			//for(unsigned int i=0;i<nbAtom1;i++){
			//	AtList_temp[i] = Grain1->getAtom(i);
			//	foraux[i] = 0.;
			//}
			for(unsigned int i=0;i<nbAt_OriSys;i++){
				AtList_temp[i+SubSys.size()] = OriSys_f1->getAtom(i);
				foraux[i+SubSys.size()] = 1.;
				//AtList_temp[i+nbAtom1] = OriSys_f1->getAtom(i);
				//foraux[i+nbAtom1] = 1.;
			}
			AtomicSystem *AtSysTemp = new AtomicSystem(AtList_temp,sizesys,_MyCrystal,this->H1,this->H2,this->H3);
			AtSysTemp->setAux(foraux,"foraux");
			AtSysTemp->printSystem_aux("FORAUX_"+to_string(f)+"_"+to_string(s_z_i)+".cfg","foraux");
			delete[] AtList_temp;
			delete[] foraux;
			delete AtSysTemp;
			// END PRINTING


		}
		// shift the OriSys using _Crystal cell vecs and choose the fully coinciding system being the closest to the surface
		hkl_search = 10;
		vector<int> hkl_shift;
		vector<double> z_shift;
		for(int i=-hkl_search;i<=hkl_search;i++){
			for(int j=-hkl_search;j<=hkl_search;j++){
				for(int k=-hkl_search;k<=hkl_search;k++){
					for(unsigned int d=0;d<3;d++) current_shift[d] = _MyCrystal->getA1()[d]*i + _MyCrystal->getA2()[d]*j + _MyCrystal->getA3()[d]*k;
				        OriSys_f1->ApplyShift(current_shift[0], current_shift[1], current_shift[2]);
					double current_sp_sum = 0.;
					for(unsigned int n1=0;n1<nbAt_OriSys;n1++){
						at_sp[n1] = 0.;
						for(unsigned int n2=0;n2<SubSys.size();n2++) at_sp[n1] += MT->gaussian(OriSys_f1->getAtomRef(n1).pos.x, OriSys_f1->getAtomRef(n1).pos.y, OriSys_f1->getAtomRef(n1).pos.z, SubSys[n2].pos.x, SubSys[n2].pos.y, SubSys[n2].pos.z, sigma);
						current_sp_sum += at_sp[n1];
					}
					if( fabs(current_sp_sum-init_sp_sum) < tol ){
						hkl_shift.push_back(i);
						hkl_shift.push_back(j);
						hkl_shift.push_back(k);
						double buffer_f = current_shift[0]*PlaneNormal[0] + current_shift[1]*PlaneNormal[1] + (current_shift[2]-shift_z_box)*PlaneNormal[2];
						z_shift.push_back(-buffer_f+FullPlaneEq[f]);
					}
				        OriSys_f1->ApplyShift(-current_shift[0], -current_shift[1], -current_shift[2]);
				}
			}
		}
		if( z_shift.size() == 0 ){
			cerr << "An issue occured when translating the oriented system to search the optimal postion with grain 1" << endl;
			exit(EXIT_FAILURE);
		}
		unsigned int opt_ind = MT->min(z_shift);
		for(unsigned int d=0;d<3;d++) current_shift[d] = _MyCrystal->getA1()[d]*hkl_shift[opt_ind*3] + _MyCrystal->getA2()[d]*hkl_shift[opt_ind*3+1] + _MyCrystal->getA3()[d]*hkl_shift[opt_ind*3+2];
		OriSys_f1->ApplyShift(current_shift[0], current_shift[1], current_shift[2]);
		for(unsigned int i=0;i<3;i++) mean_pos_orisys[i] += current_shift[i];

		

		// now the ori sys is well positionned
		// ** swap the surface of the system and remove the ions not coinciding with ori sys ** //
		//
		for(unsigned int i=0;i<3;i++) swap_shift[i] = OriSys_f1->getH2()[i];
		MT->MatDotRawVec(invTiltTrans,swap_shift,swap_shift);
		MT->MatDotRawVec(invRotMat,swap_shift,swap_shift);
		MT->MatDotRawVec(_MyCrystal->getRotMat(),swap_shift,swap_shift);
		MT->MatDotRawVec(_MyCrystal->getTiltTrans(),swap_shift,swap_shift);
		if( swap_shift[1] < 0 ) for(unsigned int d=0;d<3;d++) swap_shift[d] *= -Misfit[d];
		else for(unsigned int d=0;d<3;d++) swap_shift[d] *= Misfit[d];
		// put the ori sys at minimum y pos
		double min_y = std::numeric_limits<double>::max();
		double mean_z_pos = 0.;
		for(unsigned int i=0;i<nbAt_OriSys;i++){
			double xpos = OriSys_f1->getAtom(i).pos.x;
			double ypos = OriSys_f1->getAtom(i).pos.y;
			double zpos = OriSys_f1->getAtom(i).pos.z;
			double buffer_f = xpos*PlaneNormal[0] + ypos*PlaneNormal[1] + (zpos-shift_z_box)*PlaneNormal[2];
			mean_z_pos += -FullPlaneEq[f]-buffer_f;
			if( OriSys_f1->getAtom(i).pos.y < min_y ) min_y = OriSys_f1->getAtom(i).pos.y;
		}
		mean_z_pos /= nbAt_OriSys;
		// min_y should be at most VertPlaneEq[2*f]
		for(unsigned int i=0;i<3;i++) current_shift[i] = swap_shift[i]*((double) ceil((min_y-VertPlaneEq[2*f])/swap_shift[1]) + 1.);
		OriSys_f1->ApplyShift(-current_shift[0], -current_shift[1], -current_shift[2]);
		for(unsigned int i=0;i<3;i++) mean_pos_orisys[i] -= current_shift[i];
		unsigned int nb_shift = round((VertPlaneEq[2*f+1]-VertPlaneEq[2*f]+fabs(min_y-VertPlaneEq[2*f]))/swap_shift[1])+2; //TODO warning for not 2D


		index2rm.clear();
		double tol_rmatom = 0.9*init_sp_sum/nbAt_OriSys;
		for(unsigned int i=0;i<nbAtom;i++){
			double xpos = AtomList[i].pos.x;
			double ypos = AtomList[i].pos.y;
			double zpos = AtomList[i].pos.z;
			double buffer_v;
			if( xpar ) buffer_v = xpos;
			else buffer_v = ypos - VertPlaneSlope*xpos;
			if( ( buffer_v > VertPlaneEq[f*2] ) && ( buffer_v < VertPlaneEq[2*f+1] ) ){
				double buffer_f = xpos*PlaneNormal[0] + ypos*PlaneNormal[1] + (zpos-shift_z_box)*PlaneNormal[2];
				double dist = -FullPlaneEq[f]-buffer_f;
				if( dist < mean_z_pos ){
					double cur_sp = 0.;
					for(unsigned int j=0;j<nbAt_OriSys;j++){
						for(unsigned int k=0;k<nb_shift;k++){
							current_shift[0] = swap_shift[0]*k + OriSys_f1->getAtomRef(j).pos.x;
							current_shift[1] = swap_shift[1]*k + OriSys_f1->getAtomRef(j).pos.y;
							current_shift[2] = swap_shift[2]*k + OriSys_f1->getAtomRef(j).pos.z;
							cur_sp += MT->gaussian(xpos, ypos, zpos, current_shift[0], current_shift[1], current_shift[2], sigma);
						}
					}
					if( cur_sp < tol_rmatom ){
						index2rm.push_back(i);
					}
				}
			}
		}
		for(unsigned int i=0;i<index2rm.size();i++) opt_index2rm[s_z_i].push_back(index2rm[i]);
	} // end loop on f
	//cout << "Number of ions removed : " << index2rm.size() << endl;
	delete[] invTiltTrans;
	delete[] invRotMat;
	delete[] zaxis;
	delete[] mean_pos_orisys;
	delete[] mean_pos_subsys;
	delete[] diff_orisubsys;
	delete[] at_sp;
	delete[] current_shift;
	delete[] swap_shift;
	} // end loop on s_z
	
		unsigned int opt_ind = 0;
		for(unsigned int i=1;i<nbSample_S_z;i++)
			if( opt_index2rm[i].size() < opt_index2rm[opt_ind].size() ) opt_ind = i;
		cout << "Final number of ions removed : " << opt_index2rm[opt_ind].size() << endl;
		RemoveAtoms(opt_index2rm[opt_ind], true);

}	
	
double AtomicSystem::MakeSurfaceNeutral(bool verbose){
	if( IsCharge == false ) return 0.;
	bool Possible = false; // at least we should have two atom types having charge with opposite sign
	bool pos(false),neg(false);
	for(unsigned int n=0;n<nbAtomType;n++){
		if( AtomCharge[n] > 0 ) pos = true;
		if( AtomCharge[n] < 0 ) neg = true;
		if( pos && neg ){
			Possible = true;
			break;
		}
	}
	if( !Possible ){
		cout << "It is not possible to construct charge neutral surfaces as there are not two types of ions having charge with opposite signs" << endl;
		return 0.;
	}
	// Compute the average distance between ions to compute the standard deviation of the charge Gaussian expansion
	double ave_dist = ComputeAverageDistance();
	// Search how many ions are in the bottom slab to compute the maximum number of ions to move
	unsigned int nbMoveMax = 0;
	unsigned int MinNbMoveMax = 20;
	double slab_width = ave_dist*2.;
	for(unsigned int i=0;i<nbAtom;i++) if( AtomList[i].pos.z < slab_width ) nbMoveMax++;
	if( nbMoveMax < MinNbMoveMax ) nbMoveMax = MinNbMoveMax;
	if( verbose ) cout << "Trying to make charge neutral surfaces by moving " << nbMoveMax << " ions (or group of ions to be not separed) from the bottom surface to the upper one.." << endl;
	// Increase H3 length to be sure to have free surfaces
	double fac_vac = 1.5;
	double shift_z = this->H3[2] * ((fac_vac-1.)/2.);
	double *InitH3 = new double[3];
	for(unsigned int i=0;i<3;i++){
		InitH3[i] = this->H3[i];
		this->H3[i] *= fac_vac;
	}
	computeInverseCellVec();
	double zero(0.);
	ApplyShift(zero,zero,shift_z);
	computeWrap();

	double *true_z = new double[nbAtom];
	for(unsigned int i=0;i<nbAtom;i++) true_z[i] = AtomList[i].pos.z;
	double zmid = ( MT->max_p(true_z,nbAtom) + MT->min_p(true_z,nbAtom) ) / 2.;
	unsigned int nbPts_g = floor(H3[2]*2./ave_dist);
	if( nbPts_g % 2 != 1 ) nbPts_g++;
	unsigned int dens_ind = Compute1dDensity("Charge","z",ave_dist,nbPts_g);
	if( verbose ) Print1dDensity("InitialChargeDensity.dat","Charge");
	unsigned int mid = floor((double) (nbPts_g) * zmid / H3[2])+1;
	double low_int(0.),up_int(0.); // integrals of charge density in the upper and lower part of the system
	for(unsigned int i=0;i<mid;i++) low_int += density_prof[dens_ind][i*2]; // no need to account for distance as step is regular and system is divided by the middle
	for(unsigned int i=mid+1;i<nbPts_g;i++) up_int += density_prof[dens_ind][i*2];
	vector<unsigned int> ind_moved;
	ind_moved.push_back(0);
	double charge_bal = fabs(low_int - up_int); // try to minize this quantity by moving ions of H3 from low to high
	double init_charge_bal = charge_bal;
	// at each step move the ions which is the most at the bottom
	double *z = new double[nbAtom];
	if( !IsNotSepTag ) for(unsigned int i=0;i<nbAtom;i++) z[i] = AtomList[i].pos.z;
	else{
		ComputeNotSepList();
		for(unsigned int i=0;i<nbAtom;i++){
			true_z[i] = AtomList[i].pos.z;
			if( NotSepTag[i][0] >= 0 ) z[i] = AtomList[i].pos.z; // store only ions which are not stored because of the NotSep list
			else z[i] = std::numeric_limits<double>::max();
		}

	}
	unsigned int ind_tomove = MT->min_p_ind(z,nbAtom);
	double init_zmin = z[ind_tomove];
	double final_zmin = init_zmin;
	unsigned int opt_ind(0);
	vector<unsigned int> nbIonsMoved;
	nbIonsMoved.push_back(0);
	for(unsigned int i=0;i<nbMoveMax;i++){
		if( verbose ) Dis->ProgressBar(nbMoveMax,i);
		AtomList[ind_tomove].pos.x += InitH3[0];
		AtomList[ind_tomove].pos.y += InitH3[1];
		AtomList[ind_tomove].pos.z += InitH3[2];
		nbIonsMoved.push_back(1);
		if( IsNotSepTag){
			for(unsigned int j=0;j<NotSepTag[ind_tomove][0];j++){
				AtomList[NotSepTag[ind_tomove][j+1]].pos.x += InitH3[0];
				AtomList[NotSepTag[ind_tomove][j+1]].pos.y += InitH3[1];
				AtomList[NotSepTag[ind_tomove][j+1]].pos.z += InitH3[2];
				nbIonsMoved[nbIonsMoved.size()-1]++;
			}
		}
		computeWrap();
		true_z[ind_tomove] += InitH3[2];
		zmid = ( MT->max_p(true_z,nbAtom) + MT->min_p(true_z,nbAtom) ) / 2.;
		mid = floor((double) (nbPts_g) * zmid / H3[2])+1;
		// recompute density
		dens_ind = Compute1dDensity("Charge","z",ave_dist,nbPts_g);
		low_int = 0.;
		up_int = 0.;
		for(unsigned int i=0;i<mid;i++) low_int += ( density_prof[dens_ind][(i+1)*2] + density_prof[dens_ind][i*2] );
		for(unsigned int i=mid;i<nbPts_g-1;i++) up_int += ( density_prof[dens_ind][(i+1)*2] + density_prof[dens_ind][i*2] );
		// store value of charge balance and ion index
		ind_moved.push_back(ind_tomove);
		if( charge_bal > fabs(low_int - up_int) ){
			charge_bal = fabs(low_int - up_int);
			final_zmin = z[ind_tomove];
			opt_ind = i+1;
		}
		z[ind_tomove] = AtomList[ind_tomove].pos.z;
		ind_tomove = MT->min_p_ind(z,nbAtom);
	}
	for(unsigned int i=ind_moved.size()-1;i>opt_ind;i--){
		AtomList[ind_moved[i]].pos.x -= InitH3[0];
		AtomList[ind_moved[i]].pos.y -= InitH3[1];
		AtomList[ind_moved[i]].pos.z -= InitH3[2];
		if( IsNotSepTag){
			for(unsigned int j=0;j<NotSepTag[ind_moved[i]][0];j++){
				AtomList[NotSepTag[ind_moved[i]][j+1]].pos.x -= InitH3[0];
				AtomList[NotSepTag[ind_moved[i]][j+1]].pos.y -= InitH3[1];
				AtomList[NotSepTag[ind_moved[i]][j+1]].pos.z -= InitH3[2];
			}
		}
	}
	// recompute density
	computeWrap();
	dens_ind = Compute1dDensity("Charge","z",ave_dist,nbPts_g);
	// Reset the system at its initial position and size
	ApplyShift(zero,zero,-shift_z-(final_zmin-init_zmin));
	//ApplyShift(zero,zero,-shift_z);
	for(unsigned int i=0;i<3;i++) this->H3[i] = InitH3[i];
	computeInverseCellVec();
	delete[] InitH3;
	delete[] z;
	delete[] true_z;
	if( verbose ){
		cout << "Initial charge balance = " << init_charge_bal << ", final one = " << charge_bal << ", number of ions moved = " << nbIonsMoved[opt_ind] << endl;
		Print1dDensity("FinalChargeDensity.dat","Charge");
	}
	return final_zmin-init_zmin;
}

double AtomicSystem::ComputeAverageDistance(){
	double rcut = 6.; // maybe try latter to have an more objective value for rcut
	searchNeighbours(rcut);
	double ave_dist = 0.;
	double xp, yp ,zp, xpos, ypos, zpos;
	for(unsigned int i=0;i<nbAtom;i++){
		xpos = WrappedPos[i].x;
		ypos = WrappedPos[i].y;
		zpos = WrappedPos[i].z;
		vector<double> dist;
		for(unsigned int j=0;j<Neighbours[i*(nbMaxN+1)];j++){
			unsigned int id = Neighbours[i*(nbMaxN+1)+j+1];
			long unsigned int id1 = i*nbMaxN*3+j*3;
			long unsigned int id2 = id1+1;
			long unsigned int id3 = id1+2;
			xp = WrappedPos[id].x+CLNeighbours[id1]*H1[0]+CLNeighbours[id2]*H2[0]+CLNeighbours[id3]*H3[0]-xpos;
			yp = WrappedPos[id].y+CLNeighbours[id1]*H1[1]+CLNeighbours[id2]*H2[1]+CLNeighbours[id3]*H3[1]-ypos;
			zp = WrappedPos[id].z+CLNeighbours[id1]*H1[2]+CLNeighbours[id2]*H2[2]+CLNeighbours[id3]*H3[2]-zpos;
			dist.push_back(xp*xp+yp*yp+zp*zp);
		}
		ave_dist += sqrt(MT->min_vec(dist));
	}
	return ave_dist/nbAtom;
}

AtomicSystem::~AtomicSystem(){
	if( AtomList && this->IsAtomListMine ){
		delete[] AtomList;
	}
	if( this->IsCellVecMine ){
		if( H1 ) delete[] H1;
		if( H2 ) delete[] H2;
		if( H3 ) delete[] H3;
	}
	if( G1 ) delete[] G1;
	if( G2 ) delete[] G2;
	if( G3 ) delete[] G3;
	delete MT;
	if( IsWrappedPos ) delete[] WrappedPos;
	if( IsNeighbours ){
		delete[] Neighbours;
		delete[] CLNeighbours;
	}
	if( IsSetAux ) for(unsigned int i=0;i<Aux.size();i++) delete[] Aux[i];
	if( IsCrystalMine ) delete _MyCrystal;
	if( AreTypeMassChargeMine ){
		delete[] AtomType;
		delete[] AtomMass;
		delete[] AtomCharge;
	}
	for(unsigned int i=0;i<density_prof.size();i++){
		delete[] density_prof[i];
		delete[] density_name[i];
	}
	if( this->IsNotSepTag ) delete[] NotSepTag;
	if( IsMolId && IsMolIdMine ) delete[] MolId;
	if( IsBond && IsBondMine ){
		delete[] Bonds;
		delete[] BondType;
	}
	if( IsAngle && IsAngleMine ){
		delete[] Angles;
		delete[] AngleType;
	}
	if( IsPeriodicArr ) delete[] PeriodicArr;
	delete Dis;
}


