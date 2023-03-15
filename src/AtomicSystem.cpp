#include "AtomicSystem.h"
#include "AtomHicConfig.h"
#include "Crystal.h"
#include <cmath>
#include <sstream>
#include <vector>
#include <iomanip>

using namespace std;

AtomicSystem::AtomicSystem(Crystal *_MyCrystal, double xhi, double yhi, double zhi, vector<int> cl_box):_MyCrystal(_MyCrystal){
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
	this->IsTilted = false;
	this->G1 = new double[3];
	this->G2 = new double[3];
	this->G3 = new double[3];
	this->IsG = true;
	computeInverseCellVec();
	// compute the atomic positions by testing linear combination of unit cell crystal
	double cla1, cla2, cla3;
       	int arr[3];
	if( fabs(_MyCrystal->getA1()[0]) >1e-3 ) arr[0] = round(xhi/fabs(_MyCrystal->getA1()[0]));
	else arr[0] = 0;
	if( fabs(_MyCrystal->getA1()[1]) >1e-3 ) arr[1] = round(yhi/fabs(_MyCrystal->getA1()[1]));
	else arr[1] = 0;
	if( fabs(_MyCrystal->getA1()[2]) >1e-3 ) arr[2] = round(zhi/fabs(_MyCrystal->getA1()[2]));
	else arr[2] = 0;
	cla1 = MT->max(arr,3)+1;
	if( fabs(_MyCrystal->getA2()[0]) >1e-3 ) arr[0] = round(xhi/fabs(_MyCrystal->getA2()[0]));
	else arr[0] = 0;
	if( fabs(_MyCrystal->getA2()[1]) >1e-3 ) arr[1] = round(yhi/fabs(_MyCrystal->getA2()[1]));
	else arr[1] = 0;
	if( fabs(_MyCrystal->getA2()[2]) >1e-3 ) arr[2] = round(zhi/fabs(_MyCrystal->getA2()[2]));
	else arr[2] = 0;
	cla2 = MT->max(arr,3)+1;
	if( fabs(_MyCrystal->getA3()[0]) >1e-3 ) arr[0] = round(xhi/fabs(_MyCrystal->getA3()[0]));
	else arr[0] = 0;
	if( fabs(_MyCrystal->getA3()[1]) >1e-3 ) arr[1] = round(yhi/fabs(_MyCrystal->getA3()[1]));
	else arr[1] = 0;
	if( fabs(_MyCrystal->getA3()[2]) >1e-3 ) arr[2] = round(zhi/fabs(_MyCrystal->getA3()[2]));
	else arr[2] = 0;
	cla3 = MT->max(arr,3)+1;
	double delta_x, delta_y, delta_z, xpos, ypos, zpos;
	unsigned int count;
	double tolIonPos = 1e-9;
	bool break_comp = false;
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
		unsigned int *cl_test = new unsigned int[3];
		double tol = 1e-2;
		// initialize the NotSepTag array
		this->NotSepTag = new vector<int>[this->nbAtom];
		this->IsNotSepTag = true;
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
								id_notsep = count;
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
									this->NotSepTag[count].push_back(-1);
									this->NotSepTag[id_notsep].push_back(count);
									this->NotSepTag[id_notsep][0] += 1;
									count += 1;
								}
							}
						}
						if( break_comp) break;
					}
					if( break_comp) break;
				}
				if( break_comp) break;
			}
			if( break_comp) break;
		}
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
}
// construct AtomicSystem giving AtomList, number of atom and cell vectors 
AtomicSystem::AtomicSystem(Atom *AtomList, unsigned int nbAtom, Crystal *_MyCrystal, double *H1, double *H2, double *H3):AtomList(AtomList), nbAtom(nbAtom), _MyCrystal(_MyCrystal), H1(H1), H2(H2), H3(H3){
	read_params_atsys();
	this->IsAtomListMine = false;
	this->IsCellVecMine = false;
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
}
// Constructor using the filename of an atomic file (.cfg, .lmp, .xsf,..)
AtomicSystem::AtomicSystem(const string& filename)
{
	read_params_atsys();
	this->AtomType = new string[this->MaxAtomType];
	this->AtomMass = new double[this->MaxAtomType];
	this->AtomCharge = new double[this->MaxAtomType];
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
	string ext=filename.substr(filename.find_last_of(".") + 1);
	if( ext == "lmp" ){
		this->read_lmp_file(filename);
	}else if( ext == "cfg" || ext == "xsf" ){
		this->read_cfg_file(filename);
	}else{
		cerr << "The \"" << ext << "\" is unknown (known extension : lmp, cfg(or xsf)), aborting the construction of atomic system" << endl;
	}
	this->MT = new MathTools;
	this->G1 = new double[3];
	this->G2 = new double[3];
	this->G3 = new double[3];
	this->IsG = true;
	computeInverseCellVec();
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

unsigned int AtomicSystem::Compute1dDensity(std::string auxname, std::string dir, double sigma, unsigned int nbPts){
	if( !IsWrappedPos ) computeWrap();
	density_prof.push_back(new double[nbPts*2]);
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
					this->density_prof[this->density_prof.size()-1][i*2] = 0;
					this->density_prof[this->density_prof.size()-1][i*2+1] = this->H1[0]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[this->density_prof.size()-1][i*2] += MT->gaussian(this->density_prof[this->density_prof.size()-1][i*2+1], this->WrappedPos[j].x, sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
				}
			}else if( dir == "y" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[this->density_prof.size()-1][i*2] = 0;
					this->density_prof[this->density_prof.size()-1][i*2+1] = this->H2[1]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[this->density_prof.size()-1][i*2] += MT->gaussian(this->density_prof[this->density_prof.size()-1][i*2+1], this->WrappedPos[j].y, sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
				}
			}else if( dir == "z" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[this->density_prof.size()-1][i*2] = 0;
					this->density_prof[this->density_prof.size()-1][i*2+1] = this->H3[2]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[this->density_prof.size()-1][i*2] += MT->gaussian(this->density_prof[this->density_prof.size()-1][i*2+1], this->WrappedPos[j].z, sigma)*(this->AtomMass[this->AtomList[j].type_uint-1]);
				}
			}else{
			        cout << "The provided direction \"" << dir << "\" for density computation has not been recognize" << endl;
				return 0;
			}
		}else if( auxname == "Charge" ){
			if( dir == "x" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[this->density_prof.size()-1][i*2] = 0;
					this->density_prof[this->density_prof.size()-1][i*2+1] = this->H1[0]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[this->density_prof.size()-1][i*2] += MT->gaussian(this->density_prof[this->density_prof.size()-1][i*2+1], this->WrappedPos[j].x, sigma)*(this->AtomCharge[this->AtomList[j].type_uint-1]);
				}
			}else if( dir == "y" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[this->density_prof.size()-1][i*2] = 0;
					this->density_prof[this->density_prof.size()-1][i*2+1] = this->H2[1]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[this->density_prof.size()-1][i*2] += MT->gaussian(this->density_prof[this->density_prof.size()-1][i*2+1], this->WrappedPos[j].y, sigma)*(this->AtomCharge[this->AtomList[j].type_uint-1]);
				}
			}else if( dir == "z" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[this->density_prof.size()-1][i*2] = 0;
					this->density_prof[this->density_prof.size()-1][i*2+1] = this->H3[2]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[this->density_prof.size()-1][i*2] += MT->gaussian(this->density_prof[this->density_prof.size()-1][i*2+1], this->WrappedPos[j].z, sigma)*(this->AtomCharge[this->AtomList[j].type_uint-1]);
				}
			}else{
			        cout << "The provided direction \"" << dir << "\" for density computation has not been recognize" << endl;
				return 0;
			}
		}else if ( auxname == "Atomic" ){
			if( dir == "x" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[this->density_prof.size()-1][i*2] = 0;
					this->density_prof[this->density_prof.size()-1][i*2+1] = this->H1[0]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[this->density_prof.size()-1][i*2] += MT->gaussian(this->density_prof[this->density_prof.size()-1][i*2+1], this->WrappedPos[j].x, sigma);
				}
			}else if( dir == "y" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[this->density_prof.size()-1][i*2] = 0;
					this->density_prof[this->density_prof.size()-1][i*2+1] = this->H2[1]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[this->density_prof.size()-1][i*2] += MT->gaussian(this->density_prof[this->density_prof.size()-1][i*2+1], this->WrappedPos[j].y, sigma);
				}
			}else if( dir == "z" ){
				for(unsigned int i=0;i<nbPts;i++){
					this->density_prof[this->density_prof.size()-1][i*2] = 0;
					this->density_prof[this->density_prof.size()-1][i*2+1] = this->H3[2]*i/(nbPts-1.);
					for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[this->density_prof.size()-1][i*2] += MT->gaussian(this->density_prof[this->density_prof.size()-1][i*2+1], this->WrappedPos[j].z, sigma);
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
				this->density_prof[this->density_prof.size()-1][i*2] = 0;
				this->density_prof[this->density_prof.size()-1][i*2+1] = this->H1[0]*i/(nbPts-1.);
				for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[this->density_prof.size()-1][i*2] += MT->gaussian(this->density_prof[this->density_prof.size()-1][i*2+1], this->WrappedPos[j].x, sigma)*(this->Aux[indexaux][j]);
			}
		}else if( dir == "y" ){
			for(unsigned int i=0;i<nbPts;i++){
				this->density_prof[this->density_prof.size()-1][i*2] = 0;
				this->density_prof[this->density_prof.size()-1][i*2+1] = this->H2[1]*i/(nbPts-1.);
				for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[this->density_prof.size()-1][i*2] += MT->gaussian(this->density_prof[this->density_prof.size()-1][i*2+1], this->WrappedPos[j].y, sigma)*(this->Aux[indexaux][j]);
			}
		}else if( dir == "z" ){
			for(unsigned int i=0;i<nbPts;i++){
				this->density_prof[this->density_prof.size()-1][i*2] = 0;
				this->density_prof[this->density_prof.size()-1][i*2+1] = this->H3[2]*i/(nbPts-1.);
				for(unsigned int j=0;j<this->nbAtom;j++) this->density_prof[this->density_prof.size()-1][i*2] += MT->gaussian(this->density_prof[this->density_prof.size()-1][i*2+1], this->WrappedPos[j].z, sigma)*(this->Aux[indexaux][j]);
			}
		}else{
		        cout << "The provided direction \"" << dir << "\" for density computation has not been recognize (directions available : x y z)" << endl;
			return 0;
		}
	}
	this->density_name.push_back(new string[2]);
	this->density_name[this->density_name.size()-1][0] = auxname;
	this->density_name[this->density_name.size()-1][1] = dir;
	this->density_nbPts.push_back(nbPts);
	return this->density_name.size()-1;
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
}

// TODO : add verification that AtomType_uint of _MyCrystal and correspond to the same atom type
void AtomicSystem::setCrystal(Crystal* MyCrystal){
	this->_MyCrystal = MyCrystal;
	this->IsCrystalDefined = true;
}

// TODO : add verification that AtomType_uint of _MyCrystal and correspond to the same atom type
void AtomicSystem::setCrystal(const std::string& CrystalName){
	cout << "Setting crystal : " << CrystalName << endl;
	this->_MyCrystal = new Crystal(CrystalName);
	this->IsCrystalDefined = true;
	this->IsCrystalMine = true;
}

void AtomicSystem::computeWrap(){
	WrappedPos = new Position[this->nbAtom];
	IsWrappedPos = true;
	double x,y,z;
	for(unsigned int i=0;i<this->nbAtom;i++){
		// compute reduced coordinates
		x = this->AtomList[i].pos.x*G1[0]+this->AtomList[i].pos.y*G2[0]+this->AtomList[i].pos.z*G3[0];
		y = this->AtomList[i].pos.x*G1[1]+this->AtomList[i].pos.y*G2[1]+this->AtomList[i].pos.z*G3[1];
		z = this->AtomList[i].pos.x*G1[2]+this->AtomList[i].pos.y*G2[2]+this->AtomList[i].pos.z*G3[2];
		if( x > 1. || x < 0. ) x = x-floor(x);
		if( y > 1. || y < 0. ) y = y-floor(y);
		if( z > 1. || z < 0. ) z = z-floor(z);
		// cartesian coordinates
		this->WrappedPos[i].x = x*H1[0]+y*H2[0]+z*H3[0];
		this->WrappedPos[i].y = x*H1[1]+y*H2[1]+z*H3[1];
		this->WrappedPos[i].z = x*H1[2]+y*H2[2]+z*H3[2];
	}
}

// Searching neighbours using cell list algorithm
void AtomicSystem::searchNeighbours(const double& rc){
	computeWrap();
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
	vector<vector<unsigned int>> Cells;
	if( this->IsTilted ){
		// compute plane equations (as everywhere the only tilts considered are xy, yz, xz)
		double H1H3_z;
	        if(fabs(this->H3[2]) > 1e-6) H1H3_z = -(this->H3[1]/this->H3[2]);
		else H1H3_z = 0.;
		double H2H3_y = -(this->H2[0]/this->H2[1]); 
		double H2H3_z = ((this->H2[0]*this->H3[1]/this->H2[1])-this->H3[0])/this->H3[2];
		double planeH1H3, planeH2H3;
		for(unsigned int i=0;i<nbCellX;i++){
        	        for(unsigned int j=0;j<nbCellY;j++){
        	                for(unsigned int k=0;k<nbCellZ;k++){
        	                        Cells.push_back(vector<unsigned int>());
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
		for(unsigned int i=0;i<nbCellX;i++){
        	        for(unsigned int j=0;j<nbCellY;j++){
        	                for(unsigned int k=0;k<nbCellZ;k++){
        	                        Cells.push_back(vector<unsigned int>());
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
	for(unsigned int i=1;i<nbCellX*nbCellY*nbCellZ;i++){
		if( Cells[i].size() > this->nbMaxN ) this->nbMaxN = Cells[i].size();
	}
	this->nbMaxN *= (int) (1.5*4.*M_PI*pow(rc,3.)/(3.*CellSizeX*CellSizeY*CellSizeZ)); // 1.5 is a security factor
	if( this->IsNeighbours ){
		delete[] this->Neighbours;
		delete[] this->CLNeighbours;
	} // TODO maybe issue here if we delete the var we may need to reclare them ?
	this->Neighbours = new unsigned int[(this->nbMaxN+1)*this->nbAtom];
	this->CLNeighbours = new int[(this->nbMaxN*3)*this->nbAtom]; // contain the periodic condition (Nclx, Ncly, Nclz) applied for atom to be a neighbour

	// Perform neighbour research
	double xpos,ypos,zpos;
	int ibx, jby, kbz;
	int Nclx, Ncly, Nclz;
	double prog=0.;
	unsigned int countN = 0;
	unsigned int currentId, currentId2;
	cout << "Performing neighbour research" << endl;
	cout << "\r[" << string(bar_length*prog,'X') << string(bar_length*(1-prog),'-') << "] " << setprecision(3) << 100*prog << "%";
	for(unsigned int i=0;i<nbCellX;i++){
		for(unsigned int j=0;j<nbCellY;j++){
			for(unsigned int k=0;k<nbCellZ;k++){
				prog = double(i*nbCellZ*nbCellY+j*nbCellZ+k)/double(nbCellX*nbCellY*nbCellZ);
				cout << "\r[" << string(floor(bar_length*prog),'X') << string(ceil(bar_length*(1-prog)),'-') << "] " << setprecision(3) << 100*prog << "%";
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
										this->Neighbours[currentId*(this->nbMaxN+1)] += 1; // add this neighbours to the neighbour count
										this->Neighbours[currentId*(this->nbMaxN+1)+countN+1] = currentId2; // add this neighbours to the neighbour list 
										this->CLNeighbours[currentId*this->nbMaxN*3+countN*3] = Nclx; // store the cl used for this neighbour
										this->CLNeighbours[currentId*this->nbMaxN*3+countN*3+1] = Ncly; // store the cl used for this neighbour
										this->CLNeighbours[currentId*this->nbMaxN*3+countN*3+2] = Nclz; // store the cl used for this neighbour
										countN += 1;
									}
								}
							}
						}
						if( countN > this->nbMaxN ) cout << "Warning the number of found neighbour exceed the maximum number of neighbour allowed" << endl;
					}
				}
			}
		}
	}
	cout << endl;
	IsNeighbours = true;
}

void AtomicSystem::setAux(const double* aux, const string& AuxName){
	IsSetAux = true;
	this->Aux.push_back(new double[this->nbAtom]);
	this->Aux_name.push_back(AuxName);
	this->Aux_size.push_back(1);
	for(unsigned int i=0;i<this->nbAtom;i++) Aux[Aux.size()-1][i] = aux[i];
}

void AtomicSystem::setAux_vec(const double* aux, const unsigned int size, const string& AuxName){
	IsSetAux = true;
	this->Aux.push_back(new double[this->nbAtom*size]);
	this->Aux_name.push_back(AuxName);
	this->Aux_size.push_back(size);
	for(unsigned int i=0;i<this->nbAtom;i++){
		for(unsigned int j=0;j<size;j++){
			Aux[Aux.size()-1][i*size+j] = aux[i*size+j];
		}
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

void AtomicSystem::read_lmp_file(const string& filename){
	cout << "Reading " << filename << " file..";
	ifstream file(filename, ios::in);
	size_t pos_at, pos_x, pos_y, pos_z, pos_tilt, pos_attype, pos_Mass, pos_At;
	unsigned int line_Mass(1000), line_At(1000), buffer_uint, buffer_uint_1, count(0);
	double buffer_1, buffer_2, buffer_3, buffer_4;
	double xlo,xhi,ylo,yhi,zlo,zhi;
	string buffer_s, buffer_s_1, buffer_s_2, line;
	if(file){
		while(file){
			getline(file,line);

			// find number of atom
			pos_at=line.find("atoms");
			if(pos_at!=string::npos){
				istringstream text(line);
				text >> buffer_uint;
				this->nbAtom = buffer_uint;
				AtomList = new Atom[this->nbAtom];
			}

			// find H1 vector
			pos_x=line.find("xlo xhi");
			if(pos_x!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				xhi = buffer_2;
				xlo = buffer_1;
			}

			// find H2 vector
			pos_y=line.find("ylo yhi");
			if(pos_y!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				yhi = buffer_2;
				ylo = buffer_1;
			}

			// find H3 vector
			pos_z=line.find("zlo zhi");
			if(pos_z!=string::npos){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				zhi = buffer_2;
				zlo = buffer_1;
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
				this->nbAtomType = buffer_1;
			}

			// get lines where are the keywords Masses and Atoms to get atom type masses and positions
			pos_Mass=line.find("Masses");
			if(pos_Mass!=string::npos) line_Mass = count;
			if( count > line_Mass+1 && count < line_Mass+2+this->nbAtomType ){
				istringstream text(line);
				text >> buffer_uint >> buffer_1 >> buffer_s_1 >> buffer_s;
				this->AtomMass[buffer_uint-1] = buffer_1;
				this->AtomType[buffer_uint-1] = buffer_s;
			}
			pos_At=line.find("Atoms #");
			if(pos_At!=string::npos){
				istringstream text(line);
				text >> buffer_s_1 >> buffer_s_2 >> buffer_s;
				if( buffer_s == "charge" ) this->IsCharge = true;
			       	line_At = count;
			}
			if( count > line_At+1 ){
				istringstream text(line);
				if( this->IsCharge ){
					text >> buffer_uint >> buffer_uint_1 >> buffer_1 >> buffer_2 >> buffer_3 >> buffer_4;
					this->AtomList[buffer_uint-1].pos.x = buffer_2;
					this->AtomList[buffer_uint-1].pos.y = buffer_3;
					this->AtomList[buffer_uint-1].pos.z = buffer_4;
					this->AtomList[buffer_uint-1].type_uint = buffer_uint_1;
					this->AtomCharge[buffer_uint_1-1] = buffer_1;
				}else{
					text >> buffer_uint >> buffer_uint_1 >> buffer_2 >> buffer_3 >> buffer_4;
					this->AtomList[buffer_uint-1].pos.x = buffer_2;
					this->AtomList[buffer_uint-1].pos.y = buffer_3;
					this->AtomList[buffer_uint-1].pos.z = buffer_4;
					this->AtomList[buffer_uint-1].type_uint = buffer_uint_1;
					this->AtomCharge[buffer_uint_1-1] = 0.;
				}
			}
			count += 1;
		}
		// compute the cell vectors
		// TEST
		//double arr[4] = {0.,this->H2[0],this->H3[0],this->H2[0]+this->H3[0]};
		//this->H1[0] = xhi-xlo+this->MT->min(arr,4)-this->MT->max(arr,4);
		//double arr_2[2] = {0.,this->H3[1]};
		//this->H2[1] = yhi-ylo+this->MT->min(arr_2,2)-this->MT->max(arr_2,2);
		//this->H3[2] = zhi-zlo; // TODO verify if this is good for tilted box
		// END TEST
		this->H1[0] = xhi-xlo;
		this->H2[1] = yhi-ylo;
		this->H3[2] = zhi-zlo; // TODO verify if this is good for tilted box
		file.close();
	}else{
		cout << "The file " << filename << " cannot be openned" << endl;
	}
	cout << " done !" << endl;
}

void AtomicSystem::read_cfg_file(const string& filename){
	cout << "Reading " << filename << " file..";
	ifstream file(filename, ios::in);
	if(file){
		unsigned int line_dt(1000), line_At(1000), line_H_tilt(1000), line_H(1000), line_at(1000), buffer_uint, buffer_uint_1, count_H(0), count(0), nbAux(0), aux_count;
		size_t pos_dt, pos_At, pos_H_tilt, pos_H, pos_charge, pos_at;
		double buffer_1, buffer_2, buffer_3, buffer_4;
		double xlo,xhi,ylo,yhi,zlo,zhi;
		string buffer_s, buffer_s_1, buffer_s_2, line;
		while(file){
			getline(file,line);

			// find timestep
			pos_dt=line.find("TIMESTEP");
			if( pos_dt!=string::npos ) line_dt = count;
			if( count == line_dt+1 ){
				istringstream text(line);
				text >> buffer_1;
			       	this->timestep = buffer_1;
			}

			// find number of atom
			pos_At=line.find("NUMBER OF ATOMS");
			if( pos_At!=string::npos ) line_At = count;
			if( count == line_At+1 ){
				istringstream text(line);
				text >> buffer_uint;
			       	this->nbAtom = buffer_uint;
				AtomList = new Atom[this->nbAtom];
			}

			// find box vectors
			pos_H_tilt=line.find("BOX BOUNDS xy xz yz pp pp pp");
			if( pos_H_tilt!=string::npos ){
				line_H_tilt = count;
				this->IsTilted = true;
			}
			pos_H=line.find("BOX BOUNDS pp pp pp");
			if( pos_H!=string::npos ) line_H = count;
			if( IsTilted && count > line_H_tilt && count < line_H_tilt+4 ){
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
				}
				count_H += 1;
			}else if( !IsTilted && count > line_H && count < line_H+4 ){
				istringstream text(line);
				text >> buffer_1 >> buffer_2;
				if( count_H == 0 ) this->H1[0] = buffer_2-buffer_1;
				else if( count_H == 1 ) this->H2[1] = buffer_2-buffer_1;
				else if( count_H == 2 ) this->H3[2] = buffer_2-buffer_1;
				count_H += 1;
			}

			// search if atom are charged
			pos_charge=line.find(" q ");
			if( pos_charge!=string::npos) this->IsCharge = true;

			// find and get atom positions
			pos_at=line.find("ITEM: ATOMS");
			if( pos_at!=string::npos ){
				istringstream text(line);
				while(text >> buffer_s){
					nbAux += 1;
					if( nbAux > 9 ){
					       this->IsSetAux = true;
					       this->Aux_name.push_back(buffer_s);
					       this->Aux.push_back(new double[this->nbAtom]);
					}	       
				}
				line_at = count;
				nbAux -= 9;
			}
			if( count > line_at ){
				istringstream text(line);
				text >> buffer_uint >> buffer_uint_1 >> buffer_s >> buffer_1 >> buffer_2 >> buffer_3 >> buffer_4;
				this->AtomList[buffer_uint-1].pos.x = buffer_1;
				this->AtomList[buffer_uint-1].pos.y = buffer_2;
				this->AtomList[buffer_uint-1].pos.z = buffer_3;
				this->AtomList[buffer_uint-1].type_uint = buffer_uint_1;
				this->AtomCharge[buffer_uint_1-1] = buffer_4;
				this->AtomType[buffer_uint_1-1] = buffer_s;
				if(nbAux>0){
					aux_count = 0;
					while(text >> buffer_4){
						Aux[aux_count][buffer_uint-1] = buffer_4;
						aux_count += 1;
					}
				}
			}
			count += 1;
		}
		file.close();
		if( IsTilted ){
			// compute the cell vectors
			// TEST
			double arr[4] = {0.,this->H2[0],this->H3[0],this->H2[0]+this->H3[0]};
			this->H1[0] = xhi-this->MT->min(arr,4)-this->MT->max(arr,4);
			double arr_2[2] = {0.,this->H3[1]};
			this->H2[1] = yhi-this->MT->min(arr_2,2)-this->MT->max(arr_2,2);
			// END TEST
		}
		// search the number of atom type
		for(unsigned int i=0;i<this->MaxAtomType;i++){
			if( this->AtomType[i] == "" ){
				this->nbAtomType = i;
				break;
			}
		}

		// read MASSES database to get the masses of ions
		vector<string> element;
		vector<double> masses;
		char *database_env = getenv("MASSES_DATABASE");
		string database;
		if (database_env) {
			database = database_env;
		} else {
			#ifdef MASSES_DATABASE
			database = MASSES_DATABASE;
			#endif
		}
		if( database.empty() ) cout << "Warning database environment for masses is empty" << endl;
		else{
			string dataname = database + "/Masses.txt";
			ifstream filedata(dataname, ios::in);
			if( filedata ){
				while( filedata ){
					filedata >> buffer_s >> buffer_1;
					if( !filedata ) break;
					element.push_back(buffer_s);
					masses.push_back(buffer_1);
				}
				filedata.close();
			}
		}

		bool element_find;
		for(unsigned int j=0;j<this->nbAtomType;j++){
			element_find= false;
			for(unsigned int i=0;i<element.size();i++){
				if( this->AtomType[j] == element[i] ){
					this->AtomMass[j] = masses[i];
					element_find = true;
					break;
				}
			}
			if( !element_find ) cout << "The mass for element " << this->AtomType[j] << " has not been found in the masses databse" << endl;
		}

	// end read cfg (xsf) file
	}else{
		cout << "The file " << filename << " cannot be openned" << endl;
	}
	cout << " done !" << endl;
}


// TODO : change format for having more prec
void AtomicSystem::print_lmp(const string& filename){
	ofstream writefile(filename);
	writefile << " # File generated using AtomHic\n";
	writefile << this->File_Heading;
        writefile << "\n\t" << this->nbAtom << "\tatoms\n\t" << this->nbAtomType << "\tatom types\n\n\t0.000000000000\t" << this->H1[0] << "\txlo xhi\n\t0.000000000000\t" << H2[1] << "\tylo yhi\n\t0.000000000000\t" << H3[2] << "\tzlo zhi\n";
	if( this->IsTilted ) writefile << "\t" << H2[0] << "\t" << H3[0] << "\t" << H3[1] << "\txy xz yz\n";
	writefile << "\nMasses\n\n";
	for(unsigned int i=0;i<this->nbAtomType;i++) writefile << "\t" << i+1 << "\t" << this->AtomMass[i] << "\t# " << this->AtomType[i] << "\n";
	if( IsCharge ){
		writefile << "\nAtoms # charge\n\n";
		for(unsigned int i=0;i<this->nbAtom;i++) writefile << i+1 << "\t" << this->AtomList[i].type_uint << "\t" << this->AtomCharge[this->AtomList[i].type_uint-1] << "\t" << this->AtomList[i].pos.x << "\t" << this->AtomList[i].pos.y << "\t" << this->AtomList[i].pos.z << "\n"; 
	}else{
		writefile << "\nAtoms # atomic\n\n";
		for(unsigned int i=0;i<this->nbAtom;i++) writefile << i+1 << "\t" << this->AtomList[i].type_uint << "\t" << this->AtomList[i].pos.x << "\t" << this->AtomList[i].pos.y << "\t" << this->AtomList[i].pos.z << "\n"; 
	}
	writefile.close();
}

void AtomicSystem::print_cfg(const string& filename){
	ofstream writefile(filename);
	writefile << "ITEM: TIMESTEP\n" << (int) this->timestep << "\nITEM: NUMBER OF ATOMS\n" << this->nbAtom << "\nITEM: ";
	if( IsTilted ){
		// compute the cell vectors
		double arr[4] = {0.,this->H2[0],this->H3[0],this->H2[0]+this->H3[0]};
		double arr_2[2] = {0.,this->H3[1]};
	       	writefile << "BOX BOUNDS xy xz yz pp pp pp\n" << this->MT->min(arr,4) << "\t" << this->H1[0]+this->MT->max(arr,4) << "\t" << H2[0] << "\n" << this->MT->min(arr_2,2) << "\t" << this->H2[1]+this->MT->max(arr_2,2) << "\t" << H3[0] << "\n0\t" << H3[2] << "\t" << H3[1] << "\n";
	}
	else writefile << "BOX BOUNDS pp pp pp\n" << "0\t" << H1[0] << "\n0\t" << H2[1] << "\n0\t" << H3[2] << "\n";
	writefile << "ITEM: ATOMS id type element xu yu zu q\n";
	for(unsigned int i=0;i<this->nbAtom;i++) writefile << i+1 << " " << this->AtomList[i].type_uint << " " << this->AtomType[this->AtomList[i].type_uint-1] << " " << this->AtomList[i].pos.x << " " << this->AtomList[i].pos.y << " " << this->AtomList[i].pos.z << " " << this->AtomCharge[this->AtomList[i].type_uint-1] << "\n";
	writefile.close();
}

void AtomicSystem::printSystem_aux(const string& filename, const string& AuxName){
	vector<string> AuxNames_v;
	string buffer;
	istringstream names(AuxName);
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
	
	ofstream writefile(filename);
	writefile << "ITEM: TIMESTEP\n" << (int) this->timestep << "\nITEM: NUMBER OF ATOMS\n" << this->nbAtom << "\nITEM: ";
	if( IsTilted ){
		// compute the cell vectors
		double arr[4] = {0.,this->H2[0],this->H3[0],this->H2[0]+this->H3[0]};
		double arr_2[2] = {0.,this->H3[1]};
	       	writefile << "BOX BOUNDS xy xz yz pp pp pp\n" << this->MT->min(arr,4) << "\t" << this->H1[0]+this->MT->max(arr,4) << "\t" << H2[0] << "\n" << this->MT->min(arr_2,2) << "\t" << this->H2[1]+this->MT->max(arr_2,2) << "\t" << H3[0] << "\n0\t" << H3[2] << "\t" << H3[1] << "\n";
	}
	else writefile << "BOX BOUNDS pp pp pp\n" << "0\t" << H1[0] << "\n0\t" << H2[1] << "\n0\t" << H3[2] << "\n";
	writefile << "ITEM: ATOMS id type element xu yu zu q";
        for(unsigned int i=0;i<AuxId.size();i++){
		if( Aux_size[i] == 1 ) writefile << " " << this->Aux_name[AuxId[i]];
		else for(unsigned int j=0;j<Aux_size[i];j++) writefile << " " << this->Aux_name[AuxId[i]] << "_" << j+1;
	}
	writefile << "\n";
	for(unsigned int i=0;i<this->nbAtom;i++){
		writefile << i+1 << " " << this->AtomList[i].type_uint << " " << this->AtomType[this->AtomList[i].type_uint-1] << " " << this->AtomList[i].pos.x << " " << this->AtomList[i].pos.y << " " << this->AtomList[i].pos.z << " " << this->AtomCharge[this->AtomList[i].type_uint-1];
        	for(unsigned int j=0;j<AuxId.size();j++){
			if( Aux_size[j] == 1 ) writefile << " " << this->Aux[AuxId[j]][i];
			else for(unsigned int k=0;k<Aux_size[j];k++) writefile << " " << this->Aux[AuxId[j]][i*Aux_size[j]+k];
		}
		writefile << "\n";
	}
	writefile.close();
}

void AtomicSystem::read_params_atsys(){
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+FixedParam_Filename;
	ifstream file(filename, ios::in);
	size_t pos_rcut, pos_lsph;
	string buffer_s, line;
	if(file){
		while(file){
			getline(file,line);
			pos_rcut=line.find("RCUT_NEIGHBOURS ");
			if(pos_rcut!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->r_cut_n;
			}
			pos_lsph=line.find("L_SPH_ST ");
			if(pos_lsph!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->l_sph_st;
			}
		}
	}else{
		cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
		exit(EXIT_FAILURE);
	}
}

vector<unsigned int> AtomicSystem::selectAtomInBox(const double x_lo,const double x_hi,const double y_lo,const double y_hi,const double z_lo,const double z_hi){
	if( !IsWrappedPos ) computeWrap();
	vector<unsigned int> AtList;
	for(unsigned int i=0;i<this->nbAtom;i++){
		if( this->WrappedPos[i].x > x_lo && this->WrappedPos[i].x < x_hi && this->WrappedPos[i].y > y_lo && this->WrappedPos[i].y < y_hi && this->WrappedPos[i].z > z_lo && this->WrappedPos[i].z < z_hi ) AtList.push_back(i);
	}
	return AtList;
}

AtomicSystem::~AtomicSystem(){
	if( this->IsAtomListMine ){
		delete[] AtomList;
	}
	if( this->IsCellVecMine ){
		delete[] H1;
		delete[] H2;
		delete[] H3;
	}
	delete[] G1;
	delete[] G2;
	delete[] G3;
	delete MT;
	if( IsWrappedPos ) delete[] WrappedPos;
	if( IsNeighbours ){
		delete[] Neighbours;
		delete[] CLNeighbours;
	}
	if( IsSetAux ) for(unsigned int i=0;i<Aux.size();i++) delete[] Aux[i];
	if( IsCrystalMine ) delete _MyCrystal;
	if( !IsCrystalDefined ){
		delete[] AtomType;
		delete[] AtomMass;
		delete[] AtomCharge;
	}
	for(unsigned int i=0;i<density_prof.size();i++){
		delete[] density_prof[i];
		delete[] density_name[i];
	}
	if( this->IsNotSepTag ) delete[] NotSepTag;
}
