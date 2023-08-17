#include "AtomHicConfig.h"
#include <iostream>
#include <vector>
#include "Bicrystal.h"
#include <omp.h>

using namespace std;

Bicrystal::Bicrystal(const string& crystalName, int h_a, int k_a, int l_a, double theta):h_a(h_a), k_a(k_a), l_a(l_a), theta(theta){
	this->MT = new MathTools;
	read_params();
	setCrystal(crystalName);

}
	
// Constructor for bicrystal with facetted GB with given misorientation and GB plane and facet type (for the moment this is implemented for 2D facets and for a GB composed 2 different facet types, i.e. with one facet direction parallel to x direction)
Bicrystal::Bicrystal(const string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p, vector<int> FacetsType, unsigned int N_facet):h_a(h_a), k_a(k_a), l_a(l_a), theta(theta), h_p(h_p), k_p(k_p), l_p(l_p){ //TODO warning use this-> for the variable in the initialization list
	this->MT = new MathTools;
	read_params();
	setOrientedCrystals(crystalName, true);
	double *Dir1_G1 = new double[3];
	double *Dir2_G1 = new double[3];
	double *Dir1_G2 = new double[3];
	double *Dir2_G2 = new double[3];
	for(unsigned int i=0;i<3;i++){
		Dir1_G1[i] = FacetsType[0]*_MyCrystal->getA1()[i]-FacetsType[2]*_MyCrystal->getA2()[i]+FacetsType[1]*_MyCrystal->getA3()[i];
		Dir2_G1[i] = FacetsType[3]*_MyCrystal->getA1()[i]+FacetsType[5]*_MyCrystal->getA2()[i]-FacetsType[4]*_MyCrystal->getA3()[i];
	}
	if( Dir1_G1[2]*Dir2_G1[2] > 0 || Dir1_G1[1]*Dir2_G1[1] < 0 ){
		if( Dir1_G1[2] < 0 ){
			for(unsigned int i=0;i<3;i++){
				Dir1_G1[i] = FacetsType[0]*_MyCrystal->getA1()[i]+FacetsType[2]*_MyCrystal->getA2()[i]-FacetsType[1]*_MyCrystal->getA3()[i];
			}
		} else {
			for(unsigned int i=0;i<3;i++){
				Dir2_G1[i] = FacetsType[3]*_MyCrystal->getA1()[i]-FacetsType[5]*_MyCrystal->getA2()[i]+FacetsType[4]*_MyCrystal->getA3()[i];
			}
		}
		if( Dir1_G1[2]*Dir2_G1[2] > 0 || Dir1_G1[1]*Dir2_G1[1] < 0 ){
			cerr << "The given facets do not permit to construct the wanted GB, aborting calculation" << endl;
			exit(EXIT_FAILURE);
		}
	}
	if( Dir1_G1[1] < 0 ){
		for(unsigned int i=0;i<3;i++){
			Dir1_G1[i] *= -1;
			Dir2_G1[i] *= -1;
		}
	}
	// search LC for which the facet gives horyzontal GB
	unsigned int n1,n2;
	double tol = 1e-1;
	unsigned int nbMaxCL = 150;
	bool breaked = false;
	for(unsigned int i1=1;i1<nbMaxCL;i1++){
		for(unsigned int i2=1;i2<nbMaxCL;i2++){
			if( fabs(i1*Dir1_G1[2] + i2*Dir2_G1[2]) < tol ){
				n1 = i1;
				n2 = i2;
				breaked = true;
				break;
			}
		}
		if( breaked ) break;
	}
	if( !breaked ){
		cerr << "We don't find linear combination of facet for this GB, consider increasing tolerance or number of CL investigated, aborting calculation" << endl;
		exit(EXIT_FAILURE);
	}
	// correct the direction accounting for the small deformation applied to crystals during construction of grains
	MT->MatDotVec(this->_MyCrystal->getTiltTrans(), Dir1_G1, Dir1_G1);
	MT->MatDotVec(this->_MyCrystal->getTiltTrans(), Dir2_G1, Dir2_G1);
	MT->MatDotVec(this->_MyCrystal2->getTiltTrans(), Dir1_G1, Dir1_G2);
	MT->MatDotVec(this->_MyCrystal2->getTiltTrans(), Dir2_G1, Dir2_G2);
	n1 *= N_facet;
	n2 *= N_facet;
	double DeltaZ = fabs(n1*Dir1_G1[2]);
	// search the number of linear combination for the two system to have the same x y length
	double GBspace = 0.5; //TODO maybe define elsewhere (in a file in the main prog for example..)
	double MaxMisfit = 0.02;
	unsigned int MaxDup = 50;
	this->xl1 = this->_MyCrystal->getOrientedSystem()->getH1()[0];
	this->xl2 = this->_MyCrystal2->getOrientedSystem()->getH1()[0];
	this->yl1 = this->_MyCrystal->getOrientedSystem()->getH2()[1];
	this->yl2 = this->_MyCrystal2->getOrientedSystem()->getH2()[1];
	bool find = false;
	for(unsigned int i=0;i<MaxDup;i++){
		for(unsigned int j=0;j<MaxDup;j++){
			if( fabs(2*(xl1*i-xl2*j)/(xl1*i+xl2*j)) < MaxMisfit ){
				dupX1 = i;
				dupX2 = j;
				find = true;
				break;
			}
		}
		if( find ) break;
	}
	find = false;
	for(unsigned int i=0;i<MaxDup;i++){
		for(unsigned int j=0;j<MaxDup;j++){
			if( fabs(2*(yl1*i-yl2*j)/(yl1*i+yl2*j)) < MaxMisfit ){
				dupY1 = i;
				dupY2 = j;
				find = true;
				break;
			}
		}
		if( find ) break;
	}
	double final_xbox = (xl1*dupX1 + xl2*dupX2)/2;
    	double final_ybox = (yl1*dupY1 + yl2*dupY2)/2;
    	Mx1 = final_xbox / (xl1*dupX1);
    	My1 = final_ybox / (yl1*dupY1);
    	Mx2 = final_xbox / (xl2*dupX2);
    	My2 = final_ybox / (yl2*dupY2);
	// initialize the variables and pointers
	unsigned int nbAtom1 = this->_MyCrystal->getOrientedSystem()->getNbAtom();
	unsigned int nbAtom2 = this->_MyCrystal2->getOrientedSystem()->getNbAtom();
	unsigned int nbAtom_temp = nbAtom1*dupX1*dupY1 + nbAtom2*dupX2*dupY2;
	Atom *AtomList_temp = new Atom[nbAtom_temp];
	this->AtomList_G1 = new Atom[nbAtom1];
	this->AtomList_G2 = new Atom[nbAtom2];
	this->H1_G1 = new double[3];
	this->H2_G1 = new double[3];
	this->H3_G1 = new double[3];
	this->H1_G2 = new double[3];
	this->H2_G2 = new double[3];
	this->H3_G2 = new double[3];
	unsigned int *TagGrain = new unsigned int[nbAtom_temp];
	this->H1 = new double[3]; 
	this->H2 = new double[3]; 
	this->H3 = new double[3]; 
	this->H1[0] = final_xbox;
	this->H1[1] = 0.;
	this->H1[2] = 0.;
	this->H2[0] = 0.;
	this->H2[1] = final_ybox;
	this->H2[2] = 0.;
	this->H3[0] = 0.;
	this->H3[1] = 0.;
	this->H3[2] = this->_MyCrystal->getOrientedSystem()->getH3()[2] + this->_MyCrystal2->getOrientedSystem()->getH3()[2]+GBspace-DeltaZ;
	this->H1_G1[0] = Mx1*xl1;
	this->H1_G1[1] = 0.;
	this->H1_G1[2] = 0.;
	this->H2_G1[0] = 0.;
	this->H2_G1[1] = My1*yl1;
	this->H2_G1[2] = 0.;
	this->H3_G1[0] = 0.;
	this->H3_G1[1] = 0.;
	this->H3_G1[2] = this->_MyCrystal->getOrientedSystem()->getH3()[2];
	this->H1_G2[0] = Mx2*xl2;
	this->H1_G2[1] = 0.;
	this->H1_G2[2] = 0.;
	this->H2_G2[0] = 0.;
	this->H2_G2[1] = My2*yl2;
	this->H2_G2[2] = 0.;
	this->H3_G2[0] = 0.;
	this->H3_G2[1] = 0.;
	this->H3_G2[2] = this->_MyCrystal2->getOrientedSystem()->getH3()[2];
	this->nbAtomType = _MyCrystal->getNbAtomType();
	this->AtomType = _MyCrystal->getAtomType();
	this->AtomMass = _MyCrystal->getAtomMass();
	this->AtomCharge = _MyCrystal->getAtomCharge();
	this->IsCharge = _MyCrystal->getIsCharge();
	this->IsTilted = false;
	computeInverseCellVec();
	// paste the two grains
	double xpos, ypos, zpos, Lin;
	double fullFaceLength_1 = My1*(n1*Dir1_G1[1]+n2*Dir2_G1[1]);
	double slope1_1 = Dir1_G1[2]/(My1*Dir1_G1[1]);
	double slope2_1 = Dir2_G1[2]/(My1*Dir2_G1[1]);
	double fullFaceLength_2 = My2*(n1*Dir1_G2[1]+n2*Dir2_G2[1]);
	double slope1_2 = Dir1_G2[2]/(My2*Dir1_G2[1]);
	double slope2_2 = Dir2_G2[2]/(My2*Dir2_G2[1]);
	double Origin;
	unsigned int trueNbAt1 = 0;
	unsigned int trueNbAt2 = 0;
	if( slope1_1 > 0 ) Origin = this->_MyCrystal->getOrientedSystem()->getH3()[2]-DeltaZ;
	else Origin = this->_MyCrystal->getOrientedSystem()->getH3()[2];
	unsigned int at_count = 0;
	if( !this->_MyCrystal->getIsDoNotSep() ){
		// for crystals where ions can be separated
		for(unsigned int i=0;i<dupX1;i++){
			for(unsigned int j=0;j<dupY1;j++){
				for(unsigned int n=0;n<nbAtom1;n++){
					ypos = My1*(this->_MyCrystal->getOrientedSystem()->getAtom(n).pos.y+j*yl1);
					zpos = this->_MyCrystal->getOrientedSystem()->getAtom(n).pos.z;
					if( i == 0 && j == 0 ){
						this->AtomList_G1[n] = this->_MyCrystal->getOrientedSystem()->getAtom(n);
						this->AtomList_G1[n].pos.y = ypos;
						this->AtomList_G1[n].pos.x *= Mx1;
					}
					if( 0. <= fmod(ypos,fullFaceLength_1) && fmod(ypos,fullFaceLength_1) < (My1*n1*Dir1_G1[1]) ) Lin = zpos-slope1_1*ypos-Origin+slope1_1*((double) floor(ypos/fullFaceLength_1))*fullFaceLength_1;
					else if( (My1*n1*Dir1_G1[1]) <= fmod(ypos,fullFaceLength_1) && fmod(ypos,fullFaceLength_1) < fullFaceLength_1 ) Lin = zpos-(slope2_1*ypos)-Origin+((((double) floor(ypos/fullFaceLength_1))+1.)*slope2_1*fullFaceLength_1);
					if( Lin <= 0. ){
						AtomList_temp[at_count] = this->_MyCrystal->getOrientedSystem()->getAtom(n);
						AtomList_temp[at_count].pos.x = Mx1*(AtomList_temp[at_count].pos.x + i*xl1);
						AtomList_temp[at_count].pos.y = ypos;
						TagGrain[at_count] = 1;
						at_count += 1;
					}
				}
			}
		}
		trueNbAt1 = at_count;
		if( slope1_2 > 0 ) Origin = 0.;
		else Origin = DeltaZ;
		for(unsigned int i=0;i<dupX2;i++){
			for(unsigned int j=0;j<dupY2;j++){
				for(unsigned int n=0;n<nbAtom2;n++){
					ypos = My2*(this->_MyCrystal2->getOrientedSystem()->getAtom(n).pos.y+j*yl2);
					zpos = this->_MyCrystal2->getOrientedSystem()->getAtom(n).pos.z;
					if( i == 0 && j == 0 ){
						this->AtomList_G2[n] = this->_MyCrystal2->getOrientedSystem()->getAtom(n);
						this->AtomList_G2[n].pos.y = ypos;
						this->AtomList_G2[n].pos.x *= Mx2;
					}
					if( 0. <= fmod(ypos,fullFaceLength_2) && fmod(ypos,fullFaceLength_2) < (My2*n1*Dir1_G2[1]) ) Lin = zpos-slope1_2*ypos-Origin+slope1_2*((double) floor(ypos/fullFaceLength_2))*fullFaceLength_2;
					else if( (My2*n1*Dir1_G2[1]) <= fmod(ypos,fullFaceLength_2) && fmod(ypos,fullFaceLength_2) < fullFaceLength_2 ) Lin = zpos-(slope2_2*ypos)-Origin+((((double) floor(ypos/fullFaceLength_2))+1.)*slope2_2*fullFaceLength_2);
					if( Lin >= 0. ){
						AtomList_temp[at_count] = this->_MyCrystal2->getOrientedSystem()->getAtom(n);
						AtomList_temp[at_count].pos.x = Mx2*(AtomList_temp[at_count].pos.x + i*xl2);
						AtomList_temp[at_count].pos.y = ypos;
						AtomList_temp[at_count].pos.z += (this->_MyCrystal->getOrientedSystem()->getH3()[2]+GBspace-DeltaZ);
						TagGrain[at_count] = 2;
						at_count += 1;
					}
				}
			}
		}
		trueNbAt2 = at_count-trueNbAt1;
	}else{
		// do not sep case
		for(unsigned int i=0;i<dupX1;i++){
			for(unsigned int j=0;j<dupY1;j++){
				for(unsigned int n=0;n<nbAtom1;n++){
					if( i == 0 && j == 0 ){
						this->AtomList_G1[n] = this->_MyCrystal->getOrientedSystem()->getAtom(n);
						this->AtomList_G1[n].pos.y *= My1;
						this->AtomList_G1[n].pos.x *= Mx1;
					}
					if( this->_MyCrystal->getOrientedSystem()->getNotSepTag()[n][0] >= 0 ){
						ypos = My1*(this->_MyCrystal->getOrientedSystem()->getAtom(n).pos.y+j*yl1);
						zpos = this->_MyCrystal->getOrientedSystem()->getAtom(n).pos.z;
						if( 0. <= fmod(ypos,fullFaceLength_1) && fmod(ypos,fullFaceLength_1) < (My1*n1*Dir1_G1[1]) ) Lin = zpos-slope1_1*ypos-Origin+slope1_1*((double) floor(ypos/fullFaceLength_1))*fullFaceLength_1;
						else if( (My1*n1*Dir1_G1[1]) <= fmod(ypos,fullFaceLength_1) && fmod(ypos,fullFaceLength_1) < fullFaceLength_1 ) Lin = zpos-(slope2_1*ypos)-Origin+((((double) floor(ypos/fullFaceLength_1))+1.)*slope2_1*fullFaceLength_1);
						if( Lin <= 0. ){
							AtomList_temp[at_count] = this->_MyCrystal->getOrientedSystem()->getAtom(n);
							AtomList_temp[at_count].pos.x = Mx1*(AtomList_temp[at_count].pos.x + i*xl1);
							AtomList_temp[at_count].pos.y = ypos;
							TagGrain[at_count] = 1;
							at_count += 1;
							if( this->_MyCrystal->getOrientedSystem()->getNotSepTag()[n][0] > 0 ){
								for(unsigned int ns=0;ns<this->_MyCrystal->getOrientedSystem()->getNotSepTag()[n][0];ns++){
									AtomList_temp[at_count] = this->_MyCrystal->getOrientedSystem()->getAtom(this->_MyCrystal->getOrientedSystem()->getNotSepTag()[n][ns+1]);
									AtomList_temp[at_count].pos.x = Mx1*(AtomList_temp[at_count].pos.x + i*xl1);
									AtomList_temp[at_count].pos.y = My1*(AtomList_temp[at_count].pos.y + j*yl1);
									TagGrain[at_count] = 1;
									at_count += 1;
								}
							}
						}
					}
				}
			}
		}
		trueNbAt1 = at_count;
		if( slope1_2 > 0 ) Origin = 0.;
		else Origin = DeltaZ;
		for(unsigned int i=0;i<dupX2;i++){
			for(unsigned int j=0;j<dupY2;j++){
				for(unsigned int n=0;n<nbAtom2;n++){
					if( i == 0 && j == 0 ){
						this->AtomList_G2[n] = this->_MyCrystal2->getOrientedSystem()->getAtom(n);
						this->AtomList_G2[n].pos.y *= My2;
						this->AtomList_G2[n].pos.x *= Mx2;
					}
					if( this->_MyCrystal2->getOrientedSystem()->getNotSepTag()[n][0] >= 0 ){
						ypos = My2*(this->_MyCrystal2->getOrientedSystem()->getAtom(n).pos.y+j*yl2);
						zpos = this->_MyCrystal2->getOrientedSystem()->getAtom(n).pos.z;
						if( 0. <= fmod(ypos,fullFaceLength_2) && fmod(ypos,fullFaceLength_2) < (My2*n1*Dir1_G2[1]) ) Lin = zpos-slope1_2*ypos-Origin+slope1_2*((double) floor(ypos/fullFaceLength_2))*fullFaceLength_2;
						else if( (My2*n1*Dir1_G2[1]) <= fmod(ypos,fullFaceLength_2) && fmod(ypos,fullFaceLength_2) < fullFaceLength_2 ) Lin = zpos-(slope2_2*ypos)-Origin+((((double) floor(ypos/fullFaceLength_2))+1.)*slope2_2*fullFaceLength_2);
						if( Lin >= 0. ){
							AtomList_temp[at_count] = this->_MyCrystal2->getOrientedSystem()->getAtom(n);
							AtomList_temp[at_count].pos.x = Mx2*(AtomList_temp[at_count].pos.x + i*xl2);
							AtomList_temp[at_count].pos.y = ypos;
							AtomList_temp[at_count].pos.z += (this->_MyCrystal->getOrientedSystem()->getH3()[2]-DeltaZ);
							TagGrain[at_count] = 2;
							at_count += 1;
							if( this->_MyCrystal2->getOrientedSystem()->getNotSepTag()[n][0] > 0 ){
								for(unsigned int ns=0;ns<this->_MyCrystal2->getOrientedSystem()->getNotSepTag()[n][0];ns++){
									AtomList_temp[at_count] = this->_MyCrystal2->getOrientedSystem()->getAtom(this->_MyCrystal2->getOrientedSystem()->getNotSepTag()[n][ns+1]);
									AtomList_temp[at_count].pos.x = Mx2*(AtomList_temp[at_count].pos.x + i*xl2);
									AtomList_temp[at_count].pos.y = My2*(AtomList_temp[at_count].pos.y + j*yl2);
									AtomList_temp[at_count].pos.z += (this->_MyCrystal->getOrientedSystem()->getH3()[2]-DeltaZ);
									TagGrain[at_count] = 2;
									at_count += 1;
								}
							}
						}
					}
				}
			}
		}
		trueNbAt2 = at_count-trueNbAt1;
	}
	// Initialize the two grains
	this->Grain1 = new AtomicSystem(this->AtomList_G1,nbAtom1,_MyCrystal,this->H1_G1,this->H2_G1,this->H3_G1);
	this->Grain2 = new AtomicSystem(this->AtomList_G2,nbAtom2,_MyCrystal2,this->H1_G2,this->H2_G2,this->H3_G2);
	this->AreGrainsDefined = true;
	// verify the stoichiometry
	unsigned int *currentStoich = new unsigned int[this->_MyCrystal->getNbAtomType()];
	for(unsigned int i=0;i<this->_MyCrystal->getNbAtomType();i++) currentStoich[i] = 0;
	for(unsigned int i=0;i<at_count;i++){
		for(unsigned int t=0;t<this->_MyCrystal->getNbAtomType();t++){
			if( AtomList_temp[i].type_uint == t+1 ){
				currentStoich[t] += 1;
				break;
			}
		}
	}
	bool stoich = true;
	for(unsigned int i=0;i<this->_MyCrystal->getNbAtomType();i++){
		if( fabs(((double) currentStoich[i]/at_count) - ((double) this->_MyCrystal->getStoich()[i]/this->_MyCrystal->getNbAtom()) ) > 1e-9 ){
			stoich = false;
			cout << "Warning the stoichiometry is not the same than parent crystal, ";
			for(unsigned int t=0;t<this->_MyCrystal->getNbAtomType();t++){
				cout << "number of " << this->_MyCrystal->getAtomType(t+1) << " : " << currentStoich[t] << ", ";
			}
			cout << endl;
			break;
		}
	}
	// adjsut stoichiometry
	if( !stoich ){
		int *adjust = new int[this->_MyCrystal->getNbAtomType()];
		int *adjust_tmp = new int[this->_MyCrystal->getNbAtomType()];
		for(unsigned int i=0;i<this->_MyCrystal->getNbAtomType();i++){
			adjust[i] = this->_MyCrystal->getNbAtom()*round((double) at_count/this->_MyCrystal->getNbAtom())*((double) this->_MyCrystal->getStoich()[i]/this->_MyCrystal->getNbAtom())-currentStoich[i];
			adjust_tmp[i] = adjust[i];
		}
		// modify to have only negative value (we only want to remove atoms)
		for(unsigned int i=0;i<this->_MyCrystal->getNbAtomType();i++){
			if( adjust[i] > 0 ){
				for(unsigned int t=0;t<this->_MyCrystal->getNbAtomType();t++){
					adjust[t] -= adjust_tmp[i]*((double) this->_MyCrystal->getStoich()[t]/this->_MyCrystal->getStoich()[i]);
				}
				adjust_tmp[i] = adjust[i];
			}
		}
		// neighboring search in the GB zone for ions type to remove
		double GBup = this->_MyCrystal->getOrientedSystem()->getH3()[2]+3.;
		double GBdown = this->_MyCrystal->getOrientedSystem()->getH3()[2]-DeltaZ-3.;
		vector<vector<double>> Neigh;
		vector<int> Type2Rm;
		for(unsigned int i=0;i<this->_MyCrystal->getNbAtomType();i++){
			if( adjust[i] < 0 ){
				Type2Rm.push_back(i+1);
				Neigh.push_back(vector<double>());
			}
		}
		double rc_squared = pow(3.,2.);
		double d_squared;
		for(unsigned int i=0;i<at_count;i++){
			for(unsigned int t=0;t<Type2Rm.size();t++){
				if( AtomList_temp[i].pos.z < GBup && AtomList_temp[i].pos.z > GBdown && AtomList_temp[i].type_uint == Type2Rm[t] ){
					xpos = AtomList_temp[i].pos.x;
					ypos = AtomList_temp[i].pos.y;
					zpos = AtomList_temp[i].pos.z;
					for(unsigned int j=0;j<at_count;j++){
						if( i != j && AtomList_temp[j].pos.z < GBup && AtomList_temp[j].pos.z > GBdown  && AtomList_temp[j].type_uint == Type2Rm[t] ){
							d_squared = pow(xpos-AtomList_temp[j].pos.x,2.) + pow(ypos-AtomList_temp[j].pos.y,2.) + pow(zpos-AtomList_temp[j].pos.z,2.); 
							if( d_squared < rc_squared ){
								Neigh[t].push_back(d_squared);
								Neigh[t].push_back(j);
							}
						}
					}
				}
			}
		}
		this->nbAtom = at_count;
		for(unsigned int t=0;t<Type2Rm.size();t++){
			this->MT->sort(Neigh[t],0,2,Neigh[t]);
			this->nbAtom += adjust[Type2Rm[t]-1];
		}
		this->AtomList = new Atom[this->nbAtom];
		unsigned int at_count2 = 0;
		bool store;
		for(unsigned int i=0;i<at_count;i++){
			store = true;
			for(unsigned int t=0;t<Type2Rm.size();t++){
				if( AtomList_temp[i].type_uint == Type2Rm[t] ){
					for(unsigned int n=0;n<abs(adjust[Type2Rm[t]-1]);n++){
						if( round(Neigh[t][n*2+1]) == i ){
							store = false;
							if( TagGrain[i] == 1 ) trueNbAt1 -= 1;
							else if( TagGrain[i] == 2 ) trueNbAt2 -= 1;
							break;
						}
					}
				}
			}
			if( store ){
				this->AtomList[at_count2] = AtomList_temp[i];
				if( TagGrain[i] == 2 ) this->AtomList[at_count2].pos.z += GBspace;
				at_count2 += 1;
			}
		}
		delete[] adjust;
		delete[] adjust_tmp;
	}else{
		this->nbAtom = at_count;
		this->AtomList = new Atom[this->nbAtom];
		for(unsigned int i=0;i<this->nbAtom;i++){
			this->AtomList[i] = AtomList_temp[i];
			if( TagGrain[i] == 2 ) this->AtomList[i].pos.z += GBspace;
		}
	}
	string h_a_str = to_string(h_a);
	string k_a_str = to_string(k_a);
	string l_a_str = to_string(l_a);
	string theta_str = to_string(theta*180/M_PI);
	string h_p_str = to_string(h_p);
	string k_p_str = to_string(k_p);
	string l_p_str = to_string(l_p);
	string h_f1_str = to_string(FacetsType[0]);
	string k_f1_str = to_string(FacetsType[1]);
	string l_f1_str = to_string(FacetsType[2]);
	string h_f2_str = to_string(FacetsType[3]);
	string k_f2_str = to_string(FacetsType[4]);
	string l_f2_str = to_string(FacetsType[5]);
	string nbAt1_str = to_string(trueNbAt1);
	string nbAt2_str = to_string(trueNbAt2);
	this->File_Heading = " # ["+h_a_str+k_a_str+l_a_str+"]"+theta_str+"°("+h_p_str+k_p_str+l_p_str+") "+crystalName+" grain boundary\n # The grain boundary is faceted with ("+h_f1_str+k_f1_str+l_f1_str+") and ("+h_f2_str+k_f2_str+l_f2_str+") planes\n # The lower grain contains "+nbAt1_str+" atoms and the upper grain contains "+nbAt2_str+" atoms\n";
	this->Grain1->set_File_Heading(" # Lower grain of the ["+h_a_str+k_a_str+l_a_str+"]"+theta_str+"°("+h_p_str+k_p_str+l_p_str+") "+crystalName+" grain boundary\n");
	this->Grain2->set_File_Heading(" # Upper grain of the ["+h_a_str+k_a_str+l_a_str+"]"+theta_str+"°("+h_p_str+k_p_str+l_p_str+") "+crystalName+" grain boundary\n");

	this->xl1 *= Mx1; 
	this->xl2 *= Mx2; 
	this->yl1 *= My1; 
	this->yl2 *= My2; 
	// delete allocated memory
	delete[] currentStoich;
	delete[] AtomList_temp;
	delete[] Dir1_G1;
	delete[] Dir2_G1;
	delete[] Dir1_G2;
	delete[] Dir2_G2;
	delete[] TagGrain;
}
	
// Constructor for bicrystal with plane GB with given misorientation and GB plane
Bicrystal::Bicrystal(const string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p, bool rationalize):h_a(h_a), k_a(k_a), l_a(l_a), theta(theta), h_p(h_p), k_p(k_p), l_p(l_p){
	read_params();
	double GBspace = 0.; //TODO maybe define elsewhere (in a file in the main prog for example..)
	double MaxMisfit = 0.02;
	unsigned int MaxDup = 50;
	bool FullGrains = true; // TODO as the variables above
	this->MT = new MathTools;
	cout << "constructing crystals, " ;
	setOrientedCrystals(crystalName, rationalize);
	cout << "done" << endl;
	// search the number of linear combination for the two system to have the same x y length
	this->xl1 = this->_MyCrystal->getOrientedSystem()->getH1()[0];
	this->xl2 = this->_MyCrystal2->getOrientedSystem()->getH1()[0];
	this->yl1 = this->_MyCrystal->getOrientedSystem()->getH2()[1];
	this->yl2 = this->_MyCrystal2->getOrientedSystem()->getH2()[1];
	bool find = false;
	for(unsigned int i=0;i<MaxDup;i++){
		for(unsigned int j=0;j<MaxDup;j++){
			if( fabs(2*(xl1*i-xl2*j)/(xl1*i+xl2*j)) < MaxMisfit ){
				dupX1 = i;
				dupX2 = j;
				find = true;
				break;
			}
		}
		if( find ) break;
	}
	find = false;
	for(unsigned int i=0;i<MaxDup;i++){
		for(unsigned int j=0;j<MaxDup;j++){
			if( fabs(2*(yl1*i-yl2*j)/(yl1*i+yl2*j)) < MaxMisfit ){
				dupY1 = i;
				dupY2 = j;
				find = true;
				break;
			}
		}
		if( find ) break;
	}
	double final_xbox = (xl1*dupX1 + xl2*dupX2)/2;
    	double final_ybox = (yl1*dupY1 + yl2*dupY2)/2;
    	Mx1 = final_xbox / (xl1*dupX1);
    	My1 = final_ybox / (yl1*dupY1);
    	Mx2 = final_xbox / (xl2*dupX2);
    	My2 = final_ybox / (yl2*dupY2);
	cout << "Duplicates : " << dupX1 << " " << dupY1 << " " << dupX2 << " " << dupY2 << endl;
	cout << "Misfits : " << Mx1 << " " << My1 << " " << Mx2 << " " << My2 << endl;
	generateCSL();
	// initialize the variables and pointers
	unsigned int nbAtom1 = this->_MyCrystal->getOrientedSystem()->getNbAtom();
	unsigned int nbAtom2 = this->_MyCrystal2->getOrientedSystem()->getNbAtom();
	unsigned int nbAtom1_G, nbAtom2_G;
	this->nbAtom = nbAtom1*dupX1*dupY1 + nbAtom2*dupX2*dupY2;
	this->AtomList = new Atom[this->nbAtom];
	this->H1 = new double[3]; 
	this->H2 = new double[3]; 
	this->H3 = new double[3]; 
	this->H1_G1 = new double[3];
	this->H2_G1 = new double[3];
	this->H3_G1 = new double[3];
	this->H1_G2 = new double[3];
	this->H2_G2 = new double[3];
	this->H3_G2 = new double[3];
	this->H1[0] = final_xbox;
	this->H1[1] = 0.;
	this->H1[2] = 0.;
	this->H2[0] = 0.;
	this->H2[1] = final_ybox;
	this->H2[2] = 0.;
	this->H3[0] = 0.;
	this->H3[1] = 0.;
	this->H3[2] = this->_MyCrystal->getOrientedSystem()->getH3()[2] + this->_MyCrystal2->getOrientedSystem()->getH3()[2]+GBspace;
	this->H1_G1[1] = 0.;
	this->H1_G1[2] = 0.;
	this->H2_G1[0] = 0.;
	this->H2_G1[2] = 0.;
	this->H3_G1[0] = 0.;
	this->H3_G1[1] = 0.;
	this->H3_G1[2] = this->_MyCrystal->getOrientedSystem()->getH3()[2];
	this->H1_G2[1] = 0.;
	this->H1_G2[2] = 0.;
	this->H2_G2[0] = 0.;
	this->H2_G2[2] = 0.;
	this->H3_G2[0] = 0.;
	this->H3_G2[1] = 0.;
	this->H3_G2[2] = this->_MyCrystal2->getOrientedSystem()->getH3()[2];
	if( FullGrains ){
		this->H1_G1[0] = Mx1*xl1*dupX1;
		this->H2_G1[1] = My1*yl1*dupY1;
		this->H1_G2[0] = Mx2*xl2*dupX2;
		this->H2_G2[1] = My2*yl2*dupY2;
		nbAtom1_G = nbAtom1*dupX1*dupY1;
		nbAtom2_G = nbAtom2*dupX2*dupY2;
	}else{
		this->H1_G1[0] = Mx1*xl1;
		this->H2_G1[1] = My1*yl1;
		this->H1_G2[0] = Mx2*xl2;
		this->H2_G2[1] = My2*yl2;
		nbAtom1_G = nbAtom1;
		nbAtom2_G = nbAtom2;
	}
	this->AtomList_G1 = new Atom[nbAtom1_G];
	this->AtomList_G2 = new Atom[nbAtom2_G];
	this->nbAtomType = _MyCrystal->getNbAtomType();
	this->AtomType = _MyCrystal->getAtomType();
	this->AtomMass = _MyCrystal->getAtomMass();
	this->AtomCharge = _MyCrystal->getAtomCharge();
	this->IsCharge = _MyCrystal->getIsCharge();
	this->IsTilted = false;
	computeInverseCellVec();
	// paste the two grains
	for(unsigned int i=0;i<dupX1;i++){
		for(unsigned int j=0;j<dupY1;j++){
			for(unsigned int n=0;n<nbAtom1;n++){
				this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n] = this->_MyCrystal->getOrientedSystem()->getAtom(n);
				this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.x = Mx1*(this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.x + i*xl1);
				this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.y = My1*(this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.y + j*yl1);
				//this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.x = (this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.x + i*xl1);
				//this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.y = (this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.y + j*yl1);
				if( !FullGrains && i == 0 && j == 0 ) this->AtomList_G1[n] = this->AtomList[n];
				else this->AtomList_G1[i*dupY1*nbAtom1+j*nbAtom1+n] = this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n];
			}
		}
	}
	for(unsigned int i=0;i<dupX2;i++){
		for(unsigned int j=0;j<dupY2;j++){
			for(unsigned int n=0;n<nbAtom2;n++){
				this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n] = this->_MyCrystal2->getOrientedSystem()->getAtom(n);
				this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.x = Mx2*(this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.x + i*xl2);
				this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.y = My2*(this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.y + j*yl2);
				//this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.x = (this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.x + i*xl2);
				//this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.y = (this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.y + j*yl2);
				if( !FullGrains && i == 0 && j == 0 ) this->AtomList_G2[n] = this->AtomList[nbAtom1*dupX1*dupY1+n];
				else this->AtomList_G2[i*dupY2*nbAtom2+j*nbAtom2+n] = this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n];
				this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.z -= this->H3_G2[2];
				//this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.z += (this->_MyCrystal->getOrientedSystem()->getH3()[2]+GBspace);
			}
		}
	}
	// construct the two grains
	this->Grain1 = new AtomicSystem(this->AtomList_G1,nbAtom1_G,this->_MyCrystal,this->H1_G1,this->H2_G1,this->H3_G1);
	this->Grain2 = new AtomicSystem(this->AtomList_G2,nbAtom2_G,this->_MyCrystal2,this->H1_G2,this->H2_G2,this->H3_G2);
	this->AreGrainsDefined = true;
	string h_a_str = to_string(this->h_a);
	string k_a_str = to_string(this->k_a);
	string l_a_str = to_string(this->l_a);
	string theta_str = to_string(this->theta*180/M_PI);
	string h_p_str = to_string(this->h_p);
	string k_p_str = to_string(this->k_p);
	string l_p_str = to_string(this->l_p);
	this->File_Heading = " # ["+h_a_str+k_a_str+l_a_str+"]"+theta_str+"°("+h_p_str+k_p_str+l_p_str+") "+crystalName+" grain boundary\n";
	this->Grain1->set_File_Heading(" # Lower grain of the ["+h_a_str+k_a_str+l_a_str+"]"+theta_str+"°("+h_p_str+k_p_str+l_p_str+") "+crystalName+" grain boundary\n");
	this->Grain2->set_File_Heading(" # Upper grain of the ["+h_a_str+k_a_str+l_a_str+"]"+theta_str+"°("+h_p_str+k_p_str+l_p_str+") "+crystalName+" grain boundary\n");
	this->xl1 *= Mx1; 
	this->xl2 *= Mx2; 
	this->yl1 *= My1; 
	this->yl2 *= My2; 
}

Bicrystal::Bicrystal(const string& filename, const string NormalDir, const string CrystalName):AtomicSystem(filename){
	read_params();
	if( NormalDir != "z" && NormalDir != "y" && NormalDir != "x" && NormalDir != "Z" && NormalDir != "Y" && NormalDir != "X" ){
		cerr << "The direction normal to the GB have to be : x,X,y,Y,z or Z" << endl;
		exit(EXIT_FAILURE);
	}
	setCrystal(CrystalName);
	this->NormalDir = NormalDir;
	this->IsMassDensity = false;
	this->CA = new ComputeAuxiliary(this);
	searchGBPos();
	ComputeExcessVolume();
}

Bicrystal::Bicrystal(const string& filename, const string CrystalName):AtomicSystem(filename){
	read_params();
	setCrystal(CrystalName);
	this->CA = new ComputeAuxiliary(this);
	this->IsCA = true;
}

Bicrystal::Bicrystal(const string& filename):AtomicSystem(filename){
	read_params();
	this->CA = new ComputeAuxiliary(this);
	this->IsCA = true;
}

// search CSL primitive lattice based on the work of Xie et al. 2020 and Bonnet and Rolland 1975
//void Bicrystal::searchCSL(int h_a_func, int k_a_func, int l_a_func, double theta_func, int *CSL_vec, unsigned int verbose){
bool Bicrystal::searchCSL(double *rot_ax_func, double theta_func, int *CSL_vec, unsigned int verbose){
	cout << "Searching CSL lattice" << endl;
	this->CSL_Basis = new double[9];
	this->IsCSL_Basis = true;
	double *a1 = new double[9]; // basis vector of lattice 1 
	double *a1_inv = new double[9]; // basis vector of lattice 1 
	double *a2 = new double[9]; // basis vector of lattice 2 in coordinate of lattice 1 
	double *a2_a = new double[9]; // auxiliary lattice for CSL base computation 
	double *R = new double[9]; 
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
	double *DSC_Basis = new double[9]; 
	double *DSC_inv = new double[9]; 
	double *U_a = new double[9];
	double *searchVec = new double[3];
	double *CSL_Basis2 = new double[9]; 
	double *CSL_Basis_temp = new double[9]; 
	double *buffer_mat = new double[9]; 
	int *k = new int[3];
	double *rot_ax = new double[3];
	double *Known_CSL = new double[3];
	double *buffer_vec = new double[3];
	int *buffer_vec_int = new int[3];
	double NormRotAx = 0.;
	double tmp, tmp_1, tmp_2, BasePrec;
	for(unsigned int i=0;i<3;i++){
		a1[i*3] = this->_MyCrystal->getA1()[i];
		a1[i*3+1] = this->_MyCrystal->getA2()[i];
		a1[i*3+2] = this->_MyCrystal->getA3()[i];
		Known_CSL[i] = CSL_vec[0]*a1[i*3] + CSL_vec[1]*a1[i*3+1] + CSL_vec[2]*a1[i*3+2];
		a2[i*3] = this->_MyCrystal->getA1()[i];
		a2[i*3+1] = this->_MyCrystal->getA2()[i];
		a2[i*3+2] = this->_MyCrystal->getA3()[i];
	}
	// renormalize rot axis (secure thing) and compute rotation matrix
	for(unsigned int i=0;i<3;i++) NormRotAx += pow(rot_ax_func[i],2.);
	NormRotAx = sqrt(NormRotAx);
	for(unsigned int i=0;i<3;i++) rot_ax_func[i] /= NormRotAx;
	MT->Vec2rotMat(rot_ax_func,theta_func,R);
	unsigned int L, M, N;
	bool BaseFound = false;
	bool IsFindIntVec, BaseAligned;
	bool FirstBase = true;
	double tol_DSC = 1e-6;
	double tol_IntVec, sp;
	MT->invert3x3(a1,a1_inv);
	MT->MatDotMat(R, a2, a2);
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
		MT->MatDotVec(a1_inv,searchVec,searchVec);
		L = MT->find_integer_vector(searchVec,tol_IntVec,this->SigmaMax,k,IsFindIntVec);
		if( !IsFindIntVec ){
			if( verbose == 2 ) cout << "Fail ! Increasing tolerance !" << endl;
			if( exp > 1.5 )	exp -= .5;
			else exp -= .1;
			continue;
		}
		solve_DSC(k,L,a1,Ei,tol_DSC);
		// search if this basis is also good for a2_a_2
		for(unsigned int i=0;i<3;i++) searchVec[i] = a2_a[i*3+1];
		MT->invert3x3(Ei,Ei_inv);
		MT->MatDotVec(Ei_inv,searchVec,searchVec);
		L = MT->find_integer_vector(searchVec,tol_IntVec,this->SigmaMax,k,IsFindIntVec);
		if( !IsFindIntVec ){
			if( verbose == 2 ) cout << "Fail ! Increasing tolerance !" << endl;
			if( exp > 1.5 )	exp -= .5;
			else exp -= .1;
			continue;
		}
		if( L == 1 ){
			// Yes, test if the basis is also good for a2_a_3
			for(unsigned int i=0;i<3;i++) searchVec[i] = a2_a[i*3+2];
			MT->MatDotVec(Ei_inv,searchVec,searchVec);
			M = MT->find_integer_vector(searchVec,tol_IntVec,this->SigmaMax,k,IsFindIntVec);
			if( !IsFindIntVec ){
				if( verbose == 2 ) cout << "Fail ! Increasing tolerance !" << endl;
				if( exp > 1.5 )	exp -= .5;
				else exp -= .1;
				continue;
			}
			if( M != 1 ) solve_DSC(k,M,Ei,DSC_Basis,tol_DSC); // No, use this vector for DSC basis
			else for(unsigned int i=0;i<9;i++) DSC_Basis[i] = Ei[i];
		}else{
			// No, use this vector for DSC basis
			solve_DSC(k,L,Ei,Fi,tol_DSC);
			// search if this basis is also good for a2_a_2
			for(unsigned int i=0;i<3;i++) searchVec[i] = a2_a[i*3+2];
			MT->invert3x3(Fi,Fi_inv);
			MT->MatDotVec(Fi_inv,searchVec,searchVec);
			M = MT->find_integer_vector(searchVec,tol_IntVec,this->SigmaMax,k,IsFindIntVec);
			if( !IsFindIntVec ){
				if( verbose == 2 ) cout << "Fail ! Increasing tolerance !" << endl;
				if( exp > 1.5 )	exp -= .5;
				else exp -= .1;
				continue;
			}
			if( M != 1 ){
				solve_DSC(k,M,Fi,DSC_Basis,tol_DSC); // No, use this vector for DSC basis
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
			this->sigma = fabs(MT->det(U_1_temp));
			if( verbose == 1 ){
				cout << "We have found CSL basis corresponding to Sigma = " << this->sigma << endl;
				MT->printMat(U_1_temp);
			}
			if( OneBaseFound && this->sigma < 1.5 ){
				BaseFound = true;
				break;
			}
			BaseFound = false;
			MT->invert3x3(CSL_Basis_temp,buffer_mat);
			MT->MatDotVec(buffer_mat,Known_CSL,buffer_vec);
			tmp = 0.;
			tmp_1 = 0.;
			tmp_2 = 0.;
			for(unsigned int i=0;i<9;i++) tmp_2 += fabs(U_1_temp[i]-round(U_1_temp[i]));
			for(unsigned int i=0;i<3;i++){
				tmp += fabs(buffer_vec[i]-round(buffer_vec[i]));
				tmp_1 += fabs(round(buffer_vec[i]));
			}
			if( ( tmp < this->tolpos_known_CSL && tmp_1 > this->tolpos_known_CSL && tmp_2 < this->tol_CSL_integer ) || tmp_2 < tol_DSC ){
				BaseFound = true;
				OneBaseFound = true;
				if( verbose == 1 ) cout << "Good base found !" << endl;
				for(unsigned int ii=0;ii<9;ii++){
					CSL_Basis[ii] = CSL_Basis_temp[ii];
					U_1[ii] = U_1_temp[ii];
				}
			}
			// one of the CSL basis vector should be aligned with the rotation axis => verify that the backup base satisfy this criterium
			for(unsigned int v=0;v<3;v++){
				sp = 0.;
				NormRotAx = 0.;
				BaseAligned = false;
				for(unsigned int i=0;i<3;i++){
					sp += CSL_Basis_temp[i*3+v]*rot_ax_func[i];
					NormRotAx += pow(CSL_Basis_temp[i*3+v],2.);
				}
				if( fabs(fabs((sp/sqrt(NormRotAx)))-1.) < this->tolAlignment_CSL ){
					BaseAligned = true;
					break;
				}
			}
			if( FirstBase && this->sigma < this->SigmaMax && BaseAligned ){
				BasePrec = tmp_2;
				for(unsigned int i=0;i<9;i++) Backup_U1[i] = U_1_temp[i];
				FirstBase = false;
			}else if( tmp_2 < BasePrec && this->sigma < this->SigmaMax && BaseAligned ){
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
		this->sigma = fabs(MT->det(U_1));
		//for(unsigned int i=0;i<9;i++){// still something TODO for sharing the CSL to the two lattices (for the moment CSL coincide with grain 2) 
		//	CSL_Basis[i] += CSL_Basis2[i];
		//	CSL_Basis[i] /= 2.;
		//}
		cout << "Success ! We find a CSL corresponding to sigma =" << this->sigma << ", using a tolerance of " << tol_IntVec << endl;
		success=true;
	}else if( !FirstBase ){
		MT->MatDotMat(a1,Backup_U1,CSL_Basis);
		this->sigma = fabs(MT->det(Backup_U1));
		MT->printMat(Backup_U1);
		cout << "We find a CSL corresponding to sigma =" << this->sigma << endl;
		success=true;
	}else{
		cout << "Fail to find CSL basis" << endl;
		success=false;
	}
	cout << "CSL" << endl;
	MT->printMat(CSL_Basis);
	delete[] a1;
	delete[] a1_inv;
	delete[] a2;
	delete[] a2_a;
	delete[] U;
	delete[] U_inv;
	delete[] U_1;
	delete[] U_2;
	delete[] Ei;
	delete[] Ei_inv;
	delete[] Fi;
	delete[] Fi_inv;
	delete[] DSC_Basis;
	delete[] DSC_inv;
	delete[] U_a;
	delete[] searchVec;
	delete[] k;
	delete[] rot_ax;
	delete[] R;
	delete[] buffer_vec;
	delete[] Known_CSL;
	delete[] CSL_Basis_temp;
	delete[] CSL_Basis2;
	delete[] U_1_temp;
	return success;
}

// solve DSC equation based on the work of Bonnet and Durand 1975
void Bicrystal::solve_DSC(const int *u, const unsigned int L, const double *B, double *DSC_Base, double tol){
	unsigned int *GCD = new unsigned int[3];
	unsigned int g_v, g_lambda, g_mu, gamma, alpha, beta;
	double buffer;
	bool found = false;
	unsigned int count = 0;
	g_v = L/MT->gcd(abs(u[2]), L);
	GCD[0] = abs(u[1]);
	GCD[1] = abs(u[2]);
	GCD[2] = L;
	g_lambda = MT->gcd_mult(GCD,3);
	g_mu = MT->gcd(abs(u[2]), L)/g_lambda;
	for(unsigned int i=0;i<g_v+1;i++){
		gamma = i;
		buffer = (double) ((int) g_mu*u[1]- (int) gamma*u[2])/ (int) L;
		if( fabs(buffer-round(buffer)) < tol ) break;
	}
	alpha = 0;
	while( !found && alpha < g_mu+1 ){
		for(unsigned int i=0;i<g_v+1;i++){
			beta = i;
			buffer = (double) ((int) g_lambda*u[0] - (int) alpha*u[1] - (int) beta*u[2]) / (int) L;
			if( fabs(buffer-round(buffer)) < tol ){
				found = true;
				break;
			}
		}
		if( found ) break;
		else alpha++;
	}
	if( !found ) cout << "Failed to find DSC basis !" << endl;
	for(unsigned int i=0;i<3;i++){
		DSC_Base[i*3] = B[i*3] / g_lambda;
		DSC_Base[i*3+1] = (alpha*B[i*3]/(g_lambda * g_mu)) + B[i*3+1]/g_mu;
		DSC_Base[i*3+2] = (B[i*3]*(alpha*gamma+beta*g_mu)/(g_mu*g_v*g_lambda)) + (B[i*3+1]*gamma/(g_mu*g_v)) + B[i*3+2]/g_v;
	}

	delete[] GCD;
}
void Bicrystal::generateCSL(){
	int CLMax = 100;
	double *pos = new double[3];
	double xhi = this->_MyCrystal->getOrientedSystem()->getH1()[0]*dupX1;
	double yhi = this->_MyCrystal->getOrientedSystem()->getH2()[1]*dupY1;
	double zhi = this->_MyCrystal->getOrientedSystem()->getH3()[2]+this->_MyCrystal2->getOrientedSystem()->getH3()[2];
	double tolpos = 1.e-1;
	for(int i=-CLMax;i<CLMax;i++){
		for(int j=-CLMax;j<CLMax;j++){
			for(int k=-CLMax;k<CLMax;k++){
				for(unsigned int ii=0;ii<3;ii++) pos[ii] = i*CSL_Basis[ii*3] + j*CSL_Basis[ii*3+1] + k*CSL_Basis[ii*3+2];
				if( pos[0] < xhi+tolpos && pos[0] > -tolpos && pos[1] < yhi+tolpos && pos[1] > -tolpos && pos[2] > -zhi-tolpos && pos[2] < zhi+tolpos ){
					CSL.push_back(new double[3]);
					for(unsigned int ii=0;ii<3;ii++) CSL[CSL.size()-1][ii] = i*CSL_Basis[ii*3] + j*CSL_Basis[ii*3+1] + k*CSL_Basis[ii*3+2];
				}
			}
		}
	}
	delete[] pos;
}

void Bicrystal::printCSL(const std::string filename){
	ofstream writefile(filename);
	writefile << " # File generated using AtomHic\n";
	writefile << this->File_Heading;
        writefile << "\n\t" << this->nbAtom+CSL.size() << "\tatoms\n\t" << this->nbAtomType+1 << "\tatom types\n\n\t0.000000000000\t" << this->H1[0] << "\txlo xhi\n\t0.000000000000\t" << H2[1] << "\tylo yhi\n\t0.000000000000\t" << H3[2] << "\tzlo zhi\n";
	if( this->IsTilted ) writefile << "\t" << H2[0] << "\t" << H3[0] << "\t" << H3[1] << "\txy xz yz\n";
	writefile << "\nMasses\n\n";
	for(unsigned int i=0;i<this->nbAtomType;i++) writefile << "\t" << i+1 << "\t" << this->AtomMass[i] << "\t# " << this->AtomType[i] << "\n";
	writefile << "\t" << this->nbAtomType+1 << "\t0\t# CSL\n";

	if( IsCharge ){
		writefile << "\nAtoms # charge\n\n";
		for(unsigned int i=0;i<this->nbAtom;i++) writefile << i+1 << "\t" << this->AtomList[i].type_uint << "\t" << this->AtomCharge[this->AtomList[i].type_uint-1] << "\t" << this->AtomList[i].pos.x << "\t" << this->AtomList[i].pos.y << "\t" << this->AtomList[i].pos.z << "\n"; 
		for(unsigned int i=0;i<CSL.size();i++) writefile << i+1+this->nbAtom << "\t" << this->nbAtomType+1 << "\t0\t" << CSL[i][0] << "\t" << CSL[i][1] << "\t" << CSL[i][2] << "\n";
	}else{
		writefile << "\nAtoms # atomic\n\n";
		for(unsigned int i=0;i<this->nbAtom;i++) writefile << i+1 << "\t" << this->AtomList[i].type_uint << "\t" << this->AtomList[i].pos.x << "\t" << this->AtomList[i].pos.y << "\t" << this->AtomList[i].pos.z << "\n"; 
		for(unsigned int i=0;i<CSL.size();i++) writefile << i+1+this->nbAtom << "\t" << this->nbAtomType+1 << "\t" << CSL[i][0] << "\t" << CSL[i][1] << "\t" << CSL[i][2] << "\n";
	}
	writefile.close();

}

// search the rotation axis and angle corresponding to the closest rational GB
double Bicrystal::RationalizeOri(int h_a_func, int k_a_func, int l_a_func, double theta_func, double *rot_ax_func, int *CSL_vec){
	cout << "Rationalizing the rotation axis and angle" << endl;
	if( !this->IsCrystalDefined ){
		cerr << "The crystal is not defined, we cannot compute a rational GB" << endl;
		exit(EXIT_FAILURE);
	}
	int hp_near, kp_near, lp_near;
	double sp, norm1, norm2, theta_test;
	vector<double> theta_list;
	vector<int> hkl_list;
	double tolScalarProd = 1e-3;
	double *vec1 = new double[3];
	double *vec2 = new double[3];
	double *vec3 = new double[3];
	double *vec4 = new double[3];
	double *vec5 = new double[3];
	bool DoNotRatAxis = false;
	
	// rationalize the rotation axis by searching the closest (with quasi the same orientation) normal to a crystal plane (with relatively low integer Miller indices)
	int arr[3] = {abs(h_a_func),abs(k_a_func),abs(l_a_func)};
	int MaxHKL_Norm;
	if( MT->max(arr,3) < 3 && ( ( h_a_func == 0 && k_a_func == 0 ) || ( h_a_func == 0 && l_a_func == 0 ) || ( k_a_func == 0 && l_a_func == 0 ) ) ) MaxHKL_Norm = MT->max(arr,3)*5;
	else if( MT->max(arr,3) < 3 ) MaxHKL_Norm = MT->max(arr,3)*100;
	else if( MT->max(arr,3) < 10 ) MaxHKL_Norm = MT->max(arr,3)*50;
	else MaxHKL_Norm = MT->max(arr,3)*10;
	if( ( h_a_func == 0 && k_a_func == 0 ) || ( h_a_func == 0 && l_a_func == 0 ) || ( k_a_func == 0 && l_a_func == 0 ) ) DoNotRatAxis = true;
	if( DoNotRatAxis || this->_MyCrystal->getCrystallo() == "Cubic" ){// for cubic crystals or simple rot axis, the scalar product between plane and direction vector is expressed simply
		hp_near = h_a_func;
		kp_near = k_a_func;
		lp_near = l_a_func;
	}else{
		norm1 = 0.;
		for(unsigned int i=0;i<3;i++){
			vec1[i] = this->_MyCrystal->getA1()[i]*h_a_func + this->_MyCrystal->getA2()[i]*k_a_func + this->_MyCrystal->getA3()[i]*l_a_func;
			norm1 += pow(vec1[i],2.);
		}
		norm1 = sqrt(norm1);
		for(int i=-MaxHKL_Norm;i<MaxHKL_Norm+1;i++){
			for(int j=-MaxHKL_Norm;j<MaxHKL_Norm+1;j++){
				for(int k=-MaxHKL_Norm;k<MaxHKL_Norm+1;k++){
					sp = 0.;
					norm2 = 0.;
					for(unsigned int ii=0;ii<3;ii++){
						vec2[ii] = this->_MyCrystal->getA1_star()[ii]*i + this->_MyCrystal->getA2_star()[ii]*j + this->_MyCrystal->getA3_star()[ii]*k;
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
	norm1 = 0.;
	for(unsigned int i=0;i<3;i++){
		rot_ax_func[i] = this->_MyCrystal->getA1_star()[i]*hp_near + this->_MyCrystal->getA2_star()[i]*kp_near + this->_MyCrystal->getA3_star()[i]*lp_near;
		norm1 += pow(rot_ax_func[i],2.);
	}
	norm1 = sqrt(norm1);
	for(unsigned int i=0;i<3;i++) rot_ax_func[i] /= norm1;
	cout << hp_near << " " << kp_near << " " << lp_near << endl;
	MT->printVec(rot_ax_func);
	// Now hp_near, kp_near and lp_near are plane indices which permits to simplify scalar product and the searching of vector normal to rotation axis	
	// get cell parameter from crystal
	double a_c = this->_MyCrystal->getALength()[0];
	double b_c = this->_MyCrystal->getALength()[1];
	double c_c = this->_MyCrystal->getALength()[2];
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
	// change or not the sign of the second vector to have a direct basis
	for(unsigned int i=0;i<3;i++){
		vec1[i] = this->_MyCrystal->getA1()[i]*Normal_ind[ind_norm*6] + this->_MyCrystal->getA2()[i]*Normal_ind[ind_norm*6+1] + this->_MyCrystal->getA3()[i]*Normal_ind[ind_norm*6+2];  
		vec2[i] = this->_MyCrystal->getA1()[i]*Normal_ind[ind_norm*6+3] + this->_MyCrystal->getA2()[i]*Normal_ind[ind_norm*6+4] + this->_MyCrystal->getA3()[i]*Normal_ind[ind_norm*6+5];
	}
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
				MT->MatDotVec(RotMat,vec4,vec4);
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
					sp = 0.;
					for(unsigned int ii=0;ii<3;ii++){
						vec1[ii] = this->_MyCrystal->getA1_star()[ii]*i + this->_MyCrystal->getA2_star()[ii]*j + this->_MyCrystal->getA3_star()[ii]*k;
						sp += vec1[ii]*rot_ax_func[ii];
					}
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
			cout << "For the rotation axis [" << h_a << k_a << l_a << "] (not the one rationalized) and this angle, the corresponding symmetrical tilt GB have a (" << MilSTGB_near[ind_STGB*3] << "_" << MilSTGB_near[ind_STGB*3+1] << "_" << MilSTGB_near[ind_STGB*3+2] << ") GB plane" << endl;
		}else{
			cout << "We don't have found a corresponding STGB !" << endl;
		}
	}else{
		cout << "We don't have search the corresponding STGB of this misorientation (too long to compute)" << endl;
	}

	delete[] vec1;
	delete[] vec2;
	delete[] int_vec;
	delete[] uint_vec;
	return theta_func+true_DeltaTheta[ind_rattheta];
}

// Based on the CSL basis, we search the size of the GB plane (in x and y dir) for which there is a CSL site
void Bicrystal::searchGBSize(const int h_p_func, const int k_p_func, const int l_p_func){
	if( !this->IsCrystalDefined ){
		cerr << "The crystal is not defined, we cannot compute the GB size" << endl;
		exit(EXIT_FAILURE);
	}
	if( !IsCSL_Basis ){
		cout << "The CSL basis have not been computed, need to perform the calculation !" << endl;
		//searchCSL(this->h_a,this->k_a,this->l_a,this->theta,CSL_vec,0);
	}
	// search two orthogonal vectors of the CSL belonging to the GB plane
	double tolScalarProd = 1e-3;
	int MaxCL_CSL = 100;
	double *vec1 = new double[3];
	double *vec2 = new double[3];
	double sp;
	// get the normal direction to the GB plane
	for(unsigned int i=0;i<3;i++) vec1[i] = h_p_func*_MyCrystal->getA1_star()[0] + k_p_func*_MyCrystal->getA2_star()[i] + l_p_func*_MyCrystal->getA3_star()[i];
	// search a first normal vector
	for(int i=-MaxCL_CSL;i<MaxCL_CSL+1;i++){
		for(int j=-MaxCL_CSL;j<MaxCL_CSL+1;j++){
			for(int k=-MaxCL_CSL;k<MaxCL_CSL+1;k++){
				for(unsigned int ii=0;ii<3;ii++) vec2[ii] = i*this->CSL_Basis[ii*3] + j*this->CSL_Basis[ii*3+1] + k*this->CSL_Basis[ii*3+2];
				sp = 0;
				for(unsigned int ii=0;ii<3;ii++) sp += vec1[ii]*vec2[ii];
				if( fabs(sp) < tolScalarProd && ( i != 0 || j != 0 || k != 0 ) ){

				}
			}
		}
	}
}

void Bicrystal::setOrientedCrystals(const string& crystalName, bool rationalize){
	// TEST
	double *testvec_d = new double[3];
	double a,b,c;
	//a = 4.749;
	//b = 10.1985;
	//c = 5.9792;
	a = 4.7;
	b = 10.;
	c = 6.;
	testvec_d[0] = -0.5814;
	testvec_d[1] = 0.3234;
	testvec_d[2] = 0.7466;
	testvec_d[0] /= a;
	testvec_d[1] /= b;
	testvec_d[2] /= c;
	int *testvec_i = new int[3];
	double tol = 1e-1;
	unsigned int sigma=100000;
	bool found;
	MT->find_integer_vector(testvec_d,tol,sigma,testvec_i,found);
	if(found) cout << "rat vec : " << testvec_i[0] << " " << testvec_i[1] << " " << testvec_i[2] << endl;
	else cout << "no rat vec found" << endl;
	// END TEST
	setCrystal(crystalName);
	this->RotCartToGrain2 = new double[9];
	this->RotGrain1ToGrain2 = new double[9];
	int *CSL_vec = new int[3];
	double *rot_ax = new double[3];
	// construct the first crystal with the wanted plane horyzontal
	this->_MyCrystal->RotateCrystal(h_p,k_p,l_p);
	if( rationalize ){
		double temp_var = RationalizeOri(this->h_a,this->k_a,this->l_a,this->theta,rot_ax,CSL_vec);
		this->theta = temp_var;
		//searchCSL(this->h_a,this->k_a,this->l_a,this->theta,CSL_vec,0);
		//searchCSL(rot_ax,this->theta,CSL_vec,0);
	}
	// compute the rotation to be applied to the second grain
	double NormRotAx = 0.;
	// first rotation in the reference frame of the first crystal
	for(unsigned int i=0;i<3;i++){
		//rot_ax[i] = this->_MyCrystal->getA1()[i]*h_a+this->_MyCrystal->getA2()[i]*k_a+this->_MyCrystal->getA3()[i]*l_a;
		NormRotAx += pow(rot_ax[i],2.);
	}
	NormRotAx = sqrt(NormRotAx);
	for(unsigned int i=0;i<3;i++) rot_ax[i] /= NormRotAx;
	MT->Vec2rotMat(rot_ax,this->theta,this->RotGrain1ToGrain2);
	// add the second rotation to be in the cartesian frame
	MT->MatDotMat(this->RotGrain1ToGrain2,this->_MyCrystal->getRotMat(),this->RotCartToGrain2);
	this->_MyCrystal2 = new Crystal(crystalName);
	this->_MyCrystal2->RotateCrystal(this->RotCartToGrain2);
		searchCSL(rot_ax,this->theta,CSL_vec,0);
	//FOR TEST
	searchGBSize(this->h_p, this->k_p, this->l_p);
	// END TEST
	cout << "constructing grains" << endl;
	cout << "first" << endl;
	this->_MyCrystal->ConstructOrthogonalCell();
	this->_MyCrystal->getOrientedSystem()->print_lmp("Crystal1.lmp");
	cout << "second" << endl;
	this->_MyCrystal2->ConstructOrthogonalCell();
	this->IsCrystal2 = true;
	this->IsRotMatDefine = true;
	delete[] rot_ax;
	delete[] CSL_vec;
}

void Bicrystal::searchGBPos(){
	// Compute bond orientationnal order parameter to know where are the GBs and stored it into aux 1
	unsigned int nbPts_i = 1000;// TODO maybe set a file reading sigma nbPts, lsph and rc for the different systems
	double sigma = 2.;
	setAux(this->CA->BondOrientationalParameter(),"Disorder");
	// Compute density profile of bond ori param
	unsigned int ind_DisoDens;
	ind_DisoDens = Compute1dDensity("Disorder", NormalDir, sigma, nbPts_i);
	// Compute atomic density profile to know if there is vacuum
	unsigned int ind_AtDens;
	ind_AtDens = Compute1dDensity("Atomic", NormalDir, sigma, nbPts_i);
	double VacuumDens = 0.;
	double SecureFac = 1e-2;
	unsigned int count_dens_ave = 0;
	for(unsigned int i=0;i<this->density_nbPts[ind_AtDens];i++){
		if( this->density_prof[ind_AtDens][i*2] > 1e-1 ){
			VacuumDens += this->density_prof[ind_AtDens][i*2];
			count_dens_ave += 1;
		}
	}
	VacuumDens /= count_dens_ave;
	VacuumDens *= 1e-1;

	
	if( NormalDir == "x" ) this->Ldir = this->H1[0];
	else if( NormalDir == "y" ) this->Ldir = this->H2[1];
	else if( NormalDir == "z" ) this->Ldir = this->H3[2];
	
	this->VacuumLo = this->Ldir;
	this->MinPos = this->Ldir;
	this->SystemLength = this->Ldir;
	this->VacuumHi = 0.;
	this->MaxPos = 0.;
	this->IsVacuum = false;
	this->IsCentered = true;
	double buffer;
	vector<double> density_red;
	vector<double> density_red_forfit;
	for(unsigned int i=0;i<this->density_nbPts[ind_AtDens];i++){
		if( this->density_prof[ind_AtDens][i*2] < VacuumDens ){
			this->IsVacuum = true;
			if( this->density_prof[ind_AtDens][i*2+1] > this->VacuumHi ) this->VacuumHi = this->density_prof[ind_AtDens][i*2+1];
			if( this->density_prof[ind_AtDens][i*2+1] < this->VacuumLo ) this->VacuumLo = this->density_prof[ind_AtDens][i*2+1];
		}
	} // TODO non vacuum case ?
	if( this->IsVacuum ){
		if( NormalDir == "x" ){
			for(unsigned int i=0;i<this->nbAtom;i++){
				if( this->WrappedPos[i].x < this->MinPos ) this->MinPos = this->WrappedPos[i].x;
				if( this->WrappedPos[i].x > this->MaxPos ) this->MaxPos = this->WrappedPos[i].x;
			}
		}else if( NormalDir == "y" ){
			for(unsigned int i=0;i<this->nbAtom;i++){
				if( this->WrappedPos[i].y < this->MinPos ) this->MinPos = this->WrappedPos[i].y;
				if( this->WrappedPos[i].y > this->MaxPos ) this->MaxPos = this->WrappedPos[i].y;
			}
		}else if( NormalDir == "z" ){
			for(unsigned int i=0;i<this->nbAtom;i++){
				if( this->WrappedPos[i].z < this->MinPos ) this->MinPos = this->WrappedPos[i].z;
				if( this->WrappedPos[i].z > this->MaxPos ) this->MaxPos = this->WrappedPos[i].z;
			}
		}
		if( fabs(this->MaxPos-this->MinPos-this->Ldir) < 1. ){
			this->MinPos = this->VacuumHi-this->Ldir;
			this->MaxPos = this->VacuumLo;
			this->SystemLength = this->MaxPos-this->MinPos;
			this->IsCentered = false;
		}
		// Search GB position by computing the mean of gaussian distrib in the center of the system
		// compute the max and the mean of the disorder density 
		unsigned int indMaxDiso = 0;
		double MeanDiso = 0.;
		double FracLength = 0.01;
		for(unsigned int i=0;i<this->density_nbPts[ind_DisoDens];i++){
			buffer = this->density_prof[ind_DisoDens][i*2+1];
			if( !this->IsCentered && buffer > this->VacuumHi ) buffer -= this->Ldir;
			// account only for position 60% around the center of system
			//if( buffer > this->MinPos+FracLength*(this->SystemLength) && buffer < this->MaxPos-FracLength*(this->SystemLength) ){
				density_red.push_back(this->density_prof[ind_DisoDens][i*2]);
				density_red.push_back(buffer);
				MeanDiso += this->density_prof[ind_DisoDens][i*2];
				if( density_red[((density_red.size()/2)-1)*2] > density_red[indMaxDiso*2] ) indMaxDiso = (density_red.size()/2)-1;
			//}
		}
		MeanDiso /= (density_red.size()/2.);
		// search the minimal values around the max to define bounds for characterizing the GB
		double rightMin = 0;
		unsigned int indRight = 0;
		unsigned int count_equal;
		unsigned int max_count_equal = 20;
		count_equal = 0;
		for(unsigned int i=indMaxDiso;i<(density_red.size()/2)-1;i++){
			if( density_red[i*2] < density_red[(i+1)*2] ){
				rightMin = density_red[i*2];
				indRight = i;
				break;
			}else if( fabs(density_red[i*2]-density_red[(i+1)*2]) < 1e-8 ) count_equal++;
			if( count_equal >= max_count_equal ){
					rightMin = density_red[i*2];
					indRight = i;
					break;
			}
		}
		double leftMin = 0;
		unsigned int indLeft = 0;
		count_equal = 0;
		for(unsigned int i=indMaxDiso;i>0;i--){
			if( density_red[i*2] < density_red[(i-1)*2] ){
				leftMin = density_red[i*2];
				indLeft = i;
				break;
			}else if( fabs(density_red[i*2]-density_red[(i-1)*2]) < 1e-8 ) count_equal++;
			if( count_equal >= max_count_equal ){
					leftMin = density_red[i*2];
					indLeft = i;
					break;
			}
		}
		// if one of the two is higher than the mean we are in local minimum inside the GB
		// => search the following minimum
		double facMean = 2.3;
		unsigned int max_min_search = 250;
		unsigned int count_min_search = 0;
		bool maxfound;
		cout << "Mean diso : " << MeanDiso << endl;
		while( rightMin > MeanDiso*facMean && count_min_search < max_min_search ){
			maxfound = false;
			count_equal = 0;
			for(unsigned int i=indRight;i<(density_red.size()/2)-1;i++){
				if( !maxfound && density_red[i*2] > density_red[(i+1)*2] ) maxfound = true;
				if( maxfound && density_red[i*2] < density_red[(i+1)*2] ){
					rightMin = density_red[i*2];
					indRight = i;
					break;
				}else if( fabs(density_red[i*2]-density_red[(i+1)*2]) < 1e-8 ) count_equal++;
				if( count_equal >= max_count_equal ){
					rightMin = density_red[i*2];
					indRight = i;
					break;
				}
			}
			count_min_search++;
		}
		cout << "right min : " << rightMin << " " << count_min_search << endl;
		count_min_search = 0;
		while( leftMin > MeanDiso*facMean && count_min_search < max_min_search ){
			maxfound = false;
			count_equal = 0;
			for(unsigned int i=indLeft;i>0;i--){
				if( !maxfound && density_red[i*2] > density_red[(i-1)*2] ) maxfound = true;
				if( maxfound && density_red[i*2] < density_red[(i-1)*2] ){
					leftMin = density_red[i*2];
					indLeft = i;
					break;
				}else if( fabs(density_red[i*2]-density_red[(i-1)*2]) < 1e-8 ) count_equal++;
				if( count_equal >= max_count_equal ){
					leftMin = density_red[i*2];
					indLeft = i;
					break;
				}
			}
			count_min_search++;
		}		
		cout << "left min : " << leftMin << " " << count_min_search << endl;
		for(unsigned int i=indLeft;i<indRight+1;i++){
			density_red_forfit.push_back(density_red[i*2]);
			density_red_forfit.push_back(density_red[i*2+1]);
		}
		double MinDiso;
		if( leftMin > rightMin*facMean ) MinDiso = rightMin;
		else MinDiso = leftMin;

		// second loop to shift the distrib of order param and compute the normalization factor
		double NormFacGauss = 0;
		for(unsigned int i=0;i<(density_red_forfit.size()/2);i++){
			density_red_forfit[i*2] -= MinDiso;
			NormFacGauss += density_red_forfit[i*2];
		}
		for(unsigned int i=0;i<density_red.size()/2;i++){
			density_red[i*2] -= MinDiso;
			density_red[i*2] /= NormFacGauss;
		}
		// normalize the distribution
		for(unsigned int i=0;i<(density_red_forfit.size()/2);i++){
			density_red_forfit[i*2] /= NormFacGauss;
		}

		// set the GB profile of AtomicSystem class
		this->density_prof.push_back(new double[density_red.size()]);
		for(unsigned int i=0;i<density_red.size();i++) this->density_prof[density_prof.size()-1][i] = density_red[i];
		this->density_name.push_back(new string[2]);
		this->density_name[density_name.size()-1][0] = "GBProfile";
		this->density_name[density_name.size()-1][1] = this->NormalDir;
		this->density_nbPts.push_back(density_red.size()/2);

		// estimate the GB position 
		this->GBPos1 = 0;
		for(unsigned int i=0;i<(density_red_forfit.size()/2);i++) this->GBPos1 += density_red_forfit[i*2]*density_red_forfit[i*2+1];

		// estimate the GB width 
		this->GBwidth1 = 0;
		for(unsigned int i=0;i<(density_red_forfit.size()/2);i++) this->GBwidth1 += density_red_forfit[i*2]*pow(density_red_forfit[i*2+1]-this->GBPos1,2.);
		this->GBwidth1 = sqrt(this->GBwidth1);
	} // end if IsVacuum
	double sigma_fit, mu_fit, prefac;
	sigma_fit = this->GBwidth1;
	mu_fit = this->GBPos1;
	prefac = 1.;
	// fit a gaussian to rafine the estimation of GB width and position
	MT->gaussian_fit(density_red_forfit, mu_fit, sigma_fit, prefac);
	this->GBwidth1 = 2.355*sigma_fit;
	this->GBPos1 = mu_fit;
	// set the gaussian GB profile of AtomicSystem class (without prefac for the distribution to be normalized and for different system to be compared to each others)
	this->density_prof.push_back(new double[density_red.size()]);
	for(unsigned int i=0;i<density_red.size()/2;i++){
		//this->density_prof[density_prof.size()-1][i*2] = MT->gaussian_prefac(density_red[i*2+1], mu_fit, sigma_fit, prefac);
		this->density_prof[density_prof.size()-1][i*2] = MT->gaussian(density_red[i*2+1], mu_fit, sigma_fit);
		this->density_prof[density_prof.size()-1][i*2+1] = density_red[i*2+1];
	}
	this->density_name.push_back(new string[2]);
	this->density_name[density_name.size()-1][0] = "GBProfile_Gauss";
	this->density_name[density_name.size()-1][1] = this->NormalDir;
	this->density_nbPts.push_back(density_red.size()/2);
}

void Bicrystal::ComputeExcessVolume(){
	// Compute mass density along GB normal direction
	double sigma = 2.; // TODO once again to see if a file can be define for those vars or static vars.. ?
	unsigned int nbPts = 500;
	unsigned int index;
	index = this->Compute1dDensity("Mass", this->NormalDir, sigma, nbPts);
	double PC_dens = 0.;
	unsigned int count = 0;
	double buffer;
	this->ExcessVol = 0.;
	// compute density of perfect crystal assuming that there is perfect crystal between MinPos+GBwidth and GBPos-GBwidth and same above GBPos
	if( IsVacuum ){
		vector<double> dens_GB;
		for(unsigned int i=0;i<this->density_nbPts[index];i++){
			buffer = this->density_prof[index][i*2+1];
			if( !this->IsCentered && buffer > this->VacuumHi ) buffer -= this->Ldir;
			if( ( buffer > MinPos+GBwidth1 && buffer < GBPos1-GBwidth1 ) || ( buffer < MaxPos-GBwidth1 && buffer > GBPos1+GBwidth1 ) ){
				PC_dens += this->density_prof[index][i*2];
				count += 1;
			}else if( buffer > GBPos1-GBwidth1 && buffer < GBPos1+GBwidth1 ){
				dens_GB.push_back(density_prof[index][i*2]);
				dens_GB.push_back(buffer);
			}
		}
		PC_dens /= count;
		// compute excess volume by integrating density in the GB region
		for(unsigned int i=0;i<(dens_GB.size()/2)-2;i++) this->ExcessVol += ((2.*PC_dens/(dens_GB[i*2]+dens_GB[(i+1)*2]))-1.)*fabs((dens_GB[(i+1)*2+1]-dens_GB[i*2+1]));

	}
}
void Bicrystal::print_Grains(){
	if( this->AreGrainsDefined ){
		this->Grain1->print_lmp("Grain1.lmp");
		this->Grain2->print_lmp("Grain2.lmp");
		cout << "The two grains have been printed (\"Grain1.lmp\" and \"Grain2.lmp\")" << endl;
	}else{
		cout << "The two grains are not defined, we cannot print them" << endl;
	}
}

void Bicrystal::read_params(){
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+FixedParam_Filename;
	ifstream file(filename, ios::in);
	size_t pos_thetamax, pos_MaxHKL, pos_toldist, pos_sigmaMax, pos_tolpos_kC, pos_tolCSLint, pos_tolAlign, pos_rcut, pos_lsph;
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
			pos_sigmaMax=line.find("SIGMA_MAX ");
			if(pos_sigmaMax!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->SigmaMax;
			}
			pos_tolpos_kC=line.find("TOL_POS_KNOWN_CSL_VEC ");
			if(pos_tolpos_kC!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->tolpos_known_CSL;
			}
			pos_tolCSLint=line.find("TOL_CSL_INTEGER ");
			if(pos_tolCSLint!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->tol_CSL_integer;
			}
			pos_tolAlign=line.find("TOL_ALIGNMENT_CSL ");
			if(pos_tolAlign!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->tolAlignment_CSL;
			}
		}
	}else{
		cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
		exit(EXIT_FAILURE);
	}
}



Bicrystal::~Bicrystal(){
	if( this->IsCrystal2 ) delete _MyCrystal2;
	if( this->IsCA ) delete CA;
	if( this->AreGrainsDefined ){
		delete Grain1;
		delete Grain2;
		delete[] AtomList_G1;
		delete[] AtomList_G2;
		delete[] H1_G1;
		delete[] H2_G1;
		delete[] H3_G1;
		delete[] H1_G2;
		delete[] H2_G2;
		delete[] H3_G2;
	}
	if( IsRotMatDefine ){
		delete[] RotGrain1ToGrain2;
		delete[] RotCartToGrain2;
	}
	if( IsCSL ) for(unsigned int i=0;i<CSL.size();i++) delete[] CSL[i];
	if( IsCSL_Basis ) delete[] CSL_Basis;
}
