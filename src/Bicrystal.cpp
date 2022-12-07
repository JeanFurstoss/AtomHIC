#include "AtomHicConfig.h"
#include <iostream>
#include <vector>
#include "Bicrystal.h"

using namespace std;

// Constructor for bicrystal with facetted GB with given misorientation and GB plane and facet type (for the moment this is implemented for 2D facets and for a GB composed 2 different facet types, i.e. with one facet direction parallel to x direction)
Bicrystal::Bicrystal(const string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p, vector<int> FacetsType, unsigned int N_facet){
	this->MT = new MathTools;
	setOrientedCrystals(crystalName, h_a, k_a, l_a, theta, h_p, k_p, l_p);
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
	double xl1, xl2, yl1, yl2;
	unsigned int dupX1, dupX2, dupY1, dupY2;
	xl1 = this->_MyCrystal->getOrientedSystem()->getH1()[0];
	xl2 = this->_MyCrystal2->getOrientedSystem()->getH1()[0];
	yl1 = this->_MyCrystal->getOrientedSystem()->getH2()[1];
	yl2 = this->_MyCrystal2->getOrientedSystem()->getH2()[1];
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
    	double Mx1 = final_xbox / (xl1*dupX1);
    	double My1 = final_ybox / (yl1*dupY1);
    	double Mx2 = final_xbox / (xl2*dupX2);
    	double My2 = final_ybox / (yl2*dupY2);
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
	this->AtomType = new string[this->nbAtomType];
	this->AtomType_uint = new unsigned int[this->nbAtomType];
	this->AtomMass = new double[this->nbAtomType];
	for(unsigned int i=0;i<this->nbAtomType;i++){
		this->AtomType[i] = _MyCrystal->getAtomType(i);
		this->AtomType_uint[i] = _MyCrystal->getAtomType_uint(i);
		this->AtomMass[i] = _MyCrystal->getAtomMass(i);
	}
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
	this->Grain1 = new AtomicSystem(this->AtomList_G1,nbAtom1,this->H1_G1,this->H2_G1,this->H3_G1);
	this->Grain2 = new AtomicSystem(this->AtomList_G2,nbAtom2,this->H1_G2,this->H2_G2,this->H3_G2);
	this->AreGrainsDefined = true;
	// verify the stoichiometry
	unsigned int *currentStoich = new unsigned int[this->_MyCrystal->getNbAtomType()];
	for(unsigned int i=0;i<this->_MyCrystal->getNbAtomType();i++) currentStoich[i] = 0;
	for(unsigned int i=0;i<at_count;i++){
		for(unsigned int t=0;t<this->_MyCrystal->getNbAtomType();t++){
			if( AtomList_temp[i].type_uint == this->_MyCrystal->getAtomType_uint(t) ){
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
				cout << "number of " << this->_MyCrystal->getAtomType(t) << " : " << currentStoich[t] << ", ";
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
				Type2Rm.push_back(this->_MyCrystal->getAtomType_uint(i));
				Type2Rm.push_back(i);
				Neigh.push_back(vector<double>());
			}
		}
		double rc_squared = pow(3.,2.);
		double d_squared;
		for(unsigned int i=0;i<at_count;i++){
			for(unsigned int t=0;t<Type2Rm.size()/2;t++){
				if( AtomList_temp[i].pos.z < GBup && AtomList_temp[i].pos.z > GBdown && AtomList_temp[i].type_uint == Type2Rm[t*2] ){
					xpos = AtomList_temp[i].pos.x;
					ypos = AtomList_temp[i].pos.y;
					zpos = AtomList_temp[i].pos.z;
					for(unsigned int j=0;j<at_count;j++){
						if( i != j && AtomList_temp[j].pos.z < GBup && AtomList_temp[j].pos.z > GBdown  && AtomList_temp[j].type_uint == Type2Rm[t*2] ){
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
		for(unsigned int t=0;t<Type2Rm.size()/2;t++){
			this->MT->sort(Neigh[t],0,2,Neigh[t]);
			this->nbAtom += adjust[Type2Rm[t*2+1]];
		}
		this->AtomList = new Atom[this->nbAtom];
		unsigned int at_count2 = 0;
		bool store;
		for(unsigned int i=0;i<at_count;i++){
			store = true;
			for(unsigned int t=0;t<Type2Rm.size()/2;t++){
				if( AtomList_temp[i].type_uint == Type2Rm[t*2] ){
					for(unsigned int n=0;n<abs(adjust[Type2Rm[t*2+1]]);n++){
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
Bicrystal::Bicrystal(const string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p){
	double GBspace = 1.; //TODO maybe define elsewhere (in a file in the main prog for example..)
	double MaxMisfit = 0.02;
	unsigned int MaxDup = 50;
	bool FullGrains = true; // TODO as the variables above
	this->MT = new MathTools;
	setOrientedCrystals(crystalName, h_a, k_a, l_a, theta, h_p, k_p, l_p);
	// search the number of linear combination for the two system to have the same x y length
	double xl1, xl2, yl1, yl2;
	unsigned int dupX1, dupX2, dupY1, dupY2;
	xl1 = this->_MyCrystal->getOrientedSystem()->getH1()[0];
	xl2 = this->_MyCrystal2->getOrientedSystem()->getH1()[0];
	yl1 = this->_MyCrystal->getOrientedSystem()->getH2()[1];
	yl2 = this->_MyCrystal2->getOrientedSystem()->getH2()[1];
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
    	double Mx1 = final_xbox / (xl1*dupX1);
    	double My1 = final_ybox / (yl1*dupY1);
    	double Mx2 = final_xbox / (xl2*dupX2);
    	double My2 = final_ybox / (yl2*dupY2);
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
	this->AtomType = new string[this->nbAtomType];
	this->AtomType_uint = new unsigned int[this->nbAtomType];
	this->AtomMass = new double[this->nbAtomType];
	for(unsigned int i=0;i<this->nbAtomType;i++){
		this->AtomType[i] = _MyCrystal->getAtomType(i);
		this->AtomType_uint[i] = _MyCrystal->getAtomType_uint(i);
		this->AtomMass[i] = _MyCrystal->getAtomMass(i);
	}
	this->IsCharge = _MyCrystal->getIsCharge();
	this->IsTilted = false;
	computeInverseCellVec();
	// paste the two grains
	for(unsigned int i=0;i<dupX1;i++){
		for(unsigned int j=0;j<dupY1;j++){
			for(unsigned int n=0;n<nbAtom1;n++){
				this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n] = this->_MyCrystal->getOrientedSystem()->getAtom(n);
				this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.x = Mx1*(this->AtomList[n].pos.x + i*xl1);
				this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.y = My1*(this->AtomList[n].pos.y + j*yl1);
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
				if( !FullGrains && i == 0 && j == 0 ) this->AtomList_G2[n] = this->AtomList[nbAtom1*dupX1*dupY1+n];
				else this->AtomList_G2[i*dupY2*nbAtom2+j*nbAtom2+n] = this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n];
				this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.z += (this->_MyCrystal->getOrientedSystem()->getH3()[2]+GBspace);
			}
		}
	}
	// construct the two grains
	this->Grain1 = new AtomicSystem(this->AtomList_G1,nbAtom1_G,this->H1_G1,this->H2_G1,this->H3_G1);
	this->Grain2 = new AtomicSystem(this->AtomList_G2,nbAtom2_G,this->H1_G2,this->H2_G2,this->H3_G2);
	this->AreGrainsDefined = true;
	string h_a_str = to_string(h_a);
	string k_a_str = to_string(k_a);
	string l_a_str = to_string(l_a);
	string theta_str = to_string(theta*180/M_PI);
	string h_p_str = to_string(h_p);
	string k_p_str = to_string(k_p);
	string l_p_str = to_string(l_p);
	this->File_Heading = " # ["+h_a_str+k_a_str+l_a_str+"]"+theta_str+"°("+h_p_str+k_p_str+l_p_str+") "+crystalName+" grain boundary\n";
	this->Grain1->set_File_Heading(" # Lower grain of the ["+h_a_str+k_a_str+l_a_str+"]"+theta_str+"°("+h_p_str+k_p_str+l_p_str+") "+crystalName+" grain boundary\n");
	this->Grain2->set_File_Heading(" # Upper grain of the ["+h_a_str+k_a_str+l_a_str+"]"+theta_str+"°("+h_p_str+k_p_str+l_p_str+") "+crystalName+" grain boundary\n");
}

Bicrystal::Bicrystal(const string& filename, const string NormalDir, const string CrystalName):AtomicSystem(filename){
	if( NormalDir != "z" || NormalDir != "y" || NormalDir != "x" || NormalDir != "Z" || NormalDir != "Y" || NormalDir != "X" ){
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

Bicrystal::Bicrystal(const string& filename, const string NormalDir):AtomicSystem(filename){
	if( NormalDir != "z" || NormalDir != "y" || NormalDir != "x" || NormalDir != "Z" || NormalDir != "Y" || NormalDir != "X" ){
		cerr << "The direction normal to the GB have to be : x,X,y,Y,z or Z" << endl;
		exit(EXIT_FAILURE);
	}
	this->NormalDir = NormalDir;
	this->IsMassDensity = false;
	this->CA = new ComputeAuxiliary(this);
	this->IsCA = true;
	searchGBPos();
}

void Bicrystal::setOrientedCrystals(const string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p){
	setCrystal(crystalName);
	// construct the first crystal with the wanted plane horyzontal
	this->_MyCrystal->ConstructOrientedSystem(h_p,k_p,l_p);
	// compute the rotation to be applied to the second grain
	double *rot_ax = new double[3];
	double NormRotAx = 0.;
	// first rotation in the reference frame of the first crystal
	for(unsigned int i=0;i<3;i++){
		rot_ax[i] = this->_MyCrystal->getA1()[i]*h_a+this->_MyCrystal->getA2()[i]*k_a+this->_MyCrystal->getA3()[i]*l_a;
		NormRotAx += pow(rot_ax[i],2.);
	}
	NormRotAx = sqrt(NormRotAx);
	for(unsigned int i=0;i<3;i++) rot_ax[i] /= NormRotAx;
	this->RotCartToGrain2 = new double[9];
	this->RotGrain1ToGrain2 = new double[9];
	MT->Vec2rotMat(rot_ax,theta,this->RotGrain1ToGrain2);
	// add the second rotation to be in the cartesian frame
	MT->MatDotMat(this->_MyCrystal->getRotMat(),this->RotGrain1ToGrain2,this->RotCartToGrain2);
	this->_MyCrystal2 = new Crystal(crystalName);
	this->_MyCrystal2->ConstructOrientedSystem(this->RotCartToGrain2);
	this->IsCrystal2 = true;
	this->IsRotMatDefine = true;
	delete[] rot_ax;
}

void Bicrystal::searchGBPos(){
	// Compute bond orientationnal order parameter to know where are the GBs and stored it into aux 1
	int lsph = 20; // TODO maybe set a file reading sigma nbPts, lsph and rc for the different systems
	double rc = 5.;
	unsigned int nbPts_i = 1000;
	double sigma = 2.;
	setAux(this->CA->BondOrientationalParameter(lsph, rc),"Disorder");
	// Compute density profile of bond ori param
	unsigned int ind_DisoDens;
	ind_DisoDens = Compute1dDensity("Disorder", NormalDir, sigma, nbPts_i);
	// Compute atomic density profile to know if there is vacuum
	unsigned int ind_AtDens;
	ind_AtDens = Compute1dDensity("Atomic", NormalDir, sigma, nbPts_i);
	double VacuumDens = 1e-2;
	
	if( NormalDir == "x" ) this->Ldir = this->H1[0];
	else if( NormalDir == "y" ) this->Ldir = this->H2[1];
	else if( NormalDir == "z" ) this->Ldir = this->H3[2];
	
	this->VacuumLo = this->Ldir;
	this->MinPos = this->Ldir;
	this->SystemLength = this->Ldir;
	this->VacuumHi = 0.;
	this->MaxPos = 0.;
	this->IsVacuum = false;
	this->IsVacuum = true;
	double buffer;
	vector<double> density_red;
	vector<double> density_red_forfit;
	for(unsigned int i=0;i<this->density_nbPts[ind_AtDens];i++){
		if( this->density_prof[ind_AtDens][i*2] < VacuumDens ){
			this->IsVacuum = true;
			if( this->density_prof[ind_AtDens][i*2+1] > this->VacuumHi ) this->VacuumHi = this->density_prof[ind_AtDens][i*2+1];
			if( this->density_prof[ind_AtDens][i*2+1] < this->VacuumLo ) this->VacuumLo = this->density_prof[ind_AtDens][i*2+1];
		}
	}
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
		if( fabs(this->MinPos-(this->VacuumLo)) > 1. ){
			this->MinPos = this->VacuumHi-(this->Ldir);
			this->MaxPos = this->VacuumLo;
			this->SystemLength = this->MaxPos-(this->MinPos);
			this->IsCentered = false;
		}
		// Search GB position by computing the mean of gaussian distrib in the center of the system
		// compute the max and the mean of the disorder density 
		unsigned int indMaxDiso = 0;
		double MeanDiso = 0.;
		for(unsigned int i=0;i<this->density_nbPts[ind_DisoDens];i++){
			buffer = this->density_prof[ind_DisoDens][i*2+1];
			if( !this->IsCentered && buffer > this->VacuumHi ) buffer -= this->Ldir;
			// account only for position 60% around the center of system
			if( buffer > this->MinPos+0.2*(this->SystemLength) && buffer < this->MaxPos-0.2*(this->SystemLength) ){
				density_red.push_back(this->density_prof[ind_DisoDens][i*2]);
				density_red.push_back(buffer);
				MeanDiso += this->density_prof[ind_DisoDens][i*2];
				if( density_red[((density_red.size()/2)-1)*2] > density_red[indMaxDiso*2] ) indMaxDiso = (density_red.size()/2)-1;
			}
		}
		MeanDiso /= (density_red.size()/2.);
		// search the minimal values around the max to define bounds for characterizing the GB
		double rightMin = 0;
		unsigned int indRight = 0;
		for(unsigned int i=indMaxDiso;i<(density_red.size()/2)-1;i++){
			if( density_red[i*2] < density_red[(i+1)*2] ){
				rightMin = density_red[i*2];
				indRight = i;
				break;
			}
		}
		double leftMin = 0;
		unsigned int indLeft = 0;
		for(unsigned int i=indMaxDiso;i>=0;i--){
			if( density_red[i*2] < density_red[(i-1)*2] ){
				leftMin = density_red[i*2];
				indLeft = i;
				break;
			}
		}
		// if one of the two is higher than the mean we are in local minimum inside the GB
		// => search the second minimum
		cout << MeanDiso << endl;
		cout << leftMin << " " << rightMin << endl;
		if( rightMin > MeanDiso ){
			bool maxfound = false;
			for(unsigned int i=indRight;i<(density_red.size()/2)-1;i++){
				if( !maxfound && density_red[i*2] > density_red[(i+1)*2] ) maxfound = true;
				if( maxfound && density_red[i*2] < density_red[(i+1)*2] ){
					rightMin = density_red[i*2];
					indRight = i;
					break;
				}
			}
		}
		double facMean = 1.3;
		if( leftMin > MeanDiso*facMean ){
			bool maxfound = false;
			for(unsigned int i=indLeft;i>=0;i--){
				if( !maxfound && density_red[i*2] > density_red[(i-1)*2] ) maxfound = true;
				if( maxfound && density_red[i*2] < density_red[(i-1)*2] ){
					leftMin = density_red[i*2];
					indLeft = i;
					break;
				}
			}
		}		
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
	}
	double sigma_fit, mu_fit, prefac;
	sigma_fit = this->GBwidth1;
	mu_fit = this->GBPos1;
	prefac = 1.;
	// fit a gaussian to rafine the estimation of GB width and position
	MT->gaussian_fit(density_red_forfit, mu_fit, sigma_fit, prefac);
	this->GBwidth1 = 2.*sigma_fit;
	this->GBPos1 = mu_fit;
	// set the gaussian GB profile of AtomicSystem class
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
	// compute density of perfect crystal assuming that there is perfect crystal between MinPos+((GBPos-2*GBwidth)-MinPos)*0.5 and GBPos-2*GBwidth and same above GBPos
	if( IsVacuum ){
		vector<double> dens_GB;
		for(unsigned int i=0;i<this->density_nbPts[index];i++){
			buffer = this->density_prof[index][i*2+1];
			if( !this->IsCentered && buffer > this->VacuumHi ) buffer -= this->Ldir;
			if( ( buffer > .5*(MinPos+GBPos1-2.*GBwidth1) && buffer < GBPos1-2.*GBwidth1 ) || ( buffer < .5*(MaxPos+GBPos1)+GBwidth1 && buffer > GBPos1+2.*GBwidth1 ) ){
				PC_dens += this->density_prof[index][i*2];
				count += 1;
			}else if( buffer > GBPos1-2*GBwidth1 && buffer < GBPos1+2*GBwidth1 ){
				dens_GB.push_back(density_prof[index][i*2]);
				dens_GB.push_back(buffer);
			}
		}
		PC_dens /= count;
		// compute excess volume by integrating density in the GB region
		for(unsigned int i=0;i<(dens_GB.size()/2)-1;i++) this->ExcessVol += ((2.*PC_dens/(dens_GB[i*2]+dens_GB[(i+1)*2]))-1.)*fabs((dens_GB[(i+1)*2+1]-dens_GB[i*2+1]));

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
}
