#include "AtomHicConfig.h"
#include <iostream>
#include <vector>
#include "Bicrystal.h"

using namespace std;

Bicrystal::Bicrystal(const string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p){
	double GBspace = 1.; //TODO maybe define elsewhere (in a file in the main prog for example..)
	setCrystal(crystalName);
	this->MT = new MathTools;
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
	double *rot_mat = new double[9];
	MT->Vec2rotMat(rot_ax,theta,rot_mat);
	// add the second rotation to be in the cartesian frame
	MT->MatDotMat(this->_MyCrystal->getRotMat(),rot_mat,rot_mat);
	this->_MyCrystal2 = new Crystal(crystalName);
	this->_MyCrystal2->ConstructOrientedSystem(rot_mat);
	this->IsCrystal2 = true;
	// search the number of linear combination for the two system to have the same x y length
	double xl1, xl2, yl1, yl2;
	double MaxMisfit = 0.02;
	unsigned int MaxDup = 50;
	unsigned int dupX1, dupX2, dupY1, dupY2;
	xl1 = this->_MyCrystal->getOrientedSystem()->getH1()[0];
	xl2 = this->_MyCrystal2->getOrientedSystem()->getH1()[0];
	yl1 = this->_MyCrystal->getOrientedSystem()->getH2()[1];
	yl2 = this->_MyCrystal2->getOrientedSystem()->getH2()[1];
	this->_MyCrystal2->getOrientedSystem()->print_lmp("test.lmp");
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
	this->nbAtom = nbAtom1*dupX1*dupY1 + nbAtom2*dupX2*dupY2;
	this->AtomList = new Atom[this->nbAtom];
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
	this->H3[2] = this->_MyCrystal->getOrientedSystem()->getH3()[2] + this->_MyCrystal2->getOrientedSystem()->getH3()[2]+GBspace;
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
			}
		}
	}
	for(unsigned int i=0;i<dupX2;i++){
		for(unsigned int j=0;j<dupY2;j++){
			for(unsigned int n=0;n<nbAtom2;n++){
				this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n] = this->_MyCrystal2->getOrientedSystem()->getAtom(n);
				this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.x = Mx2*(this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.x + i*xl2);
				this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.y = My2*(this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.y + j*yl2);
				this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.z += (this->_MyCrystal->getOrientedSystem()->getH3()[2]+GBspace);
			}
		}
	}

	delete[] rot_ax;
	delete[] rot_mat;
}

Bicrystal::Bicrystal(const string& filename, const string NormalDir, const string CrystalName):AtomicSystem(filename){
	setCrystal(CrystalName);
	this->NormalDir = NormalDir;
	this->IsMassDensity = false;
	CA = new ComputeAuxiliary(this);
	searchGBPos();
	ComputeExcessVolume();
}

Bicrystal::Bicrystal(const string& filename, const string NormalDir):AtomicSystem(filename){
	this->NormalDir = NormalDir;
	this->IsMassDensity = false;
	CA = new ComputeAuxiliary(this);
	this->IsCA = true;
	searchGBPos();
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

Bicrystal::~Bicrystal(){
	if( IsCrystal2 ) delete _MyCrystal2;
	if( IsCA ) delete CA;
}
