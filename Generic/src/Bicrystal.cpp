#include "AtomHicConfig.h"
#include <iostream>
#include <vector>
#include "Bicrystal.h"

using namespace std;

Bicrystal::Bicrystal(const string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p){
	setCrystal(crystalName);
}

Bicrystal::Bicrystal(const string& filename, const string NormalDir, const string CrystalName):AtomicSystem(filename){
	setCrystal(CrystalName);
	this->NormalDir = NormalDir;
	CA = new ComputeAuxiliary(this);
	searchGBPos();
}

Bicrystal::Bicrystal(const string& filename, const string NormalDir):AtomicSystem(filename){
	this->NormalDir = NormalDir;
	CA = new ComputeAuxiliary(this);
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
	Compute1dDensity("Disorder", NormalDir, sigma, nbPts_i);
	// Compute atomic density profile to know if there is vacuum
	Compute1dDensity("Atomic", NormalDir, sigma, nbPts_i);
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
	for(unsigned int i=0;i<this->density_nbPts[1];i++){
		if( this->density_prof[1][i*2] < VacuumDens ){
			this->IsVacuum = true;
			if( this->density_prof[1][i*2+1] > this->VacuumHi ) this->VacuumHi = this->density_prof[1][i*2+1];
			if( this->density_prof[1][i*2+1] < this->VacuumLo ) this->VacuumLo = this->density_prof[1][i*2+1];
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
		// first loop to search the min order param
		double MinDiso = this->density_prof[0][0];
		for(unsigned int i=0;i<this->density_nbPts[0];i++){
			buffer = this->density_prof[0][i*2+1];
			if( !IsCentered && buffer > this->VacuumHi ) buffer -= this->Ldir;
			// account only for position 60% around the center of system
			if( buffer > this->MinPos+0.2*(this->SystemLength) && buffer < this->MaxPos-0.2*(this->SystemLength) ){
				density_red.push_back(this->density_prof[0][i*2]);
				density_red.push_back(buffer);
				if( this->density_prof[0][i*2] < MinDiso ) MinDiso = this->density_prof[0][i*2];
			}
		}
		// second loop to shift the distrib of order param and compute the normalization factor
		double NormFacGauss = 0;
		for(unsigned int i=0;i<(density_red.size()/2)-1;i++){
			density_red[i*2] -= MinDiso;
			NormFacGauss += density_red[i*2];
		}
		// normalize the distribution
		for(unsigned int i=0;i<(density_red.size()/2);i++) density_red[i*2] /= NormFacGauss;

		// set the GB profile of AtomicSystem class
		this->density_prof.push_back(new double[density_red.size()]);
		for(unsigned int i=0;i<density_red.size();i++) this->density_prof[density_prof.size()-1][i] = density_red[i];
		this->density_name.push_back(new string[2]);
		this->density_name[density_name.size()-1][0] = "GBProfile";
		this->density_name[density_name.size()-1][1] = this->NormalDir;
		this->density_nbPts.push_back(density_red.size()/2);

		// estimate the GB position 
		this->GBPos1 = 0;
		for(unsigned int i=0;i<(density_red.size()/2);i++) this->GBPos1 += density_red[i*2]*density_red[i*2+1];

		// estimate the GB width 
		this->GBwidth1 = 0;
		for(unsigned int i=0;i<(density_red.size()/2);i++) this->GBwidth1 += density_red[i*2]*pow(density_red[i*2+1]-this->GBPos1,2.);
		this->GBwidth1 = sqrt(this->GBwidth1);
	}
	double sigma_fit, mu_fit, prefac;
	sigma_fit = this->GBwidth1;
	mu_fit = this->GBPos1;
	prefac = 1.;
	// fit a gaussian to rafine the estimation of GB width and position
	MT->gaussian_fit(density_red, mu_fit, sigma_fit, prefac);
	this->GBwidth1 = sigma_fit;
	this->GBPos1 = mu_fit;
	// set the gaussian GB profile of AtomicSystem class
	this->density_prof.push_back(new double[density_red.size()]);
	for(unsigned int i=0;i<density_red.size()/2;i++){
		this->density_prof[density_prof.size()-1][i*2] = MT->gaussian_prefac(density_red[i*2+1], mu_fit, sigma_fit, prefac);
		this->density_prof[density_prof.size()-1][i*2+1] = density_red[i*2+1];
	}
	this->density_name.push_back(new string[2]);
	this->density_name[density_name.size()-1][0] = "GBProfile_Gauss";
	this->density_name[density_name.size()-1][1] = this->NormalDir;
	this->density_nbPts.push_back(density_red.size()/2);
}

void Bicrystal::ComputeExcessVolume(){
	cout << "Compute EV" << endl;
}

Bicrystal::~Bicrystal(){
	delete CA;
}
