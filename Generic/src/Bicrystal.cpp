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
	// Compute bond orientationnal order parameter to know where are the GBs and stored it into aux 1
	int lsph = 20; // TODO maybe set a file reading sigma nbPts, lsph and rc for the different systems
	double rc = 5.;
	setAux(this->CA->BondOrientationalParameter(lsph, rc),"Disorder");
	// Compute density profile of bond ori param
	Compute1dDensity("Disorder", NormalDir, 2., 100);
	// Compute atomic density profile to know if there is vacuum
	Compute1dDensity("Atomic", NormalDir, 2., 100);
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
		double MinDiso = 1;
		for(unsigned int i=0;i<this->density_nbPts[0];i++){
			buffer = this->density_prof[0][i*2+1];
			if( !IsCentered && buffer > this->VacuumHi ) buffer -= this->Ldir;
			// account only for position 60% around the center of system
			if( buffer > this->MinPos+0.2*(this->SystemLength) && buffer < this->MaxPos-0.2*(this->SystemLength) ){
				if( this->density_prof[0][i*2] < MinDiso ) MinDiso = this->density_prof[0][i*2];
			}
		}
		// second loop to normalize and shift the distrib of order param
		double NormFacGauss = 0;
		for(unsigned int i=0;i<this->density_nbPts[0];i++){
			buffer = this->density_prof[0][i*2+1];
			if( !IsCentered && buffer > this->VacuumHi ) buffer -= this->Ldir;
			// account only for position 60% around the center of system
			if( buffer > this->MinPos+0.2*(this->SystemLength) && buffer < this->MaxPos-0.2*(this->SystemLength) ) NormFacGauss += (this->density_prof[0][i*2]-MinDiso);
		}
		// third loop to compute the gaussian mean (GBPos)
		this->GBPos1 = 0;
		for(unsigned int i=0;i<this->density_nbPts[0];i++){
			buffer = this->density_prof[0][i*2+1];
			if( !IsCentered && buffer > this->VacuumHi ) buffer -= this->Ldir;
			// account only for position 60% around the center of system
			if( buffer > this->MinPos+0.2*(this->SystemLength) && buffer < this->MaxPos-0.2*(this->SystemLength) ) this->GBPos1 += buffer*(this->density_prof[0][i*2]-MinDiso)/NormFacGauss;
		}
		// fourth loop to compute the gaussian stdev (GBwidth)
		this->GBwidth1 = 0;
		for(unsigned int i=0;i<this->density_nbPts[0];i++){
			buffer = this->density_prof[0][i*2+1];
			if( !IsCentered && buffer > this->VacuumHi ) buffer -= this->Ldir;
			// account only for position 60% around the center of system
			if( buffer > this->MinPos+0.2*(this->SystemLength) && buffer < this->MaxPos-0.2*(this->SystemLength) ) this->GBwidth1 += pow((buffer*(this->density_prof[0][i*2]-MinDiso)/NormFacGauss)-GBPos1,2.);
		}
		this->GBwidth1 = sqrt(GBwidth1)/(this->density_nbPts[0]);
	}
	cout << "GB pos and width : " << GBPos1 << " " << GBwidth1 << endl;

}
void Bicrystal::ComputeExcessVolume(){
	cout << "Compute EV" << endl;
}
