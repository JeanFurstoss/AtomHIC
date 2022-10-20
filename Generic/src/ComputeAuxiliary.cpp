#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include "ComputeAuxiliary.h"
#include "AtomHicConfig.h"
#include <complex>
#include <iostream>

using namespace std;

// Compute bond orientational parameter, based on the work of Steinhardt, P. J. et al. 1983, modified by Chua et al. 2010 and modified by me for accounting for crystal with different sites 
double* ComputeAuxiliary::BondOrientationalParameter(const int& l_sph, double& rc){
	IsBondOriParam = true;
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	BondOriParam = new double[nbAt];
	unsigned int *Malpha = new unsigned int[nbAt*(nbNMax+1)]; // array containing the index of neighbours of the same species with the first line corresponding to the number of neighbours
	complex<double> *Qalpha = new complex<double>[nbAt*l_sph*2]; // complex array containing the spherical harmonic for the different modes
	for(unsigned int i=0;i<nbAt*l_sph*2;i++) Qalpha[i] = (0.,0.); // initialize it to zero
	double tolSites = 5e-0;
	double zeronum = 1e-8;
	const int bar_length = 30;
	double prog=0;
	double xpos,ypos,zpos,xp,yp,zp, colat, longit;
	unsigned int id;

	// if neighbours have not been searched perform the research
	if( !_MySystem->getIsNeighbours() ) _MySystem->searchNeighbours(rc);

	// loop on all atoms and neighbours to compute Qalpha and store neighbour of the same species
	for(unsigned int i=0;i<nbAt;i++){
		xpos = _MySystem->getAtom(i).pos.x;
		ypos = _MySystem->getAtom(i).pos.y;
		zpos = _MySystem->getAtom(i).pos.z;
		Malpha[i*(nbNMax+1)] = 0; 
		for(unsigned int j=0;j<_MySystem->getNeighbours(i*(nbNMax+1));j++){
			id = _MySystem->getNeighbours(i*(nbNMax+1)+j+1);
			// get distance vector
			xp = _MySystem->getAtom(id).pos.x+_MySystem->getCLNeighbours(i*nbNMax*3+j*3)*_MySystem->getH1()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+1)*_MySystem->getH2()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+2)*_MySystem->getH3()[0]-xpos;
			yp = _MySystem->getAtom(id).pos.y+_MySystem->getCLNeighbours(i*nbNMax*3+j*3)*_MySystem->getH1()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+1)*_MySystem->getH2()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+2)*_MySystem->getH3()[1]-ypos;
			zp = _MySystem->getAtom(id).pos.z+_MySystem->getCLNeighbours(i*nbNMax*3+j*3)*_MySystem->getH1()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+1)*_MySystem->getH2()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+j*3+2)*_MySystem->getH3()[2]-zpos;
			// compute colatitude and longitudinal angles
			colat = acos(zp)/sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.));
			if( xp > 0 ) longit = atan(yp/xp);
	                else if( ( xp < 0 ) && ( yp >= 0 ) ) longit = atan(yp/xp) + M_PI;
	                else if( ( xp < 0 ) and ( yp < 0 ) ) longit = atan(yp/xp) - M_PI;
	                else if( ( fabs(xp) < zeronum ) and ( yp > 0 ) ) longit = M_PI/2.;
	                else if( ( fabs(xp) < zeronum ) and ( yp < 0 ) ) longit = -M_PI/2.;
	                else if( ( fabs(xp) < zeronum ) and ( fabs(yp) < zeronum ) ) longit = 0.;
			// compute spherical harmonics
			for(int l=-l_sph;l<l_sph+1;l++)	Qalpha[i*l_sph*2+l+l_sph] += spherical_harmonics((unsigned int) l_sph, l, colat, longit);
			// Store the neighbour index into Malpha if it is of the same specy
			if( _MySystem->getAtom(i).type == _MySystem->getAtom(id).type ){
				Malpha[i*(nbNMax+1)] += 1;
				Malpha[i*(nbNMax+1)+Malpha[i*(nbNMax+1)]] = id;
			}
		}
	}
	return BondOriParam;

}

complex<double> ComputeAuxiliary::spherical_harmonics(const unsigned int& l, int& m, double& theta, double& phi){
	int mabs;
	if( m < 0 ) mabs = -m;
	else mabs = m;
	double leg = std::sph_legendre(l, mabs, theta);
	if( m < 0 ) leg *= pow(-1., .5*(m-mabs));
	std::complex<double> sph_harm(leg*cos(phi*m),leg*sin(phi*m));
	return sph_harm;
}

ComputeAuxiliary::~ComputeAuxiliary(){
	if( IsBondOriParam ){
		delete[] BondOriParam;

	}
}

