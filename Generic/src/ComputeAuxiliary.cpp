#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include "ComputeAuxiliary.h"
#include "AtomHicConfig.h"
#include "Crystal.h"
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <omp.h>

using namespace std;

// Compute bond orientational parameter, based on the work of Steinhardt, P. J. et al. 1983, modified by Chua et al. 2010 and modified by me for accounting for multisite crystals 
double* ComputeAuxiliary::BondOrientationalParameter(const int& l_sph, double& rc){
	// if neighbours have not been searched perform the research
	if( !_MySystem->getIsNeighbours() ) _MySystem->searchNeighbours(rc);
	cout << "Computing bond orientation parameter.. ";
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	BondOriParam = new double[nbAt];
	unsigned int *Malpha = new unsigned int[nbAt*(nbNMax+1)]; // array containing the index of neighbours of the same species (or same site in case of multisite crystal) with the first line corresponding to the number of neighbours, i.e. Malpha[i*(nbNMax+1)] = nb of neighbour of atom i, Malpha[i*(nbNMax+1)+j+1] = id of the jth neighbour of atom i
	complex<double> *Qalpha = new complex<double>[nbAt*(l_sph*2+1)]; // complex array containing the spherical harmonic for the different modes
	double *Calpha = new double[nbAt]; // normalization factor 
	for(unsigned int i=0;i<nbAt*(l_sph*2+1);i++) Qalpha[i] = (0.,0.); // initialize it to zero
	double zeronum = 1e-8;
	const int bar_length = 30;
	double prog=0;
	double xpos,ypos,zpos,xp,yp,zp, colat, longit;
	unsigned int id;
	// if the crystal has been defined, search if it is a multisite one, wuthout crystal definition it is assumed as monosite
	bool multiSite = false;
	if( _MySystem->getIsCrystalDefined() ){
		for(unsigned int i=0;i<_MySystem->getCrystal()->getNbAtomType();i++){
			if(_MySystem->getCrystal()->getNbAtomSite(i) > 1){

				multiSite = true;
				break;
			}
		}
	}

	// loop on all atoms and neighbours to compute Qalpha and store neighbour of the same species
	// Here is the most time consuming loop of the function, use parallel computation
	unsigned int j_loop;
	int l_loop;
	#pragma omp parallel for private(xpos,ypos,zpos,j_loop,id,xp,yp,zp,colat,longit,l_loop)
	for(unsigned int i=0;i<nbAt;i++){
		xpos = _MySystem->getWrappedPos(i).x;
		ypos = _MySystem->getWrappedPos(i).y;
		zpos = _MySystem->getWrappedPos(i).z;
		Malpha[i*(nbNMax+1)] = 0; 
		Calpha[i] = 0; 
		for(j_loop=0;j_loop<_MySystem->getNeighbours(i*(nbNMax+1));j_loop++){
			id = _MySystem->getNeighbours(i*(nbNMax+1)+j_loop+1);
			// get distance vector
			xp = _MySystem->getWrappedPos(id).x+_MySystem->getCLNeighbours(i*nbNMax*3+j_loop*3)*_MySystem->getH1()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+j_loop*3+1)*_MySystem->getH2()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+j_loop*3+2)*_MySystem->getH3()[0]-xpos;
			yp = _MySystem->getWrappedPos(id).y+_MySystem->getCLNeighbours(i*nbNMax*3+j_loop*3)*_MySystem->getH1()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+j_loop*3+1)*_MySystem->getH2()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+j_loop*3+2)*_MySystem->getH3()[1]-ypos;
			zp = _MySystem->getWrappedPos(id).z+_MySystem->getCLNeighbours(i*nbNMax*3+j_loop*3)*_MySystem->getH1()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+j_loop*3+1)*_MySystem->getH2()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+j_loop*3+2)*_MySystem->getH3()[2]-zpos;
			// compute colatitude and longitudinal angles
			colat = acos(zp/sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)));
			if( xp > 0 ) longit = atan(yp/xp);
	                else if( ( xp < 0 ) && ( yp >= 0 ) ) longit = atan(yp/xp) + M_PI;
	                else if( ( xp < 0 ) and ( yp < 0 ) ) longit = atan(yp/xp) - M_PI;
	                else if( ( fabs(xp) < zeronum ) and ( yp > 0 ) ) longit = M_PI/2.;
	                else if( ( fabs(xp) < zeronum ) and ( yp < 0 ) ) longit = -M_PI/2.;
	                else if( ( fabs(xp) < zeronum ) and ( fabs(yp) < zeronum ) ) longit = 0.;
			// compute spherical harmonics
			for(l_loop=-l_sph;l_loop<l_sph+1;l_loop++)	Qalpha[i*(l_sph*2+1)+l_loop+l_sph] += spherical_harmonics((unsigned int) l_sph, l_loop, colat, longit);
			// Store the neighbour index into Malpha if it is of the same specy
			if( _MySystem->getAtom(i).type == _MySystem->getAtom(id).type ){
				Malpha[i*(nbNMax+1)] += 1;
				Malpha[i*(nbNMax+1)+Malpha[i*(nbNMax+1)]] = id;
			}
		}
		// compute normalization factors
		for(int l=-l_sph;l<l_sph+1;l++)	Calpha[i] += (pow(Qalpha[i*(l_sph*2+1)+l+l_sph].real(), 2.) + pow(Qalpha[i*(l_sph*2+1)+l+l_sph].imag(), 2.));
	}
	// compute the order parameter using the formulation presented in Chua et al. 2010
	unsigned int NId, nbN;
	for(unsigned int i=0;i<nbAt;i++){
		BondOriParam[i] = 0;
		nbN = Malpha[i*(nbNMax+1)];
		for(unsigned int j=0;j<nbN;j++){
			NId = Malpha[i*(nbNMax+1)+j+1];
			for(unsigned int l=0;l<(l_sph*2+1);l++){
				BondOriParam[i] += ((Qalpha[i*(l_sph*2+1)+l].real()*Qalpha[NId*(l_sph*2+1)+l].real())+Qalpha[i*(l_sph*2+1)+l].imag()*Qalpha[NId*(l_sph*2+1)+l].imag())/(pow(Calpha[i],.5)*pow(Calpha[NId],.5)); 
			}
		}
		if( nbN == 0 ) BondOriParam[i] = 0;
		else BondOriParam[i] /= nbN;
	}

	// in the case of mutlisite crystal and/or non-centrosymmetric crystal, each atom belonging to a site will have a given BondOriParam
	// we then search the most represented BondOriParam and renormalize by the closest one 
	// this method implies to consider that most of the atoms in the system are in perfect environment
	if( multiSite ){
		// Search the most represented normalization factors
		vector<vector<double>> NormFactors; // array containing the normalization factors and the number of atom having the normalization factor for a given element, i.e. : NormFactors[i][0] = chemical element (type_uint), NormFactors[i][j*2+1] = jth normalization factor for specy i, NormFactors[i][j*2+2] = number of atom having this normalization factor
		bool ElemStored, NormFacStored;
		double tolSites = 5e-2; // !!! one of the critical values !!!
		NormFactors.push_back(vector<double>());
		NormFactors[0].push_back(_MySystem->getAtom(0).type_uint);
		NormFactors[0].push_back(BondOriParam[0]);
		NormFactors[0].push_back(1);
		for(unsigned int i=1;i<nbAt;i++){
			ElemStored = false;
			for(unsigned int j=0;j<NormFactors.size();j++){
				if( _MySystem->getAtom(i).type_uint == (unsigned int) round(NormFactors[j][0]) ){ // chemical specy already stored, use the same vector column
					NormFacStored = false;
					for(unsigned int k=0;k<(NormFactors[j].size()-1)/2;k++){
						if( fabs(BondOriParam[i]-NormFactors[j][k*2+1]) < tolSites ){ // norm factor already store, increment the counter and average the normalization factor according to the current one (to test if it is better or not)
							NormFactors[j][k*2+1] *= NormFactors[j][k*2+2];
							NormFactors[j][k*2+1] += BondOriParam[i];
							NormFactors[j][k*2+2] += 1; // maybe just keep this line..
							NormFactors[j][k*2+1] /= NormFactors[j][k*2+2];
							NormFacStored = true;
							break;
						}
					}
					if( !NormFacStored ){ // this norm factor has not been found before : create a new line and initialize the counter to 1
						NormFactors[j].push_back(BondOriParam[i]);
						NormFactors[j].push_back(1);
					}
					ElemStored = true;
					break;
				}
			}
			if( !ElemStored ){ // this specy has not been found before : create a new column in NormFactors
				NormFactors.push_back(vector<double>());
				NormFactors[NormFactors.size()-1].push_back(_MySystem->getAtom(i).type_uint);
				NormFactors[NormFactors.size()-1].push_back(BondOriParam[i]);
				NormFactors[NormFactors.size()-1].push_back(1);
			}
		}
		// Initialize the norm fac site array
		vector<vector<double>> FinalNormFactors;
		unsigned int Index;
		for(unsigned int i=0;i<NormFactors.size();i++){
			FinalNormFactors.push_back(vector<double>());
			Index = (unsigned int) round(NormFactors[i][0]);
			if( _MySystem->getCrystal()->getNbAtomSite(Index) > (NormFactors[i].size()-1)/2 ) cout << "Not enough sites have been found, consider decreasing tolSite" << endl;
			for(unsigned int j=0;j<_MySystem->getCrystal()->getNbAtomSite(Index);j++) FinalNormFactors[i].push_back(0.);
		}
		// search the most represented norm factors
		unsigned int kmax;
		for(unsigned int i=0;i<FinalNormFactors.size();i++){
			for(unsigned int j=0;j<FinalNormFactors[i].size();j++){
				kmax = 0;
				for(unsigned int k=1;k<(NormFactors[i].size()-1)/2;k++)	if( NormFactors[i][(k+1)*2] > NormFactors[i][(kmax+1)*2] ) kmax = k;
				FinalNormFactors[i][j] = NormFactors[i][(kmax+1)*2-1];
				NormFactors[i].erase(NormFactors[i].begin()+(kmax+1)*2);
				NormFactors[i].erase(NormFactors[i].begin()+(kmax+1)*2-1);
			}
		}

		// normalize order parameters according to the closest but higher most represented bondoriparam
		double tolRenorm = 1.15; // critical parameter
		double buffer_d, buffer_d_2, buffer_d_3;
		for(unsigned int i=0;i<nbAt;i++){
			for(unsigned int j=0;j<FinalNormFactors.size();j++){
				if( _MySystem->getAtom(i).type_uint == (unsigned int) round(NormFactors[j][0]) ){
					buffer_d = fabs(BondOriParam[i]-FinalNormFactors[j][0]);
					buffer_d_2 = (FinalNormFactors[j][0]*tolRenorm)-BondOriParam[i];
					kmax = 0;
					for(unsigned int k=1;k<FinalNormFactors[j].size();k++){
						buffer_d_3 = (FinalNormFactors[j][k]*tolRenorm)-BondOriParam[i];
						if( ( (buffer_d_3 > 0.) && (fabs(BondOriParam[i]-FinalNormFactors[j][k]) < buffer_d) ) || ( ( buffer_d_2 < 0. ) && ( buffer_d_3 > 0. ) ) ){
								buffer_d = fabs(BondOriParam[i]-FinalNormFactors[j][k]);
								buffer_d_2 = buffer_d_3;
								kmax = k;
						}
					}
					BondOriParam[i] /= FinalNormFactors[j][kmax];
					break;
				}
			}
		}
	} // end if multisite
	for(unsigned int i=0;i<nbAt;i++){
		if( BondOriParam[i] > 1. || BondOriParam[i] < -1. ) BondOriParam[i] = 1.;
		BondOriParam[i] = 1.-fabs(BondOriParam[i]);
	}
	cout << " Done !" << endl;
	// free allocated memory
	delete[] Malpha;
	delete[] Qalpha;
	delete[] Calpha;
	IsBondOriParam = true;
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
	if( IsBondOriParam ) delete[] BondOriParam;
}

