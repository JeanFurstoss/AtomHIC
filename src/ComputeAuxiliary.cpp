#define __STDCPP_WANT_MATH_SPEC_FUNCS__ 1
#include "ComputeAuxiliary.h"
#include "AtomHicConfig.h"
#include "Crystal.h"
#include <complex>
#include <iostream>
#include <vector>
#include <iomanip>
#include <omp.h>
#include <cstring>
#include <sstream>
#include <fstream>
#include <dirent.h>

using namespace std;

// Compute Steinhart parameters which correspond to a vector of dimension l_sph*2+1 for each atom representing the complex norm of Qalpha components
double* ComputeAuxiliary::ComputeSteinhardtParameters(const double rc, const int l_sph){
	// if neighbours have not been searched perform the research
	if( !_MySystem->getIsNeighbours() ){
		_MySystem->searchNeighbours(rc);
	}
	cout << "Computing Steinhart parameters.. ";
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	this->BondOriParam = new double[nbAt];
	this->Malpha = new unsigned int[nbAt*(nbNMax+1)]; // array containing the index of neighbours of the same species (or same site in case of multisite crystal) with the first line corresponding to the number of neighbours, i.e. Malpha[i*(nbNMax+1)] = nb of neighbour of atom i, Malpha[i*(nbNMax+1)+j+1] = id of the jth neighbour of atom i
	this->Qalpha = new complex<double>[nbAt*(l_sph*2+1)]; // complex array containing the spherical harmonic for the different modes
	complex<double> *Qlm = new complex<double>[nbAt*(l_sph*2+1)*(l_sph+1)]; // complex array containing the spherical harmonic for the different modes Qlm[i*(l_sph*2+1)*(l_sph+1)+l*(l_sph*2+1)+m] gives the spherical harmonic for atom i and degree l and m
	unsigned int lsph2 = (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1 = l_sph*2+1;
	complex<double> *buffer_complex = new complex<double>[lsph2]; // complex array containing the spherical harmonic for the different modes Qlm[i*(l_sph*2+1)*(l_sph+1)+l*(l_sph*2+1)+m] gives the spherical harmonic for atom i and degree l and m
	this->SteinhardtParams = new double[nbAt*(l_sph+1)];
	this->SteinhardtParams_ave_cutoff = new double[nbAt*(l_sph+1)];
	this->Calpha = new double[nbAt]; // normalization factor 
	for(unsigned int i=0;i<nbAt*(l_sph*2+1);i++){
		Qalpha[i] = (0.,0.); // initialize it to zero
		for(unsigned int j=0;j<l_sph+1;j++) Qlm[i*j+j] = (0.,0.);
	}
	double zeronum = 1e-8;
	const int bar_length = 30;
	double prog=0;
	double xpos,ypos,zpos,xp,yp,zp, colat, longit;
	unsigned int id;
	// loop on all atoms and neighbours to compute Qalpha and store neighbour of the same species
	// Here is the most time consuming loop of the function, use parallel computation
	unsigned int j_loop, l_loop_st, l_loop_st2, neigh, l_neigh;
	int l_loop, m_loop_st, m_loop_st0, m_loop_st1, m_loop_st2, m_loop_st3;
	#pragma omp parallel for private(xpos,ypos,zpos,j_loop,id,xp,yp,zp,colat,longit,l_loop,l_loop_st,m_loop_st,l_loop_st2,m_loop_st2,m_loop_st3,neigh,l_neigh,m_loop_st0,m_loop_st1)
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
			//for(l_loop=-l_sph;l_loop<l_sph+1;l_loop++) Qalpha[i*(l_sph*2+1)+l_loop+l_sph] += spherical_harmonics((unsigned int) l_sph, l_loop, colat, longit); // maybe false here
			for(l_loop_st=0;l_loop_st<l_sph+1;l_loop_st++){
				for(m_loop_st=-l_loop_st;m_loop_st<(int) l_loop_st+1;m_loop_st++){
					Qlm[i*lsph2 + l_loop_st*lsph1 + (unsigned int) (m_loop_st + (int) l_sph)] += spherical_harmonics(l_loop_st, m_loop_st, colat, longit);
				}
			}
			// Store the neighbour index into Malpha if it is of the same specy
			if( _MySystem->getAtom(i).type_uint == _MySystem->getAtom(id).type_uint ){
				Malpha[i*(nbNMax+1)] += 1;
				Malpha[i*(nbNMax+1)+Malpha[i*(nbNMax+1)]] = id;
			}
		}
		for(l_loop_st2=0;l_loop_st2<l_sph+1;l_loop_st2++){
			SteinhardtParams[i*(l_sph+1)+l_loop_st2] = 0.; 
			for(m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
				SteinhardtParams[i*(l_sph+1)+l_loop_st2] += norm(Qlm[i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
			}
			SteinhardtParams[i*(l_sph+1)+l_loop_st2] *= 4.*M_PI/(2.*l_loop_st2+1.); 
			SteinhardtParams[i*(l_sph+1)+l_loop_st2] /= pow(_MySystem->getNeighbours(i*(nbNMax+1)),2.);
			SteinhardtParams[i*(l_sph+1)+l_loop_st2] = sqrt(SteinhardtParams[i*(l_sph+1)+l_loop_st2]);
		}
		for(m_loop_st3=-l_sph;m_loop_st3<(int) l_sph+1;m_loop_st3++) Qalpha[i*(l_sph*2+1) + (unsigned int) (m_loop_st3 + (int) l_sph)] += Qlm[i*lsph2 + l_sph*lsph1+ (unsigned int) (m_loop_st3 + (int) l_sph)]; 
		// compute normalization factors
		for(int l=-l_sph;l<l_sph+1;l++)	Calpha[i] += (pow(Qalpha[i*(l_sph*2+1)+l+l_sph].real(), 2.) + pow(Qalpha[i*(l_sph*2+1)+l+l_sph].imag(), 2.)); // warning : l is not protected during parallel calc
	}
	//#pragma omp parallel for private(l_loop_st2,m_loop_st0,buffer_complex,neigh,m_loop_st1,m_loop_st2)
	for(unsigned int i=0;i<nbAt;i++){
		for(l_loop_st2=0;l_loop_st2<l_sph+1;l_loop_st2++){
			SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2] = 0.; 
			for(m_loop_st0=-l_loop_st2;m_loop_st0<(int) l_loop_st2+1;m_loop_st0++){
				buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] = Qlm[i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] / (double) _MySystem->getNeighbours(i*(nbNMax+1));
			}
			for(neigh=0;neigh<Malpha[i*(nbNMax+1)];neigh++){
			//for(neigh=0;neigh<_MySystem->getNeighbours(i*(nbNMax+1));neigh++){
			       for(m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++){
				       buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[Malpha[i*(nbNMax+1)+neigh+1]*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) _MySystem->getNeighbours(Malpha[i*(nbNMax+1)+neigh+1]*(nbNMax+1)));
				       //buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) _MySystem->getNeighbours(_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)*(nbNMax+1)));
			       }
			}
			for(m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
				SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2] += norm(buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
			}
			SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2] *= 4.*M_PI/(2.*l_loop_st2+1.); 
			SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2] /= pow(Malpha[i*(nbNMax+1)],2.);
			//SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2] /= pow(_MySystem->getNeighbours(i*(nbNMax+1)),2.);
			SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2] = sqrt(SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2]);
		}
	}

	//for(unsigned int i=0;i<nbAt;i++) cout << pow(SteinhardtParams[i*(l_sph+1)+l_sph]*_MySystem->getNeighbours(i*(nbNMax+1)),2.)*((2.*l_sph+1.)/(4.*M_PI)) << " " << Calpha[i] << endl;
	////cout << i << " "<< _MySystem->getNeighbours(i*(nbNMax+1)) << " " << Malpha[i*(nbNMax+1)] << " ";
	//for(unsigned int l=0;l<l_sph+1;l++) cout << SteinhardtParams[i*(l_sph+1)+l] << " ";
	//cout << endl;
	//}
	cout << "Done !" << endl;
	return this->SteinhardtParams;
}

// Same as above but with the average version where the Steinhardt params are averaged using all neighbours (in the case above only the ions of same type inside the cutoff radius are used)
double* ComputeAuxiliary::ComputeSteinhardtParameters_Multi(const double rc, const int l_sph){
	// if neighbours have not been searched perform the research
	if( !_MySystem->getIsNeighbours() ){
		_MySystem->searchNeighbours(rc);
	}
	cout << "Computing Steinhart parameters.. ";
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	this->BondOriParam = new double[nbAt];
	this->Malpha = new unsigned int[nbAt*(nbNMax+1)]; // array containing the index of neighbours of the same species (or same site in case of multisite crystal) with the first line corresponding to the number of neighbours, i.e. Malpha[i*(nbNMax+1)] = nb of neighbour of atom i, Malpha[i*(nbNMax+1)+j+1] = id of the jth neighbour of atom i
	this->Qalpha = new complex<double>[nbAt*(l_sph*2+1)]; // complex array containing the spherical harmonic for the different modes
	complex<double> *Qlm = new complex<double>[nbAt*(l_sph*2+1)*(l_sph+1)]; // complex array containing the spherical harmonic for the different modes Qlm[i*(l_sph*2+1)*(l_sph+1)+l*(l_sph*2+1)+m] gives the spherical harmonic for atom i and degree l and m
	unsigned int lsph2 = (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1 = l_sph*2+1;
	complex<double> *buffer_complex = new complex<double>[lsph2]; // complex array containing the spherical harmonic for the different modes Qlm[i*(l_sph*2+1)*(l_sph+1)+l*(l_sph*2+1)+m] gives the spherical harmonic for atom i and degree l and m
	this->SteinhardtParams = new double[nbAt*(l_sph+1)];
	this->SteinhardtParams_ave_cutoff = new double[nbAt*(l_sph+1)];
	this->Calpha = new double[nbAt]; // normalization factor 
	for(unsigned int i=0;i<nbAt*(l_sph*2+1);i++){
		Qalpha[i] = (0.,0.); // initialize it to zero
		for(unsigned int j=0;j<l_sph+1;j++) Qlm[i*j+j] = (0.,0.);
	}
	double zeronum = 1e-8;
	const int bar_length = 30;
	double prog=0;
	double xpos,ypos,zpos,xp,yp,zp, colat, longit;
	unsigned int id;
	// loop on all atoms and neighbours to compute Qalpha and store neighbour of the same species
	// Here is the most time consuming loop of the function, use parallel computation
	unsigned int j_loop, l_loop_st, l_loop_st2, neigh, l_neigh;
	int l_loop, m_loop_st, m_loop_st0, m_loop_st1, m_loop_st2, m_loop_st3;
	#pragma omp parallel for private(xpos,ypos,zpos,j_loop,id,xp,yp,zp,colat,longit,l_loop,l_loop_st,m_loop_st,l_loop_st2,m_loop_st2,m_loop_st3,neigh,l_neigh,m_loop_st0,m_loop_st1)
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
			//for(l_loop=-l_sph;l_loop<l_sph+1;l_loop++) Qalpha[i*(l_sph*2+1)+l_loop+l_sph] += spherical_harmonics((unsigned int) l_sph, l_loop, colat, longit); // maybe false here
			for(l_loop_st=0;l_loop_st<l_sph+1;l_loop_st++){
				for(m_loop_st=-l_loop_st;m_loop_st<(int) l_loop_st+1;m_loop_st++){
					Qlm[i*lsph2 + l_loop_st*lsph1 + (unsigned int) (m_loop_st + (int) l_sph)] += spherical_harmonics(l_loop_st, m_loop_st, colat, longit);
				}
			}
			// Store the neighbour index into Malpha if it is of the same specy
			if( _MySystem->getAtom(i).type_uint == _MySystem->getAtom(id).type_uint ){
				Malpha[i*(nbNMax+1)] += 1;
				Malpha[i*(nbNMax+1)+Malpha[i*(nbNMax+1)]] = id;
			}
		}
		for(l_loop_st2=0;l_loop_st2<l_sph+1;l_loop_st2++){
			SteinhardtParams[i*(l_sph+1)+l_loop_st2] = 0.; 
			for(m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
				SteinhardtParams[i*(l_sph+1)+l_loop_st2] += norm(Qlm[i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
			}
			SteinhardtParams[i*(l_sph+1)+l_loop_st2] *= 4.*M_PI/(2.*l_loop_st2+1.); 
			SteinhardtParams[i*(l_sph+1)+l_loop_st2] /= pow(_MySystem->getNeighbours(i*(nbNMax+1)),2.);
			SteinhardtParams[i*(l_sph+1)+l_loop_st2] = sqrt(SteinhardtParams[i*(l_sph+1)+l_loop_st2]);
		}
		for(m_loop_st3=-l_sph;m_loop_st3<(int) l_sph+1;m_loop_st3++) Qalpha[i*(l_sph*2+1) + (unsigned int) (m_loop_st3 + (int) l_sph)] += Qlm[i*lsph2 + l_sph*lsph1+ (unsigned int) (m_loop_st3 + (int) l_sph)]; 
		// compute normalization factors
		for(int l=-l_sph;l<l_sph+1;l++)	Calpha[i] += (pow(Qalpha[i*(l_sph*2+1)+l+l_sph].real(), 2.) + pow(Qalpha[i*(l_sph*2+1)+l+l_sph].imag(), 2.)); // warning : l is not protected during parallel calc
	}
	//#pragma omp parallel for private(l_loop_st2,m_loop_st0,buffer_complex,neigh,m_loop_st1,m_loop_st2)
	for(unsigned int i=0;i<nbAt;i++){
		for(l_loop_st2=0;l_loop_st2<l_sph+1;l_loop_st2++){
			SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2] = 0.; 
			for(m_loop_st0=-l_loop_st2;m_loop_st0<(int) l_loop_st2+1;m_loop_st0++){
				buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] = Qlm[i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] / (double) _MySystem->getNeighbours(i*(nbNMax+1));
			}
			//for(neigh=0;neigh<Malpha[i*(nbNMax+1)];neigh++){
			for(neigh=0;neigh<_MySystem->getNeighbours(i*(nbNMax+1));neigh++){
			       for(m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++){
				       //buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[Malpha[i*(nbNMax+1)+neigh+1]*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) _MySystem->getNeighbours(Malpha[i*(nbNMax+1)+neigh+1]*(nbNMax+1)));
				       buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) _MySystem->getNeighbours(_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)*(nbNMax+1)));
			       }
			}
			for(m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
				SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2] += norm(buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
			}
			SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2] *= 4.*M_PI/(2.*l_loop_st2+1.); 
			//SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2] /= pow(Malpha[i*(nbNMax+1)],2.);
			SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2] /= pow(_MySystem->getNeighbours(i*(nbNMax+1)),2.);
			SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2] = sqrt(SteinhardtParams_ave_cutoff[i*(l_sph+1)+l_loop_st2]);
		}
	}

	//for(unsigned int i=0;i<nbAt;i++) cout << pow(SteinhardtParams[i*(l_sph+1)+l_sph]*_MySystem->getNeighbours(i*(nbNMax+1)),2.)*((2.*l_sph+1.)/(4.*M_PI)) << " " << Calpha[i] << endl;
	////cout << i << " "<< _MySystem->getNeighbours(i*(nbNMax+1)) << " " << Malpha[i*(nbNMax+1)] << " ";
	//for(unsigned int l=0;l<l_sph+1;l++) cout << SteinhardtParams[i*(l_sph+1)+l] << " ";
	//cout << endl;
	//}
	cout << "Done !" << endl;
	return this->SteinhardtParams;
}

double* ComputeAuxiliary::ComputeSteinhardtParameters_OneL(const double rc, const int l_sph){ //TODO Rename this function
	// if neighbours have not been searched perform the research
	if( !_MySystem->getIsNeighbours() ){
		_MySystem->searchNeighbours(rc);
	}
	cout << "Computing Steinhart parameters.. ";
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	this->Malpha = new unsigned int[nbAt*(nbNMax+1)]; // array containing the index of neighbours of the same species (or same site in case of multisite crystal) with the first line corresponding to the number of neighbours, i.e. Malpha[i*(nbNMax+1)] = nb of neighbour of atom i, Malpha[i*(nbNMax+1)+j+1] = id of the jth neighbour of atom i
	this->Qalpha = new complex<double>[nbAt*(l_sph*2+1)]; // complex array containing the spherical harmonic for the different modes
	unsigned int lsph2 = (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1 = l_sph*2+1;
	this->Calpha = new double[nbAt]; // normalization factor 
	for(unsigned int i=0;i<nbAt*(l_sph*2+1);i++){
		Qalpha[i] = (0.,0.); // initialize it to zero
	}
	double zeronum = 1e-8;
	const int bar_length = 30;
	double prog=0;
	double xpos,ypos,zpos,xp,yp,zp, colat, longit;
	unsigned int id;
	// loop on all atoms and neighbours to compute Qalpha and store neighbour of the same species
	// Here is the most time consuming loop of the function, use parallel computation
	unsigned int j_loop, l_loop_st, l_loop_st2;
	int l_loop, m_loop_st, m_loop_st2, m_loop_st3;
	#pragma omp parallel for private(xpos,ypos,zpos,j_loop,id,xp,yp,zp,colat,longit,l_loop,l_loop_st,m_loop_st,l_loop_st2,m_loop_st2,m_loop_st3)
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
			for(l_loop=-l_sph;l_loop<l_sph+1;l_loop++) Qalpha[i*(l_sph*2+1)+l_loop+l_sph] += spherical_harmonics((unsigned int) l_sph, l_loop, colat, longit); // maybe false here
			// Store the neighbour index into Malpha if it is of the same specy
			if( _MySystem->getAtom(i).type_uint == _MySystem->getAtom(id).type_uint ){
				Malpha[i*(nbNMax+1)] += 1;
				Malpha[i*(nbNMax+1)+Malpha[i*(nbNMax+1)]] = id;
			}
		}
		// compute normalization factors
		for(int l=-l_sph;l<l_sph+1;l++)	Calpha[i] += (pow(Qalpha[i*(l_sph*2+1)+l+l_sph].real(), 2.) + pow(Qalpha[i*(l_sph*2+1)+l+l_sph].imag(), 2.)); // warning : l is not protected during parallel calc
	}
	// Test averaged within cutoff radius (not concluant I think because contrarily to the true averaged Steinhardt params where the complex are normed, the spherical harmonic are not rotationally invariant)
	//complex<double> *Qalpha_ave_cutoff = new complex<double>[nbAt*(l_sph*2+1)]; // complex array containing the spherical harmonic for the different modes
	//for(unsigned int i=0;i<nbAt;i++){
	//	for(int m_loop_st0=-l_loop_st2;m_loop_st0<(int) l_loop_st2+1;m_loop_st0++){
	//		Qalpha_ave_cutoff[i*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] = Qalpha[i*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] / (double) _MySystem->getNeighbours(i*(nbNMax+1));
	//	}
	//	for(unsigned int neigh=0;neigh<Malpha[i*(nbNMax+1)];neigh++){
	//	       for(int m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++){
	//		       Qalpha_ave_cutoff[i*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qalpha[Malpha[i*(nbNMax+1)+neigh+1]*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) _MySystem->getNeighbours(Malpha[i*(nbNMax+1)+neigh+1]*(nbNMax+1)));
	//	       }
	//	}
	//}

	//for(unsigned int i=0;i<nbAt;i++){
	//	for(int m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++) Qalpha[i*lsph1 +  (unsigned int) (m_loop_st1 + (int) l_sph)] = Qalpha_ave_cutoff[i*lsph1+ (unsigned int) (m_loop_st1 + (int) l_sph)];
	//	for(int l=-l_sph;l<l_sph+1;l++)	Calpha[i] += (pow(Qalpha[i*(l_sph*2+1)+l+l_sph].real(), 2.) + pow(Qalpha[i*(l_sph*2+1)+l+l_sph].imag(), 2.)); // warning : l is not protected during parallel calc
	//}
	//cout << "Done !" << endl;
}
double* ComputeAuxiliary::BondOriParam_SteinhardtBased(){
	// read database
	SteinhardtDatabase_read(_MySystem->getCrystal()->getName());
	// compute Steinhardt params
	ComputeSteinhardtParameters(this->rcut_ref, this->l_sph_ref);
	// compute Calpha and normalization factor for references
	double *NormSteinhardt_ref = new double[AtomTypeUINTRefPC.size()];
	for(unsigned int i=0;i<AtomTypeUINTRefPC.size();i++){
		NormSteinhardt_ref[i] = 0.;
		for(unsigned int l=0;l<this->l_sph_ref+1;l++){
			NormSteinhardt_ref[i] += pow(this->SteinhardtParams_REF_PC[i][l], 2.);
		}
		NormSteinhardt_ref[i] = sqrt(NormSteinhardt_ref[i]);
	}
	const unsigned int nbAt = _MySystem->getNbAtom();
	this->BondOriParam_Steinhardt = new double[nbAt];
	//test
	vector<double> *BondOri = new vector<double>[AtomTypeUINTRefPC.size()];
	double *NormFac = new double[AtomTypeUINTRefPC.size()];
	vector<double> BondOriFromRef;
	if( _MySystem->getCrystal()->getIsMultisite() ){ // search from the ref which is the closest to each atom
		unsigned int *closest_ref = new unsigned int[nbAt];
		vector<double> dist;
		vector<unsigned int> ind;
		double norm_dist;
		double cur_norm;
		for(unsigned int i=0;i<nbAt;i++){
			dist.clear();
			ind.clear();
			cur_norm = 0.;
			for(unsigned int l=0;l<this->l_sph_ref+1;l++)	cur_norm += pow(this->SteinhardtParams[i*(this->l_sph_ref+1)+l],2.);
			for(unsigned int j=0;j<AtomTypeUINTRefPC.size();j++){
				if( _MySystem->getAtom(i).type_uint == AtomTypeUINTRefPC[j] ){
					dist.push_back(0.);
					ind.push_back(j);
					for(unsigned int l=0;l<this->l_sph_ref+1;l++)	dist[dist.size()-1] += pow(this->SteinhardtParams[i*(this->l_sph_ref+1)+l]-this->SteinhardtParams_REF_PC[j][l],2.);
					dist[dist.size()-1] = sqrt(dist[dist.size()-1]);
					dist[dist.size()-1] /= (sqrt(cur_norm)+NormSteinhardt_ref[j])/2.;
				}
			}
			closest_ref[i] = ind[MT->min(dist)];
			this->BondOriParam_Steinhardt[i] = MT->min_vec(dist);
		}
	}//end multisite
	//ofstream writefile("closest_ref.dat");
	//for(unsigned int i=0;i<nbAt;i++) writefile << closest_ref[i] << endl;
	//writefile.close();
	
	return this->BondOriParam_Steinhardt;
}

// Compute bond orientational parameter, based on the work of Steinhardt, P. J. et al. 1983, modified by Chua et al. 2010 and modified by me for accounting for multisite crystals 
double* ComputeAuxiliary::BondOrientationalParameter(){
	if( _MySystem->getCrystal()->getIsMultisite() )	BondOriParam_MultisiteNewVersion();
	else BondOriParam_NoMultisite();
	return this->BondOriParam;
}

void ComputeAuxiliary::BondOriParam_NoMultisite(){
	double rc = _MySystem->get_rcut();
	int l_sph = _MySystem->get_lsph();
	ComputeSteinhardtParameters_OneL(rc,l_sph);
	cout << "Computing bond orientationnal parameter" << endl;

	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	if( !this->IsBondOriParam ){
		this->BondOriParam = new double[nbAt];
		this->IsBondOriParam = true;
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
}

void ComputeAuxiliary::BondOriParam_MultisiteNewVersion(){
	double rc = _MySystem->get_rcut();
	int l_sph = _MySystem->get_lsph();
	ComputeSteinhardtParameters_OneL(rc,l_sph);
	if( _MySystem->getCrystal()->getName() != "" ) SteinhardtDatabase_read(_MySystem->getCrystal()->getName());
	// search site index based on the value of Steinhardt parameter (which is sqrt(pow(Calpha*N,2.)*4pi/(2l+1)))
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	this->Atom_SiteIndex = new unsigned int[nbAt];
	vector<double> diff_St;
	double buffer_d, buffer_d1;
	cout << "steinhardt" << endl;
	for(unsigned int i=0;i<nbAt;i++){
		if( this->AtomSiteRefPC[_MySystem->getAtom(i).type_uint-1].size() == 1 ) this->Atom_SiteIndex[i] = 0;
		else{
			diff_St.clear();
			buffer_d = sqrt(Calpha[i]*4.*M_PI/((2*l_sph+1)*pow(_MySystem->getNeighbours(i*(nbNMax+1)),2.)));
			for(unsigned int j=0;j<this->AtomSiteRefPC[_MySystem->getAtom(i).type_uint-1].size();j++){
				buffer_d1 = SteinhardtParams_REF_PC[_MySystem->getAtom(i).type_uint-1][j*(l_sph+1)+l_sph]-buffer_d;
				diff_St.push_back(fabs(buffer_d1));
				//if( buffer_d1 < 0 ) diff_St.push_back(1e8);
				//else diff_St.push_back(buffer_d1);
			}
			Atom_SiteIndex[i] = MT->min(diff_St);
		}
	}
	// Compute the unnormalized order parameter	
	if( !this->IsBondOriParam ){
		this->BondOriParam = new double[nbAt];
		this->IsBondOriParam = true;
	}
	unsigned int NId, nbN, trueN;
	for(unsigned int i=0;i<nbAt;i++){
		BondOriParam[i] = 0;
		nbN = Malpha[i*(nbNMax+1)];
		trueN = 0;
		for(unsigned int j=0;j<nbN;j++){
			NId = Malpha[i*(nbNMax+1)+j+1];
			if( Atom_SiteIndex[NId] == Atom_SiteIndex[i] ){
				for(unsigned int l=0;l<(l_sph*2+1);l++){
					BondOriParam[i] += ((Qalpha[i*(l_sph*2+1)+l].real()*Qalpha[NId*(l_sph*2+1)+l].real())+Qalpha[i*(l_sph*2+1)+l].imag()*Qalpha[NId*(l_sph*2+1)+l].imag())/(pow(Calpha[i],.5)*pow(Calpha[NId],.5)); 
				}
				trueN++;
			}
		}
		if( trueN == 0 ) BondOriParam[i] = 0;
		else BondOriParam[i] /= trueN;
		//cout << i+1 << " " << trueN << " " << BondOriParam[i] << endl;
	}
	// Normalized the order parameter
	for(unsigned int i=0;i<nbAt;i++){
		BondOriParam[i] /= this->BondOriParam_REF_PC[_MySystem->getAtom(i).type_uint-1][Atom_SiteIndex[i]];
		//if( BondOriParam[i] > 1. ) BondOriParam[i] = 2.-BondOriParam[i];
		if( BondOriParam[i] > 1. || BondOriParam[i] < -1. ) BondOriParam[i] = 1.;
		BondOriParam[i] = 1.-fabs(BondOriParam[i]);
	}
}
	
void ComputeAuxiliary::BondOriParam_Multisite(){ // this version does work (it worked before => have a look if the new version is not satisfactory)
	BondOriParam_NoMultisite();
	if( _MySystem->getCrystal()->getName() != "" ) SteinhardtDatabase_read(_MySystem->getCrystal()->getName());
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	// in the case of mutlisite crystal and/or non-centrosymmetric crystal, each atom belonging to a site will have a given BondOriParam
	// we then search the most represented BondOriParam and renormalize by the closest one 
	// this method implies to consider that most of the atoms in the system are in perfect environment
	// Search the most represented normalization factors
	vector<vector<double>> NormFactors; // array containing the normalization factors and the number of atom having the normalization factor for a given element, i.e. : NormFactors[i][0] = chemical element (type_uint), NormFactors[i][j*2+1] = jth normalization factor for specy i, NormFactors[i][j*2+2] = number of atom having this normalization factor
	vector<vector<unsigned int>> SiteIndex; // array containing the site index 
	this->Atom_SiteIndex = new unsigned int[nbAt];
	// compute the tolerance for the different atom type
	vector<double> tolSites_v(_MySystem->getCrystal()->getNbAtomType());
	vector<vector<double>> BondOri_type(_MySystem->getCrystal()->getNbAtomType());
	for(unsigned int i=0;i<nbAt;i++) BondOri_type[_MySystem->getAtom(i).type_uint-1].push_back(BondOriParam[i]);
	for(unsigned int t=0;t<_MySystem->getCrystal()->getNbAtomType();t++){
		if( _MySystem->getCrystal()->getNbAtomSite(t+1) == 1 ) tolSites_v[t] = MT->max_vec(BondOri_type[t])*tolSites;
		else tolSites_v[t] = (MT->max_vec(BondOri_type[t]) - MT->min_vec(BondOri_type[t]))*tolSites;
	}
	bool ElemStored, NormFacStored;
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
					//if( fabs(BondOriParam[i]-NormFactors[j][k*2+1]) < tolSites_v[_MySystem->getAtom(i).type_uint-1] ){ // norm factor already store, increment the counter and average the normalization factor according to the current one (to test if it is better or not)
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
	unsigned int count_site = 0;
	for(unsigned int i=0;i<NormFactors.size();i++){
		FinalNormFactors.push_back(vector<double>());
		SiteIndex.push_back(vector<unsigned int>());
		Index = (unsigned int) round(NormFactors[i][0]);
		if( _MySystem->getCrystal()->getNbAtomSite(Index) > (NormFactors[i].size()-1)/2 ){
			cerr << "Not enough sites have been found, consider decreasing tolSite" << endl;
			exit(EXIT_FAILURE);
		}
		for(unsigned int j=0;j<_MySystem->getCrystal()->getNbAtomSite(Index);j++){
			FinalNormFactors[i].push_back(0.);
			SiteIndex[i].push_back(count_site);
			count_site++;
		}
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
	double tolRenorm = 1.05; // critical parameter
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
				Atom_SiteIndex[i] = SiteIndex[j][kmax];
				break;
			}
		}
	}
	for(unsigned int i=0;i<nbAt;i++){
		if( BondOriParam[i] > 1. || BondOriParam[i] < -1. ) BondOriParam[i] = 1.;
		BondOriParam[i] = 1.-fabs(BondOriParam[i]);
	}
	cout << " Done !" << endl;
	// free allocated memory
	delete[] Malpha;
	delete[] Qalpha;
	delete[] Calpha;
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

double* ComputeAuxiliary::Compute_StrainTensor(){
	if( !_MySystem->getIsCrystalDefined() ){
		cerr << "The crystal need to be defined for the strain tensor calculation" << endl;
		exit(EXIT_FAILURE);
	}
	if( !IsBondOriParam ){
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
	if( FromNum == 0 ) cout << "Computing strain tensor based on atom type" << endl;
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
	return this->StrainTensor;
}	
	
double* ComputeAuxiliary::Compute_StrainTensor_invII(){
	if( !IsStrainTensor ) Compute_StrainTensor(0);
	const unsigned int nbAt = _MySystem->getNbAtom();
	this->Strain_invII = new double[nbAt];
	this->IsStrainInvII = true;
	for(unsigned int i=0;i<nbAt;i++) this->Strain_invII[i] = this->StrainTensor[i*6+0]*this->StrainTensor[i*6+1] + this->StrainTensor[i*6+1]*this->StrainTensor[i*6+2] + this->StrainTensor[i*6+0]*this->StrainTensor[i*6+2] + this->StrainTensor[i*6+3]*this->StrainTensor[i*6+3] + this->StrainTensor[i*6+4]*this->StrainTensor[i*6+4] + this->StrainTensor[i*6+5]*this->StrainTensor[i*6+5]; 
}

// Save averaged (for each atom type) Steinhardt params
void ComputeAuxiliary::SaveAveSteinhardtParamToDatabase_PerfectCrystal(string CrystalName){
	int l_sph = _MySystem->get_lsph();
	double rc = _MySystem->get_rcut();
	ComputeSteinhardtParameters(rc, l_sph);
	string filename;
	// read the database environment
	if( _MySystem->getCrystal()->getIsMultisite() ){
		filename = "PerfectCrystal_ave.dat";
		cout << "Saving averaged Steinhardt parameters for perfect crystal of " << CrystalName << endl;
	}else{
		filename = "PerfectCrystal.dat";
		cout << "Saving Steinhardt parameters for perfect crystal of " << CrystalName << endl;
	}

	// Save the average parameters over sites
	string fullpathname = getSteinhardtDatabase(CrystalName)+filename;
	const unsigned int nbAt = _MySystem->getNbAtom();

	ofstream writefile(fullpathname);
	if( _MySystem->getCrystal()->getIsMultisite() )	writefile << "Averaged Steinhardt parameters for " << CrystalName << " perfect crystal computed using" << endl;
	else writefile << "Steinhardt parameters for " << CrystalName << " perfect crystal computed using" << endl;
	writefile << "L_SPH " << l_sph << endl;
	writefile << "RCUT " << rc << endl;

	bool Already = false;
	vector<vector<double>> St2Print;
	vector<unsigned int> type_printed;
	vector<unsigned int> count_ave;
	// average the steinhardt parameters
	for(unsigned int i=0;i<nbAt;i++){
		Already = false;
		for(unsigned int t=0;t<type_printed.size();t++){
			if( _MySystem->getAtom(i).type_uint == type_printed[t] ){
				for(unsigned int l=0;l<l_sph+1;l++) St2Print[t][l] += SteinhardtParams[i*(l_sph+1)+l];
				count_ave[t] += 1;
				Already = true;
				break;
			}
		}
		if( !Already ){
			St2Print.push_back(vector<double>());
			count_ave.push_back(1);
			type_printed.push_back(_MySystem->getAtom(i).type_uint);
			for(unsigned int l=0;l<l_sph+1;l++) St2Print[St2Print.size()-1].push_back(SteinhardtParams[i*(l_sph+1)+l]);
		}
	}
	for(unsigned int t=0;t<type_printed.size();t++){
		for(unsigned int l=0;l<l_sph+1;l++) St2Print[t][l] /= count_ave[t];
	}
	writefile << "NB_REF " << type_printed.size() << endl;
	for(unsigned int t=0;t<type_printed.size();t++){
		writefile << _MySystem->getCrystal()->getAtomType(type_printed[t]) << " " << type_printed[t] << " ";
		for(unsigned int l=0;l<l_sph+1;l++) writefile << St2Print[t][l] << " ";
		writefile << endl;
	}
	writefile.close();

	// Save average parameters over atom of same type inside the cutoff radius
	filename = "PerfectCrystal_ave_cutoff.dat";
	fullpathname = getSteinhardtDatabase(CrystalName)+filename;

	ofstream writefile2(fullpathname);
	if( _MySystem->getCrystal()->getIsMultisite() )	writefile2 << "Averaged Steinhardt parameters for atom of same type inside the cutoff radius for " << CrystalName << " perfect crystal computed using" << endl;
	else writefile2 << "Steinhardt parameters averaged over atom of same type inside the cutoff radius for " << CrystalName << " perfect crystal computed using" << endl;
	writefile2 << "L_SPH " << l_sph << endl;
	writefile2 << "RCUT " << rc << endl;

	Already = false;
	for(unsigned int i=0;i<St2Print.size();i++) St2Print.clear();
	type_printed.clear();
	count_ave.clear();
	// average the steinhardt parameters
	for(unsigned int i=0;i<nbAt;i++){
		Already = false;
		for(unsigned int t=0;t<type_printed.size();t++){
			if( _MySystem->getAtom(i).type_uint == type_printed[t] ){
				for(unsigned int l=0;l<l_sph+1;l++) St2Print[t][l] += SteinhardtParams_ave_cutoff[i*(l_sph+1)+l];
				count_ave[t] += 1;
				Already = true;
				break;
			}
		}
		if( !Already ){
			St2Print.push_back(vector<double>());
			count_ave.push_back(1);
			type_printed.push_back(_MySystem->getAtom(i).type_uint);
			for(unsigned int l=0;l<l_sph+1;l++) St2Print[St2Print.size()-1].push_back(SteinhardtParams_ave_cutoff[i*(l_sph+1)+l]);
		}
	}
	for(unsigned int t=0;t<type_printed.size();t++){
		for(unsigned int l=0;l<l_sph+1;l++) St2Print[t][l] /= count_ave[t];
	}
	writefile2 << "NB_REF " << type_printed.size() << endl;
	for(unsigned int t=0;t<type_printed.size();t++){
		writefile2 << _MySystem->getCrystal()->getAtomType(type_printed[t]) << " " << type_printed[t] << " ";
		for(unsigned int l=0;l<l_sph+1;l++) writefile2 << St2Print[t][l] << " ";
		writefile2 << endl;
	}
	writefile2.close();

}

void ComputeAuxiliary::BondOriParam_MultisiteCrystal(){
	double rc = _MySystem->get_rcut();
	int l_sph = _MySystem->get_lsph();
	ComputeSteinhardtParameters_OneL(rc,l_sph);
	cout << "Computing bond orientationnal parameter" << endl;

	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	if( !this->IsBondOriParam ){
		this->BondOriParam = new double[nbAt];
		this->Atom_SiteIndex = new unsigned int[nbAt];
		this->IsBondOriParam = true;
	}

// compute the order parameter using the formulation presented in Chua et al. 2010
	unsigned int NId, nbN, trueNeigh;
	for(unsigned int i=0;i<nbAt;i++){
		BondOriParam[i] = 0;
		nbN = Malpha[i*(nbNMax+1)];
		trueNeigh = 0;
		this->Atom_SiteIndex[i] = _MySystem->getCrystal()->getAtomSite(i);
		for(unsigned int j=0;j<nbN;j++){
			NId = Malpha[i*(nbNMax+1)+j+1];
			if( this->Atom_SiteIndex[i] == _MySystem->getCrystal()->getAtomSite(NId) ){ // restrict the sum to ions of the same site
				for(unsigned int l=0;l<(l_sph*2+1);l++){
					BondOriParam[i] += ((Qalpha[i*(l_sph*2+1)+l].real()*Qalpha[NId*(l_sph*2+1)+l].real())+Qalpha[i*(l_sph*2+1)+l].imag()*Qalpha[NId*(l_sph*2+1)+l].imag())/(pow(Calpha[i],.5)*pow(Calpha[NId],.5)); 
				}
				trueNeigh++;
			}
		}
		if( trueNeigh == 0 ) BondOriParam[i] = 0;
		else BondOriParam[i] /= trueNeigh;
	}
}

// save Steinhardt params for the different sites (defined if multisite crystal)
void ComputeAuxiliary::SaveSteinhardtParamToDatabase_PerfectCrystal(string CrystalName){
	if( _MySystem->getCrystal()->getIsMultisite() ){
		int l_sph = _MySystem->get_lsph();
		double rc = _MySystem->get_rcut();
		ComputeSteinhardtParameters(rc, l_sph);
		// read the database environment
		string filename = "PerfectCrystal.dat";
		string fullpathname = getSteinhardtDatabase(CrystalName)+filename;
		cout << "Saving Steinhardt parameters for perfect crystal of " << CrystalName << endl;
		const unsigned int nbAt = _MySystem->getNbAtom();
		ofstream writefile(fullpathname);
		writefile << "Steinhardt parameters for " << CrystalName << " perfect crystal computed using" << endl;
		writefile << "L_SPH " << l_sph << endl;
		writefile << "RCUT " << rc << endl;
		// compute the bond orientational param and store them for the different sites
		BondOriParam_MultisiteCrystal(); // warning multisite does not work every time
	//for(unsigned int i=0;i<nbAt;i++){
	//cout << i << " ";
	//	cout << BondOriParam[i] << " " << endl;
	////for(unsigned int l=0;l<l_sph+1;l++) cout << this->SteinhardtParams[i*(l_sph+1)+l] << " ";
	////cout << endl;
	//}
		bool Already = false;
		bool AlreadyType = false;
		unsigned int typeok;
		vector<vector<double>> St2Print;
		vector<vector<double>> BO2Print;
		vector<unsigned int> type_printed;
		vector<vector<unsigned int>> site_printed;
		vector<vector<unsigned int>> count_ave;
		// average the steinhardt parameters
		for(unsigned int i=0;i<nbAt;i++){
			Already = false;
			AlreadyType = false;
			for(unsigned int t=0;t<type_printed.size();t++){
				if( _MySystem->getAtom(i).type_uint == type_printed[t] ){
					AlreadyType = true;
					for(unsigned int s=0;s<site_printed[t].size();s++){
						if( _MySystem->getCrystal()->getAtomSite(i) == site_printed[t][s] ){
							BO2Print[t][s] += BondOriParam[i];
							for(unsigned int l=0;l<l_sph+1;l++){
								St2Print[t][s*(l_sph+1)+l] += SteinhardtParams[i*(l_sph+1)+l];
							}
							count_ave[t][s] += 1;
							Already = true;
							break;
						}
					}
					typeok = t;
					break;
				}
			}
			if( !Already ){
				if( !AlreadyType ){ // new type and new site
					St2Print.push_back(vector<double>());
					BO2Print.push_back(vector<double>());
					count_ave.push_back(vector<unsigned int>());
					site_printed.push_back(vector<unsigned int>());
					count_ave[count_ave.size()-1].push_back(1);
					site_printed[site_printed.size()-1].push_back(_MySystem->getCrystal()->getAtomSite(i));
					type_printed.push_back(_MySystem->getAtom(i).type_uint);
					BO2Print[BO2Print.size()-1].push_back(BondOriParam[i]);
					for(unsigned int l=0;l<l_sph+1;l++){
						St2Print[(St2Print.size()-1)].push_back(SteinhardtParams[i*(l_sph+1)+l]);
					}
				}else{ // not a new type, just a new site
					count_ave[typeok].push_back(1);
					site_printed[typeok].push_back(_MySystem->getCrystal()->getAtomSite(i));
					BO2Print[typeok].push_back(BondOriParam[i]);
					for(unsigned int l=0;l<l_sph+1;l++){
						St2Print[typeok].push_back(SteinhardtParams[i*(l_sph+1)+l]);
					}
				}
			}
		}
		unsigned int nbref = 0;
		for(unsigned int t=0;t<type_printed.size();t++){
			for(unsigned int s=0;s<site_printed[t].size();s++){
				nbref += 1;
				BO2Print[t][s] /= count_ave[t][s];
				for(unsigned int l=0;l<l_sph+1;l++){
					St2Print[t][s*(l_sph+1)+l] /= count_ave[t][s];
				}
			}
		}
		writefile << "NB_REF " << nbref << endl;
		for(unsigned int t=0;t<type_printed.size();t++){
			for(unsigned int s=0;s<site_printed[t].size();s++){
				writefile << _MySystem->getCrystal()->getAtomType(type_printed[t]) << " " << type_printed[t] << " " << site_printed[t][s]+1 << " ";
				for(unsigned int l=0;l<l_sph+1;l++) writefile << St2Print[t][s*(l_sph+1)+l] << " ";
				writefile << BO2Print[t][s] << endl;
			}
		}
		writefile.close();
	}else{ // no multisite case
		SaveAveSteinhardtParamToDatabase_PerfectCrystal(CrystalName);
	}

}

void ComputeAuxiliary::SaveSteinhardtParamToDatabase_Defect(string CrystalName, string defect_name, vector<unsigned int> At_index){
	int l_sph = _MySystem->get_lsph();
	double rc = _MySystem->get_rcut();
	ComputeSteinhardtParameters(rc, l_sph);
	string database_ext=".dat";

	// normal Steinhardt params
	string fullpathname = getSteinhardtDatabase(CrystalName)+defect_name+database_ext;
	const unsigned int nbAt = _MySystem->getNbAtom();

	ofstream writefile(fullpathname);
	writefile << "Steinhardt parameters for " << defect_name << " defect in " << CrystalName << " crystal computed using" << endl;
	writefile << "L_SPH " << l_sph << endl;
	writefile << "RCUT " << rc << endl;

	bool Already = false;
	vector<vector<double>> St2Print;
	vector<unsigned int> type_printed;
	vector<unsigned int> count_ave;
	// average the steinhardt parameters
	for(unsigned int i=0;i<At_index.size();i++){
		Already = false;
		for(unsigned int t=0;t<type_printed.size();t++){
			if( _MySystem->getAtom(At_index[i]).type_uint == type_printed[t] ){
				for(unsigned int l=0;l<l_sph+1;l++) St2Print[t][l] += SteinhardtParams[At_index[i]*(l_sph+1)+l];
				count_ave[t] += 1;
				Already = true;
				break;
			}
		}
		if( !Already ){
			St2Print.push_back(vector<double>());
			count_ave.push_back(1);
			type_printed.push_back(_MySystem->getAtom(At_index[i]).type_uint);
			for(unsigned int l=0;l<l_sph+1;l++) St2Print[St2Print.size()-1].push_back(SteinhardtParams[At_index[i]*(l_sph+1)+l]);
		}
	}
	for(unsigned int t=0;t<type_printed.size();t++){
		for(unsigned int l=0;l<l_sph+1;l++) St2Print[t][l] /= count_ave[t];
	}
	writefile << "NB_REF " << type_printed.size() << endl;
	for(unsigned int t=0;t<type_printed.size();t++){
		writefile << _MySystem->getCrystal()->getAtomType(type_printed[t]) << " " << type_printed[t] << " ";
		for(unsigned int l=0;l<l_sph+1;l++) writefile << St2Print[t][l] << " ";
		writefile << endl;
	}
	writefile.close();

	// averaged Steinhardt params over atom of same type inside cutoff radius
	string qual = "_ave_cutoff";
	fullpathname = getSteinhardtDatabase(CrystalName)+defect_name+qual+database_ext;

	ofstream writefile2(fullpathname);
	writefile2 << "Steinhardt parameters averaged over atoms of same specy inside the cutoff radius for " << defect_name << " defect in " << CrystalName << " crystal computed using" << endl;
	writefile2 << "L_SPH " << l_sph << endl;
	writefile2 << "RCUT " << rc << endl;

	Already = false;
	for(unsigned int i=0;i<St2Print.size();i++) St2Print.clear();
	type_printed.clear();
	count_ave.clear();
	// average the steinhardt parameters
	for(unsigned int i=0;i<At_index.size();i++){
		Already = false;
		for(unsigned int t=0;t<type_printed.size();t++){
			if( _MySystem->getAtom(At_index[i]).type_uint == type_printed[t] ){
				for(unsigned int l=0;l<l_sph+1;l++) St2Print[t][l] += SteinhardtParams_ave_cutoff[At_index[i]*(l_sph+1)+l];
				count_ave[t] += 1;
				Already = true;
				break;
			}
		}
		if( !Already ){
			St2Print.push_back(vector<double>());
			count_ave.push_back(1);
			type_printed.push_back(_MySystem->getAtom(At_index[i]).type_uint);
			for(unsigned int l=0;l<l_sph+1;l++) St2Print[St2Print.size()-1].push_back(SteinhardtParams_ave_cutoff[At_index[i]*(l_sph+1)+l]);
		}
	}
	for(unsigned int t=0;t<type_printed.size();t++){
		for(unsigned int l=0;l<l_sph+1;l++) St2Print[t][l] /= count_ave[t];
	}
	writefile2 << "NB_REF " << type_printed.size() << endl;
	for(unsigned int t=0;t<type_printed.size();t++){
		writefile2 << _MySystem->getCrystal()->getAtomType(type_printed[t]) << " " << type_printed[t] << " ";
		for(unsigned int l=0;l<l_sph+1;l++) writefile2 << St2Print[t][l] << " ";
		writefile2 << endl;
	}
	writefile2.close();
}

void ComputeAuxiliary::PrintSteinhardtParam(vector<unsigned int> At_index){
	int l_sph = _MySystem->get_lsph();
	double rc = _MySystem->get_rcut();
	//ComputeSteinhardtParameters_Multi(rc, l_sph);
	ComputeSteinhardtParameters_Multi(rc, l_sph);
	string at_type, filename;
	for(unsigned int t=0;t<_MySystem->getCrystal()->getNbAtomType();t++){
		at_type = _MySystem->getCrystal()->getAtomType()[t];
		filename = at_type+"_Steinhardt.dat";
		ofstream writefile(filename);
		for(unsigned int i=0;i<At_index.size();i++){
			if(_MySystem->getAtom(At_index[i]).type_uint != t+1) continue;
			else{
				writefile << 0 << " " << SteinhardtParams_ave_cutoff[At_index[i]*(l_sph+1)]-1 << " " << At_index[i] << endl;
				for(unsigned int l=1;l<l_sph+1;l++) writefile << l << " " << SteinhardtParams_ave_cutoff[At_index[i]*(l_sph+1)+l] << " " << At_index[i] << endl;
			}
		}
		writefile.close();
	}
}

void ComputeAuxiliary::SteinhardtDatabase_read(string CrystalName){
	if( !this->IsSteinhardtDatabaseRead ){
		string database = getSteinhardtDatabase(CrystalName);	
		string database_extension=".dat";	
		DIR *dir;
		struct dirent *diread;
		const char *env = database.c_str();
		string buffer_s, beg, end_qual;
		string ave_cutoff = "_ave_cutoff";
		size_t pos;
		if( (dir = opendir(env) ) != nullptr ){
			while( (diread = readdir(dir)) != nullptr ){
				buffer_s = diread->d_name;
				pos = buffer_s.find(database_extension);
				if(pos!=string::npos){
					buffer_s.erase(buffer_s.size()-database_extension.size());
					beg = buffer_s.substr(0,1); 
					if( buffer_s.size() < ave_cutoff.size() ) end_qual = "";
					else end_qual = buffer_s.substr(buffer_s.size()-ave_cutoff.size(),buffer_s.size());
					if( buffer_s != "PerfectCrystal" && buffer_s != "PerfectCrystal_ave" && buffer_s != "PerfectCrystal_ave_cutoff" && beg != "." && end_qual != ave_cutoff ) this->Ref_Def_Names.push_back(buffer_s);
				}
			}
			closedir(dir);
		}else{
			perror("opendir");
		}
		// Read perfect crystal Steinhardt params
		string PC_str="PerfectCrystal";
		string fullpath2data = database.c_str()+PC_str+database_extension;
		ifstream file(fullpath2data, ios::in);
		size_t pos_nbref, pos_lsph, pos_rc;
		unsigned int line_rc(1000), count, buffer_ui, buffer_ui2;
		double buffer_d1, buffer_d2;
		for(unsigned int i=0;i<_MySystem->getCrystal()->getNbAtomType();i++){
			this->AtomSiteRefPC.push_back(vector<unsigned int>());
			this->SteinhardtParams_REF_PC.push_back(vector<double>());
			if( _MySystem->getCrystal()->getIsMultisite() ){ // TODO change here and add IsCentrosymmetric to crystal database because we need the reference bond ori param only when it is not centrosymmetric and not when it is multisite
				this->BondOriParam_REF_PC.push_back(vector<double>());
			}
		}
		string line;
		count = 0;
		if(file){
			do{
				getline(file,line);
				// find number of ref
				pos_nbref=line.find("NB_REF ");
				if( pos_nbref!=string::npos){
					istringstream text(line);
					text >> buffer_s >> buffer_ui;
				}
				// find l_sph
				pos_lsph=line.find("L_SPH ");
				if( pos_lsph!=string::npos){
					istringstream text(line);
					text >> buffer_s >> this->l_sph_ref;
				}
				// find rcut
				pos_rc=line.find("RCUT ");
				if( pos_rc!=string::npos){
					line_rc = count;
					istringstream text(line);
					text >> buffer_s >> this->rcut_ref;
				}
				// read the parameters
				if( count > line_rc+1 && count < line_rc+buffer_ui+2 ){
					//SteinhardtParams_REF_PC.push_back(new double[(this->l_sph_ref+1)]);
					AtomTypeUINTRefPC.push_back(0);
istringstream text(line);
					text >> buffer_s >> AtomTypeUINTRefPC[AtomTypeUINTRefPC.size()-1];
					if( _MySystem->getCrystal()->getIsMultisite() ){
						text >> buffer_ui2;
						this->AtomSiteRefPC[AtomTypeUINTRefPC[AtomTypeUINTRefPC.size()-1]-1].push_back(buffer_ui2);
					}
					// get norm of steinhardt params
					for(unsigned int l=0;l<this->l_sph_ref+1;l++){
						text >> buffer_d1;
						SteinhardtParams_REF_PC[AtomTypeUINTRefPC[AtomTypeUINTRefPC.size()-1]-1].push_back(buffer_d1);
					}
					if( _MySystem->getCrystal()->getIsMultisite() ){ // TODO change here and add IsCentrosymmetric to crystal database because we need the reference bond ori param only when it is not centrosymmetric and not when it is multisite
						text >> buffer_d2;
						BondOriParam_REF_PC[AtomTypeUINTRefPC[AtomTypeUINTRefPC.size()-1]-1].push_back(buffer_d2);
					}
				}
				count += 1;
			}while(file);
		}else{
			cout << "We can't read the file containing the Steinhardt parameter for perfect crystal (/data/Steinhardt/CrystalName/PerfectCrystal.dat)" << endl;
		}
		for(unsigned int i=0;i<_MySystem->getCrystal()->getNbAtomType();i++){
			cout << "Atom type : " << i+1 << endl;
			for(unsigned int j=0;j<AtomSiteRefPC[i].size();j++){
				cout << "site " << AtomSiteRefPC[i][j] << " Steinhardt ";
			       	for(unsigned int l=0;l<this->l_sph_ref+1;l++) cout << SteinhardtParams_REF_PC[i][j*(this->l_sph_ref+1)+l] << " ";
				cout << ", Bond Ori Ref " << BondOriParam_REF_PC[i][j] << endl;
			}
		}
		unsigned int lref;
		double rcref;
		if( _MySystem->getCrystal()->getIsMultisite() ){
			PC_str = "PerfectCrystal_ave";
			fullpath2data = database.c_str()+PC_str+database_extension;
			ifstream file1(fullpath2data, ios::in);
			count = 0;
			if(file1){
				do{
					getline(file1,line);
					// find number of ref
					pos_nbref=line.find("NB_REF ");
					if( pos_nbref!=string::npos){
						istringstream text(line);
						text >> buffer_s >> this->nbref;
					}
					// find l_sph
					pos_lsph=line.find("L_SPH ");
					if( pos_lsph!=string::npos){
						istringstream text(line);
						text >> buffer_s >> lref;
						if( lref != this->l_sph_ref) cerr << "WARNING !!! The spherical harmonic degrees used for the calculation of averaged Steinhard parameters is not the same as the one for not averaged Steinhardt parameters" << endl;
					}
					// find rcut
					pos_rc=line.find("RCUT ");
					if( pos_rc!=string::npos){
						line_rc = count;
						istringstream text(line);
						text >> buffer_s >> rcref;
						if( rcref != this->rcut_ref) cerr << "WARNING !!! The cutoff radius used for the calculation of averaged Steinhard parameters is not the same as the one for not averaged Steinhardt parameters" << endl;
					}
					// read the parameters
					if( count > line_rc+1 && count < line_rc+nbref+2 ){
						SteinhardtParams_REF_PC_ave.push_back(new double[(this->l_sph_ref+1)]);
						AtomTypeUINTRefPC_ave.push_back(0);
						istringstream text(line);
						text >> buffer_s >> AtomTypeUINTRefPC_ave[AtomTypeUINTRefPC_ave.size()-1];
						// get norm of steinhardt params
						for(unsigned int l=0;l<this->l_sph_ref+1;l++) text >> SteinhardtParams_REF_PC_ave[SteinhardtParams_REF_PC_ave.size()-1][l];
					}
					count += 1;
				}while(file1);
			}else{
				cout << "We can't read the file containing the averaged Steinhardt parameters for perfect crystal (/data/Steinhardt/CrystalName/PerfectCrystal_ave.dat)" << endl;
			}
		}else{
			for(unsigned int i=0;i<AtomTypeUINTRefPC.size();i++){
				SteinhardtParams_REF_PC_ave.push_back(new double[(this->l_sph_ref+1)]);
				AtomTypeUINTRefPC_ave.push_back(AtomTypeUINTRefPC[i]);
				for(unsigned int l=0;l<this->l_sph_ref+1;l++) SteinhardtParams_REF_PC_ave[i][l] = SteinhardtParams_REF_PC[i][l];
			}
		}
		// read Steinhardt param averaged over atom of same specy inside the cutoff radius
		PC_str = "PerfectCrystal_ave_cutoff";
		fullpath2data = database.c_str()+PC_str+database_extension;
		ifstream file2(fullpath2data, ios::in);
		count = 0;
		if(file2){
			do{
				getline(file2,line);
				// find number of ref
				pos_nbref=line.find("NB_REF ");
				if( pos_nbref!=string::npos){
					istringstream text(line);
					text >> buffer_s >> buffer_ui;
					if( buffer_ui != this->nbref) cerr << "WARNING !!! The number of Q vectors for the Steinhardt parameters averaged over atom of same specy inside the cutoff radius is not the same as the one for not averaged Steinhardt parameters" << endl;
				}
				// find l_sph
				pos_lsph=line.find("L_SPH ");
				if( pos_lsph!=string::npos){
					istringstream text(line);
					text >> buffer_s >> lref;
					if( lref != this->l_sph_ref) cerr << "WARNING !!! The spherical harmonic degrees used for the calculation of averaged Steinhard parameters is not the same as the one for not averaged Steinhardt parameters" << endl;
				}
				// find rcut
				pos_rc=line.find("RCUT ");
				if( pos_rc!=string::npos){
					line_rc = count;
					istringstream text(line);
					text >> buffer_s >> rcref;
					if( rcref != this->rcut_ref) cerr << "WARNING !!! The cutoff radius used for the calculation of averaged Steinhard parameters is not the same as the one for not averaged Steinhardt parameters" << endl;
				}
				// read the parameters
				if( count > line_rc+1 && count < line_rc+nbref+2 ){
					SteinhardtParams_REF_PC_ave_cutoff.push_back(new double[(this->l_sph_ref+1)]);
					AtomTypeUINTRefPC_ave.push_back(0);
					istringstream text(line);
					text >> buffer_s >> AtomTypeUINTRefPC_ave[AtomTypeUINTRefPC_ave.size()-1];
					// get norm of steinhardt params
					for(unsigned int l=0;l<this->l_sph_ref+1;l++) text >> SteinhardtParams_REF_PC_ave_cutoff[SteinhardtParams_REF_PC_ave_cutoff.size()-1][l];
				}
				count += 1;
			}while(file2);
		}else{
			cout << "We can't read the file containing the Steinhardt parameters averaged over the atom of same specy inside the cutoff radius for perfect crystal (/data/Steinhardt/CrystalName/PerfectCrystal_ave_cutoff.dat)" << endl;
		}
		// Read the defect Steinhardt parameters
		this->nbRefDef = this->Ref_Def_Names.size();
		unsigned int count_ref;
		for(unsigned int i=0;i<this->nbRefDef;i++){
			fullpath2data = database.c_str()+this->Ref_Def_Names[i]+database_extension;
			ifstream file3(fullpath2data, ios::in);
			count = 0;
			count_ref = 0;
			if(file3){
				do{
					getline(file3,line);
					// find number of ref
					pos_nbref=line.find("NB_REF ");
					if( pos_nbref!=string::npos){
						istringstream text(line);
						text >> buffer_s >> buffer_ui;
						if( buffer_ui != this->nbref) cerr << "WARNING !!! The number of Q vectors for defect : " << this->Ref_Def_Names[i] << " is not the same as the one for not averaged Steinhardt parameters" << endl;
						this->SteinhardtParams_REF_Def.push_back(new double[this->nbref*(this->l_sph_ref+1)]);
						this->AtomTypeUINTRefDef.push_back(new unsigned int[this->nbref+1]);
						this->AtomTypeUINTRefDef[i][0]=this->nbref;
					}
					// find l_sph
					pos_lsph=line.find("L_SPH ");
					if( pos_lsph!=string::npos){
						istringstream text(line);
						text >> buffer_s >> lref;
						if( lref != this->l_sph_ref) cerr << "WARNING !!! The spherical harmonic degrees used for the calculation of defect : " << this->Ref_Def_Names[i] << " Steinhardt parameters is not the same as the one for not averaged Steinhardt parameters" << endl;
					}
					// find rcut
					pos_rc=line.find("RCUT ");
					if( pos_rc!=string::npos){
						line_rc = count;
						istringstream text(line);
						text >> buffer_s >> rcref;
						if( rcref != this->rcut_ref) cerr << "WARNING !!! The cutoff radius used for the calculation of defect : " << this->Ref_Def_Names[i] << " Steinhard parameters is not the same as the one for not averaged Steinhardt parameters" << endl;
					}
					// read the parameters
					if( count > line_rc+1 && count < line_rc+nbref+2 ){
						istringstream text(line);
						text >> buffer_s >> AtomTypeUINTRefDef[i][count_ref+1];
						// get norm of steinhardt params
						for(unsigned int l=0;l<this->l_sph_ref+1;l++) text >> SteinhardtParams_REF_Def[i][count_ref*(this->l_sph_ref+1)+l];
						count_ref++;
					}
					count += 1;
				}while(file3);
			}else{
				cout << "We can't read the file containing the averaged Steinhardt parameters for perfect crystal (/data/Steinhardt/CrystalName/" << this->Ref_Def_Names[i] << ".dat)" << endl;
			}
			fullpath2data = database.c_str()+this->Ref_Def_Names[i]+ave_cutoff+database_extension;
			ifstream file4(fullpath2data, ios::in);
			count = 0;
			count_ref = 0;
			if(file4){
				do{
					getline(file4,line);
					// find number of ref
					pos_nbref=line.find("NB_REF ");
					if( pos_nbref!=string::npos){
						istringstream text(line);
						text >> buffer_s >> buffer_ui;
						if( buffer_ui != this->nbref) cerr << "WARNING !!! The number of Q vectors averaged over atom of same type inside the cutoff radius for defect : " << this->Ref_Def_Names[i] << " is not the same as the one for not averaged Steinhardt parameters" << endl;
						this->SteinhardtParams_REF_Def_ave_cutoff.push_back(new double[this->nbref*(this->l_sph_ref+1)]);
						this->AtomTypeUINTRefDef_ave_cutoff.push_back(new unsigned int[this->nbref+1]);
						this->AtomTypeUINTRefDef_ave_cutoff[i][0]=this->nbref;
					}
					// find l_sph
					pos_lsph=line.find("L_SPH ");
					if( pos_lsph!=string::npos){
						istringstream text(line);
						text >> buffer_s >> lref;
						if( lref != this->l_sph_ref) cerr << "WARNING !!! The spherical harmonic degrees used for the calculation of defect : " << this->Ref_Def_Names[i] << " averaged over atom of same type inside the cutoff radius Steinhardt parameters is not the same as the one for not averaged Steinhardt parameters" << endl;
					}
					// find rcut
					pos_rc=line.find("RCUT ");
					if( pos_rc!=string::npos){
						line_rc = count;
						istringstream text(line);
						text >> buffer_s >> rcref;
						if( rcref != this->rcut_ref) cerr << "WARNING !!! The cutoff radius used for the calculation of defect : " << this->Ref_Def_Names[i] << " Steinhard parameters averaged over atom of same type inside the cutoff radius is not the same as the one for not averaged Steinhardt parameters" << endl;
					}
					// read the parameters
					if( count > line_rc+1 && count < line_rc+nbref+2 ){
						istringstream text(line);
						text >> buffer_s >> AtomTypeUINTRefDef_ave_cutoff[i][count_ref+1];
						// get norm of steinhardt params
						for(unsigned int l=0;l<this->l_sph_ref+1;l++) text >> SteinhardtParams_REF_Def_ave_cutoff[i][count_ref*(this->l_sph_ref+1)+l];
						count_ref++;
					}
					count += 1;
				}while(file4);
			}else{
				cout << "We can't read the file containing the averaged over atom of same type inside the cutoff radius Steinhardt parameters for perfect crystal (/data/Steinhardt/CrystalName/" << this->Ref_Def_Names[i] << ".dat)" << endl;
			}
		}
		this->IsSteinhardtDatabaseRead = true;
		//for(unsigned int i=0;i<this->Ref_Def_Names.size();i++){
		//	cout << this->Ref_Def_Names[i] << endl;
		//	for(unsigned int n=0;n<nbref;n++){
		//		cout << AtomTypeUINTRefDef_ave_cutoff[i][n+1] << " ";
		//		for(unsigned int l=0;l<this->l_sph_ref+1;l++) cout << SteinhardtParams_REF_Def_ave_cutoff[i][n*(this->l_sph_ref+1)+l] << " ";
		//		cout << endl;
		//	}
		//}
		//cout << "Perfect crystal : " << endl;
		//	for(unsigned int n=0;n<nbref;n++){
		//		cout << AtomTypeUINTRefPC_ave[n] << " ";
		//		for(unsigned int l=0;l<this->l_sph_ref+1;l++) cout << SteinhardtParams_REF_PC_ave_cutoff[n][l] << " ";
		//		cout << endl;
		//	}
	}
}

double* ComputeAuxiliary::StructuralAnalysis_Steinhardt(){
	// read database
	SteinhardtDatabase_read(_MySystem->getCrystal()->getName());
	// compute Steinhardt params
	ComputeSteinhardtParameters(this->rcut_ref, this->l_sph_ref);
	const unsigned int nbAt = _MySystem->getNbAtom();
	double *D_q = new double[nbAt*(this->Ref_Def_Names.size()+1)]; // contains the Steinhardt params Euclidian distances between each atoms and perfect crystal and all defects in the database
	// compute distances
	for(unsigned int i=0;i<nbAt;i++){
		// from perfect crystal
		D_q[i*(this->Ref_Def_Names.size()+1)] = 0.;
		for(unsigned int j=0;j<AtomTypeUINTRefPC_ave.size();j++){
			if( _MySystem->getAtom(i).type_uint == AtomTypeUINTRefPC_ave[j] ){
				for(unsigned int l=0;l<this->l_sph_ref+1;l++){ // maybe exclude l=0 ?
					//D_q[i*(this->Ref_Def_Names.size()+1)] += pow(this->SteinhardtParams[i*(this->l_sph_ref+1)+l]-this->SteinhardtParams_REF_PC_ave[j][l],2.);
					D_q[i*(this->Ref_Def_Names.size()+1)] += pow(this->SteinhardtParams_ave_cutoff[i*(this->l_sph_ref+1)+l]-this->SteinhardtParams_REF_PC_ave_cutoff[j][l],2.);
				}
				D_q[i*(this->Ref_Def_Names.size()+1)] = sqrt(D_q[i*(this->Ref_Def_Names.size()+1)]);
				break;
			}
		}
		// from all defect
		for(unsigned int j=0;j<this->Ref_Def_Names.size();j++){
			D_q[i*(this->Ref_Def_Names.size()+1)+j+1] = 0.;
			//for(unsigned int k=0;k<this->AtomTypeUINTRefDef[j][0];k++){
			for(unsigned int k=0;k<this->AtomTypeUINTRefDef_ave_cutoff[j][0];k++){
				//if( _MySystem->getAtom(i).type_uint == AtomTypeUINTRefDef[j][k+1] ){
				if( _MySystem->getAtom(i).type_uint == AtomTypeUINTRefDef_ave_cutoff[j][k+1] ){
					for(unsigned int l=0;l<this->l_sph_ref+1;l++){ // maybe exclude l=0 ?
						//D_q[i*(this->Ref_Def_Names.size()+1)+j+1] += pow(this->SteinhardtParams[i*(this->l_sph_ref+1)+l]-this->SteinhardtParams_REF_Def[j][k*(this->l_sph_ref+1)+l],2.);
						D_q[i*(this->Ref_Def_Names.size()+1)+j+1] += pow(this->SteinhardtParams_ave_cutoff[i*(this->l_sph_ref+1)+l]-this->SteinhardtParams_REF_Def_ave_cutoff[j][k*(this->l_sph_ref+1)+l],2.);
					}
					D_q[i*(this->Ref_Def_Names.size()+1)+j+1] = sqrt(D_q[i*(this->Ref_Def_Names.size()+1)+j+1]);
					break;
				}
			}
		}
	}
	// compute correlation factors such as sAB = ( exp(-lambdaAB*DAx) + 2exp(-lambdaAB*DBx) ) / ( exp(-lambdaAB*DAx) + exp(-lambdaAB*DBx) )
	// 1.a. Compute lambdaAB (to begin A is the perfect crystal and B are the defects
	double *lambdaAB = new double[this->Ref_Def_Names.size()*this->nbref];
	unsigned int *type_uint_lambda = new unsigned int[this->Ref_Def_Names.size()*this->nbref];
	for(unsigned int i=0;i<this->Ref_Def_Names.size();i++){
		for(unsigned int t=0;t<this->nbref;t++){
			lambdaAB[i*this->nbref+t] = 0.;
			for(unsigned int j=0;j<this->AtomTypeUINTRefPC_ave.size();j++){
				//if( AtomTypeUINTRefPC_ave[j] == AtomTypeUINTRefDef_[i][t+1] ){
				if( AtomTypeUINTRefPC_ave[j] == AtomTypeUINTRefDef_ave_cutoff[i][t+1] ){
					for(unsigned int l=0;l<this->l_sph_ref+1;l++){ // maybe exclude l=0 ?
						//lambdaAB[i*this->nbref+t] += pow(this->SteinhardtParams_REF_PC_ave[j][l]-this->SteinhardtParams_REF_Def[i][t*(this->l_sph_ref+1)+l],2.);
						lambdaAB[i*this->nbref+t] += pow(this->SteinhardtParams_REF_PC_ave_cutoff[j][l]-this->SteinhardtParams_REF_Def_ave_cutoff[i][t*(this->l_sph_ref+1)+l],2.);
					}
					lambdaAB[i*this->nbref+t] = 2.3/sqrt(lambdaAB[i*this->nbref+t]);
					//type_uint_lambda[i*this->nbref+t] = AtomTypeUINTRefDef[i][t+1];
					type_uint_lambda[i*this->nbref+t] = AtomTypeUINTRefDef_ave_cutoff[i][t+1];
					break;
				}
			}
		}
	}
	// 1.b Compute lambdaBC (where B and C are the different defects)
	vector<vector<double>> lambdaBC; // lambdaBC[i][j*nbref+t] gives lambda between defect i and j for type t  
	vector<vector<unsigned int>> type_uint_lambdaBC;
	for(unsigned int i=0;i<this->Ref_Def_Names.size();i++){
		lambdaBC.push_back(vector<double>());
		type_uint_lambdaBC.push_back(vector<unsigned int>());
		for(unsigned int j=0;j<this->Ref_Def_Names.size();j++){
			if( i == j ){
				for(unsigned int d=0;d<this->nbref;d++){
					lambdaBC[i].push_back(0.);
					type_uint_lambdaBC[i].push_back(0);
				}
			}else if( j > i ){
				for(unsigned int d1=0;d1<this->nbref;d1++){
					lambdaBC[i].push_back(0.);
					for(unsigned int d2=0;d2<this->nbref;d2++){
						if( AtomTypeUINTRefDef_ave_cutoff[i][d1+1] == AtomTypeUINTRefDef_ave_cutoff[j][d2+1] ){
							for(unsigned int l=0;l<this->l_sph_ref+1;l++){ // maybe exclude l=0 ?
								lambdaBC[i][lambdaBC[i].size()-1] += pow(this->SteinhardtParams_REF_Def_ave_cutoff[j][d2*(this->l_sph_ref+1)+l]-this->SteinhardtParams_REF_Def_ave_cutoff[i][d1*(this->l_sph_ref+1)+l],2.);
							}
							lambdaBC[i][lambdaBC[i].size()-1] = 2.3/sqrt(lambdaBC[i][lambdaBC[i].size()-1]);
							type_uint_lambdaBC[i].push_back(AtomTypeUINTRefDef_ave_cutoff[i][d1+1]);
							break;
						}
					}
				}
			}else{
				for(unsigned int d=0;d<this->nbref;d++){
					lambdaBC[i].push_back(lambdaBC[j][i*this->nbref+d]);
					type_uint_lambdaBC[i].push_back(type_uint_lambdaBC[j][i*this->nbref+d]);
				}
			}
		}
	}
				
	// cout << "Euclidian distance from perfect crystal to =>" << endl;
	//for(unsigned int i=0;i<this->Ref_Def_Names.size();i++){
	//	cout << this->Ref_Def_Names[i] << endl;
	//	for(unsigned int j=0;j<this->nbref;j++){
	//		cout << lambdaAB[i*this->nbref+j] << " " ;
	//	}
	//	cout << endl;
	//}
	 cout << "Euclidian distance from defect " << endl;
	for(unsigned int i=0;i<this->Ref_Def_Names.size();i++){
		for(unsigned int j=0;j<this->Ref_Def_Names.size();j++){
			cout << this->Ref_Def_Names[i] << " to " << this->Ref_Def_Names[j] << endl;
			for(unsigned int t=0;t<this->nbref;t++){
				cout << lambdaBC[i][j*this->nbref+t] << " " ;
			}
			cout << "	" << endl;
		}
		cout << endl;
	}


	// 2. Compute sAB
	double *sAB = new double[nbAt*this->Ref_Def_Names.size()];
	double buffer_exp_1, buffer_exp_2;
	for(unsigned int i=0;i<nbAt;i++){
		for(unsigned int d=0;d<this->Ref_Def_Names.size();d++){
			for(unsigned int t=0;t<this->nbref;t++){
				if( _MySystem->getAtom(i).type_uint == type_uint_lambda[d*this->nbref+t] ){
					//cout << lambdaAB[d*this->nbref+t] << " " << D_q[i*(this->Ref_Def_Names.size()+1)] << " " << 
					buffer_exp_1 = exp(-lambdaAB[d*this->nbref+t]*D_q[i*(this->Ref_Def_Names.size()+1)]);
					buffer_exp_2 = exp(-lambdaAB[d*this->nbref+t]*D_q[i*(this->Ref_Def_Names.size()+1)+d+1]);
					sAB[i*this->Ref_Def_Names.size()+d] = ( buffer_exp_1 + 2*buffer_exp_2 ) / ( buffer_exp_1 + buffer_exp_2 );
					//cout << sAB[i*this->Ref_Def_Names.size()+d] << " ";
					break;
				}
			}
		}
		//cout << endl;
	}
	// 3. when all sAB < 1.5 => we attribute the perfect crystal, if just one sAB > 1.5 we attribute the associated defect, when more than one sAB > 1.5 we compute sBC with B and C all other defect with sAB > 1.5
	double correl_fac = 1.5;
	double buffer_d;
	double *Struct = new double[nbAt*2]; // first column is integer which corresponds to the identified structure (i.e. 0 => perfect crystal (sAB < 1.5), i => index of defect (printed in DefectIndex.txt file)), the second column corresponds to the "error" in the identification (i.e. sum(sAB-1)/nbDef when perfect crystal is identified and ? when a defect is identifyed)
	vector<unsigned int> *DefCor = new vector<unsigned int>[nbAt];
	for(unsigned int i=0;i<nbAt;i++){
		for(unsigned int d=0;d<this->Ref_Def_Names.size();d++) if( sAB[i*this->Ref_Def_Names.size()+d] > correl_fac ) DefCor[i].push_back(d);
		if( DefCor[i].size() == 0 ){ // perfect crystal case
			Struct[i*2] = 0.;
			Struct[i*2+1] = 0.;
			for(unsigned int d=0;d<this->Ref_Def_Names.size();d++) Struct[i*2+1] += sAB[i*this->Ref_Def_Names.size()+d]-1;
			Struct[i*2+1] /= this->Ref_Def_Names.size();
		}else if( DefCor[i].size() == 1 ){ // defect cases
			Struct[i*2] = DefCor[i][0]+1;
			Struct[i*2+1] = 2.-sAB[i*this->Ref_Def_Names.size()+DefCor[i][0]];
		}else{
			while( DefCor[i].size() > 1 ){
				for(unsigned int j2=1;j2<DefCor[i].size();j2++){
					for(unsigned int t=0;t<this->nbref;t++){
						if( _MySystem->getAtom(i).type_uint == type_uint_lambdaBC[DefCor[i][0]][DefCor[i][j2]*this->nbref+t] ){
							buffer_exp_1 = exp(-lambdaBC[DefCor[i][0]][DefCor[i][j2]*this->nbref+t]*D_q[i*(this->Ref_Def_Names.size()+1)+DefCor[i][0]+1]);
							buffer_exp_2 = exp(-lambdaBC[DefCor[i][0]][DefCor[i][j2]*this->nbref+t]*D_q[i*(this->Ref_Def_Names.size()+1)+DefCor[i][j2]+1]);
							buffer_d = ( buffer_exp_1 + 2*buffer_exp_2 ) / ( buffer_exp_1 + buffer_exp_2 );
							if( buffer_d < correl_fac ){ // the defect correls better with DefCor[i][0] => remove DefCor[j2]
								DefCor[i].erase(DefCor[i].begin()+j2);
								Struct[i*2+1] = buffer_d-1.;
							}else{
								DefCor[i].erase(DefCor[i].begin());
								Struct[i*2+1] = 2.-buffer_d;
							}
							break;
						}
					}
					if( DefCor[i].size() == 1 ) break;
				}
			}
			Struct[i*2] = DefCor[i][0]+1;
		}
	}

	// write the DefectIndex.txt file
	ofstream writefile("DefectIndex.txt");
	writefile << "DefectIndex used for structural analysis" << endl;
	writefile << "0 PerfectCrystal" << endl;
	for(unsigned int i=0;i<this->Ref_Def_Names.size();i++) writefile << i+1 << " " << this->Ref_Def_Names[i] << endl;
		
	delete[] D_q;
	delete[] lambdaAB;
	delete[] type_uint_lambda;

	return Struct;	

}

string ComputeAuxiliary::getSteinhardtDatabase(string CrystalName){
	string database;	
	string backslash="/";
	if( CrystalName == "Forsterite" ){
		#ifdef STEINHARDT_FORSTERITE_DATABASE
		database = STEINHARDT_FORSTERITE_DATABASE;
		#endif
	}else if( CrystalName == "Periclase" ){
		#ifdef STEINHARDT_PERICLASE_DATABASE
		database = STEINHARDT_PERICLASE_DATABASE;
		#endif
	}
	if( database.empty() ){
		cerr << "Warning database environment for crystal is empty" << endl;
		exit(EXIT_FAILURE);
	} else {
		return database.c_str()+backslash;
	}

}

void ComputeAuxiliary::read_params(){
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
	}else{
		cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
		exit(EXIT_FAILURE);
	}
}

ComputeAuxiliary::~ComputeAuxiliary(){
	delete MT;
	if( IsBondOriParam ){
	cout << 1 << endl;
		delete[] BondOriParam;
	cout << 1 << endl;
		delete[] Atom_SiteIndex;
	cout << 1 << endl;
	}
	if( IsStrainTensor ){
		delete[] StrainTensor;
	}
	if( IsStrainInvII ){
		delete[] Strain_invII;
	}
}

