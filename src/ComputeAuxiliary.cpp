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
void ComputeAuxiliary::ComputeSteinhardtParameters_Mono(const double rc, const int l_sph){
	// if neighbours have not been searched perform the research
	if( !_MySystem->getIsNeighbours() || _MySystem->get_current_rc() != rc ){
		_MySystem->searchNeighbours(rc);
	}
	cout << "Computing Steinhart parameters using only ions of identical type (mono).. " << endl;
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	this->Malpha = new unsigned int[nbAt*(nbNMax+1)]; // array containing the index of neighbours of the same species (or same site in case of multisite crystal) with the first line corresponding to the number of neighbours, i.e. Malpha[i*(nbNMax+1)] = nb of neighbour of atom i, Malpha[i*(nbNMax+1)+j+1] = id of the jth neighbour of atom i
	this->Qlm = new complex<double>[nbAt*(l_sph*2+1)*(l_sph+1)]; // complex array containing the spherical harmonic for the different modes Qlm[i*(l_sph*2+1)*(l_sph+1)+l*(l_sph*2+1)+m] gives the spherical harmonic for atom i and degree l and m
	unsigned int lsph2 = (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1 = l_sph*2+1;
	this->SteinhardtParams = new double[nbAt*(l_sph+1)];
	for(unsigned int i=0;i<nbAt*(l_sph*2+1);i++){
		for(unsigned int j=0;j<l_sph+1;j++) Qlm[i*j+j] = (0.,0.);
	}
	double zeronum = 1e-8;
	const int bar_length = 30;
	double prog=0;
	double xpos,ypos,zpos,xp,yp,zp, colat, longit;
	// loop on all atoms and neighbours to compute Qalpha and store neighbour of the same species
	// Here is the most time consuming loop of the function, use parallel computation
	unsigned int j_loop, l_loop_st, l_loop_st2, id;
	int m_loop_st, m_loop_st2, m_loop_st3;
	unsigned int count_t = 0;
	#pragma omp parallel for private(xpos,ypos,zpos,j_loop,id,xp,yp,zp,colat,longit,l_loop_st,m_loop_st,l_loop_st2,m_loop_st2,m_loop_st3)
	for(unsigned int i=0;i<nbAt;i++){
		if( omp_get_thread_num() == 0 ){
			prog = double(count_t*omp_get_num_threads())/double(nbAt);
			cout << "\r[" << string(bar_length*prog,'X') << string(bar_length*(1-prog),'-') << "] " << setprecision(3) << 100*prog << "%";
			count_t++;
		}
		xpos = _MySystem->getWrappedPos(i).x;
		ypos = _MySystem->getWrappedPos(i).y;
		zpos = _MySystem->getWrappedPos(i).z;
		Malpha[i*(nbNMax+1)] = 0; 
		for(j_loop=0;j_loop<_MySystem->getNeighbours(i*(nbNMax+1));j_loop++){
			id = _MySystem->getNeighbours(i*(nbNMax+1)+j_loop+1);
			if( _MySystem->getAtom(i).type_uint == _MySystem->getAtom(id).type_uint ){
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
				for(l_loop_st=0;l_loop_st<l_sph+1;l_loop_st++){
					for(m_loop_st=-l_loop_st;m_loop_st<(int) l_loop_st+1;m_loop_st++){
						Qlm[i*lsph2 + l_loop_st*lsph1 + (unsigned int) (m_loop_st + (int) l_sph)] += spherical_harmonics(l_loop_st, m_loop_st, colat, longit);
					}
				}
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
	}
	cout << "Done !" << endl;	
	this->IsSteinhardt_Mono = true;
}

// Same as above but with the average version where the Steinhardt params are averaged using all neighbours (in the case above only the ions of same type inside the cutoff radius are used)
void ComputeAuxiliary::ComputeSteinhardtParameters_Multi(const double rc, const int l_sph){
	// if neighbours have not been searched perform the research
	if( !_MySystem->getIsNeighbours() || _MySystem->get_current_rc() != rc ){
		_MySystem->searchNeighbours(rc);
	}
	cout << "Computing Steinhart parameters using all ion type (multi).. ";
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	this->Malpha = new unsigned int[nbAt*(nbNMax+1)]; // array containing the index of neighbours of the same species (or same site in case of multisite crystal) with the first line corresponding to the number of neighbours, i.e. Malpha[i*(nbNMax+1)] = nb of neighbour of atom i, Malpha[i*(nbNMax+1)+j+1] = id of the jth neighbour of atom i
	this->Qlm = new complex<double>[nbAt*(l_sph*2+1)*(l_sph+1)]; // complex array containing the spherical harmonic for the different modes Qlm[i*(l_sph*2+1)*(l_sph+1)+l*(l_sph*2+1)+m] gives the spherical harmonic for atom i and degree l and m
	unsigned int lsph2 = (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1 = l_sph*2+1;
	this->SteinhardtParams = new double[nbAt*(l_sph+1)];
	for(unsigned int i=0;i<nbAt*(l_sph*2+1);i++){
		for(unsigned int j=0;j<l_sph+1;j++) Qlm[i*j+j] = (0.,0.);
	}
	double zeronum = 1e-8;
	const int bar_length = 30;
	double prog=0;
	double xpos,ypos,zpos,xp,yp,zp, colat, longit;
	// loop on all atoms and neighbours to compute Qalpha and store neighbour of the same species
	// Here is the most time consuming loop of the function, use parallel computation
	unsigned int j_loop, l_loop_st, l_loop_st2, id;
	int m_loop_st, m_loop_st2, m_loop_st3;
	unsigned int count_t = 0;
	#pragma omp parallel for private(xpos,ypos,zpos,j_loop,id,xp,yp,zp,colat,longit,l_loop_st,m_loop_st,l_loop_st2,m_loop_st2,m_loop_st3)
	for(unsigned int i=0;i<nbAt;i++){
		if( omp_get_thread_num() == 0 ){
			prog = double(count_t*omp_get_num_threads())/double(nbAt);
			cout << "\r[" << string(bar_length*prog,'X') << string(bar_length*(1-prog),'-') << "] " << setprecision(3) << 100*prog << "%";
			count_t++;
		}
		xpos = _MySystem->getWrappedPos(i).x;
		ypos = _MySystem->getWrappedPos(i).y;
		zpos = _MySystem->getWrappedPos(i).z;
		Malpha[i*(nbNMax+1)] = 0; 
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
	}
	cout << "Done !" << endl;
	this->IsSteinhardt_Multi = true;
}

// Compute Steinhardt parameters but filter the neighbour used for the calculation with their types (i.e. length of vector equal l_sph*(nbAtomType+1))
void ComputeAuxiliary::ComputeSteinhardtParameters_FilteredNeigh(const double rc, const int l_sph){
	// if neighbours have not been searched perform the research
	if( !_MySystem->getIsNeighbours() || _MySystem->get_current_rc() != rc ){
		_MySystem->searchNeighbours(rc);
	}
	cout << "Computing filtered Steinhart parameters (l*(n+1) dimension vector corresponding to the degree l considering all neighboring ion types (n=0) and ion of types n).. ";
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	const unsigned int kp_max = _MySystem->getCrystal()->getNbAtomType()+1;

	this->Qlm = new complex<double>[kp_max*nbAt*(l_sph*2+1)*(l_sph+1)]; // complex array containing the spherical harmonic for the different modes Qlm[i*(l_sph*2+1)*(l_sph+1)+l*(l_sph*2+1)+m] gives the spherical harmonic for atom i and degree l and m
	complex<double> *buffer_Qlm = new complex<double>[nbAt*(l_sph*2+1)*(l_sph+1)]; // complex array containing the spherical harmonic for the different modes Qlm[i*(l_sph*2+1)*(l_sph+1)+l*(l_sph*2+1)+m] gives the spherical harmonic for atom i and degree l and m
	unsigned int lsph2 = (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1 = l_sph*2+1;
	this->SteinhardtParams = new double[kp_max*nbAt*(l_sph+1)];
	for(unsigned int i=0;i<nbAt*(l_sph*2+1);i++){
		for(unsigned int j=0;j<l_sph+1;j++) Qlm[i*j+j] = (0.,0.);
	}
	this->nbNeigh_FiltNeigh = new unsigned int[kp_max*nbAt];
	double zeronum = 1e-8;
	const int bar_length = 30;
	double prog=0;
	double xpos,ypos,zpos,xp,yp,zp, colat, longit;
	// loop on all atoms and neighbours to compute Qalpha and store neighbour of the same species
	// Here is the most time consuming loop of the function, use parallel computation
	unsigned int j_loop, l_loop_st, l_loop_st1, l_loop_st2, l_loop_st3, current_type, kp_loop, nbIni, id;
	int m_loop_st, m_loop_st1, m_loop_st2, m_loop_st3;
	unsigned int count_t = 0;
	#pragma omp parallel for private(xpos,ypos,zpos,j_loop,id,xp,yp,zp,colat,longit,l_loop_st,m_loop_st,l_loop_st1,l_loop_st2,l_loop_st3,m_loop_st2,m_loop_st3,m_loop_st1, current_type, kp_loop,nbIni)
	for(unsigned int i=0;i<nbAt;i++){
		if( omp_get_thread_num() == 0 ){
			prog = double(count_t*omp_get_num_threads())/double(nbAt);
			cout << "\r[" << string(bar_length*prog,'X') << string(bar_length*(1-prog),'-') << "] " << setprecision(3) << 100*prog << "%";
			count_t++;
		}
		xpos = _MySystem->getWrappedPos(i).x;
		ypos = _MySystem->getWrappedPos(i).y;
		zpos = _MySystem->getWrappedPos(i).z;
		for(nbIni=1;nbIni<kp_max;nbIni++) nbNeigh_FiltNeigh[i*kp_max+nbIni] = 0;
		nbNeigh_FiltNeigh[i*kp_max] = _MySystem->getNeighbours(i*(nbNMax+1));
		for(j_loop=0;j_loop<nbNeigh_FiltNeigh[i*kp_max];j_loop++){
			
			id = _MySystem->getNeighbours(i*(nbNMax+1)+j_loop+1);
			current_type = _MySystem->getAtom(id).type_uint;

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
			for(l_loop_st=0;l_loop_st<l_sph+1;l_loop_st++){
				for(m_loop_st=-l_loop_st;m_loop_st<(int) l_loop_st+1;m_loop_st++) buffer_Qlm[i*lsph2 + l_loop_st*lsph1 + (unsigned int) (m_loop_st + (int) l_sph)] = spherical_harmonics(l_loop_st, m_loop_st, colat, longit);
			}
			
			// store for kp=0 (i.e. all neighbour types considered)
			for(l_loop_st1=0;l_loop_st1<l_sph+1;l_loop_st1++){
				for(m_loop_st1=-l_loop_st1;m_loop_st1<(int) l_loop_st1+1;m_loop_st1++){
					Qlm[i*lsph2 + l_loop_st1*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += buffer_Qlm[i*lsph2 + l_loop_st1*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)];
				}
			}

			// continue with the other atom types
			nbNeigh_FiltNeigh[i*kp_max+current_type] += 1;
			for(l_loop_st2=0;l_loop_st2<l_sph+1;l_loop_st2++){
				for(m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
					Qlm[current_type*nbAt*lsph2 + i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)] += buffer_Qlm[i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)];
				}
			}
		}
		for(kp_loop=0;kp_loop<kp_max;kp_loop++){
			for(l_loop_st3=0;l_loop_st3<l_sph+1;l_loop_st3++){
				SteinhardtParams[kp_loop*(l_sph+1)*nbAt+i*(l_sph+1)+l_loop_st3] = 0.; 
				for(m_loop_st3=-l_loop_st3;m_loop_st3<(int) l_loop_st3+1;m_loop_st3++){
					SteinhardtParams[kp_loop*(l_sph+1)*nbAt+i*(l_sph+1)+l_loop_st3] += norm(Qlm[kp_loop*lsph2*nbAt + i*lsph2 + l_loop_st3*lsph1 + (unsigned int) (m_loop_st3 + (int) l_sph)]);
				}
				SteinhardtParams[kp_loop*(l_sph+1)*nbAt+i*(l_sph+1)+l_loop_st3] *= 4.*M_PI/(2.*l_loop_st3+1.); 
				SteinhardtParams[kp_loop*(l_sph+1)*nbAt+i*(l_sph+1)+l_loop_st3] /= pow(nbNeigh_FiltNeigh[i*kp_max+kp_loop],2.);
				SteinhardtParams[kp_loop*(l_sph+1)*nbAt+i*(l_sph+1)+l_loop_st3] = sqrt(SteinhardtParams[kp_loop*(l_sph+1)*nbAt+i*(l_sph+1)+l_loop_st3]);
			}
		}
	}
	delete[] buffer_Qlm;
	cout << "Done !" << endl;
	this->IsSteinhardt_Filtered = true;
}

void ComputeAuxiliary::AverageSteinhardtParameters_Mono(const double rc, const int l_sph){
	cout << "Averaging Steinhardt parameters using ions of the same type (mono).." << endl;
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	const unsigned int lsph2 = (l_sph+1)*(l_sph*2+1.);
	const unsigned int lsph1 = l_sph*2+1;
	this->SteinhardtParams_ave = new double[nbAt*(l_sph+1)];
	vector<complex<double>> buffer_complex(lsph2);
	unsigned int l_loop_st2, neigh;
	int m_loop_st0, m_loop_st1, m_loop_st2;
	
	#pragma omp parallel for private(l_loop_st2,m_loop_st0,neigh,m_loop_st1,m_loop_st2) firstprivate(buffer_complex)
	for(unsigned int i=0;i<nbAt;i++){
		for(l_loop_st2=0;l_loop_st2<l_sph+1;l_loop_st2++){
			SteinhardtParams_ave[i*(l_sph+1)+l_loop_st2] = 0.; 
			for(m_loop_st0=-l_loop_st2;m_loop_st0<(int) l_loop_st2+1;m_loop_st0++){
				buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] = Qlm[i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] / (double) _MySystem->getNeighbours(i*(nbNMax+1));
			}
			for(neigh=0;neigh<Malpha[i*(nbNMax+1)];neigh++){
			       for(m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++){
				       buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[Malpha[i*(nbNMax+1)+neigh+1]*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) _MySystem->getNeighbours(Malpha[i*(nbNMax+1)+neigh+1]*(nbNMax+1)));
			       }
			}
			for(m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
				SteinhardtParams_ave[i*(l_sph+1)+l_loop_st2] += norm(buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
			}
			SteinhardtParams_ave[i*(l_sph+1)+l_loop_st2] *= 4.*M_PI/(2.*l_loop_st2+1.); 
			SteinhardtParams_ave[i*(l_sph+1)+l_loop_st2] /= pow(1+Malpha[i*(nbNMax+1)],2.);
			SteinhardtParams_ave[i*(l_sph+1)+l_loop_st2] = sqrt(SteinhardtParams_ave[i*(l_sph+1)+l_loop_st2]);
		}
	}
	cout << "Done !" << endl;
	this->IsSteinhardt_AveMono = true;
}

void ComputeAuxiliary::AverageSteinhardtParameters_Multi(const double rc, const int l_sph){
	cout << "Averaging Steinhardt parameters using all ion types (multi).." << endl;
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	unsigned int lsph2 = (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1 = l_sph*2+1;
	this->SteinhardtParams_ave = new double[nbAt*(l_sph+1)];
	vector<complex<double>> buffer_complex(lsph2);
	unsigned int l_loop_st2, neigh;
	int m_loop_st0, m_loop_st1, m_loop_st2;
	
	#pragma omp parallel for private(l_loop_st2,m_loop_st0,neigh,m_loop_st1,m_loop_st2) firstprivate(buffer_complex)
	for(unsigned int i=0;i<nbAt;i++){
		for(l_loop_st2=0;l_loop_st2<l_sph+1;l_loop_st2++){
			SteinhardtParams_ave[i*(l_sph+1)+l_loop_st2] = 0.; 
			for(m_loop_st0=-l_loop_st2;m_loop_st0<(int) l_loop_st2+1;m_loop_st0++){
				buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] = Qlm[i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] / (double) _MySystem->getNeighbours(i*(nbNMax+1));
			}
			for(neigh=0;neigh<_MySystem->getNeighbours(i*(nbNMax+1));neigh++){
			       for(m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++){
				       buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) _MySystem->getNeighbours(_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)*(nbNMax+1)));
			       }
			}
			for(m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
				SteinhardtParams_ave[i*(l_sph+1)+l_loop_st2] += norm(buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
			}
			SteinhardtParams_ave[i*(l_sph+1)+l_loop_st2] *= 4.*M_PI/((2.*l_loop_st2)+1.); 
			SteinhardtParams_ave[i*(l_sph+1)+l_loop_st2] /= pow(1+_MySystem->getNeighbours(i*(nbNMax+1)),2.);
			SteinhardtParams_ave[i*(l_sph+1)+l_loop_st2] = sqrt(SteinhardtParams_ave[i*(l_sph+1)+l_loop_st2]);
		}
	}
	cout << "Done !" << endl;
	this->IsSteinhardt_AveMulti = true;
}

void ComputeAuxiliary::AverageSteinhardtParameters_Filtered(const double rc, const int l_sph){
	cout << "Averaging Steinhardt parameters filtering ion types.." << endl;
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	unsigned int lsph2 = (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1 = l_sph*2+1;
	const unsigned int kp_max = _MySystem->getCrystal()->getNbAtomType()+1;
	this->SteinhardtParams_ave = new double[nbAt*(l_sph+1)*kp_max]; // averaged using all neighbors, Si neigh, O neigh and Mg neigh
	vector<complex<double>> buffer_complex(lsph2);
	unsigned int l_loop_st2, neigh;
	int m_loop_st0, m_loop_st1, m_loop_st2;
	unsigned int nbneightype = 0;

	// Average steinhardt parameters by filtering also the neighbours during the averaging procedure 
	#pragma omp parallel for private(l_loop_st2,m_loop_st0,neigh,m_loop_st1,m_loop_st2,nbneightype) firstprivate(buffer_complex)
	for(unsigned int i=0;i<nbAt;i++){
		for(unsigned int kp=0;kp<kp_max;kp++){
			for(l_loop_st2=0;l_loop_st2<l_sph+1;l_loop_st2++){
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] = 0.; 
				for(m_loop_st0=-l_loop_st2;m_loop_st0<(int) l_loop_st2+1;m_loop_st0++){
					buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] = Qlm[i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] / (double) Malpha[i*(nbNMax+1)];
				}
				nbneightype = 0;
				for(neigh=0;neigh<_MySystem->getNeighbours(i*(nbNMax+1));neigh++){
					if( kp > 0 ){
						if( _MySystem->getAtom(_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)).type_uint == kp ){
							nbneightype += 1;
				       			for(m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++){
				       			        buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) Malpha[_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)*(nbNMax+1)]);
				       			}
						}
					}else{
				       		for(m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++){
				       		        buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) Malpha[_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)*(nbNMax+1)]);
				       		}
					}
				}
				for(m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
					SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] += norm(buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
				}
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] *= 4.*M_PI/((2.*l_loop_st2)+1.); 
				if( kp == 0 ) nbneightype = _MySystem->getNeighbours(i*(nbNMax+1));
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] /= pow(1+nbneightype,2.);
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] = sqrt(SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2]);
			}
		}
	}	
	cout << "Done !" << endl;
	this->IsSteinhardt_AveFiltered = true;
}

void ComputeAuxiliary::AverageFilteredSteinhardtParameters_Filtered(const double rc, const int l_sph){
	cout << "Averaging filtered Steinhardt parameters filtering ion types.." << endl;
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	unsigned int lsph2 = (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1 = l_sph*2+1;
	const unsigned int kp_max = _MySystem->getCrystal()->getNbAtomType()+1;
	this->SteinhardtParams_ave = new double[nbAt*(l_sph+1)*kp_max]; // averaged using all neighbors, Si neigh, O neigh and Mg neigh
	vector<complex<double>> buffer_complex(lsph2*kp_max);
	unsigned int id, current_type, l_loop_st2, neigh;
	int m_loop_st0, m_loop_st1, m_loop_st2;

	// Average steinhardt parameters by filtering also the neighbours during the averaging procedure 
	#pragma omp parallel for private(l_loop_st2,m_loop_st0,neigh,m_loop_st1,m_loop_st2,id,current_type) firstprivate(buffer_complex)
	for(unsigned int i=0;i<nbAt;i++){
		for(l_loop_st2=0;l_loop_st2<l_sph+1;l_loop_st2++){
			for(unsigned int kp=0;kp<kp_max;kp++){
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] = 0.; 
				for(m_loop_st0=-l_loop_st2;m_loop_st0<(int) l_loop_st2+1;m_loop_st0++){
					buffer_complex[kp*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] = Qlm[kp*lsph2*nbAt + i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] / (double) nbNeigh_FiltNeigh[i*kp_max+kp];
				}
			}
			for(neigh=0;neigh<_MySystem->getNeighbours(i*(nbNMax+1));neigh++){
				id = _MySystem->getNeighbours(i*(nbNMax+1)+neigh+1);
				current_type = _MySystem->getAtom(id).type_uint;
			       	for(m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++){
			       	        buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[id*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) nbNeigh_FiltNeigh[kp_max*id]); // stored for k=0 (all ion types)
			       	        buffer_complex[current_type*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[current_type*lsph2*nbAt + id*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) nbNeigh_FiltNeigh[current_type+kp_max*id]); // store for ions of type current type)
			       	}
			}
			for(unsigned int kp=0;kp<kp_max;kp++){
				for(m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
					SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] += norm(buffer_complex[kp*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
				}
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] *= 4.*M_PI/((2.*l_loop_st2)+1.); 
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] /= pow(1+nbNeigh_FiltNeigh[i*kp_max+kp],2.);
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] = sqrt(SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2]);
			}
		}
	}
	cout << "Done !" << endl;
	this->IsSteinhardt_FilteredAveFiltered = true;
}

void ComputeAuxiliary::AverageFilteredSteinhardtParameters_Mono(const double rc, const int l_sph){
	cout << "Averaging filtered Steinhardt parameters using ion of same type.." << endl;
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	unsigned int lsph2 = (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1 = l_sph*2+1;
	const unsigned int kp_max = _MySystem->getCrystal()->getNbAtomType()+1;
	this->SteinhardtParams_ave = new double[nbAt*(l_sph+1)*kp_max]; // averaged using all neighbors, Si neigh, O neigh and Mg neigh
	vector<complex<double>> buffer_complex(lsph2*kp_max);
	unsigned int id, current_type, l_loop_st2, neigh;
	int m_loop_st0, m_loop_st1, m_loop_st2;

	// Average steinhardt parameters by filtering also the neighbours during the averaging procedure 
	#pragma omp parallel for private(l_loop_st2,m_loop_st0,neigh,m_loop_st1,m_loop_st2,id,current_type) firstprivate(buffer_complex)
	for(unsigned int i=0;i<nbAt;i++){
		current_type = _MySystem->getAtom(i).type_uint;
		for(l_loop_st2=0;l_loop_st2<l_sph+1;l_loop_st2++){
			for(unsigned int kp=0;kp<kp_max;kp++){
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] = 0.; 
				for(m_loop_st0=-l_loop_st2;m_loop_st0<(int) l_loop_st2+1;m_loop_st0++){
					buffer_complex[kp*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] = Qlm[kp*lsph2*nbAt + i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] / (double) nbNeigh_FiltNeigh[i*kp_max+current_type];
				}
			}
			for(neigh=0;neigh<_MySystem->getNeighbours(i*(nbNMax+1));neigh++){
				id = _MySystem->getNeighbours(i*(nbNMax+1)+neigh+1);
				if( _MySystem->getAtom(id).type_uint == current_type ){
					for(unsigned int kp=0;kp<kp_max;kp++){
			       			for(m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++){
			       			        buffer_complex[kp*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[kp*lsph2*nbAt + id*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) nbNeigh_FiltNeigh[current_type+kp_max*id]);
			       			}
					}
				}
			}
			for(unsigned int kp=0;kp<kp_max;kp++){
				for(m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
					SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] += norm(buffer_complex[kp*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
				}
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] *= 4.*M_PI/((2.*l_loop_st2)+1.); 
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] /= pow(1+nbNeigh_FiltNeigh[i*kp_max+current_type],2.);
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] = sqrt(SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2]);
			}
		}
	}
	cout << "Done !" << endl;
	this->IsSteinhardt_AveFilteredMono = true;
}

void ComputeAuxiliary::AverageFilteredSteinhardtParameters_Multi(const double rc, const int l_sph){
	cout << "Averaging filtered Steinhardt parameters using all ion types.." << endl;
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	unsigned int lsph2 = (l_sph+1)*(l_sph*2+1.);
	unsigned int lsph1 = l_sph*2+1;
	const unsigned int kp_max = _MySystem->getCrystal()->getNbAtomType()+1;
	this->SteinhardtParams_ave = new double[nbAt*(l_sph+1)*kp_max]; // averaged using all neighbors, Si neigh, O neigh and Mg neigh
	vector<complex<double>> buffer_complex(lsph2*kp_max);
	unsigned int id, l_loop_st2, neigh;
	int m_loop_st0, m_loop_st1, m_loop_st2;

	// Average steinhardt parameters by filtering also the neighbours during the averaging procedure 
	#pragma omp parallel for private(l_loop_st2,m_loop_st0,neigh,m_loop_st1,m_loop_st2,id) firstprivate(buffer_complex)
	for(unsigned int i=0;i<nbAt;i++){
		for(l_loop_st2=0;l_loop_st2<l_sph+1;l_loop_st2++){
			for(unsigned int kp=0;kp<kp_max;kp++){
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] = 0.; 
				for(m_loop_st0=-l_loop_st2;m_loop_st0<(int) l_loop_st2+1;m_loop_st0++){
					buffer_complex[kp*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] = Qlm[kp*lsph2*nbAt + i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] / (double) nbNeigh_FiltNeigh[i*kp_max];
				}
				for(neigh=0;neigh<_MySystem->getNeighbours(i*(nbNMax+1));neigh++){
					id = _MySystem->getNeighbours(i*(nbNMax+1)+neigh+1);
				       	for(m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++){
				       	        buffer_complex[kp*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[kp*lsph2*nbAt + id*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) nbNeigh_FiltNeigh[kp_max*id]);
				       	}
				}
				for(m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
					SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] += norm(buffer_complex[kp*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
				}
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] *= 4.*M_PI/((2.*l_loop_st2)+1.); 
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] /= pow(1+nbNeigh_FiltNeigh[i*kp_max],2.);
				SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2] = sqrt(SteinhardtParams_ave[kp*nbAt*(l_sph+1)+i*(l_sph+1)+l_loop_st2]);
			}
		}
	}
	cout << "Done !" << endl;
	this->IsSteinhardt_AveFilteredMulti = true;
}

double* ComputeAuxiliary::ComputeSteinhardtParameters(const double rc, const int l_sph, std::string SteinhardtStyle, std::string AveStyle){
	bool filtered;
	// Steinhardt style
	if( SteinhardtStyle == "Mono" ){
		if( !IsSteinhardt_Mono ) ComputeSteinhardtParameters_Mono(rc, l_sph);
		filtered = false;
	}else if( SteinhardtStyle == "Multi" ){
		if( !IsSteinhardt_Multi ) ComputeSteinhardtParameters_Multi(rc, l_sph);
		filtered = false;
	}else if( SteinhardtStyle == "Filtered" ){
		if( !IsSteinhardt_Filtered ) ComputeSteinhardtParameters_FilteredNeigh(rc, l_sph);
		filtered = true;
	}else{
		cerr << "The Steinhardt parameter style to print has not been recognized" << endl;
	       	cerr << "Possible styles : " << endl;
		cerr << "Mono : Classic Steinhardt parameters using ions ofthe same species" << endl;
		cerr << "Multi : Classic Steinhardt parameters using all ion types" << endl;
		cerr << "Filtered : Steinhardt parameters filtered during their computation using ion type (dim = l*(nbAtomType+1))" << endl;
		return 0;
	}

	// Averaging
	if( AveStyle == "none" ){
		return this->SteinhardtParams;
	}else if( AveStyle == "Mono" ){
		if( !filtered && !IsSteinhardt_AveMono ) AverageSteinhardtParameters_Mono(rc,l_sph);
		else if( filtered && !IsSteinhardt_AveFilteredMono ) AverageFilteredSteinhardtParameters_Mono(rc,l_sph);
		return this->SteinhardtParams_ave;
	}else if( AveStyle == "Multi" ){
		if( !filtered && !IsSteinhardt_AveMulti ) AverageSteinhardtParameters_Multi(rc,l_sph);
		else if( filtered && !IsSteinhardt_AveFilteredMulti ) AverageFilteredSteinhardtParameters_Multi(rc,l_sph);
		return this->SteinhardtParams_ave;
	}else if( AveStyle == "Filtered" ){
		if( ( SteinhardtStyle == "Mono" || SteinhardtStyle == "Multi" ) && !IsSteinhardt_AveFiltered ) AverageSteinhardtParameters_Filtered(rc,l_sph);
		else if( SteinhardtStyle == "Filtered" && !IsSteinhardt_FilteredAveFiltered ) AverageFilteredSteinhardtParameters_Filtered(rc,l_sph); 
		return this->SteinhardtParams_ave;
	}else{
		cerr << "The style of average for the Steinhardt parameters to print has not been recognized" << endl;
	       	cerr << "Possible styles : " << endl;
		cerr << "none : no averaging" << endl;
		cerr << "Mono : averaging over ions ofthe same species" << endl;
		cerr << "Multi : averaging over all ion types" << endl;
		cerr << "Filtered : filter during the averaging using ion type (dim = l*(nbAtomType+1))" << endl;
		return 0;
	}
}

void ComputeAuxiliary::PrintSteinhardtParam(vector<unsigned int> At_index, string ext_filename, string SteinhardtStyle, string AveStyle){
	int l_sph = _MySystem->get_lsph();
	double rc = _MySystem->get_rcut();

	ComputeSteinhardtParameters(rc, l_sph, SteinhardtStyle, AveStyle);
	
	string at_type, filename;
	unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int kp_max = _MySystem->getCrystal()->getNbAtomType();
	bool filtered;
	unsigned int dim_s;
	if( SteinhardtStyle == "Filtered" || AveStyle == "Filtered" ) dim_s = (l_sph+1)*kp_max;
	else dim_s = l_sph+1;

	if( AveStyle == "none" ){
		for(unsigned int t=0;t<_MySystem->getNbAtomType();t++){
			at_type = _MySystem->getAtomType(t);
			filename = at_type+"_Steinhardt_"+ext_filename+".dat";
			ofstream writefile(filename);
			for(unsigned int i=0;i<At_index.size();i++){
				if(_MySystem->getAtom(At_index[i]).type_uint != t+1) continue;
				else{
					for(unsigned int d=0;d<dim_s;d++) writefile << dim_s << " " << SteinhardtParams[(d/(l_sph+1))*nbAt*(l_sph+1)+At_index[i]*(l_sph+1)+(d%(l_sph+1))] << " " << At_index[i] << endl;
				}
			}
			writefile.close();
		}
	}else{
		for(unsigned int t=0;t<_MySystem->getNbAtomType();t++){
			at_type = _MySystem->getAtomType(t);
			filename = at_type+"_Steinhardt_"+ext_filename+".dat";
			ofstream writefile(filename);
			for(unsigned int i=0;i<At_index.size();i++){
				if(_MySystem->getAtom(At_index[i]).type_uint != t+1) continue;
				else{
					for(unsigned int d=0;d<dim_s;d++) writefile << dim_s << " " << SteinhardtParams_ave[(d/(l_sph+1))*nbAt*(l_sph+1)+At_index[i]*(l_sph+1)+(d%(l_sph+1))] << " " << At_index[i] << endl;
				}
			}
			writefile.close();
		}
	}
}

void ComputeAuxiliary::ComputeSteinhardtParameters_OneL(const double rc, const int l_sph){ //TODO Rename this function
	// if neighbours have not been searched perform the research
	if( !_MySystem->getIsNeighbours() || _MySystem->get_current_rc() != rc ){
		_MySystem->searchNeighbours(rc);
	}
	cout << "Computing Steinhart parameter.. ";
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
	unsigned int j_loop;
	int l_loop, l_loop_2;
	#pragma omp parallel for private(xpos,ypos,zpos,j_loop,id,xp,yp,zp,colat,longit,l_loop,l_loop_2)
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
			for(l_loop=-l_sph;l_loop<l_sph+1;l_loop++) Qalpha[i*(l_sph*2+1)+l_loop+l_sph] += spherical_harmonics((unsigned int) l_sph, l_loop, colat, longit);
			// Store the neighbour index into Malpha if it is of the same specy
			if( _MySystem->getAtom(i).type_uint == _MySystem->getAtom(id).type_uint ){
				Malpha[i*(nbNMax+1)] += 1;
				Malpha[i*(nbNMax+1)+Malpha[i*(nbNMax+1)]] = id;
			}
		}
		// compute normalization factors
		for(l_loop_2=-l_sph;l_loop_2<l_sph+1;l_loop_2++) Calpha[i] += (pow(Qalpha[i*(l_sph*2+1)+l_loop_2+l_sph].real(), 2.) + pow(Qalpha[i*(l_sph*2+1)+l_loop_2+l_sph].imag(), 2.));
	}
	cout << 1 << endl;
	cout << " Done !" << endl;
}

// Compute bond orientational parameter, based on the work of Steinhardt, P. J. et al. 1983, modified by Chua et al. 2010 and modified by me for accounting for multisite crystals 
double* ComputeAuxiliary::BondOrientationalParameter(){
	if( _MySystem->getCrystal()->getIsMultisite() )	BondOriParam_Multisite();
	else BondOriParam_NoMultisite();
	return this->BondOriParam;
}

void ComputeAuxiliary::BondOriParam_NoMultisite(){
	double rc = _MySystem->get_rcut();
	int l_sph = _MySystem->get_lsph();
	ComputeSteinhardtParameters_OneL(rc,l_sph);
	cout << "Computing bond orientationnal parameter.. " << endl;

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
		if( BondOriParam[i] > 1. || BondOriParam[i] < -1. ) BondOriParam[i] = 1.;
		BondOriParam[i] = 1.-fabs(BondOriParam[i]);
	}
	cout << "Done !" << endl;
}

void ComputeAuxiliary::BondOriParam_Multisite(){
	if( _MySystem->getCrystal()->getName() != "" ) SteinhardtDatabase_read(_MySystem->getCrystal()->getName());
	cout << "Computing bond orientationnal parameter (multisite version).. " << endl;
	ComputeSteinhardtParameters_OneL(this->rcut_ref,this->l_sph_ref);
	// search site index based on the value of Steinhardt parameter (which is sqrt(pow(Calpha*N,2.)*4pi/(2l+1)))
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	this->Atom_SiteIndex = new unsigned int[nbAt];
	vector<double> diff_St;
	double buffer_d, buffer_d1;
	cout << 1 << endl;
	for(unsigned int i=0;i<nbAt;i++){
		if( this->AtomSiteRefPC[_MySystem->getAtom(i).type_uint-1].size() == 1 ) this->Atom_SiteIndex[i] = 0;
		else{
			diff_St.clear();
			buffer_d = sqrt(Calpha[i]*4.*M_PI/((2*this->l_sph_ref+1)*pow(_MySystem->getNeighbours(i*(nbNMax+1)),2.)));
			for(unsigned int j=0;j<this->AtomSiteRefPC[_MySystem->getAtom(i).type_uint-1].size();j++){
				buffer_d1 = SteinhardtParams_REF_PC[_MySystem->getAtom(i).type_uint-1][j*(this->l_sph_ref+1)+this->l_sph_ref]-buffer_d;
				diff_St.push_back(fabs(buffer_d1));
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
				for(unsigned int l=0;l<(this->l_sph_ref*2+1);l++){
					BondOriParam[i] += ((Qalpha[i*(this->l_sph_ref*2+1)+l].real()*Qalpha[NId*(this->l_sph_ref*2+1)+l].real())+Qalpha[i*(this->l_sph_ref*2+1)+l].imag()*Qalpha[NId*(this->l_sph_ref*2+1)+l].imag())/(pow(Calpha[i],.5)*pow(Calpha[NId],.5)); 
				}
				trueN++;
			}
		}
		if( trueN == 0 ) BondOriParam[i] = 0;
		else BondOriParam[i] /= trueN;
	}
	// Normalized the order parameter
	for(unsigned int i=0;i<nbAt;i++){
		BondOriParam[i] /= this->BondOriParam_REF_PC[_MySystem->getAtom(i).type_uint-1][Atom_SiteIndex[i]];
		//if( BondOriParam[i] > 1. ) BondOriParam[i] = 2.-BondOriParam[i];
		if( BondOriParam[i] > 1. || BondOriParam[i] < -1. ) BondOriParam[i] = 1.;
		BondOriParam[i] = 1.-fabs(BondOriParam[i]);
	}
	cout << "Done !" << endl;
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
	ComputeSteinhardtParameters(rc, l_sph, "Multi", "Multi");
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
				for(unsigned int l=0;l<l_sph+1;l++) St2Print[t][l] += SteinhardtParams_ave[i*(l_sph+1)+l];
				count_ave[t] += 1;
				Already = true;
				break;
			}
		}
		if( !Already ){
			St2Print.push_back(vector<double>());
			count_ave.push_back(1);
			type_printed.push_back(_MySystem->getAtom(i).type_uint);
			for(unsigned int l=0;l<l_sph+1;l++) St2Print[St2Print.size()-1].push_back(SteinhardtParams_ave[i*(l_sph+1)+l]);
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

// Compute the bond orientational parameter knowing the site of each ion (used to store reference bondOri parameter for given crystal)
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
		ComputeSteinhardtParameters(rc, l_sph, "Multi", "Multi");//TODO verify consistency with other
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
	ComputeSteinhardtParameters(rc, l_sph, "Multi", "Multi");//TODO verify consistency with other
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
				for(unsigned int l=0;l<l_sph+1;l++) St2Print[t][l] += SteinhardtParams_ave[At_index[i]*(l_sph+1)+l];
				count_ave[t] += 1;
				Already = true;
				break;
			}
		}
		if( !Already ){
			St2Print.push_back(vector<double>());
			count_ave.push_back(1);
			type_printed.push_back(_MySystem->getAtom(At_index[i]).type_uint);
			for(unsigned int l=0;l<l_sph+1;l++) St2Print[St2Print.size()-1].push_back(SteinhardtParams_ave[At_index[i]*(l_sph+1)+l]);
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



void ComputeAuxiliary::SteinhardtDatabase_read_GMM(string CrystalName){
	if( !this->IsGMMSteinhardtDatabaseRead ){
		string database_tmp = getSteinhardtDatabase(CrystalName);	
		string database_extension=".dat";
		string prefix_GMM = "/GaussianMixtureModel/";	
		string database = database_tmp+prefix_GMM;
		DIR *dir;
		struct dirent *diread;
		const char *env = database.c_str();
		string buffer_s, beg;
		size_t pos;
		if( (dir = opendir(env) ) != nullptr ){
			while( (diread = readdir(dir)) != nullptr ){
				buffer_s = diread->d_name;
				pos = buffer_s.find(database_extension);
				if(pos!=string::npos){
					buffer_s.erase(buffer_s.size()-database_extension.size());
					beg = buffer_s.substr(0,1); 
					if( beg != "." ) this->Struct_GMM_Names.push_back(buffer_s);
				}
			}
			closedir(dir);
		}else{
			perror("opendir");
		}

		unsigned int nbG=Struct_GMM_Names.size();
		unsigned int nbAt_type = _MySystem->getCrystal()->getNbAtomType();
		this->ICovs_GMM = new vector<vector<double>>[nbG*nbAt_type*nbClusterMax_GMM];
		this->Mus_GMM = new vector<double>[nbG*nbAt_type*nbClusterMax_GMM];
		this->Det_GMM = new long double[nbG*nbAt_type*nbClusterMax_GMM];
		this->weight_GMM = new double[nbG*nbAt_type*nbClusterMax_GMM];
		this->nbCluster_GMM = new unsigned int[nbG*nbAt_type];
		vector<string> AtType_GMM;
		unsigned int *AtType_uint_GMM = new unsigned int[nbAt_type];
		unsigned int *lines_attype = new unsigned int[nbAt_type];

		string fullpath2data;
		size_t pos_dim, pos_attype, pos_mu, pos_c, pos_nbC, pos_weight;
		unsigned int line_rc(1000), count, counter_attype, buffer_ui, buffer_ui2, dim;
		double buffer_d1, buffer_d2;
		string line, buffer_s_2;

		for(unsigned int i=0;i<nbG;i++){
			fullpath2data = database.c_str()+Struct_GMM_Names[i]+database_extension;
			ifstream file(fullpath2data, ios::in);
			if(file){
				count = 0;
				counter_attype = 0;
				do{
					getline(file,line);
					// find number of dim 
					pos_dim=line.find("NUMBER_OF_DIMENSION");
					if( pos_dim!=string::npos){
						istringstream text(line);
						text >> buffer_s >> dim;
						for(unsigned int t=0;t<nbAt_type;t++){
							for(unsigned int c=0;c<nbClusterMax_GMM;c++){
								for(unsigned int d1=0;d1<dim;d1++){
									this->Mus_GMM[i*nbAt_type*nbClusterMax_GMM+t*nbClusterMax_GMM+c].push_back(0.);
									this->ICovs_GMM[i*nbAt_type*nbClusterMax_GMM+t*nbClusterMax_GMM+c].push_back(vector<double>());
									for(unsigned int d2=0;d2<dim;d2++){
										this->ICovs_GMM[i*nbAt_type*nbClusterMax_GMM+t*nbClusterMax_GMM+c][this->ICovs_GMM[i*nbAt_type*nbClusterMax_GMM+t*nbClusterMax_GMM+c].size()-1].push_back(0.);
									}
								}
							}
						}
					}
					// find atom type 
					pos_attype=line.find("ATOM_TYPE");
					if( pos_attype!=string::npos){
						istringstream text(line);
						text >> buffer_s >> buffer_s_2;
						AtType_GMM.push_back(buffer_s_2);
						for(unsigned int t=0;t<nbAt_type;t++){
							if( _MySystem->getCrystal()->getAtomType(t+1) == buffer_s_2 ){
								AtType_uint_GMM[counter_attype] = t;
								break;
							}
						}
						lines_attype[counter_attype] = count;
						counter_attype += 1;
					}
					// find atom type 
					pos_nbC=line.find("NUMBER_OF_CLUSTER");
					if( pos_nbC!=string::npos){
						istringstream text(line);
						text >> buffer_s >> nbCluster_GMM[i*nbAt_type+AtType_uint_GMM[counter_attype-1]];
					}


					if( counter_attype > 0 && count == lines_attype[counter_attype-1]+2 ){
						for(unsigned int c=0;c<nbCluster_GMM[i*nbAt_type+AtType_uint_GMM[counter_attype-1]];c++){
							if( c != 0 ) getline(file,line);
							istringstream text(line);
							text >> buffer_s >> this->weight_GMM[i*nbAt_type*nbClusterMax_GMM+AtType_uint_GMM[counter_attype-1]*nbClusterMax_GMM+c];
							getline(file,line);
							istringstream text1(line);
							text1 >> buffer_s >> this->Det_GMM[i*nbAt_type*nbClusterMax_GMM+AtType_uint_GMM[counter_attype-1]*nbClusterMax_GMM+c];
							getline(file,line);
							istringstream text2(line);
							text2 >> buffer_s;
							for(unsigned int d=0;d<dim;d++) text2 >> this->Mus_GMM[i*nbAt_type*nbClusterMax_GMM+AtType_uint_GMM[counter_attype-1]*nbClusterMax_GMM+c][d];
							getline(file,line);
							for(unsigned int d1=0;d1<dim;d1++){
								getline(file,line);
								istringstream text3(line);
								for(unsigned int d2=0;d2<dim;d2++) text3 >> this->ICovs_GMM[i*nbAt_type*nbClusterMax_GMM+AtType_uint_GMM[counter_attype-1]*nbClusterMax_GMM+c][d1][d2];
							}
						}
					}
					count++;
				}while(file);
			}
			if( counter_attype != nbAt_type ) cout << "Warning ! The number of atom type in the Gaussian mixture model of : " << Struct_GMM_Names[i] << " is different than the one of the crystal, it may cause some issues in the following" << endl;
		}
		//for(unsigned int i=0;i<nbG;i++){
		//	cout << Struct_GMM_Names[i] << endl;
		//	for(unsigned int t=0;t<nbAt_type;t++){
		//		cout << "For " << _MySystem->getCrystal()->getAtomType(t+1) << " ions :" << endl;
		//		for(unsigned int c=0;c<nbCluster_GMM[i*nbAt_type+t];c++){
		//			cout << "Custer " << c+1 << endl;
		//			cout << endl;
		//			cout << "weight : " << weight_GMM[i*nbAt_type*nbClusterMax_GMM+t*nbClusterMax_GMM+c] << endl;
		//			cout << "determinant : " << Det_GMM[i*nbAt_type*nbClusterMax_GMM+t*nbClusterMax_GMM+c] << endl;
		//			cout << "esperance : ";
 		//	        	for(unsigned int d=0;d<dim;d++) cout << Mus_GMM[i*nbAt_type*nbClusterMax_GMM+t*nbClusterMax_GMM+c][d] << " ";
		//			cout << endl;
		//			cout << "comat : " << endl;
		//			for(unsigned int d1=0;d1<dim;d1++){
		//				for(unsigned int d2=0;d2<dim;d2++) cout << ICovs_GMM[i*nbAt_type*nbClusterMax_GMM+t*nbClusterMax_GMM+c][d1][d2] << " ";
		//				cout << endl;
		//			}
		//		}
		//	}
		//}
		this->IsGMMSteinhardtDatabaseRead = true;
		delete[] AtType_uint_GMM;
		delete[] lines_attype;
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
	}
}

double* ComputeAuxiliary::StructuralAnalysis_Steinhardt_GMM(const string aux_name){
	SteinhardtDatabase_read_GMM(_MySystem->getCrystal()->getName());

	int l_sph = _MySystem->get_lsph();
	double rc = _MySystem->get_rcut();
	unsigned int aux_ind;
	unsigned int aux_size;
	if( aux_name == "none" ) ComputeSteinhardtParameters(rc, l_sph, "Multi", "Multi");// TODO : add these infos to the GMM database
	else{
		aux_ind = _MySystem->getAuxIdAndSize(aux_name,aux_size);
		if( aux_size != l_sph+1 ) cout << "Warning: the dimension of the provided Steinhardt parameters (" << aux_size << ") is different from the one in /data/FixedParameters/FixedParameters.dat (" << l_sph+1 << " => l_sph+1)" << endl;
	}

	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbAt_type = _MySystem->getCrystal()->getNbAtomType();
	const unsigned int nbStruct = this->Struct_GMM_Names.size();
	const unsigned int dim = this->Mus_GMM[0].size();
	vector<double> Filtered_St;//(20);
	for(unsigned int d=0;d<dim;d++) Filtered_St.push_back(0.);
	vector<long double> MLC(nbStruct); // Maximum likelihood classifier
	long double *N_s = new long double[nbStruct+1]; // prob associated to each structure plus in the law row the sum of these probs
	double *Struct = new double[nbAt*2];
	unsigned int current_type;
	for(unsigned int i=0;i<nbAt;i++){
		if( aux_name == "none" ) for(unsigned int l=0;l<l_sph;l++) Filtered_St[l] = SteinhardtParams_ave[i*(l_sph+1)+l+1];
		else for(unsigned int l=0;l<l_sph;l++) Filtered_St[l] = _MySystem->getAux(aux_ind)[i*(l_sph+1)+l+1];
		current_type = _MySystem->getAtom(i).type_uint;
		N_s[nbStruct] = 0.;
		for(unsigned int s=0;s<nbStruct;s++){
			N_s[s] = 0.;
			for(unsigned int c=0;c<nbCluster_GMM[s*nbAt_type+(current_type-1)];c++)	N_s[s] += MT->Prob_MultidimGaussian(this->ICovs_GMM[s*nbAt_type*nbClusterMax_GMM+(current_type-1)*nbClusterMax_GMM+c], this->Mus_GMM[s*nbAt_type*nbClusterMax_GMM+(current_type-1)*nbClusterMax_GMM+c], this->Det_GMM[s*nbAt_type*nbClusterMax_GMM+(current_type-1)*nbClusterMax_GMM+c], Filtered_St)*weight_GMM[s*nbAt_type*nbClusterMax_GMM+(current_type-1)*nbClusterMax_GMM+c];
			N_s[nbStruct] += N_s[s];
		}
		if( N_s[nbStruct] != 0 ){
			for(unsigned int s=0;s<nbStruct;s++) MLC[s] = N_s[s] / N_s[nbStruct];
			Struct[i*2] = MT->max(MLC);
			Struct[i*2+1] = MLC[Struct[i*2]];
		}else{
			Struct[i*2] = nbStruct;
			Struct[i*2+1] = 0.;
		}
	}

	// write the StructureIndex.txt file
	ofstream writefile("StructureIndex.txt");
	writefile << "Structures index used for the Gaussian Mixture Model structural analysis" << endl;
	for(unsigned int i=0;i<nbStruct;i++) writefile << i << " " << this->Struct_GMM_Names[i] << endl;
	writefile << nbStruct << " Not identified" << endl;
	writefile.close();
	delete[] N_s;
	MLC.clear();
	Filtered_St.clear();
	return Struct;
}

double* ComputeAuxiliary::StructuralAnalysis_Steinhardt(){
	// read database
	SteinhardtDatabase_read(_MySystem->getCrystal()->getName());
	// compute Steinhardt params
	ComputeSteinhardtParameters(this->rcut_ref, this->l_sph_ref, "Multi", "Multi");
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
					D_q[i*(this->Ref_Def_Names.size()+1)] += pow(this->SteinhardtParams_ave[i*(this->l_sph_ref+1)+l]-this->SteinhardtParams_REF_PC_ave_cutoff[j][l],2.);
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
						D_q[i*(this->Ref_Def_Names.size()+1)+j+1] += pow(this->SteinhardtParams_ave[i*(this->l_sph_ref+1)+l]-this->SteinhardtParams_REF_Def_ave_cutoff[j][k*(this->l_sph_ref+1)+l],2.);
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
	}else if( CrystalName == "Alumine" ){
		#ifdef STEINHARDT_ALUMINE_DATABASE
		database = STEINHARDT_ALUMINE_DATABASE;
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

// atomic strain as defined in ovito (i.e. Shimizu, Ogata, Li: Mater. Trans. 48 (2007), 2923)
double* ComputeAuxiliary::Compute_AffineTransfoMatrix(AtomicSystem &AnalyzedSystem, double rc){
	if( !_MySystem->getIsNeighbours() || _MySystem->get_current_rc() != rc ){
		_MySystem->searchNeighbours(rc);
	}

	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	const unsigned int nbAt_ref = AnalyzedSystem.getNbAtom();

	double xpos, ypos, zpos, xpos_r, ypos_r, zpos_r, id;
	double *buffer_mat = new double[9];

	if( nbAt != nbAt_ref ){
		cout << "The number of atom in the reference system and system to compute atomic strain is not the same !" << endl;
		return 0;
	}
	
	// Compute d0 and Vi if not already computed
	if( !this->Reference_AtomicStrain_Computed ){
		this->Vi_inv = new double[nbAt*9];
		this->d0 = new double[nbAt*nbNMax*3];
		
		for(unsigned int i=0;i<nbAt;i++){
			xpos_r = _MySystem->getWrappedPos(i).x;
			ypos_r = _MySystem->getWrappedPos(i).y;
			zpos_r = _MySystem->getWrappedPos(i).z;
			for(unsigned int j=0;j<9;j++) buffer_mat[j] = 0.;
			for(unsigned int n=0;n<_MySystem->getNeighbours(i*(nbNMax+1));n++){
				id = _MySystem->getNeighbours(i*(nbNMax+1)+n+1);
				this->d0[i*nbNMax*3+n*3] = _MySystem->getWrappedPos(id).x+_MySystem->getCLNeighbours(i*nbNMax*3+n*3)*_MySystem->getH1()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+1)*_MySystem->getH2()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+2)*_MySystem->getH3()[0]-xpos_r;
				this->d0[i*nbNMax*3+n*3+1] = _MySystem->getWrappedPos(id).y+_MySystem->getCLNeighbours(i*nbNMax*3+n*3)*_MySystem->getH1()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+1)*_MySystem->getH2()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+2)*_MySystem->getH3()[1]-ypos_r;
				this->d0[i*nbNMax*3+n*3+2] = _MySystem->getWrappedPos(id).z+_MySystem->getCLNeighbours(i*nbNMax*3+n*3)*_MySystem->getH1()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+1)*_MySystem->getH2()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+2)*_MySystem->getH3()[2]-zpos_r;
				for(unsigned int d1=0;d1<3;d1++){
					for(unsigned int d2=0;d2<3;d2++) buffer_mat[d1*3+d2] += this->d0[i*nbNMax*3+n*3+d1]*this->d0[i*nbNMax*3+n*3+d2];
				}
			}
			MT->invert3x3(buffer_mat,buffer_mat);
			for(unsigned int d=0;d<9;d++) Vi_inv[i*9+d] = buffer_mat[d];
		}
		this->Reference_AtomicStrain_Computed = true;
	}

	double *W = new double[9];
	double *buffer_vec = new double[3];
	double *delta = new double[nbAt*3];
	if( !this->IsJi ){
		this->Ji = new double[nbAt*9];
		this->current_d = new double[nbAt*nbNMax*3];
		this->IsJi = true;
	}

	for(unsigned int i=0;i<nbAt;i++){
		for(unsigned int d=0;d<3;d++) delta[i*3+d] = 0.;
		xpos = AnalyzedSystem.getWrappedPos(i).x;
		ypos = AnalyzedSystem.getWrappedPos(i).y;
		zpos = AnalyzedSystem.getWrappedPos(i).z;
		xpos_r = _MySystem->getWrappedPos(i).x;
		ypos_r = _MySystem->getWrappedPos(i).y;
		zpos_r = _MySystem->getWrappedPos(i).z;
		if( (xpos-xpos_r) > AnalyzedSystem.getH1()[0]/2. ) for(unsigned int d=0;d<3;d++) delta[i*3+d] -= AnalyzedSystem.getH1()[d];
		else if( (xpos-xpos_r) < -AnalyzedSystem.getH1()[0]/2. ) for(unsigned int d=0;d<3;d++) delta[i*3+d] += AnalyzedSystem.getH1()[d];
		if( (ypos-ypos_r) > AnalyzedSystem.getH2()[1]/2. ) for(unsigned int d=0;d<3;d++) delta[i*3+d] -= AnalyzedSystem.getH2()[d];
		else if( (ypos-ypos_r) < -AnalyzedSystem.getH2()[1]/2. ) for(unsigned int d=0;d<3;d++) delta[i*3+d] += AnalyzedSystem.getH2()[d];
		if( (zpos-zpos_r) > AnalyzedSystem.getH3()[2]/2. ) for(unsigned int d=0;d<3;d++) delta[i*3+d] -= AnalyzedSystem.getH3()[d];
		else if( (zpos-zpos_r) < -AnalyzedSystem.getH3()[2]/2. ) for(unsigned int d=0;d<3;d++) delta[i*3+d] += AnalyzedSystem.getH3()[d];
	}

	for(unsigned int i=0;i<nbAt;i++){
		xpos = AnalyzedSystem.getWrappedPos(i).x;
		ypos = AnalyzedSystem.getWrappedPos(i).y;
		zpos = AnalyzedSystem.getWrappedPos(i).z;
		for(unsigned int j=0;j<9;j++){
			W[j] = 0.;
			buffer_mat[j] = this->Vi_inv[i*9+j];
		}
		// Compute d and W
		for(unsigned int n=0;n<_MySystem->getNeighbours(i*(nbNMax+1));n++){
			id = _MySystem->getNeighbours(i*(nbNMax+1)+n+1);
			 this->current_d[i*nbNMax*3+n*3] = AnalyzedSystem.getWrappedPos(id).x+_MySystem->getCLNeighbours(i*nbNMax*3+n*3)*AnalyzedSystem.getH1()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+1)*AnalyzedSystem.getH2()[0]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+2)*AnalyzedSystem.getH3()[0]-xpos-delta[i*3];
			 this->current_d[i*nbNMax*3+n*3+1] = AnalyzedSystem.getWrappedPos(id).y+_MySystem->getCLNeighbours(i*nbNMax*3+n*3)*AnalyzedSystem.getH1()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+1)*AnalyzedSystem.getH2()[1]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+2)*AnalyzedSystem.getH3()[1]-ypos-delta[i*3+1];
			 this->current_d[i*nbNMax*3+n*3+2] = AnalyzedSystem.getWrappedPos(id).z+_MySystem->getCLNeighbours(i*nbNMax*3+n*3)*AnalyzedSystem.getH1()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+1)*AnalyzedSystem.getH2()[2]+_MySystem->getCLNeighbours(i*nbNMax*3+n*3+2)*AnalyzedSystem.getH3()[2]-zpos-delta[i*3+2];
			if( this->current_d[i*nbNMax*3+n*3] > AnalyzedSystem.getH1()[0]/2. ) for(unsigned int d=0;d<3;d++) this->current_d[i*nbNMax*3+n*3+d] -= AnalyzedSystem.getH1()[d];
			else if( this->current_d[i*nbNMax*3+n*3] < -AnalyzedSystem.getH1()[0]/2. ) for(unsigned int d=0;d<3;d++) this->current_d[i*nbNMax*3+n*3+d] += AnalyzedSystem.getH1()[d];
			if( this->current_d[i*nbNMax*3+n*3+1] > AnalyzedSystem.getH2()[1]/2. ) for(unsigned int d=0;d<3;d++) this->current_d[i*nbNMax*3+n*3+d] -= AnalyzedSystem.getH2()[d];
			else if( this->current_d[i*nbNMax*3+n*3+1] < -AnalyzedSystem.getH2()[1]/2. ) for(unsigned int d=0;d<3;d++) this->current_d[i*nbNMax*3+n*3+d] += AnalyzedSystem.getH2()[d];
			if( this->current_d[i*nbNMax*3+n*3+2] > AnalyzedSystem.getH3()[2]/2. ) for(unsigned int d=0;d<3;d++) this->current_d[i*nbNMax*3+n*3+d] -= AnalyzedSystem.getH3()[d];
			else if( this->current_d[i*nbNMax*3+n*3+2] < -AnalyzedSystem.getH3()[2]/2. ) for(unsigned int d=0;d<3;d++) this->current_d[i*nbNMax*3+n*3+d] += AnalyzedSystem.getH3()[d];
			for(unsigned int d1=0;d1<3;d1++){
				for(unsigned int d2=0;d2<3;d2++) W[d1*3+d2] += this->d0[i*nbNMax*3+n*3+d1]*(this->current_d[i*nbNMax*3+n*3+d2]);
			}
		}
		// Compute deformation gradient tensor
		MT->MatDotMat(buffer_mat,W,W);
		for(unsigned int d=0;d<9;d++) this->Ji[i*9+d] = W[d];
	}
	
	delete[] W;
	delete[] buffer_vec;
	delete[] buffer_mat;
	delete[] delta;
	
	return this->Ji;
}
		
double* ComputeAuxiliary::Compute_AtomicStrain(AtomicSystem &AnalyzedSystem, double rc){

	Compute_AffineTransfoMatrix(AnalyzedSystem,rc);
	
	const unsigned int nbAt = _MySystem->getNbAtom();
	
	double *buffer_mat = new double[9];
	double *buffer_mat_1 = new double[9];
	
	if( !this->IsAtomicStrain ){
		this->AtomicStrain = new double[nbAt*8]; // sorted as eps_xx,eps_yy,eps_zz,_eps_xy,eps_xz,eps_yz,shear_inv,hydro_inv
		this->IsAtomicStrain = true;
	}

	for(unsigned int i=0;i<nbAt;i++){
		for(unsigned int d=0;d<9;d++) buffer_mat_1[d] = this->Ji[i*9+d];
		// Compute Lagrangian strain tensor (.5*(J.Jt-Id))	
		MT->MatDotTMat(buffer_mat_1,buffer_mat_1,buffer_mat);

		for(unsigned int d=0;d<3;d++){
			buffer_mat[d*3+d] -= 1.;
			buffer_mat[d*3+d] *= .5; //originally all strain component must be divided by two but regarding analytical tests we don't do that here
		}
		// sort the tensor 
		AtomicStrain[i*8] = buffer_mat[0];
		AtomicStrain[i*8+1] = buffer_mat[4];
		AtomicStrain[i*8+2] = buffer_mat[8];
		AtomicStrain[i*8+3] = buffer_mat[1];
		AtomicStrain[i*8+4] = buffer_mat[2];
		AtomicStrain[i*8+5] = buffer_mat[5];
		// compute shear invariant
		AtomicStrain[i*8+6] = 0.;
		for(unsigned int d1=0;d1<3;d1++) AtomicStrain[i*8+6] += pow(AtomicStrain[i*8+3+d1],2.);
		AtomicStrain[i*8+6] += (pow(buffer_mat[4]-buffer_mat[8],2.)+pow(buffer_mat[0]-buffer_mat[8],2.)+pow(buffer_mat[0]-buffer_mat[4],2.))/8.;
		AtomicStrain[i*8+6] = sqrt(AtomicStrain[i*8+6]);
		// compute hydrostatic invariant

		AtomicStrain[i*8+7] = 0.;
		for(unsigned int d1=0;d1<3;d1++) AtomicStrain[i*8+7] += buffer_mat[d1*3+d1];
		AtomicStrain[i*8+7] /= 3.;
	}
	
	delete[] buffer_mat;
	delete[] buffer_mat_1;

	return AtomicStrain;
}

// D2Min as defined in Delbecq et al. 2023
double* ComputeAuxiliary::Compute_D2Min(AtomicSystem &AnalyzedSystem, double rc){ 

	Compute_AffineTransfoMatrix(AnalyzedSystem,rc);
	
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	if( !this->IsD2Min ){
		this->D2Min = new double[nbAt];
	}

	double *buffer_vec = new double[3];
	double *buffer_mat = new double[9];

	for(unsigned int i=0;i<nbAt;i++){
		this->D2Min[i] = 0.;
		for(unsigned int d=0;d<9;d++) buffer_mat[d] = this->Ji[i*9+d];
		for(unsigned int n=0;n<_MySystem->getNeighbours(i*(nbNMax+1));n++){
			// Compute Ji*delta0
			for(unsigned int d1=0;d1<3;d1++) buffer_vec[d1] = this->d0[i*nbNMax*3+n*3+d1];
			MT->MatDotRawVec(buffer_mat,buffer_vec,buffer_vec);
			// Compute sum( (delta - Ji*delta0)^2
			for(unsigned int d=0;d<3;d++) this->D2Min[i] += pow(this->current_d[i*nbNMax*3+n*3+d]-buffer_vec[d],2.);
		}
		// Normalize D2Min with respect to the number of neighbours
		// this->D2Min[i] /= _MySystem->getNeighbours(i*(nbNMax+1));
	}

	delete[] buffer_vec;
	delete[] buffer_mat;

	return this->D2Min;
}

ComputeAuxiliary::~ComputeAuxiliary(){
	delete MT;
	if( IsBondOriParam ){
		delete[] BondOriParam;
		delete[] Atom_SiteIndex;
	}
	if( this->Reference_AtomicStrain_Computed ){
		delete[] Vi_inv;
		delete[] d0;
		delete[] current_d;
	}
	if( this->IsJi ){
		delete[] Ji;
	}
	if( this->IsAtomicStrain ){
		delete[] AtomicStrain;
	}
	if( this->IsD2Min ){
		delete[] D2Min;
	}
	if( IsStrainTensor ){
		delete[] StrainTensor;
	}
	if( IsStrainInvII ){
		delete[] Strain_invII;
	}
	if( IsSteinhardt_Mono || IsSteinhardt_Multi || IsSteinhardt_Filtered ){
		delete[] Qlm;
		delete[] SteinhardtParams;
	}
	if( IsSteinhardt_Mono || IsSteinhardt_Multi ){
		delete[] Malpha;
	}
	if( IsSteinhardt_AveMono || IsSteinhardt_AveMulti || IsSteinhardt_AveFiltered || IsSteinhardt_FilteredAveFiltered || IsSteinhardt_AveFilteredMono || IsSteinhardt_AveFilteredMulti ){
		delete[] SteinhardtParams_ave;
	}
	if( IsSteinhardt_Filtered || IsSteinhardt_FilteredAveFiltered || IsSteinhardt_AveFilteredMono || IsSteinhardt_AveFilteredMulti ){
		delete[] nbNeigh_FiltNeigh;
	}
	if( IsGMMSteinhardtDatabaseRead ){
		delete[] ICovs_GMM;
		delete[] Mus_GMM;
		delete[] Det_GMM;
		delete[] nbCluster_GMM;
	        delete[] weight_GMM;	
	}
	//TODO structural analysis
}

