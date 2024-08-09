#include "SteinhardtDescriptors.h"
#include "AtomHicConfig.h"
#include <cmath>
#include <complex>
#include <omp.h>
#include <iomanip>

using namespace std; 

SteinhardtDescriptors::SteinhardtDescriptors(AtomicSystem *_MySystem):Descriptors(_MySystem){
	this->name = "Steinhardt";
	readFixedParams();
	dim = l_sph;
	// if neighbours have not been searched perform the research
	if( !_MySystem->getIsNeighbours() || _MySystem->get_current_rc() != rc ){
		_MySystem->searchNeighbours(rc);
	}
	nbNMax = _MySystem->getNbMaxN();
	_Descriptors = new double[nbDatTot*l_sph]; // TODO case with one L or filtered
	lsph2 = (l_sph+1)*(l_sph*2+1.);
	lsph1 = l_sph*2+1;
	Malpha = new unsigned int [nbDatTot*(nbNMax+1)];
	Qlm = new complex<double>[nbDatTot*lsph2];
	for(unsigned int i=0;i<nbDatTot*(l_sph*2+1);i++){
		for(unsigned int j=0;j<l_sph+1;j++) Qlm[i*j+j] = (0.,0.);
	}
	ComputeDescriptors();
}

void SteinhardtDescriptors::readFixedParams(){
	string fp;
	#ifdef FIXEDPARAMETERS
	fp = FIXEDPARAMETERS;
	#endif
	string backslash="/";
	string filename=fp+backslash+FixedParam_Filename;
	ifstream file(filename, ios::in);
	size_t pos_rcut, pos_lsph, pos_StStyle, pos_AveStyle;
	string buffer_s, line;
	if(file){
		while(file){
			getline(file,line);
			pos_rcut=line.find("STEINHARDT_DESCRIPTORS_RCUT ");
			if(pos_rcut!=string::npos){
				istringstream text(line);
				text >> buffer_s >> rc;
			}
			pos_lsph=line.find("STEINHARDT_DESCRIPTORS_L_SPH ");
			if(pos_lsph!=string::npos){
				istringstream text(line);
				text >> buffer_s >> l_sph;
			}
			pos_StStyle=line.find("STEINHARDT_DESCRIPTORS_STYLE ");
			if(pos_StStyle!=string::npos){
				istringstream text(line);
				text >> buffer_s >> SteinhardtStyle;
			}
			pos_AveStyle=line.find("STEINHARDT_DESCRIPTORS_AVERAGE_STYLE ");
			if(pos_AveStyle!=string::npos){
				istringstream text(line);
				text >> buffer_s >> AverageStyle;
			}
		}
	}else{
		cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
		exit(EXIT_FAILURE);
	}
	file.close();
	cout << "From /data/FixedParameters/FixedParameters.dat we will compute Steinhardt descriptors using cutoff radius of " << rc << " and degree " << l_sph << " with style \"" << SteinhardtStyle << "\" and average style \"" << AverageStyle << "\"" << endl;
}

void SteinhardtDescriptors::ComputeDescriptors(){
	// Steinhardt style
	cout << "Computing Steinhardt parameters" << endl;
	if( SteinhardtStyle == "Mono" ) ComputeSteinhardtParameters_Mono();
	else if( SteinhardtStyle == "Multi" ) ComputeSteinhardtParameters_Multi();
	else{
		cerr << "The Steinhardt parameter style to print has not been recognized" << endl;
	       	cerr << "Possible styles : " << endl;
		cerr << "Mono : Classic Steinhardt parameters using ions ofthe same species" << endl;
		cerr << "Multi : Classic Steinhardt parameters using all ion types" << endl;
		exit(EXIT_FAILURE);
	}
	cout << "Done" << endl;
	// Averaging
	if( AverageStyle == "Mono" ) AverageSteinhardtParameters_Mono();
	else if( AverageStyle == "Multi" ) AverageSteinhardtParameters_Multi();
	else if( AverageStyle == "none" ) cout << "Steinhardt descriptors successfully computed" << endl;
	else{
		cerr << "The style of average for the Steinhardt parameters to print has not been recognized" << endl;
	       	cerr << "Possible styles : " << endl;
		cerr << "none : no averaging" << endl;
		cerr << "Mono : averaging over ions ofthe same species" << endl;
		cerr << "Multi : averaging over all ion types" << endl;
		exit(EXIT_FAILURE);
	}
}

// Compute Steinhart parameters which correspond to a vector of dimension l_sph*2+1 for each atom representing the complex norm of Qalpha components
void SteinhardtDescriptors::ComputeSteinhardtParameters_Mono(){
	prog=0;
	count_t = 0;
	// loop on all atoms and neighbours to compute Qalpha and store neighbour of the same species
	#pragma omp parallel for
	for(unsigned int i=0;i<nbDatTot;i++){
		if( omp_get_thread_num() == 0 ){
			prog = double(count_t*omp_get_num_threads())/double(nbDatTot);
			cout << "\r[" << string(bar_length*prog,'X') << string(bar_length*(1-prog),'-') << "] " << setprecision(3) << 100*prog << "%";
			count_t++;
		}
		double xpos = _MySystem->getWrappedPos(i).x;
		double ypos = _MySystem->getWrappedPos(i).y;
		double zpos = _MySystem->getWrappedPos(i).z;
		Malpha[i*(nbNMax+1)] = 0; 
		for(unsigned int j_loop=0;j_loop<_MySystem->getNeighbours(i*(nbNMax+1));j_loop++){
			unsigned int id = _MySystem->getNeighbours(i*(nbNMax+1)+j_loop+1);
			if( _MySystem->getAtom(i).type_uint == _MySystem->getAtom(id).type_uint ){
				// get distance vector
				unsigned int ind1 = i*nbNMax*3+j_loop*3;
				unsigned int ind2 = ind1+1;
				unsigned int ind3 = ind2+1;
				double xp = _MySystem->getWrappedPos(id).x+_MySystem->getCLNeighbours(ind1)*_MySystem->getH1()[0]+_MySystem->getCLNeighbours(ind2)*_MySystem->getH2()[0]+_MySystem->getCLNeighbours(ind3)*_MySystem->getH3()[0]-xpos;
				double yp = _MySystem->getWrappedPos(id).y+_MySystem->getCLNeighbours(ind1)*_MySystem->getH1()[1]+_MySystem->getCLNeighbours(ind2)*_MySystem->getH2()[1]+_MySystem->getCLNeighbours(ind3)*_MySystem->getH3()[1]-ypos;
				double zp = _MySystem->getWrappedPos(id).z+_MySystem->getCLNeighbours(ind1)*_MySystem->getH1()[2]+_MySystem->getCLNeighbours(ind2)*_MySystem->getH2()[2]+_MySystem->getCLNeighbours(ind3)*_MySystem->getH3()[2]-zpos;
				// compute colatitude and longitudinal angles
				double colat = acos(zp/sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)));
				double longit;
				if( xp > 0 ) longit = atan(yp/xp);
	                	else if( ( xp < 0 ) && ( yp >= 0 ) ) longit = atan(yp/xp) + M_PI;
	                	else if( ( xp < 0 ) and ( yp < 0 ) ) longit = atan(yp/xp) - M_PI;
	                	else if( ( fabs(xp) < zeronum ) and ( yp > 0 ) ) longit = M_PI/2.;
	                	else if( ( fabs(xp) < zeronum ) and ( yp < 0 ) ) longit = -M_PI/2.;
	                	else if( ( fabs(xp) < zeronum ) and ( fabs(yp) < zeronum ) ) longit = 0.;
				// compute spherical harmonics
				for(unsigned int l_loop_st=0;l_loop_st<l_sph+1;l_loop_st++){
					for(int m_loop_st=-l_loop_st;m_loop_st<(int) l_loop_st+1;m_loop_st++){
						Qlm[i*lsph2 + l_loop_st*lsph1 + (unsigned int) (m_loop_st + (int) l_sph)] += MT->spherical_harmonics(l_loop_st, m_loop_st, colat, longit);
					}
				}
				Malpha[i*(nbNMax+1)] += 1;
				Malpha[i*(nbNMax+1)+Malpha[i*(nbNMax+1)]] = id;
			}
		}
	}

	if( AverageStyle == "none"){
		#pragma omp parallel for
		for(unsigned int i=0;i<nbDatTot;i++){
			for(unsigned int l_loop_st2=1;l_loop_st2<l_sph+1;l_loop_st2++){
				_Descriptors[i*l_sph+l_loop_st2-1] = 0.; 
				for(int m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
					_Descriptors[i*l_sph+l_loop_st2-1] += norm(Qlm[i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
				}
				_Descriptors[i*l_sph+l_loop_st2-1] *= 4.*M_PI/(2.*l_loop_st2+1.); 
				_Descriptors[i*l_sph+l_loop_st2-1] /= pow(_MySystem->getNeighbours(i*(nbNMax+1)),2.);
				_Descriptors[i*l_sph+l_loop_st2-1] = sqrt(_Descriptors[i*l_sph+l_loop_st2-1]);
			}
		}
	}
	cout << endl;
}

// Same as above but with the average version where the Steinhardt params are averaged using all neighbours (in the case above only the ions of same type inside the cutoff radius are used)
void SteinhardtDescriptors::ComputeSteinhardtParameters_Multi(){
	count_t = 0;
	prog=0;
	// loop on all atoms and neighbours to compute Qalpha and store neighbour of the same species
	#pragma omp parallel for
	for(unsigned int i=0;i<nbDatTot;i++){
		if( omp_get_thread_num() == 0 ){
			prog = double(count_t*omp_get_num_threads())/double(nbDatTot);
			cout << "\r[" << string(bar_length*prog,'X') << string(bar_length*(1-prog),'-') << "] " << setprecision(3) << 100*prog << "%";
			count_t++;
		}
		double xpos = _MySystem->getWrappedPos(i).x;
		double ypos = _MySystem->getWrappedPos(i).y;
		double zpos = _MySystem->getWrappedPos(i).z;
		Malpha[i*(nbNMax+1)] = 0; 
		for(unsigned int j_loop=0;j_loop<_MySystem->getNeighbours(i*(nbNMax+1));j_loop++){
			// get distance vector
			unsigned int id = _MySystem->getNeighbours(i*(nbNMax+1)+j_loop+1);
			unsigned int ind1 = i*nbNMax*3+j_loop*3;
			unsigned int ind2 = ind1+1;
			unsigned int ind3 = ind2+1;
			double xp = _MySystem->getWrappedPos(id).x+_MySystem->getCLNeighbours(ind1)*_MySystem->getH1()[0]+_MySystem->getCLNeighbours(ind2)*_MySystem->getH2()[0]+_MySystem->getCLNeighbours(ind3)*_MySystem->getH3()[0]-xpos;
			double yp = _MySystem->getWrappedPos(id).y+_MySystem->getCLNeighbours(ind1)*_MySystem->getH1()[1]+_MySystem->getCLNeighbours(ind2)*_MySystem->getH2()[1]+_MySystem->getCLNeighbours(ind3)*_MySystem->getH3()[1]-ypos;
			double zp = _MySystem->getWrappedPos(id).z+_MySystem->getCLNeighbours(ind1)*_MySystem->getH1()[2]+_MySystem->getCLNeighbours(ind2)*_MySystem->getH2()[2]+_MySystem->getCLNeighbours(ind3)*_MySystem->getH3()[2]-zpos;
			// compute colatitude and longitudinal angles
			double colat = acos(zp/sqrt(pow(xp,2.)+pow(yp,2.)+pow(zp,2.)));
			double longit;
			if( xp > 0 ) longit = atan(yp/xp);
	                else if( ( xp < 0 ) && ( yp >= 0 ) ) longit = atan(yp/xp) + M_PI;
	                else if( ( xp < 0 ) and ( yp < 0 ) ) longit = atan(yp/xp) - M_PI;
	                else if( ( fabs(xp) < zeronum ) and ( yp > 0 ) ) longit = M_PI/2.;
	                else if( ( fabs(xp) < zeronum ) and ( yp < 0 ) ) longit = -M_PI/2.;
	                else if( ( fabs(xp) < zeronum ) and ( fabs(yp) < zeronum ) ) longit = 0.;
			// compute spherical harmonics
			for(unsigned int l_loop_st=0;l_loop_st<l_sph+1;l_loop_st++){
				for(int m_loop_st=-l_loop_st;m_loop_st<(int) l_loop_st+1;m_loop_st++){
					Qlm[i*lsph2 + l_loop_st*lsph1 + (unsigned int) (m_loop_st + (int) l_sph)] += MT->spherical_harmonics(l_loop_st, m_loop_st, colat, longit);
				}
			}
			// Store the neighbour index into Malpha if it is of the same specy
			if( _MySystem->getAtom(i).type_uint == _MySystem->getAtom(id).type_uint ){
				Malpha[i*(nbNMax+1)] += 1;
				Malpha[i*(nbNMax+1)+Malpha[i*(nbNMax+1)]] = id;
			}
		}
	}
	if( AverageStyle == "none"){
		#pragma omp parallel for
		for(unsigned int i=0;i<nbDatTot;i++){
			for(unsigned int l_loop_st2=1;l_loop_st2<l_sph+1;l_loop_st2++){
				_Descriptors[i*l_sph+l_loop_st2-1] = 0.; 
				for(int m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
					_Descriptors[i*l_sph+l_loop_st2-1] += norm(Qlm[i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
				}
				_Descriptors[i*l_sph+l_loop_st2-1] *= 4.*M_PI/(2.*l_loop_st2+1.); 
				_Descriptors[i*l_sph+l_loop_st2-1] /= pow(_MySystem->getNeighbours(i*(nbNMax+1)),2.);
				_Descriptors[i*l_sph+l_loop_st2-1] = sqrt(_Descriptors[i*l_sph+l_loop_st2-1]);
			}
		}
	}
	cout << endl;
}

void SteinhardtDescriptors::AverageSteinhardtParameters_Mono(){
	cout << "Averaging Steinhardt parameters" << endl;
	vector<complex<double>> buffer_complex(lsph2);
	#pragma omp parallel for firstprivate(buffer_complex)
	for(unsigned int i=0;i<nbDatTot;i++){
		for(unsigned int l_loop_st2=1;l_loop_st2<l_sph+1;l_loop_st2++){
			_Descriptors[i*l_sph+l_loop_st2-1] = 0.; 
			for(int m_loop_st0=-l_loop_st2;m_loop_st0<(int) l_loop_st2+1;m_loop_st0++){
				buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] = Qlm[i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] / (double) _MySystem->getNeighbours(i*(nbNMax+1));
			}
			for(unsigned int neigh=0;neigh<Malpha[i*(nbNMax+1)];neigh++){
			       for(int m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++){
				       buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[Malpha[i*(nbNMax+1)+neigh+1]*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) _MySystem->getNeighbours(Malpha[i*(nbNMax+1)+neigh+1]*(nbNMax+1)));
			       }
			}
			for(int m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
				_Descriptors[i*l_sph+l_loop_st2-1] += norm(buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
			}
			_Descriptors[i*l_sph+l_loop_st2-1] *= 4.*M_PI/(2.*l_loop_st2+1.); 
			_Descriptors[i*l_sph+l_loop_st2-1] /= pow(1+Malpha[i*(nbNMax+1)],2.);
			_Descriptors[i*l_sph+l_loop_st2-1] = sqrt(_Descriptors[i*l_sph+l_loop_st2-1]);
		}
	}
	cout << "Steinhardt descriptors successfully computed" << endl;
}

void SteinhardtDescriptors::AverageSteinhardtParameters_Multi(){
	cout << "Averaging Steinhardt parameters" << endl;
	vector<complex<double>> buffer_complex(lsph2);
	#pragma omp parallel for firstprivate(buffer_complex)
	for(unsigned int i=0;i<nbDatTot;i++){
		for(unsigned int l_loop_st2=1;l_loop_st2<l_sph+1;l_loop_st2++){
			_Descriptors[i*l_sph+l_loop_st2-1] = 0.; 
			for(int m_loop_st0=-l_loop_st2;m_loop_st0<(int) l_loop_st2+1;m_loop_st0++){
				buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] = Qlm[i*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st0 + (int) l_sph)] / (double) _MySystem->getNeighbours(i*(nbNMax+1));
			}
			for(unsigned int neigh=0;neigh<_MySystem->getNeighbours(i*(nbNMax+1));neigh++){
			       for(int m_loop_st1=-l_loop_st2;m_loop_st1<(int) l_loop_st2+1;m_loop_st1++){
				       buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] += Qlm[_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)*lsph2 + l_loop_st2*lsph1 + (unsigned int) (m_loop_st1 + (int) l_sph)] / ((double) _MySystem->getNeighbours(_MySystem->getNeighbours(i*(nbNMax+1)+neigh+1)*(nbNMax+1)));
			       }
			}
			for(int m_loop_st2=-l_loop_st2;m_loop_st2<(int) l_loop_st2+1;m_loop_st2++){
				_Descriptors[i*l_sph+l_loop_st2-1] += norm(buffer_complex[l_loop_st2*lsph1 + (unsigned int) (m_loop_st2 + (int) l_sph)]);
			}
			_Descriptors[i*l_sph+l_loop_st2-1] *= 4.*M_PI/((2.*l_loop_st2)+1.); 
			_Descriptors[i*l_sph+l_loop_st2-1] /= pow(1+_MySystem->getNeighbours(i*(nbNMax+1)),2.);
			_Descriptors[i*l_sph+l_loop_st2-1] = sqrt(_Descriptors[i*l_sph+l_loop_st2-1]);
		}
	}
	cout << "Steinhardt descriptors successfully computed" << endl;
}

void SteinhardtDescriptors::setProperties(){
	Properties.push_back("NUMBER_OF_DIMENSION "+to_string(l_sph));
	Properties.push_back("CUTOFF_RADIUS "+to_string(rc));
	Properties.push_back("STEINHARDT_STYLE "+SteinhardtStyle);
	Properties.push_back("AVE_STYLE "+AverageStyle);
}
	
SteinhardtDescriptors::~SteinhardtDescriptors(){
}
