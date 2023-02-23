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
#include <dirent.h>

using namespace std;

// Compute Steinhart parameters which correspond to a vector of dimension l_sph*2+1 for each atom representing the complex norm of Qalpha components
double* ComputeAuxiliary::ComputeSteinhardtParameters(const double rc, const int l_sph){
	// if neighbours have not been searched perform the research
	if( !_MySystem->getIsNeighbours() ) _MySystem->searchNeighbours(rc);
	cout << "Computing Steinhart parameters.. ";
	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
	this->BondOriParam = new double[nbAt];
	this->Atom_SiteIndex = new unsigned int[nbAt];
	this->Malpha = new unsigned int[nbAt*(nbNMax+1)]; // array containing the index of neighbours of the same species (or same site in case of multisite crystal) with the first line corresponding to the number of neighbours, i.e. Malpha[i*(nbNMax+1)] = nb of neighbour of atom i, Malpha[i*(nbNMax+1)+j+1] = id of the jth neighbour of atom i
	this->Qalpha = new complex<double>[nbAt*(l_sph*2+1)]; // complex array containing the spherical harmonic for the different modes
	this->SteinhardtParams = new double[nbAt*(l_sph*2+1)];
	this->Calpha = new double[nbAt]; // normalization factor 
	for(unsigned int i=0;i<nbAt*(l_sph*2+1);i++) Qalpha[i] = (0.,0.); // initialize it to zero
	double zeronum = 1e-8;
	const int bar_length = 30;
	double prog=0;
	double xpos,ypos,zpos,xp,yp,zp, colat, longit;
	unsigned int id;
	// loop on all atoms and neighbours to compute Qalpha and store neighbour of the same species
	// Here is the most time consuming loop of the function, use parallel computation
	unsigned int j_loop, st;
	int l_loop;
	#pragma omp parallel for private(xpos,ypos,zpos,j_loop,id,xp,yp,zp,colat,longit,l_loop,st)
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
		for(st=0;st<(l_sph*2)+1;st++) SteinhardtParams[i*(l_sph*2+1)+st] = sqrt(pow(Qalpha[i*(l_sph*2+1)+st].real(),2.)+pow(Qalpha[i*(l_sph*2+1)+st].imag(),2.)); 
		// compute normalization factors
		for(int l=-l_sph;l<l_sph+1;l++)	Calpha[i] += (pow(Qalpha[i*(l_sph*2+1)+l+l_sph].real(), 2.) + pow(Qalpha[i*(l_sph*2+1)+l+l_sph].imag(), 2.));
	}
	cout << "Done !" << endl;
}

double* ComputeAuxiliary::BondOriParam_SteinhardtBased(){
	// read database
	SteinhardtDatabase_read(_MySystem->getCrystal()->getName());
	// compute Steinhardt params
	ComputeSteinhardtParameters(this->rcut_ref, this->l_sph_ref);
	// search if multisite 
	bool multiSite = false;
	if( _MySystem->getIsCrystalDefined() ){
		for(unsigned int i=0;i<_MySystem->getCrystal()->getNbAtomType();i++){
			if(_MySystem->getCrystal()->getNbAtomSite(i) > 1){
				multiSite = true;
				break;
			}
		}
	}
	// compute Calpha and normalization factor for references
	double *Calpha_ref = new double[AtomTypeUINTRefPC.size()];
	for(unsigned int i=0;i<AtomTypeUINTRefPC.size();i++){
		Calpha_ref[i] = 0;
		for(unsigned int l=0;l<this->l_sph_ref*2+1;l++)	Calpha_ref[i] += (pow(this->SteinhardtParams_REF_PC[i][(this->l_sph_ref*2+1)+l], 2.) + pow(this->SteinhardtParams_REF_PC[i][2*(this->l_sph_ref*2+1)+l], 2.));
	}
	const unsigned int nbAt = _MySystem->getNbAtom();
	this->BondOriParam_Steinhardt = new double[nbAt];
	//test
	vector<double> *BondOri = new vector<double>[AtomTypeUINTRefPC.size()];
	double *NormFac = new double[AtomTypeUINTRefPC.size()];
	vector<double> BondOriFromRef;
	//end test
	if( multiSite ){ // search from the ref which is the closest to each atom
		unsigned int *closest_ref = new unsigned int[nbAt];
		vector<double> dist;
		vector<unsigned int> ind;
		double norm_dist;
		for(unsigned int i=0;i<nbAt;i++){
			//dist.clear();
			//ind.clear();
			//for(unsigned int j=0;j<AtomTypeUINTRefPC.size();j++){
			//	if( _MySystem->getAtom(i).type_uint == AtomTypeUINTRefPC[j] ){
			//		dist.push_back(0.);
			//		ind.push_back(j);
			//		for(unsigned int l=0;l<this->l_sph_ref*2+1;l++) dist[dist.size()-1] += pow(this->SteinhardtParams[i*(this->l_sph_ref*2+1)+l]-this->SteinhardtParams_REF_PC[j][l],2.);
			//	}
			//}
			//closest_ref[i] = ind[MT->min(dist)];
			//// compute the order param
			//this->BondOriParam_Steinhardt[i] = 0.;
			//for(unsigned int l=0;l<this->l_sph_ref*2+1;l++)	BondOriParam_Steinhardt[i] += ( (fabs(this->Qalpha[i*(this->l_sph_ref*2+1)+l].real())*this->SteinhardtParams_REF_PC[closest_ref[i]][(this->l_sph_ref*2+1)+l]) + fabs((this->Qalpha[i*(this->l_sph_ref*2+1)+l].imag())*this->SteinhardtParams_REF_PC[closest_ref[i]][2*(this->l_sph_ref*2+1)+l]) );
			//this->BondOriParam_Steinhardt[i] /= ( sqrt(this->Calpha[i])*sqrt(Calpha_ref[closest_ref[i]]) );
			//BondOri[closest_ref[i]].push_back(this->BondOriParam_Steinhardt[i]);

			//test search the ref which maximize the order param
			BondOriFromRef.clear();
			for(unsigned int j=0;j<AtomTypeUINTRefPC.size();j++){
				BondOriFromRef.push_back(0.);
				if( _MySystem->getAtom(i).type_uint == AtomTypeUINTRefPC[j] ){
					for(unsigned int l=0;l<this->l_sph_ref*2+1;l++)	BondOriFromRef[BondOriFromRef.size()-1] += ( (fabs(this->Qalpha[i*(this->l_sph_ref*2+1)+l].real())*this->SteinhardtParams_REF_PC[j][(this->l_sph_ref*2+1)+l]) + fabs((this->Qalpha[i*(this->l_sph_ref*2+1)+l].imag())*this->SteinhardtParams_REF_PC[j][2*(this->l_sph_ref*2+1)+l]) );
					BondOriFromRef[BondOriFromRef.size()-1] /= ( sqrt(this->Calpha[i])*sqrt(Calpha_ref[j]) );
				}
			}
			this->BondOriParam_Steinhardt[i] = MT->max_vec(BondOriFromRef);
			closest_ref[i] = MT->max(BondOriFromRef);
			BondOri[closest_ref[i]].push_back(this->BondOriParam_Steinhardt[i]);
		}
		// search the maximum order param for each ref
		for(unsigned int i=0;i<AtomTypeUINTRefPC.size();i++){
			NormFac[i] = MT->max_vec(BondOri[i]);
			cout << NormFac[i] << " " << endl;
		}
		// renormalize according to this
		for(unsigned int i=0;i<nbAt;i++) this->BondOriParam_Steinhardt[i] /= NormFac[closest_ref[i]];
	}//end multisite
	//ofstream writefile("closest_ref.dat");
	//for(unsigned int i=0;i<nbAt;i++) writefile << closest_ref[i] << endl;
	//writefile.close();
	
	return this->BondOriParam_Steinhardt;
}

// Compute bond orientational parameter, based on the work of Steinhardt, P. J. et al. 1983, modified by Chua et al. 2010 and modified by me for accounting for multisite crystals 
double* ComputeAuxiliary::BondOrientationalParameter(){
	double rc = _MySystem->get_rcut();
	int l_sph = _MySystem->get_lsph();
	ComputeSteinhardtParameters(rc,l_sph);
	cout << "Computing bond orientationnal parameter" << endl;

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


	const unsigned int nbAt = _MySystem->getNbAtom();
	const unsigned int nbNMax = _MySystem->getNbMaxN();
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
		vector<vector<unsigned int>> SiteIndex; // array containing the site index 
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

void ComputeAuxiliary::SaveSteinhardtParamToDatabase_PerfectCrystal(string CrystalName){
	int l_sph = _MySystem->get_lsph();
	double rc = _MySystem->get_rcut();
	ComputeSteinhardtParameters(rc, l_sph);
	// read the database environment
	string filename = "PerfectCrystal.dat";
	string fullpathname = SteinhardtDatabase_write(CrystalName)+filename;
	const unsigned int nbAt = _MySystem->getNbAtom();
	bool multiSite = false;
	if( _MySystem->getIsCrystalDefined() ){
		for(unsigned int i=0;i<_MySystem->getCrystal()->getNbAtomType();i++){
			if(_MySystem->getCrystal()->getNbAtomSite(i) > 1){

				multiSite = true;
				break;
			}
		}
	}

	ofstream writefile(fullpathname);
	writefile << "Steinhardt parameters for " << CrystalName << " perfect crystal computed using" << endl;
	writefile << "L_SPH " << l_sph << endl;
	writefile << "RCUT " << rc << endl;
	
	if( multiSite ){
		// compute the bond orientational param in order to have the different sites
		BondOrientationalParameter();
		bool Already = false;
		bool AlreadyType = false;
		unsigned int typeok;
		vector<vector<double>> St2Print;
		vector<unsigned int> type_printed;
		vector<vector<unsigned int>> site_printed;
		vector<vector<unsigned int>> count_ave;
		// average the steinhardt parameters TODO maybe do the same thing with real and imag part of Qalpha
		for(unsigned int i=0;i<nbAt;i++){
			Already = false;
			AlreadyType = false;
			for(unsigned int t=0;t<type_printed.size();t++){
				if( _MySystem->getAtom(i).type_uint == type_printed[t] ){
					AlreadyType = true;
					for(unsigned int s=0;s<site_printed[t].size();s++){
						if( Atom_SiteIndex[i] == site_printed[t][s] ){
							for(unsigned int l=0;l<l_sph*2+1;l++){
								St2Print[t*3][s*(l_sph*2+1)+l] += SteinhardtParams[i*(l_sph*2+1)+l];
								St2Print[t*3+1][s*(l_sph*2+1)+l] += fabs(this->Qalpha[i*(l_sph*2+1)+l].real());
								St2Print[t*3+2][s*(l_sph*2+1)+l] += fabs(this->Qalpha[i*(l_sph*2+1)+l].imag());
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
					for(unsigned int st=0;st<3;st++) St2Print.push_back(vector<double>());
					count_ave.push_back(vector<unsigned int>());
					site_printed.push_back(vector<unsigned int>());
					count_ave[count_ave.size()-1].push_back(1);
					site_printed[site_printed.size()-1].push_back(Atom_SiteIndex[i]);
					type_printed.push_back(_MySystem->getAtom(i).type_uint);
					for(unsigned int l=0;l<l_sph*2+1;l++){
						St2Print[(St2Print.size()/3-1)*3].push_back(SteinhardtParams[i*(l_sph*2+1)+l]);
						St2Print[(St2Print.size()/3-1)*3+1].push_back(fabs(this->Qalpha[i*(l_sph*2+1)+l].real()));
						St2Print[(St2Print.size()/3-1)*3+2].push_back(fabs(this->Qalpha[i*(l_sph*2+1)+l].imag()));
					}
				}else{ // not a new type, just a new site
					count_ave[typeok].push_back(1);
					site_printed[typeok].push_back(Atom_SiteIndex[i]);
					for(unsigned int l=0;l<l_sph*2+1;l++){
						St2Print[typeok*3].push_back(SteinhardtParams[i*(l_sph*2+1)+l]);
						St2Print[typeok*3+1].push_back(fabs(this->Qalpha[i*(l_sph*2+1)+l].real()));
						St2Print[typeok*3+2].push_back(fabs(this->Qalpha[i*(l_sph*2+1)+l].imag()));
					}
				}
			}
		}
		unsigned int nbref = 0;
		for(unsigned int t=0;t<type_printed.size();t++){
			for(unsigned int s=0;s<site_printed[t].size();s++){
				nbref += 1;
				for(unsigned int l=0;l<l_sph*2+1;l++){
					St2Print[t*3][s*(l_sph*2+1)+l] /= count_ave[t][s];
					St2Print[t*3+1][s*(l_sph*2+1)+l] /= count_ave[t][s];
					St2Print[t*3+2][s*(l_sph*2+1)+l] /= count_ave[t][s];
				}
			}
		}
		writefile << "NB_REF " << nbref << endl;
		for(unsigned int t=0;t<type_printed.size();t++){
			for(unsigned int s=0;s<site_printed[t].size();s++){
				writefile << _MySystem->getCrystal()->getAtomType(type_printed[t]-1) << " " << type_printed[t] << " S" << s+1 << " ";
				for(unsigned int l=0;l<l_sph*2+1;l++) writefile << St2Print[t*3][s*(l_sph*2+1)+l] << " ";
				for(unsigned int l=0;l<l_sph*2+1;l++) writefile << St2Print[t*3+1][s*(l_sph*2+1)+l] << " " << St2Print[t*3+2][s*(l_sph*2+1)+l] << " ";
				writefile << endl;
			}
		}
		writefile.close();
	}else{ // no multisite case
		bool Already = false;
		vector<vector<double>> St2Print;
		vector<unsigned int> type_printed;
		vector<unsigned int> count_ave;
		// average the steinhardt parameters TODO maybe do the same thing with real and imag part of Qalpha
		for(unsigned int i=0;i<nbAt;i++){
			Already = false;
			for(unsigned int t=0;t<type_printed.size();t++){
				if( _MySystem->getAtom(i).type_uint == type_printed[t] ){
					for(unsigned int l=0;l<l_sph*2+1;l++){
						St2Print[t*3][l] += SteinhardtParams[i*(l_sph*2+1)+l];
						St2Print[t*3+1][l] += fabs(this->Qalpha[i*(l_sph*2+1)+l].real());
						St2Print[t*3+2][l] += fabs(this->Qalpha[i*(l_sph*2+1)+l].imag());
					}
					count_ave[t] += 1;
					Already = true;
					break;
				}
			}
			if( !Already ){
				for(unsigned int st=0;st<3;st++) St2Print.push_back(vector<double>());
				count_ave.push_back(1);
				type_printed.push_back(_MySystem->getAtom(i).type_uint);
				for(unsigned int l=0;l<l_sph*2+1;l++){
					St2Print[(St2Print.size()/3-1)*3].push_back(SteinhardtParams[i*(l_sph*2+1)+l]);
					St2Print[(St2Print.size()/3-1)*3+1].push_back(fabs(this->Qalpha[i*(l_sph*2+1)+l].real()));
					St2Print[(St2Print.size()/3-1)*3+2].push_back(fabs(this->Qalpha[i*(l_sph*2+1)+l].imag()));
				}
			}
		}
		for(unsigned int t=0;t<type_printed.size();t++){
			for(unsigned int l=0;l<l_sph*2+1;l++){
				St2Print[t*3][l] /= count_ave[t];
				St2Print[t*3+1][l] /= count_ave[t];
				St2Print[t*3+2][l] /= count_ave[t];
			}
		}
		writefile << "NB_REF " << type_printed.size() << endl;
		for(unsigned int t=0;t<type_printed.size();t++){
			writefile << _MySystem->getCrystal()->getAtomType(type_printed[t]-1) << " " << type_printed[t] << " ";
			for(unsigned int l=0;l<l_sph*2+1;l++) writefile << St2Print[t*3][l] << " ";
			for(unsigned int l=0;l<l_sph*2+1;l++) writefile << St2Print[t*3+1][l] << " " << St2Print[t*3+2][l] << " ";
			writefile << endl;
		}
		writefile.close();
	}

}

void ComputeAuxiliary::SaveSteinhardtParamToDatabase_Defect(string CrystalName, string filename){
//TODO
}

void ComputeAuxiliary::SteinhardtDatabase_read(string CrystalName){
	string database;	
	string database_extension=".dat";	
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
		DIR *dir;
		struct dirent *diread;
		const char *env = database.c_str();
		string buffer_s;
		size_t pos;
		if( (dir = opendir(env) ) != nullptr ){
			while( (diread = readdir(dir)) != nullptr ){
				buffer_s = diread->d_name;
				pos = buffer_s.find(database_extension);
				if(pos!=string::npos){
					if( buffer_s.erase(buffer_s.size()-database_extension.size()) != "PerfectCrystal" ) this->Ref_Def_Names.push_back(buffer_s);
				}
			}
			closedir(dir);
		}else{
			perror("opendir");
		}
	}
	// Search if crystal is multisite
	bool multiSite = false;
	if( _MySystem->getIsCrystalDefined() ){
		for(unsigned int i=0;i<_MySystem->getCrystal()->getNbAtomType();i++){
			if(_MySystem->getCrystal()->getNbAtomSite(i) > 1){

				multiSite = true;
				break;
			}
		}
	}
	// Read perfect crystal Steinhardt params
	string PC_str="PerfectCrystal";
	string fullpath2data = database.c_str()+backslash+PC_str+database_extension;
	ifstream file(fullpath2data, ios::in);
	size_t pos_nbref, pos_lsph, pos_rc;
	unsigned int line_rc(1000), count, nbref;
	string buffer_s, line;
	count = 0;
	if(file){
		do{
			getline(file,line);
			// find number of ref
			pos_nbref=line.find("NB_REF ");
			if( pos_nbref!=string::npos){
				istringstream text(line);
				text >> buffer_s >> nbref;
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
			if( count > line_rc+1 && count < line_rc+nbref+2 ){
				SteinhardtParams_REF_PC.push_back(new double[(this->l_sph_ref*2+1)*3]);
				AtomTypeUINTRefPC.push_back(0);
				istringstream text(line);
				text >> buffer_s >> AtomTypeUINTRefPC[AtomTypeUINTRefPC.size()-1] >> buffer_s;
				// get norm of steinhardt params
				for(unsigned int l=0;l<this->l_sph_ref*2+1;l++) text >> SteinhardtParams_REF_PC[SteinhardtParams_REF_PC.size()-1][l];
				// get real and imag part
				for(unsigned int l=0;l<this->l_sph_ref*2+1;l++) text >> SteinhardtParams_REF_PC[SteinhardtParams_REF_PC.size()-1][(this->l_sph_ref*2+1)+l] >> SteinhardtParams_REF_PC[SteinhardtParams_REF_PC.size()-1][2*(this->l_sph_ref*2+1)+l];
			}
			count += 1;
		}while(file);
	}

	this->nbRefDef = this->Ref_Def_Names.size();
	// TODO defect reading

}

string ComputeAuxiliary::SteinhardtDatabase_write(string CrystalName){
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

ComputeAuxiliary::~ComputeAuxiliary(){
	delete MT;
	if( IsBondOriParam ){
		delete[] BondOriParam;
		delete[] Atom_SiteIndex;
	}
	if( IsStrainTensor ){
		delete[] StrainTensor;
	}
	if( IsStrainInvII ){
		delete[] Strain_invII;
	}
}

