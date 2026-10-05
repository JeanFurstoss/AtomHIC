//**********************************************************************************
//*   Bicrystal.cpp                                                                *
//**********************************************************************************
//* This file contains the implementation of the Bicrystal class                   *
//**********************************************************************************
//* (C) Jan 2025 - Jean Furstoss                                                   *
//*     Université de Poitiers, Institut PPRIME                                    *
//*     UPR CNRS 3346, 86360 Chasseuneuil-du-Poitou, France                        *
//*     jean.furstoss@univ-poitiers.fr                                             *
//* Last modification: J. Furstoss - 28 Janv 2025                                  *
//**********************************************************************************
//* This program is free software: you can redistribute it and/or modify           *
//* it under the terms of the GNU General Public License as published by           *
//* the Free Software Foundation, either version 3 of the License, or              *
//* (at your option) any later version.                                            *
//*                                                                                *
//* This program is distributed in the hope that it will be useful,                *
//* but WITHOUT ANY WARRANTY; without even the implied warranty of                 *
//* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                  *
//* GNU General Public License for more details.                                   *
//*                                                                                *
//* You should have received a copy of the GNU General Public License              *
//* along with this program.  If not, see <http://www.gnu.org/licenses/>.          *
//**********************************************************************************

#include "AtomHicConfig.h"
#include <iostream>
#include <vector>
#include "Bicrystal.h"
#include <omp.h>

using namespace std;

Bicrystal::Bicrystal(const string& crystalName, int h_a, int k_a, int l_a, double theta):h_a(h_a), k_a(k_a), l_a(l_a), theta(theta){
	this->MT = new MathTools;
	BT = new BicrystalloTools;
	read_params();
	setCrystal(crystalName);
	InitCSL_DSC();
}

// Constructor for bicrystal with facetted GB with given misorientation and GB plane and facet type (for the moment this is implemented for 2D facets and for a GB composed 2 different facet types, i.e. with one facet direction parallel to x direction)
Bicrystal::Bicrystal(const string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p, vector<int> FacetsType, double L_facet, vector<string> Properties, int h_px, int k_px, int l_px):h_a(h_a), k_a(k_a), l_a(l_a), theta(theta), h_p(h_p), k_p(k_p), l_p(l_p){ // TODO will not work with cubic primitive cell
	double zero_num = 1.e-9;
	this->MT = new MathTools;
	BT = new BicrystalloTools;
	read_params();
	setCrystal(crystalName);
	InitCSL_DSC();
	string i_a_str(""), i_p_str("");
	if( _MyCrystal->getCrystallo() == "Hexagonal" ){
		i_a_str = "_"+to_string(-(this->h_a)-(this->k_a));
		i_p_str = "_"+to_string(-(this->h_p)-(this->k_p));
	}
	cout << "Constructing a GB misoriented by " << theta*180./M_PI << "° around [" << h_a << "_" << k_a << i_a_str << "_" << l_a << "] axis and one GB plane (" << h_p << "_" << k_p << i_p_str << "_" << l_p << ")" << endl;
	setOrientedCrystals(crystalName, false, Properties, h_px, k_px, l_px);
	_MyCrystal->ComputeOrthogonalPlanesAndDirections();
	_MyCrystal2->ComputeOrthogonalPlanesAndDirections();
	double *Dir1_G1 = new double[3];
	double *Dir2_G1 = new double[3];
	double *Dir1_G2 = new double[3];
	double *Dir2_G2 = new double[3];
	double facnorm1(0.), facnorm2(0.);
	double sp_t = 0.;
	for(unsigned int i=0;i<3;i++){
		Dir1_G1[i] = FacetsType[0]*_MyCrystal->getA1_star()[i]+FacetsType[1]*_MyCrystal->getA2_star()[i]+FacetsType[2]*_MyCrystal->getA3_star()[i];
		Dir2_G1[i] = FacetsType[3]*_MyCrystal->getA1_star()[i]+FacetsType[4]*_MyCrystal->getA2_star()[i]+FacetsType[5]*_MyCrystal->getA3_star()[i];
		facnorm1 += Dir1_G1[i]*Dir1_G1[i];
		facnorm2 += Dir2_G1[i]*Dir2_G1[i];
	}
	facnorm1 = sqrt(facnorm1);
	facnorm2 = sqrt(facnorm2);
	for(unsigned int d=0;d<2;d++) sp_t += Dir1_G1[d]*Dir2_G1[d];
	if( sp_t < 0. ){
		for(unsigned int i=0;i<3;i++) FacetsType[3+i] *= -1;
		facnorm2 *= -1.;
	}
	for(unsigned int i=0;i<3;i++){
		Dir1_G1[i] /= facnorm1;
		Dir2_G1[i] /= facnorm2;
	}
	double *FacetJctDir = new double[3];
	MT->mixedProd(Dir1_G1,Dir2_G1,FacetJctDir);
	double abs_xf = fabs(FacetJctDir[0]);
	double abs_yf = fabs(FacetJctDir[1]);
	double abs_zf = fabs(FacetJctDir[2]);
	bool good_facet = true;
 	if( abs_xf < zero_num && abs_yf < zero_num && abs_zf < zero_num ){
		cerr << "The two provided facet planes are colinear, the faceted GB cannot be constructed, aborting calculation" << endl;
		good_facet = false;
	}
	if( fabs(Dir1_G1[0]) < zero_num && fabs(Dir1_G1[1]) < zero_num ){
		cerr << "The first facet is horyzontal, cannot construct a faceted GB with this facet, aborting calculation" << endl;
		good_facet = false;
	}
	if( fabs(Dir2_G1[0]) < zero_num && fabs(Dir2_G1[1]) < zero_num ){
		cerr << "The second facet is horyzontal, cannot construct a faceted GB with this facet, aborting calculation" << endl;
		good_facet = false;
	}
	if( !(( Dir1_G1[2] < 0. && Dir2_G1[2] >= 0 ) || ( Dir1_G1[2] <= 0. && Dir2_G1[2] > 0 ) || ( Dir1_G1[2] > 0. && Dir2_G1[2] <= 0 ) || ( Dir1_G1[2] >= 0. && Dir2_G1[2] < 0 )) ){
		cerr << "The two provided facets planes have the same z-slope sign, the periodic boundary conditions cannot then be respected, aborting calculation" << endl; 
		good_facet = false;
	}
	if( abs_zf > zero_num ){
		cerr << "The line junction between the two provided facets is not parallel to z-plane, the periodic boundary conditions cannot then be respected, aborting calculation" << endl;
		good_facet = false;
	}
	if( !good_facet ){
		cerr << "Grain 1 (i.e. crystal reference frame for the facets) is oriented following:" << endl;
		Dis->DisplayOrthogonalCell(_MyCrystal);
		unsigned int max_hkl_search_facet = 0;
		cerr << "Enter the maximum hkl value for searching admissible facet for this system: ";
		cin >> max_hkl_search_facet;
		cout << endl;
		ShowPossibleFacets(max_hkl_search_facet);
		exit(EXIT_FAILURE);
	}
	Dis->DisplayFacetedGB(_MyCrystal,_MyCrystal2,Dir1_G1,Dir2_G1);
	// search the number of linear combination for the two system to have the same x y length (same as for classic GB construction TODO put in a function)
	this->xl1 = this->_MyCrystal->getOrientedSystem()->getH1()[0];
	this->xl2 = this->_MyCrystal2->getOrientedSystem()->getH1()[0];
	this->yl1 = this->_MyCrystal->getOrientedSystem()->getH2()[1];
	this->yl2 = this->_MyCrystal2->getOrientedSystem()->getH2()[1];
	this->zl1 = this->_MyCrystal->getOrientedSystem()->getH3()[2];
	this->zl2 = this->_MyCrystal2->getOrientedSystem()->getH3()[2];
	//  make a four imbricated loop to search in the same time what combination of G1 and G2 respect misfit and the PBC for Facets
	double final_xbox, final_ybox;
	bool found_dup = false;
	for(unsigned int i1=1;i1<MaxDup;i1++){
		for(unsigned int i2=1;i2<MaxDup;i2++){
			final_xbox = (xl1*i1 + xl2*i2)/2.;
    			double temp_MisfitX = fabs(1. - (final_xbox / (xl1*i1)));
			if( temp_MisfitX < MaxMisfit ){//&& final_xbox > InitMinBoxAside ){
				for(unsigned int j1=1;j1<MaxDup;j1++){
					for(unsigned int j2=1;j2<MaxDup;j2++){
						final_ybox = (yl1*j1 + yl2*j2)/2.;
    						double temp_MisfitY = fabs(1. - (final_ybox / (yl1*j1)));
						if( temp_MisfitY < MaxMisfit ){//&& final_ybox > InitMinBoxAside ){
    							Mx1 = final_xbox / (xl1*i1);
    							My1 = final_ybox / (yl1*j1);
							_MyCrystal->getA1()[0] *= Mx1;
							_MyCrystal->getA1()[1] *= My1;
							_MyCrystal->getA2()[0] *= Mx1;
							_MyCrystal->getA2()[1] *= My1;
							_MyCrystal->getA3()[0] *= Mx1;
							_MyCrystal->getA3()[1] *= My1;
							_MyCrystal->computeReciproqual();
							facnorm1 = 0.;
							facnorm2 = 0.;
							for(unsigned int i=0;i<3;i++){
								Dir1_G1[i] = FacetsType[0]*_MyCrystal->getA1_star()[i]+FacetsType[1]*_MyCrystal->getA2_star()[i]+FacetsType[2]*_MyCrystal->getA3_star()[i];
								Dir2_G1[i] = FacetsType[3]*_MyCrystal->getA1_star()[i]+FacetsType[4]*_MyCrystal->getA2_star()[i]+FacetsType[5]*_MyCrystal->getA3_star()[i];
								facnorm1 += Dir1_G1[i]*Dir1_G1[i];
								facnorm2 += Dir2_G1[i]*Dir2_G1[i];
							}
							facnorm1 = sqrt(facnorm1);
							facnorm2 = sqrt(facnorm2);
							for(unsigned int i=0;i<3;i++){
								Dir1_G1[i] /= facnorm1;
								Dir2_G1[i] /= facnorm2;
							}

							MT->mixedProd(Dir1_G1,Dir2_G1,FacetJctDir);
							double facnorm3 = 0.;
							for(unsigned int i=0;i<3;i++) facnorm3 += FacetJctDir[i]*FacetJctDir[i];
							facnorm3 = sqrt(facnorm3);
							for(unsigned int i=0;i<3;i++) FacetJctDir[i] /= facnorm3;
							abs_xf = fabs(FacetJctDir[0]);
							abs_yf = fabs(FacetJctDir[1]);
							double tol_bc = 1e-1; // in A TODO FixedParameters ?
							if( fabs( (final_xbox*abs_yf) - (final_ybox*abs_xf) ) < tol_bc || abs_xf < zero_num || abs_yf < zero_num ){
								found_dup = true;
								dupX1 = i1;
								dupY1 = j1;
								dupX2 = i2;
								dupY2 = j2;
    								Mx2 = final_xbox / (xl2*i2);
    								My2 = final_ybox / (yl2*j2);
								break;
							}else{
								_MyCrystal->getA1()[0] /= Mx1;
								_MyCrystal->getA1()[1] /= My1;
								_MyCrystal->getA2()[0] /= Mx1;
								_MyCrystal->getA2()[1] /= My1;
								_MyCrystal->getA3()[0] /= Mx1;
								_MyCrystal->getA3()[1] /= My1;
							}
						}
					} // end loop j2
					if( found_dup ) break;
				} // end loop j1
				if( found_dup ) break;
			} // end if i1,i2
		} // end loop i2
		if( found_dup ) break;
	} // end loop i1
	if( !found_dup ){
		cerr << "We do not find any combination of the two grains to have a misfit lower than provided by the fixed parameters and the facet respecting periodic boundary conditions" << endl;
		cerr << "Increase MAX_MISFIT or MAX_DUP parameters in a FixedParameters.ath file" << endl;
		exit(EXIT_FAILURE);
	}	
							
	cout << "Pasting the two grains over each others leads to:" << endl;
    	vector<vector<string>> arr_elements = {
		{"", "Grain 1", "Grain 2"}, 
		{"Direction", "X", "Y", "X", "Y"},
		{"Duplicates", to_string(dupX1), to_string(dupY1), to_string(dupX2), to_string(dupY2)},
		{"Misfits", to_string((Mx1-1.)*100.)+" %", to_string((My1-1.)*100.)+" %", to_string((Mx2-1.)*100.)+" %", to_string((My2-1.)*100.)+" %"}
	};
	vector<vector<unsigned int>> arr_fusion = {
		{1,  2, 2},
		{1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1}
    	};
	Dis->DisplayArray(arr_elements, arr_fusion);
	cout << "Misfit values could be reduced by decreasing MAX_MISFIT in a FixedParameters.ath file" << endl;
	cout << endl;

	// compute the corrected facet length and the number of facet to be put in the box
	// Facet length should be a mutliple of the reticular distance
	unsigned int dupSys = 1;
	double *l1 = new double[3]; // normal of facet junction comprised in facet 1
	double *l2 = new double[3]; // normal of facet junction comprised in facet 2
	double buffer;
	int buffer_int;
	MT->mixedProd(Dir1_G1,FacetJctDir,l1);
	// if the first facet is vertical exchange facets
	if( fabs(Dir1_G1[2]) < zero_num  ){
		MT->mixedProd(FacetJctDir,Dir2_G1,l1);
		MT->mixedProd(Dir1_G1,FacetJctDir,l2);
		for(unsigned int i=0;i<3;i++){
			buffer = Dir1_G1[i];
			Dir1_G1[i] = Dir2_G1[i];
			Dir2_G1[i] = buffer;
			buffer_int = FacetsType[i];
			FacetsType[i] = FacetsType[3+i];
			FacetsType[3+i] = buffer_int;
		}
	}else
	MT->mixedProd(FacetJctDir,Dir2_G1,l2);
	double norm_l1 = sqrt( pow(l1[0],2) + pow(l1[1],2) + pow(l1[2],2) );
	double norm_l2 = sqrt( pow(l2[0],2) + pow(l2[1],2) + pow(l2[2],2) );
	for(unsigned int i=0;i<3;i++){
		l1[i] /= norm_l1;
		l1[i] /= norm_l2;
	}
	double norm_xy_l1 = sqrt( pow(l1[0],2) + pow(l1[1],2) );
	double norm_xy_l2 = sqrt( pow(l2[0],2) + pow(l2[1],2) );
	double fac_L_facet_xy;
	fac_L_facet_xy = (norm_xy_l1/norm_l1) - ( norm_xy_l2*l1[2] / (norm_l2*l2[2]) );
	double alpha = acos( (norm_xy_l1*norm_xy_l2) + (l1[2]*l2[2]) );
	double dfac1 = _MyCrystal->ComputeD_hkl(FacetsType[0],FacetsType[1],FacetsType[2]) / sin(alpha);
	double dfac2 = _MyCrystal->ComputeD_hkl(FacetsType[3],FacetsType[4],FacetsType[5]) / sin(alpha);
	unsigned int Nfacets;

	// Facet length should be a mutliple of the reticular distance
	double n2On1 = -l1[2]*dfac2/(l2[2]*dfac1);
	unsigned int nLcor = round(L_facet / dfac2 );
	unsigned int nd_max = round(L_facet * maxVarL / dfac2);
	double new_L_facet = nLcor * dfac2;
	vector<unsigned int> PossibleN;
	// search possible facets for the two length being a multiple of reticular distance
	for(unsigned int i=0;i<nd_max+1;i++){
	        double n2 = double (nLcor + i) * n2On1;
	        if( fabs(round(n2)-n2) < 1e-6 ) PossibleN.push_back(nLcor+i);
	        n2 = double (nLcor - i) * n2On1;
	        if( fabs(round(n2)-n2) < 1e-6 ) PossibleN.push_back(nLcor-i);
	}
	if( PossibleN.size() == 0 ){
		cerr << "Cannot find a combination of facets length respecting reticular distances (increase MAX_VAR_FACET_LENGTH parameter in a FixedParameters.ath file)" << endl;
		exit(EXIT_FAILURE);
	}
	double denom = ((norm_xy_l1/norm_l1) - ( norm_xy_l2*l1[2] / (norm_l2*l2[2]) ))*norm_xy_l1;
	double nume = (fabs(l1[0])*final_xbox) + (fabs(l1[1])*final_ybox);
	double tol = 1e-5;
	bool found = false;
	bool ShouldBeEven = true;
	if( fabs(FacetJctDir[0]) < zero_num || fabs(FacetJctDir[1]) < zero_num ) ShouldBeEven = false;
	for(unsigned int k=1;k<dupMaxFac;k++){
		unsigned int current_N = round(k*nume/(nLcor*dfac2*denom));
		if( ShouldBeEven ){
			if( current_N < 2 ) continue;
			else if( current_N % 2 != 0 ) current_N++;
		}else if( current_N < 1 ) continue;
		for(unsigned int i_N=0;i_N<2;i_N++){
			unsigned int N = 0;
			if( ShouldBeEven ) N = current_N + i_N*2;
			else N = current_N + i_N;
			for(unsigned int i_n=0;i_n<PossibleN.size();i_n++){
				unsigned int n = PossibleN[i_n];
				if( fabs( (N*n*dfac2*denom) - (k*nume) ) < tol ){
					Nfacets = N;
					dupSys = k;
					nLcor = n;
					found = true;
					break;
				}
			}
			if( found ) break;
			if( ShouldBeEven ) N = current_N - i_N*2;
			else N = current_N - i_N;
			for(unsigned int i_n=0;i_n<PossibleN.size();i_n++){
				unsigned int n = PossibleN[i_n];
				if( fabs( (N*n*dfac2*denom) - (k*nume) ) < tol ){
					Nfacets = N;
					dupSys = k;
					nLcor = n;
					found = true;
					break;
				}
			}
			if( found ) break;
		}
		if( found ) break;
	}
	if( !found ){
		cerr << "Cannot find a facet length allowing to respect periodic boundary conditions (increase MAX_VAR_FACET_LENGTH or MAX_DUP_FACET parameters in a FixedParameters.ath file)" << endl;
		exit(EXIT_FAILURE);
	}
	new_L_facet = nLcor * dfac2;
	if( fabs(FacetJctDir[1]) > zero_num ){
		final_xbox *= dupSys;
		dupX1 *= dupSys;
		dupX2 *= dupSys;
	}
	if( fabs(FacetJctDir[0]) > zero_num ){
		final_ybox *= dupSys;
		dupY1 *= dupSys;
		dupY2 *= dupSys;
	}

	cout << "The provided facet length (" << L_facet << " A) has been replaced by " << new_L_facet << " A in order to respect PBC in GB plane and intereticular distances" << endl;
	if( Nfacets < 2 ) cout << "If the difference between the provided facet length and the new one is too high, you can try to change the MIN_BOX_ASIDE parameter in a FixedParameters.ath file, to provide a system dimension closer to the one of the provided facet length" << endl;
	unsigned int final_dupx = ceil(MinBoxSizes[0]/final_xbox);
	unsigned int final_dupy = ceil(MinBoxSizes[1]/final_ybox);
	cout << "The system will be composed of facets such that:" << endl;
	L_facet = new_L_facet;
	double L_facet2 = -L_facet*l1[2]/l2[2];
	double FacetLength;
        if( fabs(FacetJctDir[0]) < zero_num ) FacetLength = final_ybox;
	else FacetLength = final_xbox/cos(atan(FacetJctDir[1]/FacetJctDir[0]));
	if( fabs(round(L_facet/dfac2)-(L_facet/dfac2)) > 1e-5 ) cout << "WARNING !! The length of facet 1 is not a multiple of the reticular distance of facet 2, some issues could occur !" << endl;
	if( fabs(round(L_facet2/dfac1)-(L_facet2/dfac1)) > 1e-5 ) cout << "WARNING !! The length of facet 2 is not a multiple of the reticular distance of facet 1, some issues could occur !" << endl;
	cout << "\t facet type (" << FacetsType[0] << "_" << FacetsType[1] << "_" << FacetsType[2] << ") has length : " << ((double) final_dupx*final_dupy)*L_facet << " A, for total area of " << ((double) final_dupx*final_dupy)*double (Nfacets)*FacetLength*L_facet << " A2" << endl;
	cout << "\t facet type (" << FacetsType[3] << "_" << FacetsType[4] << "_" << FacetsType[5] << ") has length : " << ((double) final_dupx*final_dupy)*L_facet2 << " A, for total area of " << ((double) final_dupx*final_dupy)*double (Nfacets)*FacetLength*L_facet2 << " A2" << endl;
	double L1xy = L_facet*norm_xy_l1/norm_l1;
	double L2xy = L_facet2*norm_xy_l2/norm_l2;

	// Compute the facet height see if the grains should be duplicated in z direction
	double DeltaZ = L_facet*sin(acos(norm_xy_l1/norm_l1))/2.; // divied by two as the 2 grains share equally the facet
	unsigned int zdup1 = ceil((this->MinBoxSizes[2]+2.*DeltaZ)/_MyCrystal->getOrientedSystem()->getH3()[2])+1; 
	unsigned int zdup2 = ceil((this->MinBoxSizes[2]+2.*DeltaZ)/_MyCrystal2->getOrientedSystem()->getH3()[2])+1; 

	// initialize the variables and pointers
	unsigned int nbAtom1 = this->_MyCrystal->getOrientedSystem()->getNbAtom();
	unsigned int nbAtom2 = this->_MyCrystal2->getOrientedSystem()->getNbAtom();
	unsigned int nbAtom_temp = nbAtom1*dupX1*dupY1*zdup1 + nbAtom2*dupX2*dupY2*zdup2;
	Atom *AtomList_temp = new Atom[nbAtom_temp]; // init this array with the maximum number of atom
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
	this->H3[2] = (this->_MyCrystal->getOrientedSystem()->getH3()[2])*zdup1 + ((this->_MyCrystal2->getOrientedSystem()->getH3()[2]+GBspace)*zdup2) + GBspace;
	this->H1_G1[0] = Mx1*xl1*dupX1;
	this->H1_G1[1] = 0.;
	this->H1_G1[2] = 0.;
	this->H2_G1[0] = 0.;
	this->H2_G1[1] = My1*yl1*dupY1;
	this->H2_G1[2] = 0.;
	this->H3_G1[0] = 0.;
	this->H3_G1[1] = 0.;
	double zboxG1 = (this->_MyCrystal->getOrientedSystem()->getH3()[2])*zdup1;
	this->H3_G1[2] = zboxG1;
	this->H1_G2[0] = Mx2*xl2*dupX2;
	this->H1_G2[1] = 0.;
	this->H1_G2[2] = 0.;
	this->H2_G2[0] = 0.;
	this->H2_G2[1] = My2*yl2*dupY2;
	this->H2_G2[2] = 0.;
	this->H3_G2[0] = 0.;
	this->H3_G2[1] = 0.;
	double zboxG2 = (this->_MyCrystal2->getOrientedSystem()->getH3()[2])*zdup2;
	this->H3_G2[2] = zboxG2;
	this->nbAtomType = _MyCrystal->getNbAtomType();
	this->AtomType = _MyCrystal->getAtomType();
	this->AtomMass = _MyCrystal->getAtomMass();
	this->AtomCharge = _MyCrystal->getAtomCharge();
	this->IsCharge = _MyCrystal->getIsCharge();
	this->IsTilted = false;
	computeInverseCellVec();

	// Compute the equation of plane of the different facets (full facet planes + vertical planes)
	unsigned int NbPlanes = 2*(Nfacets+1);// times 2 because there is 2 planes per facets and just need to store 1 value per plane as the other ones are the normal
	double *FullPlaneEq = new double[NbPlanes]; 
	double *NormalPlaneEq = new double[NbPlanes*3]; // just store the normal of facet plane of each plane (redondant but allows to save time in the atom loop
	double *VertPlaneEq = new double[NbPlanes];
	double VertPlaneSlope;
	bool xpar = false;
	if( fabs(FacetJctDir[0]) > zero_num ) VertPlaneSlope = FacetJctDir[1]/FacetJctDir[0];
	else{
		VertPlaneSlope = 0.;
		xpar = true;
	}
	double invcosphi = 1./cos(atan(VertPlaneSlope));

	if( Dir1_G1[2] < 0 ) for(unsigned int i=0;i<3;i++) Dir1_G1[i] *= -1.;
	else for(unsigned int i=0;i<3;i++) Dir2_G1[i] *= -1.;

	if( VertPlaneSlope < 0. ){
		VertPlaneEq[0] = -1.; // values of y for x = 0
		FullPlaneEq[0] = 0.;
	}else{
		VertPlaneEq[0] = - VertPlaneSlope*final_xbox - 1.; 
		FullPlaneEq[0] = - Dir1_G1[1]*VertPlaneSlope*final_xbox;
	}
	for(unsigned int i=0;i<3;i++) NormalPlaneEq[i] = Dir1_G1[i];

	for(unsigned int i=0;i<Nfacets;i++){
		if( VertPlaneSlope < 0. ){
			VertPlaneEq[i*2+1] = ((L1xy/2.) + i*(L2xy+L1xy))*invcosphi; // values of y for x = 0
			VertPlaneEq[i*2+2] = ((L1xy/2.) + L2xy + i*(L2xy+L1xy))*invcosphi;
			FullPlaneEq[i*2+1] = Dir2_G1[1]*(VertPlaneEq[i*2+1]+(L2xy*invcosphi/2.)); // values of y for x = 0 and z = 0
			FullPlaneEq[i*2+2] = Dir1_G1[1]*(VertPlaneEq[i*2+2]+(L1xy*invcosphi/2.));
		}else{
			VertPlaneEq[i*2+1] = -(VertPlaneSlope*final_xbox) + ((L1xy/2.) + i*(L2xy+L1xy))*invcosphi; // values of y for x = 0 
			VertPlaneEq[i*2+2] = -(VertPlaneSlope*final_xbox) + ((L1xy/2.) + L2xy + i*(L2xy+L1xy))*invcosphi;
			if( xpar ){
				FullPlaneEq[i*2+1] = Dir2_G1[0]*(VertPlaneEq[i*2+1]+(L2xy*invcosphi/2.));
				FullPlaneEq[i*2+2] = Dir1_G1[0]*(VertPlaneEq[i*2+2]+(L1xy*invcosphi/2.));
			}else{
				FullPlaneEq[i*2+1] = Dir2_G1[1]*(VertPlaneEq[i*2+1]+(L2xy*invcosphi/2.));
				FullPlaneEq[i*2+2] = Dir1_G1[1]*(VertPlaneEq[i*2+2]+(L1xy*invcosphi/2.));
			}
		}
		for(unsigned int d=0;d<3;d++){
			NormalPlaneEq[i*6+3+d] = Dir2_G1[d];
			NormalPlaneEq[i*6+6+d] = Dir1_G1[d];
		}
	}
	if( VertPlaneSlope < 0. ){
		if( xpar ){
			VertPlaneEq[NbPlanes-1] = final_xbox + 1.; // values of y for x = 0
			FullPlaneEq[NbPlanes-1] = final_xbox;
		}else{
			VertPlaneEq[NbPlanes-1] = final_ybox - VertPlaneSlope*final_xbox + 1.; // values of y for x = 0
			FullPlaneEq[NbPlanes-1] = final_ybox - VertPlaneSlope*final_xbox;
		}
	}else{
		if( xpar ){
			VertPlaneEq[NbPlanes-1] = final_xbox + 1.; 
			FullPlaneEq[NbPlanes-1] = Dir1_G1[0]*final_xbox;
		}else{
			VertPlaneEq[NbPlanes-1] = final_ybox + 1.; 
			FullPlaneEq[NbPlanes-1] = Dir1_G1[1]*final_ybox;
		}
	}
	for(unsigned int i=0;i<3;i++) NormalPlaneEq[(NbPlanes-1)*3+i] = Dir1_G1[i];
	
	// Compute planes eq for grain 2 according to the corresponding crystallographic planes
	// 1. search Dir1_G2, Dir2_G2 and new FacetJctDir 
	int h_c, k_c, l_c, u_c_1, v_c_1, w_c_1, u_c_2, v_c_2, w_c_2;
	_MyCrystal2->SearchCrystallographicPlanesAndDirections(Dir1_G1,h_c,k_c,l_c,u_c_1,v_c_1,w_c_1);
	_MyCrystal2->SearchCrystallographicPlanesAndDirections(Dir2_G1,h_c,k_c,l_c,u_c_2,v_c_2,w_c_2);
	_MyCrystal2->getA1()[0] *= Mx2;
	_MyCrystal2->getA1()[1] *= My2;
	_MyCrystal2->getA2()[0] *= Mx2;
	_MyCrystal2->getA2()[1] *= My2;
	_MyCrystal2->getA3()[0] *= Mx2;
	_MyCrystal2->getA3()[1] *= My2;
	_MyCrystal->computeReciproqual();
	facnorm1 = 0.;
	facnorm2 = 0.;
	for(unsigned int i=0;i<3;i++){
		Dir1_G2[i] = u_c_1*_MyCrystal2->getA1_star()[i]+v_c_1*_MyCrystal2->getA2_star()[i]+w_c_1*_MyCrystal2->getA3_star()[i];
		Dir2_G2[i] = u_c_2*_MyCrystal2->getA1_star()[i]+v_c_2*_MyCrystal2->getA2_star()[i]+w_c_2*_MyCrystal2->getA3_star()[i];
		facnorm1 += Dir1_G2[i]*Dir1_G2[i];
		facnorm2 += Dir2_G2[i]*Dir2_G2[i];
	}
	double dfac1_G2 = _MyCrystal2->ComputeD_hkl(u_c_1,v_c_1,w_c_1) / sin(alpha);
	double dfac2_G2 = _MyCrystal2->ComputeD_hkl(u_c_2,v_c_2,w_c_2) / sin(alpha);
	if( fabs(round(L_facet/dfac2_G2)-(L_facet/dfac2_G2)) > 1e-5 ) cout << "WARNING !! The length of facet 1 (" << u_c_1 << "_" << v_c_1 << "_" << w_c_1 << ") is not a multiple of the reticular distance of facet 2, some issues could occur ! (diff = " << fabs(round(L_facet/dfac2_G2)-(L_facet/dfac2_G2)) << ")" << endl;
	if( fabs(round(L_facet2/dfac1_G2)-(L_facet2/dfac1_G2)) > 1e-5 ) cout << "WARNING !! The length of facet 2 (" << u_c_2 << "_" << v_c_2 << "_" << w_c_2 << ") is not a multiple of the reticular distance of facet 1, some issues could occur ! (diff = " << fabs(round(L_facet2/dfac1_G2)-(L_facet2/dfac1_G2)) << ")" << endl;
	facnorm1 = sqrt(facnorm1);
	facnorm2 = sqrt(facnorm2);
	if( Dir1_G1[2] >= 0 ) for(unsigned int i=0;i<3;i++) Dir1_G1[i] *= -1.;
	else for(unsigned int i=0;i<3;i++) Dir2_G1[i] *= -1.;
	double sp1(0.), sp2(0.);
	for(unsigned int i=0;i<3;i++){
		sp1 += Dir1_G1[i]*Dir1_G2[i];
		sp2 += Dir2_G1[i]*Dir2_G2[i];
	}
	if( sp1 < 0. ) facnorm1 *= -1.;
	if( sp2 < 0. ) facnorm2 *= -1.;
	for(unsigned int i=0;i<3;i++){
		Dir1_G2[i] /= facnorm1;
		Dir2_G2[i] /= facnorm2;
	}
	double *FacetJctDir2 = new double[3];
	MT->mixedProd(Dir1_G2,Dir2_G2,FacetJctDir2);
	double facnorm3 = 0.;
	for(unsigned int i=0;i<3;i++) facnorm3 += FacetJctDir2[i]*FacetJctDir2[i];
	facnorm3 = sqrt(facnorm3);
	for(unsigned int i=0;i<3;i++) FacetJctDir2[i] /= facnorm3;

	// 2. Compute plane eq for grain 2
	double *FullPlaneEq_2 = new double[NbPlanes]; 
	double *NormalPlaneEq_2 = new double[NbPlanes*3]; // just store the normal of facet plane of each plane (redondant but allows to save time in the atom loop
	double *VertPlaneEq_2 = new double[NbPlanes];
	double VertPlaneSlope_2;
	bool xpar_2 = false;
	if( fabs(FacetJctDir2[0]) > zero_num ) VertPlaneSlope_2 = FacetJctDir2[1]/FacetJctDir2[0];
	else{
		VertPlaneSlope_2 = 0.;
		xpar = true;
	}
	double invcosphi_2 = 1./cos(atan(VertPlaneSlope_2));

	if( Dir1_G2[2] < 0 ) for(unsigned int i=0;i<3;i++) Dir1_G2[i] *= -1.;
	else for(unsigned int i=0;i<3;i++) Dir2_G2[i] *= -1.;

	if( VertPlaneSlope_2 < 0. ){
		VertPlaneEq_2[0] = -1.; // values of y for x = 0
		FullPlaneEq_2[0] = 0.;
	}else{
		VertPlaneEq_2[0] = - VertPlaneSlope_2*final_xbox - 1.; 
		FullPlaneEq_2[0] = - Dir1_G2[1]*VertPlaneSlope_2*final_xbox;
	}
	for(unsigned int i=0;i<3;i++) NormalPlaneEq_2[i] = Dir1_G2[i];

	for(unsigned int i=0;i<Nfacets;i++){
		if( VertPlaneSlope_2 < 0. ){
			VertPlaneEq_2[i*2+1] = ((L1xy/2.) + i*(L2xy+L1xy))*invcosphi_2; // values of y for x = 0
			VertPlaneEq_2[i*2+2] = ((L1xy/2.) + L2xy + i*(L2xy+L1xy))*invcosphi_2;
			FullPlaneEq_2[i*2+1] = Dir2_G2[1]*(VertPlaneEq_2[i*2+1]+(L2xy*invcosphi_2/2.)); // values of y for x = 0 and z = 0
			FullPlaneEq_2[i*2+2] = Dir1_G2[1]*(VertPlaneEq_2[i*2+2]+(L1xy*invcosphi_2/2.));
		}else{
			VertPlaneEq_2[i*2+1] = -(VertPlaneSlope_2*final_xbox) + ((L1xy/2.) + i*(L2xy+L1xy))*invcosphi_2; // values of y for x = 0 
			VertPlaneEq_2[i*2+2] = -(VertPlaneSlope_2*final_xbox) + ((L1xy/2.) + L2xy + i*(L2xy+L1xy))*invcosphi_2;
			if( xpar_2 ){
				FullPlaneEq_2[i*2+1] = Dir2_G2[0]*(VertPlaneEq_2[i*2+1]+(L2xy*invcosphi_2/2.));
				FullPlaneEq_2[i*2+2] = Dir1_G2[0]*(VertPlaneEq_2[i*2+2]+(L1xy*invcosphi_2/2.));
			}else{
				FullPlaneEq_2[i*2+1] = Dir2_G2[1]*(VertPlaneEq_2[i*2+1]+(L2xy*invcosphi_2/2.));
				FullPlaneEq_2[i*2+2] = Dir1_G2[1]*(VertPlaneEq_2[i*2+2]+(L1xy*invcosphi_2/2.));
			}
		}
		for(unsigned int d=0;d<3;d++){
			NormalPlaneEq_2[i*6+3+d] = Dir2_G2[d];
			NormalPlaneEq_2[i*6+6+d] = Dir1_G2[d];
		}
	}
	if( VertPlaneSlope_2 < 0. ){
		if( xpar_2 ){
			VertPlaneEq_2[NbPlanes-1] = final_xbox + 1.; // values of y for x = 0
			FullPlaneEq_2[NbPlanes-1] = final_xbox;
		}else{
			VertPlaneEq_2[NbPlanes-1] = final_ybox - VertPlaneSlope_2*final_xbox + 1.; // values of y for x = 0
			FullPlaneEq_2[NbPlanes-1] = final_ybox - VertPlaneSlope_2*final_xbox;
		}
	}else{
		if( xpar ){
			VertPlaneEq_2[NbPlanes-1] = final_xbox + 1.; 
			FullPlaneEq_2[NbPlanes-1] = Dir1_G2[0]*final_xbox;
		}else{
			VertPlaneEq_2[NbPlanes-1] = final_ybox + 1.; 
			FullPlaneEq_2[NbPlanes-1] = Dir1_G2[1]*final_ybox;
		}
	}
	for(unsigned int i=0;i<3;i++) NormalPlaneEq_2[(NbPlanes-1)*3+i] = Dir1_G2[i];
	
	// paste the two grains
	double xpos, ypos, zpos, Lin;

	unsigned int trueNbAt1 = 0;
	unsigned int trueNbAt2 = 0;
	unsigned int at_count = 0;
	double hG1_inf= zboxG1-(2.*DeltaZ);
	double hG2_sup= 2.*DeltaZ;
	double eps_pos = 5e-2;
	double eps_pos_z = 3e-1;
	vector<vector<int>> cur_NotSepTag1;
	vector<vector<int>> cur_NotSepTag2;
	// TODO manage bonds and angle
	if( !this->_MyCrystal->getIsDoNotSep() ){
		// for crystals where ions can be separated
		for(unsigned int i=0;i<dupX1;i++){
			for(unsigned int j=0;j<dupY1;j++){
				for(unsigned int k=0;k<zdup1;k++){
					for(unsigned int n=0;n<nbAtom1;n++){
						xpos = Mx1*((this->_MyCrystal->getOrientedSystem()->getAtom(n).pos.x)+i*xl1); // no need to account for tilted cell vectors as the oriented system in always orthogonal
						ypos = My1*((this->_MyCrystal->getOrientedSystem()->getAtom(n).pos.y)+j*yl1);
						zpos = (this->_MyCrystal->getOrientedSystem()->getAtom(n).pos.z)+k*zl1;
						if( zpos < hG1_inf ){
							AtomList_temp[at_count] = this->_MyCrystal->getOrientedSystem()->getAtom(n);
							AtomList_temp[at_count].pos.x = xpos;
							AtomList_temp[at_count].pos.y = ypos;
							AtomList_temp[at_count].pos.z = zpos;
							TagGrain[at_count] = 1;
							at_count += 1;
						}else{
							for(unsigned int f=0;f<NbPlanes-1;f++){
								// First test if atom is between two vertical planes
								double buffer_v;
							        if( xpar ) buffer_v = eps_pos + xpos;
								else buffer_v = eps_pos + ypos - VertPlaneSlope*xpos;
								if( ( buffer_v > VertPlaneEq[f] ) && ( buffer_v < VertPlaneEq[f+1] ) ){
									// Then test if atom is bellow the real plane
									double buffer_f = eps_pos_z + xpos*NormalPlaneEq[f*3] + ypos*NormalPlaneEq[f*3+1] + (zpos-zboxG1+DeltaZ)*NormalPlaneEq[f*3+2];
									if( buffer_f < FullPlaneEq[f] ){
										AtomList_temp[at_count] = this->_MyCrystal->getOrientedSystem()->getAtom(n);
										AtomList_temp[at_count].pos.x = xpos;
										AtomList_temp[at_count].pos.y = ypos;
										AtomList_temp[at_count].pos.z = zpos;
										TagGrain[at_count] = 1;
										at_count += 1;
									}
									break;
								}
							}
						}
					}
				}
			}
		}
		trueNbAt1 = at_count;
		for(unsigned int i=0;i<dupX2;i++){
			for(unsigned int j=0;j<dupY2;j++){
				for(unsigned int k=0;k<zdup2;k++){
					for(unsigned int n=0;n<nbAtom2;n++){
						xpos = Mx2*((this->_MyCrystal2->getOrientedSystem()->getAtom(n).pos.x)+i*xl2); // no need to account for tilted cell vectors as the oriented system in always orthogonal
						ypos = My2*((this->_MyCrystal2->getOrientedSystem()->getAtom(n).pos.y)+j*yl2);
						zpos = (this->_MyCrystal2->getOrientedSystem()->getAtom(n).pos.z)+k*zl2;
						if( zpos > hG2_sup ){
							AtomList_temp[at_count] = this->_MyCrystal2->getOrientedSystem()->getAtom(n);
							AtomList_temp[at_count].pos.x = xpos;
							AtomList_temp[at_count].pos.y = ypos;
							AtomList_temp[at_count].pos.z = zpos;
							TagGrain[at_count] = 2;
							at_count += 1;
						}else{
							for(unsigned int f=0;f<NbPlanes-1;f++){
								// First test if atom is between two vertical planes
								double buffer_v;
							        if( xpar ) buffer_v = eps_pos + xpos;
								else buffer_v = eps_pos + ypos - VertPlaneSlope*xpos;
								if( ( buffer_v > VertPlaneEq[f] ) && ( buffer_v < VertPlaneEq[f+1] ) ){
									// Then test if atom is above the real plane
									double buffer_f = eps_pos_z + xpos*NormalPlaneEq[f*3] + ypos*NormalPlaneEq[f*3+1] + (zpos-DeltaZ)*NormalPlaneEq[f*3+2];
									if( buffer_f > FullPlaneEq[f] ){
										AtomList_temp[at_count] = this->_MyCrystal2->getOrientedSystem()->getAtom(n);
										AtomList_temp[at_count].pos.x = xpos;
										AtomList_temp[at_count].pos.y = ypos;
										AtomList_temp[at_count].pos.z = zpos;
										TagGrain[at_count] = 2;
										at_count += 1;
									}
									break;
								}
							}
						}
					}
				}
			}
		}
		trueNbAt2 = at_count-trueNbAt1;
	}else{
		// do not sep case
		//_MyCrystal->getOrientedSystem()->ComputeNotSepList(); // no need to compute the list because the creation of the GB has been done
		//_MyCrystal2->getOrientedSystem()->ComputeNotSepList();
		vector<int>* NotSepTag1 = _MyCrystal->getOrientedSystem()->getNotSepTag();
		vector<int>* NotSepTag2 = _MyCrystal2->getOrientedSystem()->getNotSepTag();
		for(unsigned int i=0;i<dupX1;i++){
			for(unsigned int j=0;j<dupY1;j++){
				for(unsigned int k=0;k<zdup1;k++){
					for(unsigned int n=0;n<nbAtom1;n++){
						if( NotSepTag1[n][0] >= 0 ){
							xpos = Mx1*((this->_MyCrystal->getOrientedSystem()->getAtom(n).pos.x)+i*xl1); // no need to account for tilted cell vectors as the oriented system in always orthogonal
							ypos = My1*((this->_MyCrystal->getOrientedSystem()->getAtom(n).pos.y)+j*yl1);
							zpos = (this->_MyCrystal->getOrientedSystem()->getAtom(n).pos.z)+k*zl1;
							if( zpos < hG1_inf ){
								AtomList_temp[at_count] = this->_MyCrystal->getOrientedSystem()->getAtom(n);
								AtomList_temp[at_count].pos.x = xpos;
								AtomList_temp[at_count].pos.y = ypos;
								AtomList_temp[at_count].pos.z = zpos;
								TagGrain[at_count] = 1;
								cur_NotSepTag1.push_back(vector<int>());
								int cur_at_count = at_count;
								cur_NotSepTag1[at_count].push_back(NotSepTag1[n][0]);
								at_count += 1;
								for(unsigned int neigh=0;neigh<NotSepTag1[n][0];neigh++){
									cur_NotSepTag1[(unsigned int) cur_at_count].push_back(at_count);
									cur_NotSepTag1.push_back(vector<int>());
									cur_NotSepTag1[at_count].push_back(-1-cur_at_count);
									AtomList_temp[at_count] = this->_MyCrystal->getOrientedSystem()->getAtom(NotSepTag1[n][neigh+1]);
									AtomList_temp[at_count].pos.x += i*xl1;
									AtomList_temp[at_count].pos.y += j*yl1;
									AtomList_temp[at_count].pos.z += k*zl1;
									AtomList_temp[at_count].pos.x *= Mx1;
									AtomList_temp[at_count].pos.y *= My1;
									TagGrain[at_count] = 1;
									at_count += 1;
								}
							}else{
								for(unsigned int f=0;f<NbPlanes-1;f++){
									// First test if atom is between two vertical planes
									double buffer_v;
								        if( xpar ) buffer_v = eps_pos + xpos;
									else buffer_v = eps_pos + ypos - VertPlaneSlope*xpos;
									if( ( buffer_v > VertPlaneEq[f] ) && ( buffer_v < VertPlaneEq[f+1] ) ){
										// Then test if atom is bellow the real plane
										double buffer_f = eps_pos_z + xpos*NormalPlaneEq[f*3] + ypos*NormalPlaneEq[f*3+1] + (zpos-zboxG1+DeltaZ)*NormalPlaneEq[f*3+2];
										if( buffer_f < FullPlaneEq[f] ){
											cur_NotSepTag1.push_back(vector<int>());
											int cur_at_count = at_count;
											cur_NotSepTag1[at_count].push_back(NotSepTag1[n][0]);
											AtomList_temp[at_count] = this->_MyCrystal->getOrientedSystem()->getAtom(n);
											AtomList_temp[at_count].pos.x = xpos;
											AtomList_temp[at_count].pos.y = ypos;
											AtomList_temp[at_count].pos.z = zpos;
											TagGrain[at_count] = 1;
											at_count += 1;
											for(unsigned int neigh=0;neigh<NotSepTag1[n][0];neigh++){
												cur_NotSepTag1[(unsigned int) cur_at_count].push_back(at_count);
												cur_NotSepTag1.push_back(vector<int>());
												cur_NotSepTag1[at_count].push_back(-1-cur_at_count);
												AtomList_temp[at_count] = this->_MyCrystal->getOrientedSystem()->getAtom(NotSepTag1[n][neigh+1]);
												AtomList_temp[at_count].pos.x += i*xl1;
												AtomList_temp[at_count].pos.y += j*yl1;
												AtomList_temp[at_count].pos.z += k*zl1;
												AtomList_temp[at_count].pos.x *= Mx1;
												AtomList_temp[at_count].pos.y *= My1;
												TagGrain[at_count] = 1;
												at_count += 1;
											}
										}
										break;
									}
								}
							}
						} // end if NotSepTag >= 0
					} // end loop nbAtom1
				} // end loop z grain1
			} // end loop y grain1
		} // end loop x grain1
		trueNbAt1 = at_count;
		for(unsigned int i=0;i<dupX2;i++){
			for(unsigned int j=0;j<dupY2;j++){
				for(unsigned int k=0;k<zdup2;k++){
					for(unsigned int n=0;n<nbAtom2;n++){
						if( NotSepTag2[n][0] >= 0 ){
							xpos = Mx2*((this->_MyCrystal2->getOrientedSystem()->getAtom(n).pos.x)+i*xl2); // no need to account for tilted cell vectors as the oriented system in always orthogonal
							ypos = My2*((this->_MyCrystal2->getOrientedSystem()->getAtom(n).pos.y)+j*yl2);
							zpos = (this->_MyCrystal2->getOrientedSystem()->getAtom(n).pos.z)+k*zl2;
							if( zpos > hG2_sup ){
								AtomList_temp[at_count] = this->_MyCrystal2->getOrientedSystem()->getAtom(n);
								AtomList_temp[at_count].pos.x = xpos;
								AtomList_temp[at_count].pos.y = ypos;
								AtomList_temp[at_count].pos.z = zpos;
								TagGrain[at_count] = 2;
								cur_NotSepTag2.push_back(vector<int>());
								int cur_at_count = at_count - trueNbAt1;
								cur_NotSepTag2[(unsigned int) cur_at_count].push_back(NotSepTag2[n][0]);
								at_count += 1;
								for(unsigned int neigh=0;neigh<NotSepTag2[n][0];neigh++){
									unsigned int cur_at_count_n = at_count - trueNbAt1;
									cur_NotSepTag2[(unsigned int) cur_at_count].push_back(cur_at_count_n);
									cur_NotSepTag2.push_back(vector<int>());
									cur_NotSepTag2[cur_at_count_n].push_back(-1-cur_at_count);
									AtomList_temp[at_count] = this->_MyCrystal2->getOrientedSystem()->getAtom(NotSepTag2[n][neigh+1]);
									AtomList_temp[at_count].pos.x += i*xl2;
									AtomList_temp[at_count].pos.y += j*yl2;
									AtomList_temp[at_count].pos.z += k*zl2;
									AtomList_temp[at_count].pos.x *= Mx2;
									AtomList_temp[at_count].pos.y *= My2;
									TagGrain[at_count] = 2;
									at_count += 1;
								}
							}else{
								for(unsigned int f=0;f<NbPlanes-1;f++){
									// First test if atom is between two vertical planes
									double buffer_v;
								        if( xpar ) buffer_v = eps_pos + xpos;
									else buffer_v = eps_pos + ypos - VertPlaneSlope_2*xpos;
									if( ( buffer_v > VertPlaneEq_2[f] ) && ( buffer_v < VertPlaneEq_2[f+1] ) ){
										// Then test if atom is above the real plane
										double buffer_f = eps_pos_z + xpos*NormalPlaneEq_2[f*3] + ypos*NormalPlaneEq_2[f*3+1] + (zpos-DeltaZ)*NormalPlaneEq_2[f*3+2];
										if( buffer_f > FullPlaneEq_2[f] ){
											cur_NotSepTag2.push_back(vector<int>());
											int cur_at_count = at_count - trueNbAt1;
											cur_NotSepTag2[(unsigned int) cur_at_count].push_back(NotSepTag2[n][0]);
											AtomList_temp[at_count] = this->_MyCrystal2->getOrientedSystem()->getAtom(n);
											AtomList_temp[at_count].pos.x = xpos;
											AtomList_temp[at_count].pos.y = ypos;
											AtomList_temp[at_count].pos.z = zpos;
											TagGrain[at_count] = 2;
											at_count += 1;
											for(unsigned int neigh=0;neigh<NotSepTag2[n][0];neigh++){
												unsigned int cur_at_count_n = at_count - trueNbAt1;
												cur_NotSepTag2[(unsigned int) cur_at_count].push_back(cur_at_count_n);
												cur_NotSepTag2.push_back(vector<int>());
												cur_NotSepTag2[cur_at_count_n].push_back(-1-cur_at_count);
												AtomList_temp[at_count] = this->_MyCrystal2->getOrientedSystem()->getAtom(NotSepTag2[n][neigh+1]);
												AtomList_temp[at_count].pos.x += i*xl2;
												AtomList_temp[at_count].pos.y += j*yl2;
												AtomList_temp[at_count].pos.z += k*zl2;
												AtomList_temp[at_count].pos.x *= Mx2;
												AtomList_temp[at_count].pos.y *= My2;
												TagGrain[at_count] = 2;
												at_count += 1;
											}
										}
										break;
									}
								}
							}
						} // end if NotSepTag >= 0
					} // end loop nbAtom2
				} // end z grain2
			} // end y grain2
		} // end x grain2
		trueNbAt2 = at_count-trueNbAt1;
	} // end if DoNotSep case
	nbAtom1 = trueNbAt1;
	nbAtom2 = trueNbAt2;
	this->AtomList_G1 = new Atom[trueNbAt1];
	this->AtomList_G2 = new Atom[trueNbAt2];
	unsigned int countG1=0;
	unsigned int countG2=0;
	for(unsigned int i=0;i<at_count;i++){
		if( TagGrain[i] == 1 ){
			AtomList_G1[countG1] = AtomList_temp[i];
			countG1++;
		}else{
			AtomList_G2[countG2] = AtomList_temp[i];
			countG2++;
		}
	}

	// Initialize the two grains
	if( !this->_MyCrystal->getIsDoNotSep() ){
		this->Grain1 = new AtomicSystem(this->AtomList_G1,nbAtom1,_MyCrystal,this->H1_G1,this->H2_G1,this->H3_G1,true);
		this->Grain2 = new AtomicSystem(this->AtomList_G2,nbAtom2,_MyCrystal2,this->H1_G2,this->H2_G2,this->H3_G2,true);
	}else{
		this->Grain1 = new AtomicSystem(this->AtomList_G1,nbAtom1,_MyCrystal,this->H1_G1,this->H2_G1,this->H3_G1,cur_NotSepTag1,true);
		this->Grain2 = new AtomicSystem(this->AtomList_G2,nbAtom2,_MyCrystal2,this->H1_G2,this->H2_G2,this->H3_G2,cur_NotSepTag2,true);
	}
	this->AreGrainsDefined = true;
	if( IsCharge ){ 
		// idea here is to create oriented systems with the different free surfaces in presence
		// adjust their surfaces with the MakeSurfaceNeutral function (1d-based the only one working well)
		// rotate the systems and apply the given crystal def and misfits to align with Grain1 or 2
		// make coinciding the systems at the surfaces and delete atoms which not coinciding
		
		// ** Grain1 ** //
		// Bottom surface
		//Grain1->print_lmp("Grain1_1.lmp");
		vector<int> Oris;
		for(unsigned int i=0;i<3;i++) Oris.push_back(_MyCrystal->getOrthogonalPlanes()[6+i]);
		for(unsigned int i=0;i<3;i++) Oris.push_back(_MyCrystal->getOrthogonalPlanes()[i]);
		vector<double> shifts = {lowgrain_shift_x, lowgrain_shift_y, lowgrain_shift_z};
		vector<double> Misfits = {Mx1, My1, 1.};
		vector<double> currentPlaneNormal = {0., 0., -1.};
		double currentVertPlaneSlope = 0.;
		vector<double> currentVertPlaneEq = {-5.,this->H2_G1[1]+5.};
		vector<double> currentFullPlaneEq = {0.};
		double z_shift_box = 0.;
		cout << "Adjusting electroneutrality of grain 1 bottom surface" << endl;
		Grain1->MakeSurfaceNeutral(Oris, shifts, Misfits, currentPlaneNormal, currentVertPlaneSlope, currentVertPlaneEq, currentFullPlaneEq, z_shift_box, "Grain1_neutral.lmp");
		// Facet 1
		for(unsigned int i=0;i<3;i++) Oris[i] = FacetsType[i];
		if( Dir1_G1[2] < 0 ) for(unsigned int i=0;i<3;i++) currentPlaneNormal[i] = -1.*Dir1_G1[i];
		else for(unsigned int i=0;i<3;i++) currentPlaneNormal[i] = Dir1_G1[i];
		currentVertPlaneSlope = VertPlaneSlope;
		currentVertPlaneEq.clear();
		currentFullPlaneEq.clear();
		currentFullPlaneEq.push_back(FullPlaneEq[0]);
		currentVertPlaneEq.push_back(VertPlaneEq[0]-5.);
		currentVertPlaneEq.push_back(VertPlaneEq[1]);
		for(unsigned int i=0;i<(NbPlanes/2)-1;i++){
			currentVertPlaneEq.push_back(VertPlaneEq[i*2+2]);
			currentVertPlaneEq.push_back(VertPlaneEq[i*2+3]);
			currentFullPlaneEq.push_back(FullPlaneEq[i*2+2]);
		}
		currentVertPlaneEq[currentVertPlaneEq.size()-1] += 5.;
		//cout << "G1 facet1" << endl;
		//for(unsigned int i=0;i<currentFullPlaneEq.size();i++) cout << currentVertPlaneEq[i*2] << " " << currentVertPlaneEq[i*2+1] << " " << currentFullPlaneEq[i] << endl;
		z_shift_box = zboxG1-DeltaZ;
		cout << "Adjusting electroneutrality of grain 1 facet 1" << endl;
		Grain1->MakeSurfaceNeutral(Oris, shifts, Misfits, currentPlaneNormal, currentVertPlaneSlope, currentVertPlaneEq, currentFullPlaneEq, z_shift_box);
		
		// Facet 2
		for(unsigned int i=0;i<3;i++) Oris[i] = FacetsType[i+3];
		if( Dir2_G1[2] < 0 ) for(unsigned int i=0;i<3;i++) currentPlaneNormal[i] = -1.*Dir2_G1[i];
		else for(unsigned int i=0;i<3;i++) currentPlaneNormal[i] = Dir2_G1[i];
		currentVertPlaneEq.clear();
		currentFullPlaneEq.clear();
		currentFullPlaneEq.push_back(FullPlaneEq[1]);
		currentVertPlaneEq.push_back(VertPlaneEq[1]);
		currentVertPlaneEq.push_back(VertPlaneEq[2]);
		for(unsigned int i=0;i<(NbPlanes/2)-2;i++){
			currentVertPlaneEq.push_back(VertPlaneEq[i*2+3]);
			currentVertPlaneEq.push_back(VertPlaneEq[i*2+4]);
			currentFullPlaneEq.push_back(FullPlaneEq[i*2+3]);
		}
		//cout << "G1 facet2" << endl;
		//for(unsigned int i=0;i<currentFullPlaneEq.size();i++) cout << currentVertPlaneEq[i*2] << " " << currentVertPlaneEq[i*2+1] << " " << currentFullPlaneEq[i] << endl;
		z_shift_box = zboxG1-DeltaZ;
		cout << "Adjusting electroneutrality of grain 1 facet 2" << endl;
		Grain1->MakeSurfaceNeutral(Oris, shifts, Misfits, currentPlaneNormal, currentVertPlaneSlope, currentVertPlaneEq, currentFullPlaneEq, z_shift_box);

		nbAtom1 = Grain1->getNbAtom();
		
		// check if the system if charge neutral
		double total_charge_1 = 0.;
		for(unsigned int i=0;i<nbAtom1;i++)
			total_charge_1 += AtomCharge[Grain1->getAtom(i).type_uint-1];
		//cout << "Total charge 1 = " << total_charge_1 << endl;
		//Grain1->print_lmp("Grain1_2.lmp");
		if( fabs(total_charge_1) > 1e-8 ){
			cout << "Adjusting total electroneutrality of grain 1" << endl;
			unsigned int nbInters = (NbPlanes/2);
			vector<double> dirPoint(3*nbInters);
			double n1n1(0.), n2n2(0.), n1n2(0.);
			if( Dir1_G1[2] < 0 ) for(unsigned int i=0;i<3;i++) Dir1_G1[i] *= -1.;
			else for(unsigned int i=0;i<3;i++) Dir2_G1[i] *= -1.;
			for(unsigned int i=0;i<3;i++){
				n1n1 += Dir1_G1[i]*Dir1_G1[i];
				n2n2 += Dir2_G1[i]*Dir2_G1[i];
				n1n2 += Dir1_G1[i]*Dir2_G1[i];
			}
			for(unsigned int i=0;i<nbInters;i++){
				double d1 = FullPlaneEq[i];
				double d2 = FullPlaneEq[i+1];
				double *current_d1, *current_d2;
				double firstTerm, secondTerm;
				if( i%2 == 0 ){
					current_d1 = Dir1_G1;
					current_d2 = Dir2_G1;
					firstTerm = ( d1*n2n2 - d2*n1n2 ) / ( n1n1*n2n2 - n1n2*n1n2 );
					secondTerm = ( d2*n1n1 - d1*n1n2 ) / ( n1n1*n2n2 - n1n2*n1n2 );
				}else{
					current_d1 = Dir2_G1;
					current_d2 = Dir1_G1;
					firstTerm = ( d1*n1n1 - d2*n1n2 ) / ( n1n1*n2n2 - n1n2*n1n2 );
					secondTerm = ( d2*n2n2 - d1*n1n2 ) / ( n1n1*n2n2 - n1n2*n1n2 );
					firstTerm = ( FullPlaneEq[i]*n1n1 - FullPlaneEq[i+1]*n1n2 ) / ( n1n1*n2n2 - n1n2*n1n2 );
				}
				for(unsigned int d=0;d<3;d++) dirPoint[i*3+d] = firstTerm * current_d1[d] + secondTerm * current_d2[d];
				dirPoint[i*3+2] += zboxG1-DeltaZ;
			}
			vector<double> DirVec;
			for(unsigned int i=0;i<nbInters;i++)
				for(unsigned int d=0;d<3;d++) DirVec.push_back(FacetJctDir[d]);

			double mindist2jct = 10.;
			Grain1->getH3()[2] *= 2.;
                	Grain1->computeInverseCellVec();
			Grain1->MakeSystemNeutral(dirPoint,DirVec,mindist2jct);
                	Grain1->getH3()[2] /= 2.;
                	Grain1->computeInverseCellVec();
						
		}
		
		nbAtom1 = Grain1->getNbAtom();
		total_charge_1 = 0.;
		for(unsigned int i=0;i<nbAtom1;i++)
			total_charge_1 += AtomCharge[Grain1->getAtom(i).type_uint-1];
		cout << "Total charge grain 1 = " << total_charge_1 << endl;


		// ** Grain2 ** //
		// Bottom surface
		for(unsigned int i=0;i<3;i++) Oris[i] = _MyCrystal2->getOrthogonalPlanes()[6+i];
		for(unsigned int i=0;i<3;i++) Oris[i+3] = _MyCrystal2->getOrthogonalPlanes()[i];
		shifts = {upgrain_shift_x, upgrain_shift_y, upgrain_shift_z};
		Misfits = {Mx2, My2, 1.};
		currentPlaneNormal = {0., 0., 1.};
		currentVertPlaneSlope = 0.;
		currentVertPlaneEq = {-5.,this->H2_G2[1]+5.};
		currentFullPlaneEq = {this->H3_G2[2]-20};
		//cout << "G2 surface" << endl;
		//for(unsigned int i=0;i<currentFullPlaneEq.size();i++) cout << currentVertPlaneEq[i*2] << " " << currentVertPlaneEq[i*2+1] << " " << currentFullPlaneEq[i] << endl;
		z_shift_box = 0.;
		cout << "Adjusting electroneutrality of grain 2 bottom surface" << endl;
		Grain2->MakeSurfaceNeutral(Oris, shifts, Misfits, currentPlaneNormal, currentVertPlaneSlope, currentVertPlaneEq, currentFullPlaneEq, z_shift_box, "Grain2_neutral.lmp");
		// Facet 1
		Oris[0] = u_c_1;
		Oris[1] = v_c_1;
		Oris[2] = w_c_1;
		if( Dir1_G2[2] > 0. ) for(unsigned int i=0;i<3;i++) currentPlaneNormal[i] = -1.*Dir1_G2[i];
		else for(unsigned int i=0;i<3;i++) currentPlaneNormal[i] = Dir1_G2[i];
		currentVertPlaneSlope = VertPlaneSlope_2;
		currentVertPlaneEq.clear();
		currentFullPlaneEq.clear();
		currentFullPlaneEq.push_back(-FullPlaneEq_2[0]);
		currentVertPlaneEq.push_back(VertPlaneEq_2[0]-5.);
		currentVertPlaneEq.push_back(VertPlaneEq_2[1]);
		for(unsigned int i=0;i<(NbPlanes/2)-1;i++){
			currentVertPlaneEq.push_back(VertPlaneEq_2[i*2+2]);
			currentVertPlaneEq.push_back(VertPlaneEq_2[i*2+3]);
			currentFullPlaneEq.push_back(-FullPlaneEq_2[i*2+2]);
		}
		currentVertPlaneEq[currentVertPlaneEq.size()-1] += 5.;
		//cout << "G2 facet1" << endl;
		//for(unsigned int i=0;i<currentFullPlaneEq.size();i++) cout << currentVertPlaneEq[i*2] << " " << currentVertPlaneEq[i*2+1] << " " << currentFullPlaneEq[i] << endl;
		z_shift_box = DeltaZ;
		cout << "Adjusting electroneutrality of grain 2 facet 1" << endl;
		Grain2->MakeSurfaceNeutral(Oris, shifts, Misfits, currentPlaneNormal, currentVertPlaneSlope, currentVertPlaneEq, currentFullPlaneEq, z_shift_box);

		// Facet 2
		Oris[0] = u_c_2;
		Oris[1] = v_c_2;
		Oris[2] = w_c_2;
		if( Dir2_G2[2] > 0 ) for(unsigned int i=0;i<3;i++) currentPlaneNormal[i] = -1.*Dir2_G2[i];
		else for(unsigned int i=0;i<3;i++) currentPlaneNormal[i] = Dir2_G2[i];
		currentVertPlaneEq.clear();
		currentFullPlaneEq.clear();
		currentFullPlaneEq.push_back(-FullPlaneEq_2[1]);
		currentVertPlaneEq.push_back(VertPlaneEq_2[1]);
		currentVertPlaneEq.push_back(VertPlaneEq_2[2]);
		for(unsigned int i=0;i<(NbPlanes/2)-2;i++){
			currentVertPlaneEq.push_back(VertPlaneEq_2[i*2+3]);
			currentVertPlaneEq.push_back(VertPlaneEq_2[i*2+4]);
			currentFullPlaneEq.push_back(-FullPlaneEq_2[i*2+3]);
		}
		//cout << "G2 facet2" << endl;
		//for(unsigned int i=0;i<currentFullPlaneEq.size();i++) cout << currentVertPlaneEq[i*2] << " " << currentVertPlaneEq[i*2+1] << " " << currentFullPlaneEq[i] << endl;
		cout << "Adjusting electroneutrality of grain 2 facet 2" << endl;
		Grain2->MakeSurfaceNeutral(Oris, shifts, Misfits, currentPlaneNormal, currentVertPlaneSlope, currentVertPlaneEq, currentFullPlaneEq, z_shift_box);
		nbAtom2 = Grain2->getNbAtom();

		// check if the system if charge neutral
		nbAtom2 = Grain2->getNbAtom();
		double total_charge_2 = 0.;
		for(unsigned int i=0;i<nbAtom2;i++)
			total_charge_2 += AtomCharge[Grain2->getAtom(i).type_uint-1];
		//cout << "Total charge 2 = " << total_charge_2 << endl;
		if( fabs(total_charge_2) > 1e-8 ){
			cout << "Adjusting total electroneutrality of grain 2" << endl;
			unsigned int nbInters = (NbPlanes/2);
			vector<double> dirPoint(3*nbInters);
			double n1n1(0.), n2n2(0.), n1n2(0.);
			for(unsigned int i=0;i<3;i++){
				n1n1 += Dir1_G2[i]*Dir1_G2[i];
				n2n2 += Dir2_G2[i]*Dir2_G2[i];
				n1n2 += Dir1_G2[i]*Dir2_G2[i];
			}
			for(unsigned int i=0;i<nbInters;i++){
				double d1 = FullPlaneEq[i];
				double d2 = FullPlaneEq[i+1];
				double *current_d1, *current_d2;
				double firstTerm, secondTerm;
				if( i%2 == 0 ){
					current_d1 = Dir1_G2;
					current_d2 = Dir2_G2;
					firstTerm = ( d1*n2n2 - d2*n1n2 ) / ( n1n1*n2n2 - n1n2*n1n2 );
					secondTerm = ( d2*n1n1 - d1*n1n2 ) / ( n1n1*n2n2 - n1n2*n1n2 );
				}else{
					current_d1 = Dir2_G2;
					current_d2 = Dir1_G2;
					firstTerm = ( d1*n1n1 - d2*n1n2 ) / ( n1n1*n2n2 - n1n2*n1n2 );
					secondTerm = ( d2*n2n2 - d1*n1n2 ) / ( n1n1*n2n2 - n1n2*n1n2 );
					firstTerm = ( FullPlaneEq[i]*n1n1 - FullPlaneEq[i+1]*n1n2 ) / ( n1n1*n2n2 - n1n2*n1n2 );
				}
				for(unsigned int d=0;d<3;d++) dirPoint[i*3+d] = firstTerm * current_d1[d] + secondTerm * current_d2[d];
				dirPoint[i*3+2] += DeltaZ;
			}
			vector<double> DirVec;
			for(unsigned int i=0;i<nbInters;i++)
				for(unsigned int d=0;d<3;d++) DirVec.push_back(FacetJctDir2[d]);

			double mindist2jct = 10.;
			Grain2->getH3()[2] *= 2.;
                	Grain2->computeInverseCellVec();
			Grain2->MakeSystemNeutral(dirPoint,DirVec,mindist2jct);
                	Grain2->getH3()[2] /= 2.;
                	Grain2->computeInverseCellVec();
			
			nbAtom2 = Grain2->getNbAtom();
			total_charge_2 = 0.;
			for(unsigned int i=0;i<nbAtom2;i++)
				total_charge_2 += AtomCharge[Grain2->getAtom(i).type_uint-1];
			cout << "Total charge grain 2 = " << total_charge_2 << endl;
		
		}
		cout << "Sum = " << total_charge_1+total_charge_2 << endl;


	} // end if adjust surfaces

	// merge the two grains 
	this->nbAtom = Grain1->getNbAtom() + Grain2->getNbAtom();
	this->AtomList = new Atom[this->nbAtom];
	for(unsigned int i=0;i<Grain1->getNbAtom();i++) this->AtomList[i] = Grain1->getAtom(i);
	for(unsigned int i=0;i<Grain2->getNbAtom();i++){
		this->AtomList[i+Grain1->getNbAtom()] = Grain2->getAtom(i);
		this->AtomList[i+Grain1->getNbAtom()].pos.z += (GBspace/2.)-(2.*DeltaZ)+zboxG1;
	}

	SetSystemSizeAccordingToFixedParam(true,true,false);	

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
	string nbAt1_str = to_string(nbAtom1);
	string nbAt2_str = to_string(nbAtom2);
	this->File_Heading = " # ["+h_a_str+k_a_str+l_a_str+"]"+theta_str+"°("+h_p_str+k_p_str+l_p_str+") "+crystalName+" grain boundary\n # The grain boundary is faceted with ("+to_string(FacetsType[0])+"_"+to_string(FacetsType[1])+"_"+to_string(FacetsType[2])+") and ("+to_string(FacetsType[3])+"_"+to_string(FacetsType[4])+"_"+to_string(FacetsType[5])+") GB planes\n";
	this->File_Heading += " # The lower grain (Grain1) contains "+nbAt1_str+" atoms\n";
	this->File_Heading += " # The upper grain (Grain2) contains "+nbAt2_str+" atoms\n";
	this->File_Heading += " # ("+to_string(FacetsType[0])+"_"+to_string(FacetsType[1])+"_"+to_string(FacetsType[2])+") facet has length of "+to_string(L_facet)+" A, for total area of "+to_string(double (Nfacets)*FacetLength*L_facet)+" A2\n";
	this->File_Heading += " # ("+to_string(FacetsType[3])+"_"+to_string(FacetsType[4])+"_"+to_string(FacetsType[5])+") facet has length of "+to_string(L_facet2)+" A, for total area of "+to_string(double (Nfacets)*FacetLength*L_facet2)+" A2\n";
	this->Grain1->set_File_Heading(" # Lower grain of the ["+h_a_str+k_a_str+l_a_str+"]"+theta_str+"°("+h_p_str+k_p_str+l_p_str+") "+crystalName+" grain boundary\n");
	this->Grain2->set_File_Heading(" # Upper grain of the ["+h_a_str+k_a_str+l_a_str+"]"+theta_str+"°("+h_p_str+k_p_str+l_p_str+") "+crystalName+" grain boundary\n");

	this->xl1 *= Mx1; 
	this->xl2 *= Mx2; 
	this->yl1 *= My1; 
	this->yl2 *= My2; 
	// delete allocated memory
	delete[] AtomList_temp;
	delete[] Dir1_G1;
	delete[] Dir2_G1;
	delete[] Dir1_G2;
	delete[] Dir2_G2;
	delete[] TagGrain;
	delete[] l1;
	delete[] l2;
	delete[] FacetJctDir;
	delete[] FullPlaneEq;
	delete[] NormalPlaneEq;
	delete[] VertPlaneEq;
	delete[] FacetJctDir2;
	delete[] FullPlaneEq_2;
	delete[] NormalPlaneEq_2;
	delete[] VertPlaneEq_2;
}

void Bicrystal::ShowPossibleFacets(unsigned int max_hkl_u){ // TODO will not work with cubic primitive cell
	if( max_hkl_u == 0 ) return;
	double zero_num = 1.e-9;
	int max_hkl = (int) max_hkl_u;
	vector<vector<int>> hklcombs = MT->GenerateHKLCombinations(-max_hkl,max_hkl);
	vector<int> AdmissibleFacets;
	double *Dir1_G1 = new double[3];
	double *Dir2_G1 = new double[3];
	double facnorm1(0.), facnorm2(0.), sp_t(0.);
	double *FacetJctDir = new double[3];
	for(unsigned int i=0;i<hklcombs.size();i++){
		for(unsigned int j=i;j<hklcombs.size();j++){
			facnorm1 = 0.;
			facnorm2 = 0.;
			sp_t = 0.;
			for(unsigned int d=0;d<3;d++){
				Dir1_G1[d] = hklcombs[i][0]*_MyCrystal->getA1_star()[d]+hklcombs[i][1]*_MyCrystal->getA2_star()[d]+hklcombs[i][2]*_MyCrystal->getA3_star()[d];
				Dir2_G1[d] = hklcombs[j][0]*_MyCrystal->getA1_star()[d]+hklcombs[j][1]*_MyCrystal->getA2_star()[d]+hklcombs[j][2]*_MyCrystal->getA3_star()[d];
				facnorm1 += Dir1_G1[d]*Dir1_G1[d];
				facnorm2 += Dir2_G1[d]*Dir2_G1[d];
			}
			for(unsigned int d=0;d<2;d++) sp_t += Dir1_G1[d]*Dir2_G1[d];
			facnorm1 = sqrt(facnorm1);
			facnorm2 = sqrt(facnorm2);
			if( sp_t < 0. ){
				for(unsigned int d=0;d<3;d++) hklcombs[j][d] *= -1;
				facnorm2 *= -1.;
			}
			for(unsigned int d=0;d<3;d++){
				Dir1_G1[d] /= facnorm1;
				Dir2_G1[d] /= facnorm2;
			}
			MT->mixedProd(Dir1_G1,Dir2_G1,FacetJctDir);
			double abs_xf = fabs(FacetJctDir[0]);
			double abs_yf = fabs(FacetJctDir[1]);
			double abs_zf = fabs(FacetJctDir[2]);
 			if( ( Dir2_G1[0] != 0. || Dir2_G1[1] != 0. ) && ( Dir1_G1[0] != 0. || Dir1_G1[1] != 0. ) && ( ( Dir1_G1[2] < 0. && Dir2_G1[2] >= 0 ) || ( Dir1_G1[2] <= 0. && Dir2_G1[2] > 0 ) || ( Dir1_G1[2] > 0. && Dir2_G1[2] <= 0 ) || ( Dir1_G1[2] >= 0. && Dir2_G1[2] < 0 ) ) && ( abs_xf > zero_num || abs_yf > zero_num ) && abs_zf < zero_num ){
				AdmissibleFacets.push_back(hklcombs[i][0]);
				AdmissibleFacets.push_back(hklcombs[i][1]);
				AdmissibleFacets.push_back(hklcombs[i][2]);
				AdmissibleFacets.push_back(hklcombs[j][0]);
				AdmissibleFacets.push_back(hklcombs[j][1]);
				AdmissibleFacets.push_back(hklcombs[j][2]);
			}
		}
	}
	if( AdmissibleFacets.size() == 0 ){
		cout << "No facets are possible for this GB (computed for h, k and l between -" << max_hkl << " and " << max_hkl << ")" << endl;
		return;
	}
	cout << "Admissible facets (computed for h, k and l between -" << max_hkl << " and " << max_hkl << ") for this GB are:" << endl;
	//cout << "(" << AdmissibleFacets[0] << "_" << AdmissibleFacets[1] << "_" << AdmissibleFacets[2] << ")/(" << AdmissibleFacets[3] << "_" << AdmissibleFacets[4] << "_" << AdmissibleFacets[5] << ")";
	//for(unsigned int i=1;i<AdmissibleFacets.size()/6;i++)
	//	cout << ", or (" << AdmissibleFacets[i*6] << "_" << AdmissibleFacets[i*6+1] << "_" << AdmissibleFacets[i*6+2] << ")/(" << AdmissibleFacets[i*6+3] << "_" << AdmissibleFacets[i*6+4] << "_" << AdmissibleFacets[i*6+5] << ")";
	//cout << endl;
	cout << AdmissibleFacets[0] << " " << AdmissibleFacets[1] << " " << AdmissibleFacets[2] << "  " << AdmissibleFacets[3] << " " << AdmissibleFacets[4] << " " << AdmissibleFacets[5] << endl;
	for(unsigned int i=1;i<AdmissibleFacets.size()/6;i++)
		cout << "" << AdmissibleFacets[i*6] << " " << AdmissibleFacets[i*6+1] << " " << AdmissibleFacets[i*6+2] << "  " << AdmissibleFacets[i*6+3] << " " << AdmissibleFacets[i*6+4] << " " << AdmissibleFacets[i*6+5] << endl;

	delete[] Dir1_G1;
	delete[] Dir2_G1;
	delete[] FacetJctDir;
}

// Constructor for bicrystal with plane GB with given misorientation and GB plane
Bicrystal::Bicrystal(const string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p, bool rationalize, vector<string> Properties, int h_px, int k_px, int l_px, bool vacuum):h_a(h_a), k_a(k_a), l_a(l_a), theta(theta), h_p(h_p), k_p(k_p), l_p(l_p){
	read_params();
	this->MT = new MathTools;
	BT = new BicrystalloTools;
	InitCSL_DSC();
	setCrystal(crystalName);
	ReadProperties(Properties);
	this->_MyCrystal->ReadProperties(Properties);
	string i_a_str(""), i_p_str("");
	if( _MyCrystal->getCrystallo() == "Hexagonal" ){
		i_a_str = "_"+to_string(-(this->h_a)-(this->k_a));
		i_p_str = "_"+to_string(-(this->h_p)-(this->k_p));
	}
	cout << "Constructing a GB misoriented by " << theta*180./M_PI << "° around [" << h_a << "_" << k_a << i_a_str << "_" << l_a << "] axis and one GB plane (" << h_p << "_" << k_p << i_p_str << "_" << l_p << ")" << endl;
	setOrientedCrystals(crystalName, rationalize, Properties, h_px, k_px, l_px);
	if( vacuum && _MyCrystal->getIsCharge() ){
		cout << "Making surface charge neutral" << endl;
		_MyCrystal->getOrientedSystem()->MakeSurfaceNeutral(false);
		_MyCrystal2->getOrientedSystem()->MakeSurfaceNeutral(false);
	}
	_MyCrystal->ComputeOrthogonalPlanesAndDirections();
	_MyCrystal2->ComputeOrthogonalPlanesAndDirections();
	Dis->DisplayGB(_MyCrystal,_MyCrystal2);
	cout << endl;
	// search the number of linear combination for the two system to have the same x y length
	this->xl1 = this->_MyCrystal->getOrientedSystem()->getH1()[0];
	this->xl2 = this->_MyCrystal2->getOrientedSystem()->getH1()[0];
	this->yl1 = this->_MyCrystal->getOrientedSystem()->getH2()[1];
	this->yl2 = this->_MyCrystal2->getOrientedSystem()->getH2()[1];
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
	if( !find ){
		cerr << "We do not find any combination of the two grains in X direction to have a misfit lower than provided in the fixed parameters" << endl;
		cerr << "Increase MAX_MISFIT or MAX_DUP parameters in a FixedParameters.ath file" << endl;
		exit(EXIT_FAILURE);
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
	if( !find ){
		cerr << "We do not find any combination of the two grains in Y direction to have a misfit lower than provided by the fixed parameters" << endl;
		cerr << "Increase MAX_MISFIT or MAX_DUP parameters in a FixedParameters.ath file" << endl;
		exit(EXIT_FAILURE);
	}

	double final_xbox = (xl1*dupX1 + xl2*dupX2)/2;
    	double final_ybox = (yl1*dupY1 + yl2*dupY2)/2;
    	Mx1 = final_xbox / (xl1*dupX1);
    	My1 = final_ybox / (yl1*dupY1);
    	Mx2 = final_xbox / (xl2*dupX2);
    	My2 = final_ybox / (yl2*dupY2);
	cout << "Pasting the two grains over each others leads to:" << endl;
    	vector<vector<string>> arr_elements = {
		{"", "Grain 1", "Grain 2"}, 
		{"Direction", "X", "Y", "X", "Y"},
		{"Duplicates", to_string(dupX1), to_string(dupY1), to_string(dupX2), to_string(dupY2)},
		{"Misfits", to_string((Mx1-1.)*100.)+" %", to_string((My1-1.)*100.)+" %", to_string((Mx2-1.)*100.)+" %", to_string((My2-1.)*100.)+" %"}
	};
	vector<vector<unsigned int>> arr_fusion = {
		{1,  2, 2},
		{1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1},
		{1, 1, 1, 1, 1}
    	};
	Dis->DisplayArray(arr_elements, arr_fusion);
	cout << "Misfit values could be reduced by decreasing MAX_MISFIT in a FixedParameters.ath file" << endl;
	cout << endl;

	// initialize the variables and pointers
	unsigned int nbAtom1 = this->_MyCrystal->getOrientedSystem()->getNbAtom();
	unsigned int nbAtom2 = this->_MyCrystal2->getOrientedSystem()->getNbAtom();
	unsigned int nbAtom1_G, nbAtom2_G;
	this->nbAtom = nbAtom1*dupX1*dupY1 + nbAtom2*dupX2*dupY2;
	if( nbAtom > 1e6 && FullGrains ) cout << "WARNING: the number of atoms for this GB is high (>1e6), to save memory if needed you can switch the FULL_GRAINS parameter to 0 in a FixedParameters.ath file" << endl;
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
	this->AtomType = _MyCrystal->getAtomType();
	this->AtomMass = _MyCrystal->getAtomMass();
	this->AtomCharge = _MyCrystal->getAtomCharge();
	this->IsCharge = _MyCrystal->getIsCharge();
	this->IsTilted = false;
	unsigned int nbMol_G1, nbMol_G2;
	if( _MyCrystal->getIsMolId() ){
		IsMolId = true;
		MolId = new unsigned int[nbAtom];
		MolId_G1 = new unsigned int[nbAtom1_G];
		MolId_G2 = new unsigned int[nbAtom2_G];
		nbMol_G1 = MT->max_p(_MyCrystal->getOrientedSystem()->getMolId(),nbAtom1);
		nbMol_G2 = MT->max_p(_MyCrystal2->getOrientedSystem()->getMolId(),nbAtom2);
	}
	if( _MyCrystal->getIsBond() ){
		IsBond = true;
		nbBondType = _MyCrystal->getNbBondType();
		nbBondType_G1 = _MyCrystal->getNbBondType();
		nbBondType_G2 = _MyCrystal2->getNbBondType();
		nbBonds_G1 = _MyCrystal->getOrientedSystem()->getNbBonds();
		nbBonds_G2 = _MyCrystal2->getOrientedSystem()->getNbBonds();
		nbBonds = (nbBonds_G1*dupX1*dupY1) + (nbBonds_G2*dupX2*dupY2);
		Bonds = new unsigned int[nbBonds*2];
		BondType = new unsigned int[nbBonds];
		if( FullGrains ){
			Bonds_G1 = new unsigned int[nbBonds_G1*2*dupX1*dupY1];
			BondType_G1 = new unsigned int[nbBonds_G1*dupX1*dupY1];
			Bonds_G2 = new unsigned int[nbBonds_G2*2*dupX2*dupY2];
			BondType_G2 = new unsigned int[nbBonds_G2*dupX2*dupY2];
		}else{
			Bonds_G1 = new unsigned int[nbBonds_G1*2];
			BondType_G1 = new unsigned int[nbBonds_G1];
			Bonds_G2 = new unsigned int[nbBonds_G2*2];
			BondType_G2 = new unsigned int[nbBonds_G2];
		}
	}
	if( _MyCrystal->getIsAngle() ){
		IsAngle = true;
		nbAngleType = _MyCrystal->getNbAngleType();
		nbAngleType_G1 = _MyCrystal->getNbAngleType();
		nbAngleType_G2 = _MyCrystal2->getNbAngleType();
		nbAngles_G1 = _MyCrystal->getOrientedSystem()->getNbAngles();
		nbAngles_G2 = _MyCrystal2->getOrientedSystem()->getNbAngles();
		nbAngles = (nbAngles_G1*dupX1*dupY1) + (nbAngles_G2*dupX2*dupY2);
		Angles = new unsigned int[nbAngles*3];
		AngleType = new unsigned int[nbAngles];
		if( FullGrains ){
			Angles_G1 = new unsigned int[nbAngles_G1*3*dupX1*dupY1];
			AngleType_G1 = new unsigned int[nbAngles_G1*dupX1*dupY1];
			Angles_G2 = new unsigned int[nbAngles_G2*3*dupX2*dupY2];
			AngleType_G2 = new unsigned int[nbAngles_G2*dupX2*dupY2];
		}else{
			Angles_G1 = new unsigned int[nbAngles_G1*3];
			AngleType_G1 = new unsigned int[nbAngles_G1];
			Angles_G2 = new unsigned int[nbAngles_G2*3];
			AngleType_G2 = new unsigned int[nbAngles_G2];
		}
	}

	computeInverseCellVec();
	// paste the two grains
	for(unsigned int i=0;i<dupX1;i++){
		for(unsigned int j=0;j<dupY1;j++){
			for(unsigned int n=0;n<nbAtom1;n++){
				this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n] = this->_MyCrystal->getOrientedSystem()->getAtom(n);
				this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.x = Mx1*(this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.x + i*xl1);
				this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.y = My1*(this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.y + j*yl1);
				if( !FullGrains && i == 0 && j == 0 ){
					this->AtomList_G1[n] = this->AtomList[n];
					if( IsMolId ){
						MolId_G1[n] = (this->_MyCrystal->getOrientedSystem()->getMolId(n));
						MolId[n] = MolId_G1[n];
					}
				}else if( FullGrains ){
					this->AtomList_G1[i*dupY1*nbAtom1+j*nbAtom1+n] = this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n];
					if( IsMolId ){
						MolId_G1[i*dupY1*nbAtom1+j*nbAtom1+n] = (this->_MyCrystal->getOrientedSystem()->getMolId(n)) + i*dupY1*nbMol_G1 + j*nbMol_G1;
						MolId[i*dupY1*nbAtom1+j*nbAtom1+n] = MolId_G1[i*dupY1*nbAtom1+j*nbAtom1+n];
					}
				}else if( IsMolId ) MolId[i*dupY1*nbAtom1+j*nbAtom1+n] = (this->_MyCrystal->getOrientedSystem()->getMolId(n)) + i*dupY1*nbMol_G1 + j*nbMol_G1;
				this->AtomList[i*dupY1*nbAtom1+j*nbAtom1+n].pos.z += this->H3_G2[2]+(GBspace/2.);
				
			}
			if( IsBond ){
				if( !FullGrains && i == 0 && j == 0 ){
					for(unsigned int n=0;n<nbBonds_G1;n++){
						for(unsigned int d=0;d<2;d++){
							Bonds_G1[n*2+d] = _MyCrystal->getOrientedSystem()->getBonds()[n*2+d];
							Bonds[n*2+d] = Bonds_G1[n*2+d];
						}
						BondType_G1[n] = _MyCrystal->getOrientedSystem()->getBondType(n);
						BondType[n] = BondType_G1[n];
					}
				}else if( FullGrains ){
					unsigned int c_ind = i*dupY1*nbBonds_G1+j*nbBonds_G1;
					unsigned int c_ind_a = i*dupY1*nbAtom1+j*nbAtom1;
					for(unsigned int n=0;n<nbBonds_G1;n++){
						unsigned int c_ind_n = c_ind+n;
						for(unsigned int d=0;d<2;d++){
							Bonds_G1[c_ind_n*2+d] = _MyCrystal->getOrientedSystem()->getBonds()[n*2+d] + c_ind_a;
							Bonds[c_ind_n*2+d] = Bonds_G1[c_ind_n*2+d];
						}
						BondType_G1[c_ind_n] = _MyCrystal->getOrientedSystem()->getBondType(n);
						BondType[c_ind_n] = BondType_G1[c_ind_n];
					}
				}else{
					unsigned int c_ind = i*dupY1*nbBonds_G1+j*nbBonds_G1;
					unsigned int c_ind_a = i*dupY1*nbAtom1+j*nbAtom1;
					for(unsigned int n=0;n<nbBonds_G1;n++){
						unsigned int c_ind_n = c_ind+n;
						for(unsigned int d=0;d<2;d++) Bonds[c_ind_n*2+d] = _MyCrystal->getOrientedSystem()->getBonds()[n*2+d] + c_ind_a;
						BondType[c_ind_n] = _MyCrystal->getOrientedSystem()->getBondType(n);
					}
				}
			}
			if( IsAngle ){
				if( !FullGrains && i == 0 && j == 0 ){
					for(unsigned int n=0;n<nbAngles_G1;n++){
						for(unsigned int d=0;d<3;d++){
							Angles_G1[n*3+d] = _MyCrystal->getOrientedSystem()->getAngles()[n*3+d];
							Angles[n*3+d] = Angles_G1[n*3+d];
						}
						AngleType_G1[n] = _MyCrystal->getOrientedSystem()->getAngleType(n);
						AngleType[n] = AngleType_G1[n];
					}
				}else if( FullGrains ){
					unsigned int c_ind = i*dupY1*nbAngles_G1+j*nbAngles_G1;
					unsigned int c_ind_a = i*dupY1*nbAtom1+j*nbAtom1;
					for(unsigned int n=0;n<nbAngles_G1;n++){
						unsigned int c_ind_n = c_ind+n;
						for(unsigned int d=0;d<3;d++){
							Angles_G1[c_ind_n*3+d] = _MyCrystal->getOrientedSystem()->getAngles()[n*3+d] + c_ind_a;
							Angles[c_ind_n*3+d] = Angles_G1[c_ind_n*3+d];
						}
						AngleType_G1[c_ind_n] = _MyCrystal->getOrientedSystem()->getAngleType(n);
						AngleType[c_ind_n] = AngleType_G1[c_ind_n];
					}
				}else{
					unsigned int c_ind = i*dupY1*nbAngles_G1+j*nbAngles_G1;
					unsigned int c_ind_a = i*dupY1*nbAtom1+j*nbAtom1;
					for(unsigned int n=0;n<nbAngles_G1;n++){
						unsigned int c_ind_n = c_ind+n;
						for(unsigned int d=0;d<3;d++) Angles[c_ind_n*3+d] = _MyCrystal->getOrientedSystem()->getAngles()[n*3+d] + c_ind_a;
						AngleType[c_ind_n] = _MyCrystal->getOrientedSystem()->getAngleType(n);
					}
				}
			}
		}
	}
	for(unsigned int i=0;i<dupX2;i++){
		for(unsigned int j=0;j<dupY2;j++){
			for(unsigned int n=0;n<nbAtom2;n++){
				this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n] = this->_MyCrystal2->getOrientedSystem()->getAtom(n);
				this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.x = Mx2*(this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.x + i*xl2);
				this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.y = My2*(this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n].pos.y + j*yl2);
				if( !FullGrains && i == 0 && j == 0 ){
					this->AtomList_G2[n] = this->AtomList[nbAtom1*dupX1*dupY1+n];
					if( IsMolId ){
						MolId_G2[n] = (this->_MyCrystal2->getOrientedSystem()->getMolId(n));
						MolId[nbAtom1*dupX1*dupY1+n] = MolId_G2[n];
					}
				}else if( FullGrains ){
					this->AtomList_G2[i*dupY2*nbAtom2+j*nbAtom2+n] = this->AtomList[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n];
					if( IsMolId ){
						MolId_G2[i*dupY2*nbAtom2+j*nbAtom2+n] = (this->_MyCrystal2->getOrientedSystem()->getMolId(n)) + i*dupY2*nbMol_G2 + j*nbMol_G2;
						MolId[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n] = MolId_G2[i*dupY2*nbAtom2+j*nbAtom2+n] + nbMol_G1*dupX1*dupY1;
					}
				}else if( IsMolId ) MolId[nbAtom1*dupX1*dupY1+i*dupY2*nbAtom2+j*nbAtom2+n] = (this->_MyCrystal2->getOrientedSystem()->getMolId(n)) + i*dupY2*nbMol_G2 + j*nbMol_G2 + nbMol_G1*dupX1*dupY1;
			}
			if( IsBond ){
				if( !FullGrains && i == 0 && j == 0 ){
					for(unsigned int n=0;n<nbBonds_G2;n++){
						for(unsigned int d=0;d<2;d++){
							Bonds_G2[n*2+d] = _MyCrystal2->getOrientedSystem()->getBonds()[n*2+d];
							Bonds[(n+nbBonds_G1*dupX1*dupY1)*2+d] = Bonds_G2[n*2+d] + nbAtom1*dupX1*dupY1;
						}
						BondType_G2[n] = _MyCrystal2->getOrientedSystem()->getBondType(n);
						BondType[n+nbBonds_G1*dupX1*dupY1] = BondType_G2[n];
					}
				}else if( FullGrains ){
					unsigned int c_ind = i*dupY2*nbBonds_G2+j*nbBonds_G2;
					unsigned int c_ind_a = i*dupY2*nbAtom2+j*nbAtom2;
					for(unsigned int n=0;n<nbBonds_G2;n++){
						unsigned int c_ind_n = c_ind+n;
						for(unsigned int d=0;d<2;d++){
							Bonds_G2[c_ind_n*2+d] = _MyCrystal2->getOrientedSystem()->getBonds()[n*2+d] + c_ind_a;
							Bonds[(c_ind_n+nbBonds_G1*dupX1*dupY1)*2+d] = Bonds_G2[c_ind_n*2+d] + nbAtom1*dupX1*dupY1;
						}
						BondType_G2[c_ind_n] = _MyCrystal2->getOrientedSystem()->getBondType(n);
						BondType[c_ind_n + nbBonds_G1*dupX1*dupY1] = BondType_G2[c_ind_n];
					}
				}else{
					unsigned int c_ind = i*dupY2*nbBonds_G2+j*nbBonds_G2 + nbBonds_G1*dupX1*dupY1;
					unsigned int c_ind_a = i*dupY2*nbAtom2+j*nbAtom2 + nbAtom1*dupX1*dupY1;
					for(unsigned int n=0;n<nbBonds_G2;n++){
						unsigned int c_ind_n = c_ind+n;
						for(unsigned int d=0;d<2;d++) Bonds[c_ind_n*2+d] = _MyCrystal2->getOrientedSystem()->getBonds()[n*2+d] + c_ind_a;
						BondType[c_ind_n] = _MyCrystal2->getOrientedSystem()->getBondType(n);
					}
				}
			}
			if( IsAngle ){
				if( !FullGrains && i == 0 && j == 0 ){
					for(unsigned int n=0;n<nbAngles_G2;n++){
						for(unsigned int d=0;d<3;d++){
							Angles_G2[n*3+d] = _MyCrystal2->getOrientedSystem()->getAngles()[n*3+d];
							Angles[(n + nbAngles_G1*dupX1*dupY1)*3+d] = Angles_G2[n*3+d] + nbAtom1*dupX1*dupY1;
						}
						AngleType_G2[n] = _MyCrystal2->getOrientedSystem()->getAngleType(n);
						AngleType[n + nbAngles_G1*dupX1*dupY1] = AngleType_G2[n];
					}
				}else if( FullGrains ){
					unsigned int c_ind = i*dupY2*nbAngles_G2+j*nbAngles_G2;
					unsigned int c_ind_a = i*dupY2*nbAtom2+j*nbAtom2;
					for(unsigned int n=0;n<nbAngles_G2;n++){
						unsigned int c_ind_n = c_ind+n;
						for(unsigned int d=0;d<3;d++){
							Angles_G2[c_ind_n*3+d] = _MyCrystal2->getOrientedSystem()->getAngles()[n*3+d] + c_ind_a;
							Angles[(c_ind_n + nbAngles_G1*dupX1*dupY1)*3+d] = Angles_G2[c_ind_n*3+d] + nbAtom1*dupX1*dupY1;
						}
						AngleType_G2[c_ind_n] = _MyCrystal2->getOrientedSystem()->getAngleType(n);
						AngleType[c_ind_n + nbAngles_G1*dupX1*dupY1] = AngleType_G2[c_ind_n];
					}
				}else{
					unsigned int c_ind = i*dupY2*nbAngles_G2+j*nbAngles_G2 + nbAngles_G1*dupX1*dupY1;
					unsigned int c_ind_a = i*dupY2*nbAtom2+j*nbAtom2 + nbAtom1*dupX1*dupY1;
					for(unsigned int n=0;n<nbAngles_G2;n++){
						unsigned int c_ind_n = c_ind+n;
						for(unsigned int d=0;d<3;d++) Angles[c_ind_n*3+d] = _MyCrystal2->getOrientedSystem()->getAngles()[n*3+d] + c_ind_a;
						AngleType[c_ind_n] = _MyCrystal2->getOrientedSystem()->getAngleType(n);
					}
				}
			}

		}
	}

	if( FullGrains ){
		nbBonds_G1 *= dupX1*dupY1;	
		nbBonds_G2 *= dupX2*dupY2;	
		nbAngles_G1 *= dupX1*dupY1;	
		nbAngles_G2 *= dupX2*dupY2;	
	}

	// construct the two grains
	if( !IsMolId ){
		this->Grain1 = new AtomicSystem(this->AtomList_G1, nbAtom1_G, this->_MyCrystal, this->H1_G1, this->H2_G1, this->H3_G1);
		this->Grain2 = new AtomicSystem(this->AtomList_G2, nbAtom2_G, this->_MyCrystal2, this->H1_G2, this->H2_G2, this->H3_G2);
	}else if( !IsBond ){
		this->Grain1 = new AtomicSystem(this->AtomList_G1, nbAtom1_G, this->_MyCrystal, this->H1_G1, this->H2_G1, this->H3_G1, this->MolId_G1);
		this->Grain2 = new AtomicSystem(this->AtomList_G2, nbAtom2_G, this->_MyCrystal2, this->H1_G2, this->H2_G2, this->H3_G2, this->MolId_G2);
	}else if( !IsAngle ){
		this->Grain1 = new AtomicSystem(this->AtomList_G1, nbAtom1_G, this->_MyCrystal, this->H1_G1, this->H2_G1, this->H3_G1, this->MolId_G1, nbBonds_G1, nbBondType_G1, this->Bonds_G1,this->BondType_G1);
		this->Grain2 = new AtomicSystem(this->AtomList_G2, nbAtom2_G, this->_MyCrystal2, this->H1_G2, this->H2_G2, this->H3_G2, this->MolId_G2, nbBonds_G2, nbBondType_G2, this->Bonds_G2,this->BondType_G2);
	}else{
		this->Grain1 = new AtomicSystem(this->AtomList_G1, nbAtom1_G, this->_MyCrystal, this->H1_G1, this->H2_G1, this->H3_G1, this->MolId_G1, nbBonds_G1, nbBondType_G1, this->Bonds_G1, this->BondType_G1, nbAngles_G1, nbAngleType_G1, this->Angles_G1, this->AngleType_G1);
		this->Grain2 = new AtomicSystem(this->AtomList_G2, nbAtom2_G, this->_MyCrystal2, this->H1_G2, this->H2_G2, this->H3_G2, this->MolId_G2, nbBonds_G2, nbBondType_G2, this->Bonds_G2, this->BondType_G2, nbAngles_G2, nbAngleType_G2, this->Angles_G2, this->AngleType_G2);
	}

	this->AreGrainsDefined = true;
	string h_a_str = to_string(this->h_a);
	string k_a_str = to_string(this->k_a);
	string l_a_str = to_string(this->l_a);
	string theta_str = to_string(this->theta*180/M_PI);
	string h_p_str = to_string(this->h_p);
	string k_p_str = to_string(this->k_p);
	string l_p_str = to_string(this->l_p);
	string h_p_1_z_str = to_string(_MyCrystal->getOrthogonalPlanes()[6]);
	string k_p_1_z_str = to_string(_MyCrystal->getOrthogonalPlanes()[7]);
	string l_p_1_z_str = to_string(_MyCrystal->getOrthogonalPlanes()[8]);
	string h_p_2_z_str = to_string(_MyCrystal2->getOrthogonalPlanes()[6]);
	string k_p_2_z_str = to_string(_MyCrystal2->getOrthogonalPlanes()[7]);
	string l_p_2_z_str = to_string(_MyCrystal2->getOrthogonalPlanes()[8]);
	string h_p_1_x_str = to_string(_MyCrystal->getOrthogonalPlanes()[0]);
	string k_p_1_x_str = to_string(_MyCrystal->getOrthogonalPlanes()[1]);
	string l_p_1_x_str = to_string(_MyCrystal->getOrthogonalPlanes()[2]);
	string h_p_2_x_str = to_string(_MyCrystal2->getOrthogonalPlanes()[0]);
	string k_p_2_x_str = to_string(_MyCrystal2->getOrthogonalPlanes()[1]);
	string l_p_2_x_str = to_string(_MyCrystal2->getOrthogonalPlanes()[2]);
	string h_p_1_y_str = to_string(_MyCrystal->getOrthogonalPlanes()[3]);
	string k_p_1_y_str = to_string(_MyCrystal->getOrthogonalPlanes()[4]);
	string l_p_1_y_str = to_string(_MyCrystal->getOrthogonalPlanes()[5]);
	string h_p_2_y_str = to_string(_MyCrystal2->getOrthogonalPlanes()[3]);
	string k_p_2_y_str = to_string(_MyCrystal2->getOrthogonalPlanes()[4]);
	string l_p_2_y_str = to_string(_MyCrystal2->getOrthogonalPlanes()[5]);
	string h_d_1_z_str = to_string(_MyCrystal->getOrthogonalDirs()[6]);
	string k_d_1_z_str = to_string(_MyCrystal->getOrthogonalDirs()[7]);
	string l_d_1_z_str = to_string(_MyCrystal->getOrthogonalDirs()[8]);
	string h_d_2_z_str = to_string(_MyCrystal2->getOrthogonalDirs()[6]);
	string k_d_2_z_str = to_string(_MyCrystal2->getOrthogonalDirs()[7]);
	string l_d_2_z_str = to_string(_MyCrystal2->getOrthogonalDirs()[8]);
	string h_d_1_x_str = to_string(_MyCrystal->getOrthogonalDirs()[0]);
	string k_d_1_x_str = to_string(_MyCrystal->getOrthogonalDirs()[1]);
	string l_d_1_x_str = to_string(_MyCrystal->getOrthogonalDirs()[2]);
	string h_d_2_x_str = to_string(_MyCrystal2->getOrthogonalDirs()[0]);
	string k_d_2_x_str = to_string(_MyCrystal2->getOrthogonalDirs()[1]);
	string l_d_2_x_str = to_string(_MyCrystal2->getOrthogonalDirs()[2]);
	string h_d_1_y_str = to_string(_MyCrystal->getOrthogonalDirs()[3]);
	string k_d_1_y_str = to_string(_MyCrystal->getOrthogonalDirs()[4]);
	string l_d_1_y_str = to_string(_MyCrystal->getOrthogonalDirs()[5]);
	string h_d_2_y_str = to_string(_MyCrystal2->getOrthogonalDirs()[3]);
	string k_d_2_y_str = to_string(_MyCrystal2->getOrthogonalDirs()[4]);
	string l_d_2_y_str = to_string(_MyCrystal2->getOrthogonalDirs()[5]);
	string i_d_1_x_str(""), i_d_1_y_str(""), i_d_1_z_str(""), i_d_2_x_str(""), i_d_2_y_str(""), i_d_2_z_str("");
	string i_p_1_x_str(""), i_p_1_y_str(""), i_p_1_z_str(""), i_p_2_x_str(""), i_p_2_y_str(""), i_p_2_z_str("");
	if( _MyCrystal->getCrystallo() == "Hexagonal" ){
		i_d_1_x_str = "_"+to_string(-_MyCrystal->getOrthogonalDirs()[0]-_MyCrystal->getOrthogonalDirs()[1]);
		i_d_1_y_str = "_"+to_string(-_MyCrystal->getOrthogonalDirs()[3]-_MyCrystal->getOrthogonalDirs()[4]);
		i_d_1_z_str = "_"+to_string(-_MyCrystal->getOrthogonalDirs()[6]-_MyCrystal->getOrthogonalDirs()[7]);
		i_p_1_x_str = "_"+to_string(-_MyCrystal->getOrthogonalPlanes()[0]-_MyCrystal->getOrthogonalPlanes()[1]);
		i_p_1_y_str = "_"+to_string(-_MyCrystal->getOrthogonalPlanes()[3]-_MyCrystal->getOrthogonalPlanes()[4]);
		i_p_1_z_str = "_"+to_string(-_MyCrystal->getOrthogonalPlanes()[6]-_MyCrystal->getOrthogonalPlanes()[7]);
	}
	if( _MyCrystal2->getCrystallo() == "Hexagonal" ){
		i_d_2_x_str = "_"+to_string(-_MyCrystal2->getOrthogonalDirs()[0]-_MyCrystal2->getOrthogonalDirs()[1]);
		i_d_2_y_str = "_"+to_string(-_MyCrystal2->getOrthogonalDirs()[3]-_MyCrystal2->getOrthogonalDirs()[4]);
		i_d_2_z_str = "_"+to_string(-_MyCrystal2->getOrthogonalDirs()[6]-_MyCrystal2->getOrthogonalDirs()[7]);
		i_p_2_x_str = "_"+to_string(-_MyCrystal2->getOrthogonalPlanes()[0]-_MyCrystal2->getOrthogonalPlanes()[1]);
		i_p_2_y_str = "_"+to_string(-_MyCrystal2->getOrthogonalPlanes()[3]-_MyCrystal2->getOrthogonalPlanes()[4]);
		i_p_2_z_str = "_"+to_string(-_MyCrystal2->getOrthogonalPlanes()[6]-_MyCrystal2->getOrthogonalPlanes()[7]);
	}

	SetSystemSizeAccordingToFixedParam(true,true,false);	

        if( vacuum ){
                double fac_vacuum = 1.5;
                double zero = 0.;
                double shift_z = H3[2] * ((fac_vacuum-1.)/2.);
                H3[2] *= fac_vacuum;
                ApplyShift(zero,zero,shift_z);
        }


	this->File_Heading = " # ["+h_a_str+"_"+k_a_str+i_a_str+"_"+l_a_str+"]"+theta_str+"°("+h_p_str+"_"+k_p_str+i_p_str+"_"+l_p_str+") "+crystalName+" grain boundary\n # The present GB have plane ("+h_p_1_z_str+"_"+k_p_1_z_str+i_p_1_z_str+"_"+l_p_1_z_str+") for upper grain and ("+h_p_2_z_str+"_"+k_p_2_z_str+i_p_2_z_str+"_"+l_p_2_z_str+") for lower grain\n";
	this->Grain1->set_File_Heading(" # Upper grain of the ["+h_a_str+"_"+k_a_str+i_a_str+"_"+l_a_str+"]"+theta_str+"°("+h_p_str+"_"+k_p_str+i_p_str+"_"+l_p_str+") "+crystalName+" grain boundary\n # This system has x <=> ["+h_d_1_x_str+"_"+k_d_1_x_str+i_d_1_x_str+"_"+l_d_1_x_str+"], y <=> ["+h_d_1_y_str+"_"+k_d_1_y_str+i_d_1_y_str+"_"+l_d_1_y_str+"], z <=> ["+h_d_1_z_str+"_"+k_d_1_z_str+i_d_1_z_str+"_"+l_d_1_z_str+"] and x <=> ("+h_p_1_x_str+"_"+k_p_1_x_str+i_p_1_x_str+"_"+l_p_1_x_str+"), y <=> ("+h_p_1_y_str+"_"+k_p_1_y_str+i_p_1_y_str+"_"+l_p_1_y_str+"), z <=> ("+h_p_1_z_str+"_"+k_p_1_z_str+i_p_1_z_str+"_"+l_p_1_z_str+")\n");
	this->Grain2->set_File_Heading(" # Lower grain of the ["+h_a_str+"_"+k_a_str+i_a_str+"_"+l_a_str+"]"+theta_str+"°("+h_p_str+"_"+k_p_str+i_p_str+"_"+l_p_str+") "+crystalName+" grain boundary\n # This system has x <=> ["+h_d_2_x_str+"_"+k_d_2_x_str+i_d_2_x_str+"_"+l_d_2_x_str+"], y <=> ["+h_d_2_y_str+"_"+k_d_2_y_str+i_d_2_y_str+"_"+l_d_2_y_str+"], z <=> ["+h_d_2_z_str+"_"+k_d_2_z_str+i_d_2_z_str+"_"+l_d_2_z_str+"] and x <=> ("+h_p_2_x_str+"_"+k_p_2_x_str+i_p_2_x_str+"_"+l_p_2_x_str+"), y <=> ("+h_p_2_y_str+"_"+k_p_2_y_str+i_p_2_y_str+"_"+l_p_2_y_str+"), z <=> ("+h_p_2_z_str+"_"+k_p_2_z_str+i_p_2_z_str+"_"+l_p_2_z_str+")\n");
	this->xl1 *= Mx1; 
	this->xl2 *= Mx2; 
	this->yl1 *= My1; 
	this->yl2 *= My2; 
}
//
void Bicrystal::PasteGrains(AtomicSystem* _Grain1, AtomicSystem* _Grain2) {
	this->Grain1->SetSystemSizeAccordingToFixedParam(true,true,false);
	this->Grain2->SetSystemSizeAccordingToFixedParam(true,true,false);
	unsigned int nbAtom1 = _Grain1->getNbAtom();
	unsigned int nbAtom2 = _Grain2->getNbAtom();
	unsigned int new_nbAtom = nbAtom1 + nbAtom2;
	if( nbAtom != new_nbAtom ){
		if( IsAtomListMine ) delete[] AtomList;
		this->nbAtom = new_nbAtom;
		this->AtomList = new Atom[this->nbAtom];
	}
	
	for (unsigned int i = 0; i < nbAtom1; ++i) {
		Atom A = _Grain1->getAtom(i);
		A.pos.z += _Grain2->getH3()[2] + (GBspace / 2.0);
		this->AtomList[i] = A;
	}
	for (unsigned int i = 0; i < nbAtom2; ++i) this->AtomList[nbAtom1 + i] = _Grain2->getAtom(i);
}

void Bicrystal::ShiftGrainsAlongUCInPlane(unsigned int n1, unsigned int n2, bool vacuum) {
	double shift_x, shift_y;
	double shift_z = 0.;
	double zero(0.);
	char filename[128];
	std::string File_Heading_save = this->File_Heading; 
	if( vacuum ){
		shift_z = this->H3[2] * ((fac_vacuum-1.)/2.);
		for(unsigned int i=0;i<3;i++) this->H3[i] *= fac_vacuum;
		computeInverseCellVec();

	}
	this->Grain2->ApplyShift(zero, zero, shift_z);
	for (unsigned int i = 0; i < n1; ++i) {
		for (unsigned int j = 0; j < n2; ++j) {
			// Compute the components of the shift vector by dividing i, j, and k by their respective n values
			shift_x = ((double) (i)/n1) * this->H1[0] + ((double) (j)/n2) * this->H2[0];
			shift_y = ((double) (i)/n1) * this->H1[1] + ((double) (j)/n2) * this->H2[1];
			 // Apply the shift
			this->Grain1->ApplyShift(shift_x, shift_y, shift_z);
			this->PasteGrains(this->Grain1, this->Grain2);
			// Output file
			snprintf(filename, sizeof(filename), "GB_Shift_%u_%u.lmp", i, j);
			this->File_Heading = File_Heading_save+" # Applied shift: [" + to_string(shift_x) + ", " + to_string(shift_y) + "]\n";
			this->print_lmp(filename);
			this->Grain1->ApplyShift(-shift_x, -shift_y, -shift_z);
		}
	}
	this->Grain2->ApplyShift(zero, zero, -shift_z);
}


void Bicrystal::ShiftGrainsAlongDSC(unsigned int n1, unsigned int n2, unsigned int n3, bool vacuum) { // TODO for the moment it does not account for the DoNotSep cases
	if( !this->DSC_Basis ){
		cerr << "DSC basis is not defined. Aborting shift." << endl;
		return;
	}
	
	double shift_x, shift_y, shift_z, shift_2rm;
	double shift_z_vac = 0.;
	double zero(0.);
	char filename[128];
	string File_Heading_save = this->File_Heading; 

	double fac_vacuum = 1.5;
	if( vacuum ){
		shift_z_vac = this->H3[2] * ((fac_vacuum-1.)/2.);
		for(unsigned int i=0;i<3;i++) this->H3[i] *= fac_vacuum;
		computeInverseCellVec();
	}
	this->Grain2->ApplyShift(zero, zero, shift_z_vac);

	for(unsigned int i=0;i<n1;i++){
		for(unsigned int j=0;j<n2;j++){
			for(unsigned int k=0;k<n3;k++){
				shift_x = ((double) (i)/n1) * this->DSC_Basis[0] + ((double) (j)/n2) * this->DSC_Basis[1] + ((double) (k)/n3)* this->DSC_Basis[2];
				shift_y = ((double) (i)/n1) * this->DSC_Basis[3] + ((double) (j)/n2)  * this->DSC_Basis[4] +  ((double) (k)/n3)* this->DSC_Basis[5];
				shift_z = ((double) (i)/n1)* this->DSC_Basis[6] + ((double) (j)/n2) * this->DSC_Basis[7] +  ((double) (k)/n3) * this->DSC_Basis[8];
				
				this->Grain1->ApplyShift(shift_x, shift_y, shift_z);
				this->Grain1->computeWrap();
				
				// Replace positions with the wrapped positions
				unsigned int nbAtom = this->Grain1->getNbAtom();
				for(unsigned int n=0;n<nbAtom;n++){
					Position P = this->Grain1->getWrappedPosition(n);
					this->Grain1->getAtomRef(n).pos = P;
				}
				 
				// the shift for vacuum
				shift_2rm = -shift_z;
				if( vacuum ){
					this->Grain1->ApplyShift(zero, zero, shift_z_vac);
					shift_2rm -= shift_z_vac;
				}
				this->PasteGrains(this->Grain1, this->Grain2);

				// Output file
				snprintf(filename, sizeof(filename), "GB_DSC_Shift_%u_%u_%u.lmp", i, j, k);
				this->File_Heading = File_Heading_save+"# Applied DSC shift: [" + to_string(shift_x) + ", " + to_string(shift_y) + ", " + to_string(shift_z) + "]\n";
				this->print_lmp(filename);
				this->Grain1->ApplyShift(-shift_x, -shift_y, shift_2rm);
			}
		}
	}
}

void Bicrystal::ShiftGrainsAlongCSL(unsigned int n1, unsigned int n2, unsigned int n3, bool vacuum) {// TODO for the moment it does not account for the DoNotSep cases
	if( !this->CSL_Basis ){
		cerr << "CSL basis is not defined. Aborting shift." << endl;
		return;
	}
	
	double shift_x, shift_y, shift_z, shift_2rm;
	double shift_z_vac = 0.;
	double zero(0.);
	char filename[128];
	string File_Heading_save = this->File_Heading;
	
	if( vacuum ){
		shift_z_vac = this->H3[2] * ((fac_vacuum-1.)/2.);
		for(unsigned int i=0;i<3;i++) this->H3[i] *= fac_vacuum;
		computeInverseCellVec();
	}
	this->Grain2->ApplyShift(zero, zero, shift_z_vac);
	
	for(unsigned int i=0;i<n1;i++){
		for(unsigned int j=0;j<n2;j++){
			for(unsigned int k=0;k<n3;k++){
				// Compute the translation vector based on the CSL basis
				shift_x = ((double) (i) / n1) * CSL_Basis[0] + ((double) (j) / n2) * CSL_Basis[1] + ((double) (k) / n3) * CSL_Basis[2];
				shift_y = ((double) (i) / n1) * CSL_Basis[3] + ((double) (j) / n2) * CSL_Basis[4] + ((double) (k) / n3) * CSL_Basis[5];
				shift_z = ((double) (i) / n1) * CSL_Basis[6] + ((double) (j) / n2) * CSL_Basis[7] + ((double) (k) / n3) * CSL_Basis[8];
				
				this->Grain1->ApplyShift(shift_x, shift_y, shift_z);
				this->Grain1->computeWrap();
				
				// Replace positions with the wrapped positions
				unsigned int nbAtom = this->Grain1->getNbAtom();
				for (unsigned int n = 0; n < nbAtom; ++n) {
					Position P = this->Grain1->getWrappedPosition(n);
					this->Grain1->getAtomRef(n).pos = P;
				}
	
				// the shift for vacuum
				shift_2rm = -shift_z;
				if( vacuum ){
					this->Grain1->ApplyShift(zero, zero, shift_z_vac);
					shift_2rm -= shift_z_vac;
				}
				this->PasteGrains(this->Grain1, this->Grain2);
				
				// Save the output file
				snprintf(filename, sizeof(filename), "GB_CSL_Shift_%u_%u_%u.lmp", i, j, k);
				this->File_Heading = File_Heading_save+"# Applied CSL shift: [" + to_string(shift_x) + ", " + to_string(shift_y) + ", " + to_string(shift_z) + "]\n";
				this->print_lmp(filename);
				this->Grain1->ApplyShift(-shift_x, -shift_y, shift_2rm);
			}
		}
	}
}


Bicrystal::Bicrystal(const string& filename, const string NormalDir, const string CrystalName):AtomicSystem(filename){
	InitCSL_DSC();
	read_params();
	if( NormalDir != "z" && NormalDir != "y" && NormalDir != "x" && NormalDir != "Z" && NormalDir != "Y" && NormalDir != "X" ){
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

Bicrystal::Bicrystal(const string& filename, const string CrystalName):AtomicSystem(filename){
	InitCSL_DSC();
	read_params();
	setCrystal(CrystalName);
	this->CA = new ComputeAuxiliary(this);
	this->IsCA = true;
}

Bicrystal::Bicrystal(const string& filename):AtomicSystem(filename){
	InitCSL_DSC();
	read_params();
	this->CA = new ComputeAuxiliary(this);
	this->IsCA = true;
}

void Bicrystal::InitCSL_DSC(){
	if( !CSL_Basis ){
		CSL_Basis = new double[9];
		for(unsigned int i=0;i<9;i++) CSL_Basis[i] = 0.;
	}
	if( !DSC_Basis ){
		DSC_Basis = new double[9];
		for(unsigned int i=0;i<9;i++) DSC_Basis[i] = 0.;
	}
}

void Bicrystal::setOrientedCrystals(const string& crystalName, bool rationalize, vector<string> Properties, int h_px, int k_px, int l_px){
	if( !this->RotCartToGrain2 ) this->RotCartToGrain2 = new double[9];
	for(unsigned int i=0;i<9;i++) RotCartToGrain2[i] = 0.;
	// construct the first crystal with the wanted plane horyzontal
	this->_MyCrystal->RotateCrystal(h_p,k_p,l_p,h_px,k_px,l_px);
	// shift it relative the lower grain shift
	this->_MyCrystal->ShiftMotif(lowgrain_shift_x,lowgrain_shift_y,lowgrain_shift_z);

	// compute the rotation to be applied to the second grain
	if( rationalize ){
		double *rot_ax = new double[3];
		int *CSL_vec = new int[3];
		int h_p,k_p,l_p;
		double temp_var = BT->RationalizeOri(_MyCrystal,this->h_a,this->k_a,this->l_a,this->theta,h_p,k_p,l_p,CSL_vec); // TODO
		this->theta = temp_var;
		this->_MyCrystal->ComputeRotMat(h_p,k_p,l_p,theta,RotCartToGrain2,false);
		delete[] rot_ax;
		delete[] CSL_vec;
	}else this->_MyCrystal->ComputeRotMat(h_a,k_a,l_a,theta,RotCartToGrain2);
	
	this->_MyCrystal2 = new Crystal(crystalName);
	this->_MyCrystal2->ReadProperties(Properties);
	this->_MyCrystal2->RotateCrystal(this->RotCartToGrain2);
	this->_MyCrystal2->ShiftMotif(upgrain_shift_x,upgrain_shift_y,upgrain_shift_z);
	
	BT->ComputeCSL_DSC(_MyCrystal,_MyCrystal2,CSL_Basis,DSC_Basis,0);
	BT->generateCSL();
	this->_MyCrystal->ConstructOrthogonalCell();
	this->_MyCrystal->getOrientedSystem()->getMinBoxSizes()[2] = this->MinBoxSizes[2]/2.;
	this->_MyCrystal->getOrientedSystem()->SetSystemSizeAccordingToFixedParam(false,false,true);	
	cout << "First grain constructed (" << this->_MyCrystal->getOrientedSystem()->getNbAtom() << " atoms)" << endl;
	this->_MyCrystal2->ConstructOrthogonalCell();
	this->_MyCrystal2->getOrientedSystem()->getMinBoxSizes()[2] = this->MinBoxSizes[2]/2.;
	this->_MyCrystal2->getOrientedSystem()->SetSystemSizeAccordingToFixedParam(false,false,true);	
	cout << "Second grain constructed (" << this->_MyCrystal2->getOrientedSystem()->getNbAtom() << " atoms)" << endl;

	cout << "Building orthogonal cells required deforming the crystals such that:" << endl;
	vector<vector<string>> arr_element = {
		{"","X","Y","Z","Shear"},
		{"Grain 1",to_string(_MyCrystal->GetCrystalDef()[0])+" %",to_string(_MyCrystal->GetCrystalDef()[1])+" %",to_string(_MyCrystal->GetCrystalDef()[2])+" %",to_string(_MyCrystal->GetCrystalDef()[3])+" %"},
		{"Grain 2",to_string(_MyCrystal2->GetCrystalDef()[0])+" %",to_string(_MyCrystal2->GetCrystalDef()[1])+" %",to_string(_MyCrystal2->GetCrystalDef()[2])+" %",to_string(_MyCrystal2->GetCrystalDef()[3])+" %"}
	};
	vector<vector<unsigned int>> arr_fusion = {{1,1,1,1,1},{1,1,1,1,1},{1,1,1,1,1}};
	Dis->DisplayArray(arr_element,arr_fusion);
	cout << "If these values are too large it may cause issues during the relaxation of the GB" << endl;
	cout << "These values can be reduced by decreasing TOL_ORTHO_BOX and TOL_ORTHO_BOX_Z in a FixedParameters.ath file" << endl;
	cout << endl;
	
	this->IsCrystal2 = true;
	this->IsRotMatDefine = true;
}

void Bicrystal::searchGBPos(){
	// Compute bond orientationnal order parameter to know where are the GBs and stored it into aux 1
	unsigned int nbPts_i = 1000;// TODO maybe set a file reading sigma nbPts, lsph and rc for the different systems
	double sigma = 2.;
	setAux(this->CA->BondOrientationalParameter(),"Disorder");
	// Compute density profile of bond ori param
	unsigned int ind_DisoDens;
	ind_DisoDens = Compute1dDensity("Disorder", NormalDir, sigma, nbPts_i);
	// Compute atomic density profile to know if there is vacuum
	unsigned int ind_AtDens;
	ind_AtDens = Compute1dDensity("Atomic", NormalDir, sigma, nbPts_i);
	double VacuumDens = 0.;
	double SecureFac = 1e-2;
	unsigned int count_dens_ave = 0;
	for(unsigned int i=0;i<this->density_nbPts[ind_AtDens];i++){
		if( this->density_prof[ind_AtDens][i*2] > 1e-1 ){
			VacuumDens += this->density_prof[ind_AtDens][i*2];
			count_dens_ave += 1;
		}
	}
	VacuumDens /= count_dens_ave;
	VacuumDens *= 1e-1;

	
	if( NormalDir == "x" ) this->Ldir = this->H1[0];
	else if( NormalDir == "y" ) this->Ldir = this->H2[1];
	else if( NormalDir == "z" ) this->Ldir = this->H3[2];
	
	this->VacuumLo = this->Ldir;
	this->MinPos = this->Ldir;
	this->SystemLength = this->Ldir;
	this->VacuumHi = 0.;
	this->MaxPos = 0.;
	this->IsVacuum = false;
	this->IsCentered = true;
	double buffer;
	vector<double> density_red;
	vector<double> density_red_forfit;
	for(unsigned int i=0;i<this->density_nbPts[ind_AtDens];i++){
		if( this->density_prof[ind_AtDens][i*2] < VacuumDens ){
			this->IsVacuum = true;
			if( this->density_prof[ind_AtDens][i*2+1] > this->VacuumHi ) this->VacuumHi = this->density_prof[ind_AtDens][i*2+1];
			if( this->density_prof[ind_AtDens][i*2+1] < this->VacuumLo ) this->VacuumLo = this->density_prof[ind_AtDens][i*2+1];
		}
	} // TODO non vacuum case ?
	if( this->IsVacuum ){ // set to zero the diso dens near the vacuum to not impact the gaussian fit
		double slab_width = 30.;
		//TODO correct because it is not good for all cases
		for(unsigned int i=0;i<this->density_nbPts[ind_DisoDens];i++){
			if( this->density_prof[ind_DisoDens][i*2+1] < slab_width ) this->density_prof[ind_DisoDens][i*2] = 0.;
			if( this->density_prof[ind_DisoDens][i*2+1] > VacuumLo-slab_width ) this->density_prof[ind_DisoDens][i*2] = 0.;
		}
		//Print1dDensity("dens","Disorder");
	}
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
		if( fabs(this->MaxPos-this->MinPos-this->Ldir) < 1. ){
			this->MinPos = this->VacuumHi-this->Ldir;
			this->MaxPos = this->VacuumLo;
			this->SystemLength = this->MaxPos-this->MinPos;
			this->IsCentered = false;
		}
		cout << "CENTER = " << IsCentered << ", VAC " << IsVacuum << endl;
		cout << "MaxPos = " << MaxPos << ", MinPos = " << MinPos << endl;
		cout << "VacHi = " << VacuumHi << ", VacLo = " << VacuumLo << endl;
		// Search GB position by computing the mean of gaussian distrib in the center of the system
		// compute the max and the mean of the disorder density 
		unsigned int indMaxDiso = 0;
		double MeanDiso = 0.;
		double FracLength = 0.01;
		for(unsigned int i=0;i<this->density_nbPts[ind_DisoDens];i++){
			buffer = this->density_prof[ind_DisoDens][i*2+1];
			if( this->IsVacuum && !this->IsCentered && buffer > this->VacuumHi ) buffer -= this->Ldir;
			// account only for position 60% around the center of system
			//if( buffer > this->MinPos+FracLength*(this->SystemLength) && buffer < this->MaxPos-FracLength*(this->SystemLength) ){
				density_red.push_back(this->density_prof[ind_DisoDens][i*2]);
				density_red.push_back(buffer);
				MeanDiso += this->density_prof[ind_DisoDens][i*2];
				if( density_red[((density_red.size()/2)-1)*2] > density_red[indMaxDiso*2] ) indMaxDiso = (density_red.size()/2)-1;
			//}
		}
		MeanDiso /= (density_red.size()/2.);
		// search the minimal values around the max to define bounds for characterizing the GB
		double rightMin = 0;
		unsigned int indRight = 0;
		unsigned int count_equal;
		unsigned int max_count_equal = 20;
		count_equal = 0;
		for(unsigned int i=indMaxDiso;i<(density_red.size()/2)-1;i++){
			if( density_red[i*2] < density_red[(i+1)*2] ){
				rightMin = density_red[i*2];
				indRight = i;
				break;
			}else if( fabs(density_red[i*2]-density_red[(i+1)*2]) < 1e-8 ) count_equal++;
			if( count_equal >= max_count_equal ){
					rightMin = density_red[i*2];
					indRight = i;
					break;
			}
		}
		double leftMin = 0;
		unsigned int indLeft = 0;
		count_equal = 0;
		for(unsigned int i=indMaxDiso;i>0;i--){
			if( density_red[i*2] < density_red[(i-1)*2] ){
				leftMin = density_red[i*2];
				indLeft = i;
				break;
			}else if( fabs(density_red[i*2]-density_red[(i-1)*2]) < 1e-8 ) count_equal++;
			if( count_equal >= max_count_equal ){
					leftMin = density_red[i*2];
					indLeft = i;
					break;
			}
		}
		// if one of the two is higher than the mean we are in local minimum inside the GB
		// => search the following minimum
		double facMean = 1.5;
		unsigned int max_min_search = 250;
		unsigned int count_min_search = 0;
		bool maxfound;
		cout << "Mean diso : " << MeanDiso << endl;
		while( rightMin > MeanDiso*facMean && count_min_search < max_min_search ){
			maxfound = false;
			count_equal = 0;
			for(unsigned int i=indRight;i<(density_red.size()/2)-1;i++){
				if( !maxfound && density_red[i*2] > density_red[(i+1)*2] ) maxfound = true;
				if( maxfound && density_red[i*2] < density_red[(i+1)*2] ){
					rightMin = density_red[i*2];
					indRight = i;
					break;
				}else if( fabs(density_red[i*2]-density_red[(i+1)*2]) < 1e-8 ) count_equal++;
				if( count_equal >= max_count_equal ){
					rightMin = density_red[i*2];
					indRight = i;
					break;
				}
			}
			count_min_search++;
		}
		cout << "right min : " << rightMin << " " << count_min_search << endl;
		count_min_search = 0;
		while( leftMin > MeanDiso*facMean && count_min_search < max_min_search ){
			maxfound = false;
			count_equal = 0;
			for(unsigned int i=indLeft;i>0;i--){
				if( !maxfound && density_red[i*2] > density_red[(i-1)*2] ) maxfound = true;
				if( maxfound && density_red[i*2] < density_red[(i-1)*2] ){
					leftMin = density_red[i*2];
					indLeft = i;
					break;
				}else if( fabs(density_red[i*2]-density_red[(i-1)*2]) < 1e-8 ) count_equal++;
				if( count_equal >= max_count_equal ){
					leftMin = density_red[i*2];
					indLeft = i;
					break;
				}
			}
			count_min_search++;
		}		
		cout << "left min : " << leftMin << " " << count_min_search << endl;
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
	//} // end if IsVacuum
	double sigma_fit, mu_fit, prefac;
	sigma_fit = this->GBwidth1;
	mu_fit = this->GBPos1;
	prefac = 1.;
	// fit a gaussian to rafine the estimation of GB width and position
	MT->gaussian_fit(density_red_forfit, mu_fit, sigma_fit, prefac);
	this->GBwidth1 = 2.355*sigma_fit;
	this->GBPos1 = mu_fit;
	// set the gaussian GB profile of AtomicSystem class (without prefac for the distribution to be normalized and for different system to be compared to each others)
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
	index = this->Compute1dDensity("Atomic", this->NormalDir, sigma, nbPts);
	double PC_dens = 0.;
	unsigned int count = 0;
	double buffer;
	this->ExcessVol = 0.;
	// compute density of perfect crystal assuming that there is perfect crystal between MinPos+GBwidth and GBPos-GBwidth and same above GBPos
	if( IsVacuum ){
		vector<double> dens_GB;
		for(unsigned int i=0;i<this->density_nbPts[index];i++){
			buffer = this->density_prof[index][i*2+1];
			if( !this->IsCentered && buffer > this->VacuumHi ) buffer -= this->Ldir;
			if( ( buffer > MinPos+GBwidth1 && buffer < GBPos1-GBwidth1 ) || ( buffer < MaxPos-GBwidth1 && buffer > GBPos1+GBwidth1 ) ){
				PC_dens += this->density_prof[index][i*2];
				count += 1;
			}else if( buffer > GBPos1-GBwidth1 && buffer < GBPos1+GBwidth1 ){
				dens_GB.push_back(density_prof[index][i*2]);
				dens_GB.push_back(buffer);
			}
		}
		PC_dens /= count;
		// compute excess volume by integrating density in the GB region
		for(unsigned int i=0;i<(dens_GB.size()/2)-2;i++) this->ExcessVol += ((2.*PC_dens/(dens_GB[i*2]+dens_GB[(i+1)*2]))-1.)*fabs((dens_GB[(i+1)*2+1]-dens_GB[i*2+1]));

	}
}

void Bicrystal::print_Grains(bool vacuum, string fnameG1, string fnameG2){
	if( this->AreGrainsDefined ){
		this->Grain1->print_lmp(fnameG1);
		this->Grain2->print_lmp(fnameG2);
		if( vacuum ){
			double sz1(0.), sz2(0.), zero(0.);
			sz1 = this->Grain1->getH3()[2] * ((fac_vacuum-1.)/2);
			sz2 = this->Grain2->getH3()[2] * ((fac_vacuum-1.)/2);
			for(unsigned int i=0;i<3;i++){
				this->Grain1->getH3()[i] *= fac_vacuum;
				this->Grain2->getH3()[i] *= fac_vacuum;
			}
			this->Grain1->computeInverseCellVec();
			this->Grain1->ApplyShift(zero,zero,sz1);
			this->Grain2->computeInverseCellVec();
			this->Grain2->ApplyShift(zero,zero,sz2);
			string ext1=fnameG1.substr(fnameG1.find_last_of(".") + 1);
			string ext2=fnameG2.substr(fnameG2.find_last_of(".") + 1);
			string pref1=fnameG1.substr(0,fnameG1.size()-ext1.size()-1);
			string pref2=fnameG2.substr(0,fnameG2.size()-ext2.size()-1);
			this->Grain1->print_lmp(pref1+"_vacuum."+ext1);
			this->Grain2->print_lmp(pref2+"_vacuum."+ext2);
			for(unsigned int i=0;i<3;i++){
				this->Grain1->getH3()[i] /= fac_vacuum;
				this->Grain2->getH3()[i] /= fac_vacuum;
			}
			this->Grain1->ApplyShift(zero,zero,-sz1);
			this->Grain2->ApplyShift(zero,zero,-sz2);
		}
	}else{
		cout << "The two grains are not defined, we cannot print them" << endl;
	}
}

void Bicrystal::read_params(){
	AtomicSystem::read_params();
	//string fp;
	//#ifdef FIXEDPARAMETERS
	//fp = FIXEDPARAMETERS;
	//#endif
	//string backslash="/";
	//string filename=fp+backslash+FixedParam_Filename;
	//ifstream file(filename, ios::in);
	ifstream file(FixedParam_Filename, ios::in);
	size_t pos_MaxMisfit, pos_GBSpace, pos_MaxDup, pos_FullGrains, pos_maxvarfl, pos_maxdupf, pos_shift_grain, pos_shift_upgrain;
	string buffer_s, line;
	if(file){
		while(file){
			getline(file,line);
			pos_maxvarfl=line.find("MAX_VAR_FACET_LENGTH ");
			if(pos_maxvarfl!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->maxVarL;
			}
			pos_maxdupf=line.find("MAX_DUP_FACET ");
			if(pos_maxdupf!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->dupMaxFac;
			}
			pos_shift_grain=line.find("LOWER_GRAIN_SHIFT ");
			if(pos_shift_grain!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->lowgrain_shift_x >> this->lowgrain_shift_y >> this->lowgrain_shift_z;
			}
			pos_shift_upgrain=line.find("UPPER_GRAIN_SHIFT ");
			if(pos_shift_upgrain!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->upgrain_shift_x >> this->upgrain_shift_y >> this->upgrain_shift_z;
			}
			pos_GBSpace=line.find("GB_SPACE ");
			if(pos_GBSpace!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->GBspace;
			}
			pos_MaxMisfit=line.find("MAX_MISFIT ");
			if(pos_MaxMisfit!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->MaxMisfit;
			}
			pos_MaxDup=line.find("MAX_DUP ");
			if(pos_MaxDup!=string::npos){
				istringstream text(line);
				text >> buffer_s >> this->MaxDup;
			}
			pos_FullGrains=line.find("FULL_GRAINS ");
			if(pos_FullGrains!=string::npos){
				istringstream text(line);
				unsigned int buf;
				text >> buffer_s >> buf;
				if( buf == 0 ) FullGrains = false;
				else FullGrains = true;
			}
		}
		file.close();
	}//else{
	//	cerr << "Can't read /data/FixedParameters/Fixed_Parameters.dat file !" << endl;
	//	exit(EXIT_FAILURE);
	//}
}

void Bicrystal::ReadProperties(vector<string> Properties){
	AtomicSystem::ReadProperties(Properties);
	size_t pos_MaxMisfit, pos_GBSpace, pos_MaxDup, pos_FullGrains, pos_maxvarfl, pos_maxdupf, pos_shift_grain, pos_shift_upgrain;
	string buffer_s, line;
	for(unsigned int i=0;i<Properties.size();i++){
		pos_maxvarfl=Properties[i].find("MAX_VAR_FACET_LENGTH ");
		if(pos_maxvarfl!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->maxVarL;
		}
		pos_maxdupf=Properties[i].find("MAX_DUP_FACET ");
		if(pos_maxdupf!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->dupMaxFac;
		}
		pos_shift_grain=Properties[i].find("LOWER_GRAIN_SHIFT ");
		if(pos_shift_grain!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->lowgrain_shift_x >> this->lowgrain_shift_y >> this->lowgrain_shift_z;
		}
		pos_shift_upgrain=Properties[i].find("UPPER_GRAIN_SHIFT ");
		if(pos_shift_upgrain!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->upgrain_shift_x >> this->upgrain_shift_y >> this->upgrain_shift_z;
		}
		pos_GBSpace=Properties[i].find("GB_SPACE ");
		if(pos_GBSpace!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->GBspace;
		}
		pos_MaxMisfit=Properties[i].find("MAX_MISFIT ");
		if(pos_MaxMisfit!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->MaxMisfit;
		}
		pos_MaxDup=Properties[i].find("MAX_DUP ");
		if(pos_MaxDup!=string::npos){
			istringstream text(Properties[i]);
			text >> buffer_s >> this->MaxDup;
		}
		pos_FullGrains=Properties[i].find("FULL_GRAINS ");
		if(pos_FullGrains!=string::npos){
			istringstream text(Properties[i]);
			unsigned int buf;
			text >> buffer_s >> buf;
			if( buf == 0 ) FullGrains = false;
			else FullGrains = true;
		}
	}
}

//void TestFuncAdjustStoich(){ // this function was initially in the Facet GB constructor (to see if it can be useful at some point)
//	// verify the stoichiometry
//	unsigned int *currentStoich = new unsigned int[this->_MyCrystal->getNbAtomType()];
//	for(unsigned int i=0;i<this->_MyCrystal->getNbAtomType();i++) currentStoich[i] = 0;
//	for(unsigned int i=0;i<at_count;i++){
//		for(unsigned int t=0;t<this->_MyCrystal->getNbAtomType();t++){
//			if( AtomList_temp[i].type_uint == t+1 ){
//				currentStoich[t] += 1;
//				break;
//			}
//		}
//	}
//	bool stoich = true;
//	for(unsigned int i=0;i<this->_MyCrystal->getNbAtomType();i++){
//		if( fabs(((double) currentStoich[i]/at_count) - ((double) this->_MyCrystal->getStoich()[i]/this->_MyCrystal->getNbAtom()) ) > 1e-9 ){
//			stoich = false;
//			cout << "The stoichiometry is not the same than parent crystal, ";
//			for(unsigned int t=0;t<this->_MyCrystal->getNbAtomType();t++){
//				cout << "number of " << this->_MyCrystal->getAtomType(t+1) << " : " << currentStoich[t] << ", ";
//			}
//			cout << endl;
//			break;
//		}
//	}
//	stoich = true;
//		cout << "Adjusting stoichiometry.. (this may take a while)" << endl;
//		int *adjust = new int[this->_MyCrystal->getNbAtomType()];
//		int *adjust_tmp = new int[this->_MyCrystal->getNbAtomType()];
//		for(unsigned int i=0;i<this->_MyCrystal->getNbAtomType();i++){
//			adjust[i] = this->_MyCrystal->getNbAtom()*round((double) at_count/this->_MyCrystal->getNbAtom())*((double) this->_MyCrystal->getStoich()[i]/this->_MyCrystal->getNbAtom())-currentStoich[i];
//			adjust_tmp[i] = adjust[i];
//		}
//		// modify to have only negative value (we only want to remove atoms)
//		for(unsigned int i=0;i<this->_MyCrystal->getNbAtomType();i++){
//			if( adjust[i] > 0 ){
//				for(unsigned int t=0;t<this->_MyCrystal->getNbAtomType();t++){
//					adjust[t] -= adjust_tmp[i]*((double) this->_MyCrystal->getStoich()[t]/this->_MyCrystal->getStoich()[i]);
//				}
//				adjust_tmp[i] = adjust[i];
//			}
//		}
//		// neighboring search in the GB zone for ions type to remove
//		double GBup = this->_MyCrystal->getOrientedSystem()->getH3()[2]+3.;
//		double GBdown = this->_MyCrystal->getOrientedSystem()->getH3()[2]-DeltaZ-3.;
//		vector<vector<double>> Neigh;
//		vector<int> Type2Rm;
//		for(unsigned int i=0;i<this->_MyCrystal->getNbAtomType();i++){
//			if( adjust[i] < 0 ){
//				Type2Rm.push_back(i+1);
//				Neigh.push_back(vector<double>());
//			}
//		}
//		double rc_squared = pow(3.,2.);
//		double d_squared;
//		for(unsigned int i=0;i<at_count;i++){
//			for(unsigned int t=0;t<Type2Rm.size();t++){
//				if( AtomList_temp[i].pos.z < GBup && AtomList_temp[i].pos.z > GBdown && AtomList_temp[i].type_uint == Type2Rm[t] ){
//					xpos = AtomList_temp[i].pos.x;
//					ypos = AtomList_temp[i].pos.y;
//					zpos = AtomList_temp[i].pos.z;
//					for(unsigned int j=0;j<at_count;j++){
//						if( i != j && AtomList_temp[j].pos.z < GBup && AtomList_temp[j].pos.z > GBdown  && AtomList_temp[j].type_uint == Type2Rm[t] ){
//							d_squared = pow(xpos-AtomList_temp[j].pos.x,2.) + pow(ypos-AtomList_temp[j].pos.y,2.) + pow(zpos-AtomList_temp[j].pos.z,2.); 
//							if( d_squared < rc_squared ){
//								Neigh[t].push_back(d_squared);
//								Neigh[t].push_back(j);
//							}
//						}
//					}
//				}
//			}
//		}
//		this->nbAtom = at_count;
//		for(unsigned int t=0;t<Type2Rm.size();t++){
//			this->MT->sort(Neigh[t],0,2,Neigh[t]);
//			this->nbAtom += adjust[Type2Rm[t]-1];
//		}
//		this->AtomList = new Atom[this->nbAtom];
//		unsigned int at_count2 = 0;
//		bool store;
//		for(unsigned int i=0;i<at_count;i++){
//			store = true;
//			for(unsigned int t=0;t<Type2Rm.size();t++){
//				if( AtomList_temp[i].type_uint == Type2Rm[t] ){
//					for(unsigned int n=0;n<abs(adjust[Type2Rm[t]-1]);n++){
//						if( round(Neigh[t][n*2+1]) == i ){
//							store = false;
//							if( TagGrain[i] == 1 ) trueNbAt1 -= 1;
//							else if( TagGrain[i] == 2 ) trueNbAt2 -= 1;
//							break;
//						}
//					}
//				}
//			}
//			if( store ){
//				this->AtomList[at_count2] = AtomList_temp[i];
//				if( TagGrain[i] == 2 ) this->AtomList[at_count2].pos.z += (GBspace/2.);
//				at_count2 += 1;
//			}
//		}
//		delete[] adjust;
//		delete[] adjust_tmp;
//}

Bicrystal::~Bicrystal(){
	delete BT;
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
		if( IsMolId ){
			delete[] MolId_G1;
			delete[] MolId_G2;
		}
		if( IsBond ){
			delete[] Bonds_G1;
			delete[] Bonds_G2;
			delete[] BondType_G1;
			delete[] BondType_G2;
		}
		if( IsAngle ){
			delete[] Angles_G1;
			delete[] Angles_G2;
			delete[] AngleType_G1;
			delete[] AngleType_G2;
		}
	}
	if( IsRotMatDefine ){
		delete[] RotCartToGrain2;
	}
	if( CSL_Basis ) delete[] CSL_Basis;
	//
	if( DSC_Basis ) delete[] this->DSC_Basis;
	//
}
