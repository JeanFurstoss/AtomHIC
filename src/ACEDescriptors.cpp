//**********************************************************************************
//*   ACEDescriptors.cpp                                                    *
//**********************************************************************************
//* This file contains the implementation of the ACEDescriptors class       *
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


#include "ACEDescriptors.h"
#include "AtomHicConfig.h"
#include <cmath>
#include <complex>
#include <omp.h>
#include <iomanip>
#include "ML-PACE/ace/ace_b_basis.h"
#include "ML-PACE/ace-evaluator/ace_radial.cpp"

using namespace std; 

ACEDescriptors::ACEDescriptors(AtomicSystem *_MySystem, string _yaml_filename):Descriptors(_MySystem), yaml_filename(_yaml_filename){
	this->name = "ACE";
	ComputeDescriptors();
}

void ACEDescriptors::ComputeDescriptors(){
	BBasisConfiguration MyBBasisConf;
	MyBBasisConf.load(yaml_filename,true);
	MyACEBBase = new ACEBBasisSet(MyBBasisConf);

	// Init pace arrays	
	A.init(MyACEBBase->nelements, MyACEBBase->nradmax + 1, MyACEBBase->lmax + 1, "A");
	A_rank1.init(MyACEBBase->nelements, MyACEBBase->nradbase, "A_rank1");
	
	unsigned int nbAtomType = _MySystem->getNbAtomType();
	unsigned int *elem2pace = new unsigned int[nbAtomType];
	for(unsigned int i=0;i<nbAtomType;i++){
		bool found = false;
		for(unsigned int p=0;p<MyACEBBase->nelements;p++){
			if( MyACEBBase->elements_name[p] == _MySystem->getAtomType(i) ){
				elem2pace[i] = p;
				found = true;
				break;
			}
		}
		if( !found ){
			cerr << "The element \"" << _MySystem->getAtomType(i) << "\" of the provided atomic systems was not found either in the ACE yaml file" << endl;
			exit(EXIT_FAILURE);
		}
	}

	// Compute the maximum dimension (among species) for allocation of the Descriptors array
	unsigned int max_dim = 0;
	// Search neighbours using the ceil list algo of AtomticSystem for the highest cutoff of the ACE specification (then is cutoff is different we will research the true neighbours in this function) maybe not a good idea if one of the cutoff is very large because it will uses to much memory
	// => maybe define a new function in AtomicSystem for neighbour research restricted to two different elements..
	// for the moment the first idea is implemented
	// search maximum cutoff
	double max_cutoff = 0.;
	const unsigned int lmaxi = MyACEBBase->lmax;
	const unsigned int nradiali = MyACEBBase->nradmax;
	const unsigned int nradbase = MyACEBBase->nradbase;
	for(unsigned int i=0;i<nbAtomType;i++){
		unsigned int mu_i = elem2pace[i];
		unsigned int current_dim = MyACEBBase->total_basis_size_rank1[mu_i];
    		const SHORT_INT_TYPE total_basis_size = MyACEBBase->total_basis_size[mu_i];
		cout << "For " << _MySystem->getAtomType(i) << " rank1 size = " << current_dim << ", rankN size = " << total_basis_size;
		auto basis = MyACEBBase->basis[mu_i];
		for (unsigned int func_ind = 0; func_ind < total_basis_size; ++func_ind) {
        		auto func = &basis[func_ind];
			current_dim += func->num_ms_combs;
		}
		if( max_dim < current_dim )
			max_dim = current_dim;
		cout << ", leading to " << current_dim << "-dimension ACE descriptors";
		for(unsigned int j=i;j<nbAtomType;j++){
			unsigned int mu_j = elem2pace[j];
			cout << ", cutoff with " << _MySystem->getAtomType(j) << " = " << MyACEBBase->radial_functions->cut(mu_i,mu_j);
			if( max_cutoff < MyACEBBase->radial_functions->cut(mu_i,mu_j) )
				max_cutoff = MyACEBBase->radial_functions->cut(mu_i,mu_j);
		}
		cout << endl;
	}

	// Allocate descriptor array and initialize to 0
	dim = max_dim;
	_Descriptors = new double[nbDatTot*dim];
	for(unsigned int i=0;i<nbDatTot;i++)
		for(unsigned int d=0;d<dim;d++) _Descriptors[i*dim+d] = 0.;

	// search neighbours
	if( !_MySystem->getIsNeighbours() || _MySystem->get_current_rc() != max_cutoff ){
		_MySystem->searchNeighbours(max_cutoff);
	}

	unsigned int nbNMax = _MySystem->getNbMaxN();
   
        // The following mainly comes from the compute_atom() function of the ACEBEvaluator (l.102 ace/ace_b_evaluator.cpp)	
	for(unsigned int i=0;i<nbDatTot;i++){ // parallelize ?
		unsigned int count_dim = 0;
    		Array1D<ACEComplex> A_cache(MyACEBBase->rankmax);
    		ACEComplex B{0.};
    		Array1D<ACEComplex> A_forward_prod(MyACEBBase->rankmax + 1);
    		//Array1D<ACEComplex> A_backward_prod(MyACEBBase->rankmax + 1);
		A.fill({0});
		A_rank1.fill(0);
		const Array2DLM<ACEComplex> &ylm = MyACEBBase->spherical_harmonics.ylm;
		const Array2D<DOUBLE_TYPE> &fr = MyACEBBase->radial_functions->fr;
		const Array1D<DOUBLE_TYPE> &gr = MyACEBBase->radial_functions->gr;
		ACEComplex Y{0};
		double R;

		double xpos = _MySystem->getWrappedPos(i).x;
		double ypos = _MySystem->getWrappedPos(i).y;
		double zpos = _MySystem->getWrappedPos(i).z;
		unsigned int mu_i = elem2pace[_MySystem->getAtom(i).type_uint-1];
		const SHORT_INT_TYPE total_basis_size = MyACEBBase->total_basis_size[mu_i];
		const SHORT_INT_TYPE total_basis_size_rank1 = MyACEBBase->total_basis_size_rank1[mu_i];
		//ALGORITHM 1: Atomic base A
		for(unsigned int j_loop=0;j_loop<_MySystem->getNeighbours(i*(nbNMax+1));j_loop++){
			unsigned int id = _MySystem->getNeighbours(i*(nbNMax+1)+j_loop+1);
			unsigned int ind1 = i*nbNMax*3+j_loop*3;
			unsigned int ind2 = ind1+1;
			unsigned int ind3 = ind2+1;
			// get distance vector
			double xp = _MySystem->getWrappedPos(id).x+_MySystem->getCLNeighbours(ind1)*_MySystem->getH1()[0]+_MySystem->getCLNeighbours(ind2)*_MySystem->getH2()[0]+_MySystem->getCLNeighbours(ind3)*_MySystem->getH3()[0]-xpos;
			double yp = _MySystem->getWrappedPos(id).y+_MySystem->getCLNeighbours(ind1)*_MySystem->getH1()[1]+_MySystem->getCLNeighbours(ind2)*_MySystem->getH2()[1]+_MySystem->getCLNeighbours(ind3)*_MySystem->getH3()[1]-ypos;
			double zp = _MySystem->getWrappedPos(id).z+_MySystem->getCLNeighbours(ind1)*_MySystem->getH1()[2]+_MySystem->getCLNeighbours(ind2)*_MySystem->getH2()[2]+_MySystem->getCLNeighbours(ind3)*_MySystem->getH3()[2]-zpos;

			double r_norm = sqrt(xp*xp + yp*yp + zp*zp);
			unsigned int mu_j = elem2pace[_MySystem->getAtom(id).type_uint-1];
			double current_cutoff = MyACEBBase->radial_functions->cut(mu_i,mu_j);
			if( r_norm < current_cutoff ){
				xp /= r_norm;
				yp /= r_norm;
				zp /= r_norm;
				
				MyACEBBase->radial_functions->evaluate(r_norm, MyACEBBase->nradbase, nradiali, mu_i, mu_j);
				MyACEBBase->spherical_harmonics.compute_ylm(xp, yp, zp, lmaxi);
				//loop for computing A's
				//rank = 1
				for(unsigned int n=0;n<MyACEBBase->nradbase;n++){
					//GR = gr(n);
					//A_rank1(mu_j, n) += GR * Y00;
					A_rank1(mu_j, n) += gr(n) * Y00;
				}
				//loop for computing A's
				// for rank > 1
				for(unsigned int n=0;n<nradiali;n++){
					auto &A_lm = A(mu_j, n);
					for(unsigned int l=0;l<=lmaxi;l++){
						R = fr(n, l);
						
						for(unsigned int m=0;m<=l;m++){
							Y = ylm(l, m);
							A_lm(l, m) += R * Y; //accumulation sum over neighbours
						}
					}
				}
			}
		} //end loop over neighbours
		for(unsigned int mu_j=0;mu_j<nbAtomType;mu_j++){
			for(unsigned int n=0;n<nradiali;n++){
				auto &A_lm = A(mu_j, n);
				for(unsigned int l=0;l<=lmaxi;l++){
					//fill in -m part in the outer loop using the same m <-> -m symmetry as for Ylm
					for(unsigned int m=1;m<=l;m++){
						int factor = m % 2 == 0 ? 1 : -1;
						A_lm(l, -m) = A_lm(l, m).conjugated() * factor;
					}
				}
			}
		}    //now A's are constructed
		// For rank = 1, A are still invariants by rotation => directly store them in Descriptor array
		auto basis_rank1 = MyACEBBase->basis_rank1[mu_i];
		for(unsigned int n=0;n<total_basis_size_rank1;n++){
			auto func = &basis_rank1[n];
			_Descriptors[i*dim+count_dim] = A_rank1(func->mus[0], func->ns[0] - 1);
			if( count_dim == dim ) cout << "ISSUE !" << endl;
			count_dim++;
		}

		// rank > 1
		auto basis = MyACEBBase->basis[mu_i];
		for(unsigned int func_ind=0;func_ind<total_basis_size;++func_ind){
			auto func = &basis[func_ind];
			int rank = func->rank; // to see with type of vars
			int r = rank - 1;
			int *mus = func->mus;
			short int *ns = func->ns;
			short int *ls = func->ls;
			
			//loop over {ms} combinations in sum
			for(unsigned int ms_ind=0;ms_ind<func->num_ms_combs;++ms_ind){
				short int *ms = &func->ms_combs[ms_ind * rank]; // current ms-combination (of length = rank)
				
				//loop over m, collect B  = product of A with given ms
				A_forward_prod(0) = 1;
				
				//fill forward A-product triangle
				unsigned int t;
				//cout << "for func_ind " << func_ind << " ( => " << func->num_ms_combs << " number of comb" << endl;
				//cout << "COMB = " << ms_ind << endl;
				for(t=0;t<rank;t++){
					//cout << ns[t]-1 << " " << ls[t] << " " << ms[t] << endl; 
					A_cache(t) = A(mus[t], ns[t] - 1, ls[t], ms[t]);
					A_forward_prod(t + 1) = A_forward_prod(t) * A_cache(t);
				}
				
				B = A_forward_prod(t);
				_Descriptors[i*dim+count_dim] = B.real;
				if( count_dim == dim ) cout << "ISSUE !" << endl;
				count_dim++; 
			}
		} // end loop on total basis size
	} // end loop on atoms

	delete MyACEBBase;
	delete[] elem2pace;
}


ACEDescriptors::~ACEDescriptors(){
}
