#ifndef COMPUTEAUXILIARY_H
#define COMPUTEAUXILIARY_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <complex>
#include <string>

class ATOMHIC_EXPORT ComputeAuxiliary {
protected:
	AtomicSystem *_MySystem;
	MathTools *MT;
	double *BondOriParam;
	double *StrainTensor;
	unsigned int *Atom_SiteIndex;
	bool IsBondOriParam = false;
	bool IsStrainTensor = false;
public:
	// constructors
	ComputeAuxiliary(){};
	ComputeAuxiliary(AtomicSystem *_MySystem): _MySystem(_MySystem){
		this->MT = new MathTools;
	};
	// getters
	// methods
	double* BondOrientationalParameter(const int& l, double& rc);
	double* Compute_StrainTensor();
	double* Compute_StrainTensor(unsigned int FromNum);
	unsigned int *get_AtomSiteIndex(){ return this->Atom_SiteIndex; }
	double *get_StrainTensor(){ return this->StrainTensor; }
	std::complex<double> spherical_harmonics(const unsigned int& l, int& m, double& theta, double& phi);
	// destructor
	~ComputeAuxiliary();
	
};

#endif
