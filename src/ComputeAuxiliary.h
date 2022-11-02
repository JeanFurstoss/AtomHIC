#ifndef COMPUTEAUXILIARY_H
#define COMPUTEAUXILIARY_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include <complex>
#include <string>

class ATOMHIC_EXPORT ComputeAuxiliary {
protected:
	AtomicSystem *_MySystem;
	double *BondOriParam;
	bool IsBondOriParam = false;
public:
	// constructors
	ComputeAuxiliary(){};
	ComputeAuxiliary(AtomicSystem *_MySystem): _MySystem(_MySystem){};
	// getters
	// methods
	double* BondOrientationalParameter(const int& l, double& rc);
	std::complex<double> spherical_harmonics(const unsigned int& l, int& m, double& theta, double& phi);
	// destructor
	~ComputeAuxiliary();
	
};

#endif
