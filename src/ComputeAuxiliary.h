#ifndef COMPUTEAUXILIARY_H
#define COMPUTEAUXILIARY_H

#include "SteinhardtDescriptors.h"
#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>

class ATOMHIC_EXPORT ComputeAuxiliary {
protected:
	AtomicSystem *_MySystem;
	MathTools *MT;
	SteinhardtDescriptors *StDes;
	bool IsSteinhardtDescriptor = false;
	unsigned int *Atom_SiteIndex; // determined using comparison of Steinhardt param for l=l_sph (using the database)
	bool IsAtomSiteIndex = false;
	double *StrainTensor;
	double *Strain_invII;

	bool IsStrainTensor = false;
	bool IsStrainInvII = false;
	// For atomic strain and D2Min
	double *Ji; // affine transformation matrix
	double *d0; // reference state variable for affine transfo matrix
	double *current_d; // current delta for affine transfo matrix
	double *Vi_inv; // reference state variable for affine transfo matrix
	double *AtomicStrain;
	double *D2Min;
	bool Reference_AtomicStrain_Computed = false;
	bool IsAtomicStrain = false;
	bool IsJi = false;
	bool IsD2Min = false;
	// Parameters to read
	std::string FixedParam_Filename = "Fixed_Parameters.dat";
	double tolSites;
public:
	// constructors
	ComputeAuxiliary(){};
	ComputeAuxiliary(AtomicSystem *_MySystem): _MySystem(_MySystem){
		this->MT = new MathTools;
		read_params();
	};
	// getters
	// methods
	void ComputeAtomSiteIndex();
	double* BondOrientationalParameter();
	double* Compute_AffineTransfoMatrix(AtomicSystem &ReferenceSystem, double rc); // Ji matrix as defined in ovito (i.e. Shimizu, Ogata, Li: Mater. Trans. 48 (2007), 2923)
	double* Compute_AtomicStrain(AtomicSystem &ReferenceSystem, double rc); // atomic strain as defined in ovito (i.e. Shimizu, Ogata, Li: Mater. Trans. 48 (2007), 2923)
	double* Compute_D2Min(AtomicSystem &ReferenceSystem, double rc); // D2Min as defined in Delbecq et al. 2023
	double* Compute_StrainTensor();
	double* Compute_StrainTensor(unsigned int FromNum);
	double* Compute_StrainTensor_invII();
	unsigned int *get_AtomSiteIndex(){ return this->Atom_SiteIndex; }
	double *get_StrainTensor(){ return this->StrainTensor; }
	double *get_StrainInvII(){ return this->Strain_invII; }
	void read_params();
	// destructor
	~ComputeAuxiliary();
	
};

#endif
