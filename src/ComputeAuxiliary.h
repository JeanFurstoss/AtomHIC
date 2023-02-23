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
	double *BondOriParam_Steinhardt;
	double *StrainTensor;
	double *Strain_invII;
	unsigned int *Atom_SiteIndex;
	unsigned int *Malpha;// array containing the index of neighbours of the same species (or same site in case of multisite crystal) with the first line corresponding to the number of neighbours, i.e. Malpha[i*(nbNMax+1)] = nb of neighbour of atom i, Malpha[i*(nbNMax+1)+j+1] = id of the jth neighbour of atom i
	std::complex<double> *Qalpha; // complex array containing the spherical harmonic for the different modes
	double *SteinhardtParams; // Steinhardt parameters
	std::vector<double*> SteinhardtParams_REF_PC; // Steinhardt parameters for reference perfect crystal
	double *SteinhardtParams_REF_Def; // Steinhardt parameters for reference defects => tab[i][n*(l_sph_ref*2+1)+l] gives for ref i, the lth degree of steinhard parameter : norm (for n=0), real part (n=1), imag part (n=2)
	std::vector<std::string> Ref_Def_Names; // Names of the defect present in the database
	std::vector<int> AtomTypeUINTRefPC; // array containing the atom specy => i.e. AtomTypeUINTRefPC[i] = type uint of the params in SteinhardtParams_REF_DEF[i]
	int nbRefDef;
	int l_sph_ref;
	double rcut_ref;
	double *Calpha; // normalization factor for bond orientational parameter
	bool IsBondOriParam = false;
	bool IsStrainTensor = false;
	bool IsStrainInvII = false;
public:
	// constructors
	ComputeAuxiliary(){};
	ComputeAuxiliary(AtomicSystem *_MySystem): _MySystem(_MySystem){
		this->MT = new MathTools;
	};
	// getters
	// methods
	double* BondOrientationalParameter();
	double* ComputeSteinhardtParameters(const double rc, const int l_sph);
	double* Compute_StrainTensor();
	double* Compute_StrainTensor(unsigned int FromNum);
	double* Compute_StrainTensor_invII();
	unsigned int *get_AtomSiteIndex(){ return this->Atom_SiteIndex; }
	double *get_StrainTensor(){ return this->StrainTensor; }
	double *get_StrainInvII(){ return this->Strain_invII; }
	std::complex<double> spherical_harmonics(const unsigned int& l, int& m, double& theta, double& phi);
	void SaveSteinhardtParamToDatabase_PerfectCrystal(std::string CrystalName);
	void SaveSteinhardtParamToDatabase_Defect(std::string CrystalName, std::string filename);
	std::string SteinhardtDatabase_write(std::string CrystalName);
	void SteinhardtDatabase_read(std::string CrystalName);
	double* BondOriParam_SteinhardtBased();
	// destructor
	~ComputeAuxiliary();
	
};

#endif
