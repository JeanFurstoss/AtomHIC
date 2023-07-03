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
	unsigned int *Atom_SiteIndex; // determined using comparison of Steinhardt param for l=l_sph (using the database)
	unsigned int *Malpha;// array containing the index of neighbours of the same species (or same site in case of multisite crystal) with the first line corresponding to the number of neighbours, i.e. Malpha[i*(nbNMax+1)] = nb of neighbour of atom i, Malpha[i*(nbNMax+1)+j+1] = id of the jth neighbour of atom i
	std::complex<double> *Qalpha; // complex array containing the spherical harmonic for the different modes
	double *SteinhardtParams; // Steinhardt parameters
	double *SteinhardtParams_ave_cutoff; // Steinhardt parameters
	std::vector<std::vector<double>> SteinhardtParams_REF_PC; // Steinhardt parameters for reference perfect crystal
	std::vector<std::vector<double>> BondOriParam_REF_PC; // Bond ori param for reference perfect crystal (non centrosymmetric crystals)
	std::vector<double*> SteinhardtParams_REF_PC_ave; // Steinhardt parameters for reference perfect crystal averaged over sites
	std::vector<double*> SteinhardtParams_REF_PC_ave_cutoff; // Steinhardt parameters for reference perfect crystal average over the atom of the same species within the cutoff radius
	std::vector<double*> SteinhardtParams_REF_Def; // Steinhardt parameters for reference defects => tab[i*(l_sph_ref+1)+l] gives for ref i, the lth degree of steinhard parameter
	std::vector<double*> SteinhardtParams_REF_Def_ave_cutoff; // Steinhardt parameters for reference defects => tab[i*(l_sph_ref+1)+l] gives for ref i, the lth degree of steinhard parameter
	std::vector<unsigned int*> AtomTypeUINTRefDef; // array containing the number of ref and atom specy => i.e. AtomTypeUINTRefDef[i][0] = number of ref of defect i, AtomTypeUINTRefDef[i][j+1] = type uint of the jth ref of the defect i
	std::vector<unsigned int*> AtomTypeUINTRefDef_ave_cutoff; // array containing the number of ref and atom specy => i.e. AtomTypeUINTRefDef[i][0] = number of ref of defect i, AtomTypeUINTRefDef[i][j+1] = type uint of the jth ref of the defect i
	std::vector<std::string> Ref_Def_Names; // Names of the defect present in the database
	std::vector<int> AtomTypeUINTRefPC; // array containing the atom specy => i.e. AtomTypeUINTRefPC[i] = type uint of the params in SteinhardtParams_REF_DEF[i]
	std::vector<std::vector<unsigned int>> AtomSiteRefPC; // array containing the crystallographic site => i.e. AtomSiteRefPC[i][j] = crystallographic site of the params in SteinhardtParams_REF_PC[i]
	std::vector<int> AtomTypeUINTRefPC_ave; // array containing the atom specy => i.e. AtomTypeUINTRefPC[i] = type uint of the params in SteinhardtParams_REF_DEF[i]
	unsigned int nbRefDef; // number of defect in the database
	unsigned int nbref; // number of q vectors for each defect or perfect crystal average
	int l_sph_ref;
	double rcut_ref;
	double *Calpha; // normalization factor for bond orientational parameter
	bool IsBondOriParam = false;
	bool IsStrainTensor = false;
	bool IsStrainInvII = false;
	bool IsSteinhardtDatabaseRead = false;
	bool IsSteinhardt = false;
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
	double* BondOrientationalParameter();
	void BondOriParam_MultisiteCrystal();
	void BondOriParam_Multisite();
	void BondOriParam_MultisiteNewVersion();
	void BondOriParam_NoMultisite();
	double* ComputeSteinhardtParameters(const double rc, const int l_sph);
	double* ComputeSteinhardtParameters_Multi(const double rc, const int l_sph);
	double* ComputeSteinhardtParameters_testMono(const double rc, const int l_sph);
	double* ComputeSteinhardtParameters_testMulti(const double rc, const int l_sph);
	double* ComputeSteinhardtParameters_OneL(const double rc, const int l_sph);
	double* Compute_StrainTensor();
	double* Compute_StrainTensor(unsigned int FromNum);
	double* Compute_StrainTensor_invII();
	unsigned int *get_AtomSiteIndex(){ return this->Atom_SiteIndex; }
	unsigned int getNumberRefDef_Steinhardt(){ return this->Ref_Def_Names.size(); }
	double *get_StrainTensor(){ return this->StrainTensor; }
	double *get_StrainInvII(){ return this->Strain_invII; }
	std::complex<double> spherical_harmonics(const unsigned int& l, int& m, double& theta, double& phi);
	void SaveSteinhardtParamToDatabase_PerfectCrystal(std::string CrystalName);
	void SaveAveSteinhardtParamToDatabase_PerfectCrystal(std::string CrystalName);
	void SaveSteinhardtParamToDatabase_Defect(std::string CrystalName, std::string filename, std::vector<unsigned int> At_index);
	//std::string SteinhardtDatabase_write(std::string CrystalName);
	std::string getSteinhardtDatabase(std::string CrystalName);
	void SteinhardtDatabase_read(std::string CrystalName);
	void PrintSteinhardtParam(std::vector<unsigned int> At_index, std::string ext_filename);
	double* BondOriParam_SteinhardtBased();
	double* StructuralAnalysis_Steinhardt();
	void read_params();
	// destructor
	~ComputeAuxiliary();
	
};

#endif
