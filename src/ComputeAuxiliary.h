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
	unsigned int *nbNeigh_FiltNeigh; // array containing the number of neighoring ions of particular type (i.e. nbNeigh_FiltNeigh[i*t] => nb of neighbor having type t of atom i) this array is used for the calculation of filtered steinhardt params
	std::complex<double> *Qalpha; // complex array containing the spherical harmonic for the different modes
	std::complex<double> *Qlm; // complex array containing the spherical harmonic for the different modes Qlm[i*(l_sph*2+1)*(l_sph+1)+l*(l_sph*2+1)+m] gives the spherical harmonic for atom i and degree l and m
	double *SteinhardtParams; // Steinhardt parameters
	double *SteinhardtParams_ave; // Steinhardt parameters
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
	
	// FOR Gaussian mixture model
	std::vector<std::string> Struct_GMM_Names; // name of the structures present in the database
	std::vector<std::vector<double>> *ICovs_GMM; // inverse covariant matrix of structures present in the database ( Covs_GMM[s*nbAtomType+t][i][j] => covariant matrix of structure s and atom type t )
	std::vector<double> *Mus_GMM; // esperance of structures present in the database ( Mus_GMM[s*nbAtomType+t][i] => esperance of structure s and atom type t )
	long double *Det_GMM; // determinant of the covariant matrix for structures present in the database ( Det_GMM[s*nbAtomType+t] => determinant of structure s and atom type t )

	// end for gaussian mixture model
	
	unsigned int nbRefDef; // number of defect in the database
	unsigned int nbref; // number of q vectors for each defect or perfect crystal average
	int l_sph_ref;
	double rcut_ref;
	double *Calpha; // normalization factor for bond orientational parameter
	bool IsBondOriParam = false;
	bool IsStrainTensor = false;
	bool IsStrainInvII = false;
	bool IsSteinhardtDatabaseRead = false;
	bool IsGMMSteinhardtDatabaseRead = false;
	bool IsSteinhardt_Mono = false;
	bool IsSteinhardt_Multi = false;
	bool IsSteinhardt_Filtered = false;
	bool IsSteinhardt_AveMono = false;
	bool IsSteinhardt_AveMulti = false;
	bool IsSteinhardt_AveFiltered = false;
	bool IsSteinhardt_FilteredAveFiltered = false;
	bool IsSteinhardt_AveFilteredMono = false;
	bool IsSteinhardt_AveFilteredMulti = false;
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
	void BondOriParam_NoMultisite();
	// The different Steinhardt parameters
	void ComputeSteinhardtParameters_Mono(const double rc, const int l_sph);
	void ComputeSteinhardtParameters_Multi(const double rc, const int l_sph);
	void ComputeSteinhardtParameters_FilteredNeigh(const double rc, const int l_sph);
	double* ComputeSteinhardtParameters_OneL(const double rc, const int l_sph);
	// The different averaging methods
	void AverageSteinhardtParameters_Mono(const double rc, const int l_sph);
	void AverageSteinhardtParameters_Multi(const double rc, const int l_sph);
	void AverageSteinhardtParameters_Filtered(const double rc, const int l_sph);
	void AverageFilteredSteinhardtParameters_Mono(const double rc, const int l_sph);
	void AverageFilteredSteinhardtParameters_Multi(const double rc, const int l_sph);
	void AverageFilteredSteinhardtParameters_Filtered(const double rc, const int l_sph);
	// computation of the different averaged Steinhardt parameters
	double *ComputeSteinhardtParameters(const double rc, const int l_sph, std::string SteinhardtStyle, std::string AveStyle);

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
	void SteinhardtDatabase_read_GMM(std::string CrystalName);
	void PrintSteinhardtParam(std::vector<unsigned int> At_index, std::string ext_filename, std::string SteinhardtStyle, std::string AveStyle);
	double* StructuralAnalysis_Steinhardt();
	double* StructuralAnalysis_Steinhardt_GMM();
	void read_params();
	// destructor
	~ComputeAuxiliary();
	
};

#endif
