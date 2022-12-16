#ifndef BICRYSTAL_H
#define BICRYSTAL_H

#include "AtomHicExport.h"
#include "Crystal.h"
#include "AtomicSystem.h"
#include "ComputeAuxiliary.h"
#include <string>

class ATOMHIC_EXPORT Bicrystal : public AtomicSystem {
protected:
	std::string NormalDir;
	double Ldir; // box length in the GB normal direction
	double ExcessVol;
	bool IsVacuum;
	bool IsCentered; // in the case where there is vacuum, is the system centered
	double VacuumLo; // low bound coordinate of the vacuum
	double VacuumHi; // high bound coordinate of the vacuum
	double GBPos1; // position of the first GB
	double GBPos2; // position of the second GB (if there is no vacuum)
	double GBwidth1; // width of the first GB
	double GBwidth2; // width of the second GB (if there is no vacuum)
	double MaxPos; // maximum pos of atoms
	double MinPos; // min --
	double SystemLength;
	ComputeAuxiliary *CA;
	bool IsCA = false;
	double prefac_test;
	bool IsMassDensity;
	Crystal *_MyCrystal2;
	AtomicSystem *Grain1;
	AtomicSystem *Grain2;
	Atom *AtomList_G1;
	double *H1_G1, *H2_G1, *H3_G1;
	Atom *AtomList_G2;
	double *H1_G2, *H2_G2, *H3_G2;
	double *RotGrain1ToGrain2;
	double *RotCartToGrain2;
	bool IsRotMatDefine = false;
	bool AreGrainsDefined = false;
	bool IsCrystal2 = false;
	double xl1, xl2, yl1, yl2;
	std::vector<double*> CSL;
	double *CSL_Basis;
	bool IsCSL = false;
	bool IsCSL_Basis = false;
	unsigned int dupX1, dupX2, dupY1, dupY2; // number of duplicate in X and Y dir for crystal 1 and 2 for constructing the bicrystal
	double Mx1, Mx2, My1, My2; // misfit applied to crystal 1 and 2 in dir X and Y for constructing the bicrystal
public:
	// constructors
	Bicrystal(){};
	Bicrystal(const std::string& filename, const std::string NormalDir);
	Bicrystal(const std::string& filename, const std::string NormalDir, const std::string CrystalName);
	Bicrystal(const std::string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p);// Constructor for bicrystal with plane GB with given misorientation and GB plane
	Bicrystal(const std::string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p, std::vector<int> FacetsType, unsigned int N_facet);// Constructor for bicrystal with facetted GB with given misorientation and GB plane and facet type
	// getters
	double getGBPos1(){ return this->GBPos1; }
	double getGBwidth1(){ return this->GBwidth1; }
	double getPrefac1(){ return this->prefac_test; }
	double getExcessVol() { return this->ExcessVol; }
	Crystal* getCrystal2(){ return this->_MyCrystal2; }
	double getxl1(){ return this->xl1; }
	double getyl1(){ return this->yl1; }
	double getxl2(){ return this->xl2; }
	double getyl2(){ return this->yl2; }
	// methods
	void print_Grains();
	void searchGBPos();
	void ComputeExcessVolume();
        void searchCSL(unsigned int verbose=0);
	void generateCSL();
	void solve_DSC(const int *u, const unsigned int L, const double *B, double *DSC_Base, double tol);
	void setOrientedCrystals(const std::string& crystalName, int h_a, int k_a, int l_a, double theta, int h_p, int k_p, int l_p); // method initializing 2 crystals with a given misorientation relationship and plane
	void printCSL(const std::string filename);
	// destructor
	~Bicrystal();
};

#endif
