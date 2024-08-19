#ifndef CRYSTAL_H
#define CRYSTAL_H

#include "AtomHicExport.h"
#include "MyStructs.h"
#include <string>
#include "MathTools.h"

class AtomicSystem;

class ATOMHIC_EXPORT Crystal {
private:
	std::string name;
	std::string crystallo; // cubic, orthorhombic, etc.
	std::string path2database;
	unsigned int nbAtom;
	unsigned int nbAtomType;
	std::string *AtomType; // the different atom species present in the system
	unsigned int *AtomType_uint;
	unsigned int *AtomSite;
	unsigned int *NbAtomSite; // array such as NbAtomSite[i] = number of site for atom type AtomType[i]
	bool IsCharge = false;
	double *AtomMass;
	double *AtomCharge;
	// Cell vectors
	double *a1;
	double *a2;
	double *a3;
	double *alength;
	// reciproqual lattice
	double *a1_star;
	double *a2_star;
	double *a3_star;
	double V; // elementary volum
	Atom *Motif;
	bool IsOrientedSystem = false;
	bool IsMultisite = false;
	AtomicSystem *OrientedSystem;
	MathTools *MT;
	std::vector<std::vector<unsigned int>> DoNotSep; // contain the number neighbour which should not be separe from a given atom type <=> do not separe atom type DoNotSep[i][0] from its DoNotSep[i][1] first neighbors of atom type DoNotSep[i][2]
	std::vector<std::vector<int>> NotSepList; // contain the id of neighbors which should not be separe from the given atom, NotSepList[i][j*4] = id of neighbor j which should not be separe from atom i, NotSepList[i][j*4+k] = multiplicative coefficient with ak(a1, a2, a3) for this atom to be the neighbot
	bool IsDoNotSep = false;
	double *rot_mat_total; // rotation matrix used to pass from the database oriented crystal to the current orientation
	double *TiltTrans_xyz; // transformation matrix used to construct oriented systems
        unsigned int *Stoichiometry;
	bool IsReferenceBondOriParam = false;
	std::vector<std::vector<double>> ReferenceBondOriParam; // ReferenceBondOriParam[t][s] = reference bond orientational parameter (used for searching atom site and computing order parameter of multisite and non centrosymmetric crystal) of atom type t and crystallographic site s
	std::vector<std::string> BondOriParamProperties;
	// Parameters to read
	std::string FixedParam_Filename = "Fixed_Parameters.dat";
	double TolOrthoBox;
	double TolOrthoBoxZ;
	double MinBoxHeight;
	double MinBoxAside;
public:
	// constructors
	Crystal(){};
	Crystal(const std::string& crystalName); // constructor searching in the database the different parameter for the construction
	// getters
	const std::string& getName(){ return this->name; }
	const std::string& getCrystallo(){ return this->crystallo; }
	Atom* getMotif(){ return this->Motif; } // TODO warning the const has been removed here
	double* getA1(){ return this->a1; };
	double* getA2(){ return this->a2; };
	double* getA3(){ return this->a3; };
	const double* getA1_star(){ return this->a1_star; };
	const double* getA2_star(){ return this->a2_star; };
	const double* getA3_star(){ return this->a3_star; };
	const double* getALength(){ return this->alength; };
	const double getVol(){ return this->V; };
	unsigned int getNbAtom(){ return this->nbAtom; }
	unsigned int getNbAtomType(){ return this->nbAtomType; }
	unsigned int getAtomSite(const unsigned int i){ return this->AtomSite[i]; }
	const bool getIsMultisite(){ return this->IsMultisite; }
	unsigned int getNbAtomSite(const unsigned int typeuint){ return this->NbAtomSite[typeuint-1]; }
	unsigned int *getNbAtomSite(){ return this->NbAtomSite; }
	const double* getTiltTrans(){ return this->TiltTrans_xyz; }
	const std::string getAtomType(const unsigned int typeuint){ return this->AtomType[typeuint-1]; }
	std::string* getAtomType(){ return this->AtomType; }
	const double getAtomMass(const unsigned int typeuint){ return this->AtomMass[typeuint-1]; }
	double* getAtomMass(){ return this->AtomMass; }
	const double getAtomCharge(const unsigned int typeuint){ return this->AtomCharge[typeuint-1]; }
	double* getAtomCharge(){ return this->AtomCharge; }
	const bool getIsCharge(){ return this->IsCharge; }
	AtomicSystem* getOrientedSystem(){ return this->OrientedSystem; } 
	const double* getRotMat(){ return this->rot_mat_total; };
	bool getIsDoNotSep(){ return this->IsDoNotSep; }
	unsigned int getNotSepList_size(const unsigned int i){ return this->NotSepList[i].size()/4; }
	int getNotSepList(const unsigned int i, const unsigned int j){ return this->NotSepList[i][j]; }
	std::vector<std::vector<unsigned int>> getDoNotSep(){ return this->DoNotSep; }
	unsigned int *getStoich(){ return this->Stoichiometry; }
	std::string getDatabasePath(std::string crystalName);
	std::vector<std::string> getBondOriParamProperties(){ return BondOriParamProperties; }
	std::vector<std::vector<double>> getReferenceBondOriParam(){ return ReferenceBondOriParam; }
	bool getIsReferenceBondOriParam(){ return IsReferenceBondOriParam; }
	// methods
	void read_params();
	void read_database();
	void computeReciproqual();
	void RotateCrystal(const int& h_p, const int& k_p, const int& l_p);
	void RotateCrystal(const double *RotMat);
	void ConstructNotSepList();
	void ConstructOrthogonalCell();
	void computeStoich();
	// destructor
	~Crystal();
	
};

#include "AtomicSystem.h"

#endif
