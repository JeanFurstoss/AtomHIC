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
	std::string crystallo;
	std::string path2database;
	unsigned int nbAtom;
	unsigned int nbAtomType;
	std::string *AtomType; // the different atom species present in the system
	unsigned int *AtomType_uint;
	unsigned int *NbAtomSite; // array such as NbAtomSite[i] = number of site for atom type AtomType[i]
	bool IsCharge;
	double *AtomMass;
	// Cell vectors
	double *a1;
	double *a2;
	double *a3;
	double a1length, a2length, a3length;
	// reciproqual lattice
	double *a1_star;
	double *a2_star;
	double *a3_star;
	double V; // elementary volum
	Atom *Motif;
	const std::string database_extension=".dat";
	bool IsOrientedSystem;
	AtomicSystem *OrientedSystem;
	MathTools *MT;
	std::vector<std::vector<unsigned int>> DoNotSep; // contain the number neighbour which should not be separe from a given atom type <=> do not separe atom type DoNotSep[i][0] from its DoNotSep[i][1] first neighbors of atom type DoNotSep[i][2]
	std::vector<std::vector<int>> NotSepList; // contain the id of neighbors which should not be separe from the given atom, NotSepList[i][j*4] = id of neighbor j which should not be separe from atom i, NotSepList[i][j*4+k] = multiplicative coefficient with ak(a1, a2, a3) for this atom to be the neighbot
	bool IsDoNotSep = false;
	double *rot_mat_total; // rotation matrix used to pass from the database oriented crystal to the current orientation
	double *TiltTrans_xyz; // transformation matrix used to construct oriented systems 
public:
	// constructors
	Crystal(){};
	Crystal(const std::string& crystalName); // constructor searching in the database the different parameter for the construction
	// getters
	const std::string& getName(){ return this->name; }
	const Atom* getMotif(){ return this->Motif; }
	const double* getA1(){ return this->a1; };
	const double* getA2(){ return this->a2; };
	const double* getA3(){ return this->a3; };
	const double getVol(){ return this->V; };
	const unsigned int getNbAtom(){ return this->nbAtom; }
	const unsigned int getNbAtomType(){ return this->nbAtomType; }
	const unsigned int getAtomType_uint(const unsigned int Id){ return this->AtomType_uint[Id]; }
	inline const unsigned int getNbAtomSite(const unsigned int type_uint){
		unsigned int index=0;
		for(unsigned int i=0;i<this->nbAtomType;i++){
		       if( type_uint == this->AtomType_uint[i] ){
			       index = i;
			       break;
		       }
		}
		return this->NbAtomSite[index];
	}
	const std::string getAtomType(const unsigned int Id){ return this->AtomType[Id]; }
	const double getAtomMass(const unsigned int Id){ return this->AtomMass[Id]; }
	const bool getIsCharge(){ return this->IsCharge; }
	AtomicSystem* getOrientedSystem(){ return this->OrientedSystem; } 
	const double* getRotMat(){ return this->rot_mat_total; };
	bool getIsDoNotSep(){ return this->IsDoNotSep; }
	unsigned int getNotSepList_size(const unsigned int i){ return this->NotSepList[i].size()/4; }
	int getNotSepList(const unsigned int i, const unsigned int j){ return this->NotSepList[i][j]; }
	std::vector<std::vector<unsigned int>> getDoNotSep(){ return this->DoNotSep; }
	// methods
	void read_database();
	void computeReciproqual();
	void ConstructOrientedSystem(const int& h_p, const int& k_p, const int& l_p);
	void ConstructOrientedSystem(const double *RotMat);
	void ConstructNotSepList();
	void RotateAndConstructOrthogonalCell(const double *RotMat, double &xbox, double &ybox, double &zbox, std::vector<int> &cl_box);
	// destructor
	~Crystal();
	
};

#include "AtomicSystem.h"

#endif
