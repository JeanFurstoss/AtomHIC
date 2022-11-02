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
	double V; // elementary volum
	// reciproqual lattice
	double *a1_star;
	double *a2_star;
	double *a3_star;
	Atom *Motif;
	const std::string database_extension=".dat";
	bool IsOrientedPlane;
	AtomicSystem *OrientedPlane;
	MathTools *MT;
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
	const std::string* getAtomType(){ return this->AtomType; }
	// methods
	void read_database();
	void computeReciproqual();
	void ConstructOrientedPlane(const int& h_p, const int& k_p, const int& l_p);
	// destructor
	~Crystal();
	
};

#include "AtomicSystem.h"

#endif
