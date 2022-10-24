#ifndef ATOMHIC_CRYSTAL_H
#define ATOMHIC_CRYSTAL_H

#include "AtomHicExport.h"
#include "MyStructs.h"

#include <string>

class ATOMHIC_EXPORT Crystal {
protected:
	std::string name;
	std::string path2database;
	unsigned int nbAtom;
	unsigned int nbAtomType;
	std::string *AtomType; // the different atom species present in the system
	unsigned int *AtomType_uint;
	unsigned int *NbAtomSite; // array such as NbAtomSite[i] = number of site for atom type AtomType[i]
	bool IsCharge;
	double *AtomMass;
	double *CellParameters;
	Atom *Motif;
	const std::string database_extension=".dat";

public:
	// constructors
	Crystal(){};
	Crystal(const std::string& crystalName); // constructor searching in the database the different parameter for the construction
	// getters
	const std::string& getName(){ return this->name; }
	const double* getCellParameters(){ return this->CellParameters; }
	const Atom* getMotif(){ return this->Motif; }
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
	// destructor
	~Crystal();
	
};

#endif
