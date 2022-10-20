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
	// methods
	void read_database();
	// destructor
	~Crystal();
	
};

#endif
