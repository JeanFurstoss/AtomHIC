#ifndef DESCRIPTORS_H
#define DESCRIPTORS_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>

class ATOMHIC_EXPORT Descriptors {
private:
	std::string name;
	AtomicSystem *_MySystem;
	MathTools *MT;
	unsigned int dim; // dimension of descriptors
	unsigned int nbDat; // number of data
	double *_Descriptors; // [i*dim+d] component d of descriptor i

public:
	// constructors
	Descriptors(){};
	Descriptors(const std::string& Filename); // constructor from file reading 
	// methods
	double *getDescriptors(){ return _Descriptors; }
	unsigned int getDim(){ return dim; }
	unsigned int getNbDat(){ return nbDat; }
	// destructor
	~Descriptors();
	
};

#endif
