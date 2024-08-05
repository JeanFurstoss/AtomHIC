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
	std::vector<std::string> Labels;
	unsigned int *Labels_uint;

public:
	// constructors
	Descriptors(){};
	Descriptors(const std::string& FilenameOrDir); // constructor from file reading 
	Descriptors(const std::string& FilenameOrDir, const std::string& DescriptorName); // constructor from file reading 
	// methods
	double *getDescriptors(){ return _Descriptors; }
	unsigned int *getLabels_uint(){ return Labels_uint; }
	std::string getLabels(unsigned int &l){ return Labels[l]; }
	unsigned int getDim(){ return dim; }
	unsigned int getNbDat(){ return nbDat; }
	// destructor
	~Descriptors();
	
};

#endif
