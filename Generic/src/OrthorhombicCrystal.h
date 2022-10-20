#ifndef ATOMHIC_ORTHORHOMBICCRYSTAL_H
#define ATOMHIC_ORTHORHOMBICCRYSTAL_H

#include "AtomHicExport.h"
#include "Crystal.h"

#include <string>

class ATOMHIC_EXPORT OrthorhombicCrystal : public Crystal {
public:
	// constructors
	OrthorhombicCrystal(){};
	explicit OrthorhombicCrystal(const std::string& crystalName); // constructor searching in the database the different parameter for the construction
	// destructor
	~OrthorhombicCrystal(){};
};

#endif
