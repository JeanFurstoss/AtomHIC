#ifndef ATOMHIC_MATHTOOLS_H
#define ATOMHIC_MATHTOOLS_H

#include "AtomHicExport.h"

class ATOMHIC_EXPORT MathTools {
public:
	MathTools(){};
	const double max(const double arr [], unsigned int size);
	const double min(const double arr [], unsigned int size);
	~MathTools(){};
};

#endif
