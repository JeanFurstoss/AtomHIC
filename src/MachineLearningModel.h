#ifndef MACHINELEARNINGMODEL_H
#define MACHINELEARNINGMODEL_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "Descriptors.h"

class ATOMHIC_EXPORT MachineLearningModel {
protected:
	std::string name;
	Descriptors *_MyDescriptors;
	bool IsDescriptor = false;
	MathTools *MT;

	double *buffer_vec_1_dim;
	double *buffer_vec_2_dim;
	
	unsigned int dim;
	unsigned int dim2;
	unsigned int nbDat;

public:
	// constructors
	MachineLearningModel();
	// methods
	void setDescriptors(Descriptors *D);
	void TrainModel();
	//void ReadModelParamFromDatabase();
	//void Predict();
	// destructor
	~MachineLearningModel();
	
};

#endif
