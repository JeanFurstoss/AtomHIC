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
	unsigned int nbDatMax;
	unsigned int *nbDat; // [f] number of data with filter f
	unsigned int nbFilter; // number of filter (e.g. if filtered by element, the number of element)
	// Variables for labelling the model
	bool IsLabelled = false;
	unsigned int nbLabel;
	
	// Variables for reading the database parameters
	std::vector<std::string> Labels; // Labels[l] = name of the label l
	std::string FilteringType;
	std::vector<std::string> FilterValue; // FilterValue[f] = name of the filter f
	std::string DescriptorName; // to be used for initializing descriptors
	std::vector<std::string> DescriptorProperties; // to be used for initializing descriptors
	unsigned int *FilterIndexToModify;
	bool IsFilterIndexModified = false;
	bool IsRead = false;

	// Variables for classification
	double *Classificator; // [n*2] = index of *the label with highest probability for descriptors n (index in _Descriptor array), [n*2+1] = probability
	bool IsClassified = false;	
	
	std::string FixedParam_Filename = "Fixed_Parameters.dat";

public:
	// constructors
	MachineLearningModel();
	// methods
	void setDescriptors(Descriptors *D);
	void TrainModel();
	std::string getMLDatabasePath();
	std::string getDatabasePath(const std::string &name_of_database);
	void ChangeFilterIndex();
	//void ReadModelParamFromDatabase(); // TODO implement base fonction here
	void Classify();
	void PrintClassifiedData(std::string filename);
	void readFixedParams();
	// getters
	double *getClassificator(){ return Classificator; }
	std::string getFilteringType(){ return FilteringType; }
	std::vector<std::string> getDescriptorProperties(){ return DescriptorProperties; }
	// destructor
	~MachineLearningModel();
	
};

#endif
