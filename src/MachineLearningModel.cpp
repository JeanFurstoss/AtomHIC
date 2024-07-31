#include "MachineLearningModel.h"
#include "AtomHicConfig.h"

MachineLearningModel::MachineLearningModel(){
	MT = new MathTools;
}

void MachineLearningModel::setDescriptors(Descriptors *D){ 
	_MyDescriptors = D;
	this->IsDescriptor = true;
}

MachineLearningModel::~MachineLearningModel(){
	delete MT;
}
