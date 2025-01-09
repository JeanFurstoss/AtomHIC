// Jean Furstoss 
// https://www.gnu.org/licenses/
// DBScan.h

/* This class provides many .. 
   
   email.. 
   25 August 2024
 */


#ifndef DBSCAN_H
#define DBSCAN_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "MachineLearningModel.h"

class ATOMHIC_EXPORT DBScan : public MachineLearningModel {
private:
	bool IsDescriptor = false;
	unsigned int *nbClust; // [f]
	double *centroids; // centroids[k*dim*nbFilter+d*nbFilter+f] d component of the centroid of cluster k with filter f
	AtomicSystem *_MySystem;

	// Fixed parameters
	double eps = 3.;
	unsigned int minPts = 5;

public:
	// constructors
	DBScan();
	// methods
	void readFixedParams();
	void setDescriptors(Descriptors *D);
	void TrainModel(std::string filter_name);
	void Classify(){};
	void ReadProperties(std::vector<std::string> Properties);
	//void ReadModelParamFromDatabase();
	// destructor
	~DBScan();
	
};

#endif
