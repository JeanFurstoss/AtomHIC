#ifndef KMEANS_H
#define KMEANS_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "MachineLearningModel.h"

class ATOMHIC_EXPORT KMeans : public MachineLearningModel {
private:
	bool IsDescriptor = false;
	unsigned int *nbClust; // [f]
	double *centroids; // centroids[k*dim*nbFilter+d*nbFilter+f] d component of the centroid of cluster k with filter f
	long double *V; // V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+f] d1,d2 component of variance of cluster k with filter f
	long double *V_inv; 
	long double *det_V; 
	unsigned int *Data2Cluster; // Data2Cluster[i*nbFilter+f] = id of the cluster to which the data point i belongs to with filter f
	unsigned int *nbDat2Cluster; // nbDat2Cluster[k*nbFilter+f] = number of points belonging to cluster k with filter f
	long double *LogLikelihood;
	double *BIC;
	double *centroids_old; // buffer variable for training
	bool IsInitialized = false;

	// Saved variables for the different initialization
	long double *saved_LogLikelihood;
	unsigned int *saved_nbDat2Cluster; // [c*nbClustMax*nbFilter+k*nbFilter+f]
	long double *saved_V; // [c*nbClustMax*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+f]
	long double *saved_V_inv; 
	long double *saved_det_V; 
	double *saved_centroids;
	bool SavedVariablesInitialized = false;

	double *buffer_k;
	double *buffer_dat;

	// FixedParameters
	unsigned int nbClustMax;
	double tol_KMeans;
	unsigned int MaxIter_KMeans;
	unsigned int nbInit; // number of random initialization (we keep at the end the one with the highest likelihood)

public:
	// constructors
	KMeans();
	// methods
	void readFixedParams();
	void setDescriptors(Descriptors *D);
	double SquareEuclidianDistance(const unsigned int &DescriptorIndex, const unsigned int &ClusterIndex, unsigned int &filter_value);
	void TrainModel(unsigned int &_nbClust, unsigned int &filter_value);
	void AffectData2Cluster(unsigned int &filter_value);
	void ComputeLogLikelihood(unsigned int &filter_value);
	void ComputeBIC(unsigned int &filter_value);
	void PrintModelParams(unsigned int &filter_value);
	void KMeansPPInitialization(unsigned int &_nbClust, unsigned int &filter_value);
	void SaveVariables(unsigned int &current, unsigned int &filter_value);
	void ComputeFullVariances(unsigned int &filter_value);
	double ComputeGaussianProb(unsigned int &DescriptorIndex, unsigned int &filter_value);
	double *getBIC(){ return BIC; }
	long double *getLogLikelihood(){ return LogLikelihood; }
	unsigned int *getData2Cluster(){ return Data2Cluster; }
	unsigned int *getNbDat2Cluster(){ return nbDat2Cluster; }
	double *getCentroids(){ return centroids; }
	long double *getV(){ return V; }
	long double *getV_inv(){ return V_inv; }
	long double *getDet_V(){ return det_V; }
	void ReadProperties(std::vector<std::string> Properties);
	//void ReadModelParamFromDatabase();
	// destructor
	~KMeans();
	
};

#endif
