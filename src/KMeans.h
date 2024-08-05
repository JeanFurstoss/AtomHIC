#ifndef KMEANS_H
#define KMEANS_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "MachineLearningModel.h"

class ATOMHIC_EXPORT KMeans : public MachineLearningModel {
private:
	unsigned int nbClust;
	double *centroids; // centroids[k*dim+d] d component of the centroid of cluster k
	long double *V; // V[k*dim2+d1*dim+d2] d1,d2 component of variance of cluster k
	long double *V_inv; 
	long double *det_V; 
	unsigned int *Data2Cluster; // Data2Cluster[i] = id of the cluster to which the data point i belongs to
	unsigned int *nbDat2Cluster; // nbDat2Cluster[k] = number of points belonging to cluster k
	long double LogLikelihood;
	double BIC;
	double *centroids_old; // buffer variable for training
	bool IsInitialized = false;

	// Saved variables for the different initialization
	long double *saved_LogLikelihood;
	unsigned int *saved_nbDat2Cluster; // [c*nbClustMax+k]
	long double *saved_V; // [c*nbClustMax*dim2+k*dim2+d1*dim+d2]
	long double *saved_V_inv; 
	long double *saved_det_V; 
	double *saved_centroids;
	bool SavedVariablesInitialized = false;

	double *buffer_k;
	double *buffer_dat;

	// Parameters to put in FixedParameters
	unsigned int nbClustMax = 100;
	double tol_KMeans = 1e-5;
	unsigned int MaxIter_KMeans = 1000;
	unsigned int nbInit = 100; // number of random initialization (we keep at the end the one with the highest likelihood)

public:
	// constructors
	KMeans();
	// methods
	void setDescriptors(Descriptors *D);
	double SquareEuclidianDistance(const unsigned int &DescriptorIndex, const unsigned int &ClusterIndex);
	void TrainModel(unsigned int &_nbClust);
	void AffectData2Cluster();
	void ComputeLogLikelihood();
	void ComputeBIC();
	void PrintModelParams();
	void KMeansPPInitialization(unsigned int &_nbClust);
	void SaveVariables(unsigned int &current);
	void ComputeFullVariances();
	double ComputeGaussianProb(unsigned int &DescriptorIndex);
	double getBIC(){ return BIC; }
	double getLogLikelihood(){ return LogLikelihood; }
	unsigned int *getData2Cluster(){ return Data2Cluster; }
	unsigned int *getNbDat2Cluster(){ return nbDat2Cluster; }
	double *getCentroids(){ return centroids; }
	long double *getV(){ return V; }
	long double *getV_inv(){ return V_inv; }
	long double *getDet_V(){ return det_V; }
	//void ReadModelParamFromDatabase();
	// destructor
	~KMeans();
	
};

#endif
