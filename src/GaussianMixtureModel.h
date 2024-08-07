#ifndef GAUSSIANMIXTUREMODEL_H
#define GAUSSIANMIXTUREMODEL_H

#include "AtomHicExport.h"
#include "AtomicSystem.h"
#include "MathTools.h"
#include <string>
#include "MachineLearningModel.h"
#include "KMeans.h"

class ATOMHIC_EXPORT GaussianMixtureModel : public MachineLearningModel {
private:
	unsigned int *nbClust; // [f] number of cluster in the GMM with filter f
	unsigned int *nbDat; // [f] number of data with filter f
	double *weights; // weights[k*nbFilter+f] weight of cluster k with filter f
	double *mu; // mu[k*dim*nbFilter+d*nbFilter+f] d component of the mean of cluster k with filter f
	long double *V; // V[k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+f] d1,d2 component of variance of cluster k with filter f
	long double *V_inv;
	long double *det_V;
	long double *LogLikelihood;
	double *BIC;
	
	// old variables for EM iterations
	double *weights_old;
	double *mu_old;
	long double *V_old;
	long double *V_inv_old;
	long double *det_V_old;
	
	// saved variables for training with multiple number of clusters
	unsigned int *saved_nbClust; // [i*nbFilter+f] number of cluster for run i and filter f
	double *saved_bic;
	double *saved_weights; // saved_weight[c*nbMaxClusters*nbFilter+k*nbFilter+f] = weight of cluster k for run number c with filter f
	double *saved_mu; // saved_mu[c*nbMaxClusters*dim*nbFilter+k*dim*nbFilter+d*nbFilter+f]
	long double *saved_V; // saved_V[c*nbMaxClusters*dim2*nbFilter+k*dim2*nbFilter+d1*dim*nbFilter+d2*nbFilter+f]
	long double *saved_V_inv;
	long double *saved_det_V;
	bool SavedVariablesInitialized = false;

	// buffers for EM	
	double *D_i;
	double *C_di;
	double *buffer_di;
	double *E_d;

	// Parameters to put in FixedParameters
	unsigned int nbMaxClusters = 50;
	double tol_Lkh_EM = 1e-3;
	unsigned int MaxIter_EM = 500;
	double fac_elbow = .1; // reduction factor for considering that its is a real elbow

	bool IsKMeans = false;
	KMeans *MyKM;

public:
	// constructors
	GaussianMixtureModel();
	// methods
	void setDescriptors(Descriptors *D);
	void TrainModel(unsigned int &_nbClust, unsigned int &filter_value);
	void fitOptimalGMM(unsigned int &_nbClust_min, unsigned int &_nbClust_max);
	void InitFromKMeans(unsigned int &_nbClust, unsigned int &filter_value);
	void SaveVariables(unsigned int &current_nbClust, unsigned int &filter_value);
	void SetOptimalModel(unsigned int &opt_index, unsigned int &filter_value);
	void UpdateParams(unsigned int &filter_value);
	void ComputeLogLikelihood(unsigned int &filter_value);
	void ComputeBIC(unsigned int &filter_value);
	void PrintModelParams(std::string filename, unsigned int &filter_value);
	void EM(unsigned int &filter_value);
	long double Prob_Cluster(unsigned int index_cluster, unsigned int DescriptorIndex, unsigned int &filter_value);
	//void ReadModelParamFromDatabase();
	// getters
	unsigned int getNbClust(unsigned int &filter_value){ return nbClust[filter_value]; }
	double getBIC(unsigned int &filter_value){ return BIC[filter_value]; }
	double getLogLikelihood(unsigned int &filter_value){ return LogLikelihood[filter_value]; }
	KMeans *getKMeans(){ return MyKM; }
	long double *getCov(){ return V_inv; }
	double *getMu(){ return mu; }
	// destructor
	~GaussianMixtureModel();
};

#endif
